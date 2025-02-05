import openmm as mm
import openmm.app as app
import openmm.unit as unit
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.lib.distances import distance_array
from scipy.spatial.transform import Rotation as R
import MDAnalysis as mda
import pandas as pd
import math
import os
import numpy as np
import itertools
import sys
import HPSUrry
import string
import numpy.linalg as lin
from itertools import combinations

class Disordered_Life:
    def __init__(self, experiment_name, proteins, num, style, runtime, sampling, temperature, boxsize, pbc=True,special_rule = False):
        global styles
        styles={0: 'compressed', 1: 'uncompressed', 2: 'slab'}
        self.experiment_name = experiment_name
        self.proteins = proteins
        self.num = num
        self.style = style
        self.runtime = int(runtime)
        self.boxsize = boxsize
        self.sampling = int(sampling)
        self.temperature = temperature
        self.pbc = pbc
        self.status = 'Generating'
        self.special = special_rule
        
    def execute(self):
        ### CREATE PDB FILE FOR SYSTEM ###
        if os.path.exists(f'{self.experiment_name}.pdb'):
            self.status = ['Created']
        else:
            self.create_system()
            self.pdbfix()
            self.status = ['Created']
            
        if os.path.exists(f'{self.experiment_name}_prev.xml'):
            self.status.append('Re-Running')
            self.rerun()

        elif os.path.exists(f'{self.experiment_name}_equil.xml'):
            self.status.append('Compressed')
            self.status.append('Running')
            if self.style == 2:
                self.run_slab()
            elif self.style == 1:
                self.run()
            else:
                self.run_box()
        else:
            if self.style == 1:
                self.status.append('Running')
                self.run()
            elif self.style == 0:
                self.compress_system()
            elif self.style == 2:
                self.compress_system()
        
        self.status.append('Completed')
        self.readout()
                
    def readout(self):
        print('\n')
        print(f'###   {self.experiment_name}   ###')
        print(f'Status:                {self.status}')
        print(f'Proteins of study:     {self.proteins}')
        print(f'Number of each:        {self.num}')
        print(f'System style:          {styles[self.style]}')
        print(f'PBC:                   {self.pbc}')
        print(f'Box Size:              {self.boxsize} nm')
        print(f'Runtime:               {self.runtime}')
        print('\n')
        
    def pdbfix(self):
        alphabet = string.ascii_uppercase
        combinations = []
        for i in range(50):
            for r in range(1):
                combinations.extend([''.join(p) for p in itertools.product(alphabet, repeat=1)])
        topology = app.PDBFile(f"{self.experiment_name}.pdb").getTopology()

        def format_pdb_line(record_name, atom_number, res_name, atom_name, chain_id, res_number, x, y, z, occupancy, temp_factor, element):
            return (f"{record_name:6}{atom_number:5d} {res_name:4} {atom_name} {chain_id}{res_number:4d}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}          {element:>2}")

        def generate_pdb_file(crystal_dimensions, atoms, filename):
            pdb_content = [crystal_dimensions + "\n"]
            atom_number = 1  # Start atom numbering from 1
            for atom in atoms:
                if atom_number > 99999:
                    atom_number = 1  # Restart numbering if exceeds 99999
                pdb_line = format_pdb_line("ATOM ", atom_number, atom[1], "CA-", atom[2], atom[3], atom[4], atom[5], atom[6], 1.00, 0.00, 'C')
                pdb_content.append(pdb_line + "\n")
                atom_number += 1  # Increment atom number

            with open(filename, "w") as f:
                f.write("".join(pdb_content))

        crystal_dimensions = "CRYST1 1870.000 1870.000 1870.000  90.00  90.00  90.00 P 1           1"
        data = []
        c = 0
        for chain in topology.chains():
            a = 0
            for atom in chain.atoms():
                # Update coordinates appropriately
                data.append([atom.index + 1, atom.name, combinations[c % len(combinations)], a + 1, 0.0, 0.0, 0.0])
                a += 1
            c += 1
        generate_pdb_file(crystal_dimensions, data, f"{self.experiment_name}_mda.pdb")

    def create_system(self):
        success = False
        while success == False:
            success = self.add_proteins()
            self.boxsize += 50
        self.boxsize /= 10
        data = np.loadtxt(f'{self.experiment_name}.txt',skiprows=1)
        topology = app.Topology()
        c = np.inf
        for i in range(len(data)):
            if data[i,1] != c: # START A NEW CHAIN
                c = data[i,1]
                chain = topology.addChain()
                residue = topology.addResidue(name="CG-residue", chain=chain)
                a1_name = HPSUrry.amino_acids[data[i,2]]
                atom1 = topology.addAtom(name = a1_name, element=None, residue=residue)
            else:  #CONTINUE A CHAIN
                a2_name = HPSUrry.amino_acids[data[i,2]]
                residue = topology.addResidue(name="CG-residue", chain=chain)
                atom2 = topology.addAtom(name = a2_name, element=None, residue=residue)
                topology.addBond(atom1, atom2)
                atom1 = atom2
                a1_name = a2_name
                
        ### CREATE SYSTEM ###
        print(topology)
        positions = (data[:,3:6]/10 * unit.nanometer)
        system = mm.System()
        for atom in topology.atoms():
            system.addParticle(HPSUrry.mass[atom.name])

        ### CREATE RIGID BODIES ###
        in_rigid = False
        rigid_bodies = []
        i = 0
        cs = True
        for atom in topology.atoms():
            if cs == True:
                c_old = atom.residue.chain.index
                cs = False
            if len(atom.name) > 3:
                if in_rigid == False:
                    start = i
                    in_rigid = True
            else:
                if in_rigid == True:
                    rigid_bodies.append(np.arange(start,i-1))
                    in_rigid = False
            c_new = atom.residue.chain.index
            if c_new != c_old:
                if in_rigid == True:
                    rigid_bodies.append(np.arange(start,i))
                    start = i
            c_old = c_new
            i+=1
        if in_rigid == True:
            last = np.arange(start,i)
            if rigid_bodies[-1][-1] != last[-1]:
                rigid_bodies.append(last)

        ### MAKE HYDROPHOBIC ###
        aa = len(HPSUrry.amino_acids)
        sigma_ah_map, lambda_ah_map = np.zeros((aa, aa)), np.zeros((aa, aa))
        for i in range(aa):
            for j in range(aa):
                sigma_ah_map[i,j] = np.mean([HPSUrry.size[HPSUrry.amino_acids[i+1]],HPSUrry.size[HPSUrry.amino_acids[j+1]]])/10
                lambda_ah_map[i,j] = np.mean([HPSUrry.hydropathy[HPSUrry.amino_acids[i+1]],HPSUrry.hydropathy[HPSUrry.amino_acids[j+1]]])
        lambda_ah_map = HPSUrry.mu*lambda_ah_map - HPSUrry.delta
        lj_at_cutoff = 4*HPSUrry.epsilon*((1/4)**12 - (1/4)**6)
        contacts = mm.CustomNonbondedForce(f'''energy;
                   energy=(f1+f2-offset)*step(4*sigma_ah-r);
                   offset=lambda_ah*{lj_at_cutoff};
                   f1=(lj+(1-lambda_ah)*{HPSUrry.epsilon})*step(2^(1/6)*sigma_ah-r);
                   f2=lambda_ah*lj*step(r-2^(1/6)*sigma_ah);
                   lj=4*{HPSUrry.epsilon}*((sigma_ah/r)^12-(sigma_ah/r)^6);
                   sigma_ah=sigma_ah_map(atom_type1, atom_type2);
                   lambda_ah=lambda_ah_map(atom_type1, atom_type2);
                   ''')
        n_atom_types = sigma_ah_map.shape[0]
        discrete_2d_sigma_ah_map = mm.Discrete2DFunction(n_atom_types, n_atom_types, sigma_ah_map.ravel().tolist())
        discrete_2d_lambda_ah_map = mm.Discrete2DFunction(n_atom_types, n_atom_types, lambda_ah_map.ravel().tolist())
        contacts.addTabulatedFunction('sigma_ah_map', discrete_2d_sigma_ah_map)
        contacts.addTabulatedFunction('lambda_ah_map', discrete_2d_lambda_ah_map)
        contacts.addPerParticleParameter('atom_type')
        for atom in topology.atoms():
            contacts.addParticle([HPSUrry.reverse[atom.name]-1])
        
        if self.pbc:
            contacts.setNonbondedMethod(contacts.CutoffPeriodic)
        else:
            contacts.setNonbondedMethod(contacts.CutoffNonPeriodic)
        contacts.setCutoffDistance(4*np.amax(sigma_ah_map))
        contacts.setForceGroup(2)

        ### MAKE ELECTROSTATIC ###
        alpha = HPSUrry.NA*HPSUrry.EC**2/(4*np.pi*HPSUrry.VEP)
        ldby_value = HPSUrry.ldby.value_in_unit(unit.nanometer)
        alpha_value = alpha.value_in_unit(unit.kilojoule_per_mole*unit.nanometer)
        cutoff=3.5*unit.nanometer
        cutoff_value = cutoff.value_in_unit(unit.nanometer)
        elec = mm.CustomNonbondedForce(f'''energy;
               energy=q1*q2*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{HPSUrry.dielectric_water};
               offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
               ''')
        elec.addPerParticleParameter('q')
        for atom in topology.atoms():
            elec.addParticle([HPSUrry.charge[atom.name]])
        if self.pbc:
            elec.setNonbondedMethod(elec.CutoffPeriodic)
        else:
            elec.setNonbondedMethod(elec.CutoffNonPeriodic)
        elec.setCutoffDistance(cutoff_value)
        elec.setForceGroup(3)

        if self.special == True:
            print('special in')
            elec_sri = mm.CustomNonbondedForce(f'''energy;
            energy=p1*p2*step(0-p1*p2)*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{HPSUrry.dielectric_water};
            offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
            ''')
            elec_sri.addPerParticleParameter('p')
            sri = [920,938,939,943,947]
            r5_kix = 529
            rna_kix = 1247
            connect = False
            for chain in topology.chains():
                sri_charge = []
                for atom in chain.atoms():
                    if HPSUrry.charge[atom.name] == -2:
                        sri_charge.append(-2)
                    else:
                        sri_charge.append(0)
                if len(sri_charge) == 991:
                    for bead in sri:
                        sri_charge[bead] = 2.75
                    if connect == True:
                        connect = False
                        count = 0
                        for atom in chain.atoms():
                            if count == r5_kix:
                                topology.addBond(atom1, atom)
                                print('bond formed')
                            count += 1
                if len(sri_charge) == 4380:
                    connect = True
                    count = 0
                    for atom in chain.atoms():
                        if count == rna_kix:
                            atom1 = atom
                        count += 1
                print(len(sri_charge))
                for q in sri_charge:
                    elec_sri.addParticle([q])
            elec_sri.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
            if HPSUrry.use_pbc:
                elec_sri.setNonbondedMethod(elec.CutoffPeriodic)
            else:
                elec_sri.setNonbondedMethod(elec.CutoffNonPeriodic)
            elec_sri.setCutoffDistance(cutoff_value)
            elec_sri.setForceGroup(4)
            system.addForce(elec_sri)

        ### REMOVE CMM MOTION ###
        force = mm.CMMotionRemover()
        force.setForceGroup(5)

        ### MAKE BONDS ###
        hbond = mm.HarmonicBondForce()
        for bond in topology.bonds():
            hbond.addBond(bond.atom1.index, bond.atom2.index, 0.38, 2000*HPSUrry._kcal_to_kj)
        hbond.setForceGroup(1)

        ### GENERATE SYSTEM ###
        contacts.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
        elec.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
        system.addForce(hbond)
        system.addForce(contacts)
        system.addForce(elec)
        system.addForce(force)

        box_vec_a = np.array([self.boxsize, 0, 0])*unit.nanometer
        box_vec_b = np.array([0, self.boxsize, 0])*unit.nanometer
        box_vec_c = np.array([0, 0, self.boxsize])*unit.nanometer
        system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)
        topology.setUnitCellDimensions([self.boxsize, self.boxsize, self.boxsize])
        createRigidBodies(system, positions, rigid_bodies)
        with open(f'{self.experiment_name}_system.xml', 'w') as output:
            output.write(mm.XmlSerializer.serialize(system))
        with open(f'{self.experiment_name}.pdb','w') as f:
            app.PDBFile.writeFile(topology, positions, f)

    def add_proteins(self,atoms=None):
        total_num = 0
        for i in range(len(self.proteins)):
            method = 'FastNS'
            max_n_attempts=10000
            try:
                if atoms==None:
                    atoms = pd.DataFrame()
            except:
                pass
            data = np.loadtxt(self.proteins[i]) #residue is on index 2
            new_atoms = pd.DataFrame(data)
            new_atoms.columns =['chain','residue','x','y','z','rigid']
            for index, row in new_atoms.iterrows():
                if row['rigid'] == 1:
                    new_atoms.at[index, 'residue'] += 21
            new_coord = new_atoms[['x', 'y', 'z']].to_numpy()
            new_coord -= np.mean(new_coord, axis=0)
            cutoff = 10.0
            count_n_mol = 0
            count_n_attempts = 0
            dim = np.array([self.boxsize, self.boxsize, self.boxsize, 90.0, 90.0, 90.0])
            while (count_n_mol < self.num[i]) and (count_n_attempts < max_n_attempts):
                # get a random rotation
                rotate = R.random()
                new_coord_i = rotate.apply(new_coord)
                # get a random translation
                translate = np.random.uniform(0, 1, 3)*np.array([self.boxsize, self.boxsize, self.boxsize])
                new_coord_i += translate
                if len(atoms.index) == 0:
                    new_atoms_i = new_atoms.copy()
                    new_atoms_i[['x', 'y', 'z']] = new_coord_i
                    atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
                    count_n_mol += 1
                    if method == 'FastNS':
                        coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
                        grid_search = FastNS(cutoff, coord, dim, pbc=True)
                else:
                    flag = False
                    if method == 'distance_array':
                        coord = atoms[['x', 'y', 'z']].to_numpy()
                        d = distance_array(coord, new_coord_i, dim)
                        if np.amin(d) >= cutoff:
                            flag = True
                    elif method == 'FastNS':
                        results = grid_search.search(new_coord_i.astype(np.float32))
                        if len(results.get_pair_distances()) == 0:
                            flag = True
                    if flag:
                        new_atoms_i = new_atoms.copy()
                        new_atoms_i['chain'] = new_atoms_i['chain']+np.max(atoms['chain'])+1
                        new_atoms_i[['x', 'y', 'z']] = new_coord_i
                        atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
                        count_n_mol += 1
                        if method == 'FastNS':
                            coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
                            grid_search = FastNS(cutoff, coord, dim, pbc=True)
                count_n_attempts += 1
            total_num += count_n_mol
        if total_num == np.sum(self.num):
            print('Successfully Created.')
            atoms.to_csv(f'{self.experiment_name}.txt', sep ='\t')
            return True
        else:
            print(f'Could not successfully insert {self.num[i]} molecules in {count_n_attempts} attempts.')
            print(f'Only added {count_n_mol} molecules. Increasing the box size...')
            return False

    def compress_system(self):
        ### LOAD SYSTEM ###
        top = app.PDBFile(f'{self.experiment_name}.pdb').getTopology()
        init_coord = app.PDBFile(f'{self.experiment_name}.pdb').getPositions()
        system_xml = f'{self.experiment_name}_system.xml'
        with open(system_xml, 'r') as f:
            system = mm.XmlSerializer.deserialize(f.read())

        ### INPUTS ###
        pressure = 1*unit.bar
        temperature = 300*unit.kelvin
        output_dcd = f'{self.experiment_name}_compress.dcd'
        output_interval = 5e5
        n_steps = 1e7
        friction_coeff = 1/unit.picosecond
        timestep = 10*unit.femtosecond
        barostat = mm.MonteCarloBarostat(pressure, temperature)
        system.addForce(barostat)
        integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)

        ### COMPUTING ###
        properties = {'Precision': 'mixed'}
        platform_name = 'CUDA'
        platform = mm.Platform.getPlatformByName(platform_name)
        
        ### SETUP SIMULATION ###
        simulation = app.Simulation(top, system, integrator, platform,properties)
        simulation.context.setPositions(init_coord)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(temperature)

        ### REPORTERS ###
        dcd_reporter = app.DCDReporter(output_dcd, output_interval, enforcePeriodicBox=True)
        state_data_reporter = app.StateDataReporter(f'{self.experiment_name}_compress.csv', output_interval, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
        #simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(state_data_reporter)
        
        ### RUN SIMULATION ###
        simulation.step(n_steps)

        ### GET STATUS ###
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, enforcePeriodicBox=True)
        with open(f'{self.experiment_name}_equil.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))  # Implementation goes here

    def run_slab(self):
        ### LOAD SYSTEM ###
        system_xml = f'{self.experiment_name}_system.xml'
        with open(system_xml, 'r') as f:
            system = mm.XmlSerializer.deserialize(f.read())
        npt_final_state_xml = f'{self.experiment_name}_equil.xml'
        with open(npt_final_state_xml, 'r') as f:
            npt_final_state = mm.XmlSerializer.deserialize(f.read())
        xyz = list((npt_final_state.getPeriodicBoxVectors(asNumpy=True)))
        xyz[-1]*=10
        top = app.PDBFile(f'{self.experiment_name}.pdb').getTopology()
        init_coord = np.array(npt_final_state.getPositions().value_in_unit(unit.nanometer))
        
        ### INPUTS
        temperature = self.temperature *unit.kelvin
        friction_coeff = 1/unit.picosecond
        timestep = 10*unit.femtosecond
        
        ### COMPUTING ###
        properties = {'Precision': 'mixed'}
        platform_name = 'CUDA'
        platform = mm.Platform.getPlatformByName(platform_name)
        
        ### SETUP SIMULATION ###
        integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
        simulation = app.Simulation(top, system, integrator, platform, properties)
        simulation.context.setState(npt_final_state)
        simulation.context.setPeriodicBoxVectors(xyz[0],xyz[1],xyz[2])
        
        ### REPORTERS ###
        dcd_reporter = app.DCDReporter(f'{self.experiment_name}.dcd', self.sampling, enforcePeriodicBox=True)
        state_data_reporter = app.StateDataReporter(f'{self.experiment_name}.csv', self.sampling, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
        simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(state_data_reporter)
        
        ### RUN SIMULATION ###
        if self.runtime == 0:
            simulation.runForClockTime(20)
        else:
            simulation.step(self.runtime)
        simulation.saveCheckpoint(f'{self.experiment_name}.check')
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,enforcePeriodicBox=True)
        with open(f'{self.experiment_name}_prev.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))

    def run_box(self):
        ### LOAD SYSTEM ###
        system_xml = f'{self.experiment_name}_system.xml'
        with open(system_xml, 'r') as f:
            system = mm.XmlSerializer.deserialize(f.read())
        npt_final_state_xml = f'{self.experiment_name}_equil.xml'
        with open(npt_final_state_xml, 'r') as f:
            npt_final_state = mm.XmlSerializer.deserialize(f.read())
        xyz = list((npt_final_state.getPeriodicBoxVectors(asNumpy=True)))
        xyz[-1] *= 10
        xyz[-2] *= 10
        xyz[-3] *= 10
        top = app.PDBFile(f'{self.experiment_name}.pdb').getTopology()
        box = np.asarray([np.max(xyz[0])._value, np.max(xyz[1])._value, np.max(xyz[2])._value])*unit.nanometer
        init_coord = npt_final_state.getPositions()
        init_coord += box*0.5
        
        ### INPUTS
        temperature = self.temperature *unit.kelvin
        friction_coeff = 1/unit.picosecond
        timestep = 10*unit.femtosecond
        
        ### COMPUTING ###
        properties = {'Precision': 'mixed'}
        platform_name = 'CUDA'
        platform = mm.Platform.getPlatformByName(platform_name)
        
        ### SETUP SIMULATION ###
        integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
        simulation = app.Simulation(top, system, integrator, platform, properties)
        simulation.context.setState(npt_final_state)
        simulation.context.setPeriodicBoxVectors(xyz[0],xyz[1],xyz[2])
        simulation.context.setPositions(init_coord)
        
        ### REPORTERS ###
        dcd_reporter = app.DCDReporter(f'{self.experiment_name}.dcd', self.sampling, enforcePeriodicBox=True)
        state_data_reporter = app.StateDataReporter(f'{self.experiment_name}.csv', self.sampling, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
        simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(state_data_reporter)
        
        ### RUN SIMULATION ###
        if self.runtime == 0:
            simulation.runForClockTime(20)
        else:
            simulation.step(self.runtime)
        simulation.saveCheckpoint(f'{self.experiment_name}.check')
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,enforcePeriodicBox=True)
        with open(f'{self.experiment_name}_prev.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))
        
    def run(self):
        ### LOAD SYSTEM ###
        top = app.PDBFile(f'{self.experiment_name}.pdb').getTopology()
        init_coord = app.PDBFile(f'{self.experiment_name}.pdb').getPositions()
        system_xml = f'{self.experiment_name}_system.xml'
        with open(system_xml, 'r') as f:
            system = mm.XmlSerializer.deserialize(f.read())
        xyz = list((system.getPeriodicBoxVectors(asNumpy=True)))
        top = app.PDBFile(f'{self.experiment_name}.pdb').getTopology()
        box = np.asarray([np.max(xyz[0])._value, np.max(xyz[1])._value, np.max(xyz[2])._value])*unit.nanometer
        init_coord += 0.5*box

        ### INPUTS ###
        temperature = 300*unit.kelvin
        friction_coeff = 1/unit.picosecond
        timestep = 10*unit.femtosecond
        integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)

        ### COMPUTING ###
        properties = {'Precision': 'mixed'}
        platform_name = 'CUDA'
        platform = mm.Platform.getPlatformByName(platform_name)
        
        ### SETUP SIMULATION ###
        simulation = app.Simulation(top, system, integrator, platform,properties)
        simulation.context.setPositions(init_coord)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(temperature)

        ### REPORTERS ###
        dcd_reporter = app.DCDReporter(f'{self.experiment_name}.dcd', self.sampling, enforcePeriodicBox=True)
        state_data_reporter = app.StateDataReporter(f'{self.experiment_name}.csv', self.sampling, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
        simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(state_data_reporter)
        
        ### RUN SIMULATION ###
        if self.runtime == 0:
            simulation.runForClockTime(20)
        else:
            simulation.step(self.runtime)
        simulation.saveCheckpoint(f'{self.experiment_name}.check')
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,enforcePeriodicBox=True)
        with open(f'{self.experiment_name}_prev.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))


    def rerun(self):
        ### LOAD SYSTEM ###
        system_xml = f'{self.experiment_name}_system.xml'
        with open(system_xml, 'r') as f:
            system = mm.XmlSerializer.deserialize(f.read())
        previous_state = f'{self.experiment_name}_prev.xml'
        with open(previous_state, 'r') as f:
            prev = mm.XmlSerializer.deserialize(f.read())
        top = app.PDBFile(f'{self.experiment_name}.pdb').getTopology()
        
        ### INPUTS
        temperature = self.temperature *unit.kelvin
        friction_coeff = 1/unit.picosecond
        timestep = 10*unit.femtosecond
        
        ### COMPUTING ###
        properties = {'Precision': 'mixed'}
        platform_name = 'CUDA'
        platform = mm.Platform.getPlatformByName(platform_name)
        
        ### SETUP SIMULATION ###
        integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
        simulation = app.Simulation(top, system, integrator, platform, properties)
        simulation.context.setState(prev)

        ### REPORTERS ###
        dcd_reporter = app.DCDReporter(f'{self.experiment_name}.dcd', self.sampling, enforcePeriodicBox=True,append=True)
        state_data_reporter = app.StateDataReporter(f'{self.experiment_name}.csv', self.sampling, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, speed=True,append=True)
        simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(state_data_reporter)
        ### RUN SIMULATION ###
        if os.path.isfile(f'{self.experiment_name}.check'):
            simulation.loadCheckpoint(f'{self.experiment_name}.check')
        if self.runtime == 0:
            simulation.runForClockTime(20)
        else:
            simulation.step(self.runtime)
        simulation.saveCheckpoint(f'{self.experiment_name}.check')
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,enforcePeriodicBox=True)
        with open(f'{self.experiment_name}_prev.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))
            
def createRigidBodies(system, positions, bodies):
    __author__ = "Peter Eastman"
    __version__ = "1.0"
    """Modify a System to turn specified sets of particles into rigid bodies.

    For every rigid body, four particles are selected as "real" particles whose positions are integrated.
    Constraints are added between them to make them move as a rigid body.  All other particles in the body
    are then turned into virtual sites whose positions are computed based on the "real" particles.

    Because virtual sites are massless, the mass properties of the rigid bodies will be slightly different
    from the corresponding sets of particles in the original system.  The masses of the non-virtual particles
    are chosen to guarantee that the total mass and center of mass of each rigid body exactly match those of
    the original particles.  The moment of inertia will be similar to that of the original particles, but
    not identical.

    Care is needed when using constraints, since virtual particles cannot participate in constraints.  If the
    input system includes any constraints, this function will automatically remove ones that connect two
    particles in the same rigid body.  But if there is a constraint between a particle in a rigid body and
    another particle not in that body, it will likely lead to an exception when you try to create a context.

    Parameters:
     - system (System) the System to modify
     - positions (list) the positions of all particles in the system
     - bodies (list) each element of this list defines one rigid body.  Each element should itself be a list
       of the indices of all particles that make up that rigid body.
    """
    # Remove any constraints involving particles in rigid bodies.

    for i in range(system.getNumConstraints()-1, -1, -1):
        p1, p2, distance = system.getConstraintParameters(i)
        if (any(p1 in body and p2 in body for body in bodies)):
            system.removeConstraint(i)

    # Loop over rigid bodies and process them.

    for particles in bodies:
        if len(particles) < 5:
            # All the particles will be "real" particles.

            realParticles = particles
            realParticleMasses = [system.getParticleMass(i) for i in particles]
        else:
            # Select four particles to use as the "real" particles.  All others will be virtual sites.

            pos = [positions[i] for i in particles]
            mass = [system.getParticleMass(i) for i in particles]
            cm = unit.sum([p*m for p, m in zip(pos, mass)])/unit.sum(mass)
            r = [p-cm for p in pos]
            avgR = unit.sqrt(unit.sum([unit.dot(x, x) for x in r])/len(particles))
            rank = sorted(range(len(particles)), key=lambda i: abs(unit.norm(r[i])-avgR))
            for p in combinations(rank, 4):
                # Select masses for the "real" particles.  If any is negative, reject this set of particles
                # and keep going.

                matrix = np.zeros((4, 4))
                for i in range(4):
                    particleR = r[p[i]].value_in_unit(unit.nanometers)
                    matrix[0][i] = particleR[0]
                    matrix[1][i] = particleR[1]
                    matrix[2][i] = particleR[2]
                    matrix[3][i] = 1.0
                rhs = np.array([0.0, 0.0, 0.0, unit.sum(mass).value_in_unit(unit.amu)])
                weights = lin.solve(matrix, rhs)
                if all(w > 0.0 for w in weights):
                    # We have a good set of particles.

                    realParticles = [particles[i] for i in p]
                    realParticleMasses = [float(w) for w in weights]*unit.amu
                    break

        # Set particle masses.

        for i, m in zip(realParticles, realParticleMasses):
            system.setParticleMass(i, m)

        # Add constraints between the real particles.

        for p1, p2 in combinations(realParticles, 2):
            distance = unit.norm(positions[p1]-positions[p2])
            key = (min(p1, p2), max(p1, p2))
            system.addConstraint(p1, p2, distance)

        # Select which three particles to use for defining virtual sites.

        bestNorm = 0
        for p1, p2, p3 in combinations(realParticles, 3):
            d12 = (positions[p2]-positions[p1]).value_in_unit(unit.nanometer)
            d13 = (positions[p3]-positions[p1]).value_in_unit(unit.nanometer)
            crossNorm = unit.norm((d12[1]*d13[2]-d12[2]*d13[1], d12[2]*d13[0]-d12[0]*d13[2], d12[0]*d13[1]-d12[1]*d13[0]))
            if crossNorm > bestNorm:
                bestNorm = crossNorm
                vsiteParticles = (p1, p2, p3)

        # Create virtual sites.

        d12 = (positions[vsiteParticles[1]]-positions[vsiteParticles[0]]).value_in_unit(unit.nanometer)
        d13 = (positions[vsiteParticles[2]]-positions[vsiteParticles[0]]).value_in_unit(unit.nanometer)
        cross = mm.Vec3(d12[1]*d13[2]-d12[2]*d13[1], d12[2]*d13[0]-d12[0]*d13[2], d12[0]*d13[1]-d12[1]*d13[0])
        matrix = np.zeros((3, 3))
        for i in range(3):
            matrix[i][0] = d12[i]
            matrix[i][1] = d13[i]
            matrix[i][2] = cross[i]
        for i in particles:
            if i not in realParticles:
                system.setParticleMass(i, 0)
                rhs = np.array((positions[i]-positions[vsiteParticles[0]]).value_in_unit(unit.nanometer))
                weights = lin.solve(matrix, rhs)
                system.setVirtualSite(i, mm.OutOfPlaneSite(vsiteParticles[0], vsiteParticles[1], vsiteParticles[2], weights[0], weights[1], weights[2]))
