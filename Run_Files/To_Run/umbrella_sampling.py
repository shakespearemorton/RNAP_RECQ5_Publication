import openmm as mm
import openmm.app as app
import openmm.unit as unit
import math
import numpy as np
import sys
import HPSUrry
import HPSUrry_scaled
import subprocess
import argparse
from rigid import createRigidBodies


def prep_wham(N):
    kappa = 2  # kcal / nm
    file = 'separation.txt'
    # Load the data
    data = np.loadtxt(file)
    x = []
    c = []

    # Separate the data into windows and store the values
    for index in range(0, len(data), N):
        window = data[index:index + N]
        if len(window) > 0:
            x.append(window[0,1])  # Window minimum
            c.append(window[:, 2])  # Separation distance

    # Write the time series files and metadata file for WHAM
    metadata_lines = []

    for i, (x_val, c_val) in enumerate(zip(x, c)):
        # Define the time series file name for this window
        timeseries_filename = f'timeseries_window_{i + 1}.txt'
        np.savetxt(timeseries_filename, np.column_stack((np.arange(len(c_val)), c_val)), fmt='%d %.6f')
        
        # Prepare the line for the metadata file
        metadata_lines.append(f"{timeseries_filename} {x_val} {kappa}")

    # Write the metadata file
    with open('metadata.txt', 'w') as meta_file:
        meta_file.write('\n'.join(metadata_lines))

def generate_scaled(pser_charge,pstate):
    if pstate == 'p25':
        files = ['r5.txt','pctd.txt']
    elif pstate == 'p5':
        files = ['r5.txt','p5ctd.txt']
    elif pstate == 'none':
        files = ['r5.txt','ctd.txt']
    r5 = np.loadtxt(files[0],skiprows=68,max_rows=991)
    ctd = np.loadtxt(files[1],skiprows=68,max_rows=14)
    sri_start = 908
    sri_finish = 991
    sri = r5[sri_start:sri_finish]
    topology = app.Topology()
    positions = []
    real_sys = []

    #### Make SRI ####
    chain = topology.addChain()
    for i in range(len(sri)):
        residue = topology.addResidue(name="CG-residue", chain=chain)
        aa = sri[i,2]
        real_sys.append(sri[i,:])
        atom1 = topology.addAtom(name=HPSUrry_scaled.amino_acids[int(aa)], element=None, residue=residue)
        positions.append(sri[i,4:]/10)
    rigid_body = [np.arange(len(sri))]

    #### Make CTD ####
    chain = topology.addChain()
    atom2 = 0
    for i in range(len(ctd)):
        residue = topology.addResidue(name="CG-residue", chain=chain)
        aa = ctd[i,2]
        real_sys.append(ctd[i,:])
        atom1 = topology.addAtom(name=HPSUrry_scaled.amino_acids[int(aa)], element=None, residue=residue)
        if atom2 == 0:
            atom2 = atom1
        else:
            topology.addBond(atom1, atom2)
            atom2=atom1
        positions.append(ctd[i,4:]/10)
    positions *= unit.nanometer


    system = mm.System()
    b = 1000
    box_vec_a = np.array([b/10, 0, 0])*unit.nanometer
    box_vec_b = np.array([0, b/10, 0])*unit.nanometer
    box_vec_c = np.array([0, 0, b/10])*unit.nanometer
    system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)

    charges = []
    atom_types = []
    count = 0
    for atom in topology.atoms():
        if real_sys[count][3] == 4:
            charges.append(1.1)
        else:
            charges.append(HPSUrry_scaled.charge[atom.name])
        atom_types.append(HPSUrry_scaled.reverse[atom.name])
        system.addParticle(HPSUrry_scaled.mass[atom.name])
        count += 1

    r0 = 5 * unit.nanometer
    k = 2 * HPSUrry_scaled._kcal_to_kj
    restraint = mm.CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
    restraint.addGlobalParameter('r0', r0)
    restraint.addGlobalParameter('k',k)
    g1 = np.arange(len(sri))
    g2 = np.arange(len(ctd))+len(sri)
    restraint.addGroup(g1) # index will be 0
    restraint.addGroup(g2) # index will be 1
    restraint.addBond([0,1], []) #([g1_index,g2_index], [per-bond parameter list])
    restraint.setForceGroup(6)
    system.addForce(restraint)

    hbond = mm.HarmonicBondForce()


    for bond in topology.bonds():
        hbond.addBond(bond.atom1.index, bond.atom2.index, 0.38, 2000*HPSUrry_scaled._kcal_to_kj)
    hbond.setForceGroup(1)

    aa = len(HPSUrry_scaled.amino_acids)
    sigma_ah_map, lambda_ah_map = np.zeros((aa, aa)), np.zeros((aa, aa))
    for i in range(aa):
        for j in range(aa):
            sigma_ah_map[i,j] = np.mean([HPSUrry_scaled.size[HPSUrry_scaled.amino_acids[i+1]],HPSUrry_scaled.size[HPSUrry_scaled.amino_acids[j+1]]])/10
            lambda_ah_map[i,j] = np.mean([HPSUrry_scaled.hydropathy[HPSUrry_scaled.amino_acids[i+1]],HPSUrry_scaled.hydropathy[HPSUrry_scaled.amino_acids[j+1]]])
    lambda_ah_map = HPSUrry_scaled.mu*lambda_ah_map - HPSUrry_scaled.delta

    lj_at_cutoff = 4*HPSUrry_scaled.epsilon*((1/4)**12 - (1/4)**6)
    contacts = mm.CustomNonbondedForce(f'''energy;
            energy=(f1+f2-offset)*step(4*sigma_ah-r);
            offset=lambda_ah*{lj_at_cutoff};
            f1=(lj+(1-lambda_ah)*{HPSUrry_scaled.epsilon})*step(2^(1/6)*sigma_ah-r);
            f2=lambda_ah*lj*step(r-2^(1/6)*sigma_ah);
            lj=4*{HPSUrry_scaled.epsilon}*((sigma_ah/r)^12-(sigma_ah/r)^6);
            sigma_ah=sigma_ah_map(atom_type1, atom_type2);
            lambda_ah=lambda_ah_map(atom_type1, atom_type2);
            ''')
    n_atom_types = sigma_ah_map.shape[0]
    discrete_2d_sigma_ah_map = mm.Discrete2DFunction(n_atom_types, n_atom_types, sigma_ah_map.ravel().tolist())
    discrete_2d_lambda_ah_map = mm.Discrete2DFunction(n_atom_types, n_atom_types, lambda_ah_map.ravel().tolist())
    contacts.addTabulatedFunction('sigma_ah_map', discrete_2d_sigma_ah_map)
    contacts.addTabulatedFunction('lambda_ah_map', discrete_2d_lambda_ah_map)
    contacts.addPerParticleParameter('atom_type')

    for each in atom_types:
        contacts.addParticle([each-1])
    contacts.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    if HPSUrry_scaled.use_pbc:
        contacts.setNonbondedMethod(contacts.CutoffPeriodic)
    else:
        contacts.setNonbondedMethod(contacts.CutoffNonPeriodic)
    contacts.setCutoffDistance(4*np.amax(sigma_ah_map))
    contacts.setForceGroup(2)

    alpha = HPSUrry_scaled.NA*HPSUrry_scaled.EC**2/(4*np.pi*HPSUrry_scaled.VEP)
    ldby_value = HPSUrry_scaled.ldby.value_in_unit(unit.nanometer)
    alpha_value = alpha.value_in_unit(unit.kilojoule_per_mole*unit.nanometer)
    cutoff=3.5*unit.nanometer
    cutoff_value = cutoff.value_in_unit(unit.nanometer)
    elec = mm.CustomNonbondedForce(f'''energy;
        energy=q1*q2*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{HPSUrry_scaled.dielectric_water};
        offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
        ''')
    elec.addPerParticleParameter('q')
    for q in charges:
        if q == 1.1:
            q = 0
        elec.addParticle([q])
    elec.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    if HPSUrry_scaled.use_pbc:
        elec.setNonbondedMethod(elec.CutoffPeriodic)
    else:
        elec.setNonbondedMethod(elec.CutoffNonPeriodic)
    elec.setCutoffDistance(cutoff_value)
    elec.setForceGroup(3)

    elec_sri = mm.CustomNonbondedForce(f'''energy;
        energy=p1*p2*step(0-p1*p2)*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{HPSUrry_scaled.dielectric_water};
        offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
        ''')
    elec_sri.addPerParticleParameter('p')
    for q in charges:
        if abs(q) >= 2:
            elec_sri.addParticle([q])
        elif q == 1.1:
            elec_sri.addParticle([pser_charge])
        else:
            elec_sri.addParticle([0])
    elec_sri.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    if HPSUrry_scaled.use_pbc:
        elec_sri.setNonbondedMethod(elec.CutoffPeriodic)
    else:
        elec_sri.setNonbondedMethod(elec.CutoffNonPeriodic)
    elec_sri.setCutoffDistance(cutoff_value)
    elec_sri.setForceGroup(4)

    force = mm.CMMotionRemover()
    force.setForceGroup(5)

    system.addForce(hbond)
    system.addForce(contacts)
    system.addForce(elec)
    system.addForce(elec_sri)
    system.addForce(force)

    createRigidBodies(system, positions, rigid_body)

    with open('system.xml', 'w') as output:
        output.write(mm.XmlSerializer.serialize(system))

    with open('start.pdb','w') as f:
        app.PDBFile.writeFile(topology, positions, f)

def generate(pser_charge,pstate):
    if pstate == 'p25':
        files = ['r5.txt','pctd.txt']
    elif pstate == 'p5':
        files = ['r5.txt','p5ctd.txt']
    elif pstate == 'none':
        files = ['r5.txt','ctd.txt']
    r5 = np.loadtxt(files[0],skiprows=68,max_rows=991)
    ctd = np.loadtxt(files[1],skiprows=68,max_rows=14)
    sri_start = 908
    sri_finish = 991
    sri = r5[sri_start:sri_finish]
    topology = app.Topology()
    positions = []
    real_sys = []

    #### Make SRI ####
    chain = topology.addChain()
    for i in range(len(sri)):
        residue = topology.addResidue(name="CG-residue", chain=chain)
        aa = sri[i,2]
        real_sys.append(sri[i,:])
        atom1 = topology.addAtom(name=HPSUrry.amino_acids[int(aa)], element=None, residue=residue)
        positions.append(sri[i,4:]/10)
    rigid_body = [np.arange(len(sri))]

    #### Make CTD ####
    chain = topology.addChain()
    atom2 = 0
    for i in range(len(ctd)):
        residue = topology.addResidue(name="CG-residue", chain=chain)
        aa = ctd[i,2]
        real_sys.append(ctd[i,:])
        atom1 = topology.addAtom(name=HPSUrry.amino_acids[int(aa)], element=None, residue=residue)
        if atom2 == 0:
            atom2 = atom1
        else:
            topology.addBond(atom1, atom2)
            atom2=atom1
        positions.append(ctd[i,4:]/10)
    positions *= unit.nanometer


    system = mm.System()
    b = 1000
    box_vec_a = np.array([b/10, 0, 0])*unit.nanometer
    box_vec_b = np.array([0, b/10, 0])*unit.nanometer
    box_vec_c = np.array([0, 0, b/10])*unit.nanometer
    system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)

    charges = []
    atom_types = []
    count = 0
    for atom in topology.atoms():
        if real_sys[count][3] == 4:
            charges.append(1.1)
        else:
            charges.append(HPSUrry.charge[atom.name])
        atom_types.append(HPSUrry.reverse[atom.name])
        system.addParticle(HPSUrry.mass[atom.name])
        count += 1

    r0 = 5 * unit.nanometer
    k = 2 * HPSUrry._kcal_to_kj
    restraint = mm.CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
    restraint.addGlobalParameter('r0', r0)
    restraint.addGlobalParameter('k',k)
    g1 = np.arange(len(sri))
    g2 = np.arange(len(ctd))+len(sri)
    restraint.addGroup(g1) # index will be 0
    restraint.addGroup(g2) # index will be 1
    restraint.addBond([0,1], []) #([g1_index,g2_index], [per-bond parameter list])
    restraint.setForceGroup(6)
    system.addForce(restraint)

    hbond = mm.HarmonicBondForce()


    for bond in topology.bonds():
        hbond.addBond(bond.atom1.index, bond.atom2.index, 0.38, 2000*HPSUrry._kcal_to_kj)
    hbond.setForceGroup(1)

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

    for each in atom_types:
        contacts.addParticle([each-1])
    contacts.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    if HPSUrry.use_pbc:
        contacts.setNonbondedMethod(contacts.CutoffPeriodic)
    else:
        contacts.setNonbondedMethod(contacts.CutoffNonPeriodic)
    contacts.setCutoffDistance(4*np.amax(sigma_ah_map))
    contacts.setForceGroup(2)

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
    for q in charges:
        if q == 1.1:
            q = 0
        elec.addParticle([q])
    elec.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    if HPSUrry.use_pbc:
        elec.setNonbondedMethod(elec.CutoffPeriodic)
    else:
        elec.setNonbondedMethod(elec.CutoffNonPeriodic)
    elec.setCutoffDistance(cutoff_value)
    elec.setForceGroup(3)

    elec_sri = mm.CustomNonbondedForce(f'''energy;
        energy=p1*p2*step(0-p1*p2)*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{HPSUrry.dielectric_water};
        offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
        ''')
    elec_sri.addPerParticleParameter('p')
    for q in charges:
        if abs(q) >= 2:
            elec_sri.addParticle([q])
        elif q == 1.1:
            elec_sri.addParticle([pser_charge])
        else:
            elec_sri.addParticle([0])
    elec_sri.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    if HPSUrry.use_pbc:
        elec_sri.setNonbondedMethod(elec.CutoffPeriodic)
    else:
        elec_sri.setNonbondedMethod(elec.CutoffNonPeriodic)
    elec_sri.setCutoffDistance(cutoff_value)
    elec_sri.setForceGroup(4)

    force = mm.CMMotionRemover()
    force.setForceGroup(5)

    system.addForce(hbond)
    system.addForce(contacts)
    system.addForce(elec)
    system.addForce(elec_sri)
    system.addForce(force)

    createRigidBodies(system, positions, rigid_body)

    with open('system.xml', 'w') as output:
        output.write(mm.XmlSerializer.serialize(system))

    with open('start.pdb','w') as f:
        app.PDBFile.writeFile(topology, positions, f)

def simulate(output_interval,steps,hysterisis):
    
    class COMSeparationReporter():
        def __init__(self, file, reportInterval, g1_indices, g2_indices):
            self._out = open(file, 'w')
            self._reportInterval = reportInterval
            self.g1_indices = g1_indices
            self.g2_indices = g2_indices
            self.saved_data = None


        def __del__(self):
            self._out.close()

        def describeNextReport(self, simulation):
            if self.saved_data is None:
                steps = self._reportInterval - simulation.currentStep % self._reportInterval - 1
            else:
                steps = 1
            return (steps, True, False, False, False)

        def report(self, simulation, state):
            positions = state.getPositions(asNumpy=True)

            # Calculate the masses and centers of mass for each group
            g1_masses = self.getMasses(system, self.g1_indices)
            g2_masses = self.getMasses(system, self.g2_indices)
            g1_com = self.getCenterOfMass(positions[self.g1_indices], g1_masses)
            g2_com = self.getCenterOfMass(positions[self.g2_indices], g2_masses)

            # Calculate the separation distance
            separation = np.linalg.norm(g1_com - g2_com)
            r = simulation.context.getParameter('r0')
            self._out.write(f'{simulation.currentStep}\t{r}\t{separation}\n')
            self._out.flush()

        def getMasses(self, system, indices):
            masses = [system.getParticleMass(i).value_in_unit(unit.dalton) for i in indices]
            return np.array(masses)

        def getCenterOfMass(self, positions, masses):
            total_mass = np.sum(masses)
            com = np.sum(positions * masses[:, None], axis=0) / total_mass
            return com
        
        
    temperature = 300
    box_a = 300
    box_b = 300
    box_c = 300

    system_xml = 'system.xml'
    with open(system_xml, 'r') as f:
        system = mm.XmlSerializer.deserialize(f.read())
    box_vec_a = np.array([box_a, 0, 0])*unit.nanometer
    box_vec_b = np.array([0, box_b, 0])*unit.nanometer
    box_vec_c = np.array([0, 0, box_c])*unit.nanometer
    system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)
    top = app.PDBFile('start.pdb').getTopology()
    print(top)

    init_coord = np.array(app.PDBFile('start.pdb').getPositions().value_in_unit(unit.nanometer))
    init_coord -= np.mean(init_coord, axis=0)
    init_coord += 0.5*np.array([box_a, box_b, box_c])

    start_temperature = 150
    collision = 1/unit.picosecond
    timestep = 10*unit.femtosecond
    integrator = mm.LangevinMiddleIntegrator(start_temperature*unit.kelvin, collision, timestep)
    platform_name = 'CUDA'
    platform = mm.Platform.getPlatformByName(platform_name)
    simulation = app.Simulation(top, system, integrator, platform)
    simulation.context.setPositions(init_coord)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(start_temperature*unit.kelvin)
    n_iterations = 10
    state_data_reporter = app.StateDataReporter(sys.stdout, output_interval, step=True, time=True, potentialEnergy=True,
                                            kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
    comrepo = COMSeparationReporter('separation.txt', output_interval,np.arange(97-14),np.arange(14)+83)

    for temperature_i in np.linspace(start_temperature, temperature, n_iterations):
        integrator.setTemperature(temperature_i*unit.kelvin)
        simulation.step(1000)
    if hysterisis == 0:
        r_range = np.linspace(5,1,25)
        simulation.context.setParameter('r0',r_range[0])
    else:
        r_range = np.linspace(1,5,25)
        simulation.context.setParameter('r0',r_range[0])
    simulation.step(steps/10)
    simulation.reporters.append(comrepo)
    simulation.reporters.append(state_data_reporter)
    for i in r_range:
        simulation.context.setParameter('r0',i)
        simulation.step(steps)
        print(i)


def main():
    """Main function to set up and run the simulation."""
    
    # Parse command-line arguments for temperature and model
    parser = argparse.ArgumentParser(description="Run protein simulation with specified parameters.")
    parser.add_argument("--hysterisis", type=float, default=0, help="0 - Far to Close, 1 - Close to Far")
    parser.add_argument("--scaled", type=float, default=0, help="0 - No scaling, 1 - Folded domain scaled by 0.8")
    parser.add_argument("--out_int", type=float, default=1e4, help="Output Interval")
    parser.add_argument("--window_time", type=float, default=1e7, help="How much time to spend in each window")
    parser.add_argument("--pser_charge", type=float, default=2.75, help="Charge placed on the SRI domain")
    parser.add_argument("--pstate", type=str, default='p25', help="Phosphorylation state (p25, p5, none)")

    args = parser.parse_args()

    if args.scaled == 0:
        generate(args.pser_charge,args.pstate)
    elif args.scaled == 1:
        generate_scaled(args.pser_charge,args.pstate)
    simulate(args.out_int,args.window_time,args.hysterisis)
    print(args.window_time,args.out_int)
    N = int(args.window_time/args.out_int)
    prep_wham(N)
    direc = '/storage/brno14-ceitec/shared/softmatter/shakespearem/wham/wham/wham'
    save_file = f'{args.pstate}_{args.pser_charge}_{args.scaled}_{args.hysterisis}.txt'
    command = [direc, "0.8", "6.2", "250", "1e-6", "300", "0", "metadata.txt", save_file, "1000", "12345"]


    result = subprocess.run(command, capture_output=True, text=True)

    

if __name__ == "__main__":
    main()
