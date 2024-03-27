import openmm as mm
import openmm.app as app
import openmm.unit as unit
import math
import numpy as np
from sys import stdout
import aminos
from adder_drop import add
from rigid import createRigidBodies
print(mm.__file__)
files = ['r5.txt','ppol.txt']
rigid = [[np.arange(12,250),np.arange(260,321),np.arange(324,452),np.arange(491,613),np.arange(708,790),np.arange(908,991)],[np.arange(0,3898)]]
N = 10
n_mol = [5*N,N]
r1 = [0]*(n_mol[0]-n_mol[1])
r2 = [1,0]*n_mol[1]
order = np.concatenate((r1,r2))
print(order)
topology = app.Topology()
box = 700
add(files,n_mol,box_a=box,box_b=box,box_c=box,output='output.txt')
data = np.loadtxt('output.txt',skiprows=1)
i = 0
bonds = 0
kix = n_mol[1]
r5_kix = 529
r5_bonds = []
rna_kix = 1247
rigid_bodies = []
start = False
for protein in order:
    protein = int(protein)
    atom_num = int(np.loadtxt(files[protein],skiprows=3,max_rows=1,usecols=[0]))
    rigid_region = np.concatenate((rigid[protein]))
    r = False
    chain = topology.addChain()
    residue = topology.addResidue(name="CG-residue", chain=chain)
    if data[i,3] > 22:
        data[i,3] -= 22
    atom1 = topology.addAtom(name=aminos.amino_acids[data[i,3]], element=None, residue=residue)
    i+=1
    k_use = 0
    for current in range(1, atom_num):
        residue = topology.addResidue(name="CG-residue",chain=chain)
        if data[i,3] > 22:
            data[i,3] -= 22
        atom2 = topology.addAtom(name=aminos.amino_acids[data[i,3]], element=None, residue=residue)
        if current not in rigid_region:
            if r == True:
                r = False
            topology.addBond(atom1, atom2)
            bonds +=1
        else:
            if r == False:
                topology.addBond(atom1, atom2)
                bonds+=1
                r = True
            else:
                pass
        if atom_num == 991:
            if start == True:
                if current == r5_kix:
                    topology.addBond(atom1,kix_atom)
                    start = False
        elif atom_num == 4380:
            if current == rna_kix:
                kix_atom = atom1
                start = True
        atom1 = atom2
        i+=1
# Create the rigid bodies for this protein
    if len(rigid_region)>0:
        rigid_bodies_protein = [(np.array([i-atom_num]) + region).tolist() for region in rigid[protein]]
        rigid_bodies.append(rigid_bodies_protein)

positions = (data[:,5:]/10 * unit.nanometer)
print(topology)
print(len(positions))
system = mm.System()
box_vec_a = np.array([box, 0, 0])*unit.nanometer
box_vec_b = np.array([0, box, 0])*unit.nanometer
box_vec_c = np.array([0, 0, box])*unit.nanometer
system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)

charges = []
atom_types = []
count = 0
for atom in topology.atoms():
    if data[count,4] == 4:
        charges.append(4)
    else:
        charges.append(aminos.charge[atom.name])
    atom_types.append(aminos.reverse[atom.name])
    system.addParticle(aminos.mass[atom.name])
    count += 1


hbond = mm.HarmonicBondForce()


for bond in topology.bonds():
    hbond.addBond(bond.atom1.index, bond.atom2.index, 0.38, 2000*aminos._kcal_to_kj)
hbond.setForceGroup(1)

aa = len(aminos.amino_acids)
sigma_ah_map, lambda_ah_map = np.zeros((aa, aa)), np.zeros((aa, aa))
for i in range(aa):
    for j in range(aa):
        sigma_ah_map[i,j] = np.mean([aminos.size[aminos.amino_acids[i+1]],aminos.size[aminos.amino_acids[j+1]]])/10
        lambda_ah_map[i,j] = np.mean([aminos.hydropathy[aminos.amino_acids[i+1]],aminos.hydropathy[aminos.amino_acids[j+1]]])
lambda_ah_map = aminos.mu*lambda_ah_map - aminos.delta

lj_at_cutoff = 4*aminos.epsilon*((1/4)**12 - (1/4)**6)
contacts = mm.CustomNonbondedForce(f'''energy;
           energy=(f1+f2-offset)*step(4*sigma_ah-r);
           offset=lambda_ah*{lj_at_cutoff};
           f1=(lj+(1-lambda_ah)*{aminos.epsilon})*step(2^(1/6)*sigma_ah-r);
           f2=lambda_ah*lj*step(r-2^(1/6)*sigma_ah);
           lj=4*{aminos.epsilon}*((sigma_ah/r)^12-(sigma_ah/r)^6);
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
if aminos.use_pbc:
    contacts.setNonbondedMethod(contacts.CutoffPeriodic)
else:
    contacts.setNonbondedMethod(contacts.CutoffNonPeriodic)
contacts.setCutoffDistance(4*np.amax(sigma_ah_map))
contacts.setForceGroup(2)

alpha = aminos.NA*aminos.EC**2/(4*np.pi*aminos.VEP)
ldby_value = aminos.ldby.value_in_unit(unit.nanometer)
alpha_value = alpha.value_in_unit(unit.kilojoule_per_mole*unit.nanometer)
cutoff=3.5*unit.nanometer
cutoff_value = cutoff.value_in_unit(unit.nanometer)
elec = mm.CustomNonbondedForce(f'''energy;
       energy=q1*q2*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{aminos.dielectric_water};
       offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
       ''')
elec.addPerParticleParameter('q')
for q in charges:
    if q == 4:
        q = 0
    elec.addParticle([q])
elec.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
if aminos.use_pbc:
    elec.setNonbondedMethod(elec.CutoffPeriodic)
else:
    elec.setNonbondedMethod(elec.CutoffNonPeriodic)
elec.setCutoffDistance(cutoff_value)
elec.setForceGroup(3)

elec_sri = mm.CustomNonbondedForce(f'''energy;
       energy=p1*p2*step(-8-p1*p2)*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{aminos.dielectric_water};
       offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
       ''')
elec_sri.addPerParticleParameter('p')
num = 0
for q in charges:
    if abs(q) >= 2:
        elec_sri.addParticle([q])
        num += 1
    else:
        elec_sri.addParticle([0])
print(num)
elec_sri.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
if aminos.use_pbc:
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

for i in range(len(rigid_bodies)):
    createRigidBodies(system, positions, rigid_bodies[i])

with open('system.xml', 'w') as output:
    output.write(mm.XmlSerializer.serialize(system))

with open('start.pdb','w') as f:
    app.PDBFile.writeFile(topology, positions, f)
