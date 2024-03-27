import openmm as mm
import openmm.app as app
import openmm.unit as unit
import math
import numpy as np
from sys import stdout
import aminos


topology = app.Topology()
positions = []

#### Make SRI ####
chain = topology.addChain()
residue = topology.addResidue(name="CG-residue", chain=chain)
atom1 = topology.addAtom(name='SEP', element=None, residue=residue)
positions.append([0,0,7])

#### Make CTD ####
chain = topology.addChain()
residue = topology.addResidue(name="CG-residue", chain=chain)
atom1 = topology.addAtom(name='SEP', element=None, residue=residue)
positions.append([0,0,0])
positions *= unit.nanometer


system = mm.System()
b = 1000
box_vec_a = np.array([b/10, 0, 0])*unit.nanometer
box_vec_b = np.array([0, b/10, 0])*unit.nanometer
box_vec_c = np.array([0, 0, b/10])*unit.nanometer
system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)

charges = [1.1,-2]
atom_types = []
for atom in topology.atoms():
    atom_types.append(aminos.reverse[atom.name])
    system.addParticle(aminos.mass[atom.name])

r0 = 4 * unit.nanometer
k = 20 * aminos._kcal_to_kj
restraint = mm.CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
restraint.addGlobalParameter('r0', r0)
restraint.addGlobalParameter('k',k)
g1 = [0]
g2 = [1]
restraint.addGroup(g1) # index will be 0
restraint.addGroup(g2) # index will be 1
restraint.addBond([0,1], []) #([g1_index,g2_index], [per-bond parameter list])
restraint.setForceGroup(6)
system.addForce(restraint)

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
    if q == 1.1:
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
       energy=p1*p2*step(0-p1*p2)*{alpha_value}*((exp(-r/{ldby_value})/r)-offset)*step({cutoff_value}-r)/{aminos.dielectric_water};
       offset={math.exp(-cutoff_value/ldby_value)/cutoff_value};
       ''')
elec_sri.addPerParticleParameter('p')
for q in charges:
    if abs(q) >= 2:
        elec_sri.addParticle([q])
        print(q)
    elif q == 1.1:
        elec_sri.addParticle([-2])
        print('present')
    else:
        elec_sri.addParticle([0])
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
#system.addForce(elec_sri)
system.addForce(force)

#createRigidBodies(system, positions, rigid_body)

with open('system.xml', 'w') as output:
    output.write(mm.XmlSerializer.serialize(system))

with open('start.pdb','w') as f:
    app.PDBFile.writeFile(topology, positions, f)






