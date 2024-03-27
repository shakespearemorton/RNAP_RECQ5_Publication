import MDAnalysis as mda
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

mass = {
    'MET': 131.199997,
    'GLY': 57.049999,
    'LYS': 128.199997,
    'THR': 101.099998,
    'ARG': 156.199997,
    'ALA': 71.080002,
    'ASP': 115.099998,
    'GLU': 129.100006,
    'TYR': 163.199997,
    'VAL': 99.07,
    'LEU': 113.199997,
    'GLN': 128.100006,
    'TRP': 186.199997,
    'PHE': 147.199997,
    'SER': 87.080002,
    'HIS': 137.100006,
    'ASN': 114.099998,
    'PRO': 97.120003,
    'CYS': 103.099998,
    'ILE': 113.199997,
    'SEP': 165.03,
    'SRI': 1
}

def compute_radius_of_gyration(coords, masses):
    """Compute the radius of gyration for a set of coordinates with given masses."""
    com = np.average(coords, weights=masses, axis=0)
    sq_distances = np.sum((coords - com)**2, axis=1)
    rg = np.sqrt(np.sum(masses * sq_distances) / np.sum(masses))
    return rg

def getRh(rg,N):
    a1 = 0.216
    a2 = 4.06
    a3 = 0.821
    return rg/((a1*(rg-a2*N**(0.33)))/(N**(0.60)-N**(0.33))+a3)

# Load the Universe from the topology and trajectory files
names = ['ctd','hpctd','idr','r5','p5ctd']
for i in names:
    plt.figure()
    u = mda.Universe('rg_'+i+'.dcd')
    
    # Extract amino acid names directly from PDB file
    amino_acids = []
    with open('rg_'+i+'.pdb', 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("HETATM"):
                amino_acids.append(line.split()[2])
    
    # Determine the masses of each amino acid from the aminos.py dictionary
    masses = [mass[aa] for aa in amino_acids]
    masses = np.array(masses)
    
    
    # Compute the radius of gyration for each frame
    rg = []
    n_frames = len(u.trajectory)
    for ts in range(int(n_frames*.7)):
        u.trajectory[-1-ts]
        coords = u.atoms.positions
        rg.append(compute_radius_of_gyration(coords, masses))
    rg = np.asarray(rg)
    rg /= 10
    sb.kdeplot(rg,label='Gyration',c='k',linestyle='--')
    
    
    
    Rh = getRh(np.asarray(rg), len(masses))
    sb.kdeplot(Rh,label='Hydrodynamic',c='k',linestyle='-')
    plt.xlabel('Radius (nm)')
    plt.legend()
    plt.title(i)
    print(np.mean(rg),np.mean(Rh))
