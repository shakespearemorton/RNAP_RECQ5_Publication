import MDAnalysis as mda
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.spatial import ConvexHull
import csv
colors = plt.get_cmap('Dark2').colors
colo = [1,3,4,5]
def compute_msd(ref_positions, traj_positions):
    """Compute the mean squared displacement."""
    displacements = traj_positions - ref_positions
    squared_displacements = np.sum(displacements * displacements, axis=-1)
    return np.mean(squared_displacements, axis=0)

def average_nearest_neighbor_distance(points,extra=0):
    total_distance = 0
    test = []
    for i in range(len(points)):
        min_distance = float('inf')

        for j in range(len(points)):
            if i != j:
                distance = euclidean_distance(points[i], points[j])
                if distance < min_distance:
                    min_distance = distance
        if extra == 1:
            test.append(min_distance)
        total_distance += min_distance
    if extra == 0:
        return total_distance / len(points)
    else:
        return test

def euclidean_distance(p1, p2):
    """Compute the Euclidean distance between two 3D points."""
    return np.linalg.norm(p1 - p2)

#direc = ['ten','five','three','two']
#direc = ['three','five','ten']
#labels = {direc[0]:'1:3',direc[1]:'1:5',direc[2]:'1:10'}
vols = []
#direc = ['five','seven','nine','eight'] # ,'seventy','eighty','ninety','hundo']
direc = ['sixty']
plt.figure()
alsep = []
labels = ['1:5','1:7','1:9','1:8']
use = 0
#LOAD THE INDEX FOR THE PROTEINS 0 IS RNAPOL, 1 IS RECQL5
for name in direc:
    length_protein1 = 4380
    length_protein2 = 991
    
    protein_assignment = []
    residue_counter = 0
    previous_residue_id = -1
    
    with open('../../../Condensates/Starting_Configurations/10C_5R_2P.pdb', 'r') as pdb_file:
        for line in pdb_file:
            # Look for lines starting with 'HETATM', which indicates atom records
            if line.startswith("HETATM"):
                # Extract residue number from the PDB line
                current_residue_id = int(line[22:26].strip())
                
                # Check for a jump in residue numbering
                if current_residue_id != previous_residue_id + 1 and residue_counter != 0:
                    if residue_counter == length_protein1:
                        protein_assignment.extend([0]*residue_counter)
                    elif residue_counter == length_protein2:
                        protein_assignment.extend([1]*residue_counter)
                    else:
                        raise ValueError(f"Unexpected number of residues: {residue_counter}")
    
                    residue_counter = 0
    
                previous_residue_id = current_residue_id
                residue_counter += 1
    
    # Handle the last set of residues after iterating
    if residue_counter == length_protein1:
        protein_assignment.extend([0]*residue_counter)
    elif residue_counter == length_protein2:
        protein_assignment.extend([1]*residue_counter)
    else:
        raise ValueError(f"Unexpected number of residues: {residue_counter}")
    
    u = mda.Universe('../../../Condensates/Trajectories/10C_5R_2P.dcd')
    time = np.genfromtxt('../../../Condensates/Simulation_Information/10C_5R_2P.csv',delimiter=',')
    time = time[:,0]
    # Use the protein_assignment to get indices of beads belonging to protein 1 copies
    protein1_indices = [i for i, val in enumerate(protein_assignment) if val == 0]
    
    # Compute the start and end indices for each protein 1 copy
    length_protein1 = 4380
    protein1_starts = range(0, len(protein1_indices), length_protein1)
    protein1_ends = range(length_protein1, len(protein1_indices) + 1, length_protein1)
    
    # Compute the center of mass for each protein 1 copy in each frame
    com_coordinates = []
   
    for tsi in range(len(u.trajectory)):
        ts = u.trajectory[tsi]
    for ts in u.trajectory:
        frame_coms = []
        for start, end in zip(protein1_starts, protein1_ends):
            subset_indices = protein1_indices[start:end-462]
            subset_atoms = u.atoms[subset_indices]
            
            # Using numpy's mean function to compute the center of mass
            com = np.mean(subset_atoms.positions, axis=0)
            frame_coms.append(com)
        com_coordinates.append(frame_coms)
        
   
    com_coordinates = np.array(com_coordinates)
    msd_values = []
    for i in range(com_coordinates.shape[1]):
        ref_position = com_coordinates[0, i]  # Take the initial position as reference
        msd = compute_msd(ref_position, com_coordinates[:, i])
        msd_values.append(msd)
    
    msd_values = np.array(msd_values)
    avg_msd = np.mean(msd_values, axis=0)
    avg_separation = []
    
    for frame_coms in com_coordinates:
        pairwise_distances = average_nearest_neighbor_distance(frame_coms,extra=1)
        avg_separation.append(np.mean(pairwise_distances))
    alsep.append(pairwise_distances)

    y = np.asarray(avg_separation)/10
    x = time*1e-8
    plt.plot(y, label=name,linewidth=4,c=colors[colo[use]])
    print(name, avg_separation[0]/10)
    print(name, np.mean(y[-5:]), np.std(y[-5:]))
    use += 1
plt.axhline(12.51, c='r',linestyle = '--')
plt.rcParams.update({'font.size': 14})
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=3, fancybox=True, shadow=True)
#plt.xlim([0,10])
#plt.ylim([11,14])
plt.xlabel("Time ($\mu$s)")
plt.ylabel("Average COM Separation (nm)")
plt.text(5, 12.1,'Experimental\n(12.5)', ha='center', color='r')
plt.savefig('com.pdf',format='pdf')
np.savetxt('com.txt',y)

plt.show()
