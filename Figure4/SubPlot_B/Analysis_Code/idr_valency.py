import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from progressbar import progressbar
import pandas as pd

def save_data_to_csv(x_data, y_data, std_dev=None, filename='output.csv'):
    # Create a dictionary from the input data
    data_dict = {'X Data': x_data, 'Y Data': y_data}
    
    # If standard deviation is provided, add it to the dictionary
    if std_dev is not None:
        data_dict['Standard Deviation'] = std_dev
    
    # Create a pandas DataFrame from the dictionary
    df = pd.DataFrame(data_dict)
    
    # Save the DataFrame to a CSV file
    df.to_csv(filename, index=False)
    
    print(f'Data saved to {filename}')


rnapol = 4380
ctd = 364
r5 = 991
iri = 200
sri = [920,921,938,939,943,947]
protein_assignment = []
residue_counter = 0
previous_residue_id = -1
resnames = []
proteins = []
protein = []

index = 0
pstop = False
with open('../../../Condensates/Starting_Configurations/70R.pdb', 'r') as pdb_file:
    for line in pdb_file:
        if line.startswith("HETATM"):
            current_residue_id = int(line[22:26].strip())
            if current_residue_id != previous_residue_id + 1 and residue_counter != 0:
                pstop = True
                if residue_counter == rnapol:
                    protein_assignment.extend([0]*residue_counter)
                elif residue_counter == r5:
                    protein_assignment.extend([1]*residue_counter)
                else:
                    raise ValueError(f"Unexpected number of residues: {residue_counter}")

                residue_counter = 0
            
            resnames.append(line[12:15].strip())
            
            previous_residue_id = current_residue_id
            if pstop:
                if len(proteins) > 0:
                    protein.append(index)
                proteins.append(protein)
                protein = []
                pstop=False
            else:
                protein.append(index)
            residue_counter += 1
            index +=1
            
    protein.append(index)        
    proteins.append(protein) 
    if residue_counter == rnapol:
        protein_assignment.extend([0]*residue_counter)
    elif residue_counter == r5:
        protein_assignment.extend([1]*residue_counter)
    else:
        raise ValueError(f"Unexpected number of residues: {residue_counter}")
    
    
u = mda.Universe('../../../Condensates/Trajectories/70R.dcd')
protein1_all_indices = [i for i, val in enumerate(protein_assignment) if val == 0]
r5_indices = [i for i, val in enumerate(protein_assignment) if val == 1]
ctd_indices = []
for start_idx in range(0, len(protein1_all_indices), rnapol):
    ctd_indices.extend(protein1_all_indices[start_idx + rnapol - ctd:start_idx + rnapol])
iri_indices = []   
for start_idx in range(0, len(r5_indices), r5):
    iri_indices.extend(r5_indices[start_idx + 625:start_idx + 825])

ctd_res = [resnames[i] for i in ctd_indices[0:ctd]]
iri_res = [resnames[i] for i in iri_indices[0:iri]]
r5_res = [resnames[i] for i in r5_indices[0:r5]]

indexes = {'r5':r5_indices,'ctd':ctd_indices,'iri':iri_indices}
lengths = {'r5':r5,'ctd':ctd,'iri':iri}
residues = {'r5':r5_res,'ctd':ctd_res,'iri':iri_res}

r5_copies = [r5_indices[i:i+r5] for i in range(0, len(r5_indices), r5)]
iri_copies = [iri_indices[i:i+iri] for i in range(0, len(iri_indices), iri)]
pol_copies = [protein1_all_indices[i+1247] for i in range(0, len(protein1_all_indices), rnapol)]
ctd_copies = [ctd_indices[i:i+ctd] for i in range(0,len(ctd_indices),ctd)]

protein_numbers = {}
for i, prot in enumerate(proteins):
    for idx in prot:
        protein_numbers[idx] = i

contacts = {}
n_frames = len(u.trajectory)
user = np.arange(0,n_frames-1,1).astype(int)

dont = [7,9,15,26,34,46,48,68,69] #Specifically for this R5_DROPLET, do not use proteins not in the condensate
valency_step = []
for frame in progressbar(user):
    u.trajectory[frame]
    valency_copy = []
    num = 0
    for r in iri_copies:
        if num in dont:
            pass
        else:
            # Determine which protein number of R5 we are working with
            pnum_r = protein_numbers[r[0]]
    
            distances = distance_array(u.atoms[r].positions, 
                                       u.atoms[iri_indices].positions, 
                                       box=u.dimensions)
            mask = (distances > 30)*1e8
            min_distances = np.min(distances*mask, axis=0)
    
            # Identify all iri_indices that are in contact with this r domain
            in_contact_indices = np.where(min_distances < 1e7)[0]
            contact_copies = [iri_indices[i] for i in in_contact_indices]
            i_temp = []
            for c in contact_copies:
                # Determine which protein number is in contact with this copy of R5
                i_temp.append(protein_numbers[c])
            
            valency_copy.append(len(np.unique(i_temp))-1)
        num += 1
    valency_step.append(valency_copy)
    

v = np.asarray(valency_step)
v_m = np.mean(v,axis=1)
v_s = np.std(v,axis=1)
user = user/100
plt.figure(dpi=1200)
plt.fill_between(user,8.8+3.3,8.8-3.3,alpha = 0.5,color='lightsteelblue')
plt.plot(user,v_m+20,c='lightsteelblue', label='Homotypic',linewidth = 4)
#plt.plot(user,v_m,c='darkgoldenrod', label='Heterotypic',linewidth = 4)
#plt.scatter(user,v_m,c='lightsteelblue', label='Homotypic',s=0.01)
plt.scatter(user,v_m,c='darkgoldenrod', label='Heterotypic',s=20)
#plt.fill_between(user,v_m+v_s, v_m-v_s,alpha=0.2,color='darkgoldenrod')
plt.legend()
plt.ylim([3,13])
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Valency of $RECQ5^{625-825}$')
plt.savefig('idr_valency.png')
print('IDR-IDR',np.mean(v_m),np.mean(v_s))
save_data_to_csv(user,v_m,v_s,'homogeneous_raw_data.csv')


