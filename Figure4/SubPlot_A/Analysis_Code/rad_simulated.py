import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import analyticalrdf
import glob
import pandas as pd
from lxml import etree


def load(fname):
    tree = etree.parse(fname)
    root = tree.getroot()
    data = []

    for marker in root.iter('marker'):
        # Extract the relevant attributes
        id_ = int(marker.attrib['id'])
        x = float(marker.attrib['x'])
        y = float(marker.attrib['y'])
        z = float(marker.attrib['z'])
        coordX = float(marker.attrib['coordX'].rstrip('px'))
        coordY = float(marker.attrib['coordY'].rstrip('px'))
        coordZ = float(marker.attrib['coordZ'].rstrip('px'))
        data.append([id_, x, y, z, coordX, coordY, coordZ])

    df = pd.DataFrame(data, columns=['id', 'x', 'y', 'z', 'coordX', 'coordY', 'coordZ'])
    df = center(df)
    return df[['x', 'y', 'z']].values

def center(df):
    mean_x = df['x'].mean()
    mean_y = df['y'].mean()
    mean_z = df['z'].mean()
    df_centered = df.copy()
    df_centered['x'] = df['x'] - mean_x
    df_centered['y'] = df['y'] - mean_y
    df_centered['z'] = df['z'] - mean_z

    return df_centered

plt.figure(dpi=1200)

'''files = glob.glob('*.cmm')
exp = []
for file in files:
    points = load(file)/10
    r = np.arange(0,100,0.5)
    gr = analyticalrdf.RDF_AnalyticNorm(points,r,0.5)
    for i in range(len(r)):
        if gr[-1-i] > 0:
            cut = r[-1-i]/2
            break
    V = (4/3) * np.pi * cut**3
    rho = len(points) / V
    rsp = (4*np.pi*r**2)
    Fsphere = rsp/(rsp*(1-(r/(2*cut))**2)*(1+r/(4*cut)))
    for i in range(len(r)):
        if r[i] > cut:
            imax = i
            break
    sgr = gr*Fsphere
    exp.append(sgr)
stds = np.std(exp,axis=0)
rdfs = np.mean(exp,axis=0)
plt.plot(r[:imax],rdfs[:imax], linewidth=3,label='RNAPII - exp',c='rosybrown')
plt.fill_between(r[:imax],rdfs[:imax]+stds[:imax],rdfs[:imax]-stds[:imax],alpha=0.3,color='rosybrown')
'''
    
rnapol = 4380
r5 = 991
file = ['60C_5R_2P']
names = ['RNAPII','RECQ5']
color = ["dimgrey", "darkgoldenrod", "#2ca02c", "#d62728",'#9467bd'] 
use = 0


for name in file:
    protein_assignment = []
    residue_counter = 0
    previous_residue_id = -1
    resnames = []
    with open(name+'.pdb', 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("HETATM"):
                current_residue_id = int(line[22:26].strip())
                if current_residue_id != previous_residue_id + 1 and residue_counter != 0:
                    if residue_counter == rnapol:
                        protein_assignment.extend([0]*residue_counter)
                    elif residue_counter == r5:
                        protein_assignment.extend([1]*residue_counter)
                    else:
                        raise ValueError(f"Unexpected number of residues: {residue_counter}")
    
                    residue_counter = 0
                
                resnames.append(line[12:15].strip())
                previous_residue_id = current_residue_id
                residue_counter += 1
    
    
    if residue_counter == rnapol:
        protein_assignment.extend([0]*residue_counter)
    elif residue_counter == r5:
        protein_assignment.extend([1]*residue_counter)
    else:
        raise ValueError(f"Unexpected number of residues: {residue_counter}")
        
    rnap_indices = [i for i, val in enumerate(protein_assignment) if val == 0]    
    rnap_copies = [resnames[i] for i in rnap_indices[0:rnapol]]
    rnap_copies = [rnap_indices[i:i+3900] for i in range(0, len(rnap_indices), rnapol)]
    ctd_copies = [rnap_indices[i+rnapol-364:i+rnapol] for i in range(0, len(rnap_indices), rnapol)]
    r5_indices = [i for i, val in enumerate(protein_assignment) if val == 1]
    r5_res = [resnames[i] for i in r5_indices[0:r5]]
    r5_copies = [r5_indices[i:i+400] for i in range(0, len(r5_indices), r5)]
    u = mda.Universe(name+'.dcd')
    n_frames = len(u.trajectory)
    use_frames = int(n_frames*.05)

copies = [rnap_copies,r5_copies]
for copier in copies:    
    rdfs = []
    for frame in np.linspace(.33*n_frames,n_frames-1,5).astype(int):
        u.trajectory[frame]
        com = []
        for copy in copier:
            subset_atoms = u.atoms[copy]
            com.append(np.mean(subset_atoms.positions, axis=0))
        com = np.stack(com)/10
        r = np.arange(0,100,0.5)
        gr = analyticalrdf.RDF_AnalyticNorm(com,r,0.5)
        for i in range(len(r)):
            if gr[-1-i] > 0:
                cut = r[-1-i]/2
                break
        V = (4/3) * np.pi * cut**3
        rho = len(com) / V
        rsp = (4*np.pi*r**2)
        Fsphere = rsp/(rsp*(1-(r/(2*cut))**2)*(1+r/(4*cut)))
        for i in range(len(r)):
            if r[i] > cut:
                imax = i
                break
        sgr = gr*Fsphere
        rdfs.append(sgr)
    stds = np.std(rdfs,axis=0)
    rdfs = np.mean(rdfs,axis=0)
    r *= 10
    plt.plot(r[:imax],rdfs[:imax], linewidth=3,c=color[use],label=names[use])
    plt.fill_between(r[:imax],rdfs[:imax]+stds[:imax],rdfs[:imax]-stds[:imax],color=color[use],alpha=0.3)
    use += 1
plt.legend()
plt.xlim([1,300])
plt.ylim([0,4])
plt.xlabel('Distance ($\AA$)',fontsize=15)
plt.ylabel('RDF',fontsize = 15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.axvline(x=125.45, color='dimgrey', linestyle="--", linewidth=2)
plt.text(125.45+50, plt.ylim()[1]*0.80,'Experimental\n(125.45)', ha='center', color='dimgrey',fontsize=12)
plt.savefig('rdf.pdf',format='pdf')
