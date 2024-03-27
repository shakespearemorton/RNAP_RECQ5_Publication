import numpy as np
import pandas as pd
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.lib.distances import distance_array
from scipy.spatial.transform import Rotation as R

def add(files,n_mol,box_a=500,box_b=500,box_c=500,output='output.txt',atoms=None):
    rna = n_mol[-1]
    for i in range(len(files)):
        method = 'FastNS'
        max_n_attempts=100000
        try:
            if atoms==None:
                atoms = pd.DataFrame()
        except:
            pass
        atom_num = int(np.loadtxt(files[i],skiprows=3,max_rows=1,usecols=[0]))
        data = np.loadtxt(files[i],skiprows=68,max_rows=atom_num) #residue is on index 2
        new_atoms = pd.DataFrame(data)
        new_atoms.columns =['index','molecule','residue','charge','x','y','z']
        if atom_num == 991:
            n_mol[i] -= rna
        if atom_num == 4380:
            data = np.loadtxt('r5.txt',skiprows=68,max_rows=991)
            r5 = pd.DataFrame(data)
            r5.columns =['index','molecule','residue','charge','x','y','z']
            r5['molecule'] = 2
            new_atoms = pd.concat([new_atoms,r5])
        new_coord = new_atoms[['x', 'y', 'z']].to_numpy()
        new_coord -= np.mean(new_coord, axis=0)
        cutoff = 10.0
        count_n_mol = 0
        count_n_attempts = 0
        dim = np.array([box_a, box_b, box_c, 90.0, 90.0, 90.0])
        while (count_n_mol < n_mol[i]) and (count_n_attempts < max_n_attempts):
            # get a random rotation
            rotate = R.random()
            new_coord_i = rotate.apply(new_coord)
            # get a random translation
            translate = np.random.uniform(0, 1, 3)*np.array([box_a, box_b, box_c])
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
                    new_atoms_i[['x', 'y', 'z']] = new_coord_i
                    atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
                    count_n_mol += 1
                    if method == 'FastNS':
                        coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
                        grid_search = FastNS(cutoff, coord, dim, pbc=True)
            count_n_attempts += 1
        if count_n_mol == n_mol[i]:
            print(f'Successfully inserted {n_mol[i]} molecules.')
        else:
            print(f'Could not successfully insert {n_mol[i]} molecules in {count_n_attempts} attempts.')
            print(f'Only added {count_n_mol} molecules. Try increasing the box size or number of attempts to add more molecules.')
    atoms.to_csv(output, sep ='\t')
