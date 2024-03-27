import numpy as np
import MDAnalysis as mda
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import openmm as mm
import openmm.unit as unit

start_frame = 100
end_frame = 129
boundary = 10.0
output_density_csv = 'density.csv'
u = mda.Universe( 'slab.dcd')
n_atoms = len(u.atoms)
n_monomers = 500 # monomer number
n_frames = len(u.trajectory)
n_atoms_per_monomer = int(n_atoms/n_monomers)
box_a, box_b, box_c = 35, 35, 350

print(f'Use frames {start_frame}-{end_frame} when computing density')
print(f'The trajecotry has {n_frames} frames')
print(f'The system has {n_monomers} monomers.')

# double check if all the COM are within the main box
# load mass
system_xml = 'system.xml'
with open(system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
atom_mass = []
n_atoms = system.getNumParticles()
n_atoms_per_monomer = int(n_atoms/n_monomers)
for i in range(n_atoms):
    atom_mass.append(system.getParticleMass(i).value_in_unit(unit.dalton))
monomer_mass = []
for i in range(n_atoms_per_monomer):
    monomer_mass.append(system.getParticleMass(i).value_in_unit(unit.dalton))
monomer_mass = np.sum(np.array(monomer_mass))

bin_width = 1
n_bins = int(box_c/bin_width)
bins = np.linspace(0, box_c, n_bins + 1)
bin_width = bins[1] - bins[0] # reset bin_width
z = 0.5*(bins[1:] + bins[:-1])
rho_M, rho_g_per_L = [], []
NA = 6.02214076e+23
for i in range(start_frame, end_frame + 1):
    u.trajectory[i]
    count_i, _ = np.histogram(u.atoms.positions[:,2]/10, bins=bins)
    rho_M_i = count_i/(NA*box_a*box_b*bin_width*10**-27*10**3) # unit mol/L
    rho_g_per_L_i = rho_M_i*monomer_mass # unit g/L
    rho_M.append(rho_M_i)
    rho_g_per_L.append(rho_g_per_L_i)
rho_M = np.mean(np.array(rho_M), axis=0)
rho_g_per_L = np.mean(np.array(rho_g_per_L), axis=0)
z_shifted = z - np.mean(z)
df_density = pd.DataFrame(columns=['z (nm)', 'rho (M)', 'rho (g/L)'])

non_zero_indices = np.nonzero(rho_M)[0]

if non_zero_indices.size > 0:
    # Calculate the amount of shift needed
    start_index = non_zero_indices[0]
    end_index = non_zero_indices[-1]
    shift_amount = (len(rho_M) - (end_index - start_index)) // 2 - start_index

    # Shift the 'rho' array
    rho_M_shifted = np.roll(rho_M, shift_amount)
    # After shifting, ensure that the 'rho' values outside the original non-zero range are set to 0
    rho_M_shifted[:max(0, shift_amount)] = 0
    rho_M_shifted[min(len(rho_M), len(rho_M) + shift_amount):] = 0
else:
    # If all values are zero, just use the original rho_M array
    rho_M_shifted = rho_M

df_density['z (nm)'] = z_shifted
df_density['rho (M)'] = rho_M_shifted
df_density['rho (g/L)'] = rho_g_per_L
df_density.to_csv(output_density_csv, index=False)


# compute the density of two phases
print(f'Dense phase regime is {-1*boundary} <= z <= {boundary}')
rho_dilute_phase = df_density.loc[(df_density['z (nm)'] < -50) | (df_density['z (nm)'] > 50), 'rho (g/L)'].mean()
rho_concentrated_phase = df_density.loc[(df_density['z (nm)'] >= -1*boundary) & (df_density['z (nm)'] <= boundary), 'rho (g/L)'].mean()
print(f'Dilute phase concentration is {rho_dilute_phase:.2f} g/L (mg/mL), and concentrated phase concentration is {rho_concentrated_phase:.2f} g/L (mg/mL)')
z = df_density['z (nm)']
rho = df_density['rho (M)'] # g/L is equivalent to mg/mL
plt.plot(z, rho)
plt.xlabel('z (nm)')
plt.ylabel('Density (M)')
plt.show()
