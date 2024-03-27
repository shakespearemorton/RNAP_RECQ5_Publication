import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import matplotlib.pyplot as plt
import openmm.app as app
from progressbar import progressbar

rnapol = 4380
ctd_len = 364
r5 = 991
sri_residues = [920,921,938,939,943,947]

file = 'start.pdb'
top = app.PDBFile(file).getTopology()
ctd = []
idr = []
sri = []

for chain in top.chains():
    c = []
    for atom in chain.atoms():
        c.append(atom.index)
    if len(c) == rnapol:
        ctd.append(c[rnapol-ctd_len:])
    if len(c) == r5:
        idr.append(c[625:825])
        s = [c[resi] for resi in sri_residues]
        sri.append(s)

u = mda.Universe('slab.dcd')
n_frames = len(u.trajectory)
u = mda.Universe('slab.dcd')

plt.figure(dpi = 1200)
user = np.arange(len(u.trajectory))
valency_step = []
for frame in progressbar(user):
    u.trajectory[frame]
    valency_copy = []
    for c_cop in ctd:
        i_temp = 0
        for sri_copy in sri:
            distances = distance_array(u.atoms[sri_copy].positions, 
                                       u.atoms[c_cop].positions, 
                                       box=u.dimensions)
            if np.min(distances) < 30:
                i_temp += 1
        valency_copy.append(i_temp)
    valency_step.append(valency_copy)
v = np.asarray(valency_step)
v_m = np.mean(v,axis=1)
v_s = np.std(v,axis=1)
plt.scatter(user/100,v_m,c='b',label='SRI-CTD')
plt.fill_between(user/100,v_m+v_s, v_m-v_s,alpha=0.3,color='b')
print('SRI-CTD',np.mean(v_m),np.mean(v_s))
np.savetxt('sixty_sri_valency.txt',v_m)
np.savetxt('sixty_sri_valency_std.txt',v_s)

plt.ylim([1,20])
plt.xlabel('Time ($\mu$s)')
plt.ylabel('Valency of CTD')
plt.savefig('sri_interact.pdf',format='pdf')
plt.show()
