from sim_urry import *
import os
import numpy as np

experiment_name = 'ratio_5'
proteins = ['../../Protein_Files/combined.txt','../../Protein_Files/r5.txt']                                         # File(s) with starting configurations of the carbon alphas
num = [20,80]                                                                      # How many of each protein you want to add in the system. They'll be randomly orientated
styles = {0:'compressed',1:'uncompressed',2:'slab'}
style = 0
runtime = 0                                                                   # Femtoseconds (timestep = 10 fs)
sampling = 1e6
temperature = 300                                                               # Kelvin
boxsize = 1000                                                                   # Angstroms
pbc = True                                                                      # Use periodic boundary conditions?

exp = Disordered_Life(experiment_name, proteins, num, style, runtime,sampling, temperature, boxsize, pbc,special_rule=True)
exp.execute()