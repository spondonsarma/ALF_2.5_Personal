#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 05:20:44 2020

@author: fassaad
"""

from py_alf import Simulation            # Interface with ALF
import numpy as np                       # Numerical library
sims = []                                # Vector of Simulation instances
for Ham_Uf in [0.0, 0.1, 0.2]:           # Values of dtau
    print(Ham_Uf)
    sim_dict = {"Model": "Kondo", 
                "Lattice_type": "Bilayer_square", 
                "L1": 4 , "L2": 4, 
                "Ham_Uf": Ham_Uf,  
                "Beta": 5.0, 
                "Nsweep": 100, 
                "NBin": 5}
    sim = Simulation('Kondo', sim_dict,
                     alf_dir = '/Users/fassaad/Programs/ALF/Work/',
                     branch = '151-introduce-sun-kondo-hamiltonian-for-bilayer-lattices')
    sims.append(sim)
#sims[0].compile(target = "Kondo")
Con = np.empty((len(sims), 2))          # Matrix for storing energy values
Uf  = np.empty((len(sims),))           # Matrix for Dtau values, for plotting
for i, sim in enumerate(sims):
    #sim.run()
    #sim.analysis() 
    Uf[i] = sim.sim_dict['Ham_Uf']                             # Store Uf value
    Con[i] = sim.get_obs(['Constraint_scalJ'])['Constraint_scalJ']['obs']  # Store constraint

print(Uf)
print(Con)
print()
with open('Constraint.dat', 'w') as file:
    for i in range(len(sims) ):
        file.write('%6.6f\t' % (Uf[i])) 
        file.write('%6.6f\t' % (Con[i,0])) 
        file.write('%6.6f\t' % (Con[i,1])) 
        file.write('\n') 
