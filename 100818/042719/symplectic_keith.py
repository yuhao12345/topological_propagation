# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 20:21:01 2019

@author: ykang
"""

import kwant
import tinyarray
import math
from matplotlib import pyplot
import numpy as np
import cmath
from kwant.digest import uniform

sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])

def make_system(a, t, W, L, salt):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity).

    def hopping(site1, site2):
        
        alpha=uniform(repr(site1),salt)*2*np.pi   
        beta=np.arccos(1-2*uniform(repr(site1.pos),salt))/2
        gamma=uniform(repr(site2),salt)*2*np.pi
#        np.random.seed(seed=salt)
#        alpha=np.random.rand()*2*np.pi            
#        beta=np.arccos(1-2*np.random.rand())/2
#        gamma=np.random.rand()*2*np.pi   
        return -tinyarray.array([[cmath.exp(1j*alpha)*math.cos(beta), cmath.exp(1j*gamma)*math.sin(beta)],
                          [-cmath.exp(-1j*gamma)*math.sin(beta), cmath.exp(-1j*alpha)*math.cos(beta)]])
#        return np.array([[1,0],[0,1]])
    lat = kwant.lattice.square(a)

    sys = kwant.Builder()

    #### Define the scattering region. ####
    sys[(lat(x, y) for x in range(L) for y in range(W))] = tinyarray.array([[0, 0], [0, 0]])
        
    # hoppings in x-direction
    sys[lat.neighbors()]=hopping

    #### Define the left lead. ####
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))

    lead[(lat(0, j) for j in range(W))] = tinyarray.array([[0, 0], [0, 0]])
    # hoppings in x-direction
    lead[lat.neighbors()]=-sigma_0
    #### Attach the leads and return the finalized system. ####
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys


def plot_conductance(sys, energies):
    # Compute conductance
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(sys, energy)
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()



sys = make_system(1, 1.0, 8, 200, '55')

# Check that the system looks as intended.
#kwant.plot(sys)

# Finalize the system.
sys = sys.finalized()

# We should see non-monotonic conductance steps.
plot_conductance(sys, energies=[0.01 * i  for i in range(50)])
    
#wf=kwant.wave_function(sys,3.7)(0)  #,check_hermiticity=False
#kwant.plotter.map(sys, (abs(wf[0][0:600])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
#g=kwant.smatrix(sys, 0.01).transmission(1, 0)
