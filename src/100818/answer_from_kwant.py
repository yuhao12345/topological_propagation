# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 15:28:32 2018

@author: user
"""

import kwant
from matplotlib import pyplot
import numpy as np

def make_system(a=1, W=30, L=20, barrier=1., barrierpos=(3, 8),
                mu=0.4, Delta=0.1, Deltapos=4, t=1.0):
    # Start with an empty tight-binding system and two square lattices,
    # corresponding to electron and hole degree of freedom
    lat_e = kwant.lattice.square(a, name='e')
    lat_h = kwant.lattice.square(a, name='h')

    sys = kwant.Builder()
    sys1 = kwant.Builder()
    sys2 = kwant.Builder()
    #### Define the scattering region. ####
    sys1[(lat_e(x, y) for x in range(L) for y in range(W))] = 4 * t - mu
    sys2[(lat_h(x, y) for x in range(L) for y in range(W))] = mu - 4 * t
    
    # the tunnel barrier
    sys1[(lat_e(x, y) for x in range(barrierpos[0], barrierpos[1])
         for y in range(W))] = 4 * t + barrier - mu
    sys2[(lat_h(x, y) for x in range(barrierpos[0], barrierpos[1])
         for y in range(W))] = mu - 4 * t - barrier

    # hoppings for both electrons and holes
    sys1[lat_e.neighbors()] = -t
    sys2[lat_h.neighbors()] = t
    sys+=sys1
    sys+=sys2
    # Superconducting order parameter enters as hopping between
    # electrons and holes
    sys[((lat_e(x, y), lat_h(x, y)) for x in range(Deltapos, L)
         for y in range(W))] = Delta

    #### Define the leads. ####
    # Symmetry for the left leads.
    sym_left = kwant.TranslationalSymmetry((-a, 0))

    # left electron lead
    lead0 = kwant.Builder(sym_left)
    lead0[(lat_e(0, j) for j in range(W))] = 4 * t - mu
    lead0[lat_e.neighbors()] = -t

    # left hole lead
    lead1 = kwant.Builder(sym_left)
    lead1[(lat_h(0, j) for j in range(W))] = mu - 4 * t
    lead1[lat_h.neighbors()] = t

    # Then the lead to the right
    # this one is superconducting and thus is comprised of electrons
    # AND holes
    sym_right = kwant.TranslationalSymmetry((a, 0))
    lead2 = kwant.Builder(sym_right)
    lead2 += lead0
    lead2 += lead1
    lead2[((lat_e(0, j), lat_h(0, j)) for j in range(W))] = Delta

    #### Attach the leads and return the system. ####
    sys.attach_lead(lead0)
    sys.attach_lead(lead1)
    sys.attach_lead(lead2)

    return sys,sys1,sys2


    
#I used sys1 and sys2 to help in the 2D-plot
sys,sys1,sys2 = make_system()

# Check that the system looks as intended.
kwant.plot(sys)

# Finalize the system.
sys = sys.finalized()


 #sys and sys1 have the same positions of sites

wf=kwant.wave_function(sys,energy=1)

wavefunction=wf(0)[1].reshape(2,-1)# lead number 0, mode number 2

#plot of wf for electrons:
kwant.plotter.map(sys1.finalized(),abs(wavefunction[0])**2 )


#plot of wf for holes:
kwant.plotter.map(sys1.finalized(),abs(wavefunction[1])**2 )