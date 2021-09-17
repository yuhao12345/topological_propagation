# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 16:35:58 2018

@author: user
"""
import kwant

import numpy as np
import tinyarray

# For plotting
from matplotlib import pyplot

tau_x = tinyarray.array([[0, 1], [1, 0]])
tau_z = tinyarray.array([[1, 0], [0, -1]])

def make_system(a=1, W=10, L=10, barrier=1.5, barrierpos=(3, 4),
                mu=0.4, Delta=0.1, Deltapos=4, t=1.0):
    # Start with an empty tight-binding system and two square lattices,
    # corresponding to electron and hole degree of freedom
    lat_e = kwant.lattice.square(a, name='e')
    lat_h = kwant.lattice.square(a, name='h')
    sys = kwant.Builder()

    #### Define the scattering region. ####
    sys[(lat_e(x, y) for x in range(L) for y in range(W))] = 4 * t - mu
    sys[(lat_h(x, y) for x in range(L) for y in range(W))] = mu - 4 * t

    # the tunnel barrier
    sys[(lat_e(x, y) for x in range(barrierpos[0], barrierpos[1])
         for y in range(W))] = 4 * t + barrier - mu
    sys[(lat_h(x, y) for x in range(barrierpos[0], barrierpos[1])
         for y in range(W))] = mu - 4 * t - barrier

    # hoppings for both electrons and holes
    sys[lat_e.neighbors()] = -t
    sys[lat_h.neighbors()] = t

    # Superconducting order parameter enters as hopping between
    # electrons and holes
    sys[((lat_e(x, y), lat_h(x, y)) for x in range(Deltapos, L)
         for y in range(W))] = Delta

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
    sym_right = kwant.TranslationalSymmetry((a, 0))
    lead2 = kwant.Builder(sym_right)
    lead2 += lead0
    lead2 += lead1
    lead2[((lat_e(0, j), lat_h(0, j)) for j in range(W))] = Delta
    sys.attach_lead(lead0)
    sys.attach_lead(lead1)
    sys.attach_lead(lead2)

    return sys

def plot_conductance(sys, energies):
    # Compute conductance
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(sys, energy)
        # Conductance is N - R_ee + R_he
        data.append(smatrix.submatrix(0, 0).shape[0] -
                    smatrix.transmission(0, 0) +
                    smatrix.transmission(1, 0))
    pyplot.plot(data)
sys=make_system().finalized()       
plot_conductance(sys, np.linspace(0,0.2,10))