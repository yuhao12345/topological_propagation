# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 11:09:28 2019

@author: ykang
"""
import kwant
import tinyarray
from matplotlib import pyplot

sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])

def make_system(a=1, t=1.0, alpha=0.5, W=10, L=30):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity).
    lat = kwant.lattice.square(a)

    sys = kwant.Builder()

    #### Define the scattering region. ####
    sys[(lat(x, y) for x in range(L) for y in range(W))] = \
        4 * t * sigma_0 
    # hoppings in x-direction
    sys[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_y
    # hoppings in y-directions
    sys[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_x

    #### Define the left lead. ####
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))

    lead[(lat(0, j) for j in range(W))] = 4 * t * sigma_0 
    # hoppings in x-direction
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_y
    # hoppings in y-directions
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_x

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



sys = make_system()

# Check that the system looks as intended.
kwant.plot(sys)

# Finalize the system.
sys = sys.finalized()

# We should see non-monotonic conductance steps.
#plot_conductance(sys, energies=[0.01 * i - 0.3 for i in range(100)])
    
wf=kwant.wave_function(sys,0)(0)  #,check_hermiticity=False
#kwant.plotter.map(sys, (abs(wf[0][0:600])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
