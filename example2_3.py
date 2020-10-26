# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 15:08:10 2017

@author: user
"""

import kwant
import tinyarray
# For plotting
from matplotlib import pyplot

sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])
e_z=1
alpha=1.8
def make_system(a=1, t=1.0, W=10, L=30):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity.
    lat = kwant.lattice.square(a)

    syst = kwant.Builder()
    
    syst[(lat(x, y) for x in range(L) for y in range(W))] = \
        4 * t * sigma_0 + e_z * sigma_z
    # hoppings in x-direction
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_y
    # hoppings in y-directions
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_x
        
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
    lead[(lat(0, j) for j in range(W))] = 4 * t * sigma_0 + e_z * sigma_z
    # hoppings in x-direction
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_y
    # hoppings in y-directions
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_x
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst

def plot_conductance(syst, energies):
    # Compute conductance
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()
    
def main():
    syst = make_system()

    # Check that the system looks as intended.
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # We should see conductance steps.
    plot_conductance(syst, energies=[0.01 * i for i in range(100)])
    
if __name__ == '__main__':
    main()