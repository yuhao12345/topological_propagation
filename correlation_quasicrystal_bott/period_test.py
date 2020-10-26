# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 13:36:16 2018

@author: ykang
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants
 
import kwant
import sys
sys.path.append('./modules/')
import wraparound
 
# -------------------------------------------------------------------
# Constants
qe = sp.constants.value("elementary charge") # unit: C
me = sp.constants.value("electron mass")/qe*1e-18 #unit: eV*s^2/nm^2
hP = sp.constants.value("Planck constant in eV s") #unit: eV*s
hbar = hP/(2*sp.pi) #unit: eV*s
# -------------------------------------------------------------------
 
a = 0.03
V0 = 0.0
 
W = 20
L = 50
 
t = hbar**2/(2*me*a**2) # units: eV
 
lat = kwant.lattice.square(a)
 
# Infinite potential plane in y direction
sys = kwant.Builder(kwant.TranslationalSymmetry(lat.vec((0, W))))
sys[(lat(i,j) for i in range(L) for j in range(W))] = lambda p: 4*t
sys[lat.neighbors(1)] = -t
 
sys = wraparound.wraparound(sys)
 
lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0), lat.vec((0, W))))
lead[(lat(0,j) for j in range(W))] = 4*t
lead[lat.neighbors(1)] = -t
 
lead = wraparound.wraparound(lead, keep=0)
 
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())
 
kwant.plot(sys)
 
sys = sys.finalized()
 
# -------------------------------------------------------
# Calculation
ky = 0.0
 
energies = np.arange(0.0, 5.0, 0.05)
transmission = []
num_prop = []
for energy in energies:
    smatrix = kwant.smatrix(sys, energy, [ky])
    transmission.append(smatrix.transmission(1, 0))
    num_prop.append(smatrix.num_propagating(0))
# -------------------------------------------------------
 
# Plot transmission and propagating modes
plt.plot(energies, transmission, '.')
plt.show()
 
plt.plot(energies, num_prop, '.')
plt.show()
 
# Plot wave function squared for the first mode for a specified energy and ky
wf = kwant.solvers.default.wave_function(sys, energy=2.0, args=[0.0])
kwant.plotter.map(sys, np.abs(wf(0)[0])**2, fig_size=(8,5));