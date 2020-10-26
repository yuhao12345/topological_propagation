# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 16:25:38 2018

@author: ykang
"""

# according to the replying of kwant, create a system with lead and periodic bc 
# can plot a good figure, but cannot calculate scattering matrix and wavefunction
import kwant
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin, cos

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

def onsite(site, mu):
    return -mu

def hopping_1(site1, site2, t1):
    return t1

def hopping_2(site1, site2, t2):
    return t2

# create a supercell with periodic bc in both direction, define the width of sample
mod_graphene = kwant.Builder(kwant.TranslationalSymmetry(*graphene.prim_vecs))
mod_graphene[graphene.shape(lambda pos: True, (0, 0))] = onsite
mod_graphene[graphene.neighbors(n=1)] = hopping_1
mod_graphene[graphene.neighbors(n=2)] = hopping_2

kwant.plot(mod_graphene);


# define the length of sample
#def shape(site): return abs(site.pos[0]) < 5
#
##the horizontal direction is not periodic
#lead = kwant.wraparound.wraparound(mod_graphene, keep=0)  
#
#syst = kwant.Builder()
#syst.fill(lead, shape, (0, 0))
#syst.attach_lead(lead)
#syst.attach_lead(lead.reversed())
#kwant.plot(syst, fig_size=(10, 5));
#syst=syst.finalized()
#gf_mode=kwant.smatrix(syst,1,[0])
#T=gf_mode.transmission(1,0)
