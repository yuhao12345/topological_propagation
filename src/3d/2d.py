# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 15:30:53 2018

@author: user
"""

import kwant
#import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform

t_ini=time()

lat = kwant.lattice.square()
salt='5'
dis=1
hop=0.5
def make_cuboid(a):
    def cuboid_shape(pos):
        x, y= pos
        return 0 <= x < a and 0 <= y < a
    
    def onsite(site):
        return dis*(uniform(repr(site),salt)-.5)

    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = hop
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    def lead_shape(pos):
        x, y= pos
        return 0 <= y < a 

    lead[lat.shape(lead_shape, (0, 0))] = 0
    lead[lat.neighbors()] = hop

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

en=1


sys = make_cuboid(16)
gf=kwant.greens_function(sys,en)
gfm=gf.submatrix(1,0)
gg=np.zeros(shape=(16,16),dtype=complex)
for m in range(16):
    for n in range(16):
        gg[n,m]=np.sin((1+np.arange(16))/17*(1+n)*np.pi) @ gfm @ np.transpose(np.sin((1+np.arange(16))/17*(1+m)*np.pi))/2
gf_mode=kwant.smatrix(sys, en)
t=gf.transmission(1,0)
gfmm=gf_mode.submatrix(1,0)
t_mode=gf_mode.transmission(1,0)

sum(np.linalg.svd(gg,compute_uv=False)**2)
sum(np.linalg.svd(gfmm,compute_uv=False)**2)
ind=sys.lead_interfaces
sys.pos(ind[1][10])
