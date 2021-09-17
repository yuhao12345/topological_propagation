# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:11:47 2019

@author: user
"""

import kwant
import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform
import scipy.io as sio
import scipy
import os.path
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()

length=30  # toi calculate loc length, sample length=1e4
width=10

dis=4

lat = kwant.lattice.square()

def make_system(length,width,salt):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<=length and abs(y)<width

    def lead_shape(pos):
        x, y= pos
        return abs(y)<width #2
    def onsite(site):
        x,y=site.pos
        if abs(x)==length:
            return 4-1j
        else:
            return dis*(uniform(repr(site),salt)-.5)+4
    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = 1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 4
    lead[lat.neighbors()] = 1

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

sys=make_system(length,width,'5')
kwant.plot(sys)
ha=sys.hamiltonian_submatrix()

#w,v=np.linalg.eig(ha)
#w,v=scipy.sparse.linalg.eigs(ha, k=15, sigma=4.5, which='LR',  
#                         return_eigenvectors=True)
#kwant.plotter.map(sys, (abs(v[:,0])**2),num_lead_cells=2,fig_size=(15, 10))
#
#print([min(np.real(w)),max(np.real(w))])

elapsed=time()-t_ini