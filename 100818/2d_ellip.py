# -*- coding: utf-8 -*-
"""
Created on Mon May 20 13:48:45 2019

@author: ykang
"""

import kwant
import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform
import scipy.io as sio
import os.path
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()

length=50  # toi calculate loc length, sample length=1e4
width=80

dis=3

lat = kwant.lattice.square()

def make_system(length,width,salt):
    def cuboid_shape(pos):
        x, y= pos
        return x**2/length**2+y**2/width**2<1
    #    return abs(x) + abs(y)*np.sqrt(3) < 100
    #    return abs(x)+10>abs(y)*3 and abs(x)<100
    def lead_shape(pos):
        x, y= pos
        return abs(y)<40 #2
    def onsite(site):
        x,y=site.pos
        if abs(y)<20:
            return 4
        else:
            return dis*(uniform(repr(site),salt)-.5)+4
    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = -1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 4
    lead[lat.neighbors()] = -1

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

#sys=make_system(length,width,'5')
#kwant.plot(sys)


def tau_n(cishu):
    en=0.4
    salt=cishu+random.random()*1000
    sys=make_system(length,width,str(salt))
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#    s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False) #u, s, vh
    
    wf=kwant.wave_function(sys,en)(0)
#    kwant.plotter.map(sys, (abs(wf[14])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

    ldos = kwant.ldos(sys,en)
    if cishu==0:
        coord=np.array([sys.pos(i) for i in range(ldos.shape[0])])
        sio.savemat('E:/dwell3/172/coord.mat', {'coord':coord},oned_as='row') 
    myDict = {'gf_mode':gf_mode, 'wf':wf, 'ld':ldos}  #'s':s, 
    completeName = os.path.join('E:/dwell3/172/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row')    
#    return s

Parallel(n_jobs=8)(delayed(tau_n)(cishu=j) for j in np.arange(0,500,1)) 

#s=Parallel(n_jobs=8)(delayed(tau_n)(cishu=j) for j in np.arange(0,1000,1)) 
#myDict = {'s':s}
#completeName = 'E:/dwell3/142.mat'
#sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini