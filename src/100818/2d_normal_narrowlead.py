# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 16:51:41 2019

@author: ykang
"""
# barrow lead, so input only one channel

import kwant
from random import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform
import scipy.io as sio
import os.path
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()

length=40  # toi calculate loc length, sample length=1e4
width=10
lead_width=10
dis=.5

lat = kwant.lattice.square()

def make_system(length,width,lead_width,salt):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    def lead_shape(pos):
        x, y= pos
        return abs(y)<lead_width
    def onsite(site):
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

def lnt_x1(cishu):
    en=0.4
    salt=cishu+random()
    sys=make_system(length,width, lead_width,str(salt))
    sm=kwant.smatrix(sys,en)
    gf_mode=sm.submatrix(1,0) 
    t=sm.transmission(1,0) 
    wff=kwant.wave_function(sys,en)(0)
    w=np.abs(wff.T)
    coord=np.array([sys.pos(i) for i in range(w.shape[0])])
    u=np.concatenate((coord,w),axis=1)

    myDict = {'u':u,'t':t,'g':gf_mode}
    completeName = os.path.join('E:/narrowlead/11/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
    
#Parallel(n_jobs=10)(delayed(lnt_x1)(cishu=j) for j in np.arange(0,5000,1)) 

sys=make_system(length,width,lead_width, 'ed')  
en=0.4
wf=kwant.wave_function(sys,en)(0)
kwant.plotter.map(sys, abs(wf[0]),num_lead_cells=5,fig_size=(15, 10))  
t=kwant.smatrix(sys,en).transmission(1,0)

elapsed=time()-t_ini
