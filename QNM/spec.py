# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:00:35 2019

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
            return 4
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

#sys=make_system(length,width,'5')
#
#num_f=4500
#energies=np.linspace(4.5,5.4,num_f)
#gf=np.zeros((num_f,2*width-1,2*width-1),dtype=complex)
#for k in range(num_f):
#    gf[k,:,:]=kwant.greens_function(sys,energies[k]).submatrix(1,0)
#    
#plt.figure(figsize=(20,10))
#plt.plot(energies,abs(gf[:,5,9])**2,'.')
#plt.plot(np.real(w),np.zeros(w.shape),'o')
#
#myDict = {'gf':gf,'en':energies,'w':w,'v':v}
#completeName = 'E:/QNM/0.mat'
#sio.savemat(completeName,myDict,oned_as='row') 
#
wff=kwant.wave_function(sys,4.873)
wf=wff(0)
kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=2,fig_size=(15, 8),colorbar=False)

kwant.plotter.map(sys, (abs(v[:,61])**2),num_lead_cells=2,fig_size=(15, 8),colorbar=False)
#kwant.plotter.map(sys, (abs(v[:,51])**2),num_lead_cells=2,fig_size=(15, 8),colorbar=False)

elapsed=time()-t_ini