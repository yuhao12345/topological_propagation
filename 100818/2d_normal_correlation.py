# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 20:48:56 2018

@author: user
"""
#https://kwant-project.org/doc/1/reference/generated/kwant.physics.StabilizedModes#kwant.physics.StabilizedModes
#all the propagating modes are normalized to carry unit current

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
import long_range_correlation

t_ini=time()

length=115  # toi calculate loc length, sample length=1e4
width=350

w=1.1   #root mean squar, equivalent to dis in homogeneous case
v=10  #correlation length

lat = kwant.lattice.square()

def make_system(length,width,salt):
    #correlation_pattern_guassian(t1,t2,v,w,salt)
    # t:size   v:correlation length   w: root mean square, control the magnitude of disorder  
    pattern=long_range_correlation.correlation_pattern_guassian(2*length-1,2*width-1,v,w,salt)+4
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    #    return abs(x) + abs(y)*np.sqrt(3) < 100
    #    return abs(x)+10>abs(y)*3 and abs(x)<100
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width #2
    def onsite(site):
        x,y=site.pos
        return pattern[int(x+length-1),int(y+width-1)]
#        if abs(y)<10:
#            return 4
#        else:
#            return dis*(uniform(repr(site),salt)-.5)+4
    
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


def tau_n(cishu):
    en=0.4
    salt=cishu+random.randint(1,2**31)
    sys=make_system(length,width,salt)
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
    s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False) #u, s, vh
    
#    wf=kwant.wave_function(sys,en)(0)
#    kwant.plotter.map(sys, (abs(wf[14])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

#    ldos = kwant.ldos(sys,en)
#    if cishu==0:
#        coord=np.array([sys.pos(i) for i in range(ldos.shape[0])])
#        sio.savemat('E:/dwell3/162/coord.mat', {'coord':coord},oned_as='row') 
#    myDict = {'s':s }  #'s':s, 'gf_mode':gf_mode, 'wf':wf, 'ld':ldos
#    completeName = os.path.join('E:/dwell3/162/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row')    
    return s

#Parallel(n_jobs=8)(delayed(tau_n)(cishu=j) for j in np.arange(0,1000,1)) 

s=Parallel(n_jobs=10)(delayed(tau_n)(cishu=j) for j in np.arange(0,100,1)) 
myDict = {'s':s}
completeName = 'E:/dwell3/220.mat'
sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini
#
#sys=make_system(length,width, 'ed')  
#en=0.4
#t=kwant.smatrix(sys,en).transmission(1,0)  
#wff=kwant.wave_function(sys,en)
#wf=wff(0)
#kwant.plotter.map(sys, abs(wf[0]),num_lead_cells=5,fig_size=(15, 10))  
#
#plt.figure()
##plt.imshow(abs(wf[0]).reshape(39,9))
#tt=abs(wf[0]).reshape(39,9)
#plt.plot(tt[0,])
#sum(tt[0,]**2)
#y=[abs(np.sin(3*np.pi/10*i))*np.sqrt(2/14.5) for i in np.arange(1,10)]
#plt.plot(y)

