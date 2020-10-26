# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 10:29:07 2018

@author: ykang
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import kwant
#from kwant.digest import uniform
#from time import time
import scipy.io as sio
import os.path
import long_range_correlation_square as corr

#salt=1
#t_ini=time()
width=50
length=100

def onsite(site):
    x,y=site.pos
    return 4+pattern[int(x),int(y)]

def disk(pos):
    x,y=pos
    return 0<=x<length and 0<=y<width

def edge(pos):
    x,y=pos
    return 0<=y<width

def make_system():
    sys=kwant.Builder()
    lat=kwant.lattice.square()
    sys[lat.shape(disk,(0,0))]=onsite
    sys[lat.neighbors()]=-1
    
    lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
    lead[lat.shape(edge,(0,1))]=4
    lead[lat.neighbors()]=-1
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    return sys.finalized()
en=1

#save_path = 'C:/Users/ykang/Documents/correlation eigenchannel/0/'
for cishu in np.arange(0,1):

    pattern=corr.correlation_pattern_guassian(length,width,5,1.5,cishu+random.randint(0,10000))
#    pattern=corr.uniformrandom(length,width,1,cishu+random.randint(0,10000))
    sys=make_system()
    scat=kwant.smatrix(sys,en)
    gf_mode=scat.submatrix(1,0)
      
    u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
    wf=kwant.wave_function(sys,en)(0)
    T=scat.transmission(1,0)
    if s[0]>0.98:
#            v=sio.loadmat(completeName)['v']
        f_tau1=abs(np.array(sum(np.dot(np.matrix(np.diag(v[0,:])).H,wf))))**2
        f_tau1_reshape=f_tau1.reshape((length,width))
        plt.subplot(121)
        plt.imshow(pattern)
#        plt.imshow((pattern-.25)*(1-np.sign(pattern-.25)))
        plt.subplot(122)
        plt.imshow(f_tau1_reshape)


#elapsed=time()-t_ini