# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:19:44 2019

@author: ykang
"""

import numpy as np
from matplotlib import pyplot
import scipy.io as sio
import os.path
import random
import kwant
from kwant.digest import uniform
from random import choices
from time import time
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()


length=50  # toi calculate loc length, sample length=1e4
width=15

dis=6

lat = kwant.lattice.square()

def make_system(length,width,salt):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width 
    def onsite(site):
        x,y=site.pos
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

df=1e-10   # It is freq not omega !!!

en=0.41
en1=en+df

sys=make_system(length,width, '0')
tm=kwant.smatrix(sys,en)
gf_mode=tm.submatrix(1,0) 
number=gf_mode.shape[0]





#  gf=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!
def eigentime(u0,u1,vh0,vh1):

    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t
def delay_tm(cishu):
    salt=cishu+random.random()*100
    sys=make_system(length,width, str(salt))
    gf=kwant.greens_function(sys,en).submatrix(1,0)
#    gf=kwant.smatrix(sys,en).submatrix(1,0)
    u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
    gf1=kwant.greens_function(sys,en1).submatrix(1,0)
#    gf1=kwant.smatrix(sys,en1).submatrix(1,0)
    u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
    
    t=eigentime(u0,u1,vh0,vh1)
    return np.concatenate((np.diag(t)[0:number],s0[0:number]))

#dwell=delay_tm(0)
dwell=Parallel(n_jobs=8)(delayed(delay_tm)(cishu=j) for j in np.arange(0,10000,1))

myDict = {'dwell':dwell}
completeName = 'E:/dwell3/522.mat'
sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini