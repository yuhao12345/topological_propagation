# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 23:03:06 2019

@author: user
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

width=20
length=30

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m1 = .52
dis=1.5

def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1

    def onsite(site):
        x,y=site.pos
        if abs(y)>1.5:
            return (uniform(repr(site),salt)-0.5)*dis+m1 * (1 if site.family == a else -1)
        else:
            return 0

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width 
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym) 
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

#syst=make_system(width, length, '0')
#attach_lead(syst)
#sys = syst.finalized()
#wf=kwant.wave_function(sys,0.4)(0)
#kwant.plotter.map(sys, (abs(wf[16])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
#   
def TM(cishu):
    en=0.4
    salt=cishu+random.random()*100
    syst=make_system(width, length, str(salt))
    attach_lead(syst)
    sys = syst.finalized()
    wf=kwant.wave_function(sys,en)(0)
  
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#    s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False) #u, s, vh
    ldos = kwant.ldos(sys, en)
    if cishu==0:
        coord=np.array([sys.pos(i) for i in range(ldos.shape[0])])
        sio.savemat('E:/dwell3/703/coord.mat', {'coord':coord},oned_as='row') 
    myDict = {'gf_mode':gf_mode, 'wf':wf, 'ld':ldos}  #'s':s, 
    completeName = os.path.join('E:/dwell3/703/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row')      

#TM(0)
Parallel(n_jobs=10)(delayed(TM)(n) for n in np.arange(0,100))

#dwell=delay_gf(0)
#dwell=Parallel(n_jobs=12)(delayed(delay_gf)(cishu=j) for j in np.arange(0,1000)) 
#myDict = {'s':s}

#myDict={'dwell':dwell}
#completeName = 'E:/dwell3/106.mat'
#sio.savemat(completeName,myDict,oned_as='row') 


elapsed=time()-t_ini