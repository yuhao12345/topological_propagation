# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 19:48:46 2019

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

width=30
length=50

di=10  # disorder: w-di<y<w
alpha=0.5

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  

def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width #width-di<y<width  #
    def nnn(site1, site2):
        x,y=site1.pos
        if abs(x)<length and y>width-di:
            return 1j*m2*np.sign(uniform(repr(site1.tag),salt)-alpha)
        else:
            return 1j *m2 

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= 0
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn #1j *m2 * sz #

    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 #1
    lead[graphene.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
    return sys.finalized()

df=1e-8   # It is freq not omega !!!
en=0.3
en1=en+df

#sys=make_system(width, length, '0')
#tm=kwant.smatrix(sys,en)
#gf_mode=tm.submatrix(1,0) 
#number=gf_mode.shape[0]
def eigentime(u0,u1,vh0,vh1):
    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t
def delay_tm(cishu):
    salt=cishu+random.random()*100
    sys=make_system(width, length, str(salt))
    gf=kwant.greens_function(sys,en).submatrix(1,0)
#    gf=kwant.smatrix(sys,en).submatrix(1,0)
    u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
    gf1=kwant.greens_function(sys,en1).submatrix(1,0)
#    gf1=kwant.smatrix(sys,en1).submatrix(1,0)
    u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
    
    t=eigentime(u0,u1,vh0,vh1)
    
    g=kwant.smatrix(sys,en).transmission(1,0)
    return np.concatenate((np.diag(t)[0:5],s0[0:5],[g]))

#dwell=delay_tm(0)
dwell=Parallel(n_jobs=6)(delayed(delay_tm)(cishu=j) for j in np.arange(0,10000,1))

myDict = {'dwell':dwell}
completeName = 'E:/dwell3/527.mat'
sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini