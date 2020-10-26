# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:30:39 2018

@author: ykang
"""
# haldane model, long range correlated disorder, calculate loc length
# periodic along transverse direction

#import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d
import numpy as np
from numpy import sqrt,conj,exp,zeros,pi,angle,cos,sin,arccos,linspace,arctan
import kwant
#import scipy.sparse.linalg as sla
import wraparound
#from numpy import linalg as LA
import scipy.io as sio
import os.path
import math
import long_range_correlation
from time import time
#from kwant.digest import uniform
#from random import random
t_ini=time()

Y=10
length=2500
width=sqrt(3)/2*Y

#salt='55'
#dis=2

t1=2*int(2**(math.ceil(np.log2(length))))
t2=int(2**(math.ceil(np.log2(3*Y))))

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def disk(pos):
    x,y=pos
    return 0<=x<length 

def nnn_hopping(site1, site2):
    return .15j 

def hopping1(site1,site2):
    return -1
def onsite_clean(site):
    return 0.5*(1 if site.family == a else -1)
def onsite(site):
    x,y=site.pos
#    return 0
#    return (uniform(repr(site),str(salt))-0.5)*dis+0.5*(1 if site.family == a else -1)
#    return (random()-0.5)*dis+0.5*(1 if site.family == a else -1)
    return pattern[int(round(2*x)),int(round(2*sqrt(3)*y))]+0.5*(1 if site.family == a else -1)

 
# Infinite potential plane in y direction
def make_system():
    sys = kwant.Builder(kwant.TranslationalSymmetry(graphene.vec((-Y/2,Y))))
    sys[graphene.shape(disk,(0, 0))] = onsite
    sys[graphene.neighbors(1)] = hopping1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn_hopping
    sys = wraparound.wraparound(sys)
    #kwant.plot(sys, fig_size=(10, 5))
    
    lead = kwant.Builder(kwant.TranslationalSymmetry((1, 0), graphene.vec((-Y/2,Y))))
    lead[graphene.shape(disk,(0, 0))] = onsite_clean
    lead[graphene.neighbors(1)] = hopping1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn_hopping
    lead = wraparound.wraparound(lead, keep=0)
     
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
     
    return sys.finalized()

ky = 0.0


energies = np.linspace(0.01, 3, 10)
#transmission = []
#save_path = 'C:/Users/ykang/Documents/correlated_loc_downdstair/34/'
save_path = 'E:/haldane_correlated/g10/'

v=2   
w=4

for salt in np.arange(5):
    seed=np.random.randint(10000)
    pattern=long_range_correlation.correlation_pattern_guassian(t1,t2,v,w,seed)
#    pattern=long_range_correlation.correlation_pattern(t1,t2,v,w,seed)
    os.makedirs(os.path.join(save_path,str(salt)+'/'))

    sys=make_system()
    cishu=0 
    for energy in energies:
        try:
            gf=kwant.greens_function(sys,energy, [ky]).submatrix(1,0)
#            gf_mode = kwant.smatrix(sys, energy, [ky])
#            T=gf_mode.transmission(1,0)
#            transmission.append(T)
            
            completeName = os.path.join(save_path+str(salt)+'/',str(cishu)+".mat")
            myDict = {'energy':energy, 'gf':gf, 'length':length,'width':width}
            sio.savemat(completeName,myDict,oned_as='row')
            
            cishu=cishu+1
        except:
            continue

elapsed=time()-t_ini
#2*sqrt(3)*sys.pos(39)[1]
