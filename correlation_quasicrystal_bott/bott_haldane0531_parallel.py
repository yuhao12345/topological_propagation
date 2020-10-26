# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:12:59 2019

@author: ykang
"""

import matplotlib
from matplotlib import pyplot
from mpl_toolkits import mplot3d
import numpy as np
from numpy import sqrt,conj,exp,zeros,pi,angle,cos,sin,arccos
import kwant
import scipy.sparse.linalg as sla
from wraparound import wraparound
from numpy import linalg as LA
from kwant.digest import uniform
import scipy
from joblib import Parallel, delayed
import multiprocessing
import random
from time import time
import scipy.io as sio
import os.path

t_ini=time()

length=15

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def make_system(length,salt):
    def onsite(site):
        return (uniform(repr(site),salt)-0.5)*2.5  
    sys = kwant.Builder(kwant.TranslationalSymmetry(*graphene.prim_vecs*length)) #kwant.TranslationalSymmetry(*graphene.prim_vecs*20)
    sys[graphene.shape((lambda pos: True), (0, 0))] = onsite
    sys[graphene.neighbors(1)] = -1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = .1j
    sys=wraparound(sys).finalized()
    return sys

def coord():
    sys=make_system(length,'0')  
    p=[]
    for i in range(2*length**2):
        p.append(sys.pos(i))
    p=np.array(p)
    return p

p=coord()

def bott(cishu):
    salt=cishu+random.random()*100
    sys=make_system(length,str(salt))

    lattice_k = np.array([0,0])
    ham_mat = sys.hamiltonian_submatrix(args=(list(lattice_k)))

    vals,evecs=LA.eigh(ham_mat)
    
    l=sum(vals<-0.4)
    c01=zeros([l,l],dtype=np.complex)
    for i in range(l):
        for k in range(l):
            c01[i,k]=sum(conj(evecs[:, i])*exp(-1j*2*pi/length*p[:,0])*evecs[:, k])
    c12=zeros([l,l],dtype=np.complex)
    for i in range(l):
        for k in range(l):
            c12[i,k]=sum(conj(evecs[:, i])*exp(-1j*pi/length*(p[:,0]+sqrt(3)*p[:,1]))*evecs[:, k])
    c23=zeros([l,l],dtype=np.complex)
    for i in range(l):
        for k in range(l):
            c23[i,k]=sum(conj(evecs[:, i])*exp(1j*2*pi/length*p[:,0])*evecs[:, k])
    c30=zeros([l,l],dtype=np.complex)
    for i in range(l):
        for k in range(l):
            c30[i,k]=sum(conj(evecs[:, i])*exp(1j*pi/length*(p[:,0]+sqrt(3)*p[:,1]))*evecs[:, k])
    c = c01 @ c12 @ c23 @ c30
    chern=sum(angle(LA.eigvals(c)))/2/pi
    return chern

bo=Parallel(n_jobs=5)(delayed(bott)(n) for n in np.arange(0,200))
bott_index=np.mean(bo)

myDict = {'bott':bo}  
sio.savemat('E:/dwell3/bott2.5.mat',myDict,oned_as='row')    

elapsed=time()-t_ini