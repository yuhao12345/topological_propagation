# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 21:45:58 2018

@author: ykang
"""
# Aubry-andre model, sweep phi, topological edge state appears
import numpy as np
from numpy import pi,sqrt,linspace
import matplotlib.pyplot as plt
#from random import random
import kwant
import scipy.sparse.linalg as sla
from numpy import linalg as LA

b=2#np.sqrt(5)/2+0.5

dx=1
L=99
spec=[]
def onsite(site):
    n,=site.pos
    return .5*np.cos(2*pi*b*n+phi) 

for phi in linspace(0,2*pi,num=100):   
    lat=kwant.lattice.chain(dx)
    sys=kwant.Builder()      
    sys[(lat(i) for i in range(L))]=onsite
    sys[lat.neighbors()]= 1
    sys=sys.finalized()
    ham_mat = sys.hamiltonian_submatrix(sparse=False)
    vals,evecs= LA.eigh(ham_mat)
    spec.append(vals)
    #plt.plot(np.abs(evecs[:, 0])**2)
spec=np.array(spec)
plt.plot(spec)

#plt.plot(vals,'.')
#plt.plot(np.abs(evecs[:, 37])**2)
