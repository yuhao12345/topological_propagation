# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 15:18:33 2018

@author: ykang
"""
# close bc

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from numpy import sqrt,conj,exp,zeros,pi,angle,cos,sin,arccos,linspace,arctan
import kwant
import scipy.sparse.linalg as sla
from wraparound import wraparound
from numpy import linalg as LA
import random

p1=1
p2=0
theta=arctan(1/sqrt(3)) #arccos(11/14)
b1=1#np.sqrt(5)/2+0.5
length=10
phi=0
spec=[]
ksi=[]
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def disk(pos):
    x,y=pos
    return abs(x)<length and abs(y)<length

def nnn_hopping(site1, site2):
    return 1j * 1

def v1(x,y):
    return cos(2*pi*b1*x+phi)+2*cos(pi*b1*x-phi)*cos(sqrt(3)*pi*b1*y)
    

def potential(x,y):
    return p1*v1(x,y)+p2*v1(cos(theta)*x-sin(theta)*y,sin(theta)*x+cos(theta)*y)

def onsite(site):
    x,y=site.pos
    return (potential(x,y) if site.family == a else -potential(x,y-1/sqrt(3)))

#for phi in linspace(0,pi,num=1): 
for p1 in linspace(0,3,num=31):
    sys = kwant.Builder() #kwant.TranslationalSymmetry(*graphene.prim_vecs*20)
    sys[graphene.shape(disk, (0, 0))] = onsite
    sys[graphene.neighbors(1)] = 1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn_hopping
    sys=sys.finalized()
        
    ham_mat = sys.hamiltonian_submatrix(sparse=False)
    #vals,evecs= sla.eigsh(ham_mat,k=10, sigma=None,which='SM', return_eigenvectors=True)  #,evecs
    vals,evecs= LA.eigh(ham_mat)
    #vals=LA.eigvals(ham_mat)
    plt.figure()
    kwant.plotter.map(sys, np.abs(evecs[:,500])**2,colorbar=False, oversampling=1)
    #vals.shape[0]//2
    #plt.plot(vals,'.')
    spec.append(vals)
    ksi.append(sum(abs(evecs[:,500])**4)**0.5)
#spec=np.array(spec)
plt.plot(spec)
plt.figure()
plt.plot(ksi,'.')
#for i in range(41):
#    plt.figure()
#    plt.plot(spec[i],'.')
#plt.plot(spec[30],'.')    
#temp=np.abs(evecs[:, 500])**2
#for i in range(100):
#    temp=temp+np.abs(evecs[:, 499-i])**2
#kwant.plotter.map(sys, temp,colorbar=False, oversampling=1)    


#kwant.plot(sys, fig_size=(10, 5))

