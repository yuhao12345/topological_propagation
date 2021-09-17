# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:50:29 2019

@author: ykang
"""
# periodic along x direction

import kwant
import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform
import scipy.io as sio
import os.path
import wraparound
import scipy.sparse.linalg as sla

length=3
width=2
salt='0'
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def disk(pos):
    x,y=pos
    return 0<=y<width # and 0<=x<length
def onsite(site):
    x,y=site.pos
    return  (uniform(repr(site),salt)-0.5)*4
sys = kwant.Builder(kwant.TranslationalSymmetry(graphene.vec((length,0))))
sys[graphene.shape(disk, (0, 0))] = onsite
sys[graphene.neighbors(1)] = 1
sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = .1j
#sys = wraparound.wraparound(sys).finalized()

kwant.plot(sys, fig_size=(10, 5))
sys=sys.finalized()
#
#def momentum_to_lattice(k):
#    """Transform momentum to the basis of reciprocal lattice vectors.
#    
#    See https://en.wikipedia.org/wiki/Reciprocal_lattice#Generalization_of_a_dual_lattice
#    """
#    B = np.array(graphene.prim_vecs).T
#    A = B.dot(np.linalg.inv(B.T.dot(B)))
#    return np.linalg.solve(A, k)
#lattice_k = momentum_to_lattice([0, 0])
#ham_mat = sys.hamiltonian_submatrix(args=(list(lattice_k)))

ham_mat = sys.hamiltonian_submatrix(sparse=True)
#evecs = sla.eigsh(ham_mat, k=20, which='SM')[1]

# Plot the probability density of the 10th eigenmode.
#kwant.plotter.map(sys, np.abs(evecs[:, 15])**2,colorbar=False, oversampling=1,fig_size=(10, 5))

w,v=sla.eigs(ham_mat, k=20, sigma=0.01, which='LR',return_eigenvectors=True)

def make():
    sys = kwant.Builder()
    def disk0(pos):
        x,y=pos
        return 0<=y<width and 0<=x<length
    sys[graphene.shape(disk0, (0, 0))] = onsite
    sys[graphene.neighbors(1)] = 1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = .1j
    return sys.finalized()
sys0=make()
kwant.plotter.map(sys0, (abs(v[:,1])**2),num_lead_cells=2,fig_size=(15, 10))
