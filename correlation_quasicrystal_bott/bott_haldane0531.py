# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 13:13:37 2018

@author: user
"""
# periodic bc to calcuulate bott index, haldane model

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

salt='5555'
length=5

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def momentum_to_lattice(k):
    """Transform momentum to the basis of reciprocal lattice vectors.
    
    See https://en.wikipedia.org/wiki/Reciprocal_lattice#Generalization_of_a_dual_lattice
    """
    B = np.array(graphene.prim_vecs).T
    A = B.dot(np.linalg.inv(B.T.dot(B)))
    return np.linalg.solve(A, k)

def onsite(site):
    return (uniform(repr(site),salt)-0.5)*2


sys = kwant.Builder(kwant.TranslationalSymmetry(*graphene.prim_vecs*length)) #kwant.TranslationalSymmetry(*graphene.prim_vecs*20)
sys[graphene.shape((lambda pos: True), (0, 0))] = onsite
sys[graphene.neighbors(1)] = -1
sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = .2j
sys=wraparound(sys).finalized()
#kwant.plot(sys, fig_size=(10, 5))

lattice_k = momentum_to_lattice([0, 0])
ham_mat = sys.hamiltonian_submatrix(args=(list(lattice_k)))

##vals,evecs= sla.eigsh(ham_mat, k=1,sigma=None,return_eigenvectors=True)  #,evecs
vals,evecs=LA.eigh(ham_mat)
#kwant.plotter.map(sys, np.abs(evecs[:, 100])**2,colorbar=False, oversampling=1)

#vals,evecs=scipy.sparse.linalg.eigs(ham_mat, k=length**2, M=None, which='SR')

p=[]
for i in range(2*length**2):
    p.append(sys.pos(i))
p=np.array(p)

l=sum(vals<1)
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


#pyplot.plot(vals,'.')