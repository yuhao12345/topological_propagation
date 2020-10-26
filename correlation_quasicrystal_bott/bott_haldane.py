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
import random

p1=0   # ~0.9 is the phase transition, bandgap reopen
p2=3 # ~2.1 is phase transition, bandgap disappear since onsite potential is aperiodic
theta=arccos(11/14)
length=10

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

# only a slice of band structure
def dispersion_2D(syst, args=None, lim=4*np.pi/length, num_points=15):
    """A simple plot of 2D band structure."""
    if args is None:
        args = []
    momenta0=[0]
    momenta = np.linspace(-lim, lim, num_points)
    energies = []
    for kx in momenta:
        for ky in momenta0:
            lattice_k = momentum_to_lattice([kx, ky])
            h = syst.hamiltonian_submatrix(args=(list(args) + list(lattice_k)))
            energies.append(np.linalg.eigvalsh(h))
    
#    energies = np.array(energies).reshape(num_points, num_points, -1)
    energies = np.array(energies).reshape( num_points,1, -1)
    emin, emax = np.min(energies), np.max(energies)
    kx, ky = np.meshgrid(momenta, momenta0)
    fig = pyplot.figure()
    axes = fig.add_subplot(1, 1, 1, projection='3d')
    for band in range(energies.shape[-1]):
#        axes.plot_surface(kx, ky, energies[:, :, band], cstride=2, rstride=2,
#                          cmap=matplotlib.cm.RdBu_r, vmin=emin, vmax=emax,
#                          linewidth=0.1)
        axes.scatter(kx, ky, energies[:, :, band])
    axes.view_init(0, azim=90)
def nnn_hopping(site1, site2):
    return 1j * 1

def v1(x,y):
    return cos(2*pi*x)+2*cos(pi*x)*cos(sqrt(3)*pi*y)

def potential(x,y):
    return p1*v1(x,y)+p2*v1(cos(theta)*x-sin(theta)*y,sin(theta)*x+cos(theta)*y)

def onsite(site):
    x,y=site.pos
    return (potential(x,y) if site.family == a else -potential(x,y-1/sqrt(3)))
#    return 0 * (1 if site.family == a else -1)

sys = kwant.Builder(kwant.TranslationalSymmetry(*graphene.prim_vecs*length)) #kwant.TranslationalSymmetry(*graphene.prim_vecs*20)
sys[graphene.shape((lambda pos: True), (0, 0))] = onsite
sys[graphene.neighbors(1)] = 1
sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn_hopping
sys=wraparound(sys).finalized()
kwant.plot(sys, fig_size=(10, 5))


#dispersion_2D(sys)
#
#
#lattice_k = momentum_to_lattice([0, 0])
#ham_mat = sys.hamiltonian_submatrix(args=(list(lattice_k)))
#
###vals,evecs= sla.eigsh(ham_mat, k=1,sigma=None,return_eigenvectors=True)  #,evecs
#vals,evecs=LA.eigh(ham_mat)
##kwant.plotter.map(sys, np.abs(evecs[:, 100])**2,colorbar=False, oversampling=1)
#
#p=[]
#for i in range(2*length**2):
#    p.append(sys.pos(i))
#p=np.array(p)
#
#l=length**2
#c01=zeros([l,l],dtype=np.complex)
#for i in range(l):
#    for k in range(l):
#        c01[i,k]=sum(conj(evecs[:, i])*exp(-1j*2*pi/length*p[:,0])*evecs[:, k])
#c12=zeros([l,l],dtype=np.complex)
#for i in range(l):
#    for k in range(l):
#        c12[i,k]=sum(conj(evecs[:, i])*exp(-1j*pi/length*(p[:,0]+sqrt(3)*p[:,1]))*evecs[:, k])
#c23=zeros([l,l],dtype=np.complex)
#for i in range(l):
#    for k in range(l):
#        c23[i,k]=sum(conj(evecs[:, i])*exp(1j*2*pi/length*p[:,0])*evecs[:, k])
#c30=zeros([l,l],dtype=np.complex)
#for i in range(l):
#    for k in range(l):
#        c30[i,k]=sum(conj(evecs[:, i])*exp(1j*pi/length*(p[:,0]+sqrt(3)*p[:,1]))*evecs[:, k])
#c = c01 @ c12 @ c23 @ c30
#chern=sum(angle(LA.eigvals(c)))/2/pi


#pyplot.plot(vals,'.')