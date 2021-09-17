# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 10:39:16 2017

@author: user
"""

from types import SimpleNamespace

from ipywidgets import interact
import matplotlib
from matplotlib import pyplot
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import kwant
from wraparound import wraparound


def momentum_to_lattice(k):
    """Transform momentum to the basis of reciprocal lattice vectors.
    
    See https://en.wikipedia.org/wiki/Reciprocal_lattice#Generalization_of_a_dual_lattice
    """
    B = np.array(graphene.prim_vecs).T
    A = B.dot(np.linalg.inv(B.T.dot(B)))
    return np.linalg.solve(A, k)


def dispersion_2D(syst, args=None, lim=1.5*np.pi, num_points=200):
    """A simple plot of 2D band structure."""
    if args is None:
        args = []
    momenta = np.linspace(-lim, lim, num_points)
    energies = []
    for kx in momenta:
        for ky in momenta:
            lattice_k = momentum_to_lattice([kx, ky])
            h = syst.hamiltonian_submatrix(args=(list(args) + list(lattice_k)))
            energies.append(np.linalg.eigvalsh(h))
    
    energies = np.array(energies).reshape(num_points, num_points, -1)
    emin, emax = np.min(energies), np.max(energies)
    kx, ky = np.meshgrid(momenta, momenta)
    fig = pyplot.figure()
    axes = fig.add_subplot(1, 1, 1, projection='3d')
    for band in range(energies.shape[-1]):
        axes.plot_surface(kx, ky, energies[:, :, band], cstride=2, rstride=2,
                          cmap=matplotlib.cm.RdBu_r, vmin=emin, vmax=emax,
                          linewidth=0.1)
        
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

s0=np.identity(2)
sx=np.array([[0,1],[1,0]])
sz=np.array([[1,0],[0,-1]])
ho=np.array([.2, -.2, 0])
def onsite1(site):
    x,y=site.pos
    t=[0 if y>=0 else 1]
    return ho[t]*sz+.1*sx
def onsite2(site):
    x,y=site.pos
    t=[0 if y>=0 else 1]
    return -ho[t]*sz+.1*sx
#def onsite1a(site):
#    x,y=site.pos
#    if y>0.7:
#        t=0
#    if y<0:
#        t=1
#    else:
#        t=2
#    return ho[t]*sz+.4*sx
#def onsite2a(site):
#    x,y=site.pos
#    if y>0.7:
#        t=0
#    if y<0:
#        t=1
#    else:
#        t=2
#    return -ho[t]*sz+.4*sx
def hopp(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1)*(y2+0.1)>0 and y1+0.1>0:
        return .15*s0
    if (y1+0.1)*(y2+0.1)>0 and y1+0.1<0:
        return -.15*s0
    else:
        return np.zeros([2,2])    
def hopp2(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1)*(y2+0.1)>0 and y1+0.1>0:
        return -.15*s0
    if (y1+0.1)*(y2+0.1)>0 and y1+0.1<0:
        return .15*s0
    else:
        return np.zeros([2,2])
m, t=0.5, -0.133
def onsite3(site):
    x,y=site.pos
    if y>-0.1:
        return m*sz+.0*sx
    else:
        return np.zeros([2,2])
def onsite4(site):
    x,y=site.pos
    if y>-0.1:
        return -m*sz+.0*sx
    else:
        return np.zeros([2,2])  
def hopp3(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1)*(y2+0.1)>0 and y1+0.1<0:
        return -t*s0
    else:
        return np.zeros([2,2])
def hopp4(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1)*(y2+0.1)>0 and y1+0.1<0:
        return t*s0
    else:
        return np.zeros([2,2])
def edge(pos):
    x,y=pos
    return abs(y)<30 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

#armchair_ribbon = kwant.Builder(kwant.TranslationalSymmetry([0, np.sqrt(3)]))
#armchair_ribbon[graphene.shape((lambda pos: abs(pos[0]) < 9), (0, 0))] = 0
#armchair_ribbon[graphene.neighbors(1)] = 1

lead=kwant.Builder(kwant.TranslationalSymmetry((1,0))) 
# QSH
#lead[a.shape(edge,(0,0))]=onsite1
#lead[b.shape(edge,(0,0))]=onsite2
#lead[graphene.neighbors()]=-s0*2/3
# QVH
#lead[graphene.shape(edge,(0,0))]=np.zeros([2,2])
#lead[graphene.neighbors()]=-s0*2/3
#lead[a.neighbors()]=hopp
#lead[b.neighbors()]=hopp2
# QVH+QSH
lead[a.shape(edge,(0,0))]=onsite3
lead[b.shape(edge,(0,0))]=onsite4
lead[graphene.neighbors()]=-s0*2/3
lead[a.neighbors()]=hopp3
lead[b.neighbors()]=hopp4
#kwant.plot(lead,fig_size=(12, 30))
#plt.hold(True)
kwant.plotter.bands(lead.finalized(), fig_size=(10, 12));

#ax.add_line(mlines.Line2D([-np.pi,np.pi],[.3,.3]))
