# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 10:39:16 2017

@author: user   for QVH (onsite energy)
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

m1,m2=.5,.5  #m1 valley, m2 spin


graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
A, B = graphene.sublattices

nnn_hoppings_a = (((-1, 0), A, A), ((0, 1), A, A), ((1, -1), A, A))
nnn_hoppings_b = (((1, 0), B, B), ((0, -1), B, B), ((-1, 1), B, B))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

s0=np.identity(2)
#sx=np.array([[0,1],[1,0]])
sz=np.array([[1,0],[0,-1]])

def spin_orbit(site1, site2):
    return 1j * m2 * sz
def add_hoppings(syst):
    syst[graphene.neighbors(1)] = s0
    syst[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit
    
# for QSH AND QVH
def onsite_qsvh(site):
    x,y=site.pos
    if y<0:
        onsite_a = m1*s0
        onsite_b = -m1*s0
        return onsite_a if site.family == A else onsite_b
    else:
        return np.zeros([2,2])
def hopp_qsvh(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    hop_a = m2*1j*sz
    hop_b = -m2*1j*sz
    if y1>0:
        return hop_a if site1.family == A else hop_b
    else:
        return np.zeros([2,2])

# for QSH
def onsite_qsh(site):
#    onsite_a = m1*s0
#    onsite_b = -m1*s0
#    return onsite_a if site.family == A else onsite_b
    return np.zeros([2,2])    
def hopp_qsh(site1,site2):
    hop_a = m2*1j*sz
    hop_b = -m2*1j*sz
    return hop_a if site1.family == A else hop_b

# for QVH
def onsite_qvh(site):
    onsite_a = m1*s0
    onsite_b = -m1*s0
    return onsite_a if site.family == A else onsite_b

# for QVH-QSH-QVH
def onsite_qvsvh(site):
    x,y=site.pos
    if abs(y)>15:
        onsite_a = m1*s0
        onsite_b = -m1*s0
        return onsite_a if site.family == A else onsite_b
    else:
        return np.zeros([2,2])
def hopp_qvsvh(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    hop_a = m2*1j*sz
    hop_b = -m2*1j*sz
    if abs(y1)>15:
        return np.zeros([2,2])       
    else:
        return hop_a if site1.family == A else hop_b
    
# for QSH-QSH-QSH
def hopp_qssh(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    hop_a = m2*1j*sz
    hop_b = -m2*1j*sz
    if y1>0:
        return hop_a if site1.family == A else hop_b
    else:
        return hop_b if site1.family == A else hop_a
    
lead = kwant.Builder(kwant.TranslationalSymmetry([1, 0]))
lead[graphene.shape((lambda pos: abs(pos[1]) < 8), (0, 0))] = onsite_qsh

#lead=kwant.Builder(kwant.TranslationalSymmetry([1,0])) 
# QSH
#lead[graphene.shape(edge,(0,0))]=onsite_qsh
#lead[graphene.neighbors()]=-s0
#lead[A.neighbors()]=hopp_qsh
#lead[B.neighbors()]=hopp_qsh
#kwant.plot(lead,fig_size=(10, 12))
# QVH
#lead[graphene.shape(edge,(0,0))]=onsite_qvh
lead[graphene.neighbors()]=-s0

# QVH+QSH
#lead[graphene.shape(edge,(0,0))]=onsite_qsvh
#lead[graphene.neighbors()]=-s0
#lead[A.neighbors()]=hopp_qsvh
#lead[B.neighbors()]=hopp_qsvh

#QVH-QSH-QVH
#lead[graphene.shape(edge,(0,0))]=onsite_qvsvh
#lead[graphene.neighbors()]=-s0
#lead[A.neighbors()]=hopp_qvsvh
#lead[B.neighbors()]=hopp_qvsvh

# QSH-QSH-QSH
#lead[graphene.shape(edge,(0,0))]=np.zeros([2,2])
#lead[graphene.neighbors()]=-s0
#lead[A.neighbors()]=hopp_qssh
#lead[B.neighbors()]=hopp_qssh

#kwant.plot(lead,fig_size=(12, 30))
#plt.hold(True)
#add_hoppings(lead)
kwant.plotter.bands(lead.finalized(), fig_size=(10, 12))

#ax.add_line(mlines.Line2D([-np.pi,np.pi],[.3,.3]))
