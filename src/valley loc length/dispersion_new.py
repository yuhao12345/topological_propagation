# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 16:38:08 2018

@author: user
"""

import matplotlib
from matplotlib import pyplot
from mpl_toolkits import mplot3d
import numpy as np
import scipy.io as sio
import os.path

import kwant
from wraparound import wraparound

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m1 = .52  #valley   0.4
m2 = .1  #spin    0.077

s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def spin_orbit(site1, site2):
    return 1j * m2 * sz

def spin_orbit_v(site1, site2):
    x,y=site1.pos
    if y>-0.1:     #y>-0.1:
        return 1j * m2 * sz
    else:
        return np.zeros([2,2])

def onsite_s(site):
    return np.zeros([2,2])

def onsite_v(site):
    return s0 * m1 * (1 if site.family == a else -1)

def onsite_sv(site):
    x,y=site.pos
    if y>-0.1:     #y>-0.1:
        return np.zeros([2,2])
    else:
        return s0 * m1 * (1 if site.family == a else -1)

def free_space(syst):
    syst[graphene.neighbors(1)] = s0
    
def qsh(syst):
    syst[graphene.neighbors(1)] = s0
    syst[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit

def qvh(syst):
    syst[graphene.neighbors(1)] = -s0  
    
def qsvh(syst):
    syst[graphene.neighbors(1)] = -s0
    syst[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v

zigzag = kwant.Builder(kwant.TranslationalSymmetry([1, 0]))

zigzag[graphene.shape((lambda pos: abs(pos[1]) < 15.1), (0, 0))] = onsite_s  #8, chuxian bandgap he zhixian
qsh(zigzag)    #9, 15.1, 10,19.4  zigzag two side
#free_space(zigzag)
#zigzag[graphene.shape((lambda pos: abs(pos[1]) <10.6), (0, 0))] = onsite_v  #8, chuxian bandgap he zhixian
#qvh(zigzag)   #18.5

#zigzag[graphene.shape((lambda pos: abs(pos[1]) < 19.4), (0, 0))] = onsite_sv  #8, chuxian bandgap he zhixian
#qsvh(zigzag)   #19.4

def family_color(site):
    return 'black' if site.family == a else 'white'

def hopping_lw(site1, site2):
    return 0.04 if site1.family == site2.family else 0.1

#kwant.plot(zigzag,fig_size=(10, 12),site_lw=0.1,site_color=family_color,hop_lw=hopping_lw)
#kwant.plotter.bands(zigzag.finalized(), fig_size=(10, 12))

#bands = kwant.physics.Bands(zigzag)
#momenta = np.linspace(-np.pi, np.pi, 101)
#energies = [bands(k) for k in momenta]
#pyplot.plot(momenta, energies)
#pyplot.show()

bands = kwant.physics.Bands(zigzag.finalized())
momenta = np.linspace(-np.pi, np.pi, 601)
energies = [bands(k) for k in momenta]
myDict = {'k':momenta,'en':energies}
completeName = 'C:/Users/ykang/Dropbox/TI8/dispersion.mat'
sio.savemat(completeName,myDict,oned_as='row') 

#pyplot.plot(momenta, energies)
#pyplot.show()