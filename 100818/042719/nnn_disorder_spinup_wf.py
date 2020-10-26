# -*- coding: utf-8 -*-
"""
Created on Thu May  9 18:29:01 2019

@author: ykang
"""

import numpy as np
from matplotlib import pyplot
import scipy.io as sio
import os.path
import random
import kwant
from kwant.digest import uniform
from random import choices
from time import time
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()

width=30
length=25

#di=10  # disorder: w-di<y<w

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  

#alpha=0.5  # probability of push down 

def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1

    def nnn(site1, site2):
        return 1j*m2

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*4


    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn #1j *m2 * sz #
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width 
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    


#def gf_virtual(cishu):

syst=make_system(width, length, str(8900))
attach_lead(syst)
sys = syst.finalized()
#kwant.plot(sys, fig_size=(10, 3))
wff=kwant.wave_function(sys,.4)
wf=wff(0)
kwant.plotter.map(sys, (abs(wf[16])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)  #


elapsed=time()-t_ini