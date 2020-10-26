# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:34:19 2019

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
length=200


dis=0

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  


def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1

    def onsite(site):
        x,y=site.pos
        if y>width-10:
            return (uniform(repr(site),salt)-0.5)*dis
        else:
            return 0


    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width 
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    

#def mount_vlead(sys, vlead_interface, norb):
#    """Mounts virtual lead to interfaces provided.
#
#    :sys: kwant.builder.Builder
#        An unfinalized system to mount leads
#    :vlead_interface: sequence of kwant.builder.Site
#        Interface of lead
#    :norb: integer
#        Number of orbitals in system hamiltonian.
#    """
#    dim = len(vlead_interface)*norb
#    zero_array = np.zeros((dim, dim), dtype=float)
#    def selfenergy_func(energy, args=()):
#        return zero_array
#
#    vlead = kwant.builder.SelfEnergyLead(selfenergy_func, vlead_interface)
#    sys.leads.append(vlead)

########### for one configuration

#en=0.4
syst=make_system(width, length, '9')   # whole system as virtual lead
attach_lead(syst)
##kwant.plot(syst0,fig_size=(25, 10))
#
##greens_function_sites = syst0.sites()
##mount_vlead(syst0, greens_function_sites, 1)
sys = syst.finalized()
#kwant.plot(sys)

en=0.4008
G=kwant.smatrix(sys,en).transmission(1,0) 
wf=kwant.wave_function(sys,en)(0)
kwant.plotter.map(sys, (abs(wf[0])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
  
### step 1, gf spectrum

energies=np.linspace(0.35,0.45,1000)
def gf_01(cishu):
    en=energies[cishu]
    gf=kwant.greens_function(sys,en).submatrix(1,0)
    myDict = {'gf':gf} #,'ld':ld
    completeName = os.path.join('E:/dwell3/751/69/', str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
    return gf

Parallel(n_jobs=10)(delayed(gf_01)(cishu) for cishu in np.arange(0,1000,1))

elapsed=time()-t_ini