# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:23:41 2019

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

width=10
length=50


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
        return -3<y<width and abs(x)<length   #25.1

    def onsite(site):
        x,y=site.pos
#        return 0
        if (x+10)**2+(y-4)**2<16 or (x-20)**2+(y-5)**2<9:
            return 3
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
        return -3<y<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
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
syst=make_system(width, length, '2')   # whole system as virtual lead
attach_lead(syst)
##kwant.plot(syst0,fig_size=(25, 10))
#
##greens_function_sites = syst0.sites()
##mount_vlead(syst0, greens_function_sites, 1)
sys = syst.finalized()
#kwant.plot(sys)
### step 1, gf spectrum

ens=np.linspace(0.3,0.38,500)
def gf_01(cishu):
    en=ens[cishu]
    gf=kwant.greens_function(sys,en).submatrix(1,0)[0,0]
#    gf=kwant.smatrix(sys,en).submatrix(1,0)[0,0]
#    if gf.shape[0]==0:
#        gf=2
#    else:
#        gf=gf[0,0]
#    myDict = {'gf':gf} #,'ld':ld
#    completeName = os.path.join('E:/dwell3/756/', str(cishu)+'.mat')
#    sio.savemat(completeName,myDict,oned_as='row') 
    return gf
spec=Parallel(n_jobs=10)(delayed(gf_01)(cishu) for cishu in np.arange(0,500,1))
myDict = {'spec':spec,'ens':ens} #,'ld':ld
completeName = os.path.join('E:/dwell3/903.mat')
sio.savemat(completeName,myDict,oned_as='row') 
    
#def g_01(cishu):
#    en=ens[cishu]
#    g=kwant.smatrix(sys,en).transmission(1,0) 
#    return g

#g=Parallel(n_jobs=10)(delayed(g_01)(cishu) for cishu in np.arange(0,512,1))
#
#myDict = {'g':g} #,'ld':ld
#completeName = os.path.join('E:/dwell3/770.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#gf_list=Parallel(n_jobs=10)(delayed(gf_01)(cishu) for cishu in range(200))
#pyplot.plot(np.abs(gf_list)**2,'.')

#myDict = {'energies':energies,'gf':gf} #,'ld':ld
#completeName = os.path.join('E:/dwell3/751/1.mat')
#sio.savemat(completeName,myDict,oned_as='row') 


### step2  wf
en=ens[442]
wf=kwant.wave_function(sys,en)(0)   #.403586
kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

elapsed=time()-t_ini