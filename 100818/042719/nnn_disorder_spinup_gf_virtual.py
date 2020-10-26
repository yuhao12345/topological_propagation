# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:12:41 2019

@author: user
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

width=6.2
length=50

di=10  # disorder: w-di<y<w

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  

alpha=0.5  # probability of push down 

def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1

    def nnn(site1, site2):
        x,y=site1.pos
#        return 1j*m2
        if abs(x)<length and y>width-di:
#            return 1j * m2 * (1-2*uniform(repr(site1.tag),salt)*1)
            return 1j*m2*np.sign(uniform(repr(site1.tag),salt)-alpha)
        else:
            return 1j *m2 
    def onsite(site):
        x,y=site.pos
#        return (uniform(repr(site),salt)-0.5)*4
        return 0
#        if abs(x)<length and y>width-di:
#            return (uniform(repr(site),salt)-0.5)*1
#        else:
#            return 0

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
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    

def mount_vlead(sys, vlead_interface, norb):
    """Mounts virtual lead to interfaces provided.

    :sys: kwant.builder.Builder
        An unfinalized system to mount leads
    :vlead_interface: sequence of kwant.builder.Site
        Interface of lead
    :norb: integer
        Number of orbitals in system hamiltonian.
    """
    dim = len(vlead_interface)*norb
    zero_array = np.zeros((dim, dim), dtype=float)
    def selfenergy_func(energy, args=()):
        return zero_array

    vlead = kwant.builder.SelfEnergyLead(selfenergy_func, vlead_interface)
    sys.leads.append(vlead)

def make_system0(width, length):
    def disk(pos):
        x,y=pos
        return width-di<y<width and abs(x)<length   #25.1
    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,width-1))]= 0  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    def lead_shape(pos):
        x,y=pos
        return width-di<y<width 
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])    
    lead = kwant.Builder(sym) #,conservation_law=-sz    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 #1
    lead[graphene.neighbors()]=1       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys

syst0=make_system0(width, length)   # part of system



#syst0=make_system(width, length, '0')   # whole system as virtual lead
#attach_lead(syst0)
#kwant.plot(syst0,fig_size=(25, 10))

greens_function_sites = syst0.sites()
mount_vlead(syst0, greens_function_sites, 1)
sys0 = syst0.finalized()
coord_inject=np.array([sys0.pos(i) for i in sys0.lead_interfaces[0]]) 
coord=np.array([sys0.pos(i) for i in sys0.lead_interfaces[2]])  # sequence of virtual lead


def gf_virtual(cishu):
    salt=cishu+random.random()*100
    syst=make_system(width, length, str(salt))
    attach_lead(syst)
    greens_function_sites = syst0.sites()
    mount_vlead(syst, greens_function_sites, 1)
    sys = syst.finalized()
    #kwant.plot(sys, fig_size=(10, 3))
    
#    gf=kwant.greens_function(sys,0.4)
#    gf20=gf.submatrix(2,0)[:,0]
    gf20=kwant.greens_function(sys,0.4).submatrix(2,0)[:,0]
#    gf2=gf.submatrix(2,2)
#    gf22=gf2[:,3406]
#    ld=np.diag(gf.submatrix(2,2))
    myDict = {'coord':coord, 'coord_inject':coord_inject, 'gf20':gf20} #,'ld':ld
    completeName = os.path.join('E:/dwell3/53/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 


#gf_virtual(0)
Parallel(n_jobs=8)(delayed(gf_virtual)(n) for n in np.arange(405,500))

#syst=make_system(width, length, '0')
#attach_lead(syst)
#t=kwant.smatrix(syst.finalized(),.35).transmission(1,0)


elapsed=time()-t_ini
