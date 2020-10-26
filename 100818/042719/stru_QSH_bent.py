# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 10:08:39 2020

@author: user
"""

import numpy as np
import kwant
from kwant.digest import uniform
#from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from matplotlib import pyplot

t_ini=time()
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b 

def make_system(width, length, salt, dis):
    def disk(pos):
        x,y=pos
        return -width<y<width+(0<x and x<30)*np.sqrt(3)*x+(x>=30)*np.sqrt(3)*30 and  abs(x)<length
#        return -width<y<width and  200<=x<600 #abs(x)<length  # #   #25.1 #0<y<width 

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*dis
#        if x>=200:
#            return (uniform(repr(x*y),salt)-0.5)*dis
#        else:
#            return (uniform(repr((x+400)*y),salt)-0.5)*dis

#        if y>-30:
#            return (uniform(repr(site),salt)-0.5)*dis
#        else:
#            return 0

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def make_lead(width):
    def lead_shape(pos):
        x,y=pos
        return -width<y<width #0<y<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, 0))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
    return lead

def make_lead1(width):
    def lead_shape(pos):
        x,y=pos
        return -width<y<np.sqrt(3)*30+width
    sym = kwant.TranslationalSymmetry((1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, 50))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
    return lead

def attach_lead(sys,width):
    lead=make_lead(width)
    lead1=make_lead1(width)
    sys.attach_lead(lead)
    sys.attach_lead(lead1)
    
def make_system_all(width, length, salt, dis):    # disordered QSH system
    syst=make_system(width, length, salt,dis) 
    attach_lead(syst,width)
    return syst.finalized()

length=100
dis=1.6
width=30
salt=0
sys=make_system_all(width, length, str(salt), dis)
#kwant.plot(sys)

wf=kwant.wave_function(sys,0)(0)  #,check_hermiticity=False
kwant.plotter.map(sys, (abs(wf[0]))**0.5,num_lead_cells=5,fig_size=(12, 10),colorbar=False)

ens=np.linspace(-0.5,0.5,256)
def spec(e1):
    return np.angle(kwant.greens_function(sys,e1).submatrix(1,0)[0,0])
#    return kwant.smatrix(sys,e1).submatrix(1,0)
#t=Parallel(n_jobs=8)(delayed(spec)(e1) for e1 in ens)
#tp=kwant.greens_function(sys,0).submatrix(1,0)[:,0]
#pyplot.plot(np.abs(tp))

#completeName0 = os.path.join('C:/Users/ykang/Dropbox/TI8/data/603.mat')
#sio.savemat(completeName0,{'t':t,'ens':ens},oned_as='row') 

elapsed=time()-t_ini