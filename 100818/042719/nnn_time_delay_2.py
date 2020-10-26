# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 22:50:11 2019

@author: ykang
"""
##### many configuration, time delay
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

width=60
length=50

dis=1.6
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
        return 0<y<width and abs(x)<length   #25.1

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*dis

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return 0<y<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

e1 = 0.39 #0.39515
df=1e-9


#gf_mode=kwant.smatrix(sys,e1).submatrix(1,0)
#n=gf_mode.shape[0]

def eigentime(u0,u1,vh0,vh1):

    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t

def t_gf(sys):
    try:
        gf=kwant.greens_function(sys,e1).submatrix(1,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=kwant.greens_function(sys,e1+df).submatrix(1,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
        
        t=np.diag(eigentime(u0,u1,vh0,vh1))[0]  # time ,not dos
    except:
        t=0  
    return t

def dwell(cishu):
    salt=cishu#+random.random()*100
    
#    width=49.6
    syst=make_system(width, length, str(salt))
    attach_lead(syst)
#    kwant.plot(syst,fig_size=(25, 30))
    sys = syst.finalized()
    tt = t_gf(sys)
    wf=kwant.wave_function(sys,e1)(0)
    t_wf=sum(abs(wf[0])**2)/2/np.pi
    return [tt,t_wf]

#t_list=Parallel(n_jobs=5)(delayed(dwell)(cishu) for cishu in np.arange(0,2500,1))
#myDict = {'t_list':t_list}  #'s':s, 
#completeName = os.path.join('E:/dwell3/827/',str(width)+'w3.mat')
#sio.savemat(completeName,myDict,oned_as='row')    


def make_lead():
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width 
    sym = kwant.TranslationalSymmetry((1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
    return lead.finalized()
#flead=make_lead()
#
#prop_modes, _ = flead.modes(energy=e1)
#velo_lead=prop_modes.velocities[1]*2*np.pi
#t_ballistic=(2*length-1)/velo_lead

syst=make_system(width, length, str(1))
attach_lead(syst)
sys = syst.finalized()
wf=kwant.wave_function(sys,e1)(0)   #.403586
kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)


elapsed=time()-t_ini