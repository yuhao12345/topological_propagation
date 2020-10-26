# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 17:19:49 2019

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
        return -width<y<width and abs(x)<length   #25.1

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
        return -width<y<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, 0))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

syst=make_system(width, length, '2')   # whole system as virtual lead
attach_lead(syst)
sys = syst.finalized()

e1 = .39
df=1e-9

#def s_sys(en):
#    s=kwant.smatrix(sys,en).data
#    return s

#si0=s_sys(e1)
#sf0=s_sys(e1+df)
#dos_sys=np.trace(si0.conj().T@(sf0-si0))/df/(4*np.pi**2)*(-1j)

#gf_mode=kwant.smatrix(sys,e1).submatrix(1,0)
#n=gf_mode.shape[0]

def eigentime(u0,u1,vh0,vh1):

    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t
def t_gf(sys,e1):
    try:
        gf=kwant.greens_function(sys,e1).submatrix(1,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=kwant.greens_function(sys,e1+df).submatrix(1,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
        
        t=np.diag(eigentime(u0,u1,vh0,vh1))[0:1]  # time ,not dos
    except:
        t=0  
    return t
def r_gf(sys,e1):
    try:
        gf=kwant.greens_function(sys,e1).submatrix(0,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=kwant.greens_function(sys,e1+df).submatrix(0,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
        
        t=np.diag(eigentime(u0,u1,vh0,vh1))#[0:1]  # time ,not dos
    except:
        t=0  
    return t
### from gf
#gf=kwant.greens_function(sys,e1).submatrix(1,0)[:,0]
#u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
#gf1=kwant.greens_function(sys,e1+df).submatrix(1,0)[:,0]
#u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)

### from gf_mode
#gf_mode=kwant.smatrix(sys,e1).submatrix(1,0)
#u0, s0, vh0=np.linalg.svd(gf_mode, full_matrices=True, compute_uv=True)
#gf1_mode=kwant.smatrix(sys,e1+df).submatrix(1,0)
#u1, s1, vh1=np.linalg.svd(gf1_mode, full_matrices=True, compute_uv=True)

#t=np.diag(eigentime(u0,u1,vh0,vh1))#[0:n]  # time ,not dos
#tt=sum(t[np.abs(np.imag(t))<1e-4])
#tt=sum(t)

#t=t_gf(sys,e1)
#t_r=r_gf(sys,e1)
wf=kwant.wave_function(sys,0.389)(0)
#t_wf=sum(abs(wf[0])**2)/2/np.pi
kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10))

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
#
#prop_modes, _ = flead.modes(energy=e1)
#velo_lead=prop_modes.velocities[1]*2*np.pi
#t_ballistic=(2*length-1)/velo_lead

#wf_lead=prop_modes.wave_functions
#
#c=np.array([flead.pos(n) for n in range(138)])  # lead sequence
#c2=np.array([sys.pos(sys.lead_interfaces[0][n]) for n in range(138)]) # real system
#
#keys=c[:,1]
#v0=wf_lead[:,0]
#v1=wf_lead[:,1]
#d0 = dict(zip(keys, v0))
#d1 = dict(zip(keys, v1))
#wf_l=np.array([d0[n]  for n in c2[:,1]])  # sort based on system sequence
#wf_r=np.array([d1[n]  for n in c2[:,1]])
#
##pyplot.plot(c2[:,1],abs(wf_l),'.')
##pyplot.figure()
##pyplot.plot(c2[:,1],abs(wf_r),'o')
#
#g=kwant.greens_function(sys,en)
#g10=g.submatrix(1,0)
#g00=g.submatrix(0,0)
#g11=g.submatrix(1,1)
#g01=g.submatrix(0,1)
#t10=1j*velo_lead[1]*(wf_r.conj().T @ g10 @ wf_r)
#t01=1j*velo_lead[1]*(wf_l.conj().T @ g01 @ wf_l)
#t00= -1 + 1j*velo_lead[1]*(wf_l.conj().T @ g00 @ wf_r)
#t11= -1 + 1j*velo_lead[1]*(wf_r.conj().T @ g11 @ wf_l)
#s1=np.array([[t00[0,0],t01[0,0]],[t10[0,0],t11[0,0]]])
##s1=1/(s0.conj().T@s0)[0,0]**0.5*s0   # make sure unitary
#    
#s=kwant.smatrix(sys,en).data
##np.trace(s1.conj().T@(s12-s1))/df/(4*np.pi**2)*(-1j)
##np.trace(s.conj().T@(s2-s))/df/(4*np.pi**2)*(-1j)
#s1.conj().T@s1

elapsed=time()-t_ini