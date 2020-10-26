# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 20:40:21 2019

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


length=50  # toi calculate loc length, sample length=1e4
width=15

dis=0

lat = kwant.lattice.square()

def make_system(length,width,salt):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width 
    def onsite(site):
        x,y=site.pos
        return dis*(uniform(repr(site),salt)-.5)+4

    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = -1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 4
    lead[lat.neighbors()] = -1

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

sys=make_system(length,width, '0')


def make_lead():
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width 
    sym = kwant.TranslationalSymmetry((1,0))

    
    lead = kwant.Builder(sym)
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 4
    lead[lat.neighbors()] = -1
    return lead.finalized()
flead=make_lead()

e1 = .7
df=1e-8

#c=np.array([flead.pos(n) for n in range(29)])  # lead sequence
#c2=np.array([sys.pos(sys.lead_interfaces[0][n]) for n in range(29)]) # real system
gf_mode=kwant.smatrix(sys,e1).submatrix(1,0)
n=gf_mode.shape[0]


def s(en):
    prop_modes, _ = flead.modes(energy=en)
    velo_lead=prop_modes.velocities[n:2*n]
    
    wf_lead=np.abs(prop_modes.wave_functions[:,n:2*n]) @ np.diag(velo_lead**0.5)
    
    g=kwant.greens_function(sys,en)
    g10=g.submatrix(1,0)
    g00=g.submatrix(0,0)
    g11=g.submatrix(1,1)
    g01=g.submatrix(0,1)
    t10= 1j*(wf_lead.conj().T @ g10 @ wf_lead)
    t01= 1j*(wf_lead.conj().T @ g01 @ wf_lead)
    t00= -np.eye(n) + 1j*(wf_lead.conj().T @ g00 @ wf_lead)
    t11= -np.eye(n) + 1j*(wf_lead.conj().T @ g11 @ wf_lead)
    s1=np.bmat([[t00, t01], [t10, t11]])
    return s1

si=s(e1)
sf=s(e1+df)

uni=si.conj().T@si
dos=np.trace(si.conj().T@(sf-si))/df/(4*np.pi**2)*(-1j)

def s_sys(en):
    s=kwant.smatrix(sys,en).data
    return s
si0=s_sys(e1)
sf0=s_sys(e1+df)
dos_sys=np.trace(si0.conj().T@(sf0-si0))/df/(4*np.pi**2)*(-1j)

def eigentime(u0,u1,vh0,vh1):

    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t

gf=kwant.greens_function(sys,e1).submatrix(1,0)
u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
gf1=kwant.greens_function(sys,e1+df).submatrix(1,0)
u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
    
#u0, s0, vh0=np.linalg.svd(si[n:2*n,0:n], full_matrices=True, compute_uv=True)
#u1, s1, vh1=np.linalg.svd(sf[n:2*n,0:n], full_matrices=True, compute_uv=True)

t=np.diag(eigentime(u0,u1,vh0,vh1)/np.pi)
tt=sum(t[np.abs(np.imag(t))<1e-4])
elapsed=time()-t_ini