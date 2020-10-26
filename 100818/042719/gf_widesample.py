# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 16:40:50 2019

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
#        return -50<y<width and abs(x)<length   #25.1
        return -50<y<width and length<=x<3*length

    def onsite(site):
        x,y=site.pos
        if y>-30:
            return (uniform(repr(site),salt)-0.5)*dis
        else:
            return 0

    sys=kwant.Builder()
    sys[graphene.shape(disk,(400,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return -50<y<width #abs(y)<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

syst=make_system(width, length, '2')
attach_lead(syst)
sys = syst.finalized()

ens = np.linspace(0.3924,0.3978,8192)
#ens = np.linspace(0.3924,0.3978,1000)
#ens=np.linspace(0.3870,0.3924,8192)
df=1e-9


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

def dwell(cishu):
    e1=ens[cishu]

    tt = t_gf(sys,e1)
    return tt

#t_list=Parallel(n_jobs=10)(delayed(dwell)(cishu) for cishu in np.arange(0,1000,1))
    
def wf_integral(cishu):
    e1=ens[cishu]
    wf=kwant.wave_function(sys,e1)(0)[0]
    return sum(abs(wf)**2)

#wf=Parallel(n_jobs=10)(delayed(wf_integral)(cishu) for cishu in np.arange(0,1000,1))
def T(cishu):
    e1=ens[cishu]
    g=kwant.smatrix(sys,e1)
    t=g.transmission(1,0)
    r=g.transmission(0,0)
    return [t,r]
#condu=Parallel(n_jobs=10)(delayed(T)(cishu) for cishu in np.arange(0,1000,1))
    
#myDict = {'t_list':t_list, 'ens' :ens }  #'s':s, 
#completeName = os.path.join('E:/dwell3/807/2.mat')
#sio.savemat(completeName,myDict,oned_as='row')    
    
def gf_01(cishu):
    en=ens[cishu]
    gf=kwant.greens_function(sys,en).submatrix(1,0)[:,0]
#    gf=kwant.smatrix(sys,en).submatrix(1,0)[0,0]
    myDict = {'gf':gf} #,'ld':ld
    completeName = os.path.join('E:/dwell3/810/', str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#    return gf
#Parallel(n_jobs=8)(delayed(gf_01)(cishu) for cishu in np.arange(0,8192,1))
    
#def make_lead():
#    def lead_shape(pos):
#        x,y=pos
#        return abs(y)<width 
#    sym = kwant.TranslationalSymmetry((1,0))
#    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
#    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
#    
#    lead = kwant.Builder(sym)
#    
#    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
#    lead[graphene.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
#    return lead.finalized()
#flead=make_lead()
#
#
#prop_modes, _ = flead.modes(energy=e1)
#velo_lead=prop_modes.velocities[1]*2*np.pi
#t_ballistic=(2*length-1)/velo_lead

en=0.3916 #ens[1255] #energies[522]
#wf=kwant.wave_function(sys,en)(0)   #.403586
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10)) #,colorbar=False
gf10=kwant.greens_function(sys,en).submatrix(1,0)[0,0]
gf01=kwant.greens_function(sys,en).submatrix(0,1)[0,0]
#coord=np.array([sys.pos(i) for i in range(73600)])
#myDict = {'wf':wf,'coord':coord} #,'ld':ld
#completeName = os.path.join('E:/dwell3/807/tt.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini