# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 12:58:35 2019

@author: ykang
"""
#### spectrum of gf
#### narrow sample
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
import dwell_gf_summary as dg
t_ini=time()

width=60
length=50


dis= 1.6
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
        return 0<y<width and abs(x)<length #abs(x)<length   #25.1 #

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*dis


    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,30))]= onsite  #0 #
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
    

#en=0.4
syst=make_system(width, length, '1')   # whole system as virtual lead
attach_lead(syst)
sys = syst.finalized()
#kwant.plot(sys)
### step 1, gf spectrum

#energies=np.linspace(0.3924,0.3978,8192)
#energies=np.linspace(0.35,0.4,500)
#energies=np.linspace(0.389,0.3924,5158)
def gf_01(cishu):
    en=energies[cishu]
    gf=kwant.greens_function(sys,en).submatrix(1,0)[:,0]
#    gf=kwant.smatrix(sys,en).submatrix(1,0)[0,0]
    myDict = {'gf':gf} #,'ld':ld
    completeName = os.path.join('E:/dwell3/824/', str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#    return gf
#Parallel(n_jobs=8)(delayed(gf_01)(cishu) for cishu in np.arange(0,5158,1))

#myDict = {'en':energies} #,'ld':ld
#completeName = os.path.join('E:/dwell3/792/en.mat')
#sio.savemat(completeName,myDict,oned_as='row') 
    


### step2  wf

#en=0.3906 #energies[475]
#wf=kwant.wave_function(sys,en)(1)   #.403586
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
#coord=np.array([sys.pos(i) for i in range(55338)])
#myDict = {'coord':coord,'wf':wf}
#completeName = os.path.join('E:/dwell3/wf/sample2_0.39403.mat')
#sio.savemat(completeName,myDict,oned_as='row') 



gf=kwant.smatrix(sys,.4).submatrix(1,0)
s0=np.linalg.svd(gf, full_matrices=True, compute_uv=False)
gf01=kwant.smatrix(sys,.4).submatrix(0,1)
s1=np.linalg.svd(gf01, full_matrices=True, compute_uv=False)
t01=kwant.smatrix(sys,.4).transmission(0,1)
t10=kwant.smatrix(sys,.4).transmission(1,0)
#en_wf=np.linspace(energies[2300],energies[5300],512)
    
def wf_01(cishu):
    en=en_wf[cishu]
    wf=kwant.wave_function(sys,en)(0)[0]
#    gf=kwant.greens_function(sys,en).submatrix(1,0)[:,0]
    myDict = {'wf':wf}
    completeName = os.path.join('E:/dwell3/786/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
#    return wf

#Parallel(n_jobs=10)(delayed(wf_01)(cishu) for cishu in np.arange(0,512,1))

ens=np.linspace(0.389,0.392,500)
df=1e-9
#def eigentime(u0,u1,vh0,vh1):
#
#    vv0=vh0.conj().T
#    vv1=vh1.conj().T
#    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
#    return t

#def t_gf(sys,e1):
#    try:
#        gf=kwant.greens_function(sys,e1).submatrix(1,0)
#        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
#        gf1=kwant.greens_function(sys,e1+df).submatrix(1,0)
#        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
#        
#        t=np.diag(eigentime(u0,u1,vh0,vh1))[0]  # time ,not dos
#    except:
#        t=0  
#    return t

def dwell(cishu):
    e1=ens[cishu]

    tt = dg.t_gf(sys,e1,df)
    wf=kwant.wave_function(sys,e1)(0)[0]
    wf_integral=sum(abs(wf)**2)/2/np.pi
    return [tt,wf_integral]

def T(cishu):
    e1=ens[cishu]
    return kwant.smatrix(sys,e1).transmission(1,0)
#t=Parallel(n_jobs=5)(delayed(T)(cishu) for cishu in np.arange(0,500,1))
#t_list=Parallel(n_jobs=5)(delayed(dwell)(cishu) for cishu in np.arange(0,1500,1))
#myDict = {'en':ens,'t':t} #,'ld':ld
#completeName = os.path.join('E:/dwell3/828/3.mat')
#sio.savemat(completeName,myDict,oned_as='row') 


elapsed=time()-t_ini

#sys.pos(sys.lead_interfaces[0][0])      #right lead first point

