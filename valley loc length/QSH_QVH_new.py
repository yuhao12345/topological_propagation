# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:40:31 2017

@author: user
"""

import numpy as np
from matplotlib import pyplot
import scipy.io as sio
import os.path
import random
import kwant
from kwant.digest import uniform
from time import time

t_ini=time()
#dis = 0   # nnn disorder

en=0.25

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m1 = .5  #valley
m2 = .077  #spin

s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b
def disk(pos):
    x,y=pos
    return abs(y)<15.1 and abs(x)<20

def onsite(site):
    x,y=site.pos
    if y>0:
#        return dis2*(uniform(repr(site),salt)-.5)*s0
        return np.zeros([2,2])     
    else:
        return s0 * m1 * (1 if site.family == a else -1)

def spin_orbit_v(site1, site2):
    x,y=site1.pos
    if y>0:
        return 1j * m2 * sz
    else:
        return np.zeros([2,2])
def spin_orbit_v_lead(site1, site2):
    x,y=site1.pos
    if y>0:
        return 1j * m2 * sz
    else:
        return np.zeros([2,2])
        
def onsite_lead(site):
    x,y=site.pos
    if y>0:
        return np.zeros([2,2])
    else:
        return s0 * m1 * (1 if site.family == a else -1)

sys=kwant.Builder()
sys[graphene.shape(disk,(0,0))]=onsite_lead
sys[graphene.neighbors()]=-s0 
sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v_lead

lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
lead[graphene.shape((lambda pos: abs(pos[1]) < 10), (0, 0))]=onsite_lead
lead[graphene.neighbors()]=-s0
lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v_lead
     
#kwant.plot(lead)
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())

sys=sys.finalized()

#kwant.plot(sys, fig_size=(30, 12))
#
#matr=kwant.greens_function(sys,0).submatrix(1,0)[67,67]

energy=np.linspace(0,.4,100)  
spec=[kwant.greens_function(sys,en).submatrix(1,0)[67,67] for en in energy]  #4,4 and 18,18 spin up/down
ang=np.angle(spec) #np.unwrap(np.angle(spec))
trans=[abs(sp)**2 for sp in spec]
pyplot.figure(figsize=(10,6))
pyplot.plot(energy,ang,'o')
pyplot.title('phase')

pyplot.figure()
pyplot.figure(figsize=(10,6))
pyplot.plot(energy,trans,'o')
pyplot.title('intensity')
#save_path = 'E:/TI transmission simulation/vsv12/'
#k=1
#for dis in np.arange(.25, .55,.05): 
#    s_list=[]
#    t_list=[]
#    
#    for cishu in np.arange(0,250):
#        salt=str(cishu+random.random())
#        try:
#            sca=kwant.smatrix(sys,en)
#            gf_mode=sca.submatrix(1,0)  
##            wff=kwant.wave_function(sys,en)
#            u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#            T=sca.transmission(1,0)  
#            s_list.append(s)
#            t_list.append(T)
#        except:
#            continue
#    myDict = {'s':s_list,'t':t_list}
#    completeName = os.path.join(save_path, str(k)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row')   
#    k=k+1
#
#wf=wff(0)
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i])
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(20, 10))        
#
elapsed=time()-t_ini

#temp=sys.lead_interfaces
#sys.pos(1871)   #34<0,35>0
