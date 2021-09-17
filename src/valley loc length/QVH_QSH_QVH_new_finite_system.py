# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 19:23:53 2018

@author: user
"""

# emulate experiment, only consider disorder of spin-orbit coupling near boundary
# get 16*16 matrix of greens function 

import numpy as np
from matplotlib import pyplot
import scipy.io as sio
import os.path
#import random
import kwant
#from kwant.digest import uniform
from time import time

t_ini=time()
#dis = 0   # nnn disorder

en=0.25

l0=15.1  # system size in y direction
l=9  #boundary position
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
    return abs(y)<l0 and abs(x)<10

def onsite(site):
    x,y=site.pos
    if abs(y)<l:
#        return dis2*(uniform(repr(site),salt)-.5)*s0
        return np.zeros([2,2]) 
    else:
        return s0 * m1 * (1 if site.family == a else -1)

def spin_orbit_v(site1, site2):
    x,y=site1.pos
    if abs(y)<l:
        if x**2+(abs(y)-8)**2<25:       #point defect
            return -1j * m2 * sz
        else:
            return 1j * m2 * sz
    else:
        return np.zeros([2,2])
def spin_orbit_v_lead(site1, site2):
    x,y=site1.pos
    if abs(y)<l:
        return 1j * m2 * sz
    else:
        return np.zeros([2,2])
        
def onsite_lead(site):
    x,y=site.pos
    if abs(y)<l:
        return np.zeros([2,2])
    else:
        return s0 * m1 * (1 if site.family == a else -1)

sys=kwant.Builder()
sys[graphene.shape(disk,(0,0))]=onsite
sys[graphene.neighbors()]=-s0 
sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v

lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
lead[graphene.shape((lambda pos: abs(pos[1]) < l0), (0, 0))]=onsite_lead
lead[graphene.neighbors()]=-s0
lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v_lead
     

sys.attach_lead(lead)
sys.attach_lead(lead.reversed())

sys=sys.finalized()
#kwant.plot(sys,fig_size=(20, 10))


#
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

sca=kwant.smatrix(sys,en)
gf_mode=sca.submatrix(1,0) 
T4=sca.transmission(1,0)  
#wff=kwant.wave_function(sys,en)
#wf=wff(0)     #wave function of lead
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i])    # guess: 0 or i: positive or negative velocity; 2i or 2i+1: spin up or down
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(20, 10))        

gf=kwant.greens_function(sys, en).submatrix(1,0)

tau4=np.linalg.svd(gf,compute_uv=False)**2
#pyplot.plot(tau0,'.')
#
#sys.lead_interfaces
#sys.pos(1646)
##elapsed=time()-t_ini
#
#pyplot.imshow(np.log(abs(gf)))
