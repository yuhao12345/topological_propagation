# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 19:23:53 2018

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

t_ini=time()
#dis = 0   # nnn disorder

population = [1, -1]
weights = [0, 2]

width=10
bo=0  #boundary between qvh and qsh

dis2=1
salt='ed'
en=0

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]],norbs=2)  # Coordinates of the sites
a, b = graphene.sublattices

m1 = .4  #valley
m2 = .1  #spin    3*sqrt(3)*m2
Lambda_R=2

s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b
def disk(pos):
    x,y=pos
    return abs(y)<width and abs(x)<10   #25.1

def onsite(site):
    x,y=site.pos

    if y>bo:
        return np.zeros([2,2])     
    else:
        return s0 * m1 * (1 if site.family == a else -1)

def spin_orbit_v(site1, site2):
    x,y=site1.pos

    if y<=bo:
        return np.zeros([2,2])
#    elif x**2+y**2<6 or (x+5)**2+(y-2)**2<6 or (x+8)**2+(y-3)**2<6 \
#    or (x+15)**2+(y-3)**2<6 or (x+22)**2+(y-4)**2<4 or (x+24)**2+(y-5)**2<4:
#        return -1j * m2 * sz*2
#    elif y>bo and y<5:
##        return 1j * m2 * choices(population, weights)[0] * sz
#        return 1j * m2 * 2*((uniform(repr(site1),salt)-0.5)) * sz
    else: 
        return 1j *m2 * sz
        
    
def spin_orbit_v_lead(site1, site2):
    x,y=site1.pos
#    return np.zeros([2,2])
    if y>bo:
        return 1j * m2 * sz
    else:
        return np.zeros([2,2])
        
def onsite_lead(site):
    x,y=site.pos
#    return s0
    if y>bo:
        return np.zeros([2,2])
    else:
        return s0 * m1 * (1 if site.family == a else -1)

def hop_ras_e1(site1,site2,Lambda_R=Lambda_R,t1=1):
    x,y=site1.pos
    if abs(x)<5 and 0<y<0:
        return s0-1j*Lambda_R*sx
    else:
        return s0
    

def hop_ras_e2(site1,site2,Lambda_R=Lambda_R,t1=1):
    x,y=site1.pos
    if abs(x)<5 and 0<y<0:
        return s0+1j*Lambda_R*(0.5*sx - np.sqrt(3)/2.0*sy)
    else:
        return s0
    

def hop_ras_e3(site1,site2,Lambda_R=Lambda_R,t1=1):
    x,y=site1.pos
    if abs(x)<5 and 0<y<0:
        return s0+1j*Lambda_R*(0.5*sx + np.sqrt(3)/2.0*sy)
    else:
        return s0
    

def make_system():
    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]=onsite
    sys[graphene.neighbors()]=-s0 
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v
    
    sys[kwant.HoppingKind((0, 0), a,b)] = hop_ras_e1
    sys[kwant.HoppingKind((0, 1), a,b)] = hop_ras_e2
    sys[kwant.HoppingKind((-1, 1), a,b)] = hop_ras_e3
    
    lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
    lead[graphene.shape((lambda pos: abs(pos[1]) < 10), (0, 0))]=onsite_lead
    lead[graphene.neighbors()]=-s0
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_v_lead
         
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
    return sys.finalized()

sys=make_system()

en=.4
wff=kwant.wave_function(sys,en)
wf=wff(0)
kwant.plotter.map(sys,abs(wf[0][1::2])**2, fig_size=(10, 10))

#def plot_conductance(sys, energies):
#    data = []
#    for energy in energies:
#        smatrix = kwant.smatrix(sys, energy)
#        data.append(smatrix.transmission(1, 0))
#
#    pyplot.figure()
#    pyplot.plot(energies, data)
#    pyplot.xlabel("energy [t]")
#    pyplot.ylabel("conductance [e^2/h]")
#    pyplot.show()
#energies=[i*0.02 for i in range(20)]
#plot_conductance(sys,energies)


#for cishu in range(1):
#    sys=make_system()
#    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
##    con=[]
##    enegies=np.linspace(0,0.4,51)
##    for en in enegies:
##        t=kwant.smatrix(sys,en).transmission(1,0)
##        con.append(t)
##    pyplot.plot(enegies,con,'.')
##    print(t)
#    wf=kwant.wave_function(sys,en)(0)
#    wavef1=[]    
#    for i in range(wf.shape[1]//2):
#        wavef1.append(wf[0,2*i+1])       
#    kwant.plotter.map(sys,((abs(np.array(wavef1))**2)), fig_size=(20, 5),num_lead_cells=5)
    

#    ## ln(I(x))~x
#    wf_down=[wf[:,i] for i in np.arange(1,wf.shape[1],2)]
#    wf1=[]
#    wf1=2*np.mean(np.log(np.abs(wf_down)),1)
#
#    ans=[]
#    for k in range(len(wf1)):
#        ans.append(np.append(sys.pos(k),wf1[k]))
#    
#    myDict = {'ans':ans}
#    completeName = os.path.join('E:/qv_qs_lnI/1/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row')   

#gf=kwant.greens_function(sys,en).submatrix(1,0)
#t1=kwant.greens_function(sys,en).transmission(1,0)
#sg=np.linalg.svd(gf, compute_uv=False)
#taug1=sg**2

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


