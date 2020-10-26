# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 17:20:10 2018

@author: user
"""
# https://kwant-project.org/doc/1.0/tutorial/tutorial5#id1
# set spin up/down as two layers
import kwant
import numpy as np
from matplotlib import pyplot
from kwant.digest import uniform
import scipy.io as sio
import os.path

#salt='ed'
m1 = .4  #valley
m2 = .1
Lambda_R=2

bo=0  #boundary between qvh and qsh

width=15
length=10
delta=0#2*width+5

lat_u = kwant.lattice.honeycomb(a=1, name='u')  #spin up
a, b = lat_u.sublattices
nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

lat_d = kwant.lattice.honeycomb(a=1, name='d')  #spin down
ad, bd = lat_d.sublattices
nnn_hoppings_ad = (((-1, 0), ad, ad), ((0, 1), ad, ad), ((1, -1), ad, ad))
nnn_hoppings_bd = (((1, 0), bd, bd), ((0, -1), bd, bd), ((-1, 1), bd, bd))
nnn_hoppingsd = nnn_hoppings_ad + nnn_hoppings_bd

def disk(pos):
    x,y=pos
    return abs(y)<width and abs(x)<length   #25.1

def lead0_shape(pos):
    x,y=pos
    return abs(y)<width

def onsite(site):
    x,y=site.pos
    if y>bo:
        return 0   
#        return -m1 * (1 if site.family == a else -1)
    else:
        return m1 * (1 if site.family == a or site.family ==ad else -1)
    
def nnn(site1, site2):
    x,y=site1.pos
    if y<=bo:
        return 0
#    elif x**2+(y-3)**2<5:
#        return -1j *m2
    else: 
        return 1j *m2*(1 if site1.family == a or site1.family ==b else -1)

def hop_ras_e1(site1,site2,Lambda_R=Lambda_R):   #rashba
    x,y=site1.pos
    if x**2+(y-3)**2<25 and y>0:
        return 1j*Lambda_R#* 2*(uniform(repr(site1),salt)-0.5)
    else:
        return 0
def hop_ras_e2(site1,site2,Lambda_R=Lambda_R):
    x,y=site1.pos
    if x**2+(y-3)**2<25 and y>0:
        return 1j*Lambda_R*(0.5 + 1j*np.sqrt(3)/2.0 \
               *(1 if site1.family == a or site1.family ==b else -1))#* 2*(uniform(repr(site1),salt)-0.5)
    else:
        return 0
def hop_ras_e3(site1,site2,Lambda_R=Lambda_R):
    x,y=site1.pos
    if x**2+(y-3)**2<25 and y>0:
        return 1j*Lambda_R*(0.5 - 1j*np.sqrt(3)/2.0\
               *(1 if site1.family == a or site1.family ==b else -1))#* 2*(uniform(repr(site1),salt)-0.5)
    else:
        return 0
        
def make_system():
    sys = kwant.Builder()

    sys[lat_u.shape(disk,(0,0))] = onsite
    sys[lat_d.shape(disk,(0,0))] = onsite

    sys[lat_u.neighbors()] = 1
    sys[lat_d.neighbors()] = 1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn
    
    sys[kwant.HoppingKind((0, 0), a,bd)] = hop_ras_e1
    sys[kwant.HoppingKind((0, 1), a,bd)] = hop_ras_e2
    sys[kwant.HoppingKind((-1, 1), a,bd)] = hop_ras_e3
    
    sys[kwant.HoppingKind((0, 0), ad,b)] = hop_ras_e1
    sys[kwant.HoppingKind((0, 1), ad,b)] = hop_ras_e2
    sys[kwant.HoppingKind((-1, 1), ad,b)] = hop_ras_e3
    
#    sym_left = kwant.TranslationalSymmetry((-1, 0))
#    lead0 = kwant.Builder(sym_left)
#    lead0[lat_u.shape((lambda pos: abs(pos[1]) < width), (0, 0))]=onsite#1
#    lead0[lat_u.neighbors()] = 1
#    lead0[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn
##    # left hole lead
#    lead1 = kwant.Builder(sym_left)
#    lead1[lat_d.shape((lambda pos: abs(pos[1]) < width), (0, 0))]=onsite#1
#    lead1[lat_d.neighbors()] = 1
#    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn  
#    
#    sys.attach_lead(lead0)   #spin up
#    sys.attach_lead(lead1)   #spin down
#    sys.attach_lead(lead0.reversed())   #lead3 spin up
#    sys.attach_lead(lead1.reversed())   #lead4 spin down
    
    sym0u = kwant.TranslationalSymmetry((-1,0))
    sym0u.add_site_family(lat_u.sublattices[0], other_vectors=[(-1, 2)])
    sym0u.add_site_family(lat_u.sublattices[1], other_vectors=[(-1, 2)])
    lead0_u = kwant.Builder(sym0u)    
    lead0_u[lat_u.shape(lead0_shape,(0,1))]= onsite #  1#  
    lead0_u[lat_u.neighbors()] = 1   
    lead0_u[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn    
    sys.attach_lead(lead0_u)
    #================================left_down0=====================================   
    sym0d = kwant.TranslationalSymmetry((-1,0))
    sym0d.add_site_family(lat_d.sublattices[0], other_vectors=[(-1, 2)])
    sym0d.add_site_family(lat_d.sublattices[1], other_vectors=[(-1, 2)])
    lead0_d = kwant.Builder(sym0d)    
    lead0_d[lat_d.shape(lead0_shape,(0,1))]= onsite #    1#
    lead0_d[lat_d.neighbors()] = 1    
    lead0_d[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn    
    sys.attach_lead(lead0_d)
    sys.attach_lead(lead0_u.reversed())
    sys.attach_lead(lead0_d.reversed())
    return sys

sys=make_system().finalized()

pyplot.rcParams["figure.figsize"] = [20,15]
kwant.plot(sys)

#gf_mode=kwant.smatrix(sys,0).submatrix(2,1) 
tm=kwant.smatrix(sys,0)
t01=tm.transmission(0,1)
t11=tm.transmission(1,1) 
t21=tm.transmission(2,1)
t31=tm.transmission(3,1)   

#wff=kwant.wave_function(sys,0)
#wf=wff(1)
#
#wf1=wf[0]
#ans1=[]
#for k in range(len(wf1)):
#    ans1.append(np.append(sys.pos(k),wf1[k]))
#myDict = {'ans1':ans1}
#completeName = os.path.join('E:/temp/', "9.mat")
#sio.savemat(completeName,myDict,oned_as='row') 

#sys1 = kwant.Builder()
#sys1[lat_u.shape(disk,(0,0))] = 0
#sys1[lat_u.neighbors()] = 0
#sym_left = kwant.TranslationalSymmetry((-1, 0))
#leadn = kwant.Builder(sym_left)
#leadn[lat_u.shape((lambda pos: abs(pos[1]) < width), (0, 0))]=0
#leadn[lat_u.neighbors()] = 0
#sys1.attach_lead(leadn) 
#sys1.attach_lead(leadn.reversed()) 
#sys1=sys1.finalized()
#kwant.plot(sys1)

#wf_u=[wf[0,i] for i in np.arange(0,wf.shape[1]//2)]
#kwant.plotter.map(sys1,(np.abs(np.array(wf_u))**2), fig_size=(10, 5)) 
## 
#wf_d=[wf[0,i] for i in np.arange(wf.shape[1]//2,wf.shape[1])]
#kwant.plotter.map(sys1,(np.abs(np.array(wf_d))**2), fig_size=(10, 5)) 
#wf_d2=[wf[0,i+wf.shape[1]//2] for i in np.arange(0,wf.shape[1]//2)]

#def plot_conductance(sys, energies):
#    # Compute conductance
#    data = []
#    for energy in energies:
#        smatrix = kwant.smatrix(sys, energy)
#        # Conductance is N - R_ee + R_he
#        data.append(smatrix.submatrix(0, 0).shape[0] -
#                    smatrix.transmission(0, 0) +
#                    smatrix.transmission(1, 0))
#    pyplot.plot(data)
#sys=make_system().finalized()       
#plot_conductance(sys, np.linspace(0,0.2,10))

#pyplot.plot(np.abs(wf[0]))

#lnt=np.sum(np.log(np.abs(wf)**2),axis=0)

#kwant.plotter.map(sys, np.abs(wf[0]),num_lead_cells=5,fig_size=(20, 15))
