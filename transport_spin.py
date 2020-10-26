# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:10:22 2017

@author: Yuhao Kang  %kymwswj   QSH
"""

import numpy as np
from matplotlib import pyplot
import random
import kwant
from kwant.digest import uniform

a=1
salt='wdlxwws'
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices

sys=kwant.Builder()

s0=np.identity(2)
sx=np.array([[0,1],[1,0]])
#sy=np.array([[0,-1j],[1j,0]])
sz=np.array([[1,0],[0,-1]])

#m=np.array([-6.1360,   -4.2952,         0,    4.2952,    6.1360])
#v=np.array([3.869,   -2.247,   -5.3252,   -2.247,    3.869])
ho=np.array([.4, -.4])

def disk(pos):
    x,y=pos
    return abs(y)<30 and np.sqrt(3)*x-100<y<np.sqrt(3)*x+100

def edge(pos):
    x,y=pos
    return abs(y)<30 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

def onsite1(site):
    x,y=site.pos
    t=[0 if y>=0 else 1]
    return ho[t]*sz+0.5*sx#+.5*uniform(repr(site),salt)*sx
#    if x**2+y**2<50:
#        return ho[t]*sz+1*random.random()*sx
#    else:
#        return ho[t]*sz+0.1*sx#m[t]*sz+v[t]*sx

def onsite2(site):
    x,y=site.pos
    t=[0 if y>=0 else 1]
    return -ho[t]*sz+0.5*sx
#    return -ho[t]*sz+.5*uniform(repr(site),salt)*sx
#    return -random.random()*ho[t]*sz+.5*random.random()*sx
#    if x**2+y**2<50:
#        return -ho[t]*sz+1*random.random()*sx
#    else:
#        return -ho[t]*sz+0.1*sx#-m[t]*sz+v[t]*sx

sys[A.shape(disk,(0,0))]=onsite1
sys[B.shape(disk,(0,0))]=onsite2
sys[lat.neighbors()]=-s0*2/3
#kwant.plot(sys)
lead=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
lead[A.shape(edge,(0,0))]=onsite1
lead[B.shape(edge,(0,0))]=onsite2
lead[lat.neighbors()]=-s0*2/3
     
#kwant.plot(lead)
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())
pyplot.rcParams["figure.figsize"] = [20,10]
#kwant.plot(sys)

sys=sys.finalized()
en=0.35

#sca=kwant.smatrix(sys,en).data
#print(sca)

wf=kwant.wave_function(sys,en)(0)
wavef1=[]
for i in range(wf.shape[1]//2):
    wavef1.append(wf[0,2*i+1])
kwant.plotter.map(sys,(abs(np.array(wavef1))**2))
#wavef2=[]
#for i in range(wf.shape[1]//2):
#    wavef2.append(wf[1,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wavef2))**2))

gf_mode=kwant.smatrix(sys,en).submatrix(1,0)      
u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
print(s)
#
#print(kwant.smatrix(sys,en).transmission(0,1))
#print(kwant.smatrix(sys,en).num_propagating(0))
##
#gf=kwant.greens_function(sys,en).submatrix(1,0)
#print(kwant.greens_function(sys,en).transmission(1,0))
#
#u1, s1, v1 = np.linalg.svd(gf, full_matrices=True)
#print(s1)

#mm=kwant.physics.modes()
#mn=kwant.physics.PropagatingModes()
