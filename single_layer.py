# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 18:41:37 2017

@author: Yuhao Kang  %kymwswj
"""

import numpy as np
from matplotlib import pyplot
import kwant

a=1
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices

sys=kwant.Builder()

#s0=np.identity(2)
#sx=np.array([[0,1],[1,0]])
#sy=np.array([[0,-1j],[1j,0]])
#sz=np.array([[1,0],[0,-1]])

#pot=np.array([[1,   -6.136,     3.869],
#    [1,   -4.2952,   -2.2475],
#    [1,    0,       -5.3252],
#    [1,    4.2952,   -2.2475],
#    [1,    6.136,     3.869]])
ho=np.array([.4, -.4])
def disk(pos):
    x,y=pos
    return abs(y)<25 and np.sqrt(3)*x-50<y<np.sqrt(3)*x+50

def edge(pos):
    x,y=pos
    return abs(y)<25 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

def onsite1(site):
    x,y=site.pos
    t=[0 if y>0 else 1]
    return ho[t]

def onsite2(site):
    x,y=site.pos
    t=[0 if y>0 else 1]
    return -ho[t]

sys[A.shape(disk,(0,1))]=onsite1
sys[B.shape(disk,(0,1))]=onsite2
#sys[lat.shape(disk,(0,1))]=onsite2
sys[lat.neighbors(1)]=-2/3
#kwant.plot(sys)
lead=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
lead[A.shape(edge,(0,0))]=onsite1
lead[B.shape(edge,(0,0))]=onsite2
lead[lat.neighbors(1)]=-2/3
     
#kwant.plot(lead)
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())
pyplot.rcParams["figure.figsize"] = [15,10]
kwant.plot(sys)


sys=sys.finalized()
en=0.25
sca=kwant.smatrix(sys,en).data


wf=kwant.wave_function(sys,en)
wavef=wf(0)
d=(abs(wavef)**2).sum(axis=0)
    
kwant.plotter.map(sys,d)
mode=sca.shape[1]//2
t=np.zeros([mode,mode],dtype=np.complex)
for i in range(mode):
    for j in range(mode):
        t[i,j]=sca[i+mode,j]
        
#t2=np.dot(t,np.matrix.getH(t))
#w,v= np.linalg.eig(t2)
u, s, v = np.linalg.svd(t, full_matrices=True)
print(s)
print(kwant.smatrix(sys,en).transmission(1,0))