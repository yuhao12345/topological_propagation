# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:40:31 2017

@author: user
"""

import numpy as np
from matplotlib import pyplot
from random import random
import kwant
from kwant.digest import uniform
from time import time
#from kwant.digest import uniform
t_ini=time()
m, t=0.4, -0.15
dis=0.2
en=0.38
salt='wdlxwws'
a=1
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices

sys=kwant.Builder()

s0=np.identity(2)
sx=np.array([[0,1],[1,0]])
#sy=np.array([[0,-1j],[1j,0]])
sz=np.array([[1,0],[0,-1]])


def disk(pos):
    x,y=pos
    return abs(y)<26 and np.sqrt(3)*x-99<y<np.sqrt(3)*x+9

def edge(pos):
    x,y=pos
    return abs(y)<14 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

def hopp(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return -t*s0
    else:
        return np.zeros([2,2])
    
def hopp2(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return t*s0
    else:
        return np.zeros([2,2])
    

def onsite1a(site):
    x,y=site.pos
    if y>-0.1*a:
        return m*sz+0*sx
    else:
        return np.zeros([2,2])
    
def onsite2a(site):
    x,y=site.pos
    if y>-0.1*a:
        return -m*sz+0*sx
    else:
        return np.zeros([2,2])
def onsite1(site):
    x,y=site.pos
    global rl
    t=random()
    rl.append(t)
    if y>-0.1*a:
        return (m-(t-.5)*dis)*sz+0*sx
    else:
        return np.zeros([2,2])
    
def onsite2(site):
    global mn
    mn+=1
    x,y=site.pos
    if y>-0.1*a:
        return -(m-(rl[mn-1]-0.5)*dis)*sz+0*sx
    else:
        return np.zeros([2,2])
    
    
sys[A.shape(disk,(0,0))]=onsite1
sys[B.shape(disk,(0,0))]=onsite2
sys[lat.neighbors()]=-s0*2/3
sys[A.neighbors()]=hopp
sys[B.neighbors()]=hopp2

lead0=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
lead0[A.shape(edge,(0,0))]=onsite1a
lead0[B.shape(edge,(0,0))]=onsite2a
lead0[lat.neighbors()]=-s0*2/3
lead0[A.neighbors()]=hopp
lead0[B.neighbors()]=hopp2

sys.attach_lead(lead0)
sys.attach_lead(lead0.reversed())

#pyplot.rcParams["figure.figsize"] = [20,15]
#kwant.plot(sys)
#pyplot.savefig('testplot.png')

sys=sys.finalized()

#sca=kwant.smatrix(sys,en).data

#wf=kwant.wave_function(sys,en)(0)
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i])
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2))

#wavef2=[]
#for i in range(wf.shape[1]//2):
#    wavef2.append(wf[1,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wavef2))**2))
s_list=[]
for conf in range(2500):
    rl=[]
    mn=0
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0)      
    u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#    print(s)
    s_list.append(s)
    
import scipy.io
x = s_list
matfile = 's_list3.mat'
scipy.io.savemat(matfile, mdict={'tau': x}, oned_as='row')
matdata = scipy.io.loadmat(matfile)
assert np.all(x == matdata['tau'])

elapsed=time()-t_ini