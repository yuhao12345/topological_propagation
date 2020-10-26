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
m, t=0.4, -0.133
dis=0
#en=0.3
#salt='wdlxwws'
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
    return abs(y)<30 and np.sqrt(3)*x-69<y<np.sqrt(3)*x+9

def edge(pos):
    x,y=pos
    return abs(y)<10 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 
def edge1(pos):
    x,y=pos
    return abs(y)<10 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

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
#    global rl
#    t=random()
#    rl.append(t)
    if (x-25*a)**2+(y-3*a)**2<.5:  #(x-25*a)**2+(y-9*a)**2<19:(x-25*a)**2+(y-5*a)**2<2
        return -m*sz+0*sx
    elif y>-0.1*a:
        return m*sz+0*sx
    else:
        return np.zeros([2,2])
    
def onsite2(site):
#    global mn
#    mn+=1
    x,y=site.pos
    if (x-25*a)**2+(y-3*a)**2<.5:
        return m*sz+0*sx
    elif y>-0.1*a:
        return -m*sz+0*sx
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

lead1=kwant.Builder(kwant.TranslationalSymmetry((a,0))) 
lead1[A.shape(edge1,(0,0))]=onsite1a
lead1[B.shape(edge1,(0,0))]=onsite2a
lead1[lat.neighbors()]=-s0*2/3
lead1[A.neighbors()]=hopp
lead1[B.neighbors()]=hopp2

sys.attach_lead(lead1)
#sys.attach_lead(lead0.reversed())

#pyplot.rcParams["figure.figsize"] = [20,15]
#kwant.plot(sys, fig_size=(30, 12))
#pyplot.savefig('testplot.png')

sys=sys.finalized()

#sca=kwant.smatrix(sys,en).data

#wf=kwant.wave_function(sys,.1196)(0)
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i])
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(30, 12))

#wavef2=[]
#for i in range(wf.shape[1]//2):
#    wavef2.append(wf[1,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wavef2))**2))
#s_list=[]
#for conf in range(1):
#    rl=[]
#    mn=0
#gf_mode=kwant.smatrix(sys,1).submatrix(1,0)      
#u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#    print(s)
#    s_list.append(s)
#kwant.greens_function(sys,en).transmission(1,0)  

energy=np.linspace(0.11,.16,400)  
spec=[kwant.greens_function(sys,en).submatrix(1,0)[2,2] for en in energy]  #2,2
ang=np.unwrap(np.angle(spec))
#trans=[abs(sp)**2 for sp in spec]
pyplot.plot(energy,ang,'o')

#pyplot.figure()
#pyplot.plot(energy,trans,'o')
#temp=kwant.greens_function(sys,.3).out_block_coords(1)
#print(kwant.greens_function(sys,.3).num_propagating(0))

#condu=[kwant.greens_function(sys,en).transmission(1,0) for en in energy]

condu=[kwant.smatrix(sys,en).transmission(1,0) for en in energy]
pyplot.figure()
pyplot.plot(energy,condu,'o')

#matr=kwant.greens_function(sys,0).submatrix(1,0)
#print(matr.max())
#print(np.unravel_index(np.argmax(matr, axis=None), matr.shape))
#import scipy.io
#x = s_list
#matfile = 's_list3.mat'
#scipy.io.savemat(matfile, mdict={'tau': x}, oned_as='row')
#matdata = scipy.io.loadmat(matfile)
#assert np.all(x == matdata['tau'])

elapsed=time()-t_ini