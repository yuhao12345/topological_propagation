# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 13:17:59 2017

@author: Yuhao Kang  
"""

import numpy as np
from matplotlib import pyplot
from random import random
import kwant
from kwant.digest import uniform

a=1
d,l=40,70 
salt='wdlxwws'
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices

sys=kwant.Builder()

s0=np.identity(2)
sx=np.array([[0,1],[1,0]])
#sy=np.array([[0,-1j],[1j,0]])
sz=np.array([[1,0],[0,-1]])
ho=np.array([.4, -.4])

def disk(pos):
    x,y=pos         
    return (-0.866*d<y<0.866*d and -l<x<=0 or 
            0<=y<=0.866*l and x-d<=y/np.sqrt(3)<x+d or 
            -0.866*l<=y<=0 and x-0.9*d<-y/np.sqrt(3)<x+d)

def edge0(pos):
    x,y=pos
    return -0.866*d<y<0.866*d and np.sqrt(3)*x-0.1<=y<=np.sqrt(3)*x+0.58 

def edge1(pos):
    x,y=pos
    return abs(x)<=d and -0.1*a<=y<=1.6*a

def hopp(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if ((y1-np.sqrt(3)*x1+0.57*a)*(y2-np.sqrt(3)*x2+0.57*a)>0 
        and y1-np.sqrt(3)*x1+0.57*a<0 
        and (y1+np.sqrt(3)*x1+0.57*a)*(y2+np.sqrt(3)*x2+0.57*a)>0 
        and y1+np.sqrt(3)*x1+0.57*a>0):
        return 0.26*(random()-0.5)*s0
    else:
        return np.zeros([2,2])
    
def hopp2(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if ((y1-np.sqrt(3)*x1+0.57*a)*(y2-np.sqrt(3)*x2+0.57*a)>0 
        and y1-np.sqrt(3)*x1+0.57*a<0 
        and (y1+np.sqrt(3)*x1+0.57*a)*(y2+np.sqrt(3)*x2+0.57*a)>0 
        and y1+np.sqrt(3)*x1+0.57*a>0):
        return -0.26*(random()-0.5)*s0
    else:
        return np.zeros([2,2])
def hoppa(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if ((y1-np.sqrt(3)*x1+0.57*a)*(y2-np.sqrt(3)*x2+0.57*a)>0 
        and y1-np.sqrt(3)*x1+0.57*a<0 
        and (y1+np.sqrt(3)*x1+0.57*a)*(y2+np.sqrt(3)*x2+0.57*a)>0 
        and y1+np.sqrt(3)*x1+0.57*a>0):
        return (0.13)*s0
    else:
        return np.zeros([2,2])
    
def hopp2a(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if ((y1-np.sqrt(3)*x1+0.57*a)*(y2-np.sqrt(3)*x2+0.57*a)>0 
        and y1-np.sqrt(3)*x1+0.57*a<0 
        and (y1+np.sqrt(3)*x1+0.57*a)*(y2+np.sqrt(3)*x2+0.57*a)>0 
        and y1+np.sqrt(3)*x1+0.57*a>0):
        return -(0.13)*s0
    else:
        return np.zeros([2,2])
def onsite1(site):
    x,y=site.pos
    if y>=0 and y-np.sqrt(3)*x+0.57*a>0:
        return ho[0]*sz+0*sx
    if y<0 and y+np.sqrt(3)*x+0.57*a<0:
        return ho[1]*sz+0*sx
    else:
        return np.zeros([2,2])
    
def onsite2(site):
    x,y=site.pos
    if y>=0 and y-np.sqrt(3)*x+0.57*a>0:
        return -ho[0]*sz+0*sx
    if y<0 and y+np.sqrt(3)*x+0.57*a<0:
        return -ho[1]*sz+0*sx
    else:
        return np.zeros([2,2])
    
sys[A.shape(disk,(0,0))]=onsite1
sys[B.shape(disk,(0,0))]=onsite2
sys[lat.neighbors()]=-s0*2/3
sys[A.neighbors()]=hopp
sys[B.neighbors()]=hopp2

lead0=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
lead0[A.shape(edge0,(0,0))]=onsite1
lead0[B.shape(edge0,(0,0))]=onsite2
lead0[lat.neighbors()]=-s0*2/3
lead0[A.neighbors()]=hoppa
lead0[B.neighbors()]=hopp2a
#kwant.plot(lead0)
lead1=kwant.Builder(kwant.TranslationalSymmetry((a/2,a/2*np.sqrt(3)))) 
lead1[A.shape(edge1,(0,0))]=onsite1
lead1[B.shape(edge1,(0,0))]=onsite2
lead1[lat.neighbors()]=-s0*2/3
lead1[A.neighbors()]=hoppa
lead1[B.neighbors()]=hopp2a
#kwant.plot(lead1)
lead2=kwant.Builder(kwant.TranslationalSymmetry((a/2,-a/2*np.sqrt(3)))) 
lead2[A.shape(edge0,(0,0))]=onsite1
lead2[B.shape(edge0,(0,0))]=onsite2
lead2[lat.neighbors()]=-s0*2/3
lead2[A.neighbors()]=hoppa
lead2[B.neighbors()]=hopp2a
#kwant.plot(lead2)

sys.attach_lead(lead0)
sys.attach_lead(lead1)
sys.attach_lead(lead2)
pyplot.rcParams["figure.figsize"] = [20,15]
kwant.plot(sys)
#pyplot.savefig('testplot.png')

#sys=sys.finalized()
#en=0.1
##
#sca=kwant.smatrix(sys,en).data
#
#wf=kwant.wave_function(sys,en)(0)
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i])
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2))
#
#wavef2=[]
#for i in range(wf.shape[1]//2):
#    wavef2.append(wf[1,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wavef2))**2))
#
#gf_mode=kwant.smatrix(sys,en).submatrix(1,0)      
#u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#print(s)
##
#print(kwant.smatrix(sys,en).transmission(1,0))
##print(kwant.smatrix(sys,en).num_propagating(0))
###
#gf=kwant.greens_function(sys,en).submatrix(1,0)
#print(kwant.greens_function(sys,en).transmission(1,0))
##
#u1, s1, v1 = np.linalg.svd(gf, full_matrices=True)
#print(s1)