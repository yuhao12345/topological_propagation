# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:22:55 2018

@author: user
"""

import numpy as np
import kwant
 
m, t=0.4, -0.133
dis=0
a=1
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices
sys=kwant.Builder()
s0=np.identity(2)
sx=np.array([[0,1],[1,0]])
sz=np.array([[1,0],[0,-1]])
 
def disk(pos):
    x,y=pos
    return abs(y)<2 and np.sqrt(3)*x-1<y<np.sqrt(3)*x+5
 
def edge(pos):
    x,y=pos
    return abs(y)<2 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58
def edge1(pos):
    x,y=pos
    return abs(y)<2 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58
 
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
   
 
def onsite(site):
    x,y=site.pos
    if (x-25*a)**2+(y-9*a)**2<8:
        return -m*sz+0*sx
    elif y>-0.1*a:
        return m*sz+0*sx
    else:
        return np.zeros([2,2])
   
    
sys[lat.shape(disk,(0,0))]=onsite
sys[lat.neighbors()]=-s0*2/3
sys[A.neighbors()]=hopp
sys[B.neighbors()]=hopp2
 
lead0=kwant.Builder(kwant.TranslationalSymmetry((-a,0)))
lead0[lat.shape(edge,(0,0))]=onsite
lead0[lat.neighbors()]=-s0*2/3
lead0[A.neighbors()]=hopp
lead0[B.neighbors()]=hopp2
sys.attach_lead(lead0)
sys.attach_lead(lead0.reversed())
 
kwant.plot(sys, fig_size=(12, 12))
 
sys=sys.finalized()
 
gf=kwant.greens_function(sys,0).submatrix(1,0)