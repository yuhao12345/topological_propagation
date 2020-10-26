# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 21:22:35 2017

@author: Yuhao Kang    %kymwswj QVH  triangular lattice 
"""

import numpy as np
from matplotlib import pyplot
from random import random
import kwant
from time import time
#from kwant.digest import uniform
t_ini=time()
s_list=[]
en=0
a=1
w=0.25
# define lattice and shape
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices
#s0=np.identity(2)
def disk(pos):
    x,y=pos
    return abs(y)<30 and np.sqrt(3)*x-150<y<np.sqrt(3)*x+100

def edge(pos):
    x,y=pos
    return abs(y)<28 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

#define self and hopping energy

def hopp(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    global rl,rl2
    
#    mm+=1
    t=random()
    t2=random()
    rl.append(t)
    rl2.append(t2)
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a>0:
        return (0.15-w*t)#*s0
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return -(0.15-w*t2)#*s0
    else:
        return 0#np.zeros([2,2])
    
    
    
def hopp2(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    global mn
    mn+=1
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a>0:
        return -(0.15-w*rl[mn-1])#*s0
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return (0.15-w*rl2[mn-1])#*s0
    else:
        return 0#np.zeros([2,2])
    
def hopp_0(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a>0:
        return 0.15#*s0
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return -0.15#*s0
    else:
        return 0#np.zeros([2,2])
    
def hopp2_0(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a>0:
        return -0.15#*s0
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return 0.15#*s0
    else:
        return 0#np.zeros([2,2])
def make_system():  
    global rl,rl2,mn
    rl=[]
    rl2=[]    
#    mm=0
    mn=0 
    
    sys=kwant.Builder()    
    sys[lat.shape(disk,(0,0))]=0#np.zeros([2,2])
    sys[lat.neighbors()]=-2/3#*s0
    sys[A.neighbors()]=hopp
    sys[B.neighbors()]=hopp2
    lead=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
    lead[lat.shape(edge,(0,0))]=0#np.zeros([2,2])
    lead[lat.neighbors()]=-2/3#*s0
    lead[A.neighbors()]=hopp_0
    lead[B.neighbors()]=hopp2_0
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

#    pyplot.rcParams["figure.figsize"] = [20,10]
#kwant.plot(sys)

    return sys.finalized()   
sys=make_system()
for times in range(1):    
  
#    salt='wdlxwws'    
#    sca=kwant.smatrix(sys,en).data
#    sca2=kwant.smatrix(sys,en).data
    ##print(sca)
    wf=kwant.wave_function(sys,en)(0)
    wavef1=[]    
    for i in range(wf.shape[1]):
        wavef1.append(wf[1,i])       
    kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(30, 12))
    rl=[]
    rl2=[]    
    mn=0 
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0)      
    u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
    print(s)
    
    s_list.append(s)
#    print(mn)
#print(np.array(s_list))
#print(kwant.smatrix(sys,en).transmission(1,0))
##print(kwant.smatrix(sys,en).num_propagating(0))
###
#gf=kwant.greens_function(sys,en).submatrix(1,0)
##print(kwant.greens_function(sys,en).transmission(1,0))
##
#u1, s1, v1 = np.linalg.svd(gf, full_matrices=True)
#print(s1)
#import scipy.io
#x = s_list
#matfile = 's_list3.mat'
#scipy.io.savemat(matfile, mdict={'tau': x}, oned_as='row')
#matdata = scipy.io.loadmat(matfile)
#assert np.all(x == matdata['tau'])

elapsed=time()-t_ini
#pyplot.hist(s_list)
