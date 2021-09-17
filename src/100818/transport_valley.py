# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 21:22:35 2017

@author: Yuhao Kang    %kymwswj QVH  edge mode between up and down triangles
"""

import numpy as np
from matplotlib import pyplot
#from random import random
import kwant
#from time import time
from kwant.digest import uniform
#t_ini=time()

salt='ed4'
en=0
m1 = .2  #valley
dis=0

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

def disk(pos):
    x,y=pos
    return abs(y)<20 and abs(x)<20

def edge(pos):
    x,y=pos
    return abs(y)<10

def onsite_lead(site):
    return 4
#    x,y=site.pos
#    if y>0:
#        return m1 * (1 if site.family == a else -1)  
#    else:
#        return -m1 * (1 if site.family == a else -1)    

def onsite(site):
    x,y=site.pos
#    if x**2+(y-1)**2<5:
#        return .4 * (1 if site.family == a else -1)  

    if y>0:
        return m1* ((uniform(repr(site),salt)-0.5)*dis+1)* (1 if site.family == a else -1)  
    else:
        return -m1 * (1 if site.family == a else -1)   
    
def make_system():   
    sys=kwant.Builder()    
    sys[graphene.shape(disk,(0,0))]=onsite
    sys[graphene.neighbors()]=-1
    lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
    lead[graphene.shape(edge,(0,0))]=onsite_lead
    lead[graphene.neighbors()]=-1
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    return sys.finalized() 
#pyplot.rcParams["figure.figsize"] = [20,10]
#kwant.plot(sys)

      
sys=make_system()
wf=kwant.wave_function(sys,en)(0)
#wavef1=[]    
#for i in range(wf.shape[1]):
#    wavef1.append(wf[0,i])       
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(20, 12))

tpp=(abs(wf)**2)[0,]
kwant.plotter.map(sys,tpp, fig_size=(20, 12))

gf=kwant.greens_function(sys,en).submatrix(1,0)
t1=kwant.greens_function(sys,en).transmission(1,0)
sg=np.linalg.svd(gf, compute_uv=False)
taug1=sg**2

#gf_measure=np.zeros([14,14],dtype=complex)
#for j in range(14):
#    for jj in range(14):
#        gf_measure[j,jj]=gf[3*j+1,3*jj+1]
#sg_m=np.linalg.svd(gf_measure, compute_uv=False)
#taug1_m=sg_m**2    # measured

gf_mode=kwant.smatrix(sys,en).submatrix(1,0)    
s=np.linalg.svd(gf_mode, compute_uv=False)
tau1=s**2
#for times in range(1):    
#  
##    salt='wdlxwws'    
##    sca=kwant.smatrix(sys,en).data
##    sca2=kwant.smatrix(sys,en).data
#    ##print(sca)
#
#    rl=[]
#    rl2=[]    
#    mn=0 
#    gf_mode=kwant.smatrix(sys,en).submatrix(1,0)      
#    u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#    print(s)
#    
#    s_list.append(s)
##    print(mn)
##print(np.array(s_list))
##print(kwant.smatrix(sys,en).transmission(1,0))
###print(kwant.smatrix(sys,en).num_propagating(0))
####
##gf=kwant.greens_function(sys,en).submatrix(1,0)
###print(kwant.greens_function(sys,en).transmission(1,0))
###
##u1, s1, v1 = np.linalg.svd(gf, full_matrices=True)
##print(s1)
##import scipy.io
##x = s_list
##matfile = 's_list3.mat'
##scipy.io.savemat(matfile, mdict={'tau': x}, oned_as='row')
##matdata = scipy.io.loadmat(matfile)
##assert np.all(x == matdata['tau'])
#
#elapsed=time()-t_ini
#pyplot.hist(s_list)
