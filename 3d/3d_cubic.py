# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:02:25 2018

@author: user
"""
import kwant
import random
import numpy as np
from time import time
import matplotlib.pyplot as plt

t_ini=time()

lat = kwant.lattice.general([(0, 0, 1), (0, 1, 0), (1, 0, 0)],
                            [(0, 0, 0)])

hop=0.5
def make_cuboid(a):
    def cuboid_shape(pos):
        x, y, z = pos
        return 0 <= x < a and 0 <= y < a and 0 <= z < a
    
    def onsite(pos):
        return (random.random()-0.5)*dis

    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0, 0))] = onsite
    sys[lat.neighbors()] = hop
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0,0))
    lead = kwant.Builder(sym_lead)

    def lead_shape(pos):
        x, y, z = pos
        return x==0 and 0 <= y < a and 0 <= z < a

    lead[lat.shape(lead_shape, (0, 0, 0))] = 0
    lead[lat.neighbors()] = hop

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

en=2
tt=[]
#disorder=np.arange(2,16,2)
#for dis in disorder:
#    t_list=[]
#    for cishu in range(50):
#        sys = make_cuboid(a=4)
#    
#        gf_mode = kwant.smatrix(sys, en)
#        T=gf_mode.transmission(1,0)
#        t_list.append(T)
#    tt.append(t_list)
#        
##plt.plot(np.arange(2,16,2),np.log(tt))  
##plt.ylabel('log T')
##plt.xlabel('disorder')
##plt.rc('font', size=5)  
##plt.tick_params(labelsize=14)
#log_mean=np.mean(np.log(tt),1)   
##plt.figure()
##plt.plot(np.arange(2,16,2),log_mean) 
##plt.ylabel('<log T>')
##plt.xlabel('disorder')
##plt.rc('font', size=5)  
##plt.tick_params(labelsize=14)
#
#tmp4=log_mean
#
#plt.plot(disorder,tmp4)
##plt.plot(disorder,tmp6)
##plt.plot(disorder,tmp8)
##plt.plot(disorder,tmp10)
##plt.plot(disorder,tmp12)
#plt.ylabel('<log T>')
#plt.xlabel('disorder')
#plt.rc('font', size=5)  
#plt.tick_params(labelsize=14)

dis=6
size=np.arange(4,16,1)
for a in size:
    t_list=[]
    for cishu in range(50):
        sys = make_cuboid(a)
    
        gf_mode = kwant.smatrix(sys, en)
        T=gf_mode.transmission(1,0)
        t_list.append(T)
    tt.append(t_list)
        
#plt.plot(np.arange(2,16,2),np.log(tt))  
#plt.ylabel('log T')
#plt.xlabel('size')
#plt.rc('font', size=5)  
#plt.tick_params(labelsize=14)
log_mean=np.mean(np.log(tt),1)   
plt.figure()
plt.plot(size,log_mean,'o') 
plt.ylabel('<log T>')
plt.xlabel('size')
plt.rc('font', size=5)  
plt.tick_params(labelsize=14)

#tmp10=log_mean
#plt.plot(np.arange(2,16,2),tmp6)
#plt.plot(np.arange(2,16,2),tmp14)
#plt.plot(np.arange(2,16,2),tmp10)
elapsed=time()-t_ini
#wf=kwant.wave_function(sys,en)(0)
#kwant.plot(sys,site_size=0.05, site_lw=0.05, hop_lw=0.05)

#kwant.plotter.map(sys,(abs(np.array(wf[0,]))**2), fig_size=(20, 10))  only works in 2d

#sys.lead_interfaces
#sys.pos(14)
#kwant.greens_function(sys,en).submatrix(1,0)[67,67] 
