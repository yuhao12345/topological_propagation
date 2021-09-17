# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:49:17 2018

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import kwant
#from kwant.digest import uniform
from time import time
import scipy.io as sio
import os.path
import long_range_correlation_square as corr

#salt=1
t_ini=time()
width=50
length=100

n=1
n1=1
n2=1


#dis=0.15



def onsite(site):
    x,y=site.pos
    return 4+pattern[int(x),int(y)]

def disk(pos):
    x,y=pos
    return 0<=x<length and 0<=y<width

def edge(pos):
    x,y=pos
    return 0<=y<width

def make_system():
    sys=kwant.Builder()
    lat=kwant.lattice.square()
    sys[lat.shape(disk,(0,0))]=onsite
    sys[lat.neighbors()]=-1
    
    lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
    lead[lat.shape(edge,(0,1))]=4
    lead[lat.neighbors()]=-1
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    return sys.finalized()
en=1

save_path = 'C:/Users/ykang/Documents/correlation eigenchannel/2/'
for cishu in np.arange(500,1000):
    salt=cishu+random.randint(0,10000)
    pattern=corr.correlation_pattern_guassian(length,width,5,2,salt)
#    pattern=corr.uniformrandom(length,width,2,salt)
    sys=make_system()
    scat=kwant.smatrix(sys,en)
    gf_mode=scat.submatrix(1,0)
    
#    ggf=np.abs(gf_mode)**2
#    
#    gf_mode2=scat.submatrix(0,0)
#    ggf2=np.abs(gf_mode2)**2
    try:
        u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
        wf=kwant.wave_function(sys,en)(0)
        T=scat.transmission(1,0)
        completeName = os.path.join(save_path, str(cishu)+".mat")
        myDict = {'u':u,'v':v,'s':s,'wf':wf,'T':T,'gf_mode':gf_mode,'pattern':pattern}  #,'gf_mode':gf_mode
        
        sio.savemat(completeName,myDict,oned_as='row')
    except:
        continue

#mode_number=gf_mode.shape[0]
#mode_v=-scat.lead_info[0].velocities[0:mode_number] 
#completeName1 = os.path.join(save_path, "velocity.mat")   
#myDict = {'mode_v':mode_v}  
#sio.savemat(completeName1,myDict,oned_as='row')      
    #['E:\\kwant\\1\\',str(cishu),'.mat']
#    final_T.append(kwant.smatrix(sys,en).transmission(1,0))
#    final_tau_sqrt.append(s)
#    a_test=abs(wf[0,:].reshape((length+1,width-1)))**2
#    print(s)

#    if s[0]>0.997:
#        f_tau1=abs(np.array(sum(np.dot(np.matrix(np.diag(v[0,:])).H,wf))))**2
#        ff_tau1=f_tau1.reshape((length+1,width-1))
#        plt.figure()
#        plt.imshow(ff_tau1)
#
#        shape_tau1=np.sum(ff_tau1,axis=1)
#        shape.append(shape_tau1)
#
#
#    
#    kwant.plotter.map(sys,(abs(np.array(wf[0,:]))**2))
##    local_dos = kwant.ldos(sys, en)
##    kwant.plotter.map(sys, local_dos, num_lead_cells=0)
##point_pos=np.array([sys.pos(i) for i in range(sys.graph.num_nodes)])  #get position of sites
##temp=np.array([abs(np.array(wf[2,:]))**2])
##field=np.concatenate((point_pos,temp.T),axis=1)
#plt.figure()
#plt.plot(np.mean(np.array(shape),axis=0))
elapsed=time()-t_ini

#kwant.plotter.map(sys,(abs(np.array(wf[2,:]))**2))

#sites = sys.leads[0].sites
#mode = scat.lead_info
#psi = mode[0].wave_functions[:, 2]
#psi2 = dict(zip(sites, psi))
#print(scat.transmission(0,0))
#kwant.plotter.map(sys,(abs(np.array(wf[0,:]))**2), fig_size=(20, 10))


