# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 17:20:10 2018

@author: user
"""
# https://kwant-project.org/doc/1.0/tutorial/tutorial5#id1
# set spin up/down as two layers
#parallel calculation
import kwant
import numpy as np
import scipy.io as sio
import os.path
from matplotlib import pyplot
from kwant.digest import uniform
from joblib import Parallel, delayed
import multiprocessing
from time import time
from random import random

t_ini=time()

m1 = .4  #valley
m2 = .1
bo=0  #boundary between qvh and qsh
en=0.15
width=20

lat_d = kwant.lattice.honeycomb(a=1, name='d')  #spin down

ad, bd = lat_d.sublattices
nnn_hoppings_ad = (((-1, 0), ad, ad), ((0, 1), ad, ad), ((1, -1), ad, ad))
nnn_hoppings_bd = (((1, 0), bd, bd), ((0, -1), bd, bd), ((-1, 1), bd, bd))
nnn_hoppingsd = nnn_hoppings_ad + nnn_hoppings_bd


def onsite_d(site):
    x,y=site.pos
    if y>bo:
        return 0   
    else:
        return m1 * (1 if site.family == ad else -1)
    


def nnn_d_lead(site1, site2):
    x,y=site1.pos
    if y<=bo:
        return 0
    else: 
        return -1j *m2
    
def make_system(length):
    def disk_d(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1
    def nnn_d(site1, site2):
        x,y=site1.pos
        if y<=bo:
            return 0
    #    elif x**2+(y+delta-3)**2<9:
        elif y<10 and abs(x)<length:
            return 1j * m2 * 2*(uniform(repr(site1),salt)-0.5)
        else: 
            return -1j *m2
        
    sys = kwant.Builder()
    sys[lat_d.shape(disk_d,(0,0))] = onsite_d
    sys[lat_d.neighbors()] = 1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn_d

    sym_left = kwant.TranslationalSymmetry((-1, 0))
    lead1 = kwant.Builder(sym_left)
    lead1[lat_d.shape((lambda pos: abs(pos[1]) < width), (0, 0))]= onsite_d#1
    lead1[lat_d.neighbors()] = 1
    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn_d_lead   
    

    sys.attach_lead(lead1)   #spin down
    sys.attach_lead(lead1.reversed())   #lead4 spin down

    return sys.finalized()

num_cores = multiprocessing.cpu_count()
## parllel computing   log_I
## sometimes it calculates the previous length rather than the new one, restart the console!
def lnt_x(cishu):
    global salt
    salt=str(cishu+random())
    sys=make_system(length=100)
    
#    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#    t=kwant.smatrix(sys,en).transmission(1,0) 
    wff=kwant.wave_function(sys,en)
    wf=wff(0)    
    for num_channel in range(wf.shape[0]):
        wf1=wf[num_channel,]
        ans=[]
        for k in range(len(wf1)):
            ans.append(np.append(sys.pos(k),wf1[k]))
        myDict = {'ans':ans}
        completeName = os.path.join('E:/yuhao/62/', str(cishu*wf.shape[0]+num_channel)+".mat")
        sio.savemat(completeName,myDict,oned_as='row') 
#    wf1=np.mean(np.log(np.abs(wf)**2),axis=0)
#    ans=[]
#    for k in range(len(wf1)):
#        ans.append(np.append(sys.pos(k),wf1[k]))
#    myDict = {'ans':ans}
#    completeName = os.path.join('C:/Users/ykang/Documents/yuhao_share/17/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 
    

Parallel(n_jobs=num_cores-2)(delayed(lnt_x)(cishu=j) for j in np.arange(0,100,1)) 
     
## parllel computing    spectrum
#energies=np.linspace(0,0.8,50)
#salt='ed'
#def con(en):    
#    sys=make_system(length=100)
#    return kwant.smatrix(sys,en).transmission(1,0) 
#
#conduc = Parallel(n_jobs=num_cores-2)(delayed(con)(i) for i in energies)
#pyplot.plot(energies,conduc,'.')



#salt='r'
#sys=make_system(length=100)
## 
###gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
###t3=kwant.smatrix(sys,en).transmission(1,0) 
#wff=kwant.wave_function(sys,en)
#wf=wff(0)
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5)  #,fig_size=(25, 10)

#wf1=np.mean(np.log(np.abs(wf)**2),axis=0)
#kwant.plotter.map(sys, wf1,num_lead_cells=5)  

elapsed=time()-t_ini