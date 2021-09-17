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

m2 = .1

en=0.3

width=15
length=180

lat_d = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                              [[0, 0], [0, 1/np.sqrt(3)]],norbs=2)  # Coordinates of the sites

ad, bd = lat_d.sublattices
nnn_hoppings_ad = (((-1, 0), ad, ad), ((0, 1), ad, ad), ((1, -1), ad, ad))
nnn_hoppings_bd = (((1, 0), bd, bd), ((0, -1), bd, bd), ((-1, 1), bd, bd))
nnn_hoppingsd = nnn_hoppings_ad + nnn_hoppings_bd

ll=1
def onsite_d(site):
    return 0 
    
def make_system(width, length, salt):
    def disk_d(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width
    def nnn_d(site1, site2):
        x,y=site1.pos
        if abs(x)<length*ll and y>width-5:
            return 1j * m2 * (1-2*1*uniform(repr(site1),salt))
        else: 
            return 1j *m2
        
    sys = kwant.Builder()
    sys[lat_d.shape(disk_d,(0,0))] = onsite_d
    sys[lat_d.neighbors()] = 1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn_d

    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(lat_d.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(lat_d.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym) #,conservation_law=-sz
    
    lead[lat_d.shape(lead_shape, (0, 0))] = 1 # onsite# 
    lead[lat_d.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

num_cores = multiprocessing.cpu_count()
## parllel computing   log_I
## sometimes it calculates the previous length rather than the new one, restart the console!
def lnt_x(cishu):
    salt=cishu+random()
    sys=make_system(width, length, str(salt))
    
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
        completeName = os.path.join('E:/yuhao/65/', str(cishu*wf.shape[0]+num_channel)+".mat")
        sio.savemat(completeName,myDict,oned_as='row') 
#    wf1=np.mean(np.log(np.abs(wf)**2),axis=0)
#    ans=[]
#    for k in range(len(wf1)):
#        ans.append(np.append(sys.pos(k),wf1[k]))
#    myDict = {'ans':ans}
#    completeName = os.path.join('C:/Users/ykang/Documents/yuhao_share/17/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 
    

Parallel(n_jobs=num_cores-2)(delayed(lnt_x)(cishu=j) for j in np.arange(0,1000,1)) 
     
## parllel computing    spectrum
#energies=np.linspace(0,0.8,50)
#salt='ed'
#def con(en):    
#    sys=make_system(length=100)
#    return kwant.smatrix(sys,en).transmission(1,0) 
#
#conduc = Parallel(n_jobs=num_cores-2)(delayed(con)(i) for i in energies)
#pyplot.plot(energies,conduc,'.')


#sys=make_system(width, length, 'r')
#
####gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#t3=kwant.smatrix(sys,en).transmission(1,0) 
#wff=kwant.wave_function(sys,en)
#wf=wff(0)
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5)  #,fig_size=(25, 10)

#wf1=np.mean(np.log(np.abs(wf)**2),axis=0)
#kwant.plotter.map(sys, wf1,num_lead_cells=5)  

elapsed=time()-t_ini