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

#en=0.3

width=30
length=50

lat_d = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                              [[0, 0], [0, 1/np.sqrt(3)]],norbs=2)  # Coordinates of the sites

ad, bd = lat_d.sublattices
nnn_hoppings_ad = (((-1, 0), ad, ad), ((0, 1), ad, ad), ((1, -1), ad, ad))
nnn_hoppings_bd = (((1, 0), bd, bd), ((0, -1), bd, bd), ((-1, 1), bd, bd))
nnn_hoppingsd = nnn_hoppings_ad + nnn_hoppings_bd


def onsite_d(site):
    return 0 
    
def make_system(width, length):
    def disk_d(pos):
        x,y=pos
        return -width<y<15 and abs(x)<length   #25.1
    def lead_shape(pos):
        x,y=pos
        return -width<y<15
    def nnn_d(site1, site2):
        x,y=site1.pos
        if x**2+(y+10)**2<16:
            return -1j * m2
        elif y<0: 
            return 1j *m2
        else:
            return -1j *m2
    def nnn_lead(site1, site2):
        x,y=site1.pos
        if y<0: 
            return 1j *m2
        else:
            return -1j *m2       
    sys = kwant.Builder()
    sys[lat_d.shape(disk_d,(0,0))] = onsite_d
    sys[lat_d.neighbors()] = 1
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]] = nnn_d

    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(lat_d.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(lat_d.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[lat_d.shape(lead_shape, (0, 0))] = onsite_d
    lead[lat_d.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppingsd]]=nnn_lead
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

num_cores = multiprocessing.cpu_count()
## parllel computing   log_I
## sometimes it calculates the previous length rather than the new one, restart the console!
#def lnt_x(cishu):
#    salt=cishu+random()
#    sys=make_system(width, length, str(salt))
#    
##    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
##    t=kwant.smatrix(sys,en).transmission(1,0) 
#    wff=kwant.wave_function(sys,en)
#    wf=wff(0)    
#    for num_channel in range(wf.shape[0]):
#        wf1=wf[num_channel,]
#        ans=[]
#        for k in range(len(wf1)):
#            ans.append(np.append(sys.pos(k),wf1[k]))
#        myDict = {'ans':ans}
#        completeName = os.path.join('E:/yuhao/65/', str(cishu*wf.shape[0]+num_channel)+".mat")
#        sio.savemat(completeName,myDict,oned_as='row') 
        
#    wf1=np.mean(np.log(np.abs(wf)**2),axis=0)
#    ans=[]
#    for k in range(len(wf1)):
#        ans.append(np.append(sys.pos(k),wf1[k]))
#    myDict = {'ans':ans}
#    completeName = os.path.join('C:/Users/ykang/Documents/yuhao_share/17/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 
    

#Parallel(n_jobs=num_cores-2)(delayed(lnt_x)(cishu=j) for j in np.arange(0,1000,1)) 
     
## parllel computing    spectrum
energies=np.linspace(0.33,0.35,1000)
sys=make_system(width,length)
def con(n):    
    en=energies[n]
    gf=kwant.greens_function(sys,en).submatrix(1,0)
    t=kwant.smatrix(sys,en).transmission(1,0) 
    myDict = {'gf':gf,'t':t}
    completeName = os.path.join('E:/fano/7/', str(n)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
    
spec = Parallel(n_jobs=20)(delayed(con)(n) for n in range(1000))

#X = [x.real for x in spec]
#Y = [x.imag for x in spec]
#pyplot.scatter(X,Y, color='red')

#ang=np.angle(spec)
#r=np.abs(spec)
#pyplot.polar(ang,r,'.')

#pyplot.plot(energies,spec,'.')
#pyplot.title('Transmittance')
##pyplot.ylim((0,1.5))
#pyplot.xlabel('freq')
    
#sys=make_system(width, length)
#en=0.37#0.34025     #0.459--6    0.44--8   x**2+(y+8)**2<9:
#wff=kwant.wave_function(sys,en)
#wf=wff(0)
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10))

#t=kwant.smatrix(sys,en).transmission(1,0)

####gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#t3=kwant.smatrix(sys,en).transmission(1,0) 

#t1=kwant.greens_function(sys,en).transmission(1,0)
#gf_point=kwant.greens_function(sys,en).submatrix(1,0)[0,0]

#energy=np.linspace(0.44,.47,100)  
#spec=[kwant.greens_function(sys,en).submatrix(1,0)[0,0] for en in energy]

#ang1=np.unwrap(np.angle(spec))
#pyplot.figure(figsize=(10,6))
#pyplot.plot(energies,ang1,'.')
#pyplot.title('phase')
#pyplot.xlabel('freq')
#
#trans=[abs(sp)**2 for sp in spec]
#pyplot.figure(figsize=(10,6))
#pyplot.plot(energies,trans,'.')
#pyplot.title('intensity')
    
#sys.lead_interfaces[0][17]
#sys.pos(4909)  #17 
y=[sys.pos(sys.lead_interfaces[0][i])[1] for i in range(sys.lead_interfaces[0].size)]

elapsed=time()-t_ini

#myDict = {'s':spec}
#completeName = os.path.join('C:/Users/ykang/Desktop/', str(0)+".mat")
#sio.savemat(completeName,myDict,oned_as='row') 