# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 17:31:45 2018

@author: user
"""

import numpy as np
from matplotlib import pyplot
import random
import kwant
from kwant.digest import uniform
import scipy.io as sio
import os.path
from time import time

t_ini=time()

yl, yu = -8, 8
d = 5  #disorder depth
dis = 0   # nnn disorder
dis2 =.5  # onsite disorder    .2
m1, m2 = .5, .5
en=2.5

a=1


s0=np.identity(2)
#sx=np.array([[0,1],[1,0]])
sz=np.array([[1,0],[0,-1]])


graphene = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=graphene.sublattices
def disk(pos):
    x,y=pos
    return abs(y)<8 and np.sqrt(3)*x-300<y<np.sqrt(3)*x+10

def edge(pos):
    x,y=pos
    return abs(y)<8 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

def onsite_qvh(site):
    onsite_a = m1+dis2*(uniform(repr(site),salt)-.5)
    onsite_b = -m1+dis2*(uniform(repr(site),salt)-.5)
    return onsite_a if site.family == A else onsite_b
        

def onsite_qvh_lead(site):
    onsite_a = m1
    onsite_b = -m1
    return onsite_a if site.family == A else onsite_b      

    
sys=kwant.Builder()
sys[graphene.shape(disk,(0,0))]=onsite_qvh
sys[graphene.neighbors()]=-1

#kwant.plot(sys)
lead=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
lead[graphene.shape(disk,(0,0))]=onsite_qvh_lead
lead[graphene.neighbors()]=-1

     
#kwant.plot(lead)
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())

#pyplot.rcParams["figure.figsize"] = [20,15]
#kwant.plot(sys, fig_size=(30, 12))
#pyplot.savefig('testplot.png')

sys=sys.finalized()

#sca=kwant.smatrix(sys,en)
#
#wff=kwant.wave_function(sys,0)

#
#print(sca.transmission(1,0))
#temp=sca.submatrix(1,0)

s_list=[]
t_list=[]
save_path = 'E:/qv_qs_qv/'
for cishu in np.arange(0,300):
    salt=str(cishu/2+random.random())
    try:
        wff=kwant.wave_function(sys,en)
        sca=kwant.smatrix(sys,en)
        gf_mode=sca.submatrix(1,0)      
        u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
        T=sca.transmission(1,0)  
        s_list.append(s)
        t_list.append(T)
#        myDict = {'s':s,'t':T,'salt':salt}
#        completeName = os.path.join(save_path, str(cishu)+".mat")
#        sio.savemat(completeName,myDict,oned_as='row')
    except:
        continue
myDict = {'s':s_list,'t':t_list}
completeName = os.path.join(save_path, str(24)+".mat")
sio.savemat(completeName,myDict,oned_as='row')   

#wf=wff(0)
##wavef1=[]
##for i in range(wf.shape[1]//2):
##    wavef1.append(wf[0,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wf[2,]))**2), fig_size=(20, 10))

#energy=np.linspace(0.03,.05,300)  
#spec=[kwant.greens_function(sys,en).submatrix(0,0)[4,4] for en in energy]  #4,4 and 18,18 spin up/down
#ang=np.unwrap(np.angle(spec))
#trans=[abs(sp)**2 for sp in spec]
#pyplot.figure(figsize=(10,6))
#pyplot.plot(energy,ang,'o')
#pyplot.title('phase')
##temp=kwant.greens_function(sys,0).submatrix(1,0)
#pyplot.figure()
#pyplot.figure(figsize=(10,6))
#pyplot.plot(energy,trans,'o')
#pyplot.title('intensity')

#temp=kwant.greens_function(sys,.3).out_block_coords(1)
#print(kwant.greens_function(sys,.3).num_propagating(0))
#
#condu=[kwant.greens_function(sys,en).transmission(1,0) for en in energy]
#pyplot.figure()
#pyplot.plot(energy,condu)

#matr=kwant.greens_function(sys,0).submatrix(1,0)
#print(matr.max())
#print(np.unravel_index(np.argmax(matr, axis=None), matr.shape))
#import scipy.io
#x = s_list
#matfile = 's_list3.mat'
#scipy.io.savemat(matfile, mdict={'tau': x}, oned_as='row')
#matdata = scipy.io.loadmat(matfile)
#assert np.all(x == matdata['tau'])
#sys.pos(2373)    #23  2373
elapsed=time()-t_ini