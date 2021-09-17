# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 14:23:14 2018

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


#dis = 0   # nnn disorder
#dis2=0  # onsite disorder
#salt='kyh'
en=0.25

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m1 = .5  #valley
m2 = .077   #spin

s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def disk(pos):
    x,y=pos
    return abs(y)<9 and abs(x)<75

#def onsite_s(site):
#    x,y=site.pos
#    if abs(y)<0:
#        return np.zeros([2,2]) 
#    else:
##    return np.zeros([2,2]) 
#        return dis2*(uniform(repr(site),salt)-.5)*s0
def spin_orbit(site1, site2):
    return 1j * (m2+dis*(uniform(repr(site1),salt)-.5)) * sz
def spin_orbit_lead(site1, site2):
    return 1j * m2 * sz   
sys=kwant.Builder()
sys[graphene.shape(disk,(0,0))]=np.zeros([2,2])
sys[graphene.neighbors()]=-s0
sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit

lead=kwant.Builder(kwant.TranslationalSymmetry([-1,0])) 
lead[graphene.shape((lambda pos: abs(pos[1]) < 9), (0, 0))]=np.zeros([2,2])
lead[graphene.neighbors()]=-s0
lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = spin_orbit_lead
     
#kwant.plot(lead)
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())

#pyplot.rcParams["figure.figsize"] = [20,15]
#kwant.plot(sys, fig_size=(30, 12))
#pyplot.savefig('testplot.png')

sys=sys.finalized()

save_path = 'E:/TI transmission simulation/s12/'
k=1
#kwant.wave_function(sys,en)
for dis in np.arange(.25, .55,.05):  # onsite disorder    .2

    s_list=[]
    t_list=[]    
    for cishu in np.arange(0,250):
        salt=str(cishu+random.random())
        try:
#            wff=kwant.wave_function(sys,en)
            sca=kwant.smatrix(sys,en)
            gf_mode=sca.submatrix(1,0)      
            u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
            T=sca.transmission(1,0)  
            s_list.append(s)
            t_list.append(T)
        except:
            continue
    myDict = {'s':s_list,'t':t_list}
    completeName = os.path.join(save_path, str(k)+".mat")
    sio.savemat(completeName,myDict,oned_as='row')   
    k=k+1

#wf=wff(0)
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(20, 10))



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