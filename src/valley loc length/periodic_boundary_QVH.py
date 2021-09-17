# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 23:39:48 2018

@author: user     QVH periodic bc
"""

import numpy as np
import kwant
import matplotlib.pyplot
import random
import scipy.io as sio
import os.path
from kwant.digest import uniform
import time

start=time.time()
#def main() :
X = 10
Y = 4
t=1.0
m=.4
dis=.00
en=0.1

#salt='kyh'
graphene = kwant.lattice.honeycomb(1,'b')
A,B=graphene.sublattices

def rectangle(pos):
    x, y = pos
    return -X/2 < x < X/2

def onsite(site):
    onsite_a = m+(uniform(repr(site),salt)-.5)*dis
    onsite_b = -m+(uniform(repr(site),salt)-.5)*dis
    return onsite_a if site.family == A else onsite_b

def onsite_lead(site):
    onsite_a = m
    onsite_b = -m
    return onsite_a if site.family == A else onsite_b

sym = kwant.TranslationalSymmetry(graphene.vec((-Y/2,Y)))
anc = kwant.Builder(sym) ### 2D periodic conditions
anc[graphene.shape(rectangle,(0, 0))] = None
anc[graphene.neighbors()] = None 
sys = kwant.Builder()
  
sys[anc.sites()] = onsite
sys[((a, sym.to_fd(b)) for a, b in anc.hoppings())] = -t 
#to_fd(a, b=None)[source]
#Map a site or hopping to the fundamental domain.
#If b is None, return a site equivalent to a within the fundamental 
#domain. Otherwise, return a hopping equivalent to (a, b) 
#but where the first element belongs to the fundamental domain.
sym_anc = kwant.TranslationalSymmetry(graphene.vec((1,0)),graphene.vec((-Y/2,Y)))
anc_left = kwant.Builder(sym_anc)

sym_left = kwant.TranslationalSymmetry(graphene.vec((1,0)))
lead_left = kwant.Builder(sym_left)


anc_left[graphene.shape(lambda p: True,(0, 0))] = None
anc_left[graphene.neighbors()] = None

lead_left[anc_left.sites()] = onsite_lead
lead_left[((a, sym.to_fd(b)) for a, b in anc_left.hoppings())] = -t 

sys.attach_lead(lead_left)
sys.attach_lead(lead_left.reversed())

sys = sys.finalized() 
kwant.plot(sys)

#kwant.plot(lead_left)



#save_path = 'E:/loc scale/6/'
#for cishu in np.arange(1,2):
#    salt=str(cishu+random.random())
#    try:
##        u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
##        wf=kwant.wave_function(sys,en)(0)
#        
#        gf_mode=kwant.smatrix(sys,en)
#        T=gf_mode.transmission(1,0)
#        s00=gf_mode.submatrix(0,0)
#        s10=gf_mode.submatrix(1,0)
#        s01=gf_mode.submatrix(0,1)
#        s11=gf_mode.submatrix(1,1)
#        completeName = os.path.join(save_path, str(cishu)+".mat")
#        myDict = {'s00':s00,'s01':s01,'s10':s10,'s11':s11,
#                  't':T,'salt':salt,'length':X,'width':Y}
#        sio.savemat(completeName,myDict,oned_as='row')
#    except:
#        continue


#wf=kwant.wave_function(sys,en)(1)
#kwant.plotter.map(sys,(abs(np.array(wf[0,]))**2))

end=time.time()
print(end-start)