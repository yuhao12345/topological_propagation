# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 13:57:58 2018

@author: user
"""

#       Aubry Andre model
import numpy as np
import matplotlib.pyplot as plt
#from random import random
import kwant
import scipy.sparse.linalg as sla

la=2.2
#alpha=np.sqrt(5)/2+0.5
miu=1.2
dx=1
L=1000

def onsite(site):
    n,=site.pos
    return la*np.cos(.2*pow(n,miu)) #np.cos(np.pi*alpha*pow(n,miu))

    
def make_system():
  
   lat=kwant.lattice.chain(dx)
   sys=kwant.Builder()
  
   sys[(lat(i) for i in range(L))]=onsite
   sys[lat.neighbors()]= 1
  
#   sym = kwant.TranslationalSymmetry((-dx,)) #here i also tried to use 
#
#   l_lead=kwant.Builder(sym)
#   l_lead[(lat(0))]=onsite
#   l_lead[lat.neighbors()]= 1
#   sys.attach_lead(l_lead)
#   
#   sym1 = kwant.TranslationalSymmetry((dx,)) #here i also tried to use 
#
#   l_lead1=kwant.Builder(sym1)
#   l_lead1[(lat(L-1))]=onsite
#   l_lead1[lat.neighbors()]= 1
#   sys.attach_lead(l_lead1)
#
   return sys.finalized()

sys=make_system()
#kwant.plot(sys, fig_size=(10, 3))
#
###en=2
#wf=kwant.wave_function(sys,2.3911)(0)
#plt.plot(abs(np.array(wf[0,]))**2)

ham_mat = sys.hamiltonian_submatrix(sparse=True)
vals,evecs  = sla.eigsh(ham_mat, k=1, sigma=1,which='LA', return_eigenvectors=True)
#vals, vecs = sla.eigsh(ham_mat, k=5, which='SM')[1]
#
#    # Plot the probability density of the 10th eigenmode.
plt.plot(np.abs(evecs[:, 0])**2)
#plot_wave_function(sys)
#energy=np.linspace(.1,.5,21)  
#spec=[kwant.greens_function(sys,en).submatrix(1,0)[0,0] for en in energy]  #2,2
##ang=np.unwrap(np.angle(spec))
#trans=[abs(sp)**2 for sp in spec]
#plt.figure()
#plt.plot(energy,trans,'o')

#kwant.greens_function(sys,1.613).submatrix(0,0)
