# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 13:57:58 2018

@author: user
"""

# add virtual lead   the wavefunction cannot be obtained
import numpy as np
import matplotlib.pyplot as plt
#from random import random
import kwant

la=.4
#alpha=np.sqrt(5)/2+0.5
miu=.7
dx=1
L=10
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
   
#   sym1 = kwant.TranslationalSymmetry((dx,)) #here i also tried to use 
#
#   l_lead1=kwant.Builder(sym1)
#   l_lead1[(lat(L))]=onsite
#   l_lead1[lat.neighbors()]= 1
#   sys.attach_lead(l_lead1)
#
#   return sys.finalized()
   return sys

def mount_vlead(sys, vlead_interface, norb):
    """Mounts virtual lead to interfaces provided.

    :sys: kwant.builder.Builder
        An unfinalized system to mount leads
    :vlead_interface: sequence of kwant.builder.Site
        Interface of lead
    :norb: integer
        Number of orbitals in system hamiltonian.
    """
    dim = len(vlead_interface)*norb
    #print(dim)
    zero_array = np.zeros((dim, dim), dtype=float)
    def selfenergy_func(energy, args=()):
        return zero_array

    vlead = kwant.builder.SelfEnergyLead(selfenergy_func, vlead_interface)
    sys.leads.append(vlead)

def vshape(pos):
    x,  = pos
    return x==0 

#vlead = kwant.Builder(kwant.TranslationalSymmetry((a/2, a/2*np.sqrt(3))))
lat=kwant.lattice.chain(dx)
vlead=kwant.Builder()
vlead[lat.shape(vshape, (0,))] = 0    # just define shape, self energy is defined in kwant.builder.SelfEnergyLead
kwant.plot(vlead, fig_size=(10, 3))
#vlead[(lat(i) for i in range(L))]=0 
gf_sites = vlead.sites()    
syst = make_system()
mount_vlead(syst, gf_sites, 1)
sys = syst.finalized()

#sys=make_system()
kwant.plot(sys, fig_size=(10, 3))

kwant.greens_function(sys,1.6).submatrix(0,0)
##en=2
wf=kwant.wave_function(sys,1.613)(0)
##kwant.plotter.map(sys,(abs(np.array(wf))**2), fig_size=(10, 3))
#plt.plot(abs(np.array(wf[0,]))**2)

#energy=np.linspace(.1,.5,21)  
#spec=[kwant.greens_function(sys,en).submatrix(1,0)[0,0] for en in energy]  #2,2
##ang=np.unwrap(np.angle(spec))
#trans=[abs(sp)**2 for sp in spec]
#plt.figure()
#plt.plot(energy,trans,'o')

