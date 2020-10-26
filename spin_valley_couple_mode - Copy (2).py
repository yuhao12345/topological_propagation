# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 22:40:31 2017

@author: user
"""

import numpy as np
from matplotlib import pyplot
from random import random
import kwant
#from kwant.digest import uniform
#from time import time
#from kwant.digest import uniform
#t_ini=time()
m, t=0.4, -0.133
dis=0
#en=0.3
#salt='wdlxwws'
a=1
lat = kwant.lattice.general([[a, 0], [a/2, a*np.sqrt(3)/2]],
                            [[0, 0], [0, a/np.sqrt(3)]])
A,B=lat.sublattices

s0=np.identity(2)
sx=np.array([[0,1],[1,0]])
sz=np.array([[1,0],[0,-1]])

def disk(pos):
    x,y=pos
    return abs(y)<3 and np.sqrt(3)*x-1<y<np.sqrt(3)*x+9

def edge(pos):
    x,y=pos
    return abs(y)<3 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 
#def edge1(pos):
#    x,y=pos
#    return abs(y)<5 and np.sqrt(3)*x-0.1<y<np.sqrt(3)*x+0.58 

def hopp(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return -t*s0
    else:
        return np.zeros([2,2])
    
def hopp2(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    if (y1+0.1*a)*(y2+0.1*a)>0 and y1+0.1*a<0:
        return t*s0
    else:
        return np.zeros([2,2])
    
def onsite1a(site):
    x,y=site.pos
    if y>-0.1*a:
        return m*sz+0*sx
    else:
        return np.zeros([2,2])
    
def onsite2a(site):
    x,y=site.pos
    if y>-0.1*a:
        return -m*sz+0*sx
    else:
        return np.zeros([2,2])
def onsite1(site):
    x,y=site.pos
    if (x-25*a)**2+(y-5*a)**2<0:
        return -m*sz+0*sx
    elif y>-0.1*a:
        return m*sz+0*sx
    else:
        return np.zeros([2,2])
    
def onsite2(site):
    x,y=site.pos
    if (x-25*a)**2+(y-5*a)**2<0:
        return m*sz+0*sx
    elif y>-0.1*a:
        return -m*sz+0*sx
    else:
        return np.zeros([2,2])
     
def scatter():
    sys=kwant.Builder()   
    sys[A.shape(disk,(0,0))]=onsite1
    sys[B.shape(disk,(0,0))]=onsite2
    sys[lat.neighbors()]=-s0*2/3
    sys[A.neighbors()]=hopp
    sys[B.neighbors()]=hopp2
    
    lead0=kwant.Builder(kwant.TranslationalSymmetry((-a,0))) 
    lead0[A.shape(edge,(0,0))]=onsite1a
    lead0[B.shape(edge,(0,0))]=onsite2a
    lead0[lat.neighbors()]=-s0*2/3
    lead0[A.neighbors()]=hopp
    lead0[B.neighbors()]=hopp2
    sys.attach_lead(lead0)
    sys.attach_lead(lead0.reversed())
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
    zero_array = np.zeros((dim, dim), dtype=float)
    def selfenergy_func(energy, args=()):
        return zero_array

    vlead = kwant.builder.SelfEnergyLead(selfenergy_func, vlead_interface)
    sys.leads.append(vlead)

def vshape(pos):
    x, y = pos
    return abs(y)<0.1*a and -2*a>x>-4*a

#vlead = kwant.Builder(kwant.TranslationalSymmetry((a/2, a/2*np.sqrt(3))))
vlead=kwant.Builder()
vlead[lat.shape(vshape, (-3,0))] = 0    # just define shape, self energy is defined in kwant.builder.SelfEnergyLead
    
#vlead[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = 1
gf_sites = vlead.sites()    
syst = scatter()
mount_vlead(syst, gf_sites, 2)
sys = syst.finalized()
    
#pyplot.rcParams["figure.figsize"] = [20,15]
kwant.plot(sys, fig_size=(10, 12))
#pyplot.savefig('testplot.png')
#sca=kwant.smatrix(sys,en).data

#wf=kwant.wave_function(sys,.15)(0)
#wavef1=[]
#for i in range(wf.shape[1]//2):
#    wavef1.append(wf[0,2*i])
#kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(20, 10))

#wavef2=[]
#for i in range(wf.shape[1]//2):
#    wavef2.append(wf[1,2*i+1])
#kwant.plotter.map(sys,(abs(np.array(wavef2))**2))
#s_list=[]
#for conf in range(1):
#    rl=[]
#    mn=0
#gf_mode=kwant.smatrix(sys,1).submatrix(1,0)      
#u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#    print(s)
#    s_list.append(s)
#kwant.greens_function(sys,en).transmission(1,0)  

#energy=np.linspace(0,.1,100)  
#spec=[kwant.greens_function(sys,en).submatrix(1,0)[14,14] for en in energy]  #2,2
#ang=np.unwrap(np.angle(spec))
#trans=[abs(sp)**2 for sp in spec]
#pyplot.figure(figsize=(10,6))
#pyplot.plot(energy,ang,'o')
#
#pyplot.figure()
#pyplot.figure(figsize=(10,6))
#pyplot.plot(energy,trans,'o')
#temp=kwant.greens_function(sys,.3).out_block_coords(1)
#print(kwant.greens_function(sys,.3).num_propagating(0))
#
#condu=[kwant.greens_function(sys,en).transmission(1,0) for en in energy]
#pyplot.figure()
#pyplot.plot(energy,condu)

temp3=kwant.greens_function(sys,0).submatrix(2,2)

#matr=kwant.greens_function(sys,0).submatrix(1,0)
#print(matr.max())
#print(np.unravel_index(np.argmax(matr, axis=None), matr.shape))
#import scipy.io
#x = s_list
#matfile = 's_list3.mat'
#scipy.io.savemat(matfile, mdict={'tau': x}, oned_as='row')
#matdata = scipy.io.loadmat(matfile)
#assert np.all(x == matdata['tau'])

#elapsed=time()-t_ini
#b=kwant.plotter.sys_leads_sites(sys,num_lead_cells=2)
#kwant.plotter.sys_leads_sites(sys,num_lead_cells=1)[1]
#sys.graph.num_nodes
#sys.pos(27)
#kwant.plotter.sys_leads_pos(sys, site_lead_nr=320)
#temp=kwant.greens_function(sys,en).submatrix(1,0)
#temp1=kwant.greens_function(sys,en).out_block_coords(1)
#temp2=kwant.greens_function(sys,en).block_coords(1, 0)
#  %reset   clear variables
lead_site=sys.lead_interfaces
sys.hamiltonian_submatrix()
sys.leads[0].selfenergy(energy=0.01)