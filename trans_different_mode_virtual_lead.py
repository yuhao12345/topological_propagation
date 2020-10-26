# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:49:17 2018

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import kwant
from kwant.digest import uniform
import time
import scipy.io as sio
import os.path

start=time.time()
width=5
length=15
n1=1
n2=1
n=1
sys=kwant.Builder()
lat=kwant.lattice.square()



def onsite(site):
#    return 4/(n+.3*uniform(repr(site),salt))**2
#    return 4/(n+.0*(np.sign(uniform(repr(site),salt)-0.5)+1))**2
    return 4/(n)**2


def disk(pos):
    x,y=pos
    return abs(x)<length and abs(y)<width

def edge(pos):
    x,y=pos
    return x==0 and abs(y)<width

def scatter():
    sys[lat.shape(disk,(1,1))]=onsite
    sys[lat.neighbors()]=-1/n**2
    
    lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
    lead[lat.shape(edge,(0,1))]=4/n1**2
    lead[lat.neighbors()]=-1/n1**2
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    #lead=kwant.Builder(kwant.TranslationalSymmetry((1,0))) 
    #lead[lat.shape(edge,(0,1))]=4/n2**2
    #lead[lat.neighbors()]=-1/n2**2
    #sys.attach_lead(lead)
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
    return y==0 and abs(x)<length

#vlead = kwant.Builder(kwant.TranslationalSymmetry((a/2, a/2*np.sqrt(3))))
vlead=kwant.Builder()
vlead[lat.shape(vshape, (0,0))] = 0    # just define shape, self energy is defined in kwant.builder.SelfEnergyLead
    
gf_sites = vlead.sites()    
syst = scatter()
mount_vlead(syst, gf_sites, 1)
sys = syst.finalized()

en=1
sys.pos(148)
lead_site=sys.lead_interfaces
temp=kwant.greens_function(sys,en).submatrix(2,2)
#local_dos = kwant.ldos(sys, en)
#shape=[]
#final_T=[]
#final_tau_sqrt=[]
#save_path = 'E:/kwant/26/'
#for cishu in np.arange(1):
#    salt=str(cishu+random.random())
#    gf_mode=kwant.smatrix(sys,en).submatrix(1,0)  
#    try:
#        u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
#        wf=kwant.wave_function(sys,en)(0)
#        T=kwant.smatrix(sys,en).transmission(1,0)
#        completeName = os.path.join(save_path, str(cishu)+".mat")
#        myDict = {'u':u,'v':v,'s':s,'wf':wf,'T':T,'gf_mode':gf_mode}
#    ##    myDict['u'] = u
#        sio.savemat(completeName,myDict,oned_as='row')
#    except:
#        continue
            
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
#    kwant.plotter.map(sys, local_dos, num_lead_cells=0)
##point_pos=np.array([sys.pos(i) for i in range(sys.graph.num_nodes)])  #get position of sites
##temp=np.array([abs(np.array(wf[2,:]))**2])
##field=np.concatenate((point_pos,temp.T),axis=1)
#plt.figure()
#plt.plot(np.mean(np.array(shape),axis=0))
end=time.time()
print(end-start)
#print(np.mean(final_T))
#plt.imshow
