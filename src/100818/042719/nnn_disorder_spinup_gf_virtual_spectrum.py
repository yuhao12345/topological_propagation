# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 15:04:22 2019

@author: user
"""

import numpy as np
from matplotlib import pyplot
import scipy.io as sio
import os.path
import random
import kwant
from kwant.digest import uniform
from random import choices
from time import time
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()

width=30
length=200


dis=1.6
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  


def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and -length<x<length   #25.1

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*dis


    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    

#def mount_vlead(sys, vlead_interface, norb):
#    """Mounts virtual lead to interfaces provided.
#
#    :sys: kwant.builder.Builder
#        An unfinalized system to mount leads
#    :vlead_interface: sequence of kwant.builder.Site
#        Interface of lead
#    :norb: integer
#        Number of orbitals in system hamiltonian.
#    """
#    dim = len(vlead_interface)*norb
#    zero_array = np.zeros((dim, dim), dtype=float)
#    def selfenergy_func(energy, args=()):
#        return zero_array
#
#    vlead = kwant.builder.SelfEnergyLead(selfenergy_func, vlead_interface)
#    sys.leads.append(vlead)

########### for one configuration
#di=11
#def make_system0(width, length):
#    def disk(pos):
#        x,y=pos
#        return width-di<y<width and -length<x<3*length   #25.1
#    sys=kwant.Builder()
#    sys[graphene.shape(disk,(300,width-1))]= 0  #0 #
#    sys[graphene.neighbors()]=1      # comment it, when has rashba
#    def lead_shape(pos):
#        x,y=pos
#        return width-di<y<width 
#    sym = kwant.TranslationalSymmetry((-1,0))
#    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
#    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])    
#    lead = kwant.Builder(sym) #,conservation_law=-sz    
#    lead[graphene.shape(lead_shape, (0, width-1))] = 0 #1
#    lead[graphene.neighbors()]=1       
#    sys.attach_lead(lead)
#    sys.attach_lead(lead.reversed())
#
#    return sys

#syst0=make_system0(width, length)   # part of system

#en=0.4
syst=make_system(width, length, '2')   # whole system as virtual lead
attach_lead(syst)
##kwant.plot(syst0,fig_size=(25, 10))
#
#greens_function_sites = syst0.sites()
#mount_vlead(syst, greens_function_sites, 1)
sys = syst.finalized()

#coord=np.array([sys.pos(i) for i in sys.lead_interfaces[2]]) 
#myDict = {'coord':coord} #,'ld':ld
#completeName = os.path.join('E:/dwell3/797/coor.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#kwant.plot(sys)
### step 1, gf spectrum

#energies=np.linspace(0.3924,0.3978,8192)
#energies=np.linspace(0.3870,0.3924,8192)
#energies=np.linspace(0.3924,0.3978,1000)
energies=np.linspace(0.3978,0.39915,2048)
def gf_01(cishu):
    en=energies[cishu]
    gf=kwant.greens_function(sys,en).submatrix(1,0)[:,0]
#    gf=kwant.smatrix(sys,en).submatrix(1,0)[0,0]
    myDict = {'gf':gf} #,'ld':ld
    completeName = os.path.join('E:/dwell3/811/', str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#    return gf
#Parallel(n_jobs=5)(delayed(gf_01)(cishu) for cishu in np.arange(0,2048,1))
    
#myDict = {'spec':spec} #,'ld':ld
#completeName = os.path.join('E:/dwell3/782.mat')
#sio.savemat(completeName,myDict,oned_as='row') 
    
#def g_01(cishu):
#    en=energies[cishu]
#    g=kwant.smatrix(sys,en).transmission(1,0) 
#    return g

#g=Parallel(n_jobs=10)(delayed(g_01)(cishu) for cishu in np.arange(0,1000,1))

#myDict = {'g':g} #,'ld':ld
#completeName = os.path.join('E:/dwell3/751/9.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#gf_list=Parallel(n_jobs=10)(delayed(gf_01)(cishu) for cishu in range(200))
#pyplot.plot(np.abs(gf_list)**2,'.')

#myDict = {'energies':energies,'gf':gf} #,'ld':ld
#completeName = os.path.join('E:/dwell3/751/1.mat')
#sio.savemat(completeName,myDict,oned_as='row') 


### step2  wf
en=0.39   #energies[667] #energies[522]
wf=kwant.wave_function(sys,en)(0)   #.403586
kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10)) #,colorbar=False

#gf10=kwant.greens_function(sys,.4).submatrix(1,0)
#gf01=kwant.greens_function(sys,.4).submatrix(0,1)
    
#coord=np.array([sys.pos(i) for i in range(55200)])
#myDict = {'wf':wf,'coord':coord} #,'ld':ld
#completeName = os.path.join('E:/dwell3/807/t.mat')
#sio.savemat(completeName,myDict,oned_as='row') 
##    
#def wf_01(cishu):
#    en=energies[cishu]
#    wf=kwant.wave_function(sys,en)(0)[0]
#    myDict = {'wf':wf} #,'ld':ld
#    completeName = os.path.join('E:/dwell3/794/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 
#    return wf

#Parallel(n_jobs=10)(delayed(wf_01)(cishu) for cishu in np.arange(10,100,1))



def gf_virtual(cishu):
    en = energies[cishu]
    gf20 = kwant.greens_function(sys,en).submatrix(2,0)[:,0]
#    gf2=gf.submatrix(2,2)
#    gf22=gf2[:,3406]
#    ld=np.diag(gf.submatrix(2,2))
    myDict = {'gf20':gf20} 
    completeName = os.path.join('E:/dwell3/797/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 


#gf_virtual(1)
#Parallel(n_jobs=3)(delayed(gf_virtual)(n) for n in np.arange(0,3))

#coord_inject=np.array([sys0.pos(i) for i in sys0.lead_interfaces[0]]) 
#coord=np.array([sys0.pos(i) for i in sys0.lead_interfaces[2]])  # sequence of virtual lead

#gf20=kwant.greens_function(sys0,en).submatrix(2,0)[:,0]
#myDict = {'coord':coord, 'coord_inject':coord_inject, 'gf20':gf20} #,'ld':ld
#completeName = os.path.join('E:/dwell3/750/', str(0)+".mat")
#sio.savemat(completeName,myDict,oned_as='row') 


##########  spectrum for many conf
#energies=np.linspace(0.35,0.45,50)
#def gf_01_2(cishu):
#    salt=cishu#+random.random()*100
#    syst=make_system(width, length, str(salt))   # whole system as virtual lead
#    attach_lead(syst)
#    sys = syst.finalized()
#    gf=[kwant.greens_function(sys,en).submatrix(1,0)[0,0] for en in energies]
#    myDict = {'gf':gf} #,'ld':ld
#    completeName = os.path.join('E:/dwell3/751/', str(cishu)+'.mat')
#    sio.savemat(completeName,myDict,oned_as='row') 
##    return gf
#gf_01_2(1)
#Parallel(n_jobs=2)(delayed(gf_01_2)(cishu) for cishu in np.arange(0,4,1))


elapsed=time()-t_ini
#sys.lead_interfaces[1]
#sys.pos(sys.lead_interfaces[1][0])

## plot onsite energy

##sys.sites[1]  #return sites information
#coo=np.array([sys.pos(i) for i in range(27600)])
#potential=[(uniform(repr(sys.sites[i]),'2')-0.5)*dis for i in range(27600)]
#myDict = {'coo':coo, 'potential':potential} #,'ld':ld
#completeName = os.path.join('E:/dwell3/790_onsite.mat')
#sio.savemat(completeName,myDict,oned_as='row') 