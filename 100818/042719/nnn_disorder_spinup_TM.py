# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 20:48:07 2019

@author: ykang
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
length=1000

di=10  # disorder: w-di<y<w

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  

alpha=0.5  # probability of push down 

def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1

    def nnn(site1, site2):
        x,y=site1.pos
#        return 1j*m2
        if abs(x)<length and y>width-di:
#            return 1j * m2 * (1-2*uniform(repr(site1.tag),salt)*1)
            return 1j*m2*np.sign(uniform(repr(site1.tag),salt)-alpha)
        else:
            return 1j *m2 
    def onsite(site):
        x,y=site.pos
#        return (uniform(repr(site),salt)-0.5)*6
        return 0
#        if abs(x)<length and y>width-di:
#            return (uniform(repr(site),salt)-0.5)*4
#        else:
#            return 0

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn #1j *m2 * sz #
    return sys

def attach_lead(sys):
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width 
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 #1
    lead[graphene.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    

def TM(cishu):
    salt=cishu+random.random()*100
    syst=make_system(width, length, str(salt))
    attach_lead(syst)
    sys = syst.finalized()
    #kwant.plot(sys, fig_size=(10, 3))
#    wf=kwant.wave_function(sys,0.35)(0)
#    kwant.plotter.map(sys, (abs(wf[14])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
#    
#    t=kwant.smatrix(sys,.35).transmission(1,0)
    
    gf_mode=kwant.smatrix(sys,0.4).submatrix(1,0) 
    s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False) #u, s, vh
    ldos = kwant.ldos(sys, 0.4)
    if cishu==0:
        coord=np.array([sys.pos(i) for i in range(ldos.shape[0])])
        sio.savemat('E:/dwell3/97/coord.mat', {'coord':coord},oned_as='row') 
    myDict = {'ld':ldos,'s':s} 
    completeName = os.path.join('E:/dwell3/97/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row')      
#    kwant.plotter.map(sys, ldos)


## delay time from TM
    
f0=0.33
f1=0.35
df=1e-8   # It is freq not omega !!!
f_points=10
energies=np.linspace(f0,f1,f_points)
energies1=energies+df
    
def uv(energies,sys):  # get u,v from TM or TM_gf, TM_gf is right, but TM does not work
    u_mat=[]
    v_mat=[]
    u2_mat=[]
    v2_mat=[]
    for en in energies:
        gf=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!

        u, s, vh=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        u_mat.append(u[:,0])
        v_mat.append(vh.conj().T[:,0])
        u2_mat.append(u[:,1])
        v2_mat.append(vh.conj().T[:,1])
    return u_mat,v_mat ,u2_mat,v2_mat

def delay_tm(cishu):
    salt=cishu+random.random()*100
    syst=make_system(width, length, str(salt))
    attach_lead(syst)
    sys = syst.finalized()
#    u_mat,v_mat=uv(energies,sys)  
#    u_mat1,v_mat1=uv(energies1,sys)  
    u_mat,v_mat,u2_mat,v2_mat=uv(energies,sys)  
    u_mat1,v_mat1,u2_mat1,v2_mat1=uv(energies1,sys)   
    u0=np.column_stack(u_mat)
    v0=np.column_stack(v_mat)
    u1=np.column_stack(u_mat1)
    v1=np.column_stack(v_mat1)
    dphi_tm=np.diag((u0.conj().T @ (u1-u0)-v0.conj().T @ (v1-v0)))/1j/df/(2*np.pi)
    
    u20=np.column_stack(u2_mat)   # second eigenstate
    v20=np.column_stack(v2_mat)
    u21=np.column_stack(u2_mat1)
    v21=np.column_stack(v2_mat1)
    dphi2_tm=np.diag((u20.conj().T @ (u21-u20)-v20.conj().T @ (v21-v20)))/1j/df/(2*np.pi)
    
    return np.append(dphi_tm,dphi2_tm) #dphi_tm # 

## delay time from gf

#def phase_gf(energies,sys):
#    gf=[]
#    for en in energies:
#        gf10=kwant.greens_function(sys,en).submatrix(1,0)
#        gf10_uu=gf10[0,0]
#        gf.append(gf10_uu)
#    return gf
#
#def delay_gf(cishu):      # from gf
#    salt=cishu+random.random()*100
#    syst=make_system(width, length, str(salt))
#    attach_lead(syst)
#    sys = syst.finalized()
#    gf0=phase_gf(energies, sys)
#    gf1=phase_gf(energies1, sys)    
#    dphi_gf=(np.angle(gf1)-np.angle(gf0))/df/(2*np.pi)
#
#    return dphi_gf
    
#dwell=delay_tm(0)
dwell=Parallel(n_jobs=10)(delayed(delay_tm)(cishu=j) for j in np.arange(0,500)) 

#s=TM(0)
#s=Parallel(n_jobs=5)(delayed(TM)(n) for n in np.arange(0,2000))

#dwell=delay_gf(0)
#dwell=Parallel(n_jobs=12)(delayed(delay_gf)(cishu=j) for j in np.arange(0,1000)) 
#myDict = {'s':s}

myDict={'dwell':dwell}
completeName = 'E:/dwell3/106.mat'
sio.savemat(completeName,myDict,oned_as='row') 


elapsed=time()-t_ini



