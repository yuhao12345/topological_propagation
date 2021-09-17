# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 20:47:26 2019

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
length=50
ratio=1

di=10  # disorder: w-di<y<w

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
        return abs(y)<width and abs(x)<length   #25.1
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width #width-di<y<width  #
#    def lead_shape1(pos):
#        x,y=pos
#        return  abs(y)<width  #width-15<y<width #
    def nnn(site1, site2):
        x,y=site1.pos
        if abs(x)<length and y>width-di:
            return 1j * m2 * (1-2*uniform(repr(site1.tag),salt)*1)
        else:
            return 1j *m2 

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= 0
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn #1j *m2 * sz #

    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 0 #1
    lead[graphene.neighbors()]=1
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
#    sym1 = kwant.TranslationalSymmetry((1,0))
#    sym1.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
#    sym1.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
#    lead1 = kwant.Builder(sym1) #,conservation_law=-sz
#    
#    lead1[graphene.shape(lead_shape1, (0, width-1))] =  0#  s0 # 
#    lead1[graphene.neighbors()]=1
#    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2     
#    sys.attach_lead(lead1)
    
    return sys.finalized()

en=0.5
sys=make_system(width, length, '7527d')
wff=kwant.wave_function(sys,en)
wf=wff(0)
kwant.plotter.map(sys, (abs(wf[0])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

#tm=kwant.smatrix(sys,en)
#t=tm.transmission(1,0)
#gf_mode=tm.submatrix(1,0) 
#s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False)

#f0=0.3
#f1=0.35
#df=1e-11   # It is freq not omega !!!
#f_points=20
#energies=np.linspace(f0,f1,f_points)
#energies1=energies+df

# from TM
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
    sys=make_system(width, length, str(salt))
#    trans=[]
#    for en in energies:
#        trans.append(kwant.smatrix(sys,en).transmission(1,0))
    u_mat,v_mat,u2_mat,v2_mat=uv(energies,sys)  
    u_mat1,v_mat1,u2_mat1,v2_mat1=uv(energies1,sys)   
    u0=np.column_stack(u_mat)
    v0=np.column_stack(v_mat)
    u1=np.column_stack(u_mat1)
    v1=np.column_stack(v_mat1)
    
    u20=np.column_stack(u2_mat)   # second eigenstate
    v20=np.column_stack(v2_mat)
    u21=np.column_stack(u2_mat1)
    v21=np.column_stack(v2_mat1)
    dphi_tm=np.diag((u0.conj().T @ (u1-u0)-v0.conj().T @ (v1-v0)))/1j/df/(2*np.pi)
    dphi2_tm=np.diag((u20.conj().T @ (u21-u20)-v20.conj().T @ (v21-v20)))/1j/df/(2*np.pi)
    return np.append(dphi_tm,dphi2_tm)

#dwell=Parallel(n_jobs=10)(delayed(delay_tm)(cishu=j) for j in np.arange(0,1000,1)) 
#myDict = {'dwell':dwell}
#completeName = 'E:/dwell2/0/9.mat'
#sio.savemat(completeName,myDict,oned_as='row') 

def tau_n(cishu):
    en=0.5
    salt=cishu+random.random()*1000
    sys=make_system(width, length, str(salt))
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
    s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False) #u, s, vh
    return s

#tau_half=Parallel(n_jobs=10)(delayed(tau_n)(cishu=j) for j in np.arange(0,500,1)) 
#myDict = {'tau_half':tau_half}
#completeName = 'E:/dwell2/43.mat'
#sio.savemat(completeName,myDict,oned_as='row') 

def lni(cishu):
    en=0.5
    salt=cishu+random.random()*1000
    sys=make_system(width, length, str(salt))
    wff=kwant.wave_function(sys,en)(0)
    w=np.abs(wff.T)
    coord=np.array([sys.pos(i) for i in range(w.shape[0])])
    up=np.concatenate((coord,w),axis=1)

    myDict = {'u':up}
    completeName = os.path.join('E:/dwell2/49/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
    
#Parallel(n_jobs=10)(delayed(lni)(cishu=j) for j in np.arange(500,2000,1)) 

def mean_T(cishu):
    en=0.4
    salt=cishu+random.random()*1000
    sys=make_system(width, length, str(salt))
    t=kwant.smatrix(sys,en).transmission(1,0)
    return t

#condu=Parallel(n_jobs=10)(delayed(mean_T)(cishu=j) for j in np.arange(500,1500,1))
##np.mean(condu[1:30])
#myDict = {'condu':condu}
#completeName = 'E:/dwell2/47.mat'
#sio.savemat(completeName,myDict,oned_as='row') 

#energies=np.linspace(0.2,0.5,500)
#sys=make_system(width, length, 'e57d')
#def spec(en):
#    t=kwant.smatrix(sys,en).transmission(1,0)
#    return t
#spectrum=Parallel(n_jobs=10)(delayed(spec)(en=j) for j in energies)
#myDict = {'spec':spectrum,'en':energies}
#completeName = 'E:/dwell2/45.mat'
#sio.savemat(completeName,myDict,oned_as='row') 
#
#pyplot.plot(energies,spectrum)

elapsed=time()-t_ini