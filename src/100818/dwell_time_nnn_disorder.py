# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 17:11:39 2019

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

width=15
length=50
ratio=0.9

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]],norbs=2)  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b  

def make_system(width, length, salt):
    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width
    def nnn(site1, site2):
        x,y=site1.pos
        if abs(x)<length*ratio and y>width-5:
            return 1j * m2 * (1-2*uniform(repr(site1),salt)*1) * sz
        else:
            return 1j *m2 * sz

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= np.zeros([2,2]) 
    sys[graphene.neighbors()]=s0      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn #1j *m2 * sz #

    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym,conservation_law=-sz) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, 0))] = s0 # onsite# 
    lead[graphene.neighbors()]=s0
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 * sz
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
    return sys.finalized()

en=0.1
sys=make_system(width, length, 'ed')
wff=kwant.wave_function(sys,en)
wf=wff(0)

kwant.plotter.map(sys,abs(wf[12][0::2])**2)
    
#f0=0.4
#f1=0.45
#df=1e-11   # It is freq not omega !!!
#f_points=20
#energies=np.linspace(f0,f1,f_points)
#energies1=energies+df
#
## from TM
#def uv(energies,sys):  # get u,v from TM or TM_gf, TM_gf is right, but TM does not work
#    u_mat=[]
#    v_mat=[]
#    u2_mat=[]
#    v2_mat=[]
#    for en in energies:
#        gf=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!
#        gf1=gf[:,0::2]
#        u, s, vh=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
#        u_mat.append(u[:,0])
#        v_mat.append(vh.conj().T[:,0])
#        u2_mat.append(u[:,1])
#        v2_mat.append(vh.conj().T[:,1])
#    return u_mat,v_mat ,u2_mat,v2_mat
#
#
#def delay_tm(cishu):
#    salt=cishu+random.random()*100
#    sys=make_system(width, length, str(salt))
#    trans=[]
#    for en in energies:
#        trans.append(kwant.smatrix(sys,en).transmission((1, 0), (0, 0)))
#    u_mat,v_mat,u2_mat,v2_mat=uv(energies,sys)  
#    u_mat1,v_mat1,u2_mat1,v2_mat1=uv(energies1,sys)   
#    u0=np.column_stack(u_mat)
#    v0=np.column_stack(v_mat)
#    u1=np.column_stack(u_mat1)
#    v1=np.column_stack(v_mat1)
#    dphi_tm=np.diag((u0.conj().T @ (u1-u0)-v0.conj().T @ (v1-v0)))/1j/df/(2*np.pi)
#    return np.append(dphi_tm,trans)
#
#dwell=Parallel(n_jobs=20)(delayed(delay_tm)(cishu=j) for j in np.arange(0,2000,1)) 
#myDict = {'dwell':dwell}
#completeName = 'E:/dwell2/0/1.mat'
#sio.savemat(completeName,myDict,oned_as='row') 
#
#elapsed=time()-t_ini