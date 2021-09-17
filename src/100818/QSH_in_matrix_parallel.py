# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 21:50:27 2018

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

width=15
whole_length=100
ratio=100/whole_length  #disorder/whole length
length=100/ratio

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]],norbs=2)  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2
#Lambda_R = 2

s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b


def onsite(site):
    return np.zeros([2,2])     

#ll=1  # region with rashba or nnn disorder
      # measure dwell time, ll=0.9 


    
#dis=.7  #disordser stength of nnn hopping
def make_system(width, length, salt):

    def disk(pos):
        x,y=pos
        return abs(y)<width and abs(x)<length   #25.1
    def lead_shape(pos):
        x,y=pos
        return abs(y)<width #width-5<y<width    #
    
    def nnn(site1, site2):
        x,y=site1.pos
        if abs(x)<length*ratio and y>width-5:
#            return 1j * m2 * (uniform(repr(site1),salt)-0.5) * 2 * sz
            return 1j * m2 * (1-2*uniform(repr(site1),salt)*0) * sz
#            return 1j *m2 * sz
        else:
            return 1j *m2 * sz
        
    # center point of rashba region
#    xr= (np.divide(random.sample(range(100), 30),100)-0.5)*2*length*ratio
#    yr=width-np.divide(random.sample(range(100), 30),100)*5   #density!!!
#    def hop_ras_e1(site1,site2,Lambda_R=Lambda_R,t1=1):
#        x,y=site1.pos
#        if min(abs(x-xr)+abs(y-yr))<1.5:
#            return s0-1j*Lambda_R*sx
#        else:
#            return s0
#        
#    
#    def hop_ras_e2(site1,site2,Lambda_R=Lambda_R,t1=1):
#        x,y=site1.pos
#        if min(abs(x-xr)+abs(y-yr))<1.5:
#            return s0+1j*Lambda_R*(0.5*sx - np.sqrt(3)/2.0*sy)
#        else:
#            return s0
#        
#    
#    def hop_ras_e3(site1,site2,Lambda_R=Lambda_R,t1=1):
#        x,y=site1.pos
#        if min(abs(x-xr)+abs(y-yr))<1.5:
#            return s0+1j*Lambda_R*(0.5*sx + np.sqrt(3)/2.0*sy)
#        else:
#            return s0
    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]=onsite
    sys[graphene.neighbors()]=s0      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn #1j *m2 * sz #
#    sys[kwant.HoppingKind((0, 0), a,b)] = hop_ras_e1
#    sys[kwant.HoppingKind((0, 1), a,b)] = hop_ras_e2
#    sys[kwant.HoppingKind((-1, 1), a,b)] = hop_ras_e3

    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym,conservation_law=-sz) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, 0))] = onsite# s0 # 
    lead[graphene.neighbors()]=s0
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 * sz
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
    return sys.finalized()

#en=.5
#sys=make_system(width, length, '5')
#wff=kwant.wave_function(sys,en)
#wf=wff(0)
#
#kwant.plotter.map(sys,abs(wf[0][0::2])**2,num_lead_cells=5, fig_size=(20, 10))
#tm=kwant.smatrix(sys,.1)
#t10_uu=tm.transmission((1, 0), (0, 0))
#t10_du=tm.transmission((1, 1), (0, 0))
#t00_uu=tm.transmission((0, 0), (0, 0))
#t00_du=tm.transmission((0, 1), (0, 0))
#gf=kwant.greens_function(sys,.1).submatrix(1,0)
#smatrix.transmission((lead1, q1), (lead0, q0)) is the transmission from the
# q0 block of the lead0 into the q1 block of lead1.

#tm=kwant.smatrix(sys,0.11)
#gf_mode=tm.submatrix((1, 0), (0, 0)) 
#u, s, vh=np.linalg.svd(gf_mode, full_matrices=True, compute_uv=True)
    
f0=0.4
f1=0.45
df=1e-11   # It is freq not omega !!!
f_points=20
energies=np.linspace(f0,f1,f_points)
energies1=energies+df

# for only one freq
#def uv(energies,sys):  # get u,v from TM or TM_gf, TM_gf is right, but TM does not work
#    u_mat=[]
#    v_mat=[]
#    u2_mat=[]
#    v2_mat=[]
##    s_mat=[]
#    en= energies
#    gf=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!
#    gf1=gf[:,0::2]
#    u, s, vh=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
#    u_mat=(u[:,0])
#    v_mat=(vh.conj().T[:,0])
#    u2_mat=(u[:,1])
#    v2_mat=(vh.conj().T[:,1])
#    s_mat=(s)
#    return u_mat,v_mat ,u2_mat,v2_mat,s_mat

#def delay2(cishu):
#    salt=cishu+random.random()*100
#    sys=make_system2(width, length, str(salt))
#    tm=kwant.smatrix(sys,f0)
#    t10_uu=tm.transmission((1, 0), (0, 0))
#    t10_du=tm.transmission((1, 1), (0, 0))
##    t00_uu=tm.transmission((0, 0), (0, 0))
##    t00_du=tm.transmission((0, 1), (0, 0))
#    u_mat,v_mat,u2_mat,v2_mat,s_mat=uv(energies,sys)  
#    u_mat1,v_mat1,u2_mat1,v2_mat1,s_mat1=uv(energies1,sys)   
#    u0=(u_mat)
#    v0=(v_mat)
#    u1=(u_mat1)
#    v1=(v_mat1)
#    
#    u20=(u2_mat)   # second eigenstate
#    v20=(v2_mat)
#    u21=(u2_mat1)
#    v21=(v2_mat1)
#    dphi_tm=((u0.conj().T @ (u1-u0)-v0.conj().T @ (v1-v0)))/1j/df/(2*np.pi)
#    dphi2_tm=((u20.conj().T @ (u21-u20)-v20.conj().T @ (v21-v20)))/1j/df/(2*np.pi)
#    return np.append(np.array([t10_uu,t10_du,dphi_tm,dphi2_tm]),s_mat[0:3])


# from TM
def uv(energies,sys):  # get u,v from TM or TM_gf, TM_gf is right, but TM does not work
    u_mat=[]
    v_mat=[]
    u2_mat=[]
    v2_mat=[]
#    s_mat=[]
    for en in energies:
        gf=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!
        gf1=gf[:,0::2]
        u, s, vh=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
        u_mat.append(u[:,0])
        v_mat.append(vh.conj().T[:,0])
        u2_mat.append(u[:,1])
        v2_mat.append(vh.conj().T[:,1])
#        s_mat.append(s)
    return u_mat,v_mat ,u2_mat,v2_mat #,s_mat






## from greens function
#def phase_gf(energies,sys):
#    gf=[]
#    for en in energies:
#        gf10=kwant.greens_function(sys,en).submatrix(1,0)
#        gf10_uu=gf10[0,0]
#        gf.append(gf10_uu)
#    return gf

def delay_tm(cishu):
    salt=cishu+random.random()*100
    sys=make_system(width, length, str(salt))
    trans=[]
    for en in energies:
        trans.append(kwant.smatrix(sys,en).transmission((1, 0), (0, 0)))
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
    return np.append(np.append(dphi_tm,trans),salt)

## transmission of incident spin up
def trans_spin_up(cishu):
    en=0.5
    salt=cishu+random.random()*100
    sys=make_system2(width, length, str(salt))  
    tm=kwant.smatrix(sys,en)
    t10_uu=tm.transmission((1, 0), (0, 0))
#    t10_du=tm.transmission((1, 1), (0, 0))
#    t00_uu=tm.transmission((0, 0), (0, 0))
#    t00_du=tm.transmission((0, 1), (0, 0))
    return np.array([t10_uu,salt])   #t10_du,t00_uu,t00_du,

#def delay_gf(cishu):      # from gf
#    salt=cishu+random.random()*100
#    sys=make_system(width, length, str(salt))
#    gf0=phase_gf(energies, sys)
#    gf1=phase_gf(energies1, sys)    
#    dphi_gf=(np.angle(gf1)-np.angle(gf0))/df/(2*np.pi)
#    trans=[]
#    for en in energies:
#        trans.append(kwant.smatrix(sys,en).transmission((1, 0), (0, 0)))
#    return np.append(np.append(dphi_gf,trans),salt)
    
def ldos(cishu):
    en=0.4
    salt=cishu+random.random()*100
    sys=make_system2(width, length, str(salt))
    local_dos = kwant.ldos(sys, en)
    ld = local_dos.reshape(-1, 2)

    coord=np.array([sys.pos(i) for i in range(ld.shape[0])])
    ans=np.concatenate((coord,ld),axis=1)
    myDict = {'ans':ans}
    completeName = os.path.join('E:/dwell/11/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 

def lni(cishu):
    en=0.5
    salt=cishu+random.random()*100
    sys=make_system(width, length, str(salt))
    wff=kwant.wave_function(sys,en)(0)
    w=np.abs(wff.T)
    coord=np.array([sys.pos(i) for i in range(w.shape[0]//2)])
    up=np.concatenate((coord,w[0::2,]),axis=1)
    down=np.concatenate((coord,w[1::2,]),axis=1)
    myDict = {'u':up,'d':down}
    completeName = os.path.join('E:/dwell2/5/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
    
#energies=np.linspace(0.33,0.51,500)
#salt='ed'
#sys=make_system2(width, length, salt)
#def con(n):    
#    en=energies[n]
#    gf=kwant.greens_function(sys,en).submatrix(1,0)
#    t=kwant.smatrix(sys,en).transmission((1, 0), (0, 0))
#    myDict = {'gf':gf,'t':t}
#    completeName = os.path.join('E:/dwell/72/', str(n)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 

#spec = Parallel(n_jobs=20)(delayed(con)(n) for n in range(500))

#num_cores = multiprocessing.cpu_count()
#dwell=[]
#dwell=Parallel(n_jobs=20)(delayed(delay_gf)(cishu=j) for j in np.arange(0,500,1)) 
#condu=Parallel(n_jobs=num_cores-2)(delayed(trans_spin_up)(cishu=j) for j in np.arange(1000,5000,1)) 
#condu=Parallel(n_jobs=20)(delayed(delay2)(cishu=j) for j in np.arange(0,100,1)) 
#myDict = {'condu':condu}   # {'dwell':dwell} #
#completeName = 'E:/dwell/33.mat'
#sio.savemat(completeName,myDict,oned_as='row') 

Parallel(n_jobs=1)(delayed(lni)(cishu=j) for j in np.arange(0,4,1)) 

elapsed=time()-t_ini

#sys.pos(290)
#sys.lead_interfaces