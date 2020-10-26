# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 20:48:56 2018

@author: user
"""
#https://kwant-project.org/doc/1/reference/generated/kwant.physics.StabilizedModes#kwant.physics.StabilizedModes
#all the propagating modes are normalized to carry unit current

import kwant
import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform
import scipy.io as sio
import os.path
from joblib import Parallel, delayed
import multiprocessing

t_ini=time()

length=20  # toi calculate loc length, sample length=1e4
width=30

dis=2

lat = kwant.lattice.square()

def make_system(length,width,salt):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    #    return abs(x) + abs(y)*np.sqrt(3) < 100
    #    return abs(x)+10>abs(y)*3 and abs(x)<100
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width #2
    def onsite(site):
        x,y=site.pos
#        return dis*(uniform(repr(site),salt)-.5)+4
        if y>width-10:
            return 4
        else:
            return dis*(uniform(repr(site),salt)-.5)+4
    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = -1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 4
    lead[lat.neighbors()] = -1

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()


def loc_length(cishu):    # loc length from green function 
    en=0.4
    salt=cishu+random()
    sys=make_system(length,width, str(salt))
    gf=kwant.greens_function(sys,en).submatrix(1,0)
    abs_diag_gf=abs(np.diag(gf))
    max_d=max(abs_diag_gf)
    log_trace_d=2*np.log(max_d)+np.log(sum((abs_diag_gf/max_d)**2))
    xi_inverse=-log_trace_d/2/length
    xi=1/xi_inverse
    return xi
#loclength1=Parallel(n_jobs=20)(delayed(loc_length)(cishu=j) for j in np.arange(0,40,1)) 
    
#f0=0.4
#f1=0.45
#df=1e-10   # It is freq not omega !!!
#f_points=3
#energies=np.linspace(f0,f1,f_points)
#energies1=energies+df
##
#### from TM
#def uv(energies):  # get u,v from TM or TM_gf, TM_gf is right, but TM does not work
#    u_mat=[]
#    v_mat=[]
#    s_mat=[]
#    u2_mat=[]
#    v2_mat=[]
#    u3_mat=[]
#    v3_mat=[]
#    for en in energies:
##        tm=kwant.smatrix(sys,en)
##        gf_mode=tm.submatrix((1, 0), (0, 0)) 
##        gf_mode1=tm.submatrix((1, 1), (0, 0)) 
##        gf1=np.concatenate((gf_mode, gf_mode1))
##        gf1=kwant.smatrix(sys,en).submatrix(1,0)
#        gf1=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!
#       
#        u, s, vh=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
#        u_mat.append(u[:,0])
#        v_mat.append(vh.conj().T[:,0])
#        u2_mat.append(u[:,1])
#        v2_mat.append(vh.conj().T[:,1])
#        u3_mat.append(u[:,2])
#        v3_mat.append(vh.conj().T[:,2])
#        s_mat.append(s)
#    return u_mat,v_mat,u2_mat,v2_mat,u3_mat,v3_mat,s_mat#, 
#u_mat,v_mat,u2_mat,v2_mat,u3_mat,v3_mat,s_mat=uv(energies)   #
#u_mat1,v_mat1,u2_mat1,v2_mat1,u3_mat1,v3_mat1,s_mat1=uv(energies1)#
##
#u0=np.column_stack(u_mat)
#v0=np.column_stack(v_mat)
#u1=np.column_stack(u_mat1)
#v1=np.column_stack(v_mat1)
##
#u20=np.column_stack(u2_mat)   # second eigenstate
#v20=np.column_stack(v2_mat)
#u21=np.column_stack(u2_mat1)
#v21=np.column_stack(v2_mat1)
#
#u30=np.column_stack(u3_mat)   # second eigenstate
#v30=np.column_stack(v3_mat)
#u31=np.column_stack(u3_mat1)
#v31=np.column_stack(v3_mat1)
#
#dphi_tm=np.diag((u0.conj().T @ (u1-u0)-v0.conj().T @ (v1-v0)))/1j/df/(2*np.pi)
#dphi2_tm=np.diag((u20.conj().T @ (u21-u20)-v20.conj().T @ (v21-v20)))/1j/df/(2*np.pi)
#dphi3_tm=np.diag((u30.conj().T @ (u31-u30)-v30.conj().T @ (v31-v30)))/1j/df/(2*np.pi)

## first calculate the mean of log(I) for all channels, then do average over the cross section
#def lnt_x(cishu):
#    global salt
#    salt=str(cishu+random())
#    sys=make_system()
#    
##    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
##    t=kwant.smatrix(sys,en).transmission(1,0) 
#    wff=kwant.wave_function(sys,en)
#    wf=wff(0)    
#    wf1=np.mean(np.log(np.abs(wf)**2),axis=0)
#    ans=[]
#    for k in range(len(wf1)):
#        ans.append(np.append(sys.pos(k),wf1[k]))
#    myDict = {'ans':ans}
#    completeName = os.path.join('C:/Users/ykang/Documents/yuhao_share/23/', str(cishu)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 
    
# first do integral of cross section for each channel, then do average for all channels in all samples
en=0.4
def lnt_x1(cishu):
    salt=cishu+random()
    sys=make_system(length,width, str(salt))
    
#    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#    t=kwant.smatrix(sys,en).transmission(1,0) 
    wff=kwant.wave_function(sys,en)(0)
    w=np.abs(wff.T)
    coord=np.array([sys.pos(i) for i in range(w.shape[0])])
    u=np.concatenate((coord,w),axis=1)

    myDict = {'u':u}
    completeName = os.path.join('E:/dwell2/14/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 

#Parallel(n_jobs=10)(delayed(lnt_x1)(cishu=j) for j in np.arange(0,2000,1)) 
        
        
# transmission relative to smaple length
def trans(cishu):
    salt=str(cishu+random())
    sys=make_system(length,width, str(salt))  
    t=kwant.smatrix(sys,en).transmission(1,0)  
    return t

#tra=Parallel(n_jobs=10)(delayed(trans)(cishu=j) for j in np.arange(0,500,1)) 
#myDict = {'tra':tra}
#completeName = os.path.join('E:/dwell2/t7.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

def tau_n(cishu):
    en=0.4
    salt=cishu+random.random()*1000
    sys=make_system(length,width,str(salt))
    gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
#    s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False) #u, s, vh
    
    wf=kwant.wave_function(sys,en)(0)
#    kwant.plotter.map(sys, (abs(wf[14])**2),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

    ldos = kwant.ldos(sys,en)
    if cishu==0:
        coord=np.array([sys.pos(i) for i in range(ldos.shape[0])])
        sio.savemat('E:/dwell3/167/coord.mat', {'coord':coord},oned_as='row') 
    myDict = {'gf_mode':gf_mode, 'wf':wf, 'ld':ldos}
    completeName = os.path.join('E:/dwell3/167/', str(cishu)+".mat")
    sio.savemat(completeName,myDict,oned_as='row')    
#    return s

#Parallel(n_jobs=8)(delayed(tau_n)(cishu=j) for j in np.arange(0,200,1)) 

#s=Parallel(n_jobs=12)(delayed(tau_n)(cishu=j) for j in np.arange(0,1000,1)) 
#myDict = {'s':s}
#completeName = 'E:/dwell3/210.mat'
#sio.savemat(completeName,myDict,oned_as='row') 

sys=make_system(length,width,'5')
gf_mode=kwant.smatrix(sys,en).submatrix(1,0) 
gf_mode1=kwant.smatrix(sys,en).submatrix(0,1) 

elapsed=time()-t_ini
#
#sys=make_system(length,width, 'ed')  
#en=0.4
#t=kwant.smatrix(sys,en).transmission(1,0)  
#wff=kwant.wave_function(sys,en)
#wf=wff(0)
#kwant.plotter.map(sys, abs(wf[0]),num_lead_cells=5,fig_size=(15, 10))  
#
#plt.figure()
##plt.imshow(abs(wf[0]).reshape(39,9))
#tt=abs(wf[0]).reshape(39,9)
#plt.plot(tt[0,])
#sum(tt[0,]**2)
#y=[abs(np.sin(3*np.pi/10*i))*np.sqrt(2/14.5) for i in np.arange(1,10)]
#plt.plot(y)

