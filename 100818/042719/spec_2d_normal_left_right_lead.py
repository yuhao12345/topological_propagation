# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:37:44 2019

@author: ykang
"""

import stru_QSH as dg
import stru_2d_normal_left_right_lead as dnlr
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
from matplotlib import pyplot

t_ini=time()

width=6
length=50
width_left=6
width_right=6
dis=1.3
salt='10'

sys=dnlr.make_system_all(length,width,dis,salt,width_left,width_right)
#kwant.plot(sys,fig_size=(15, 10))


num_freq=2000
ens=np.linspace(0.45,0.46,num_freq)   # ens 
df=1e-9

flead0 = sys.leads[0]
#flead1 = sys.leads[1]
energy=0.4

# t_gf==t_mode2
t_mode=dnlr.t_mode(sys,energy,df,1,0)
t_mode2=dnlr.t_mode_DIY(sys,energy,df,1,0)
t_gf=dnlr.t_gf(sys,energy,df,1,0)

t_s=dnlr.tranmission_DIY(sys,energy)
tau=np.linalg.svd(t_s, full_matrices=True, compute_uv=False) 
trans2=np.sum(tau**2)

trans=kwant.smatrix(sys,energy).transmission(1,0)

#wf_inte=sum(dnlr.wf_integral_all(sys,energy))/2

G = kwant.greens_function(sys, energy)
Gr=G.submatrix(1,0) # r indicates retarded
l=np.linalg.svd(Gr, full_matrices=True, compute_uv=False)
#
#tr=kwant.smatrix(sys, energy).submatrix(1,0)
#l1=np.linalg.svd(tr, full_matrices=True, compute_uv=False)
#
##Ga=np.conj((Gr).T)        # a indicates advanced
#
#Sigma_Lr =flead0.selfenergy(energy)  #r indicates retarded
#Sigma_La =np.conj((Sigma_Lr).T)           #a indicates advanced
#
#Sigma_Rr =flead1.selfenergy(energy)
#Sigma_Ra =np.conj((Sigma_Rr).T)
#
#Gamma_in =1j*(Sigma_Lr-Sigma_La)
#Gamma_out=1j*(Sigma_Rr-Sigma_Ra)
#
##### TM=Gamma_out*Gr*Gamma_in*Ga
##TM=np.dot(np.dot(np.dot(Gamma_out,Gr),Gamma_in),Ga) #Fisher-Lee relation
##trans=np.trace(TM)
#
### from gf to TM
##g10=kwant.greens_function(sys,energy).submatrix(1,0)
##flead0 = sys.leads[0]
prop_modes, _ = flead0.modes(energy)
v=prop_modes.velocities#*2*np.pi  # 0: left 1:right
n=v.size//2   # number of channel
##
#wf_lead=(prop_modes.wave_functions[:,1])
#wf_lead_n=np.sum(np.abs(wf_lead)**2)
#w=np.matrix(wf_lead)
#gamma=np.conj(w).T*v@w/wf_lead_n   # correct with kwant result
##lead_vel=np.sqrt(v)*wf_lead/np.sqrt(wf_lead_n)
##
#t_s2=1j*np.sqrt(v)*np.sqrt(v)*w @ Gr @ np.conj(w).T /wf_lead_n   # trnamsission matrix

#prop_modes1, _ = flead1.modes(energy)
#v1=prop_modes1.velocities[1]#*2*np.pi  # 0: left 1:right
#wf_lead1=(prop_modes1.wave_functions[:,0])
#wf_lead1_n=np.sum(np.abs(wf_lead1)**2)
#w1=np.matrix(wf_lead1)
#lead_vel1=np.sqrt(v1)*wf_lead1/np.sqrt(wf_lead1_n)
#t_s2=1j*np.sqrt(v)*np.sqrt(v)*w @ Gr @ np.conj(w1).T /wf_lead_n
#t=1j*np.conj(lead_vel1)@(g10)@(lead_vel.T)
##
#s=kwant.smatrix(sys,energy)
#tt10=s.data[1,0]
#t10=s.transmission(1,0)
#
#wf=kwant.wave_function(sys,0.45)(0)
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

#t=dnlr.t_mode_all(sys,.45,df)
#t2=dnlr.t_gf_all(sys,.45,df)

def spec_trans(cishu):
    return kwant.smatrix(sys,ens[cishu]).transmission(1,0)
#trans=Parallel(n_jobs=8)(delayed(spec_trans)(cishu) for cishu in np.arange(0,num_freq,1))

def spec_gf(cishu):
    gf=kwant.greens_function(sys,ens[cishu]).submatrix(1,0) #,check_hermiticity=False
    myDict = {'en':ens,'gf':gf} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
    completeName = os.path.join('E:/dwell4/9/',str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=8)(delayed(spec_gf)(cishu) for cishu in np.arange(0,num_freq,1))
    
#pyplot.plot(trans,'.')

def spec_time(cishu):
    t_mode=dnlr.t_mode_all(sys,ens[cishu],df)   #t10,t01,t00,t11
    myDict = {'en':ens,'t_mode':t_mode} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
    completeName = os.path.join('E:/dwell4/5/',str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=8)(delayed(spec_time)(cishu) for cishu in np.arange(0,num_freq,1))
    
def spec_wf(cishu):
    return dnlr.wf_integral_all(sys,ens[cishu])
#wf=Parallel(n_jobs=8)(delayed(spec_wf)(cishu) for cishu in np.arange(0,num_freq,1))

#myDict = {'en':ens,'trans':trans} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
#completeName = os.path.join('E:/dwell4/6.5.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#num_freq=200
#ens=np.linspace(0.42,0.45,num_freq)   # ens 
#
##e1=0.39
#df=1e-9
#
#gf=kwant.greens_function(sys,.38).submatrix(1,0)
#u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
#
#g=kwant.smatrix(sys,.38).submatrix(1,0)
#u1, s1, vh1=np.linalg.svd(g, full_matrices=True, compute_uv=True)

#t=dg.t_mode_all(sys,e1,df) #t10,t01,t00,t11
#wf=dg.wf_integral_all(sys,e1)
#dg.plot_wf(sys,0.4513,0)
#kwant.smatrix(sys,ens[181]).transmission(1,0)

#wf=kwant.wave_function(sys,0.4513)(0)   #,check_hermiticity=False
#coord=np.array([sys.pos(i) for i in range(7581)])
#myDict = {'wf':wf,'coord':coord} #,'ld':ld
#completeName = os.path.join('C:/Users/ykang/Desktop/TI8/F4_1_wf_0.4513.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini

#energy=0.45


