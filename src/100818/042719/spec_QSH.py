# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:07:13 2019

@author: ykang
"""
### one configuration, obtain spectrum
import stru_QSH as dg
import stru_2d_normal_left_right_lead as dnlr
#import stru_QSH_loss as dg_loss
#import stru_QSH_periodic_cavity as dg_p
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
from matplotlib import pyplot

# use fisher lee relation to obtain smatrix is hard here, it is ok in trivial system
t_ini=time()

width=30
length=100
salt='2'
dis=2

sys=dg.make_system_all(width, length, salt,dis)
#sys=dg_loss.make_system_all(width, length, salt,dis)   # (width, length, salt, dis)
#sys=dg_p.make_system_all_periodic_cavity(30,100)
#kwant.plot(sys)

#l=dg.make_lead(width)
#kwant.plot(l)
num_freq=512
ens=np.linspace(.392,.398,num_freq) #(0.3,0.4,10000)   # ens 

#e1=0.39
df=1e-9

#t_gf=dg.t_gf(sys,.3935,df,1,0)

#t_gf10=dg.t_gf(sys,.398,df,1,0) 

#gf=kwant.greens_function(sys,ens[100]).submatrix(1,0)[:,0]

#t=dg.t_gf_all(sys,e1,df)  #t10,t01,t00,t11
#wf=dg.wf_integral_all(sys,e1)
#wf0=dg.wf_integral(sys,e1,0)
#wf1=dg.wf_integral(sys,e1,1)

#wf=kwant.wave_function(sys,0.35)(0)   #,check_hermiticity=False
#dg.plot_wf(sys,ens[21],0)
#dg.plot_wf(sys,0.4,0)   #0.39649189

#s=kwant.smatrix(sys,0.39992172)
#sdata=s.data
#t10=kwant.smatrix(sys,0.35).transmission(1,0)
sm=kwant.smatrix(sys,0.35).data
tm=kwant.smatrix(sys,0.35).submatrix(1,0)
#t01=s.transmission(0,1)
#t_gf=dg.t_gf(sys,0.35,1e-9,1,0)
#t_gf2=dg.t_gf(sys,0.39649243,1e-10,0,0)
#t_gf3=dg.t_gf(sys,0.39649243,1e-9,0,0)

#myDict = {'t':t_gf} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
#completeName = os.path.join('E:/dwell4/133/', str(5)+".mat")
#sio.savemat(completeName,myDict,oned_as='row') 
    
# save wf
#coord=np.array([sys.pos(i) for i in range(27600)])
#myDict = {'wf':wf,'coord':coord} #,'ld':ld
#completeName = os.path.join('C:/Users/ykang/Desktop/TI8/F2_4_wf_0.39649243.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#energy=0.35
#
##t_r=dnlr.tranmission_DIY(sys,energy)  # DIY transmission matrix does not work
#flead0 = sys.leads[0]
#flead1 = sys.leads[1]
#
#G = kwant.greens_function(sys, energy)
#Gr=G.submatrix(1,0) # r indicates retarded
#l=np.linalg.svd(Gr, full_matrices=True, compute_uv=False)
#Ga=np.conj((Gr).T)        # a indicates advanced
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
#### TM=Gamma_out*Gr*Gamma_in*Ga
#TM=np.dot(np.dot(np.dot(Gamma_out,Gr),Gamma_in),Ga) #Fisher-Lee relation
#trans=np.trace(TM)
#
#g10=kwant.greens_function(sys,energy).submatrix(1,0)
#prop_modes, _ = flead0.modes(energy)
#v=prop_modes.velocities[1]#*2*np.pi  # 0: left 1:right
#
#wf_lead=(prop_modes.wave_functions[:,0])
#
#cood_lead=np.array([flead0.pos(i) for i in range(21)])
##pyplot.plot(cood_lead[:,1],np.abs(wf_lead),'.')
##pyplot.plot(cood_lead[:,1],np.angle(wf_lead),'.')
#wf_lead_n=np.sum(np.abs(wf_lead)**2)
#lead_vel=np.sqrt(v)*wf_lead/np.sqrt(wf_lead_n)
#w=np.matrix(wf_lead)
#gamma=np.conj(w).T*v@(w)/wf_lead_n   # correct with kwant result
##
#prop_modes1, _ = flead1.modes(energy) 
#v1=prop_modes1.velocities[1]#*2*np.pi  # 0: left 1:right
#wf_lead1=(prop_modes1.wave_functions[:,0])
#wf_lead1_n=np.sum(np.abs(wf_lead1)**2)
#lead_vel1=np.sqrt(v1)*wf_lead1/np.sqrt(wf_lead1_n)
##
#t=1j*np.conj(lead_vel1)@(g10)@(lead_vel.T)
#
#s=kwant.smatrix(sys,energy)
#tt10=s.data[1,0]
#t10=s.transmission(1,0)

def spec_time(cishu):
#    return dg.t_gf_all(sys,ens[cishu],df)   #t10,t01,t00,t11
    t_gf10=dg.t_gf(sys,ens[cishu],df,1,0)   #t10,t01,t00,t11
    t_gf01=dg.t_gf(sys,ens[cishu],df,0,1) 
    myDict = {'en':ens,'t_gf10':t_gf10,'t_gf01':t_gf01} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
    completeName = os.path.join('E:/dwell4/141/',str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 

#Parallel(n_jobs=6)(delayed(spec_time)(cishu) for cishu in np.arange(0,num_freq,1))

def spec_wf(cishu):
#    return dg.wf_integral(sys,ens[cishu],0)
    return dg.wf_integral_all(sys,ens[cishu])

#wf=Parallel(n_jobs=8)(delayed(spec_wf)(cishu) for cishu in np.arange(0,num_freq,1))

def spec_trans(cishu):
    return kwant.smatrix(sys,ens[cishu]).transmission(1,0)
#    return dg.trans_all(sys,ens[cishu])   #t10,t01,t00,t11

#trans=Parallel(n_jobs=8)(delayed(spec_trans)(cishu) for cishu in np.arange(0,num_freq,1))

def spec_gf(cishu):   
#    return kwant.greens_function(sys,ens[cishu]).submatrix(1,0)[0,0] #,check_hermiticity=False
    gf=kwant.greens_function(sys,ens[cishu]).submatrix(1,0)[:,0]
    myDict = {'en':ens,'gf':gf} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
    completeName = os.path.join('E:/dwell4/123/',str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row')     
#Parallel(n_jobs=4)(delayed(spec_gf)(cishu) for cishu in np.arange(0,num_freq,1))

def spec_t(cishu):   
#    return kwant.greens_function(sys,ens[cishu]).submatrix(1,0)[0,0] #,check_hermiticity=False
    t=kwant.smatrix(sys,ens[cishu]).data
    myDict = {'en':ens,'t':t} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
    completeName = os.path.join('E:/dwell4/124/',str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=4)(delayed(spec_t)(cishu) for cishu in np.arange(0,num_freq,1))
    
    
#myDict = {'en':ens, 'wf':wf} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
#completeName = os.path.join('E:/dwell4/140_wf.mat')
#sio.savemat(completeName,myDict,oned_as='row') 


def wf_en(cishu):
    wf=kwant.wave_function(sys,ens[cishu])(0)[0]
    wf_I=np.abs(wf)**2
    myDict = {'wf_I':wf_I} #,'ld':ld ,  'g':g 't':t, 'wf':wf, 'trans':trans
    completeName = os.path.join('E:/dwell4/126/',str(cishu)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=8)(delayed(wf_en)(cishu) for cishu in np.arange(0,num_freq,1))
#v=dg.ballistic_velocity(width,.4)

#pyplot.plot(wf,'.')

elapsed=time()-t_ini

#sys.pos(sys.lead_interfaces[0][0])

#flead0 = sys.leads[0]
##flead1 = sys.leads[1]
#energy=0.35
#
#t_mode=dnlr.t_mode(sys,energy,df,1,0)
#t_gf=dnlr.t_gf(sys,energy,df,1,0)
#
#G = kwant.greens_function(sys, energy)
#Gr=G.submatrix(1,0) # r indicates retarded
#l=np.linalg.svd(Gr, full_matrices=True, compute_uv=False)
#
#tr=kwant.smatrix(sys, energy).submatrix(1,0)
#l1=np.linalg.svd(tr, full_matrices=True, compute_uv=False)
##Ga=np.conj((Gr).T)        # a indicates advanced
#prop_modes, _ = flead0.modes(energy)
#v=prop_modes.velocities[1]#*2*np.pi  # 0: left 1:right
