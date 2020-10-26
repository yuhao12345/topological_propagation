# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:00:42 2020

@author: user
"""

#import stru_QSH as dg
import stru_2d_normal_left_right_lead as dnlr
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from numpy import linalg as LA
#from matplotlib import pyplot

t_ini=time()

width=20
length=30
width_left=20
width_right=width_left
dis=2.2
salt='1'

df=1e-9

sys=dnlr.make_system_all(length,width,dis,salt,width_left,width_right)
sys=sys.finalized()
#wf=kwant.wave_function(sys,0.3495,check_hermiticity=False)(0)
#kwant.plotter.map(sys, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)
#ld=sum(kwant.ldos(sys,0.42))
c=8
en1=2-2*np.cos(c*np.pi/(2*width_left))  # lower
en2=2-2*np.cos((c+1)*np.pi/(2*width_left))  # higher
#tm=dnlr.tranmission_DIY(sys,0.2)
num_freq=256
ens=np.linspace(0.39,0.46,num_freq)   # ens 

#tt=dnlr.t_mode_DIY(sys,0.5,df,1,0)
#
ha=sys.hamiltonian_submatrix()    # h_in
## fisher lee
#en=0.75
#G = kwant.greens_function(sys, en)
#Gr=G.submatrix(1,0) # r indicates retarded
#G.submatrix(1,0)-G.submatrix(0,1).T
#Ga=np.conj((Gr).T)        # a indicates advanced
#
flead0 = sys.leads[0]
flead1 = sys.leads[1]
#Sigma_Lr =flead0.selfenergy(en)  #r indicates retarded
#Sigma_La =np.conj((Sigma_Lr).T)           #a indicates advanced
#
#Sigma_Rr =flead1.selfenergy(en)
#Sigma_Ra =np.conj((Sigma_Rr).T)
#
#Gamma_in =1j*(Sigma_Lr-Sigma_La)
#Gamma_out=1j*(Sigma_Rr-Sigma_Ra)
#
#T=np.trace(Gr @ Gamma_in @ Ga @ Gamma_out)
#T_kwant=kwant.smatrix(sys,en).transmission(1,0)
#
#n=Gamma_in.shape[0]
#m=ha.shape[0]
#coupling = np.zeros(ha.shape, dtype=np.complex128)
#coupling[0:n,0:n]=coupling[0:n,0:n]+Sigma_Lr
#coupling[m-n:m,m-n:m]=coupling[m-n:m,m-n:m]+Sigma_Rr
#Heff=ha+coupling
#I = np.eye(ha.shape[0], ha.shape[1], dtype=np.complex128)*en
#g_hinv=np.linalg.inv(I-Heff)
#G10=g_hinv[m-n:m,0:n]   # same as kwant result
#
#w,v=LA.eig(Heff)
#
#tp=g_hinv*np.linalg.det(I-Heff)
#tpp=np.trace(Gr)*np.linalg.det(I-Heff)
#a=I-Heff
##completeName = os.path.join('D:/tp/0.mat')  #7_tracet
##sio.savemat(completeName,{'a':a},oned_as='row') 
#
#Sigma_Lr+Sigma_La
#
#prop_modes, _ = flead0.modes(en)
#v=prop_modes.velocities#*2*np.pi  # 0: left 1:right
#n=v.size//2  # number of channel
#v_matrix_sqrt= np.diag([v[i]**0.5 for i in range(n,2*n)])
#wf_lead=prop_modes.wave_functions[:,n:2*n]
#wf_lead_n=np.sum(np.abs(wf_lead)**2,0)**0.5
#wf_lead_unit=wf_lead/wf_lead_n
#tppp=np.conj(wf_lead_unit).T@wf_lead_unit
##t_s=1j*v_matrix_sqrt @ (wf_lead_unit).T  @ Gr @ np.conj(wf_lead_unit)  @ v_matrix_sqrt 




#gf=kwant.greens_function(sys,.38).submatrix(1,0)
def gf(cishu):
    gf=kwant.greens_function(sys,ens[cishu]).submatrix(1,0)
    gf1=kwant.greens_function(sys,ens[cishu]+df).submatrix(1,0)
    myDict = {'en':ens,'gf':gf,'gf1':gf1} 
    completeName = os.path.join('E:/channel time/30/'+str(cishu)+'.mat')  #7_tracet
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=5)(delayed(gf)(cishu) for cishu in np.arange(0,num_freq,1))

def channeltime(cishu):
    t=dnlr.tranmission_DIY(sys,ens[cishu])
    t1=dnlr.tranmission_DIY(sys,ens[cishu]+df)
    return np.array([np.trace(t),np.trace(t1)])

#trace_t=Parallel(n_jobs=5)(delayed(channeltime)(cishu) for cishu in np.arange(0,num_freq,1))
#
def DOS(cishu):  # ONLY FOR UNITARY SYSTEM
    return sum(kwant.ldos(sys,ens[cishu]))

def invG(cishu):
    en=ens[cishu]
    Sigma_Lr =flead0.selfenergy(en)  #r indicates retarded
#    Sigma_La =np.conj((Sigma_Lr).T)           #a indicates advanced
    
    Sigma_Rr =flead1.selfenergy(en)
#    Sigma_Ra =np.conj((Sigma_Rr).T)
    
#    Gamma_in =1j*(Sigma_Lr-Sigma_La)
#    Gamma_out=1j*(Sigma_Rr-Sigma_Ra)

    n=Sigma_Lr.shape[0]
    m=ha.shape[0]
    coupling = np.zeros(ha.shape, dtype=np.complex128)
    coupling[0:n,0:n]=coupling[0:n,0:n]+Sigma_Lr
    coupling[m-n:m,m-n:m]=coupling[m-n:m,m-n:m]+Sigma_Rr
    Heff=ha+coupling
    I = np.eye(ha.shape[0], ha.shape[1], dtype=np.complex128)*en
    g_hinv=np.linalg.inv(I-Heff)
    dos=-sum(np.diag(g_hinv).imag/np.pi)
    return dos
#dos_invG=Parallel(n_jobs=3)(delayed(invG)(cishu) for cishu in np.arange(0,num_freq,1))

def inte_wf(cishu):   # =DOS when no loss
    return (np.sum(np.abs(kwant.wave_function(sys,ens[cishu],check_hermiticity=False)(0))**2)+\
    np.sum(np.abs(kwant.wave_function(sys,ens[cishu],check_hermiticity=False)(1))**2))/2/np.pi
#inteI=Parallel(n_jobs=5)(delayed(inte_wf)(cishu) for cishu in np.arange(0,num_freq,1))
def eigen(cishu):
    return np.real(sum(dnlr.t_mode_DIY(sys,ens[cishu],df,1,0)))
#eigtime=Parallel(n_jobs=5)(delayed(eigen)(cishu) for cishu in np.arange(0,num_freq,1))

def dettime(cishu):
    t1=np.linalg.det(dnlr.tranmission_DIY(sys,ens[cishu]))
    t2=np.linalg.det(dnlr.tranmission_DIY(sys,ens[cishu]+df))
    return np.array([t1,t2])
#dett=Parallel(n_jobs=5)(delayed(dettime)(cishu) for cishu in np.arange(0,num_freq,1))
#
#myDict = {'en':ens,'dos_invG':dos_invG,'dett':dett,'df':df} #,'ld':ld ,'dos':dos   'wf':wf, 'trace_t':trace_t
#completeName = os.path.join('E:/channel time/63_dos_dett_0.001.mat')  #7_tracet
#sio.savemat(completeName,myDict,oned_as='row') 

def spec(cishu):
#    t=dnlr.tranmission_DIY(sys,ens[cishu])
    t=kwant.smatrix(sys,ens[cishu]).submatrix(1,0)
    myDict = {'t':t} 
    completeName = os.path.join('E:/channel time/65/'+str(cishu)+'.mat')  #7_tracet
    sio.savemat(completeName,myDict,oned_as='row') 
Parallel(n_jobs=5)(delayed(spec)(cishu) for cishu in np.arange(0,num_freq,1))

###sys.pos(sys.lead_interfaces[1][0])

#### plot field profile
#t=kwant.smatrix(sys,0.73509286).submatrix(1,0)
#wf=kwant.wave_function(sys,0.73509286,check_hermiticity=False)(0)
#coord=np.array([sys.pos(i) for i in range(wf.shape[1])])
#myDict = {'wf':wf,'coord':coord,'t':t} 
#completeName = os.path.join('C:/Users/ykang/Dropbox/taut/Q1Dzero_0.73509286.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini