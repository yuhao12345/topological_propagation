# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 20:45:34 2020

@author: user
"""

import stru_2d_normal2 as dn2
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from matplotlib import pyplot

t_ini=time()

width=10
#width_lead=3
length=30
dis=2.2
salt='10555'
df=1e-9
sys=dn2.make_system_all(length,width,dis,salt)

c=5
en1=2-2*np.cos(c*np.pi/(2*width))  # lower
en2=2-2*np.cos((c+1)*np.pi/(2*width))  # higher

#tm=dn.tranmission_DIY(sys,0.5)
num_freq=1024
ens=np.linspace(0.6,0.73,num_freq)   # ens 

ha=sys.hamiltonian_submatrix()    # h_in
flead0 = sys.leads[0]
flead1 = sys.leads[1]

def spec(cishu):
#    t=dn2.tranmission_DIY(sys,ens[cishu])
    t=kwant.smatrix(sys,ens[cishu],check_hermiticity=False).submatrix(1,0)
    myDict = {'t':t} 
    completeName = os.path.join('E:/channel time/79/'+str(cishu)+'.mat')  #7_tracet
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=5)(delayed(spec)(cishu) for cishu in np.arange(0,num_freq,1))   

def DOS(cishu):  # ONLY FOR UNITARY SYSTEM
    return sum(kwant.ldos(sys,ens[cishu]))

def invG(cishu):
    en=ens[cishu]
    Sigma_Lr =flead0.selfenergy(en)  #r indicates retarded    
    Sigma_Rr =flead1.selfenergy(en)
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
#dos_invG=Parallel(n_jobs=5)(delayed(invG)(cishu) for cishu in np.arange(0,num_freq,1))

def dettime(cishu):
    t1=dn2.tranmission_DIY(sys,ens[cishu])
    t2=dn2.tranmission_DIY(sys,ens[cishu]+df)
    myDict = {'t1':t1,'t2':t2,'df':df} 
    completeName = os.path.join('E:/channel time/79/'+str(cishu)+'.mat')  #7_tracet
    sio.savemat(completeName,myDict,oned_as='row') 
#    return np.array([t1,t2])
#Parallel(n_jobs=5)(delayed(dettime)(cishu) for cishu in np.arange(0,num_freq,1))
#
#myDict = {'en':ens,'dos_invG':dos_invG,} #,'ld':ld ,'dos':dos   'wf':wf, 'trace_t':trace_t
#completeName = os.path.join('E:/channel time/79_dos_dett_0_2.mat')  #7_tracet
#sio.savemat(completeName,myDict,oned_as='row') 

def dettm(cishu):
    t=dn2.tranmission_DIY(sys,ens[cishu])
    de=np.linalg.det(t)
    return de #[np.angle(de),abs(de)]
dett=Parallel(n_jobs=5)(delayed(dettm)(cishu) for cishu in np.arange(0,num_freq,1))
   
myDict = {'en':ens,'dett':dett} #,'ld':ld ,'dos':dos   'wf':wf, 'trace_t':trace_t
completeName = os.path.join('E:/channel time/79phase/0.mat')  #7_tracet
sio.savemat(completeName,myDict,oned_as='row') 
 
#en=0.68728
#t=dn2.tranmission_DIY(sys,en)
#wf=kwant.wave_function(sys,en,check_hermiticity=False)(0)   #
#coord=np.array([sys.pos(i) for i in range(1121)])
#myDict = {'t':t,'wf':wf,'coord':coord} #,'ld':ld
#completeName = os.path.join('C:/Users/ykang/Dropbox/taut//wf_0.68728_loss0.mat')
#sio.savemat(completeName,myDict,oned_as='row') 
    
elapsed=time()-t_ini