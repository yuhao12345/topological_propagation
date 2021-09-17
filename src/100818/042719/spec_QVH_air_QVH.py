# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:35:02 2019

@author: ykang
"""

##### trivial system, QVH with a air gap, keep single channel
import stru_QVH_air_QVH as qa
import stru_2d_normal as dn
import stru_QSH as dg
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
t_ini=time()

width=40.7
width_lead=3
length=100
dis=1
salt='5'

sys=qa.make_system_all(width, width_lead, length, salt, dis) 

wf=kwant.wave_function(sys,0.2)(0)

ens=np.linspace(0.3,0.4,200)   # ens 

e1=0.39
df=1e-9

#t=dg.t_mode_all(sys,0.25,df) #t10,t01,t00,t11
#wf=dg.wf_integral_all(sys,0.3)
#
#dg.plot_wf(sys,ens[15],0)
#print(kwant.smatrix(sys,ens[15]).transmission(1,0))

def spec_time(cishu):
    return dg.t_mode_all(sys,ens[cishu],df)   #t10,t01,t00,t11

t=Parallel(n_jobs=5)(delayed(spec_time)(cishu) for cishu in np.arange(0,200,1))

def spec_wf(cishu):
    return dg.wf_integral_all(sys,ens[cishu])

wf=Parallel(n_jobs=5)(delayed(spec_wf)(cishu) for cishu in np.arange(0,200,1))

def spec_trans(cishu):
    return dg.trans_all(sys,ens[cishu])   #t10,t01,t00,t11

trans=Parallel(n_jobs=5)(delayed(spec_trans)(cishu) for cishu in np.arange(0,200,1))

def spec_gf(cishu):
    return kwant.greens_function(sys,ens[cishu]).submatrix(1,0)[7,3]
g=Parallel(n_jobs=5)(delayed(spec_gf)(cishu) for cishu in np.arange(0,200,1))
#
myDict = {'en':ens,'t':t, 'wf':wf, 'trans':trans,'g':g} #,'ld':ld , 
completeName = os.path.join('E:/dwell3/829/17.mat')
sio.savemat(completeName,myDict,oned_as='row') 

###sys.pos(sys.lead_interfaces[1][0])
#[sys.pos(sys.lead_interfaces[0][i]) for i in range(10)]
elapsed=time()-t_ini