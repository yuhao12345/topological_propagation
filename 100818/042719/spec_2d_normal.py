# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 10:53:48 2019

@author: ykang
"""

import stru_QSH as dg
import stru_2d_normal as dn
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
from matplotlib import pyplot

t_ini=time()

width=5
#width_lead=3
length=70
dis=2.2
salt='1'

sys=dn.make_system_all(length,width,dis,salt)
#sys=dn.make_system_crys_air_crys(length,width,width_lead,dis,salt)
#wf=kwant.wave_function(sys,0.42)(0)
#ld=sum(kwant.ldos(sys,0.42))
#tm=dn.tranmission_DIY(sys,0.5)
num_freq=500
ens=np.linspace(0.42,0.5,num_freq)   # ens 

#e1=0.39
df=1e-9

#gf=kwant.greens_function(sys,.38).submatrix(1,0)
#u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
#
#g=kwant.smatrix(sys,.55).submatrix(1,0)

#u1, s1, vh1=np.linalg.svd(g, full_matrices=True, compute_uv=True)

#t=dg.t_mode_all(sys,0.42,df) #t10,t01,t00,t11
#wf=dg.wf_integral_all(sys,e1)
#dg.plot_wf(sys,.33,0)
#kwant.smatrix(sys,ens[181]).transmission(1,0)

def spec_time(cishu):
    return dg.t_gf_all(sys,ens[cishu],df)   #t10,t01,t00,t11

#t=Parallel(n_jobs=5)(delayed(spec_time)(cishu) for cishu in np.arange(0,num_freq,1))
#
def spec_wf(cishu):
    return dg.wf_integral_all(sys,ens[cishu])
#
#wf=Parallel(n_jobs=5)(delayed(spec_wf)(cishu) for cishu in np.arange(0,num_freq,1))
#
def spec_trans(cishu):
    return kwant.smatrix(sys,ens[cishu]).transmission(1,0)
#    return dg.trans_all(sys,ens[cishu])   #t10,t01,t00,t11

#trans=Parallel(n_jobs=8)(delayed(spec_trans)(cishu) for cishu in np.arange(0,num_freq,1))
#pyplot.plot(ens,trans)

def spec_eigentrans(cishu):  # for fist eigenchannel
    s0=np.linalg.svd(kwant.smatrix(sys,.38).submatrix(1,0),full_matrices=False,compute_uv=False)
    return s0[0]**2
#trans_eigen=Parallel(n_jobs=8)(delayed(spec_eigentrans)(cishu) for cishu in np.arange(0,num_freq,1))

def spec_gf(cishu):
    return kwant.greens_function(sys,ens[cishu]).submatrix(1,0)[6,3]
##g=Parallel(n_jobs=5)(delayed(spec_gf)(cishu) for cishu in np.arange(0,num_freq,1))


def channeltime(cishu):
    t=dn.tranmission_DIY(sys,ens[cishu])
    t1=dn.tranmission_DIY(sys,ens[cishu]+df)
    
    return np.array([np.trace(t),np.trace(t1)])

trace_t=Parallel(n_jobs=5)(delayed(channeltime)(cishu) for cishu in np.arange(0,num_freq,1))
#
def DOS(cishu):
#    wf=kwant.wave_function(sys,ens[cishu])(0)
#    left=np.sum(np.abs(wf)**2,1)
#    wf=kwant.wave_function(sys,ens[cishu])(1)
#    right=np.sum(np.abs(wf)**2,1)
#    return left+right  
    return sum(kwant.ldos(sys,ens[cishu]))
#dos=Parallel(n_jobs=5)(delayed(DOS)(cishu) for cishu in np.arange(0,num_freq,1))
myDict = {'en':ens,'trace_t':trace_t} #,'ld':ld ,'dos':dos   'wf':wf, 'trace_t':trace_t
completeName = os.path.join('E:/channel time/17_tracet.mat')  #7_tracet
sio.savemat(completeName,myDict,oned_as='row') 

#completeName = os.path.join('E:/channel time/1/'+str(cishu)+'.mat')
#sio.savemat(completeName,myDict,oned_as='row')

###sys.pos(sys.lead_interfaces[1][0])
    

elapsed=time()-t_ini