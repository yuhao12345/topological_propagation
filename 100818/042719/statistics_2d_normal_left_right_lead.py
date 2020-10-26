# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 14:49:40 2019

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

width=30
length=80
width_left=width
width_right=width
#width_left=30
#width_right=30
dis=1.2

en=1.4

syst=dnlr.make_system_all(length,width,dis,'10',width_left,width_right)
sys=syst.finalized()
wf_t=kwant.wave_function(sys,en)(0)
coord=np.array([sys.pos(i) for i in range(wf_t[0].size)])
#kwant.plot(sys,fig_size=(15, 10))

myDict = {'coord':coord} 
completeName = os.path.join('E:/pt/20/coord.mat')
sio.savemat(completeName,myDict,oned_as='row') 

#trans=kwant.smatrix(sys,en).transmission(1,0)
#wf=kwant.wave_function(sys,en)(0)
#t_s=dnlr.tranmission_DIY(sys,en)
#tau=np.linalg.svd(t_s, full_matrices=True, compute_uv=False) 
#trans2=np.sum(tau**2)

def trans_stat(salt):
    syst=dnlr.make_system_all(length,width,dis,salt,width_left,width_right)
    sys=syst.finalized()
    return kwant.smatrix(sys,en).transmission(1,0)

def wf_stat(salt):
    syst=dnlr.make_system_all(length,width,dis,str(salt),width_left,width_right)
    sys=syst.finalized()
    t_s=dnlr.tranmission_DIY(sys,en)
    wf=kwant.wave_function(sys,en)(0)
    myDict = {'TM':t_s,'wf':wf} 
    completeName = os.path.join('E:/pt/20/', str(salt)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 

trans=Parallel(n_jobs=5)(delayed(wf_stat)(cishu) for cishu in np.arange(10,100,1))

elapsed=time()-t_ini