# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 11:32:26 2020

@author: ykang
"""

import stru_2d_normal_lossvariation as dn3
import numpy as np
#import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from matplotlib import pyplot

t_ini=time()

width=10
length=30
dis=2.2
salt='10555'
df=1e-9

num_loss=201
loss=np.linspace(-0.02,0.01,num_loss)

num_freq=1024
ens=np.linspace(0.6,0.73,num_freq)   # ens 


def dettm(cishu):
    sys=dn3.make_system_all(length,width,dis,salt,loss[cishu])
    t=[np.linalg.det(dn3.tranmission_DIY(sys,en)) for en in ens]
    myDict = {'t':t} 
    completeName = os.path.join('E:/channel time/79phase/'+str(cishu)+'.mat')  #7_tracet
    sio.savemat(completeName,myDict,oned_as='row') 
Parallel(n_jobs=5)(delayed(dettm)(cishu) for cishu in np.arange(0,num_loss,1))

elapsed=time()-t_ini