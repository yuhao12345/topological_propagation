# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 21:55:52 2019

@author: ykang
"""

import stru_QSH as dg
#import stru_QSH_loss as dg_loss
#import stru_QSH_periodic_cavity as dg_p
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from matplotlib import pyplot

t_ini=time()

width=50
length=60
dis=2

num_freq=500
ens=np.linspace(.38,.4,num_freq) #(0.3,0.4,10000)   # ens 

def wf_stat(salt):
    sys=dg.make_system_all(width, length, str(salt), dis)
    wf=[]
    for cishu in np.arange(0,num_freq,1):
        wf.append(dg.wf_integral(sys,ens[cishu],0))
    
    myDict = {'en':ens,'wf':wf} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
    completeName = os.path.join('E:/dwell4/51/', str(salt)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
#    del sys
Parallel(n_jobs=5)(delayed(wf_stat)(cishu) for cishu in np.arange(505,1000,1))

elapsed=time()-t_ini