# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 11:39:35 2019

@author: ykang
"""

##### many configuration statistics
import stru_QSH as dg
import numpy as np
#import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
t_ini=time()

width=60
length=50
dis=1.6

e1=0.39
df=1e-9

def stat_time(cishu):   # return time, wf_integral
    
    sys=dg.make_system_all(width, length, str(cishu),dis)   # (width, length, salt, dis)
    t = dg.t_gf_all(sys,e1,df)
    wf = dg.wf_integral_all(sys,e1)
    return np.concatenate([t,wf])
#stat_t=Parallel(n_jobs=5)(delayed(stat_time)(cishu) for cishu in np.arange(0,1000,1))
#
#myDict = {'stat_t':stat_t} #,'ld':ld
#completeName = os.path.join('E:/dwell3/828/10.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#sys=dg.make_system_all(width, length, '1',dis)   # (width, length, salt, dis)
#dg.plot_wf(sys,0.390677354709419,1)
