# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 20:48:33 2019

@author: ykang
"""
import kwant
import stru_QSH as dg
#import stru_QSH_loss as dg_loss
#import stru_QSH_periodic_cavity as dg_p
import numpy as np
#import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from matplotlib import pyplot

t_ini=time()

width=50
length=100
dis=1.6

#num_freq=500
#ens=np.linspace(.38,.4,num_freq) #(0.3,0.4,10000)   # ens 


en=0.35
df=1e-9

#sys=dg.make_system_all(width, length, '62',dis)
##wf=kwant.wave_function(sys,0.35)(0)   #,check_hermiticity=False
#dg.plot_wf(sys,en,0)

def wf_stat(salt):
    sys=dg.make_system_all(width, length, str(salt), dis)
    return dg.wf_integral(sys,en,0)
        
#wf=Parallel(n_jobs=8)(delayed(wf_stat)(cishu) for cishu in np.arange(0,200,1))

def t_stat(salt):
    sys=dg.make_system_all(width, length, str(salt), dis)
    t_gf=dg.t_gf(sys,en,df,1,0) 
    return t_gf[0]
#    myDict = {'t':t_gf} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
#    completeName = os.path.join('E:/dwell4/145/', str(salt)+".mat")
#    sio.savemat(completeName,myDict,oned_as='row') 
    
#t_list=Parallel(n_jobs=8)(delayed(t_stat)(cishu) for cishu in np.arange(0,200,1))
    
#sys0=dg.make_system_all(width, length, str(0), 0)
##wf=np.abs(kwant.wave_function(sys0,en)(0)[0])**2
#coord=np.array([sys0.pos(i) for i in range(46200)])
#completeName0 = os.path.join('E:/dwell4/147/coord.mat')
#
#sio.savemat(completeName0,{'coord':coord},oned_as='row') 

def speckle_stat(salt):
    sys=dg.make_system_all(width, length, str(salt), dis)
    wf=np.abs(kwant.wave_function(sys,en)(0)[0])**2
    myDict = {'wf':wf} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
    completeName = os.path.join('E:/dwell4/147/', str(salt)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
Parallel(n_jobs=6)(delayed(speckle_stat)(cishu) for cishu in np.arange(200,1000,1))
#myDict = {'wf':wf} #,'ld':ld 't':t,  'trans':trans 'g':g t':t_list
#completeName = os.path.join('E:/dwell4/146/70_wf.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini