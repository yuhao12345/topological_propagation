# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 20:48:33 2019

@author: ykang
"""

import stru_QSH as dg
import numpy as np
import kwant
#from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
#from matplotlib import pyplot

t_ini=time()

width=30
length=200
dis=0

en=0.35
df=1e-9

#ens=np.linspace(-0.4,0.4,256)
#sys=dg.make_system_all(width, length, '55',dis)
##wf=kwant.wave_function(sys,0.35)(0)   #,check_hermiticity=False
#dg.plot_wf(sys,0.4)
#t=kwant.smatrix(sys,en).transmission(1,0)
    
def t_gf00(sys,en,df):     # only point to point
    g1=kwant.greens_function(sys,en).submatrix(1,0)[0,0]
    g2=kwant.greens_function(sys,en+df).submatrix(1,0)[0,0]
    return (np.angle(g2)-np.angle(g1))/df

def stat(salt):
    sys=dg.make_system_all(width, length, str(salt), dis)
#    t=dg.t_gf(sys,en,df,1,0)[0]*2*np.pi
#    wf=dg.wf_integral(sys,en,0)*2*np.pi
    gf0=t_gf00(sys,en,df)
#    return [wf,gf0]
    return gf0

    
def spec(e1):
    wf=kwant.wave_function(sys,e1)(0)[0]
    return sum(abs(wf)**2)
#    return kwant.greens_function(sys,e1).submatrix(1,0)[0,0]
#    return t_gf00(sys,e1,df)
#    return [kwant.greens_function(sys,e1).submatrix(1,0)[0,0],\
#            kwant.smatrix(sys,e1).transmission(1,0)]

#tspec= Parallel(n_jobs=8)(delayed(spec)(cishu) for cishu in ens)

#def stat_gf0_spec(salt):
#    sys=dg.make_system_all(width, length, str(salt), dis)
#    t=[t_gf00(sys,en,df) for en in ens]
#    completeName = os.path.join('E:/dwell5/50/', str(salt)+".mat")
#    sio.savemat(completeName,{'t':t},oned_as='row') 

def speckle_stat(salt):
    sys=dg.make_system_all(width, length, str(salt), dis)
    wf=kwant.wave_function(sys,en)(0)[0]
    myDict = {'wf':wf} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
    completeName = os.path.join('E:/dwell5/100/', str(salt)+".mat")
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=5)(delayed(speckle_stat)(cishu) for cishu in np.arange(10,1000,1))

sys0=dg.make_system_all(width, length, str(0), 0)
wf=np.abs(kwant.wave_function(sys0,en)(0)[0])**2
coord=np.array([sys0.pos(i) for i in range(wf.size)])
completeName0 = os.path.join('C:/Users/yumao/OneDrive - The University of Chicago/502.4.mat')
sio.savemat(completeName0,{'wf':wf,'coord':coord},oned_as='row') 

#st=Parallel(n_jobs=5)(delayed(stat)(cishu) for cishu in np.arange(500,1500))
#s=Parallel(n_jobs=8)(delayed(spec)(cishu) for cishu in ens)


#myDict = {'tspec':tspec,'ens':ens}
##myDict = {'st':st} #,'ld':ld 't':t,  'trans':trans 'g':g t':t_list  'st':st
#completeName = os.path.join('E:/dwell5/118.mat')
#sio.savemat(completeName,myDict,oned_as='row') 
#c1=os.path.join('C:/Users/ykang/Dropbox/TI8/data/118.mat')
#sio.savemat(c1,myDict,oned_as='row') 

elapsed=time()-t_ini