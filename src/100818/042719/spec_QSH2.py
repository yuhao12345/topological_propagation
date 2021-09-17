# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 15:07:13 2019

@author: ykang
"""
### one configuration, obtain spectrum
import stru_QSH as dg
#import stru_2d_normal_left_right_lead as dnlr
#import stru_QSH_loss as dg_loss
#import stru_QSH_periodic_cavity as dg_p
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
from matplotlib import pyplot


#####  eigenchannel time was divided by 2pi in previous definition!!!!!
# use fisher lee relation to obtain smatrix is hard here, it is ok in trivial system
t_ini=time()

width=30
length=50
salt='25'
dis=3

sys=dg.make_system_all(width, length, salt,dis)
dg.plot_wf(sys,.395,0)
def t_wiger(sys,en,df):
    s1=kwant.smatrix(sys,en).data
    s2=kwant.smatrix(sys,en+df).data
#    tw=s2.transpose().conj()@(s2-s1)/df
#    return np.trace(tw)
    p1=np.angle(np.linalg.det(s1))
    p2=np.angle(np.linalg.det(s2))
    return (p2-p1)/df

num_freq=64
ens=np.linspace(.392,.398,num_freq) #(0.3,0.4,10000)   # ens 

#e1=0.39
df=1e-9

#def det_gf(sys,en,df):    # useless, det(gf) gives 0
#    g1=kwant.greens_function(sys,en).submatrix(1,0)
#    g2=kwant.greens_function(sys,en+df).submatrix(1,0)
#    d1=np.linalg.det(g1)
#    d2=np.linalg.det(g2)
#    return (np.angle(d2)-np.angle(d1))/df
#


def t_gf00(sys,en,df):     # only point to point
    g1=kwant.greens_function(sys,en).submatrix(1,0)[0,0]
    g2=kwant.greens_function(sys,en+df).submatrix(1,0)[0,0]
    return (np.angle(g2)-np.angle(g1))/df

def spec_wf(cishu):
    return dg.wf_integral(sys,ens[cishu],0)

#wf=Parallel(n_jobs=5)(delayed(spec_wf)(cishu) for cishu in np.arange(0,num_freq,1))



def t_gf(cishu):
    return dg.t_gf(sys,ens[cishu],df,1,0)[0]*2*np.pi

#tgf=Parallel(n_jobs=5)(delayed(t_gf)(cishu) for cishu in np.arange(0,num_freq,1))

def t_gf1(cishu):
    return t_gf00(sys,ens[cishu],df)
#tgf1=Parallel(n_jobs=5)(delayed(t_gf1)(cishu) for cishu in np.arange(0,num_freq,1))



def spec_trans(cishu):
    return kwant.smatrix(sys,ens[cishu]).transmission(1,0)
#    return dg.trans_all(sys,ens[cishu])   #t10,t01,t00,t11

#trans=Parallel(n_jobs=8)(delayed(spec_trans)(cishu) for cishu in np.arange(0,num_freq,1))

  
    
#myDict = {'en':ens, 'wf':wf,'tgf':tgf, 'tgf1':tgf1} #,'ld':ld 't':t, 'wf':wf, 'trans':trans 'g':g
#completeName = os.path.join('E:/dwell4/152.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

#sys.pos(sys.lead_interfaces[0][0])

elapsed=time()-t_ini
