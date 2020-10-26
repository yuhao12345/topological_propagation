# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 15:21:12 2018

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
#import random
import kwant
from kwant.digest import uniform
import time
import scipy.io as sio
import os.path

width=49
length=159
shape_whole=[]
shape_tau1=[]
v_dist=[]
t_list=[]
#final_tau_sqrt=[]
save_path = 'E:/index mismatch/16/'
for cishu in np.arange(0,500):#range(1000):
    completeName = os.path.join(save_path, str(cishu)+".mat")
#    myDict = {'u':u,'v':v,'s':s,'wf':wf,'T':T,'gf_mode':gf_mode}

    wf=sio.loadmat(completeName)['wf']
    wf_ave=np.mean(abs(wf)**2,axis=0)
    wf_ave_reshape=wf_ave.reshape((length,width))
    wf_shape=np.sum(wf_ave_reshape,axis=1)
    shape_whole.append(wf_shape)
    s=sio.loadmat(completeName)['s']
    if s[0,0]>0.98:
        v=sio.loadmat(completeName)['v']
        f_tau1=abs(np.array(sum(np.dot(np.matrix(np.diag(v[0,:])).H,wf))))**2
        f_tau1_reshape=f_tau1.reshape((length,width))
        shape_tau1.append(np.sum(f_tau1_reshape,axis=1))
#        shape.append(shape_tau1)
##        plt.figure()
##        plt.imshow(ff_tau1)
        v_dist.append(abs(v[0,:])**2)
        
    t=sio.loadmat(completeName)['T']
    t_list.append(t)

plt.title('ave_whole eigenchannel')
plt.plot(np.mean(np.array(shape_whole),axis=0))
plt.figure()
plt.plot(np.mean(np.array(shape_tau1),axis=0))
plt.title('ave_tau1')
plt.figure()
v_mean=np.mean(np.array(v_dist),axis=0)
plt.plot(v_mean)
plt.title('ave_v_tau1')
print(np.mean(t_list))
#%%
