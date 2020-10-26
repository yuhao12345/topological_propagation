# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 16:39:05 2018

@author: ykang
"""

# loc length from transfer matrix, conductance g
# One-Parameter Scaling of Localization Length and Conductance in Disordered Systems

import matplotlib.pyplot as plt
import scipy.io as sio
import os
from numpy.linalg import inv,svd
import numpy as np

file_path='C:/Users/ykang/Documents/correlated_loc/33/'
xi_inverse_list=[]
en=[]
for i in range(31):
    completeName = os.path.join(file_path, str(i)+".mat")
    temp=sio.loadmat(completeName)
    gf=temp['gf']
    length=temp['length']
    energy=temp['energy']
    width=temp['width'][0,0]
    # if diag number is 10^(-205), square will make it become 0 !!!
    
#    xi_inverse=-np.log(np.sum(abs(np.diag(gf))**2))/2/length
    abs_diag_gf=abs(np.diag(gf))
    max_d=max(abs_diag_gf)
    log_trace_d=2*np.log(max_d)+np.log(sum((abs_diag_gf/max_d)**2))
    xi_inverse=-log_trace_d/2/length

    xi_inverse_list.append(xi_inverse)
    en.append(energy)
    
#    s00=temp['s00']
#    s01=temp['s01']
#    s10=temp['s10']
#    s11=temp['s11']
#    t00=s10-s11@inv(s01)@s00
#    t01=s11@inv(s01)
#    t10=-inv(s01)@s00
#    t11=inv(s01)
#    transfer=np.bmat([[t00, t01], [t10, t11]])
#        
#    u, s, vh=svd(transfer)
#    s_list.append(s)
#    s_list(k,:)=s;
#    t_list(k)=t;

#s_list=np.array(s_list)
#xi_inverse=np.log(s_list[:,0])
xi=np.array(xi_inverse_list)[:,0,0]**(-1)/width

energies=np.array(en)[:,0,0]

plt.plot(energies,xi,'.')
#plt.yscale('log', basey=10)

#plt.imshow(pattern)
