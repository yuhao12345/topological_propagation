# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 16:39:05 2018

@author: ykang
"""

# loc length from transfer matrix, conductance g

import matplotlib.pyplot as plt
import scipy.io as sio
import os
from numpy.linalg import inv,svd
import numpy as np

file_path='E:/haldane_correlated/g6/'
energies = np.linspace(0.01, 3, 10)


#en=[]
xi_matrix=[]
for salt in np.arange(0,5):
    xi_inverse_list=[]
    for i in range(10):
        completeName = os.path.join(file_path+str(salt)+'/', str(i)+".mat")
        temp=sio.loadmat(completeName)
        gf=temp['gf']
        length=temp['length']
#        energy=temp['energy']
        width=temp['width'][0,0]
        # if diag number is 10^(-205), square will make it become 0 !!!
        
    #    xi_inverse=-np.log(np.sum(abs(np.diag(gf))**2))/2/length
        abs_diag_gf=abs(np.diag(gf))
        max_d=max(abs_diag_gf)
        log_trace_d=2*np.log(max_d)+np.log(sum((abs_diag_gf/max_d)**2))
        xi_inverse=-log_trace_d/2/length
    
        xi_inverse_list.append(xi_inverse)
#        en.append(energy)

    xi=np.array(xi_inverse_list)[:,0,0]**(-1)/width
    xi_matrix.append(xi)
#energies=np.array(en)[:,0,0]
plt.plot(energies,np.transpose(np.array(xi_matrix)),'s')
plt.figure()
plt.plot(energies,np.nanmean(np.array(xi_matrix),axis=0),'s')

tempg6=np.transpose(np.nanmean(np.array(xi_matrix),axis=0))

plt.figure(figsize=(20,10))
plt.plot(energies,tempg6,'r.')
plt.plot(energies,tempg10,'ro')
plt.plot(energies,tempg16,'rs')
plt.plot(energies,tempp6,'b.')
plt.plot(energies,tempp10,'bo')
plt.plot(energies,tempp16,'bs')
plt.yscale('log', basey=2)
plt.ylabel('loc length')
plt.xlabel('energy')
matplotlib.rc('font', size=5)