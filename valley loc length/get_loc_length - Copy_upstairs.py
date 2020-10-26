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

file_path='E:/correlated_loc/7/'
energies = np.linspace(0, 3, 21)


#en=[]
xi_matrix=[]
for salt in range(5):
    xi_inverse_list=[]
    for i in range(21):
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
plt.plot(energies,np.transpose(np.array(xi_matrix)))
plt.figure()
plt.plot(energies,np.transpose(np.array(xi_matrix).mean(0)))

temp7=np.transpose(np.array(xi_matrix).mean(0))
#plt.plot(energies,xi,'.')
#plt.yscale('log', basey=2)
#plt.plot(energies,temp1)
#plt.plot(energies,temp2,'.')
#plt.plot(energies,temp3,'o')
#plt.plot(energies,temp4,'s')
#plt.plot(energies,temp5,'g')