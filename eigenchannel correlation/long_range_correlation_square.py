# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 16:42:05 2018

@author: ykang
"""

# long range correlation based on fourier transfomation, 2d system

from scipy.special import kv,gamma
from numpy import pi,exp,sqrt
import numpy as np
from numpy import random
#import matplotlib.pyplot as plt
#from time import time

#t_ini=time()


#t1=2**8  # size of system
#t2=2**14



# power law  correlation(q)~q^(-v)



def f_sq(x,y,v):  
#    v=.1
    p=v/2-1  # beta
    q=(x**2+y**2)**0.5
    return 2*pi/gamma(p+1)*((q/2)**p)*kv(p,q)

def correlation_pattern(t1,t2,v,w,salt):
    # power law  correlation(q)~q^(-v)
    # w: root mean square, control the magnitude of disorder  
    l1=(t1-1)/2
    l2=(t2-1)/2
#    salt=256
    random.seed(salt)
    phase=random.rand(t1,t2)*2*pi
    fk=np.zeros((t1,t2))
    for i in range(t1):
        for j in range(t2):
            fk[i,j]=f_sq(2*pi*(i-l1)/t1/sqrt(3),2*pi*(j-l2)/t2,v)   # unit vector in k space: (1/sqrt(3),1)
    
    fp=sqrt(fk)*exp(1j*phase)  
    h=np.real(np.fft.ifft2(fp))
    #sqrt(np.mean(h**2))
    return h/sqrt(np.mean(h**2))*w   # set root mean square as a const

# guassian   ~exp(-l^2/2/v^2)   prb 94,134203(2016)
#v=.3    # correlation length
def gua_f(x,y,v):
    return exp(-(x**2+y**2)/v**2)/pi/v**2
def correlation_pattern_guassian(t1,t2,v,w,salt):
    random.seed(salt)
    l1=(t1-1)/2
    l2=(t2-1)/2
    fx=np.zeros((t1,t2))
    for i in range(t1):
        for j in range(t2):
            fx[i,j]=gua_f(i-l1,j-l2,v)
    fk=np.fft.fft2(fx)
    phase=random.rand(t1,t2)*2*pi
    fp=fk*exp(1j*phase)
    h=np.real(np.fft.ifft2(fp))
    return h/sqrt(np.mean(h**2))*w/sqrt(12)
#    return h/np.max(abs(h))*w

def uniformrandom(t1,t2,w,salt):
    np.random.seed(salt)
    return (np.random.rand(t1,t2)-0.5)*w
#plt.imshow(tt)    
#np.mean(h)
#np.max(h)
#print([np.max(tt), np.min(tt), np.std(tt)])

#elapsed=time()-t_ini