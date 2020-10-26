# -*- coding: utf-8 -*-
"""
Yuhao Kang  rgf G_whole space and ldos
"""
import numpy as np
from numpy import (pi,sqrt,sin,cos,arange,diag,squeeze,eye,imag,
                   floor,arccos,ones,zeros,concatenate,meshgrid)
from numpy.linalg import inv
import matplotlib.pyplot as plt
from matplotlib import cm
import cmath
from mpl_toolkits.mplot3d import Axes3D

# define system size and mesh
k=1
d=1
tiny=0
length, width=200, 100
zushu=1
kd=k*d
m,n =int(length/d+1), int(width/d-1) 
mode_max=int(floor(arccos(1-kd**2/2)*(n+1)/pi))  #calculate the useful mode number
h0=4*diag(ones(n))-diag(ones(n-1),1)-diag(ones(n-1),-1)  #define the hopping
wg=ones((n,m-2))

re_g=kd**2/2-2+cos((arange(mode_max)+1)*pi/(n+1))
g_mode=re_g-1j*sqrt(1-re_g**2+0j) #gf in mode space 

phi=zeros((n,mode_max),dtype=np.complex)
gl_nn=zeros((m,n,n),dtype=np.complex)
gl_n0=zeros((m,n,n),dtype=np.complex)
gl_0n=zeros((m,n,n),dtype=np.complex)
gr_nn=zeros((m,n,n),dtype=np.complex) 
gr_n1=zeros((m,n,n),dtype=np.complex) 
gr_1n=zeros((m,n,n),dtype=np.complex) 
g_n0=zeros((m,n,n),dtype=np.complex)
g_n1=zeros((m,n,n),dtype=np.complex)
g_nn=zeros((m,n,n),dtype=np.complex)
g_0n=zeros((m,n,n),dtype=np.complex)
g_1n=zeros((m,n,n),dtype=np.complex)
for y in range(mode_max):    #convert real space to mode space
    phi[:,y]=sqrt(2/(n+1))*sin(pi*(y+1)*(arange(n)+1)/(n+1))

g0=phi.dot(diag(g_mode)).dot(phi.T)  #real space nearest hopping in lead    from mode space to real space
gl_nn[0,:,:]=g0   
gl_n0[0,:,:]=g0
gl_0n[0,:,:]=g0
gr_nn[m-1,:,:]=g0
gr_n1[m-1,:,:]=g0
gr_1n[0,:,:]=g0

for conf in range(zushu):
    eps=concatenate((ones((n,1)),wg+np.random.random((n,m-2))*tiny,ones((n,1))),axis=1)  #define disorder system

    for x in range(2,m):                #iteration
        temp=inv(diag(eps[:,x-1])*kd**2-h0-squeeze(gl_nn[x-2,:,:]))
        gl_nn[x-1,:,:]=temp  
        gl_n0[x-1,:,:]=-temp.dot(squeeze(gl_n0[x-2,:,:]))
        gl_0n[x-1,:,:]=-squeeze(gl_0n[x-2,:,:]).dot(temp)
    
    for x in range(m-1,1,-1):
        temp=inv(diag(eps[:,x-1])*kd**2-h0-squeeze(gr_nn[x,:,:]))
        gr_nn[x-1,:,:]=temp;
        gr_n1[x-1,:,:]=-temp.dot(squeeze(gr_n1[x,:,:]))
        gr_1n[x-1,:,:]=-squeeze(gr_1n[x,:,:]).dot(temp)
    
    for x in range(2,m):
        glnn=squeeze(gl_nn[x-1,:,:])
        grnn=squeeze(gr_nn[x,:,:])
        gln0=squeeze(gl_n0[x-1,:,:])
        gg=eye(n)-glnn.dot(grnn)
        g_n0[x-1,:,:]=inv(gg).dot(gln0)
        temp=inv(gg).dot(glnn)
        g_nn[x-1,:,:]=temp
        g_n1[x-1,:,:]=-temp.dot(squeeze(gr_n1[x,:,:]))
        g_0n[x-1,:,:]=squeeze(gl_0n[x-1,:,:])+squeeze(gl_0n[x-1,:,:]).dot(squeeze(gr_nn[x,:,:])).dot(squeeze(g_nn[x-1,:,:]))
        g_1n[x-1,:,:]=-squeeze(gr_1n[x,:,:]).dot(squeeze(g_nn[x-1,:,:]))
    
    s11=zeros((n,m-2),dtype=np.complex)          #return probability
    for x in range(2,m):
        s11[:,x-2]=diag(squeeze(g_nn[x-1,:,:]))
    
    g=squeeze(g_n0[:,:,5])
#    plt.imshow(abs(g)**2)
    xv,yv=meshgrid(arange(n),arange(m))
    xv1,yv1=meshgrid(arange(n),arange(m-2))
    zv=abs(g)**2
    plt.rcParams["figure.figsize"] = [20,15]
    #plot speckle pattern
    
    surf = plt.figure().gca(projection='3d').plot_surface(xv, yv, zv, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    #plot LDOS
    ldos=-imag(s11).T/pi
    surf = plt.figure().gca(projection='3d').plot_surface(xv1, yv1, ldos,
                      cmap=cm.coolwarm, linewidth=0, antialiased=False)




