# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:49:17 2018

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import kwant
from kwant.digest import uniform
import time
import scipy.io as sio
import os.path

start=time.time()
width=30
length=150
n1=1
n2=1
n=1
sys=kwant.Builder()
lat=kwant.lattice.square()


en=1
def disk(pos):
    x,y=pos
    return 0<=x<=length and 0<y<width
def edge(pos):
    x,y=pos
    return x==0 and 0<y<width
def onsite(site):
    return 4/(n+.15*(np.sign(uniform(repr(site),salt)-0.5)+1))**2

sys[lat.shape(disk,(1,1))]=onsite
sys[lat.neighbors()]=-1/n**2
        
lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0))) 
lead[lat.shape(edge,(0,1))]=4/n1**2
lead[lat.neighbors()]=-1/n1**2
sys.attach_lead(lead)
sys.attach_lead(lead.reversed())
sys=sys.finalized()
#shape=[]
#final_T=[]
#final_tau_sqrt=[]
#save_path = 'E:/kwant/26/'
ans=np.zeros((200,10))
ans1=np.zeros((200,10))
ind=0



for cishu in range(200):

    ans[ind,0]=length
    ans1[ind,0]=length
    salt=str(cishu+random.random())

    gf_mode=kwant.smatrix(sys,en).submatrix(1,0)  
    try:
        u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
        wf=kwant.wave_function(sys,en)(0)
        temp=[]
        for mode in range(9):
            f_tau1=abs(np.array(sum(np.dot(np.matrix(np.diag(v[mode,:])).H,wf))))**2
            shape_tau1=np.sum(f_tau1.reshape((length+1,width-1)),axis=1)
            ans[ind,mode+1]=np.argmax(shape_tau1)/length
            ans1[ind,mode+1]=s[mode]
#        ans=np.row_stack([ans,np.array(hang)])
#        T=kwant.smatrix(sys,en).transmission(1,0)
#        completeName = os.path.join(save_path, str(cishu)+".mat")
#        myDict = {'u':u,'v':v,'s':s,'wf':wf,'T':T,'gf_mode':gf_mode}
#    ##    myDict['u'] = u
#        sio.savemat(completeName,myDict,oned_as='row')
    except:
        continue
    ind=ind+1
            
    #['E:\\kwant\\1\\',str(cishu),'.mat']
#    final_T.append(kwant.smatrix(sys,en).transmission(1,0))
#    final_tau_sqrt.append(s)
#    a_test=abs(wf[0,:].reshape((length+1,width-1)))**2
#    print(s)

#    if s[0]>0.997:
#        f_tau1=abs(np.array(sum(np.dot(np.matrix(np.diag(v[0,:])).H,wf))))**2
#        ff_tau1=f_tau1.reshape((length+1,width-1))
#        plt.figure()
#        plt.imshow(ff_tau1)
#
#        shape_tau1=np.sum(ff_tau1,axis=1)
#        shape.append(shape_tau1)
#
#
#    
#    kwant.plotter.map(sys,(abs(np.array(wf[0,:]))**2))
##    local_dos = kwant.ldos(sys, en)
##    kwant.plotter.map(sys, local_dos, num_lead_cells=0)
##point_pos=np.array([sys.pos(i) for i in range(sys.graph.num_nodes)])
##temp=np.array([abs(np.array(wf[2,:]))**2])
##field=np.concatenate((point_pos,temp.T),axis=1)
#plt.figure()
#plt.plot(np.mean(np.array(shape),axis=0))
import pandas as pd 
df = pd.DataFrame(ans)
df.to_csv("l150.csv")
df1 = pd.DataFrame(ans1)
df1.to_csv("lt150.csv")
end=time.time()
print(end-start)
#print(np.mean(final_T))
#plt.imshow