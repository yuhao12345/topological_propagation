# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:23:38 2019

@author: ykang
"""

import stru_1d as od
import numpy as np
import kwant
from matplotlib import pyplot
import scipy.io as sio
import os.path
from time import time
from joblib import Parallel, delayed
t_ini=time()

length=50
dis=1

f_num=50;
ens=np.linspace(5,6,f_num)   # ens 
#
syst=od.make_system(length,dis,'35')
sys1=syst.finalized()
##kwant.plot(sys,fig_size=(5, 10))
tp=kwant.smatrix(sys1,5.5).data
#t=[kwant.smatrix(sys,ens[cishu]).transmission(1,0) for cishu in range(50)]
#pyplot.plot(ens,t,'.')
syst0=od.make_system(length, 0,'0')   # part of system
greens_function_sites = syst0.sites()
od.mount_vlead(syst, greens_function_sites, 1)
sys=syst.finalized()
gf=np.diag(kwant.greens_function(sys,5.5,check_hermiticity=False).submatrix(2,2))
#gf=(kwant.greens_function(sys,5.5).submatrix(2,0))
ldos=-np.imag(gf)/np.pi
coord=np.array([sys.pos(i) for i in sys.lead_interfaces[2]]) 

wf_l=kwant.wave_function(sys1,5.5,check_hermiticity=False)(0)[0]
wf_r=kwant.wave_function(sys1,5.5,check_hermiticity=False)(1)[0]


#
flead0 = sys1.leads[0]
prop_modes, _ = flead0.modes(5.5)
v=prop_modes.velocities#*2*np.pi  # 0: left 1:right

#pyplot.plot(np.abs(gf)*np.sqrt(v[1]),'.')
#pyplot.plot(np.abs(wf_l))

#pyplot.plot(np.angle(gf))
#pyplot.plot(np.angle(wf_l))

#ldos_kwant = kwant.ldos(sys1, 5.5,check_hermiticity=False)

pyplot.plot(np.abs(wf_l)**2+np.abs(wf_r)**2,'.')
pyplot.plot(ldos*2*np.pi)
#pyplot.plot(ldos_kwant*2*np.pi,'.')
#kwant.smatrix(sys1,5.5).transmission(1,0)

def stat_wf(salt):
    wf_spec=np.zeros((f_num,length),dtype=complex)
    t_spec=[]
    sys=od.make_system(length,dis,str(salt))
    sys=sys.finalized()
    for i in range(f_num):
        wf_spec[i,]=kwant.wave_function(sys,ens[i])(0)
        t_spec.append(kwant.smatrix(sys,ens[i]).transmission(1,0))
    myDict = {'wf':wf_spec,'t':t_spec} 
    completeName = os.path.join('E:/pt/5/',str(salt)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
    
#Parallel(n_jobs=5)(delayed(stat_wf)(salt) for salt in np.arange(8000,9000,1))
elapsed=time()-t_ini