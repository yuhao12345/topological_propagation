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
import cmath
from numpy import linalg as LA
t_ini=time()

length=30
dis=1.5

lat = kwant.lattice.chain()

#o=np.array([4,4,4,4,4,4,4,4,4,4])
o=np.ones(length)*4
def onsite(site):
#    print(site.tag)   position of the site
    return o[site.tag]
syst = kwant.Builder()
syst[(lat(x) for x in range(length))] = onsite
syst[lat.neighbors()] = -1

lead = kwant.Builder(kwant.TranslationalSymmetry([-1]))
lead[(lat(0))] = 4
lead[lat.neighbors()] = -1

syst.attach_lead(lead)
syst.attach_lead(lead.reversed())
sys=syst.finalized()


f_num=50;
ens=np.linspace(3,5,f_num)   # ens 
##
#
#syst=od.make_system(length,dis,'35')
#sys1=syst.finalized()
#
######## for one sample
####  construct effecitive H_eff
ham=sys.hamiltonian_submatrix() 
ham[0,0]=ham[0,0]-cmath.exp(1j*2)  # add lead term
ham[-1,-1]=ham[-1,-1]-cmath.exp(1j*2)
w, v = LA.eig(ham)  # get QNM , w is complex energy
#pyplot.plot(np.abs(v[:,9])**2)
## filter w within [3.25,4.25]
#ind=(np.real(w)>3.25)*(np.real(w)<4.25)
#en_comp=w[ind]
#eigenwf=v[:,ind] # note it is normal;ized rather than biorthogonal!!! 
###kwant.plot(sys,fig_size=(5, 10))
#
#### spec
t=[kwant.smatrix(sys,ens[cishu]).transmission(1,0) for cishu in range(f_num)]
#t=[kwant.smatrix(sys1,ens[cishu]).submatrix(1,0) for cishu in range(f_num)]
#pyplot.plot(ens,t,'.')

wf=np.zeros((f_num,length),dtype=complex)
for i in range(f_num):
    wf[i,]=kwant.wave_function(sys,ens[i])(0)
    
pyplot.plot(np.abs(wf[25,])**2)

u=np.zeros((length,length),dtype=complex)
f=np.zeros(length,dtype=complex)
for ind in range(1):
    u[:,ind]=v[:,ind]/np.sqrt(np.sum(v[:,ind]**2))
    f=f+u[0,ind]*u[:,ind]/(ens[25]-w[ind])
pyplot.plot(np.abs(f)**2)

#syst0=od.make_system(length, 0,'0')   # part of system
#greens_function_sites = syst0.sites()
#od.mount_vlead(syst, greens_function_sites, 1)
#sys=syst.finalized()
#gf=np.diag(kwant.greens_function(sys,5.5,check_hermiticity=False).submatrix(2,2))
##gf=(kwant.greens_function(sys,5.5).submatrix(2,0))
#ldos=-np.imag(gf)/np.pi
#coord=np.array([sys.pos(i) for i in sys.lead_interfaces[2]]) 
#
#wf_l=kwant.wave_function(sys1,5.5,check_hermiticity=False)(0)[0]
#wf_r=kwant.wave_function(sys1,5.5,check_hermiticity=False)(1)[0]
#
#
##
#flead0 = sys1.leads[0]
#prop_modes, _ = flead0.modes(5.5)
#v=prop_modes.velocities#*2*np.pi  # 0: left 1:right

#pyplot.plot(np.abs(gf)*np.sqrt(v[1]),'.')
#pyplot.plot(np.abs(wf_l))

#pyplot.plot(np.angle(gf))
#pyplot.plot(np.angle(wf_l))

#ldos_kwant = kwant.ldos(sys1, 5.5,check_hermiticity=False)

#pyplot.plot(np.abs(wf_l)**2+np.abs(wf_r)**2,'.')
#pyplot.plot(ldos*2*np.pi)
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

def stat_QNM(salt):
    syst=od.make_system(length,dis,salt)
    sys1=syst.finalized()

    ham=sys1.hamiltonian_submatrix() 
    ham[0,0]=ham[0,0]-cmath.exp(1j*2)  # add lead term
    ham[-1,-1]=ham[-1,-1]-cmath.exp(1j*2)
    w, v = LA.eig(ham)  # get QNM , w is complex energy
    # filter w within [3.25,4.25]
    ind=(np.real(w)>3.25)*(np.real(w)<4.25)
    en_comp=w[ind]
    eigenwf=v[:,ind] # note it is normalized rather than biorthogonal!!! 
    myDict = {'en_comp':en_comp,'eigenwf':eigenwf} 
    completeName = os.path.join('E:/mode/1/',str(salt)+'.mat')
    sio.savemat(completeName,myDict,oned_as='row') 
#Parallel(n_jobs=5)(delayed(stat_QNM)(salt) for salt in np.arange(1000,5000,1))

elapsed=time()-t_ini