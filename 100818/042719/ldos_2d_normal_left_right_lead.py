# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 15:58:01 2019

@author: ykang
"""

import stru_QSH as dg
import stru_1d as od
import stru_2d_normal_left_right_lead as dnlr
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
from matplotlib import pyplot

t_ini=time()

width=15
length=50
width_left=6
width_right=6
dis=1.3
salt='10'

en=.45  # 2.1 when onsite_lead=6

syst=dnlr.make_system_all(length,width,dis,salt,width_left,width_right)
sys1=syst.finalized()  # real structure
#kwant.plot(sys,fig_size=(15, 10))
wf_kwant=kwant.wave_function(sys1,en)(0)[1]
#kwant.plotter.map(sys1, (abs(wf[0]))**0.5,num_lead_cells=5,fig_size=(10, 10),colorbar=False)

## attcah virtual lead
syst0=dnlr.make_system_all(length,width,0,'0',width_left,width_right) # virtual lead
greens_function_sites = syst0.sites()
od.mount_vlead(syst, greens_function_sites, 1)
sys=syst.finalized()

#gf=np.diag(kwant.greens_function(sys,en).submatrix(2,2))
#ldos=-np.imag(gf)/np.pi

gf20=kwant.greens_function(sys,en).submatrix(2,0)
coord_virtual=np.array([sys.pos(i) for i in sys.lead_interfaces[2]]) # coord of virtual lead

flead0 = sys.leads[0]
prop_modes, _ = flead0.modes(en)
v=prop_modes.velocities#*2*np.pi  # 0: left 1:right
wf_lead=prop_modes.wave_functions
n=v.size//2   # number of channel
wf=1j*gf20@wf_lead[:,3]*np.sqrt(v[3])/np.sqrt(sum(abs(wf_lead[:,3])**2))

#ldos_kwant = kwant.ldos(sys1, en)
coord_sys=np.array([sys1.pos(i) for i in range(sys.lead_interfaces[2].size)]) 

#kwant.plotter.map(sys1,ldos_kwant)
#myDict = {'ldos':ldos,'ldos_kwant':ldos_kwant,'coord_virtual':coord_virtual,'coord_sys':coord_sys} 
#completeName = os.path.join('E:/pt/ldos3.mat')
#sio.savemat(completeName,myDict,oned_as='row') 

myDict = {'wf':wf,'wf_kwant':wf_kwant,'coord_virtual':coord_virtual,'coord_sys':coord_sys} 
completeName = os.path.join('E:/pt/wf2.mat')
sio.savemat(completeName,myDict,oned_as='row') 

elapsed=time()-t_ini