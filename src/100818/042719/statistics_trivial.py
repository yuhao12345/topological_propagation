# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 12:48:35 2019

@author: ykang
"""

import stru_QSH as dg
import stru_QVH_air_QVH as qa
import numpy as np
import kwant
from joblib import Parallel, delayed
import scipy.io as sio
import os.path
from time import time
import pandas as pd
import matplotlib.pyplot as plt
t_ini=time()


width=30
width_lead=30
length=50
dis=3.5

sys0=qa.make_system_all(width, width_lead, length, '0', dis) 

#e1=0.19
#wf=kwant.wave_function(sys0,e1)(0)   #.403586
#kwant.plotter.map(sys0, (abs(wf[0])),num_lead_cells=5,fig_size=(15, 10),colorbar=False)

#t=[kwant.smatrix(sys0,e1).transmission(1,0) for e1 in np.linspace(0.2,0.3,100)]
#plt.plot(np.linspace(0.2,0.3,100),t)

#coord=qa.coordinate(sys0,e1)
#[sys0.pos(sys0.lead_interfaces[0][i]) for i in range(20)]
#sys0.pos(sys0.lead_interfaces[1][5])
gf=[kwant.greens_function(sys0,e1).submatrix(1,0)[5,5] for e1 in np.linspace(0.18,0.22,50)]
plt.plot(np.linspace(0.18,0.22,50),np.unwrap(np.angle(gf)))    
plt.plot(np.linspace(0.18,0.22,50),np.abs(np.angle(gf)))        
#def stat_I(cishu):
#    sys=qa.make_system_all(width, width_lead, length, str(cishu), dis) 
#    return qa.I_integral(sys,e1,coord,width_lead)
#
#I=Parallel(n_jobs=8)(delayed(stat_I)(cishu) for cishu in np.arange(0,1000,1))
#
#myDict = {'I':I} #,'ld':ld
#completeName = os.path.join('E:/dwell3/829/4.mat')
#sio.savemat(completeName,myDict,oned_as='row')    

elapsed=time()-t_ini