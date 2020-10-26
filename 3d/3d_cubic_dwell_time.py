# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 23:14:06 2018

@author: user
"""

import kwant
import random
import numpy as np
from time import time
import matplotlib.pyplot as plt
from kwant.digest import uniform
t_ini=time()

lat = kwant.lattice.general([(0, 0, 1), (0, 1, 0), (1, 0, 0)],
                            [(0, 0, 0)])

hop=0.5
def make_cuboid(a,seed):
    def cuboid_shape(pos):
        x, y, z = pos
        return 0 <= x < a and 0 <= y < a and 0 <= z < a
    
    def onsite(site):
        x,y,z=site.pos
        return (uniform(repr(site),seed)-0.5)*dis
#        return (random.random()-0.5)*dis

    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0, 0))] = onsite
    sys[lat.neighbors()] = hop
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0,0))
    lead = kwant.Builder(sym_lead)

    def lead_shape(pos):
        x, y, z = pos
        return x==0 and 0 <= y < a and 0 <= z < a

    lead[lat.shape(lead_shape, (0, 0, 0))] = 0
    lead[lat.neighbors()] = hop

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

dtime=[]
freq_num=6
energies=np.linspace(1.99,2,freq_num)
dis=6
size=np.arange(4,16,1)
for a in size:
    dwell=[]
    for cishu in range(50):
        salt=np.random.randint(10000)
        sys = make_cuboid(a,str(salt))
        gf=[]
        for en in energies:
            gf.append(kwant.greens_function(sys,en).submatrix(1,0)[0,0])
        phase=np.unwrap(np.angle(gf))
#        plt.plot(phase,'.')
        phase_de=np.mean((phase[np.arange(1,freq_num,2)]-phase[np.arange(0,freq_num,2)])/(energies[1]-energies[0])/2/np.pi)
        
        dwell.append(phase_de)
    dtime.append(dwell)

elapsed=time()-t_ini        

#plt.plot(np.mean(dtime,1),'.')
plt.plot(dtime,'.')
plt.ylabel('<dwell>')
plt.xlabel('size')
plt.rc('font', size=5)  
plt.tick_params(labelsize=14)