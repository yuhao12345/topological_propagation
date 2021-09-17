# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 10:36:48 2019

@author: ykang
"""
import numpy as np
import kwant
from kwant.digest import uniform

lat = kwant.lattice.square()

def make_system_all(length,width,dis,salt,command='null'):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width 
    def onsite(site):
        x,y=site.pos
        if command=='crystal':
            if (x%2==0) & (y%2==0):
                return 2
            else:
                return 0
        else:
            if abs(x)<length:    # add clean region
                return dis*(uniform(repr(site),salt)-.5)+4
            else:
                return 4

    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = -1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 4
    lead[lat.neighbors()] = -1

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

def make_system_crys_air_crys(length,width,width_lead,dis,salt):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    def lead_shape(pos):
        x, y= pos
        return abs(y)<width_lead
    def onsite(site):
        x,y=site.pos
        if abs(y)<width_lead:
            return 0
        elif y>=width_lead:
            return ((x%2==0) & (y%2==0))*2
        else:
            return ((x%2==0) & (y%2==0))*2
    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = -1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)

    lead[lat.shape(lead_shape, (0, 0))] = 0
    lead[lat.neighbors()] = -1

    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())

    return sys.finalized()

def tranmission_DIY(sys,energy):  # based on fisher-lee relation, different from kwant
    Gr = kwant.greens_function(sys, energy).submatrix(1,0)
    flead0 = sys.leads[0]
    prop_modes, _ = flead0.modes(energy)
    v=prop_modes.velocities#*2*np.pi  # 0: left 1:right
    n=v.size//2  # number of channel
    v_matrix_sqrt= np.diag([v[i]**0.5 for i in range(n,2*n)])
#    v_matrix_sqrt2= np.diag([i**1 for i in range(n,2*n)])
    wf_lead=prop_modes.wave_functions[:,n:2*n]
    wf_lead_n=np.sum(np.abs(wf_lead)**2,0)**0.5
    wf_lead_unit=wf_lead/wf_lead_n
    t_s=1j*v_matrix_sqrt @ (wf_lead_unit).T  @ Gr @ np.conj(wf_lead_unit)  @ v_matrix_sqrt 
    return t_s