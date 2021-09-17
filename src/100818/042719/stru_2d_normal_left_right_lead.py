# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 15:17:18 2019

@author: ykang
"""

import numpy as np
import kwant
from kwant.digest import uniform

lat = kwant.lattice.square()

def make_system_all(length,width,dis,salt,width_left,width_right):
    def cuboid_shape(pos):
        x, y= pos
        return abs(x)<length and abs(y)<width
    def lead_shape_left(pos):
        x, y= pos
        return abs(y)<width_left 
    def lead_shape_right(pos):
        x, y= pos
        return abs(y)<width_right
    def onsite(site):
        x,y=site.pos
        return dis*(uniform(repr(site),salt)-.5)+4#-0.001j

    
    sys = kwant.Builder()
    sys[lat.shape(cuboid_shape, (0, 0))] = onsite
    sys[lat.neighbors()] = -1
    
    sym_lead = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym_lead)
    lead[lat.shape(lead_shape_left, (0, 0))] = 4
    lead[lat.neighbors()] = -1
    sys.attach_lead(lead)
    
    sym_lead_r = kwant.TranslationalSymmetry((1, 0))
    lead_r = kwant.Builder(sym_lead_r)
    lead_r[lat.shape(lead_shape_right, (0, 0))] = 4
    lead_r[lat.neighbors()] = -1
    sys.attach_lead(lead_r)

    return sys

def make_lead(width):
    def lead_shape_left(pos):
        x,y=pos
        return abs(y)<width #0<y<width #width>y>width-8 #
    sym_lead = kwant.TranslationalSymmetry((1, 0))
    lead = kwant.Builder(sym_lead)
    lead[lat.shape(lead_shape_left, (0, 0))] = 4
    lead[lat.neighbors()] = -1
    return lead

def ballistic_velocity(width,en):     ## clean system, group velocity
    flead=make_lead(width).finalized()
    prop_modes, _ = flead.modes(en)
    v=prop_modes.velocities[1]#*2*np.pi  # 0: left 1:right
    wf_lead=(prop_modes.wave_functions[:,1])
    return np.sqrt(v)*wf_lead  

#flead=make_lead(3).finalized()
#prop_modes, _ = flead.modes(0.45)
#wf_lead=np.abs(prop_modes.wave_functions[:,1])
#v=prop_modes.velocities[1]*2*np.pi  # 0: left 1:right
#wf_v=np.sqrt(v)*wf_lead  

def eigentime(u0,u1,vh0,vh1,df):
    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t

def t_gf(sys,e1,df,o,i):   # output_site,incident_site, only one mode case
    g=kwant.greens_function(sys,e1)
    g1=kwant.greens_function(sys,e1+df)
    try:
        gf=g.submatrix(o,i)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(o,i)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
        
        t=np.diag(eigentime(u0,u1,vh0,vh1,df)) # time ,not dos
    except:
        t=0  
    return t

def t_gf_all(sys,e1,df):   # only one mode case
    g=kwant.greens_function(sys,e1)
    g1=kwant.greens_function(sys,e1+df)
    try:
        gf=g.submatrix(1,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(1,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t10=np.diag(eigentime(u0,u1,vh0,vh1,df))  # time ,not dos
    except:
        t10=0  
    try:
        gf=g.submatrix(0,1)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(0,1)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t01=np.diag(eigentime(u0,u1,vh0,vh1,df))  # time ,not dos
    except:
        t01=0
    try:
        gf=g.submatrix(0,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(0,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t00=np.diag(eigentime(u0,u1,vh0,vh1,df))  # time ,not dos
    except:
        t00=0  
    try:
        gf=g.submatrix(1,1)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(1,1)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t11=np.diag(eigentime(u0,u1,vh0,vh1,df)) # time ,not dos
    except:
        t11=0  
    return np.array([t10,t01,t00,t11])

def tranmission_DIY(sys,energy):  # based on fisher-lee relation, different from kwant
    Gr = kwant.greens_function(sys, energy,check_hermiticity=False).submatrix(1,0)
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

def t_mode_DIY(sys,e1,df,o,i):  
    g=tranmission_DIY(sys,e1)
    g1=tranmission_DIY(sys,e1+df)
    try:
        u0, s0, vh0=np.linalg.svd(g, full_matrices=True, compute_uv=True)
        u1, s1, vh1=np.linalg.svd(g1, full_matrices=True, compute_uv=True)        
        t=np.diag(eigentime(u0,u1,vh0,vh1,df))  # time ,not dos
    except:
        t=0  
    return t

def t_mode(sys,e1,df,o,i):   # only one mode case
    g=kwant.smatrix(sys,e1)
    g1=kwant.smatrix(sys,e1+df)
    try:
        gf=g.submatrix(o,i)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(o,i)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t=np.diag(eigentime(u0,u1,vh0,vh1,df))  # time ,not dos
    except:
        t=0  
    return t
    
def t_mode_all(sys,e1,df):   # only one mode case
    g=kwant.smatrix(sys,e1)
    g1=kwant.smatrix(sys,e1+df)
    try:
        gf=g.submatrix(1,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(1,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t10=np.diag(eigentime(u0,u1,vh0,vh1,df))[0]  # time ,not dos
    except:
        t10=0  
    try:
        gf=g.submatrix(0,1)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(0,1)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t01=np.diag(eigentime(u0,u1,vh0,vh1,df))[0]  # time ,not dos
    except:
        t01=0
    try:
        gf=g.submatrix(0,0)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(0,0)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t00=np.diag(eigentime(u0,u1,vh0,vh1,df))[0]  # time ,not dos
    except:
        t00=0  
    try:
        gf=g.submatrix(1,1)
        u0, s0, vh0=np.linalg.svd(gf, full_matrices=True, compute_uv=True)
        gf1=g1.submatrix(1,1)
        u1, s1, vh1=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)        
        t11=np.diag(eigentime(u0,u1,vh0,vh1,df))[0]  # time ,not dos
    except:
        t11=0  
    return np.array([t10,t01,t00,t11])

def wf_integral_all(sys,e1):  # incident_site=0/1 , only one mode case
    wf=kwant.wave_function(sys,e1)
    wf0=wf(0)[0]
    wf1=wf(1)[0]
    return np.array([sum(abs(wf0)**2)/2/np.pi, sum(abs(wf1)**2)/2/np.pi])