# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:59:55 2019

@author: ykang
"""
import numpy as np
import kwant
from kwant.digest import uniform

def eigentime(u0,u1,vh0,vh1,df):
    vv0=vh0.conj().T
    vv1=vh1.conj().T
    t=(u0.conj().T@(u1-u0)-vv0.conj().T@(vv1-vv0))/1j/df/(2*np.pi)
    return t

##### only for single mode case
    
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

def wf_integral(sys,e1,incident_site):  # incident_site=0/1 , only one mode case
    wf=kwant.wave_function(sys,e1)(incident_site)[0]
    return sum(abs(wf)**2)/2/np.pi

def wf_integral_all(sys,e1):  
    wf=kwant.wave_function(sys,e1)
    wf0=wf(0)[0]
    wf1=wf(1)[0]
    return np.array([sum(abs(wf0)**2)/2/np.pi, sum(abs(wf1)**2)/2/np.pi])

def trans_all(sys,e1):
    t=kwant.smatrix(sys,e1)
    t10=t.transmission(1,0)
    t01=t.transmission(0,1)
    t00=t.transmission(0,0)
    t11=t.transmission(1,1)
    return np.array([t10,t01,t00,t11])

def plot_wf(sys,e1,site):  # site=0/1
    wf=kwant.wave_function(sys,e1)(site)  #,check_hermiticity=False
    kwant.plotter.map(sys, (abs(wf[0]))**0.5,num_lead_cells=5,fig_size=(12, 10),colorbar=False)
    
#################################################
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b 

def make_system(width, length, salt, dis):
    def disk(pos):
        x,y=pos
        return -width<y<width and  abs(x)<length
#        return -width<y<width and  200<=x<600 #abs(x)<length  # #   #25.1 #0<y<width 

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*dis
#        if x>=200:
#            return (uniform(repr(x*y),salt)-0.5)*dis
#        else:
#            return (uniform(repr((x+400)*y),salt)-0.5)*dis

#        if y>-30:
#            return (uniform(repr(site),salt)-0.5)*dis
#        else:
#            return 0

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def make_lead(width):
    def lead_shape(pos):
        x,y=pos
        return -width<y<width #0<y<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, 0))] = 0 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
    return lead

def attach_lead(sys,width):
    lead=make_lead(width)
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
def make_system_all(width, length, salt, dis):    # disordered QSH system
    syst=make_system(width, length, salt,dis) 
    attach_lead(syst,width)
    return syst.finalized()

#def ballistic_velocity(width,en):     ## clean system, group velocity
#    flead=make_lead(width).finalized()
#    prop_modes, _ = flead.modes(en)
#    v=prop_modes.velocities[1]*2*np.pi  # 0: left 1:right
#    wf_lead=np.abs(prop_modes.wave_functions[:,1])
#    return np.sqrt(v)*wf_lead   

###sys.pos(sys.lead_interfaces[1][0])
