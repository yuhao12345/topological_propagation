# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 21:50:27 2018

@author: user
"""

import numpy as np
from matplotlib import pyplot
import scipy.io as sio
import os.path
#import random
import kwant
from kwant.digest import uniform
from random import choices
from time import time

t_ini=time()

width=15
ratio=100/100
length=100/ratio
salt='555'
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]],norbs=2)  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2
Lambda_R = 2
s0 = np.identity(2)
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1, -1])

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b
def disk(pos):
    x,y=pos
    return abs(y)<width and abs(x)<length   #25.1

def lead_shape(pos):
    x,y=pos
    return width-5<y<width  #abs(y)<width #

def lead_shape1(pos):
    x,y=pos
    return abs(y)<width  #width-5<y<width #

def onsite(site):
    return np.zeros([2,2])     


dis=1

def nnn(site1, site2):
    x,y=site1.pos
    if abs(x)<length*ratio and y>width-5:
##        return 1j * m2 * choices(population, weights)[0] * sz
        return 1j * m2 * (1-2*uniform(repr(site1),salt)*dis) * sz
    else:
        return 1j *m2 * sz
    
   
# center point of rashba region
#xr= (np.divide(random.sample(range(100), 30),100)-0.5)*2*length*ratio
#yr=width-np.divide(random.sample(range(100), 30),100)*5   #density  !!!
#def hop_ras_e1(site1,site2,Lambda_R=Lambda_R,t1=1):
#    x,y=site1.pos
#    if min(abs(x-xr)+abs(y-yr))<1.5:
#        return s0-1j*Lambda_R*sx
#    else:
#        return s0
#    
#
#def hop_ras_e2(site1,site2,Lambda_R=Lambda_R,t1=1):
#    x,y=site1.pos
#    if min(abs(x-xr)+abs(y-yr))<1.5:
#        return s0+1j*Lambda_R*(0.5*sx - np.sqrt(3)/2.0*sy)
#    else:
#        return s0
#    
#
#def hop_ras_e3(site1,site2,Lambda_R=Lambda_R,t1=1):
#    x,y=site1.pos
#    if min(abs(x-xr)+abs(y-yr))<1.5:
#        return s0+1j*Lambda_R*(0.5*sx + np.sqrt(3)/2.0*sy)
#    else:
#        return s0
    

def make_system2():
    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]=onsite
    sys[graphene.neighbors()]=s0      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn # 1j *m2 * sz #
#    sys[kwant.HoppingKind((0, 0), a,b)] = hop_ras_e1
#    sys[kwant.HoppingKind((0, 1), a,b)] = hop_ras_e2
#    sys[kwant.HoppingKind((-1, 1), a,b)] = hop_ras_e3

 
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym,conservation_law=-sz) #,conservation_law=-sz
    
    lead[graphene.shape(lead_shape, (0, width-1))] =  onsite# s0 #  
    lead[graphene.neighbors()]=s0
#    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 * sz
       
    sys.attach_lead(lead)
#    sys.attach_lead(lead.reversed())
    
    sym1 = kwant.TranslationalSymmetry((1,0))
    sym1.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym1.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    lead1 = kwant.Builder(sym1,conservation_law=-sz) #,conservation_law=-sz
    
    lead1[graphene.shape(lead_shape1, (0, width-1))] =  onsite#  s0 # 
    lead1[graphene.neighbors()]=s0
#    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 * sz      
    sys.attach_lead(lead1)
    
    return sys.finalized()

sys=make_system2()

en=0.5 #.519
wff=kwant.wave_function(sys,en)
wf=wff(0)
kwant.plotter.map(sys,(abs(wf[0][::2])**2),num_lead_cells=5, fig_size=(20, 10))
#t=kwant.smatrix(sys,en).transmission(1,0)
#tm=kwant.smatrix(sys,en)
#t10_uu=tm.transmission((1, 0), (0, 0))
#t10_du=tm.transmission((1, 1), (0, 0))
#t00_uu=tm.transmission((0, 0), (0, 0))
#t00_du=tm.transmission((0, 1), (0, 0))
#t10_dd=tm.transmission((1, 1), (0, 1))

##gf=kwant.greens_function(sys,.1).submatrix(1,0)
##smatrix.transmission((lead1, q1), (lead0, q0)) is the transmission from the
## q0 block of the lead0 into the q1 block of lead1.

#tm=kwant.smatrix(sys,en)
#gf_mode=tm.submatrix((1, 1), (0, 1)) 
#s=np.linalg.svd(gf_mode, full_matrices=False, compute_uv=False)  #u, s, vh
#f0=0.4
#f1=0.45
#df=1e-11   # It is freq not omega !!!
#f_points=1
#energies=np.linspace(f0,f1,f_points)
#energies1=energies+df
##
#### from TM
#def uv(energies):  # get u,v from TM or TM_gf, TM_gf is right, but TM does not work
#    u_mat=[]
#    v_mat=[]
#    s_mat=[]
#    u2_mat=[]
#    v2_mat=[]
#    u3_mat=[]
#    v3_mat=[]
#    for en in energies:
##        tm=kwant.smatrix(sys,en)
##        gf_mode=tm.submatrix((1, 0), (0, 0)) 
##        gf_mode1=tm.submatrix((1, 1), (0, 0)) 
##        gf1=np.concatenate((gf_mode, gf_mode1))
##        gf_mode=kwant.smatrix(sys,en).submatrix(1,0)
#        gf=kwant.greens_function(sys,en).submatrix(1,0)   # TM_gf not TM !!!
#        gf1=gf[:,0::2]
#        u, s, vh=np.linalg.svd(gf1, full_matrices=True, compute_uv=True)
#        u_mat.append(u[:,0])
#        v_mat.append(vh.conj().T[:,0])
#        u2_mat.append(u[:,1])
#        v2_mat.append(vh.conj().T[:,1])
#        u3_mat.append(u[:,2])
#        v3_mat.append(vh.conj().T[:,2])
#        s_mat.append(s)
#    return u_mat,v_mat,u2_mat,v2_mat,u3_mat,v3_mat,s_mat#
#
#u_mat,v_mat,u2_mat,v2_mat,u3_mat,v3_mat,s_mat=uv(energies)   #
#u_mat1,v_mat1,u2_mat1,v2_mat1,u3_mat1,v3_mat1,s_mat1=uv(energies1)#
##
#u0=np.column_stack(u_mat)
#v0=np.column_stack(v_mat)
#u1=np.column_stack(u_mat1)
#v1=np.column_stack(v_mat1)
##
#u20=np.column_stack(u2_mat)   # second eigenstate
#v20=np.column_stack(v2_mat)
#u21=np.column_stack(u2_mat1)
#v21=np.column_stack(v2_mat1)
#
##u30=np.column_stack(u3_mat)   # 3rd eigenstate
##v30=np.column_stack(v3_mat)
##u31=np.column_stack(u3_mat1)
##v31=np.column_stack(v3_mat1)
#
#dphi_tm1=np.diag((u0.conj().T @ (u1-u0)-v0.conj().T @ (v1-v0)))/1j/df/(2*np.pi)
#dphi2_tm1=np.diag((u20.conj().T @ (u21-u20)-v20.conj().T @ (v21-v20)))/1j/df/(2*np.pi)
##dphi3_tm=np.diag((u30.conj().T @ (u31-u30)-v30.conj().T @ (v31-v30)))/1j/df/(2*np.pi)
#
### from greens function
#def phase_gf(energies):
#    gf_uu=[]
#    gf_du=[]
#    for en in energies:
#        gf10=kwant.greens_function(sys,en).submatrix(1,0)
#        gf10_uu=gf10[0,0]
#        gf10_du=gf10[1,0]
#        gf_uu.append(gf10_uu)
#        gf_du.append(gf10_du)
#    return gf_uu,gf_du
#gfuu0, gfdu0=phase_gf(energies)
#gfuu1, gfdu1=phase_gf(energies1)    
#dphi_gfuu=(np.angle(gfuu1)-np.angle(gfuu0))/df/(2*np.pi)
#dphi_gfdu=(np.angle(gfdu1)-np.angle(gfdu0))/df/(2*np.pi)

## ldos 
#en=0.4
#local_dos = kwant.ldos(sys, en)
#ld = local_dos.reshape(-1, 2)
#
#coord=np.array([sys.pos(i) for i in range(ld.shape[0])])
#ans=np.concatenate((coord,ld),axis=1)
#myDict = {'ans':ans}
#completeName = os.path.join('E:/dwell/1/', str(0)+".mat")
#sio.savemat(completeName,myDict,oned_as='row') 

#en=0.4
#wff=kwant.wave_function(sys,en)(0)
#w=np.abs(wff.T)
#coord=np.array([sys.pos(i) for i in range(w.shape[0]//2)])
#up=np.concatenate((coord,w[0::2,]),axis=1)
#down=np.concatenate((coord,w[1::2,]),axis=1)
#myDict = {'u':up,'d':down}
#completeName = os.path.join('E:/dwell/13/', str(0)+".mat")
#sio.savemat(completeName,myDict,oned_as='row') 
elapsed=time()-t_ini

#sys.pos(290)
#sys.lead_interfaces
#en=0.2
#df=1e-8
#tm=kwant.smatrix(sys,en).data
#tm1=kwant.smatrix(sys,en+df).data
#wigner=tm.conj().T@(tm1-tm)/df/(2*np.pi)*(-1j)
#np.trace(wigner)

#local_dos = kwant.ldos(sys, .5)
#kwant.plotter.map(sys, local_dos[0::2] , num_lead_cells=5,fig_size=(10, 10))
#kwant.plotter.map(sys, local_dos[1::2] , num_lead_cells=5,fig_size=(10, 10))

#y=[sys.pos(sys.lead_interfaces[0][i])[1] for i in range(sys.lead_interfaces[0].size)]
