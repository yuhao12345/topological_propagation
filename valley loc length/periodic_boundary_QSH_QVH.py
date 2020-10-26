# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 23:50:20 2018

@author: user     QSH+QVH PBC
"""
import numpy as np
import kwant
import matplotlib.pyplot
import random
import scipy.io as sio
import os.path
from kwant.digest import uniform
import time
#
start=time.time()
#def main() :
X = 10
Y = 30
yl = Y/4
d=3   #depth of disorder
yu = .75*Y  # yl-yu is boundary of QSH
m1,m2=.5,.5

dis=.2
en=0.52
salt=str(random.random()+5)
s0=np.identity(2)
#sx=np.array([[0,1],[1,0]])
sz=np.array([[1,0],[0,-1]])

graphene = kwant.lattice.honeycomb(1,'b')
A,B=graphene.sublattices

def rectangle(pos):
    x, y = pos
    return -X/2 < x < X/2

def onsite_qsvh(site):
    x,y=site.pos
    if yl<y<yu:
        return np.zeros([2,2])        
    else: 
        onsite_a = m1*s0
        onsite_b = -m1*s0
        return onsite_a if site.family == A else onsite_b
 
        
def hopp_qsvh(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    hop_a_dis = (m2+(uniform(repr(site1),salt)-.5)*dis)*1j*sz
    hop_b_dis = -(m2+(uniform(repr(site1),salt)-.5)*dis)*1j*sz
    hop_a = m2*1j*sz
    hop_b = -m2*1j*sz
    if yl+d<y1<yu-d:
        return hop_a if site1.family == A else hop_b
    if yl<y1<=yl+d or yu-d<=y1<yu:
        return hop_a_dis if site1.family == A else hop_b_dis
    else:
        return np.zeros([2,2])

def onsite_qsvh_lead(site):
    x,y=site.pos
    if yl<y<yu:
        return np.zeros([2,2])        
    else: 
        onsite_a = m1*s0
        onsite_b = -m1*s0
        return onsite_a if site.family == A else onsite_b
def hopp_qsvh_lead(site1,site2):
    x1,y1=site1.pos
    x2,y2=site2.pos
    hop_a = m2*1j*sz
    hop_b = -m2*1j*sz
    if yl<y1<yu and yl<y2<yu:
        return hop_a if site1.family == A else hop_b
    else:
        return np.zeros([2,2])

sym = kwant.TranslationalSymmetry(graphene.vec((-Y/2,Y)))
anc = kwant.Builder(sym) ### 2D periodic conditions
anc[graphene.shape(rectangle,(0, 0))] = None
anc[graphene.neighbors()] = None 
sys = kwant.Builder()
  
sys[anc.sites()] = onsite_qsvh_lead
sys[((a, sym.to_fd(b)) for a, b in anc.hoppings())] = -s0 
sys[A.neighbors()]=hopp_qsvh
sys[B.neighbors()]=hopp_qsvh

sym_anc = kwant.TranslationalSymmetry(graphene.vec((1,0)),graphene.vec((-Y/2,Y)))
anc_left = kwant.Builder(sym_anc)

sym_left = kwant.TranslationalSymmetry(graphene.vec((1,0)))
lead_left = kwant.Builder(sym_left)

anc_left[graphene.shape(lambda p: True,(0, 0))] = None
anc_left[graphene.neighbors()] = None

lead_left[anc_left.sites()] = onsite_qsvh_lead
lead_left[((a, sym.to_fd(b)) for a, b in anc_left.hoppings())] = -s0
lead_left[A.neighbors()]=hopp_qsvh_lead
lead_left[B.neighbors()]=hopp_qsvh_lead

sys.attach_lead(lead_left)
sys.attach_lead(lead_left.reversed())

sys = sys.finalized() 
#kwant.plot(sys,fig_size=(20, 10))

#kwant.plot(lead_left)

#save_path = 'E:/loc scale/6/'
#for cishu in np.arange(1,2):
#    salt=str(cishu+random.random())
#    try:
##        u, s, v = np.linalg.svd(gf_mode, full_matrices=True)
##        wf=kwant.wave_function(sys,en)(0)
#        
#        gf_mode=kwant.smatrix(sys,en)
#        T=gf_mode.transmission(1,0)
#        s00=gf_mode.submatrix(0,0)
#        s10=gf_mode.submatrix(1,0)
#        s01=gf_mode.submatrix(0,1)
#        s11=gf_mode.submatrix(1,1)
#        completeName = os.path.join(save_path, str(cishu)+".mat")
#        myDict = {'s00':s00,'s01':s01,'s10':s10,'s11':s11,
#                  't':T,'salt':salt,'length':X,'width':Y}
#        sio.savemat(completeName,myDict,oned_as='row')
#    except:
#        continue
gf_mode=kwant.smatrix(sys,en)
print(gf_mode.transmission(1,0))
wff=kwant.wave_function(sys,en)
wf=wff(0)
wavef1=[]
for i in range(wf.shape[1]//2):
    wavef1.append(wf[0,2*i+1])
kwant.plotter.map(sys,(abs(np.array(wavef1))**2), fig_size=(20, 10))

end=time.time()
print(end-start)

#if __name__ == '__main__':
#    main()