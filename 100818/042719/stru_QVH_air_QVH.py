# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 14:33:55 2019

@author: ykang
"""

import numpy as np
import kwant
from kwant.digest import uniform
import pandas as pd
import stru_QSH as dg

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m1 = 0.52   # for valley

def make_system(width, width_lead, length, salt, dis):
    def disk(pos):
        x,y=pos
        return -width<y<10 and abs(x)<length #abs(x)<length   #25.1 #

    def onsite(site):
        x,y=site.pos
        if abs(y)<3:
            return 0;
        elif x**2+(y+8)**2<9:
            return -m1 if site.family == a else m1
        else:
            return m1 if site.family == a else -m1
#        if y>=2:
#            return m1 if site.family == a else -m1
#        elif abs(y)<3:
#            return 0
#        else:
#            onsite_a = m1+dis*(uniform(repr(site.tag),salt)-.5)
#            return onsite_a if site.family == a else -onsite_a
     
    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,0))]=onsite
    sys[graphene.neighbors()]=-1
    return sys

def attach_lead(sys,width_lead):
    def lead_shape(pos):
        x,y=pos
        return -width_lead<=y<10
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, 0))] = 0 
    lead[graphene.neighbors()]=-1
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
def make_system_all(width, width_lead, length, salt, dis):    # disordered QSH system
    syst=make_system(width, width_lead, length, salt, dis) 
    attach_lead(syst,width_lead)
    return syst.finalized()    

###  I
    
def coordinate(sys,e1):
    wf=kwant.wave_function(sys,e1)(0)
    return np.array([sys.pos(i) for i in range(wf.shape[1])])

def I_integral(sys,e1,coord,width_lead):   # integral over cross section of lead
    wf=kwant.wave_function(sys,e1)(0)
    df=pd.DataFrame({'x':coord[:,0],'y':coord[:,1],'wf':abs(wf[0,:])**2})  
    df1=df[abs(df.y)<width_lead]
    wf_int=df1.groupby(['x']).sum()
    return np.array(wf_int.wf)




