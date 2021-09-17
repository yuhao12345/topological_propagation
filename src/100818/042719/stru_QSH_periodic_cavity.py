# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 18:35:16 2019

@author: ykang
"""
import numpy as np
import kwant
import stru_QSH as dg

graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

m2 = .1  #spin    3*np.sqrt(3)*m2

#center=np.linspace(-80,80,11)
center=np.linspace(-50,50,1)
nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b 

def make_system(width, length):
    def disk(pos):
        x,y=pos
        return 0<y<width and abs(x)<length #abs(x)<length   #25.1 #

    def nnn(site1, site2):
        x,y=site1.pos
        for i in range(center.size):
            if (x-center[i])**2+(y-23)**2<4**2:
#        if x**2+(y-23)**2<4**2 or (x-50)**2+(y-23)**2<4**2:
                return -1j * m2
        return 1j *m2 
    
    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,width-1))]= 4-0.0001j
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = nnn
    return sys

def attach_lead(sys,width):
    def lead_shape(pos):
        x,y=pos
        return 0<y<width #width>y>width-8 #
    sym = kwant.TranslationalSymmetry((-1,0))
    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])
    
    lead = kwant.Builder(sym)
    
    lead[graphene.shape(lead_shape, (0, width-1))] = 4 
    lead[graphene.neighbors()]=1
    lead[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]]=1j *m2 
       
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
def make_system_all_periodic_cavity(width, length):    # disordered QSH system
    syst=make_system(width, length) 
    attach_lead(syst,width)
    return syst.finalized()

#sys=make_system_all_periodic_cavity(60, 100)
#dg.plot_wf(sys,.45,0)
