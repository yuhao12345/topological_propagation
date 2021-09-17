# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:34:19 2019

@author: ykang
"""
import numpy as np
import kwant
from kwant.digest import uniform
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
        return 0<y<width and abs(x)<length #abs(x)<length   #25.1 #

    def onsite(site):
        x,y=site.pos
        return (uniform(repr(site),salt)-0.5)*dis+4-0.0001j

    sys=kwant.Builder()
    sys[graphene.shape(disk,(0,width-1))]= onsite  #0 #
    sys[graphene.neighbors()]=1      # comment it, when has rashba
    sys[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 1j*m2
    return sys

def make_lead(width):
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
    return lead

def attach_lead(sys,width):
    lead=make_lead(width)
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed())
    
def make_system_all(width, length, salt, dis):    # disordered QSH system
    syst=make_system(width, length, salt,dis) 
    attach_lead(syst,width)
    return syst.finalized()