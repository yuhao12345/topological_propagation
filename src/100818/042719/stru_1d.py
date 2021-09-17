# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:26:48 2019

@author: ykang
"""

import numpy as np
import kwant
from kwant.digest import uniform
#from matplotlib import pyplot

lat = kwant.lattice.chain()
def make_system(length,dis,salt):
    def onsite(site):
        return dis*(uniform(repr(site),salt)-.5)+4 #+0.005j
    syst = kwant.Builder()
    syst[(lat(x) for x in range(length))] = onsite
    syst[lat.neighbors()] = -1

    lead = kwant.Builder(kwant.TranslationalSymmetry([-1]))
    lead[(lat(0))] = 4
    lead[lat.neighbors()] = -1

    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())
    return syst


def mount_vlead(sys, vlead_interface, norb):
    """Mounts virtual lead to interfaces provided.

    :sys: kwant.builder.Builder
        An unfinalized system to mount leads
    :vlead_interface: sequence of kwant.builder.Site
        Interface of lead
    :norb: integer
        Number of orbitals in system hamiltonian.
    """
    dim = len(vlead_interface)*norb
    zero_array = np.zeros((dim, dim), dtype=float)
    def selfenergy_func(energy, args=()):
        return zero_array

    vlead = kwant.builder.SelfEnergyLead(selfenergy_func, vlead_interface)
    sys.leads.append(vlead)

########### for one configuration
#di=11
#def make_system0(width, length):
#    def disk(pos):
#        x,y=pos
#        return width-di<y<width and -length<x<3*length   #25.1
#    sys=kwant.Builder()
#    sys[graphene.shape(disk,(300,width-1))]= 0  #0 #
#    sys[graphene.neighbors()]=1      # comment it, when has rashba
#    def lead_shape(pos):
#        x,y=pos
#        return width-di<y<width 
#    sym = kwant.TranslationalSymmetry((-1,0))
#    sym.add_site_family(graphene.sublattices[0], other_vectors=[(-1, 2)])
#    sym.add_site_family(graphene.sublattices[1], other_vectors=[(-1, 2)])    
#    lead = kwant.Builder(sym) #,conservation_law=-sz    
#    lead[graphene.shape(lead_shape, (0, width-1))] = 0 #1
#    lead[graphene.neighbors()]=1       
#    sys.attach_lead(lead)
#    sys.attach_lead(lead.reversed())
#
#    return sys





