# -*- coding: utf-8 -*-
"""
Created on Thu May 21 16:07:44 2020

@author: user
"""

import numpy as np
import kwant
#from matplotlib import pyplot

lat = kwant.lattice.chain()
def make_system(length):
    def onsite(site):
        x=site.pos[0]
        if x<0:
            return 2+0.05j
        return 2-0.05j


    syst = kwant.Builder()
    syst[(lat(x) for x in range(-length,length))] = onsite
    syst[lat.neighbors()] = -1

    lead = kwant.Builder(kwant.TranslationalSymmetry([-1]))
    lead[(lat(0))] = 1
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

sys=make_system(10)
sys=sys.finalized()
kwant.plot(sys,fig_size=(5, 10))
s=kwant.smatrix(sys,1.5,check_hermiticity=False).data
abs(np.linalg.det(s))