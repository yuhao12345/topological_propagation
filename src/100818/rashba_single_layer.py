# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 21:36:18 2018

@author: user
"""
# original code from website https://www.mail-archive.com/kwant-discuss@kwant-project.org/msg01041.html
import kwant
from matplotlib import pyplot
import tinyarray
from numpy import sqrt
from math import *

sys=kwant.Builder()
lat_u = kwant.lattice.honeycomb(1, name='up')
a,b = lat_u.sublattices

lat_d = kwant.lattice.honeycomb(1, name='down')
c,d = lat_d.sublattices

lx,ly=11,7
Lambda_R=2
Lambda=0.0
Lambda_nu=0#1.45

def Rectangle(pos,lx=lx,ly=ly):
    x,y=pos
    return 0<=x<lx and 1<=y<ly

def onsite(site,Lambda_nu=Lambda_nu):
	return Lambda_nu if (site.family == a or site.family == c) else -Lambda_nu

           
def Hop_Value_u(site1,site2,Lambda=Lambda):
    return 1j*Lambda*(2./sqrt(3))*1

def Hop_Value_d(site1,site2,Lambda=Lambda):
    return -1j*Lambda*(2./sqrt(3))

def hop_ras_e1(site1,site2,Lambda_R=Lambda_R):
    return 1j*Lambda_R

def hop_ras_e2_u(site1,site2,Lambda_R=Lambda_R):
    return 1j*Lambda_R*(0.5 + 1j*sqrt(3)/2.0)

def hop_ras_e3_u(site1,site2,Lambda_R=Lambda_R):
    return 1j*Lambda_R*(0.5 - 1j*sqrt(3)/2.0)


def hop_ras_e2_d(site1,site2,Lambda_R=Lambda_R):
    return 1j*Lambda_R*(0.5 - 1j*sqrt(3)/2.0)

def hop_ras_e3_d(site1,site2,Lambda_R=Lambda_R):
    return 1j*Lambda_R*(0.5 + 1j*sqrt(3)/2.0)




sys[lat_u.shape(Rectangle,(0,1))]=onsite   #even with 0 onsite potential, we need to use the zero matrix and not a scalar.
sys[lat_d.shape(Rectangle,(0,1))]=onsite

sys[lat_u.neighbors()] = 1
sys[lat_d.neighbors()] = 1

sys[kwant.HoppingKind((0, 0), a,d)] = hop_ras_e1
sys[kwant.HoppingKind((0, 1), a,d)] = hop_ras_e2_u
sys[kwant.HoppingKind((-1, 1), a,d)] = hop_ras_e3_u

sys[kwant.HoppingKind((0, 0), c,b)] = hop_ras_e1
sys[kwant.HoppingKind((0, 1), c,b)] = hop_ras_e2_d
sys[kwant.HoppingKind((-1, 1), c,b)] = hop_ras_e3_d

#Here we chose the clockwise Vij which are +1 (as retuned by Hop_value). The other ones (-1) will be given directly 
#by the hermeticity of the Hamiltonian
hoppings1_u = (((1, 0), a,a), ((-1, 1), a,a), ((0, -1), a,a))    
hoppings2_u = (((1, 0), b,b), ((-1, 1), b,b), ((0, -1), b,b))

hoppings1_d = (((1, 0), c,c), ((-1, 1), c,c), ((0, -1), c,c))    
hoppings2_d = (((1, 0), d,d), ((-1, 1), d,d), ((0, -1), d,d))
sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1_u]] = Hop_Value_u
sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2_u]] = Hop_Value_u

sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1_d]] = Hop_Value_d
sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2_d]] = Hop_Value_d
#================================left_up=====================================
def lead0_shape_u(pos,ly=ly):
    x,y=pos
    return    1<=y<ly

sym0u = kwant.TranslationalSymmetry((-1,0))
sym0u.add_site_family(lat_u.sublattices[0], other_vectors=[(-1, 2)])
sym0u.add_site_family(lat_u.sublattices[1], other_vectors=[(-1, 2)])
lead0_u = kwant.Builder(sym0u)

lead0_u[lat_u.shape(lead0_shape_u,(0,1))]=onsite

lead0_u[lat_u.neighbors()] = 1

lead0_u[kwant.HoppingKind((0, 0), a,d)] = hop_ras_e1
lead0_u[kwant.HoppingKind((0, 1), a,d)] = hop_ras_e2_u
lead0_u[kwant.HoppingKind((-1, 1), a,d)] = hop_ras_e3_u

lead0_u[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1_u]] = Hop_Value_u
lead0_u[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2_u]] = Hop_Value_u

sys.attach_lead(lead0_u)
#================================left_down0=====================================
def lead0_shape_d(pos,ly=ly):
    x,y=pos
    return    1<=y<ly

sym0d = kwant.TranslationalSymmetry((-1,0))
sym0d.add_site_family(lat_d.sublattices[0], other_vectors=[(-1, 2)])
sym0d.add_site_family(lat_d.sublattices[1], other_vectors=[(-1, 2)])
lead0_d = kwant.Builder(sym0d)

lead0_d[lat_d.shape(lead0_shape_d,(0,1))]=onsite

lead0_d[lat_d.neighbors()] = 1

lead0_d[kwant.HoppingKind((0, 0), c,b)] = hop_ras_e1
lead0_d[kwant.HoppingKind((0, 1), c,b)] = hop_ras_e2_d
lead0_d[kwant.HoppingKind((-1, 1), c,b)] = hop_ras_e3_d

lead0_d[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1_d]] = Hop_Value_d
lead0_d[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2_d]] = Hop_Value_d

sys.attach_lead(lead0_d)
#=================================right_up========================================
def lead1_shape_u(pos,ly=ly):
    x,y=pos
    return    1<=y<ly
    
sym1u = kwant.TranslationalSymmetry((1,0))
sym1u.add_site_family(lat_u.sublattices[0], other_vectors=[(-1, 2)])
sym1u.add_site_family(lat_u.sublattices[1], other_vectors=[(-1, 2)])
lead1_u = kwant.Builder(sym1u)
lead1_u[lat_u.shape(lead1_shape_u, (0, 1))] = onsite
lead1_u[lat_u.neighbors()] = 1

lead1_u[kwant.HoppingKind((0, 0), a,d)] = hop_ras_e1
lead1_u[kwant.HoppingKind((0, 1), a,d)] = hop_ras_e2_u
lead1_u[kwant.HoppingKind((-1, 1), a,d)] = hop_ras_e3_u

lead1_u[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1_u]] = Hop_Value_u
lead1_u[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2_u]] = Hop_Value_u

sys.attach_lead(lead1_u)
#=================================right_down========================================
def lead1_shape_d(pos,ly=ly):
    x,y=pos
    return    1<=y<ly
    
sym1d = kwant.TranslationalSymmetry((1,0))
sym1d.add_site_family(lat_d.sublattices[0], other_vectors=[(-1, 2)])
sym1d.add_site_family(lat_d.sublattices[1], other_vectors=[(-1, 2)])
lead1_d = kwant.Builder(sym1d)
lead1_d[lat_d.shape(lead1_shape_d, (0, 1))] = onsite
lead1_d[lat_d.neighbors()] = 1

lead1_d[kwant.HoppingKind((0, 0), c,b)] = hop_ras_e1
lead1_d[kwant.HoppingKind((0, 1), c,b)] = hop_ras_e2_d
lead1_d[kwant.HoppingKind((-1, 1), c,b)] = hop_ras_e3_d

lead1_d[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1_d]] = Hop_Value_d
lead1_d[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2_d]] = Hop_Value_d

sys.attach_lead(lead1_d)
#==============================================================================
def family_colors(site):
        return 0 if (site.family == a or site.family == c)  else 1
def hopping_lw(site1, site2):
        return 0.04 if site1.family == site2.family else 0.1
def hopping_color(site1,site2):
         return 'g' if site1.family==site2.family else 'g'
        

kwant.plot(sys,site_color=family_colors,site_lw=0.1, hop_lw=hopping_lw,hop_color=hopping_color,colorbar=False)

sys=sys.finalized()
#======================
def plot_conductance(sys, energies):
    
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(sys, energy)
        data.append(smatrix.transmission(2, 0)+smatrix.transmission(3, 0)+smatrix.transmission(2, 1)+smatrix.transmission(3, 1))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()
energies=[-3+i*0.02 for i in range(100)]
plot_conductance(sys,energies)