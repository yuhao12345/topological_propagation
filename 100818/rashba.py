# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 21:11:23 2018

@author: user
"""
# https://www.mail-archive.com/kwant-discuss@kwant-project.org/msg01041.html

import kwant
from matplotlib import pyplot
import tinyarray
from numpy import sqrt
from math import *
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
Sigma_z = tinyarray.array([[1, 0], [0, -1]])

sys=kwant.Builder()
lat=kwant.lattice.honeycomb()
a,b= lat.sublattices
lx,ly=11,7
Lambda_R=2
Lambda=0.0
Lambda_nu=0#1.45

def Rectangle(pos,lx=lx,ly=ly):
    x,y=pos
    return 0<=x<lx and 1<=y<ly

def onsite(site,Lambda_nu=Lambda_nu):
	return Lambda_nu*sigma_0 if site.family == a else -Lambda_nu*sigma_0
           
def Hop_Value(site1,site2,Lambda=Lambda):
    return 1j*Lambda*(2./sqrt(3))*Sigma_z

def hop_ras_e1(site1,site2,Lambda_R=Lambda_R,t1=1):
    return t1*sigma_0-1j*Lambda_R*sigma_x

def hop_ras_e2(site1,site2,Lambda_R=Lambda_R,t1=1):
    return t1*sigma_0+1j*Lambda_R*(0.5*sigma_x - sqrt(3)/2.0*sigma_y)

def hop_ras_e3(site1,site2,Lambda_R=Lambda_R,t1=1):
    return t1*sigma_0+1j*Lambda_R*(0.5*sigma_x + sqrt(3)/2.0*sigma_y)

sys[lat.shape(Rectangle,(0,1))]=onsite   #even with 0 onsite potential, we need to use the zero matrix and not a scalar.

#sys[lat.neighbors()]=t1*sigma_0
sys[kwant.HoppingKind((0, 0), a,b)] = hop_ras_e1
sys[kwant.HoppingKind((0, 1), a,b)] = hop_ras_e2
sys[kwant.HoppingKind((-1, 1), a,b)] = hop_ras_e3
#Here we chose the clockwise Vij which are +1 (as retuned by Hop_value). The other ones (-1) will be given directly 
#by the hermeticity of the Hamiltonian
hoppings1 = (((1, 0), a, a), ((-1, 1), a, a), ((0, -1), a, a))    
hoppings2 = (((1, 0), b, b), ((-1, 1), b, b), ((0, -1), b, b))
sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1]] = Hop_Value
sys[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2]] = Hop_Value
#==========================left==================================
def lead_shape(pos,ly=ly):
    x,y=pos
    return    1<=y<ly

sym = kwant.TranslationalSymmetry((-1,0))
sym.add_site_family(lat.sublattices[0], other_vectors=[(-1, 2)])
sym.add_site_family(lat.sublattices[1], other_vectors=[(-1, 2)])
lead = kwant.Builder(sym)

lead[lat.shape(lead_shape, (0, 1))] = onsite

lead[kwant.HoppingKind((0, 0), a,b)] = hop_ras_e1
lead[kwant.HoppingKind((0, 1), a,b)] = hop_ras_e2
lead[kwant.HoppingKind((-1, 1), a,b)] = hop_ras_e3
lead[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1]] = Hop_Value
lead[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2]] = Hop_Value

sys.attach_lead(lead)
sys.attach_lead(lead.reversed())
#=======================right=============================
#def lead1_shape(pos,ly=ly):
#    x,y=pos
#    return    1<=y<ly
#    
#sym1 = kwant.TranslationalSymmetry((1,0))
#sym1.add_site_family(lat.sublattices[0], other_vectors=[(-1, 2)])
#sym1.add_site_family(lat.sublattices[1], other_vectors=[(-1, 2)])
#lead1 = kwant.Builder(sym1)
#
#lead1[lat.shape(lead1_shape, (0, 1))] = onsite
#
#lead1[kwant.HoppingKind((0, 0), a,b)] = hop_ras_e1
#lead1[kwant.HoppingKind((0, 1), a,b)] = hop_ras_e2
#lead1[kwant.HoppingKind((-1, 1), a,b)] = hop_ras_e3
#lead1[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings1]] = Hop_Value
#lead1[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings2]] = Hop_Value
#
#sys.attach_lead(lead1)
#======================================================
def family_colors(site):
        return 0 if (site.family == a)  else 1
def hopping_lw(site1, site2):
        return 0.04 if site1.family == site2.family else 0.1
def hopping_color(site1,site2):
         return 'g' if site1.family==site2.family else 'g'
        

kwant.plot(sys,site_color=family_colors,site_lw=0.1, hop_lw=hopping_lw,hop_color=hopping_color,colorbar=False)

sys=sys.finalized()

#pyplot.rcParams["figure.figsize"] = [20,15]
#kwant.plot(sys)
#=================================================
#def plot_conductance(sys, energies):
#    
#    data = []
#    for energy in energies:
#        smatrix = kwant.smatrix(sys, energy)
#        data.append(smatrix.transmission(1, 0))
#
#    pyplot.figure()
#    pyplot.plot(energies, data)
#    pyplot.xlabel("energy [t]")
#    pyplot.ylabel("conductance [e^2/h]")
#    pyplot.show()
#energies=[-3+i*0.02 for i in range(100)]
#plot_conductance(sys,energies)


