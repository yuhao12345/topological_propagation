# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 15:17:07 2018

@author: user
"""

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import *
from matplotlib import rcParams
from numpy import *
from numpy.linalg import *
import pickle
import sys
import os
import string
import heapq
import kwant
import tinyarray
from matplotlib import pyplot

chiral=True
if chiral:
    p = pi/5    #phi
    t = 0.66    #theta
    a = 0.34
    x = 1.4
    e1 = 0
    e2 = 0.3
    t2=0.1
    t1=-x*t2
    t0 = 2
    lam=-0.08
    t_so1 = 0.01 #spin-orbit coupling param
    t_so2 = x*t_so1 #spin-orbit coupling param
    tl=tr=0.3
    N = 30
    sigma_0 = tinyarray.array([[1, 0], [0, 1]])
    sigma_x = tinyarray.array([[0, 1], [1, 0]])
    sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
    sigma_z = tinyarray.array([[1, 0], [0, -1]])
    no=2            #number of orbitals
    def sigma_v1(ap):         # pauli metrix along the vertical axis
        value=sigma_z*cos(t)+sin(t)*(sigma_x*sin(ap)-sigma_y*cos(ap))
        return value

    def sigma_v2(ap):         # pauli metrix along the vertical axis
        value=sigma_z*cos(t)-sin(t)*(sigma_x*sin(ap)-sigma_y*cos(ap))
        return value

    def family_color(sites):
        return 'black' #if site.family == sites

    def hopping_lw(site1, site2):
        return 0.08


    class Amorphous(kwant.builder.SiteFamily):
        def __init__(self, coords):
            self.coords = coords
            super(Amorphous, self).__init__("amorphous", "",no)

        def normalize_tag(self, tag):
            try:
                tag = int(tag[0])
            except:
                raise KeyError

            if 0 <= tag < len(coords):
                return tag
            else:
                raise KeyError

        def pos(self, tag):
            return self.coords[tag]

    coords=[(0.0000000000, 0.0000000000, 0.0000000000), (-0.1336881039,
0.4114496766, 0.3400000000), (-0.4836881039, 0.6657395614, 0.6800000000),
(-0.9163118961, 0.6657395614, 1.0200000000), (-1.2663118961, 0.4114496766,
1.3600000000), (-1.4000000000, 0.0000000000, 1.7000000000), (-1.2663118961,
-0.4114496766, 2.0400000000), (-0.9163118961, -0.6657395614, 2.3800000000),
(-0.4836881039, -0.6657395614, 2.7200000000), (-0.1336881039,
-0.4114496766, 3.0600000000), (0.0000000000, -0.0000000000, 3.4000000000),
(-0.1336881039, 0.4114496766, 3.7400000000), (-0.4836881039, 0.6657395614,
4.0800000000), (-0.9163118961, 0.6657395614, 4.4200000000), (-1.2663118961,
0.4114496766, 4.7600000000), (-1.4000000000, 0.0000000000, 5.1000000000),
(-1.2663118961, -0.4114496766, 5.4400000000), (-0.9163118961,
-0.6657395614, 5.7800000000), (-0.4836881039, -0.6657395614, 6.1200000000),
(-0.1336881039, -0.4114496766, 6.4600000000), (0.0000000000, -0.0000000000,
6.8000000000), (-0.1336881039, 0.4114496766, 7.1400000000), (-0.4836881039,
0.6657395614, 7.4800000000), (-0.9163118961, 0.6657395614, 7.8200000000),
(-1.2663118961, 0.4114496766, 8.1600000000), (-1.4000000000, 0.0000000000,
8.5000000000), (-1.2663118961, -0.4114496766, 8.8400000000),
(-0.9163118961, -0.6657395614, 9.1800000000), (-0.4836881039,
-0.6657395614, 9.5200000000), (-0.1336881039, -0.4114496766, 9.8600000000),
(-1.4000000000, 0.0000000000, 0.0000000000), (-1.2663118961, -0.4114496766,
0.3400000000), (-0.9163118961, -0.6657395614, 0.6800000000),
(-0.4836881039, -0.6657395614, 1.0200000000), (-0.1336881039,
-0.4114496766, 1.3600000000), (0.0000000000, -0.0000000000, 1.7000000000),
(-0.1336881039, 0.4114496766, 2.0400000000), (-0.4836881039, 0.6657395614,
2.3800000000), (-0.9163118961, 0.6657395614, 2.7200000000), (-1.2663118961,
0.4114496766, 3.0600000000), (-1.4000000000, 0.0000000000, 3.4000000000),
(-1.2663118961, -0.4114496766, 3.7400000000), (-0.9163118961,
-0.6657395614, 4.0800000000), (-0.4836881039, -0.6657395614, 4.4200000000),
(-0.1336881039, -0.4114496766, 4.7600000000), (0.0000000000, -0.0000000000,
5.1000000000), (-0.1336881039, 0.4114496766, 5.4400000000), (-0.4836881039,
0.6657395614, 5.7800000000), (-0.9163118961, 0.6657395614, 6.1200000000),
(-1.2663118961, 0.4114496766, 6.4600000000), (-1.4000000000, 0.0000000000,
6.8000000000), (-1.2663118961, -0.4114496766, 7.1400000000),
(-0.9163118961, -0.6657395614, 7.4800000000), (-0.4836881039,
-0.6657395614, 7.8200000000), (-0.1336881039, -0.4114496766, 8.1600000000),
(0.0000000000, -0.0000000000, 8.5000000000), (-0.1336881039, 0.4114496766,
8.8400000000), (-0.4836881039, 0.6657395614, 9.1800000000), (-0.9163118961,
0.6657395614, 9.5200000000), (-1.2663118961, 0.4114496766, 9.8600000000)]
    amorphous_lat = Amorphous(coords)

    syst = kwant.Builder()

    #adding the onsite and hopping to the system
    for i in range(N):
        syst[amorphous_lat(i)] = e1*sigma_0
        syst[amorphous_lat(N+i)] = e2*sigma_0
        syst[amorphous_lat(i), amorphous_lat(N+i)] = lam*sigma_0
        if i > 0:
            syst[amorphous_lat(i), amorphous_lat(i-1)] = t1*sigma_0 +\
1j*t_so1*(sigma_v1(i*p)+sigma_v1((i-1)*p))
            syst[amorphous_lat(N+i),amorphous_lat(N+i-1)] = t2*sigma_0 +\
1j*t_so2*(sigma_v2(i*p)+sigma_v2((i-1)*p))
    # If we want to attach to vertical 1D chains to the system
    # we first add a site of the down lead to the scattering region
    lat=kwant.lattice.cubic(a, norbs=no)

    syst[lat(0, 0, -1)] = e1*sigma_0
    syst[amorphous_lat(0), lat(0, 0, -1)] = tl*sigma_0

    # We make a regular down lead and attach it to the system
    dn_lead = kwant.Builder(kwant.TranslationalSymmetry((0, 0, -a)))
    dn_lead[lat(0, 0, -2)] = e1*sigma_0
    dn_lead[lat.neighbors()] = t0*sigma_0
    syst.attach_lead(dn_lead)

    prim_vecs=tinyarray.array([(a,0.,0.),(0.,a,0.),(0.,0.,a)])
    offset=tinyarray.array((-1.2663118961, 0.4114496766,0.0))
    lat2=kwant.lattice.Monatomic(prim_vecs, offset, norbs=no)

    syst[lat2(0, 0, N)] = e1*sigma_0
    syst[amorphous_lat(2*N-1), lat2(0, 0, N)] = tr*sigma_0

    up_lead = kwant.Builder(kwant.TranslationalSymmetry((0, 0, a)))
    up_lead[lat2(0, 0, N+1)] = e1*sigma_0
    up_lead[lat2.neighbors()] = t0*sigma_0
    syst.attach_lead(up_lead)

    system=kwant.plot(syst, site_lw=0.1, site_color=family_color,hop_lw=hopping_lw,fig_size=(20, 15))


    trans=True
    if trans:
        syst = syst.finalized()
        energies = []
        data = []

        for ie in range(-320,520):
            energy = ie * 0.001
            smatrix = kwant.smatrix(syst, energy=energy)
            energies.append(energy)
            data.append(0.5*smatrix.transmission(1, 0))
        fig = pyplot.figure(figsize=(6,2))
        pyplot.plot(energies, data)
        pyplot.xlim([-0.32,0.52])
        pyplot.ylim([-0.03,1.25])
        pyplot.xlabel("energy [eV]")
        pyplot.ylabel("conductance [e^2/h]")
        pyplot.show()


    strans1=True
    if strans1:
        #syst = syst.finalized()
        energies = []
        data = []

        def oscle(ene, lead_nr):
            wfs=kwant.wave_function(syst, ene,check_hermiticity=True)(lead_nr)
            spin_current_z = 0
            for psi in wfs:
                psi_start = psi[0 : 2]
                psi_end = psi[2 * 61: 2 * 61 + 2]
                #spin_current_z += -2 *imag(psi_end.conjugate().dot(sigma_z).dot(psi_start))
                spin_current_z +=abs(imag(psi_end.conjugate().dot(sigma_z).dot(psi_start)))
            return spin_current_z
        for ie in range(-320,520):
            energy = ie * 0.001
            energies.append(energy)
            data.append(oscle(ene=energy, lead_nr=1))
        fig = pyplot.figure(figsize=(6,2))
        pyplot.plot(energies, data)
        pyplot.show()




    strans=True
    if strans:
        #syst = syst.finalized()
        J_spin = kwant.operator.Current(syst, sigma_z,
where=[(amorphous_lat(0), amorphous_lat(N-2))], sum=True)
        all_wfs = kwant.wave_function(syst, energy=0.25)(1)
        spin_current_list = sum(J_spin(wf) for wf in all_wfs)
        print(spin_current_list)