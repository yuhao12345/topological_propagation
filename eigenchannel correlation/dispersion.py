# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 01:25:24 2017

@author: user
"""

from types import SimpleNamespace

from ipywidgets import interact
import matplotlib
from matplotlib import pyplot
from mpl_toolkits import mplot3d
import numpy as np

import kwant
from wraparound import wraparound


def momentum_to_lattice(k):
    """Transform momentum to the basis of reciprocal lattice vectors.
    
    See https://en.wikipedia.org/wiki/Reciprocal_lattice#Generalization_of_a_dual_lattice
    """
    B = np.array(lat.prim_vecs).T
    A = B.dot(np.linalg.inv(B.T.dot(B)))
    return np.linalg.solve(A, k)


def dispersion_2D(syst, args=None, lim=1.5*np.pi, num_points=50):
    """A simple plot of 2D band structure."""
    if args is None:
        args = []
    momenta = np.linspace(-lim, lim, num_points)
    energies = []
    for kx in momenta:
        for ky in momenta:
            lattice_k = momentum_to_lattice([kx, ky])
            h = syst.hamiltonian_submatrix(args=(list(args) + list(lattice_k)))
            energies.append(np.linalg.eigvalsh(h))
    
    energies = np.array(energies).reshape(num_points, num_points, -1)
    emin, emax = np.min(energies), np.max(energies)
    kx, ky = np.meshgrid(momenta, momenta)
    fig = pyplot.figure()
    axes = fig.add_subplot(1, 1, 1, projection='3d')
    for band in range(energies.shape[-1]):
        axes.plot_surface(kx, ky, energies[:, :, band], cstride=2, rstride=2,
                          cmap=matplotlib.cm.RdBu_r, vmin=emin, vmax=emax,
                          linewidth=0.1)

l=1
lat = kwant.lattice.square()

bulk_graphene = kwant.Builder(kwant.TranslationalSymmetry(*lat.prim_vecs))
bulk_graphene[lat.shape((lambda pos: True), (0, 0))] = 4
bulk_graphene[lat.neighbors()] = -1
#bulk_graphene[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = -0j
#kwant.plot(bulk_graphene)
dispersion_2D(wraparound(bulk_graphene).finalized())
        
#zigzag_haldane = kwant.Builder(kwant.TranslationalSymmetry([l, 0]))
#zigzag_haldane[graphene.shape((lambda pos: abs(pos[1]) < 20), (0, 0))] = onsite
#zigzag_haldane[graphene.neighbors(1)] = 1
#zigzag_haldane[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 0.5j
#
#kwant.plotter.bands(zigzag_haldane.finalized())
