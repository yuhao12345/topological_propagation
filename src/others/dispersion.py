# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 01:25:24 2017

@author: user
"""
# haldane model: bulk and edge spectrum

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
    B = np.array(graphene.prim_vecs).T
    A = B.dot(np.linalg.inv(B.T.dot(B)))
    return np.linalg.solve(A, k)


def dispersion_2D(syst, args=None, lim=2*np.pi, num_points=50):
    """A simple plot of 2D band structure."""
    global energies,kx
    if args is None:
        args = []
    momenta0=[0]
    momenta = np.linspace(-lim, lim, num_points)
    energies = []
    for kx in momenta:
        for ky in momenta0:
            lattice_k = momentum_to_lattice([kx, ky])
            h = syst.hamiltonian_submatrix(args=(list(args) + list(lattice_k)))
            energies.append(np.linalg.eigvalsh(h))
    
#    energies = np.array(energies).reshape(num_points, num_points, -1)
    energies = np.array(energies).reshape( num_points,1, -1)
    emin, emax = np.min(energies), np.max(energies)
    kx, ky = np.meshgrid(momenta, momenta0)
#    fig = pyplot.figure()
#    axes = fig.add_subplot(1, 1, 1, projection='3d')
#    for band in range(energies.shape[-1]):
##        axes.plot_surface(kx, ky, energies[:, :, band], cstride=2, rstride=2,
##                          cmap=matplotlib.cm.RdBu_r, vmin=emin, vmax=emax,
##                          linewidth=0.1)
#        axes.scatter(kx, ky, energies[:, :, band])
#    axes.view_init(0, azim=0)
    
    pyplot.plot(np.transpose(kx),energies[:,0,:])

l=1
graphene = kwant.lattice.general([[l, 0], [l/2, l*np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, l/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices

nnn_hoppings_a = (((-1, 0), a, a), ((0, 1), a, a), ((1, -1), a, a))
nnn_hoppings_b = (((1, 0), b, b), ((0, -1), b, b), ((-1, 1), b, b))
nnn_hoppings = nnn_hoppings_a + nnn_hoppings_b

def onsite(site):
    x,y=site.pos
    return .5*(1 if site.family == a else -1)

#bulk_graphene = kwant.Builder(kwant.TranslationalSymmetry(*graphene.prim_vecs))
#bulk_graphene[graphene.shape((lambda pos: True), (0, 0))] = 0
#bulk_graphene[graphene.neighbors()] = 1
#bulk_graphene[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = -0j
##kwant.plot(bulk_graphene)
#dispersion_2D(wraparound(bulk_graphene).finalized())
        
zigzag_haldane = kwant.Builder(kwant.TranslationalSymmetry([l, 0]))
zigzag_haldane[graphene.shape((lambda pos: abs(pos[1]) < 9), (0, 0))] = onsite  #up and down edge should be zigzag
zigzag_haldane[graphene.neighbors(1)] = -1
zigzag_haldane[[kwant.builder.HoppingKind(*hopping) for hopping in nnn_hoppings]] = 0.15j

kwant.plotter.bands(zigzag_haldane.finalized())
