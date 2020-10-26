# -*- coding: utf-8 -*-
"""
Created on Sat Jun 23 20:37:04 2018

@author: user
"""

import numpy as np
from numpy import cos,sin,pi,arctan,sqrt
import matplotlib.pyplot as plt
#from random import random
import kwant
import scipy.sparse.linalg as sla

p1=0.5
p2=1.3
the=arctan(1/sqrt(3))
length=30
def onsite(site):
    x,y=site.pos
    return p1*(cos(2*x)+cos(2*y))+p2*(cos(2*(x*cos(the)-y*sin(the))+cos(2*(x*sin(the)+y*cos(the)))))

def disk(pos):
    x,y=pos
    return abs(x)<length and abs(y)<length

sys=kwant.Builder()
lat=kwant.lattice.square()
sys[lat.shape(disk,(0,0))]=onsite
sys[lat.neighbors()]= 1
sys=sys.finalized()

ham_mat = sys.hamiltonian_submatrix(sparse=True)
vals,evecs= sla.eigsh(ham_mat, k=1, sigma=None,which='LA', return_eigenvectors=True)  #,evecs

kwant.plotter.map(sys, np.abs(evecs[:, 0])**2,colorbar=False, oversampling=1)

