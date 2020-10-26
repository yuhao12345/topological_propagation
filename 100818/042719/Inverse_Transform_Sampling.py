# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 21:54:50 2019

@author: ykang
"""
##Inverse Transform Sampling 
import numpy as np
import matplotlib.pyplot as plt


x=np.arccos(1-2*np.random.rand(50000))/2
n, bins, patches = plt.hist(x, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)