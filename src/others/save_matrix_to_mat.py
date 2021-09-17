# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 14:53:51 2017

@author: user
"""

import numpy as np
import scipy.io

# Some test data
x = ld

# Specify the filename of the .mat file
matfile = 'ldos80_0.3.mat'

# Write the array to the mat file. For this to work, the array must be the value
# corresponding to a key name of your choice in a dictionary
scipy.io.savemat(matfile, mdict={'t2': x}, oned_as='row')

# For the above line, I specified the kwarg oned_as since python (2.7 with 
# numpy 1.6.1) throws a FutureWarning.  Here, this isn't really necessary 
# since oned_as is a kwarg for dealing with 1-D arrays.

# Now load in the data from the .mat that was just saved
matdata = scipy.io.loadmat(matfile)

# And just to check if the data is the same:
assert np.all(x == matdata['t2'])