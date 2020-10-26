# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 13:38:37 2017

@author: user
"""
import kwant
import inspect
lines = inspect.getsourcelines(kwant.physics.modes)
print("".join(lines[0]))

