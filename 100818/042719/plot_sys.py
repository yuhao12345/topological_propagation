# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:10:00 2019

@author: ykang
"""
import kwant
import stru_QSH as dg
import numpy as np
width=2.3
length=4
salt='2'
dis=2.4
graphene = kwant.lattice.general([[1, 0], [1/2, np.sqrt(3)/2]],  # lattice vectors
                                 [[0, 0], [0, 1/np.sqrt(3)]])  # Coordinates of the sites
a, b = graphene.sublattices
#sys=dg.make_system_all(width, length, salt,dis)
syst=dg.make_system(width, length, salt,dis) 
dg.attach_lead(syst,width)
def family_color(site):
        return 'black' if site.family == a else 'white'

def hopping_lw(site1, site2):
    return 0.04 if site1.family == site2.family else 0.1

kwant.plot(syst, site_lw=0.1, site_color=family_color, hop_lw=hopping_lw,fig_size=(15, 10))
#kwant.plot(sys)