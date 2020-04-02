# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 17:36:32 2018

@author: prowe

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

"""

# .. Built in modules
from netCDF4 import Dataset
import numpy as np


class AtmProf:
    NoOfInstances = 0
    __slots__ = (
      'Nodsinit',
      'NLYR',
      'prof_file',
      'zs',
      'Ps',
      'Ts',
      'h2os',
      'h2o_mb',
      )
     
    
    def __init__(self, prof_file):
        
        # .. Load atmospheric profile
        self.prof_file = prof_file
        
        # .. Load in the profile
        with Dataset(prof_file, 'r', format='NETCDF4_CLASSIC') as nc:    
            self.zs = np.double(nc['z'][:].data)
            self.Ps = np.double(nc['P'][:].data)
            self.Ts = np.double(nc['T'][:].data)
            self.h2os = np.double(nc['h2o'][:].data)
                
        
        # .. Get water vapor in mb (not used yet?)
        h2o_ppp = 1e-6*self.h2os 
        h2o1_mb = h2o_ppp * self.Ps                 # mb x 10^6
        h2o1_mb = h2o_ppp * (self.Ps-h2o1_mb)       # mb x 10^6
        self.h2o_mb = h2o_ppp * (self.Ps-h2o1_mb)   # mb x 10^6
        
        
    
