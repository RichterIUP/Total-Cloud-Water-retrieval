# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 09:57:17 2018

@author: prowe

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

"""

import numpy as np
import copy


class MicroNus:
    __slots__ = (
      'nu',                # microwindo wavenumber (nu)
      'Nmicnus',           # number of micronus
      'inus',              # indices to measurement nus
      'Nnus',              # number of measurement nus per micronu
      'inu_cloudcheck',    # index to measurement nu for checking for cloud
      'Rmeas',             # measured radiance
      'Rsigma',            # std dev measured radiance
      'tsc',               # transmittance, surface-to-layer
      'od_layer',          # optical depth of layer
      'rad_above',         # radiance above troposphere x trans of trop
      'rad_at_top',        # radiance impinging on troposphere top
      'rad_clear',         # clear-sky radiance
      'rads',              # layer radiances up to troposphere
     )


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def __init__(self, mrad_nu, odMatch, microwins):
      
        #% % % %    Set up wavenumbers    % % % % % % % % %
        # .. If the wavenumbers do not change, this only needs to happen
        #    once. But if they change, it must be run again 
        
        # .. Get indices to measured/calculated radiances to use
        nu_microlist = microwins
        nu_microlist = nu_microlist.astype(float)
        self.Nmicnus = nu_microlist.shape[0]
        self.nu = np.zeros(self.Nmicnus)
        if odMatch == 'dum':
          self.inus = np.zeros(self.Nmicnus, dtype='int')
          for inu in range(self.Nmicnus):  
        
            inus0 = np.where(np.logical_and(mrad_nu>=microwins[inu,0:1], \
                                          mrad_nu<=microwins[inu,1:2]))[0]+1
            if len(inus0)==1:
              self.inus[inu] = inus0  # check this agrees with matlab
            else:
              # This is a circuitous way of copying what I did in matlab
              self.inus[inu] = inus0[int(np.ceil( (len(inus0)+0.0)/2))] \
                                      - 1 -1
           
            self.nu[inu] = mrad_nu[self.inus[inu]]
          
        elif odMatch == 'lowResTrans1pt':   
          # .. Use one point at the average nu for each microwindow
          #    this choice is probably obsolete.
          self.inus  = np.zeros(self.Nmicnus,dtype='int')
          self.Nnus = np.zeros(self.Nmicnus, dtype='int')
        
          for inu in range(self.Nmicnus):  
            inus0 = np.where(np.logical_and(mrad_nu>=microwins[inu,0:1], \
                                            mrad_nu<=microwins[inu,1:2]))[0]
            self.Nnus[inu] = int(len(inus0))
                        
            self.nu[inu] = np.mean(mrad_nu[inus0])
            try:
                self.inus[inu] = int(np.floor(np.mean(inus0)))
            except:
                raise ValueError('Something is wrong')
        elif odMatch == 'lowResMicrowin':   
          # .. We will average all points in each microwindow
          #    Here we get the indices for each microwindow
          #    for the measurement
          self.inus  = np.zeros((self.Nmicnus, 2),dtype='int')
          self.Nnus = np.zeros(self.Nmicnus, dtype='int')
          for imic in range(self.Nmicnus):          
            inus0 = np.where(np.logical_and(mrad_nu >= microwins[imic,0:1], \
                                            mrad_nu <= microwins[imic,1:2]))[0]
            self.Nnus[imic] = int(len(inus0))
            try:
                self.inus[imic,0] = inus0[0]
            except:
                raise NameError('No wavenumbers found in microwindow')
            self.inus[imic,1] = inus0[-1]+1  # Add one b/c of Python indexing
            self.nu[imic] = np.mean(mrad_nu[inus0]) 
        else:
          raise NameError('Bad value for ip.odMatch:'+odMatch)
          
        self.inu_cloudcheck = np.where(mrad_nu > 960)[0][0]
        
        

    def getMeasRadInMicrowins(self, mrad_rad, mrad_stdDev, odMatch):
    
        # .. If we are going to retrieve cloud microphysical properties,
        #    we will need the measurements in the microwindows.
        #    It is better to do this here rather than within the loop
        #    over ispec, to save time
        if odMatch == 'dum':                    # 'lowResTrans1pt':
            self.Rmeas = mrad_rad[self.inus,:]
            self.Rsigma = mrad_stdDev[self.inus]
        elif odMatch =='lowResTrans1pt':         
            self.Rmeas  = mrad_rad[self.inus,:]
            self.Rsigma = mrad_stdDev[self.inus]
        elif odMatch == 'lowResMicrowin':
            # .. We will average all points in each microwindow
            #    We already got the indices for each microwindow
            #    for the measurement (during the init, above)
            Nspecs = mrad_rad.shape[1]
            self.Rmeas = np.zeros((self.Nmicnus, Nspecs))
            self.Rsigma = np.zeros((self.Nmicnus, Nspecs))
            for inu in range(self.Nmicnus):
                self.Rmeas[inu,:] = np.mean(\
                      mrad_rad[self.inus[inu][0]:self.inus[inu][1],:], \
                      axis=0)
                self.Rsigma[inu,:] = np.mean(\
                      mrad_stdDev[self.inus[inu][0]:self.inus[inu][1],:], \
                      axis=0)
                    
        else:
            raise NameError('Bad value for odMatch: ' + odMatch)
    

    def getSimValsInMicrowins(self, clrSky, odMatch):
        #Get wavenumbers, radiances, standard deviations, and
        #    optical depths, in microwindows
        if odMatch == 'lowResMicrowin':
            nlyr = clrSky.tsc.shape[1]
            nlyr_trop = clrSky.rads.shape[1]
            # .. Get the mean transmittance in each microwindow
            self.tsc = np.zeros((self.Nmicnus, nlyr))
            self.rad_clear = np.zeros(self.Nmicnus)
            self.rad_above = np.zeros(self.Nmicnus)
            self.rad_at_top = np.zeros(self.Nmicnus)
            self.rads = np.zeros((self.Nmicnus, nlyr_trop))
            self.od_layer = np.zeros((self.Nmicnus, nlyr))
            for inu in range(self.Nmicnus):
                inus = [i for i in range(self.inus[inu][0],self.inus[inu][1])]
                self.tsc[inu,:] = np.mean(clrSky.tsc[inus,:], axis=0)
                self.rad_clear[inu] = np.mean(clrSky.rad_clear[inus])
                self.rad_above[inu] = np.mean(clrSky.rad_above[inus])
                self.rads[inu,:] = np.mean(clrSky.rads[inus,:], axis=0)

                try:
                    self.rad_at_top[inu] = np.mean(clrSky.rad_at_top[inus])
                except:
                    pass
                

            # .. Get the layer optical depths corresponding to the 
            #    transmittances and get the radiances above
            odtot_Lm1 = np.zeros(self.Nmicnus)
            for ilayer in range(nlyr):
                # .. Layer effective-resolution optical depth
                trans0 = copy.deepcopy(self.tsc[:,ilayer])
                trans0[trans0>=1] = 1
                trans0[trans0<=0] = 1e-40
                odtot = -np.log(trans0)
                self.od_layer[:,ilayer] = odtot - odtot_Lm1
                                
                # .. Set up for next time around
                odtot_Lm1 = odtot
            
        elif odMatch == 'lowResTrans1pt':   
            nlyr = clrSky.tsc.shape[1]
            nlyr_trop = clrSky.rads.shape[1]
            # .. Get the transmittance in each microwindow
            self.tsc = clrSky.tsc[self.inus,:]
            self.rad_clear = clrSky.rad_clear[self.inus]
            self.rads = clrSky.rads[self.inus,:]
            self.od_layer = np.zeros((self.Nmicnus, nlyr))
                

            # .. Get the layer optical depths corresponding to the 
            #    transmittances and get the radiances above
            odtot_Lm1 = np.zeros(self.Nmicnus)
            for ilayer in range(nlyr):
                # .. Layer effective-resolution optical depth
                trans0 = copy.deepcopy(self.tsc[:,ilayer])
                trans0[trans0>=1] = 1
                trans0[trans0<=0] = 1e-40
                odtot = -np.log(trans0)
                self.od_layer[:,ilayer] = odtot - odtot_Lm1
                                
                # .. Set up for next time around
                odtot_Lm1 = odtot
            
        else:
            raise NameError('Bad value for odMatch: ' + odMatch)
