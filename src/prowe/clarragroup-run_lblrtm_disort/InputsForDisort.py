# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 17:36:32 2018

@author: prowe

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.
"""

# .. Built-in Modules
import numpy as np

# .. My modules
from sun_position import sun_position

class InputsForDisort:
    __slots__ = (
      'nu',
      'dtau_gas',
      'iobs',
      'NLYR',
      'saa',
      'UMU0',
      'UMU',
      'albedo',
      'FBEAM',
      'delv',
      'NSTR', 
      'LAMBER',
      'TEMPER',
      )
    
    #  'Sa',
    #  'Xa',
    
      
    def __init__(self, ip, mnu):
        #self.Sa                 = ip.Sa                     # if single phase?
        #self.Xa                 = ip.Xa
        self.NSTR          = ip.disort_nstr
        self.LAMBER        = ip.disort_lamber
        
        
        # # # # # # # # # #    Variables that depend on mnu   # # # # # # # # #
        # .. Assuming mnu never changes, get things that depend on it now
        # .. Surface albedo: Load from input file
        nu_surf_albedo, surf_albedo = \
            get_surface_albedo_from_file(ip.surf_emiss_data_file)

        self.albedo = np.interp(mnu, nu_surf_albedo, surf_albedo)
        if any([any(self.albedo < 0), any(self.albedo > 1), \
           any(np.isnan(self.albedo))]):
               raise NameError('bad albedo')
            
        # .. Solar spectrum
        self.FBEAM = getsolarbeam_IR(mnu, ip.solar_source_fun_file)
        self.delv = 0 * mnu + 1
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    def getAngleStuff(self, lat, lon, alt, viewAngle, thisdatetime):
        
        # .. Solar angles
        sun = sun_position(thisdatetime, lat, lon, alt)
        sza = sun.zenith
        
        # .. UMU: cosine of scene viewing angle, defined from TOA
        #         -1=>looking up (downwelling), 1=>looking down (upwelling)
        viewAngle_TOA = 180 - viewAngle
        
        # Assign values        
        self.saa = sun.azimuth
        self.UMU0 = np.cos(np.deg2rad(sza[0]))
        self.UMU = [np.cos(np.deg2rad(viewAngle_TOA))]
        
        

    def setDisortVals(self, lat, lon, alt, viewAngle, thisdatetime,
                      nu, dtau_gas, NLYR, zs, TEMPER):

        
        # .. Set solar and viewing angles for radiative transfer (DISORT)
        self.getAngleStuff(lat, lon, alt, viewAngle, thisdatetime)
         

        # .. If we have an off-zenith viewing angle,
        #    the optical depths were for the slant view,
        #    so we need to undo that so DISORT can run at
        #    the proper viewing angle relative to straight up/down
        if self.UMU != -1 and self.UMU != 1:
            dtau_gas = dtau_gas * np.abs(self.UMU[0])

        # .. Note for future: change the following to quality control
        #    the optical depths for values below 1e-5 and return
        #    wavenumber-dependent indices, which would then be
        #    read in to run_disort.py. This would be better because
        #    for multiple runs, it only needs to be done once
        #    But for now, just comment out the following, which
        #    no longer works since we have allowed a thin optical 
        #    depth between the surface and observer, which will be
        #    removed later
        #
        ## .. Remove optical depths <1e-5 below 12 km  
        ##    Above 12 km allow them; the atmosphere will be chopped
        ##    at the first occurrence of dtau_gas<1e-5
        #idtau_gas = np.where(zs[:-1]<12)[0]
        #count = 0
        #while np.any(dtau_gas[:,idtau_gas] < 1e-5):
        #    if count>len(nu):
        #        raise NameError('While loop continuing too long')
        #    iout = np.where(dtau_gas[count,idtau_gas]<1e-5)[0]
        #    dtau_gas[count,iout] = 1e-5
        #    count+=1

        self.nu          = (nu.T).astype('float64')
        self.dtau_gas    = np.flipud(dtau_gas.T)  # deep copies?
        self.NLYR        = NLYR
        self.TEMPER      = TEMPER

        # .. Set the index to the observer layer, iobs
        #    iobs corresponds to the number of layers above the observer. 
        #    So the surface => len(DTAUC), where 
        #    DTAUC is self.dtau_gas at some wavenumber j, 
        #    self.dtau_gas[:,j]
        alt_km = alt/1000
        if (zs[0] - alt_km) < 1e-4 :
            self.iobs = np.shape(self.dtau_gas)[0]  # => surface
        elif zs[0] < alt_km:
            self.iobs = np.where(np.flipud(zs) < alt_km)[0][0] - 1
        else:
            raise ValueError('Atmospheric profile above the observer altitude')

        
# # # # # # # #      LOAD SURFACE ALBEDO    # # # # # # # # # # # #  # # 
def get_surface_albedo_from_file(surfEmissDataFile):

    albedoData = np.loadtxt(surfEmissDataFile, comments='%')
    nu_surf_albedo = albedoData[:, 1]
    surf_albedo = 1-albedoData[:, 2]
    
    return nu_surf_albedo, surf_albedo

        
def getsolarbeam_IR (wnum=None, solarbeam_IR_file=None):

    kurucz = np.loadtxt(solarbeam_IR_file) #'kurucz.dat')
    beam = np.interp(wnum,kurucz[:,0],kurucz[:,1])/1000;
    return (beam)

        
