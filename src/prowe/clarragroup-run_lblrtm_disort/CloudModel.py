'''
CloudModel.py

@author: prowe

Copyright 2017-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

'''

# .. Built-in modules
import numpy as np

# .. My modules
import ssp_stuff as cloudx

class Cloud:
    __slots__ = (
      'base_liq', 'base_ice', 'top_liq', 'top_ice',
      'tempLiq', 'tempIce', 'temp',
      'layerLiq',                      # Atmospheric layer containing liquid
      'layerIce',                      # Atmospheric layer containing ice
      'layerCloud',                    # Atmospheric layer containing any
      'cloudLiqDist', 'cloudIceDist',
      'layers_liq_from_surface',
      'layers_ice_from_surface',
      'z',                            # Height from TOA to surface
      'dz',                           # Height differentials, TOA to surface
      'TEMPER',                       # Temperature from TOA to surface
      'modelName', 
      'height_CO2slicing', 'height_mlev',
      'base_liq_file', 'top_liq_file', 'base_ice_file', 'top_ice_file', 
      'microphysical_retrieval_uses_heights_from',
      'cloudMaskFileBef', 'cloudMaskFileAft', 
      'cloudPhaseFromTemp',
      'phase',
      'ssp_liqTemps', 'liqTdependence', 
      'ssp_iceTemps', 'iceTdependence',
      'cloudht_combineNus',
      'cloudht_chooseLyr',
      'iTempLiq', 'wTempLiq',
      'iTempIce', 'wTempIce',
      'wCloudLiqDist',             # weighted cloud distribution liq/ice:
      'wCloudIceDist',             # cloud layers x ssp files
      'iLiqLayer',                 # index to layerCloud with liquid
      'iIceLayer',                 # index to layerCloud with ice
      'iscloud',                   # True if cloud exists, 0 if not
      'hasLiq',
      'hasIce',
      )
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
     
    def __init__(self, ip):

        self.ssp_liqTemps = ip.ssp_liq_temps
        self.liqTdependence = ip.liq_tdependence
        self.ssp_iceTemps = ip.ssp_ice_temps
        self.iceTdependence = ip.ice_tdependence
        self.base_liq = ip.cloudbase_liq
        self.top_liq  = ip.cloudtop_liq
        self.base_ice = ip.cloudbase_ice
        self.top_ice  = ip.cloudtop_ice
        
                    
    def setCloudHeightsFromInputs(self, base_liq, base_ice, top_liq, top_ice):
        self.base_liq = base_liq
        self.base_ice = base_ice
        self.top_liq = top_liq
        self.top_ice = top_ice

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      

    def cloudModel(self, zs, Ts, thisdatetime,
                   ssp_liqTemps, liqTdependence, 
                   ssp_iceTemps, iceTdependence): #,
       #            cloudHtsFrom):
        """
        Inputs:
          self.base_liq,self.top_liq,self.base_ice,self.top_ice,
          self.z,
          self.TEMPER,
          self.modelName
        Outputs (as part of class):
          layerLiq,layerIce,cldTemp_wat,tempIce,
          
        Notes: cloud heights should already have been set
        """
    
        # .. Disort takes temperature from TOA to surface
        TEMPER = np.flipud(Ts)
        self.TEMPER = TEMPER.round(decimals=2)
        

        # Heights: make sure heights go from top down
        self.z = np.flipud(zs)
        self.dz = np.diff(self.z)
        if np.any(self.dz >= 0):
          raise NameError('One or more height does not decrease.')

        self.layerLiq = np.where(np.logical_and(self.z > self.base_liq,
                                                self.z < self.top_liq))[0]
        self.layerIce = np.where(np.logical_and(self.z > self.base_ice,
                                               self.z < self.top_ice))[0] 

        self.layerCloud = np.unique([self.layerLiq, self.layerIce])
        Ncloud = len(self.layerCloud)

        # .. Cloud temperatures are taken as layer mean temperature
        TEMPERmean = np.diff(TEMPER)/2 + TEMPER[:-1]
        cloud_temps = TEMPERmean[ self.layerCloud]

        if np.any(self.layerIce):
            Nssp_ice = len(ssp_iceTemps)
            self.wTempIce = np.zeros((Ncloud, Nssp_ice))
            self.wCloudIceDist = np.zeros((Ncloud, Nssp_ice))
            
            # .. Ice cloud: indices and weights to ssp files
            iTempIce, wTempIce = cloudx.interp1_weights2(
                ssp_iceTemps, cloud_temps, iceTdependence)
            
            # .. Set distributions of ice in clouds and normalize
            ice_dist = -self.dz[self.layerIce]
            self.cloudIceDist = ice_dist / np.sum(ice_dist) 

            self.iIceLayer = np.zeros(len(self.layerIce)).astype('int')
            for i,layer in enumerate(self.layerIce):
                self.iIceLayer[i] = self.layerCloud.tolist().index(layer)
                j = self.iIceLayer[i] 
                self.wTempIce[j,iTempIce[j,:]] = wTempIce[j,:]
                self.wCloudIceDist[j,iTempIce[j,:]] = \
                                 wTempIce[j,:] * -self.dz[self.layerCloud[j]]            
            # Normalize
            self.wCloudIceDist = self.wCloudIceDist \
                                  / np.sum(self.wCloudIceDist)
            self.iTempIce = iTempIce


        if np.any(self.layerLiq):
            Nssp_liq = len(ssp_liqTemps)
            self.wTempLiq = np.zeros((Ncloud, Nssp_liq))
            self.wCloudLiqDist = np.zeros((Ncloud, Nssp_liq))
            
            # .. Liq cloud: indices and weights to ssp files
            iTempLiq, wTempLiq = cloudx.interp1_weights2(
                ssp_liqTemps, cloud_temps, liqTdependence)
            
            # .. Set distributions of liq in clouds and normalize
            liq_dist = -self.dz[self.layerLiq]
            self.cloudLiqDist = liq_dist / np.sum(liq_dist) 

            self.iLiqLayer = np.zeros(len(self.layerLiq)).astype('int')
            for i,layer in enumerate(self.layerLiq):
                self.iLiqLayer[i] = self.layerCloud.tolist().index(layer)
                j = self.iLiqLayer[i] 
                self.wTempLiq[j,iTempLiq[j,:]] = wTempLiq[j,:]
                self.wCloudLiqDist[j,iTempLiq[j,:]] = \
                                  wTempLiq[j,:] * -self.dz[self.layerCloud[j]]            
                # Normalize
            self.wCloudLiqDist = self.wCloudLiqDist \
                                 / np.sum(self.wCloudLiqDist)
            self.iTempLiq = iTempLiq
        
 