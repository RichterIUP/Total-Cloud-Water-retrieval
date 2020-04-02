#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:04:15 2019

@author: prowe

Purpose: Run DISORT using inputs specified in a namelist file

Instructions: 
    1) You will need to install f90nml:
       $ pip install f90nml
    2) Install rundisort_py
    3) Change the path below as needed for your setup.
    4) Make sure you have the file "sample.nml" in your path.
    5) Results should agree with those from sampleRun.
    
Copyright 2020 by Penny Rowe and NorthWest Research Associates     
"""

# .. Built-in modules
import numpy as np
#import f90nml


class Inputs(object):
    '''
    Get inputs from namelist file and put into class
    '''
    def __init__(self, namelist):
        
        # .. Load in the namelist
        nml = None#f90nml.read(namelist)
        input_nml = nml['inputs_nml']
        
        
        # .. Correct the liq ssp files
        ssp_temps = 'ssp_liq_temps'
        pmomfiles = 'pmomfiles_liq'
        input_nml[pmomfiles] = []
        if type(input_nml[ssp_temps]) == int:
            Ntemps = 1
            input_nml[ssp_temps] = [input_nml[ssp_temps]]
        else:
            Ntemps = len(input_nml[ssp_temps])
        for i in range(Ntemps):
            pmomfile_name = pmomfiles + str(i)
            input_nml[pmomfiles].append(input_nml[pmomfile_name])
            del input_nml[pmomfile_name]
        
        
        # .. Correct the ice ssp files
        ssp_temps = 'ssp_ice_temps'
        pmomfiles = 'pmomfiles_ice'
        input_nml[pmomfiles] = []
        if type(input_nml[ssp_temps]) == int:
            Ntemps = 1
            input_nml[ssp_temps] = [input_nml[ssp_temps]]
        else:
            Ntemps = len(input_nml[ssp_temps])
        for i in range(Ntemps):
            pmomfile_name = pmomfiles + str(i)
            input_nml[pmomfiles].append(input_nml[pmomfile_name])
            del input_nml[pmomfile_name]
        
                
        # .. Correct the microwindows
        dum = np.asarray(input_nml['microwins'])
        dum.resize(int(len(dum)/2),2)
        input_nml['microwins'] = dum
    
            
        # .. Turn dictionary into class
        for key in input_nml:
            setattr(self, key, input_nml[key])
            

