#!/usr/bin/python
'''@package docstring
Class to call LBLRTM and DISORT
'''

# -*- coding: utf8 -*-

import sys
import numpy             as np
sys.path.append("./src/")
sys.path.append("./src/prowe/clarragroup-run_lblrtm_disort/")
import aux2 as aux
import inp
import run_lblrtm_disort_micro
import inputs_run_lblrtm_disort_prowe as irldp

def forward_run(atmospheric_param, thread_fact, lblrtm, file_num):
    '''Start LBLRTM/DISORT
    
    @param atmospheric_param The current MCP
    @param thread_fact Factor need for the derivative. Can be +1 or -1
    @param lblrtm If true, LBLRTM will be executed
    @param file_num The number of the MCP, whose derivative shall be calculated
    '''

    '''
    Calculate the disturbance for each MCP
    '''
    drl = thread_fact[1] * np.array([0.0, aux.STEPSIZE, 0.0, 0.0, 0.0])
    dri = thread_fact[1] * np.array([0.0, 0.0, aux.STEPSIZE, 0.0, 0.0])
    dtl = thread_fact[1] * np.array([0.0, 0.0, 0.0, aux.STEPSIZE, 0.0])
    dti = thread_fact[1] * np.array([0.0, 0.0, 0.0, 0.0, aux.STEPSIZE])

    atm = [0.0, 0.0, 0.0, 0.0]
    atm[0] = atmospheric_param[0]+dtl[thread_fact[0]]
    atm[1] = atmospheric_param[1]+dti[thread_fact[0]]
    atm[2] = atmospheric_param[2]+drl[thread_fact[0]]
    atm[3] = atmospheric_param[3]+dri[thread_fact[0]]

    '''
    Set up LBLRTM and DISORT and run LBLRTM/DISORT
    '''

    #inputfile = Combine_Results(atm)
    inputfile = '/home/philipp/TCWret/src/prowe/clarragroup-run_lblrtm_disort/sample_run/inputs.nml'
    microwindows, nu_disort, Rcloudy, Rclear_disort, lblrtm_rad_clear = run_lblrtm_disort_micro.run_lblrtm_disort_micro(inputfile)
    if inp.SEARCH_INIT:
        aux.RADIANCE_SKIP[file_num] = np.array(Rcloudy)

    else:
        aux.WAVENUMBER_LBLDIS = nu_disort
        aux.RADIANCE_LBLDIS[file_num].append(np.array(Rcloudy))

    return

class Combine_Results():
    def __init__(self, atm_param):
        self.date_year_month_day = [int(aux.YEAR), int(aux.MONTH), int(aux.DAY)]
        self.time_hour_min_sec = [int(aux.HOUR), int(aux.MINUTE), 0]
        self.lat = aux.LAT
        self.lon = aux.LON     
        self.cld_od_vis_liq = atm_param[0] * (1 - atm_param[1])
        self.cld_od_vis_ice = atm_param[0] * atm_param[1]
        self.reff_wat = atm_param[2]
        self.reff_ice = atm_param[3]
        self.cloudbase_liq = aux.CLOUD_BASE[0]
        self.cloudbase_ice = aux.CLOUD_BASE[0]
        self.cloudtop_liq = aux.CLOUD_TOP[0]
        self.cloudtop_ice = aux.CLOUD_TOP[0]
        self.prof_file = aux.PROF_FILE

        
        self.out_dir = irldp.out_dir
        self.prof_dir = irldp.prof_dir
        self.od_file = irldp.od_file       
        self.lblrtm_dir = irldp.lblrtm_dir
        self.runlblrtm = irldp.runlblrtm
        self.solar_source_fun_file = irldp.solar_source_fun_file
        self.surf_emiss_data_file = irldp.surf_emiss_data_file
        self.pmomfiles_liq = irldp.pmomfiles_liq
        self.pmomfiles_ice = irldp.pmomfiles_ice
        self.alt = irldp.alt
        self.view_angle = irldp.view_angle
        self.bwn = irldp.bwn
        self.ewn = irldp.ewn
        self.dnu = irldp.dnu
        self.npts_lores_ifgram = irldp.npts_lores_ifgram
        self.dnu_lores = irldp.dnu_lores
        self.npts_padded_ifgram = irldp.npts_padded_ifgram
        self.liq_tdependence = irldp.liq_tdependence
        self.ice_tdependence = irldp.ice_tdependence   
        self.ssp_liq_temps = irldp.ssp_liq_temps
        self.ssp_ice_temps = irldp.ssp_ice_temps
        self.microwins = irldp.microwins
        self.z_trop = irldp.z_trop
        self.od_match = irldp.od_match
        self.disort_nstr = irldp.disort_nstr
        self.disort_lamber = irldp.disort_lamber
