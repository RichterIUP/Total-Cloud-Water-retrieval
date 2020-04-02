
"""
run_radtran_micro.py

Copyright 2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose:
    Perform cloudy-sky radiative transfer for microwindows, using
    LBLRTM and DISORT
      
    Inputs:
        Namelist file: inputs.nml
    
    Outputs (where RU = mW/[m2 sr cm-1]):
        Variables:
            ip.microwins: microwindow lower and upper bound, cm-1
            micro.nu: microwindow effective wavenumber, cm-1
            Rcloudy: Cloudy-sky radiance computed by DISORT, RU
            Rclear_disort,clear: Clear-sky radiance computed by DISORT, RU
            micro.rad_clear: Clear-sky radiance computed using LBLRTM, RU
            
        Text file with columns:
            nu1 nu2 nu_eff DISORT,cloudy  DISORT,clear LBLRTM,clear

            nu1: microwindow lower bound, cm-1
            nu2: microwindow lower bound, cm-1
            nu_eff: microwindow effective wavenumber, cm-1
            DISORT,cloudy: Cloudy-sky radiance computed by DISORT, RU
            DISORT,clear: Clear-sky radiance computed by DISORT, RU
            LBLRTM,clear: Clear-sky radiance computed using LBLRTM, RU
            
        Optional: Optical depth file


    Classes:
        cloud
        inputsForDisort
        atmProf
        sspLiq
        sspIce


    Other Selected outputs of interest (not returned)
        cloud.z: Up/down flipped copy of atmospheric profile heights
        izm: wavenumber-depedendent index to cloud.z giving maximum height 
             of DISORT run (above this value gas optical depths are low)
        rfldn: "Diffuse down-flux (total minus direct-beam)
               (without delta-M scaling)." See DISORT.doc
        flup: "Diffuse up-flux." See DISORT.doc
        disort_msg: String containing DISORT warnings or errors
                    Note: this does not capture all possibilities!
                    Check the output in the terminal window.
        _clr: extension indicating clear sky run
        

    Note: 
        1) You will need to change the directories after "import site"
           below to those for your installation of rundisort_py.
           (I haven't figured out if/how I want to bundle this as a package)
        
        2) To get access to other outputs, classes, this code can easily be
          converted into a script.
   


    Example commands for running
        inputfile = '/Users/prowe/Git_repos/run_LBLRTM_DISORT/sample_run/' \
                    + 'inputs.nml'

        microwindows, \
        nu_eff, \
        Rcloudy, \
        Rclear_disort, \
        Rclear_lblrtm = run_lblrtm_disort_micro(inputfile)

       
"""    

import site
#site.addsitedir("/Users/prowe/Git_repos/rundisort_py/installation/")
site.addsitedir("./src/prowe/clarragroup-run_lblrtm_disort/")

# .. Built-in modules
import numpy as np
import datetime as dt
    
# .. My Modules
from Inputs import Inputs
from AtmosphericProfile import AtmProf
from CloudModel import Cloud
from MicroNus import MicroNus
from InputsForDisort import InputsForDisort
from ssp_stuff import get_ssp
from run_clear_sky_sim import run_clear_sky_sim
from run_disort import run_disort



def run_lblrtm_disort_micro(inputfile):
        
    # # #  # # # # #       DEBUGGING and PRINTING      # # # # # # # # # # # #
    print_stdout = True                # Print messages to screen?
    rerun = False                      # Re-run LBLRTM for gas optical depths
                                       # if files already exist?
    save_results = True                # Save results to text file?
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    
    
    # .. Set up
    ip = Inputs(inputfile)                                # Inputs from file
    #ip = inputfile

    date = dt.datetime(ip.date_year_month_day[0],         # Date (datetime obj)
                       ip.date_year_month_day[1], 
                       ip.date_year_month_day[2],
                       ip.time_hour_min_sec[0],
                       ip.time_hour_min_sec[1],
                       ip.time_hour_min_sec[2])
    outfile = ip.out_dir + 'radiance' \
              + date.strftime('%Y%m%d_%H%M%S') \
              + '.txt'                                      # Output file name      
    nu = np.arange(0., ip.ewn+ip.dnu, ip.dnu)               # nu, from 0
    nu = nu[nu >= ip.bwn-1e-6]                              # nu, from bwn
    micro = MicroNus(nu, ip.od_match, ip.microwins)         # Microwindow nus
    sspLiq = get_ssp(ip.pmomfiles_liq, micro.nu)            # SSP liquid
    sspIce = get_ssp(ip.pmomfiles_ice, micro.nu)            # SSP ice
    inputsForDisort = InputsForDisort(ip, micro.nu)         # DISORT inputs
    cloud = Cloud(ip)                                       # Cloud model
    
    
    
    # .. Get the clear sky (gas) variables, from LBLRTM
    #    Load in the atmospheric profile, 
    #    create the optical depth file, if indicated,
    #    and get the clear sky radiance from LBLRTM
    atmProf = AtmProf(ip.prof_dir + ip.prof_file)             # Atmos. profile
    if print_stdout: 
        print('Running LBLRTM to get clear sky variables.')
        print(ip.view_angle)
    clrSky = run_clear_sky_sim(ip.od_file, ip.prof_dir,       # Gas opt depths
                               ip.lblrtm_dir, ip.runlblrtm,
                               ip.alt, 
                               ip.z_trop, 
                               np.round(ip.bwn)-100, 
                               np.round(ip.ewn)+50, 
                               ip.npts_lores_ifgram,
                               ip.npts_padded_ifgram,
                               date, ip.dnu_lores, ip.prof_file,
                               np.ones(1)*ip.view_angle, 
                               atmProf.zs[-1], rerun, 
                               True, print_stdout)
    clrSky.restrict_nus(nu)
    
    # .. Check if wavenumbers agree to within 1e-6
    if (len(clrSky.nu) != len(nu)) \
      or (np.allclose(clrSky.nu, nu, rtol=1e-6) == False):
        raise NameError('Measured and simulated wavenumber vectors disagree!')
    
    
    micro.getSimValsInMicrowins(clrSky, ip.od_match)        # ODs in micronus
    cloud.setCloudHeightsFromInputs(ip.cloudbase_liq, 
                                    ip.cloudbase_ice, 
                                    ip.cloudtop_liq, 
                                    ip.cloudtop_ice)        # Set cloud heights
    
    
    # .. The cloud model
    cloud.cloudModel(atmProf.zs, atmProf.Ts, date,
                     ip.ssp_liq_temps, ip.liq_tdependence, 
                     ip.ssp_ice_temps, ip.ice_tdependence)    
                 
    # .. Inputs for DISORT
    inputsForDisort.setDisortVals(
       ip.lat, ip.lon, ip.alt, ip.view_angle, 
       date, micro.nu, micro.od_layer, 
       clrSky.rads.shape[1], atmProf.zs, cloud.TEMPER)
    
    
    
    # .. Get cloudy sky radiance from DISORT
    Rcloudy, izm, \
    rfldn, flup, \
    disort_msg = run_disort(inputsForDisort.nu, inputsForDisort.dtau_gas, 
                            cloud.layerCloud, 
                            ip.reff_wat, ip.cld_od_vis_liq,
                            ip.reff_ice, ip.cld_od_vis_ice,
                            cloud.iTempLiq, 
                            cloud.wCloudLiqDist * ip.cld_od_vis_liq, 
                            cloud.iTempIce, 
                            cloud.wCloudIceDist * ip.cld_od_vis_ice,
                            cloud.iLiqLayer, cloud.iIceLayer, 
                            sspLiq, sspIce, 
                            inputsForDisort.TEMPER, inputsForDisort.NSTR, 
                            inputsForDisort.UMU0, inputsForDisort.UMU,
                            inputsForDisort.albedo, inputsForDisort.FBEAM, 
                            inputsForDisort.delv, inputsForDisort.iobs)
    if print_stdout & disort_msg:
        print('for clear-sky DISORT run at ' 
              + date.strftime('%Y%m%d %H:%M:%S'))
    
    
    # .. Get clear sky radiance from DISORT
    Rclear_disort, \
    izm_clr, \
    rfldn_clr, \
    flup_clr, \
    disort_msg_clr = run_disort(inputsForDisort.nu, inputsForDisort.dtau_gas, 
                                np.array([len(atmProf.zs)-4]), 
                                10., 1e-6,
                                10., 0.0,
                                np.array([0]), np.array([[1.]]) * 1e-6,
                                np.array([0]), np.array([[0.]]),
                                np.array([0]), [],
                                sspLiq, sspIce, 
                                inputsForDisort.TEMPER, inputsForDisort.NSTR, 
                                inputsForDisort.UMU0, inputsForDisort.UMU,
                                inputsForDisort.albedo, inputsForDisort.FBEAM, 
                                inputsForDisort.delv, inputsForDisort.iobs)
    if print_stdout:
        if disort_msg:
            print('for clear-sky DISORT run at ' 
                  + date.strftime('%Y%m%d %H:%M:%S'))
        print('DISORT, cloudy:', Rcloudy)
        print('DISORT, clear:', Rclear_disort)
        print('LBLRTM:', micro.rad_clear)
        print('LBLRTM - DISORT, clear:' + str(Rclear_disort - micro.rad_clear))  
        
    
    if save_results:
        # .. Organize output and Save the results as an ASCII file
        X = np.vstack([micro.nu, Rcloudy, Rclear_disort, micro.rad_clear]).T
        X = np.hstack([ip.microwins, X])
        hdr = 'lower bound (cm-1))     upper bound (cm-1)       ' \
              + 'effective nu (cm-1) ' \
              + 'DISORT cloudy (RU = mW/[m2 sr cm-1]) '\
              + 'DISORT clear (RU)' \
              + 'LBLRTM clear (RU)  ' 
        np.savetxt(outfile, X, header = hdr)
    
    
    
    return ip.microwins, micro.nu, Rcloudy, Rclear_disort, micro.rad_clear
