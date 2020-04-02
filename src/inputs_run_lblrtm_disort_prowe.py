# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 12:51:20 2020

@author: Philipp
"""

import inp

# Inputs:
#   out_dir: Directory where output will be save, string
#   lblrtm_dir: Directory where LBLRTM is, string
#   runlblrtm: Command for running LBLRTM, string
#   prof_dir: Directory of profile file, string
#   prof_file: File with profile, string
#   od_file: Optical depth file, including directory (optional, if already created or if it will be saved), string     
#   solar_source_fun_file: file name for solar source function, string
#   surf_emiss_data_file: file name for surface emissivity data, string
#   pmomfiles_liq0: file name of ssp parameters for liquid for first temp, string
#   pmomfiles_liq1: file name of ssp parameters for liquid for second temp, string
#   pmomfiles_liq2: etc (up to the number of temps provided in ssp_liq_temps), string
#   pmomfiles_ice0: As above, for ice, string
#   date_year_month_day: year, month, day, list of int
#   time_hour_min_sec: hour, minute, second, list of int
#   lat: latitude, float
#   lon: longittude, float
#   alt: altitude in km, , float
#   view_angle: Instrument viewing angle (0=>zenith) in degrees, float
#   cld_od_vis_liq: Optical depth (referenced to visible wavelenths) of liquid cloud, float
#   cld_od_vis_ice: Optical depth (referenced to visible wavelenths) of ice cloud, float
#   reff_wat: Effective radius of liquid droplets in micron, float
#   reff_ice: Effective radius of ice droplets in micron, float
#   cloudbase_liq: Base of liquid cloud in km, float
#   cloudbase_ice: Base of ice cloud in km, float
#   cloudtop_liq: Top of liquid cloud in km, float
#   cloudtop_ice: Top of ice cloud in km, float
#   bwn: Beginning wavenumber, float
#   ewn: Ending wavenumber, float
#   dnu: Spacing of wavenumbers, float
#   npts_lores_ifgram: Number of points in measured interferogram, for resolution-matching
#   npts_padded_ifgram: Number of points in measured interferogram after zero-padding, for resolution-matching
#   liq_tdependence:  Handling of ssp temperature dependence for liquid, e.g. 'nearest' or 'interpolate'
#   ice_tdependence: Handling of ssp temperature dependence for ice, e.g. 'nearest' or 'interpolate'
#   ssp_liq_temps: temperatures corresponding to liquid ssp files (1 per file), list
#   ssp_ice_temps: temperatures corresponding to ice ssp files (1 per file), list
#   microwins: Beginning and ending wavenumbers for microwindows, list
#              e.g. first microwindow low, first microwindow high,
#                   second microwindow low, second microwindow high, etc
#   sim_toa: Maximum altitude of simulation, int
#   z_trop: Altitude of tropopause, float
#   od_match: Option for resolution matching, string, e.g. 'lowResMicrowin'
#   disort_nstr: Number of streams to use in DISORT, int
#   disort_lamber: Boolean for lambertian for DISORT, 0 or 1

overall_path = "/home/phi.richter"
out_dir = '{}/OUTFOLDER/'.format(overall_path)
prof_dir = '{}/atmospheric_profiles'.format(overall_path)
prof_file = 'prof20170701_0621.nc'
od_file = '{}/od20170701_0621.mat'.format(inp.PATH)       
lblrtm_dir = "{}/bin".format(inp.PATH_TO_LBLRTM)
runlblrtm = 'lblrtm'
solar_source_fun_file = '{}/TCWret/src/prowe/clarragroup-run_lblrtm_disort/input_files/kurucz.dat'.format(overall_path)
surf_emiss_data_file = '{}/TCWret/src/prowe/clarragroup-run_lblrtm_disort/input_files/EmisIceSnow1.txt'.format(overall_path)
pmomfiles_liq = ['{}/ssp_clarra/pmom_water_T240_RFN_S331.nc'.format(overall_path), \
                  '{}/ssp_clarra/pmom_water_T253_RFN_S331.nc'.format(overall_path), \
                  '{}/ssp_clarra/pmom_water_T263_RFN_S331.nc'.format(overall_path), \
                  '{}/ssp_clarra/pmom_water_T273_RFN_S331.nc'.format(overall_path), \
                  '{}/pmom_water_T300_DW_S331.nc'.format(overall_path)]
pmomfiles_ice = ['{}/ssp_clarra/pmom_ice_sphere_T266_Warren_S331.nc'.format(overall_path)]
alt = 0.02
view_angle = 0.0
cloudtop_ice = 1.0
bwn = 500.046104877
ewn = 1399.98446716
dnu = 0.1205291748046875
npts_lores_ifgram = 23697
dnu_lores = 0.3333333333333333
npts_padded_ifgram = 65536
liq_tdependence = 'interpolate'
ice_tdependence = 'nearest'   
ssp_liq_temps = [240, 253, 263, 273, 300]
ssp_ice_temps = [273]
microwins =     [770.9,  774.8,  785.9, 790.7, 809.5,  813.5,  817.0, 823.5, \
                 828.6,  834.6,  843.1, 848.1, 860.1,  864.0,  872.5, 877.5, \
                 891.9,  895.8,  898.2, 904.8, 929.6,  939.7,  958.0, 964.3, \
                 985.0,  991.5,  1076.6, 1084.8, 1092.2, 1098.1, 1113.6, 1116.6, \
                 1124.4, 1132.6, 1142.2, 1148.0, 1155.2, 1163.4]
z_trop = 2.
od_match = 'lowResMicrowin'
disort_nstr = 16
disort_lamber = 1
   
