#!/usr/bin/python
'''@package docstring
Class to call LBLRTM and DISORT
'''

# -*- coding: utf8 -*-

import subprocess
import sys
import numpy             as np
import netCDF4           as nc
import scipy.signal      as sig
import datetime          as dt
sys.path.append("./src/")
import aux2 as aux
import inp
import rundecker         as rd
import lblrun

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
    Set up LBLDIS and run LBLRTM/DISORT
    '''
    
    if inp.SEARCH_INIT:
        forward = LBLDIS(atm, file_num)
    else:
        forward = LBLDIS(atm, int(thread_fact[0]*thread_fact[1]))
    radiance = forward.execute(lblrtm=lblrtm)
    
    '''
    Save the calculated radiance
    '''
    #if inp.SEARCH_INIT:
    #    #aux.RADIANCE_LBLDIS[file_num] = np.array(radiance)
    #    aux.RADIANCE_SKIP[file_num] = np.array(radiance)
    #elif inp.SKIP_SN:
    #    aux.RADIANCE_SKIP[file_num] = np.array(radiance)
    #else:
    aux.RADIANCE_LBLDIS[file_num].append(np.array(radiance))

    return

####################################################################################
    
class LBLDIS:
    '''Class to setup and run LBLDIS
    '''

    ####################################################################################

    def __init__(self, mcp, thread):
        '''Read the MCP
        
        @param self self pointer
        @param mcp The current MCP
        @param thread The number of the thread, needed for consistent execution of the derivatives
        '''
        self.__radio = {"height": aux.ATMOSPHERIC_GRID[aux.ALT], \
                        "pressure": aux.ATMOSPHERIC_GRID[aux.PRES], \
                        "temperature": aux.ATMOSPHERIC_GRID[aux.TEMP], \
                        "humidity": aux.ATMOSPHERIC_GRID[aux.HUM]}
        self.__r_eff_liq = mcp[2]
        self.__r_eff_ice = mcp[3]
        self.__tau_total = mcp[0]
        self.__f_ice = mcp[1]
        self.__thread = thread
        

    ####################################################################################

    def run_lblrtm(self):
        '''Set up and execution of LBLRTM
        
        @param self self pointer
        '''

        '''
        Use the correct unit for the humidity
        '''
        if inp.HUMIDITY == 'g/kg':
            hmd_unit = 'C'
        elif inp.HUMIDITY == 'g/m3':
            hmd_unit = 'D'
        elif inp.HUMIDITY == 'ppmv':
            hmd_unit = 'A'
        else:
            raise ValueError("Need unit for humidity")
            
        '''
        Write the configuration file for LBLRTM
        '''
        rd.rundecker(self.__radio['height'], self.__radio['pressure'], \
                     self.__radio['temperature'], self.__radio['humidity'], \
                     "{}".format(aux.LBLTP5), co2_ppm=aux.CO2_PROFILE, \
                     o3_ppm=aux.O3_PROFILE, atm=inp.PREDEF_ATM, hmd_unit=hmd_unit, \
                     sample = 4.0)

        '''
        Run LBLRTM
        '''
        print("[{}] Run LBLRTM...".format(dt.datetime.now()))
        lblrun.lblrun(aux.LBLTP5, aux.LBLDIR, aux.LBLLOG, lbl_home=inp.PATH_TO_LBLRTM)
        print("[{}] LBLRTM finished...".format(dt.datetime.now()))

        return

    ####################################################################################

    def write_lbldis_parm(self):
        '''Create the configuration file for LBLDIS
        
        @param self self pointer
        '''
        n_layer_liq = len(aux.CLOUD_GRID)
        n_layer_ice = len(aux.CLOUD_GRID)
        cloud_grid = []
        for layer in sorted(aux.CLOUD_GRID):
            cloud_grid.append(layer)
            
            
        if inp.ONLY_OD or True:
            tau_ice = self.__tau_total * self.__f_ice
            tau_liquid = self.__tau_total * (1 - self.__f_ice)
        if tau_liquid < 0.0:
            tau_liquid = 0.0
        if tau_ice < 0.0:
            tau_ice = 0.0

        wavenumber_low = aux.MICROWINDOWS[inp.WINDOWS[0]][0]-50.0
        wavenumber_high = aux.MICROWINDOWS[inp.WINDOWS[-1]][-1]+50.0

        with open("{}/.lbldis_{}.parm".format(inp.PATH, self.__thread), "w") as file_:
            file_.write("LBLDIS parameter file\n")
            file_.write("16		Number of streams\n")
            file_.write("{:04.1f} 30. 1.0	Solar ".format(aux.SOLAR_ZENITH_ANGLE))
            file_.write("zenith angle (deg), relative azimuth (deg), solar distance (a.u.)\n")
            file_.write(" 180           Zenith angle (degrees): 0 -> ")
            file_.write("upwelling, 180 -> downwelling\n")
            file_.write(" {} {} {}".format(wavenumber_low, \
                        wavenumber_high, inp.RESOLUTION))
            file_.write(" v_start, v_end, and v_delta [cm-1]\n")
            file_.write("1               Cloud parameter option flag: ")
            file_.write("0: reff and numdens, >=1:  reff and tau\n")
            if inp.BUILTIN:
                file_.write("{}               ".format(n_layer_liq+n_layer_ice))
                file_.write("Number of cloud layers\n")
            else:
                file_.write("{}".format(n_layer_liq+(int(inp.NUM_OF_SHAPES*n_layer_ice))))
                file_.write("               Number of cloud layers\n")
            for loop_liq_layer in range(n_layer_liq):
                alt = cloud_grid[loop_liq_layer]*1e-3
                for loop_radio_height in range(len(self.__radio['height'])):
                    if(self.__radio['height'][loop_radio_height] > alt - 1e-3 and \
                       self.__radio['height'][loop_radio_height] < alt + 1e-3):
                        ind = loop_radio_height
                        break
                temp_of_layer = self.__radio['temperature'][ind]
                if temp_of_layer < 240 + (253 - 240)/2.0:
                    mie_liq = 12
                elif temp_of_layer <= 253 + (263 - 253)/2.0:
                    mie_liq = 13
                elif temp_of_layer <= 263 + (273 - 263)/2.0:
                    mie_liq = 14
                else:
                    mie_liq = 0
                tau_liq_lay = tau_liquid/float(n_layer_liq)
                reff_liq = self.__r_eff_liq
                file_.write("{} {:5.3f} {:10.8f} -1 {:10.8f}\n".format(mie_liq, \
                                                                       alt, reff_liq, tau_liq_lay))
            for loop_ice_layer in range(n_layer_ice):
                alt = cloud_grid[loop_ice_layer]*1e-3
                reff_ice = self.__r_eff_ice
                tau_ice_lay = tau_ice/float(n_layer_ice)
                if inp.BUILTIN:
                    if self.__r_eff_ice <= 16:
                        mie_ice = 4
                    else:
                        mie_ice = 7
                    file_.write("{} {:5.3f} {:10.8f} -1 {:10.8f}\n".format(mie_ice, \
                                                                           alt, reff_ice, tau_ice_lay))
                else:
                    for loop_ice_shape in range(len(inp.SHAPES)):
                        if inp.SHAPES[loop_ice_shape] > 0.0:
                            mie_ice = loop_ice_shape+1
                            tau_ice_lay = inp.SHAPES[loop_ice_shape] * \
                            tau_ice/float(n_layer_ice)
                            file_.write("{} {:5.3f} ".format(mie_ice, alt))
                            file_.write("{:10.8f} -1 {:10.8f}\n".format(reff_ice, tau_ice_lay))
            file_.write("{}\n".format(aux.LBLDIR))
            file_.write("solar/solar.kurucz.rad.1cm-1binned.full_disk.asc\n")
            file_.write("15       Number of scattering property databases\n")
            file_.write("ssp/ssp_db.mie_wat.gamma_sigma_0p100\n")
            file_.write("ssp/ssp_db.mie_ice.gamma_sigma_0p100\n")#1
            file_.write("ssp/ssp_db.Aggregate.gamma.0p100\n")#2
            file_.write("ssp/ssp_db.BulletRosette.gamma.0p100\n")#3
            file_.write("ssp/ssp_db.Droxtal.gamma.0p100\n")#4
            file_.write("ssp/ssp_db.HollowCol.gamma.0p100\n")#5
            file_.write("ssp/ssp_db.Plate.gamma.0p100\n")#6
            file_.write("ssp/ssp_db.SolidCol.gamma.0p100\n")#7
            file_.write("ssp/ssp_db.Spheroid.gamma.0p100\n")#8
            file_.write("ssp/ssp_db.mie_gypsum.lognormal_sigma_0p699\n")#9
            file_.write("ssp/ssp_db.mie_kaolinite.lognormal_sigma_0p699\n")#12
            file_.write("ssp/ssp_db.mie_quartz.lognormal_sigma_0p699\n")#11
            file_.write("ssp/ssp_db.mie_wat_zasetsky240.gamma_sigma_0p100\n")#12
            file_.write("ssp/ssp_db.mie_wat_zasetsky253.gamma_sigma_0p100\n")#13
            file_.write("ssp/ssp_db.mie_wat_zasetsky263.gamma_sigma_0p100\n")#14
            file_.write("-1.	Surface temperature (specifying a negative")
            file_.write("value takes the value from profile)\n")
            file_.write("4	Number of surface spectral emissivity lines (wnum, emis)\n")
            file_.write("100 {}\n".format(inp.EMISSIVITY))
            file_.write("700 {}\n".format(inp.EMISSIVITY))
            file_.write("800 {}\n".format(inp.EMISSIVITY))
            file_.write("3000 {}\n".format(inp.EMISSIVITY))
        return

    ####################################################################################

    def run_disort(self):
        '''Start DISORT
        
        @param self self pointer
        '''
        lbldisout_file = '{}/.lbldisout_{}'.format(inp.PATH, self.__thread)
        lbldislog = '{}/.lbldislog_{}.txt'.format(inp.PATH, self.__thread)
        with open("{}/.run_disort_{}.sh".format(inp.PATH, self.__thread), "w") as file_:
            file_.write("#!/bin/bash\n")
            exec_lbldis = '({}/lbldis {}/.lbldis_{}.parm 0 {}) >& {}\n'
            file_.write(exec_lbldis.format(inp.PATH_TO_LBLDIS, inp.PATH, \
                                           self.__thread, lbldisout_file, lbldislog))
        print("[{}] Run DISORT ...".format(dt.datetime.now()))
        subprocess.call(["bash", "{}/.run_disort_{}.sh".format(inp.PATH, self.__thread)])
        print("[{}] DISORT finished ...".format(dt.datetime.now()))

        return

    ####################################################################################

    def write_data_to_file(self):
        '''Read the calculated radiances
        
        @param self self pointer
        
        @return The calculated radiances
        '''
        wavenumber = []
        radiance = []
        wn_interpolated = []
        radiance_interpolated = []
        radiance_out = []
        disort_out = nc.Dataset('{}/.lbldisout_{}.cdf'.format(inp.PATH, self.__thread))
        
        '''
        Read the radiances
        '''
        for i in range(len(disort_out.variables['wnum'])):
            #if aux.in_windows(disort_out.variables['wnum'][i], inp.WINDOWS):
            wavenumber.append(disort_out.variables['wnum'][i])
            radiance.append(disort_out.variables['radiance'][i][0])

        '''
        Convolve the radiances with a boxcar function
        '''
        if inp.CONVOLVE and inp.RESOLUTION < 2.0:
            ft_boxcar = lambda xx: 2 * (inp.OPD) * np.sinc(2 * (inp.OPD) * xx)
            radiance = np.array(radiance)
            wavenumber = np.array(wavenumber)
            x_axis = np.array([loop_count for \
                               loop_count in range(len(wavenumber))])
            convolution = sig.convolve(radiance, \
                                       ft_boxcar(x_axis), mode='full')
            convolution = np.array(convolution[:len(radiance)])
            normalisation = max(radiance)/max(convolution)
            radiance = normalisation * convolution

        '''
        If necessery, interpolate the calculated radiances to the FTIR grid
        '''
        for element in aux.WAVENUMBER_FTIR:
            try:
                if(element > wavenumber[0] and element < wavenumber[-1]):
                    radiance_interpolated.append(np.interp(element, wavenumber, radiance))
                    wn_interpolated.append(element)
            except IndexError:
                sys.stderr.write("{}\n".format(wavenumber))
                sys.exit(-1)
        '''
        Discard all values outside the microwindows
        '''
        number_of_datapoints = len(wn_interpolated)
        for i in range(number_of_datapoints):
            if aux.in_windows(wn_interpolated[i], inp.WINDOWS):
                radiance_out.append(radiance_interpolated[i])

        return radiance_out

    ####################################################################################

    def execute(self, lblrtm=False):
        '''Execute LBLDIS
        
        @param self self pointer
        @param lblrtm If true, also execute lblrtm

        @return The calculated radiances
        '''
        if lblrtm:
            self.run_lblrtm()


        self.write_lbldis_parm()
        self.run_disort()
        radiance = self.write_data_to_file()
        return radiance
