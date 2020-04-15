#!/usr/bin/python3
'''@package docstring
Read the case-file
'''

import os
import aux2 as aux
import inp
import numpy as np
import netCDF4 as nc
import convert_humidity

ALT_GRID = np.array([0.018,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.36,0.37,0.38,0.39,0.4,\
            0.45,0.5,0.52,0.57,0.62,0.72,0.78,0.81,0.86,0.92,1.1,1.25,1.4,1.6,1.8,2.0,\
            2.2,2.4,2.6,2.8,3.0,3.2,3.6,3.8,4.0,4.2,4.4,4.7,5.0,5.3,5.6,5.9,6.2,6.5,6.8,\
            7.1,7.4,7.7,8.0,8.4,8.8,9.2,9.6,10.0,10.5,11.0,11.5,12.0,12.5,13.0,14.0,15.0,\
            16.0,17.0,19.0,21.0,23.0,25.0,28.0,30.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,\
            40.0,41.0,42.0,43.0,44.0,45.0,60.0])*1000.0#m
RED_ALT = np.array([0.018,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.36,0.37,0.38,0.39,0.4,\
            0.45,0.5,0.52,0.57,0.62,0.72,0.78,0.81,0.86,0.92,1.1,1.25,1.4,1.6,1.8,2.0,\
            2.2,2.4,2.6,2.8,3.0,3.2,3.6,3.8,4.0,4.2,4.4,4.7,5.0,5.3,5.6,5.9,6.2,6.5,6.8,\
            7.1,7.4,7.7,8.0,8.4,8.8,9.2,9.6,10.0,10.5,11.0,11.5,12.0,12.5,13.0,14.0,15.0,\
            18.0,20.0,25.0,30.0])#km
CO2_PROFILE = np.array([405.3291,405.3511,405.38287,405.4892,405.5955,405.71167,405.83932,\
            405.96698,405.9925,406.01804,406.04355,406.0691,406.0946,406.17392,406.13852,406.12436,\
            406.08896,406.05356,405.98273,405.93842,405.91577,405.87802,405.8327,405.6968,405.629,\
            405.58966,405.53723,405.51263,405.51242,405.51224,405.51035,405.50354,405.49677,405.49,\
            405.476,405.42474,405.3991,405.37347,405.34763,405.32175,405.2829,405.24408,405.23068,\
            405.22003,405.2094,405.19553,405.16727,405.13904,405.1108,405.07285,404.98526,404.89764,\
            404.78082,404.51868,404.03073,403.5428,403.05484,401.9772,400.86832,399.75946,399.02756,\
            398.35522,397.6829,396.789,396.09695,395.41132,394.72803,393.70438,392.97867,392.41183,\
            391.8871,391.1372,390.78885,390.54504,390.42316,390.30127,390.24466,390.19,390.13535,\
            390.0807,390.02603,389.97137,389.93207,389.92896,389.92584,389.92276,389.91965,389.8731])#ppm
O3_PROFILE = np.array([0.028653354859483823,0.04816822763873067,0.062286941134498966, \
            0.06361594526875432,0.06562818790032524,0.06718196686573343,0.06934568150436156, \
            0.07221109493803343,0.07281664421210514,0.0734231754710233,0.0740355578454113,\
            0.0746479402197993,0.07526032259418729,0.07799265076504781,0.08023077911371435,\
            0.08096513007809727,0.08241287676725799,0.08334526155496731,0.08349884510978584,\
            0.08338812828554264,0.08327307725555631,0.08308132553891244,0.083021469655878,\
            0.08302719883536393,0.08272117649685426,0.08223139941051287,0.08107562464428096,\
            0.07999553181801221,0.07988466250603704,0.08051944724174995,0.08099382796904951,\
            0.08175633707444052,0.0836308223210225,0.08596569804296818,0.08831025869140584,\
            0.09328180570548257,0.09637165027119245,0.0971838523746458,0.09689815803950802,\
            0.09509955435961069,0.0897057292998673,0.0879370857432865,0.08773533817747944,\
            0.08744303993249522,0.08922775084966299,0.09428451110217577,0.09914140284099297,\
            0.09803966783899515,0.09257356576788885,0.09283326786816312,0.10299721665531336,\
            0.10986575614830985,0.10138482566077192,0.10646474178267536,0.11544223721793528,\
            0.107026374254647,0.10588336118258168,0.11060159906916746,0.13674266186538248,\
            0.4123691152699151,0.7865569122761463,1.0030563073798568,1.106550866520788,\
            1.331694421054253,1.6706791781524704,2.044586978803936,2.2670215619914376,\
            3.1097515340926822,4.310736265632324,5.515553476397224,5.391763629599047,\
            7.3319322799625315,8.748508668734008,9.302802233094853,9.38256177560111,\
            9.391204168001336,9.347020698354777,9.206611149299984,8.957953967903576,\
            8.61224204067852,8.168958481101429,7.627530729841363,7.010544484050797,\
            6.373037917960377,5.815065206180527,5.327697407238587,4.912090903394789,\
            2.370877704296867])#ppm

def read_cloud_boundaries(cbh_liq, cth_liq, cbh_ice, cth_ice):
    if len(cbh_liq) == 0:
        return cbh_ice, cth_ice
    elif len(cbh_ice) == 0:
        return cbh_liq, cth_liq
        
    cbh = []
    cth = []
    '''
    Test if clouds superimpose
    '''
    for ii in range(len(cbh_ice)):
        for ll in range(len(cbh_liq)):
            if cbh_ice[ii] == cbh_liq[ll]:
                cbh.append(cbh_ice[ii])
            elif cbh_ice[ii] < cbh_liq[ll] and cth_ice[ii] > cbh_liq[ll]:
                cbh.append(cbh_ice[ii])
            elif cbh_liq[ll] < cbh_ice[ii] and cth_liq[ll] > cbh_ice[ii]:
                cbh.append(cbh_liq[ll])
            if cth_ice[ii] == cth_liq[ll]:
                cth.append(cth_ice[ii])
            elif cth_ice[ii] > cth_liq[ll] and cbh_ice[ii] < cth_liq[ll]:
                cth.append(cth_ice[ii])
            elif cth_liq[ll] > cth_ice[ii] and cbh_liq[ll] < cth_ice[ii]:
                cth.append(cth_liq[ll])
    return cbh, cth
    
def read_testcase(fname):
    '''
    Read the data for the testcase
    
    @param fname The name of the netCDF4 file
    '''
    
    with nc.Dataset(fname) as dataset:
        aux.FTIR = fname.split("/")[-1]
        #month = aux.FTIR[7:9]
        #day = aux.FTIR[9:11]
        #year = aux.FTIR[3:7]
        #hour = aux.FTIR[12:14]
        #minute = aux.FTIR[14:16]
        aux.LAT = dataset.variables['lat'][0]
        aux.LON = dataset.variables['lon'][0]
        aux.SOLAR_ZENITH_ANGLE = dataset.variables['sza'][0]
        if aux.SOLAR_ZENITH_ANGLE > 90.0:
            aux.SOLAR_ZENITH_ANGLE = -1.0
        aux.WAVENUMBER_FTIR = np.array(dataset.variables['nu_cld'][:])
        aux.RADIANCE_FTIR = np.array(dataset.variables['rad_cld_down'][:])+inp.OFFSET
        #plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR)
        #plt.show()
        aux.NOISE_FTIR = np.array([0.0 for ii in range(len(aux.RADIANCE_FTIR))])
        cbh_liq = dataset.variables['cldbase_liq'][:]
        cbh_ice = dataset.variables['cldbase_ice'][:]
        cth_liq = dataset.variables['cldtop_liq'][:]
        cth_ice = dataset.variables['cldtop_ice'][:]
        
      


        if inp.CLOUD_BASE[0] == -1:
            if np.isnan(cbh_liq) and not np.isnan(cbh_ice): 
                aux.CLOUD_BASE = [np.float_(cbh_ice)*1000.0]
            else: 
                aux.CLOUD_BASE = [np.float_(cbh_liq)*1000.0]
        else:
            aux.CLOUD_BASE = inp.CLOUD_BASE
        if inp.CLOUD_TOP[0] == -1:
            if np.isnan(cth_liq) and not np.isnan(cth_ice): 
                aux.CLOUD_TOP = [np.float_(cth_ice)*1000.0]
            else: 
                aux.CLOUD_TOP = [np.float_(cth_liq)*1000.0]
        else:
            aux.CLOUD_TOP = inp.CLOUD_TOP

        if inp.TESTCASE and "multiLayer" not in fname:
            aux.CLOUD_BASE = aux.CLOUD_TOP 
            
        inp.HUMIDITY = 'g/m3'

        altitude = np.array(dataset.variables['prof_hts'][:])#m
        pressure = np.array(dataset.variables['prof_pres'][:])#hPa
        temperature = np.array(dataset.variables['prof_temp'][:])+inp.DISTURB_TEMPERATURE#K
        humidity = np.array(dataset.variables['prof_rh'][:])
        humidity = convert_humidity.convert(["-r", humidity], ["-tk", temperature], ["-p", pressure])[0]*1000.0
        aux.ATMOSPHERIC_GRID = [pressure, altitude, temperature, humidity+inp.DISTURB_HUMIDITY*humidity]
        aux.CO2_PROFILE = np.interp(altitude, ALT_GRID, CO2_PROFILE)
        aux.O3_PROFILE = np.interp(altitude, ALT_GRID, O3_PROFILE)
    return
    
def read_input(fname):
    '''Read the data from the netCDF4-file
    
    @param fname The name of the netCDF4 file
    '''
    if inp.TESTCASE:
        read_testcase(fname)
        return 
        
    dataset = nc.Dataset(fname)
    aux.FTIR = fname.split("/")[-1]
    aux.MONTH = aux.FTIR[7:9]
    aux.DAY = aux.FTIR[9:11]
    aux.YEAR = aux.FTIR[3:7]
    aux.HOUR = aux.FTIR[12:14]
    aux.MINUTE = aux.FTIR[14:16]
    aux.LAT = dataset.variables['lat'][0]
    aux.LON = dataset.variables['lon'][0]
    #f = open("/home/phi.richter/SAMPLES/samples_dship/{}{}{}_{}{}/case_{}{}{}_{}{}.0", "r")
    #cont = f.readlines()
    #f.close()
    #aux.SOLAR_ZENITH_ANGLE = float(90) - float(cont[4].split(":")[-1])
    aux.SOLAR_ZENITH_ANGLE = dataset.variables['sza'][0]
    aux.WAVENUMBER_FTIR = np.array(dataset.variables['wavenumber'][:])
    aux.RADIANCE_FTIR = np.array(dataset.variables['radiance'][:])+inp.OFFSET
    aux.NOISE_FTIR = aux.calc_noise()*2#np.array(dataset.variables['stdDev'][:])

    '''
    Test if there is an additional file with cloud height informations 
    from cloud radar
    '''
    path = ""
    for element in fname.split("/")[1:-1]:
        path = path + "/" + element
    print(path)
    for file_ in os.listdir(path):
        if inp.USE_CLOUD_FILES and ("CLOUD" in file_ and "{}{}{}.{}{}".format(aux.YEAR, aux.MONTH, aux.DAY, aux.HOUR, aux.MINUTE) in file_):
            aux.CLOUD_LAYERS = []
            cloudfile_= nc.Dataset("{}/{}".format(path, file_))
            levels = cloudfile_.variables['levels'][:]
            dz = np.mean(np.ediff1d(levels))

            lwc_lay = cloudfile_.variables['layers_liq'][:]
            iwc_lay = cloudfile_.variables['layers_ice'][:]
            lwp = np.sum(lwc_lay) * dz
            iwp = np.sum(iwc_lay) * dz

            if inp.FI_FROM_CLOUDS:
                inp.MCP[1] = iwp / (lwp+iwp)
                inp.MICROPHYSICAL_CLOUD_PARAMETERS = [inp.MCP[:], inp.MCP[:]]
                inp.MCP_APRIORI = np.array(inp.MCP[:])

            aux.CLOUD_LAYERS_ICE = levels[np.where(iwc_lay > 0.0)[0]]
            aux.CLOUD_LAYERS_LIQ = levels[np.where(lwc_lay > 0.0)[0]]
            for kk in range(len(lwc_lay)):
                if lwc_lay[kk] > 0.0 or iwc_lay[kk] > 0.0:
                    aux.CLOUD_LAYERS.append(levels[kk])
            aux.CLOUD_LAYERS = np.array(sorted(aux.CLOUD_LAYERS))
            if len(aux.CLOUD_LAYERS) != 0:
                inp.CLOUD_BASE = [-999]
                inp.CLOUD_TOP = [-999]
            break
    
    if inp.CLOUD_BASE[0] == -1:
        aux.CLOUD_BASE = [dataset.variables['cloud_base'][0]*1000.0]#m
    else:
        aux.CLOUD_BASE = inp.CLOUD_BASE
    if inp.CLOUD_TOP[0] == -1:
        aux.CLOUD_TOP = [dataset.variables['cloud_top'][0]*1000.0]#m
    else:
        aux.CLOUD_TOP = inp.CLOUD_TOP
    #if inp.TESTCASE:
    #    aux.CLOUD_BASE = aux.CLOUD_TOP

    altitude = np.array(dataset.variables['z'][:])
    if inp.HUMIDITY == "g/kg":
        key_humd = 'sh'
    elif inp.HUMIDITY == "g/m3":
        key_humd = 'ah'
    elif inp.HUMIDITY == 'ppmv':
        key_humd = 'h2o'
    elif inp.HUMIDITY == "%":
        key_humd = "rh"
    pressure_pre_int = np.array(dataset.variables['P'][:])
    temperature_pre_int = np.array(dataset.variables['T'][:])
    humidity_pre_int = np.array(dataset.variables[key_humd][:])
    if inp.HUMIDITY == "%":
        humidity_pre_int = np.array(convert_humidity.convert(["-r", humidity_pre_int], \
                                                             ["-tk", temperature_pre_int], \
                                                             ["-p", pressure_pre_int])[1]*1000.0)
        inp.HUMIDITY = "g/kg"
    if len(altitude) > 70:
        pressure = np.interp(RED_ALT, altitude, pressure_pre_int)
        temperature = np.interp(RED_ALT, altitude, temperature_pre_int+inp.DISTURB_TEMPERATURE)
        humidity = np.interp(RED_ALT, altitude, humidity_pre_int)
        altitude = RED_ALT

    aux.ATMOSPHERIC_GRID = [pressure, altitude, temperature, humidity+inp.DISTURB_HUMIDITY*humidity]
    aux.CO2_PROFILE = np.array(dataset.variables['co2'][:])
    aux.O3_PROFILE = np.array(dataset.variables['o3'][:])
    dataset.close()
    return
    
if __name__ == '__main__':
    inp.TESTCASE = True
    read_input("/home/philippr/TESTCASES_nc/simulation_2012021817_cld1_0006_ch2_singleLayer.nc")
    #read_input("/home/philippr/Arbeit/RETRIEVAL/PS.20170611_161500.nc")
