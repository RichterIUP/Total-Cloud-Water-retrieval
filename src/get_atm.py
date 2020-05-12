#!/usr/bin/python
'''@package docstring
Read and manipulate the atmospheric profiles
'''

#Standard-library
import sys
sys.path.append("./src/")

#Third party modules
import numpy as np

#Self-defined modules
import aux2 as aux
import inp
import convert_humidity

'''
Indices to avoid magic numbers
'''
PRES = 0
ALT = 1
TEMP = 2
HUM = 3

def calc_pwv(prof_pres, prof_humd, prof_temp):
    '''Calculate the precipitable water path
    
    @param prof_pres The pressure profile in hPa
    @param prof_humd The humidity profile in the specified unit (see inp.py)
    @param prof_temp The temperature profile in K
    
    @return The precipitable water vapour
    '''
    delta_pres = np.abs(np.ediff1d(prof_pres))
    rho_wat = 1000.0#kg/m3
    gravity = 9.83#m/s2
    gas_constant_water = 461.5#J/(kg*K)
    precipitable_water_vapour = 0.0#m
    
    '''
    If the humidity is given as specific humidity, then it need to
    be converted to absolute humidity
    '''
    if inp.HUMIDITY == 'g/kg':
        prof_humd = convert_humidity.convert(["-s", np.array(prof_humd)], \
                                ["-tk", np.array(prof_temp)], \
                                ["-p", np.array(prof_pres)])[0] * 1000.0
       
    '''
    Calculate the PWV for each layer, sum up and convert to cm
    '''
    for loop in range(len(prof_pres)-1):
        wat_pres = prof_humd[loop] * 1e-3 * gas_constant_water * \
        prof_temp[loop]
        numerator = 0.622 * wat_pres
        denominator = prof_pres[loop] *1e2 - wat_pres
        precipitable_water_vapour = precipitable_water_vapour + \
            numerator / denominator * delta_pres[loop] *1e2 * 1 / \
            (rho_wat * gravity) * 100.0

    return precipitable_water_vapour
    
def interpolate_to_altitude_grid(atm):
    '''Interpolate the atmosphere data to a new altitude grid with the number of layers specifiec in aux.py
    
    @param atm The atmospheric data
    
    @return The interpolated profiles
    '''
    
    '''
    Convert to km
    '''
    cbh = aux.CLOUD_BASE

    cth = aux.CLOUD_TOP
    nlayer = aux.TOTAL_NUM_OF_LAYER
    
    atm_out = {'ALT': [], 'HMD': [], \
                'PRS': [], 'TMP': []}
    
    '''
    Set the maximum to 30km
    '''
    #for element in atm[1]:
    #    if element >= inp.MAX_ALT:
    #        alt_high = element
    #        break
    alt_grid = atm[1]
    alt_grid = alt_grid[np.where(alt_grid <= inp.MAX_ALT)[0]]
    alt_high = alt_grid[-1]

    '''
    The first layer is the ground layer
    '''
    alt_low = atm[1][0]
    if alt_low == 0.0:
        alt_low = 0.05
    print(alt_low, alt_high, nlayer)
    #alt_grid = atm[1]
    #alt_grid = np.geomspace(alt_low, alt_high, nlayer)
    
    '''
    Insert the cloud boundaries into the new altitude grid
    '''
    for jj in range(len(list(cbh))):
        for ii in range(len(alt_grid)-1):
            if alt_grid[ii] < cbh[jj]*1e-3 and alt_grid[ii+1] > cbh[jj]*1e-3:
                alt_grid[ii] = cbh[jj]*1e-3
            if alt_grid[ii] < cth[jj]*1e-3 and alt_grid[ii+1] > cth[jj]*1e-3:
                alt_grid[ii] = cth[jj]*1e-3

    '''
    Interpolate all profiles to the new grid
    '''
    atm_out['PRS'] = np.interp(np.array(alt_grid), np.array(atm[1]), np.array(atm[PRES]))
    atm_out['TMP'] = np.interp(np.array(alt_grid), np.array(atm[1]), np.array(atm[TEMP]))
    atm_out['HMD'] = np.interp(np.array(alt_grid), np.array(atm[1]), np.array(atm[HUM]))
    atm_out['ALT'] = alt_grid
    print(alt_grid)
    return atm_out
    
def get_atm():
    '''
    Read the atmospheric profile
    
    @return The precipitabe water vapour
    '''

    pres = aux.ATMOSPHERIC_GRID[PRES]
    alt = aux.ATMOSPHERIC_GRID[ALT]
    temp = aux.ATMOSPHERIC_GRID[TEMP]
    humd = aux.ATMOSPHERIC_GRID[HUM]

    '''
    Interpolate the profiles to the desired altitude grid
    '''
    print(len(aux.CLOUD_LAYERS))
    if len(aux.CLOUD_LAYERS) == 0:
        atm_interpolated = interpolate_to_altitude_grid([pres, alt, temp, humd])
        
        pres = atm_interpolated['PRS']
        alt = atm_interpolated['ALT']
        temp = atm_interpolated['TMP']
        humd = atm_interpolated['HMD']
     
    if inp.HUMIDITY == 'ppmv' or inp.HUMIDITY == 'g/m3':
        pwv = -1.0
    else:
        pwv = calc_pwv(pres, humd, temp)
        
    if len(aux.CLOUD_LAYERS) == 0:
        aux.CLOUD_BASE = np.round(aux.CLOUD_BASE, 3)
        aux.CLOUD_TOP = np.round(aux.CLOUD_TOP, 3)
        
    cloud_grid = []
    cloud_temp = []
    
    '''
    Find all layers within the cloud and calculate the mean temperature
    '''

    if len(aux.CLOUD_LAYERS) == 0: 
        for jj in range(len(aux.CLOUD_BASE)):
            for loop in range(len(alt)):
                if np.round(alt[loop]*1e3, 3) >= aux.CLOUD_BASE[jj] and \
                   np.round(alt[loop]*1e3, 3) <= aux.CLOUD_TOP[jj]:
                    cloud_grid.append(alt[loop]*1e3)
                    cloud_temp.append(temp[loop])
                    #if inp.NUMBER_OF_CLOUD_LAYER != 1:
        NUMBER_OF_CLOUD_LAYER = len(cloud_grid)
    else:
        for element in aux.CLOUD_LAYERS:
            for loop in range(len(alt)-1):
                if alt[loop]*1e3 <= element and alt[loop+1]*1e3 >= element:
                    if np.abs(alt[loop]*1e3 - element) <= np.abs(alt[loop+1]*1e3 - element):
                        cloud_grid.append(alt[loop]*1e3)
                        cloud_temp.append(temp[loop])
                        loop = loop + 1
                    else:
                        cloud_grid.append(alt[loop+1]*1e3)
                        cloud_temp.append(temp[loop+1])
        #cloud_grid = list(set(cloud_grid))
        idx_cloud_grid = aux.get_unique_indices(cloud_grid)
        cloud_grid_list = []
        for ii in idx_cloud_grid:
            cloud_grid_list.append(cloud_grid[ii])
        cloud_grid = cloud_grid_list
        
        idx_cloud_temp = aux.get_unique_indices(cloud_temp)
        cloud_temp_list = []
        for ii in idx_cloud_temp:
            cloud_temp_list.append(cloud_temp[ii])
        cloud_temp = cloud_temp_list
        
        #cloud_temp = list(set(cloud_temp))
        NUMBER_OF_CLOUD_LAYER = len(cloud_grid)
 
    aux.ATMOSPHERIC_GRID = [pres, alt, temp, humd]
    aux.CLOUD_GRID = np.array(cloud_grid[0:NUMBER_OF_CLOUD_LAYER])

    aux.CLOUD_TEMP = np.mean(cloud_temp[0:NUMBER_OF_CLOUD_LAYER])
    aux.PRECIPITABLE_WATER_VAPOUR = pwv
    return
