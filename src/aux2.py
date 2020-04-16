#!/usr/bin/python
'''@package docstring
Some auxillary functions and global variables. These global variables
are designed not to be changes by the user but by L-IWP
'''

import datetime as dt
import sys

#Third-party modules
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

sys.path.append("./src/")

import inp
import log

PRES = 0
ALT = 1
TEMP = 2
HUM = 3

'''
Read the current time. This should avoid problems if several spectra are
retrieved parallel
'''
NOW = dt.datetime.now()
TIME_INDEX = "{}_{}_{}_{}_{}_{}".format(NOW.month, NOW.day, \
              NOW.hour, NOW.minute, NOW.second, NOW.microsecond)
YEAR = 2017
MONTH = 0
DAY = 0
HOUR = 0
MINUTE = 0


'''
Wildcards for several paths
'''
FTIR = ''
LBLTP5 = ''
LBLTMP = ''
LBLDIR = ''
LBLLOG = ''

'''
Wildcards for trace gas profiles
'''
CO2_PROFILE = None
O3_PROFILE = None

#Parameter der Iteration
'''
Stepsize for the calculation of the derivative
'''
STEPSIZE_TAU = 1e-3#5e-4
STEPSIZE_RADIUS = 1e-2
INCREASE_LM = 4.0
DECREASE_LM = 2.0
#STEPSIZE_TT = 1e-3
#STEPSIZE_FI = 1e-3
#STEPSIZE_RL = 1e-3
#STEPSIZE_RI = 1e-3

'''
Maximum number of iteration
'''
MAX_ITER = 15

'''
Wildcards for atmospheric grid and cloud grid
'''
ATMOSPHERIC_GRID = [None, None, None, None]
CLOUD_GRID = [None]
CLOUD_TEMP = 0.0

'''
Wildcard for the solar zenith angle
'''
SOLAR_ZENITH_ANGLE = 0.0

'''
Wildcards for the cloud boundaries
'''
CLOUD_BASE = []
CLOUD_TOP = []
CLOUD_LAYERS = []
CLOUD_LAYERS_LIQ = []
CLOUD_LAYERS_ICE = []
'''
Total number of atmospheric layers. Must be less then 70
'''
TOTAL_NUM_OF_LAYER = 69

'''
Wildcard for the variance of the spectral radiances
'''
NOISE_FTIR = False

'''
Wildcards for the geoposition
''' 
LAT = 0.0
LON = 0.0

'''
Wildcards for the iteration matrices
'''
S_A_INV_MATRIX = None
S_Y_INV_MATRIX = None

'''
Wildcards for the radiances and wavenumber
'''
RADIANCE_LBLDIS = [[], [], [], [], [], [], [], [], []]
RADIANCE_SKIP = [[], [], [], [], [], [], [], [], []]
RADIANCE_FTIR = []
WAVENUMBER_FTIR = []

TOTAL_OPTICAL_DEPTH = []
ICE_FRACTION = []
RADIUS_LIQUID = []
RADIUS_ICE = []
CHI2 = []
RESIDUUM = []
T_MATRIX = []

VARIANCE_RA = 0.0
SLOPE_RETR = False

PROF_FILE = "prof20170611_1651.nc"
'''
Set this true if convergence reached
'''
CONVERGED = False

'''
Microwindows in cm-1
'''

MICROWINDOWS = [[] for ii in range(18)]
MICROWINDOWS[0] = [529.9, 533.0]
#MICROWINDOWS[1] = [538.0, 539.0]
#MICROWINDOWS[2] = [542.0, 544.0]
#MICROWINDOWS[3] = [551.0, 554.0]
MICROWINDOWS[1] = [558.5, 562.0]
MICROWINDOWS[2] = [572.0, 575.0]
#MICROWINDOWS[6] = [578.0, 579.0]
MICROWINDOWS[3] = [770.9, 774.8]
MICROWINDOWS[4] = [785.9, 790.7]
MICROWINDOWS[5] = [809.5, 813.5]
MICROWINDOWS[6] = [817.0, 823.5]
#MICROWINDOWS[8] = [820.5, 824.4]
MICROWINDOWS[7] = [828.3, 834.6]
#MICROWINDOWS[10] = [835.8, 838.7]
#MICROWINDOWS[8] = [843.1, 848.1]
#MICROWINDOWS[12] = [855.5, 858.0]
#MICROWINDOWS[8] = [860.1, 864.0]
#MICROWINDOWS[14] = [866.0, 870.0]
MICROWINDOWS[8] = [872.2, 877.5]
#MICROWINDOWS[16] = [884.5, 886.0]
MICROWINDOWS[9] = [891.9, 895.8]
MICROWINDOWS[10] = [898.2, 904.8]
#MICROWINDOWS[19] = [911.0, 913.0]
#MICROWINDOWS[20] = [914.5, 917.0]
MICROWINDOWS[11] = [929.6, 939.7]
MICROWINDOWS[12] = [958.0, 964.3]
#MICROWINDOWS[14] = [985.0, 991.5]
#MICROWINDOWS[22] = [931.5, 932.0]
#MICROWINDOWS[16] = [1076.6, 1084.8]
MICROWINDOWS[13] = [1092.2, 1098.1]
#MICROWINDOWS[24] = [1102.5, 1105.0]
MICROWINDOWS[14] = [1113.3, 1116.6]
MICROWINDOWS[15] = [1124.4, 1132.6]
MICROWINDOWS[16] = [1142.2, 1148.0]
MICROWINDOWS[17] = [1155.2, 1163.4]

TEMP_OF_CLOUD = 0.0
ENABLE_LM_DURING_ITER = True
TEST_STARTPARAM = []
CHI_ADJ = 10000000.0
CHI = 1e20
PRECIPITABLE_WATER_VAPOUR = 0.0


####################################################################################

def in_windows(wavenumber, num_windows):
    '''Determine if a wavenumber is part of a microwindow
    
    @param wavenumber Wavenumber to test
    @param num_windows Numbers of the microwindow
    @return True if inside of a microwindow, false if not
    '''
    for i in num_windows:
        if wavenumber > MICROWINDOWS[i][0] and wavenumber < MICROWINDOWS[i][1]:
            return True
    return False

####################################################################################

def change_resolution():
    '''Interpolate to a new spectral resolution and remove microwindows from
    the retrieval, if no datapoints are left
    '''

    global WAVENUMBER_FTIR
    global RADIANCE_FTIR
    wavenumber = []
    radiance = []

    wavenumber_raw = WAVENUMBER_FTIR
    radiance_raw = RADIANCE_FTIR
    num_of_lines = len(WAVENUMBER_FTIR)
    resolution = inp.RESOLUTION
    '''
    If RESOLUTION is negative, calculate the spectral resolution from the 
    spectrum
    '''
    if inp.RESOLUTION < 0.0:
        inp.RESOLUTION = np.mean(np.ediff1d(WAVENUMBER_FTIR))
    '''
    Remove microwindows with no datapoints
    '''
    #__redefine_microwindows()        
    
    '''
    Remove all datapoints outside the microwindows
    '''
    for loop in range(num_of_lines):
        '''
        Interpolate the spectral radiance to the chosen resolution. If
        the chosen resolution is negative, apply no interpolation
        '''
        if in_windows(wavenumber_raw[0]+loop*inp.RESOLUTION, inp.WINDOWS) and resolution > 0.0:
            wavenumber.append(wavenumber_raw[0]+loop*inp.RESOLUTION)
            radiance.append(np.interp(wavenumber[-1], wavenumber_raw, \
                                        radiance_raw))
        elif in_windows(wavenumber_raw[loop], inp.WINDOWS) and resolution < 0.0:
            wavenumber.append(wavenumber_raw[loop])
            radiance.append(radiance_raw[loop])

    WAVENUMBER_FTIR = wavenumber
    RADIANCE_FTIR = radiance

    return 

####################################################################################

def __redefine_microwindows():
    '''Remove microwindows, if no datapoints are inside the window
    '''
    global WAVENUMBER_FTIR
        
    diff = []
    for loop in inp.WINDOWS:
        diff.append([loop, MICROWINDOWS[loop][-1] - MICROWINDOWS[loop][-2]])

    inp.WINDOWS = []
    num_of_windows = len(diff)
    for loop in range(num_of_windows):
        if diff[loop][-1] > inp.RESOLUTION:
            inp.WINDOWS.append(diff[loop][-2])

    return
    
####################################################################################

def average(wavenumber, radiance):
    
    global SLOPE_RETR
    
    windows = inp.WINDOWS
    
    num_of_lines = len(wavenumber)
    '''
    Remove microwindows with no datapoints
    '''
    func = lambda x, a, b: a * x + b 
    new_radiance = []
    new_wavenumber = []
    '''
    Remove all datapoints outside the microwindows
    '''

    for win in inp.WINDOWS:
        dummy_wn = []
        dummy_ra = []
        for loop in range(num_of_lines):
            if in_windows(wavenumber[loop], [win]):
                dummy_wn.append(wavenumber[loop])
                dummy_ra.append(radiance[loop])
        if SLOPE_RETR:
            new_radiance.append((dummy_ra[0] - dummy_ra[-1])/(dummy_wn[0]-dummy_wn[-1]))
        else:
            new_radiance.append(np.median(dummy_ra))
        new_wavenumber.append(np.median(dummy_wn))
    
    '''
    Speichere nur den Slope und den RMS ab
    '''

    #if SLOPE_RETR: 
    #    new_radiance = [(new_radiance[0] - new_radiance[7])/(new_wavenumber[0] - new_wavenumber[7]), (new_radiance[8] - new_radiance[-1])/(new_wavenumber[8] - new_wavenumber[-1]), np.sum(new_radiance[0:8]), np.sum(new_radiance[8:-1])]
    #    new_wavenumber = range(len(new_radiance))#wavenumber[0:4]
    #else:
    #    new_radiance = radiance
    #    new_wavenumber = wavenumber

    return [np.array(new_wavenumber), np.array(new_radiance)]
        
    

####################################################################################

def calc_noise():
    '''
    Calculate the noise of the spectrum
    
    @return The calculated variance of the spectrum and the inverse of the covariance matrix S_y
    '''
    global NOISE_FTIR
    global WAVENUMBER_FTIR
    global RADIANCE_LBLDIS

    RADIANCE_LBLDIS = [[], [], [], [], [], [], [], [], []]

    
    if inp.STDDEV >= 0.0:
        '''
        Apply a non-weighted least squares algorithm
        '''
        stdDev = [inp.STDDEV for ii in range(len(WAVENUMBER_FTIR))]

    else:
        '''
        Use the standard deviation retrieved from the spectrum
        '''
        wavenumber = []
        spectral_radiance = []
        histogram = []
        wn_window = []
        ra_window = []
        stdDev = []
        dp = len(WAVENUMBER_FTIR)
        for line in range(dp):
            wavenumber.append(WAVENUMBER_FTIR[line])#line[0])
            spectral_radiance.append(RADIANCE_FTIR[line])#line[1])
        #   stdDev.append(0.16)
            for number in inp.WINDOWS:
                for window in [number]:
                    if in_windows(WAVENUMBER_FTIR[line], [window]):
                        wn_window.append(WAVENUMBER_FTIR[line])
                        ra_window.append(RADIANCE_FTIR[line])
                if len(ra_window) > 3:
                    func = lambda x, a, b, c: a * x**2 + b * x + c
                    popt, pcov = opt.curve_fit(func, wn_window, ra_window)
                    for ii in range(len(ra_window)):
                        histogram.append(np.array(ra_window[ii])-(popt[0]*np.array(wn_window[ii])**2+popt[1]*np.array(wn_window[ii])+popt[2]))
                    wn_window = []
                    ra_window = []
        for element in spectral_radiance:
            stdDev.append(np.std(histogram))


    return stdDev


####################################################################################

def s_y_inv():
    '''
    Calculate the noise of the spectrum
    
    @return The calculated variance of the spectrum and the inverse of the covariance matrix S_y
    '''
    global NOISE_FTIR
    global WAVENUMBER_FTIR
    global RADIANCE_LBLDIS

    RADIANCE_LBLDIS = [[], [], [], [], [], [], [], [], []]

    if inp.STDDEV < 0.0 and not inp.TESTCASE:
        variance_ra = np.mean(np.array(NOISE_FTIR))**2
    elif inp.STDDEV == 0.0:
        variance_ra = 1.0
    else:
        variance_ra = inp.STDDEV**2
    vec_error = np.array([variance_ra for ii in range(len(WAVENUMBER_FTIR))])
    s_y_inv_matrix = np.reciprocal(vec_error) * np.identity(len(vec_error))

    return [variance_ra, s_y_inv_matrix]


####################################################################################

def add_noise():
    '''Add noise to an artifical spectrum
    
    @return The noisy spectral radiance
    '''
    global RADIANCE_FTIR

    radiance = RADIANCE_FTIR
    return np.array(radiance) + np.random.normal(0.0, \
                   inp.STDDEV, len(radiance))

class MCPOutOfBoundsError(Exception):
    def __init__(self, value):
        self.value = value

class TooHighLMError(Exception):
    def __init__(self, value):
        self.value = value      

def get_unique_indices(array):
    idx = [0]
    for ii in range(1, len(array)):
        if array[ii-1] != array[ii]:
            idx.append(ii)
    return idx
    
if __name__ == '__main__':
    array = [150.0, 200.0, 250.0, 300.0, 300.0, 350.0, 380.0, 400.0, 2600.0, 2600.0, 2800.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3200.0, 3200.0, 3200.0, 3200.0, 3200.0, 3200.0, 3200.0, 3200.0, 3200.0, 3600.0, 3600.0, 3600.0, 3600.0, 3600.0, 3600.0, 3600.0, 3600.0, 3600.0, 3600.0, 3800.0, 3800.0, 3800.0, 4700.0, 4700.0, 4700.0, 4700.0, 4700.0, 4700.0, 4700.0, 4700.0, 4700.0, 4700.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5000.0, 5300.0, 5300.0, 5300.0, 5300.0, 5300.0, 5300.0, 5300.0, 5300.0, 5300.0, 5300.0, 5600.0, 5600.0, 5600.0, 5600.0, 5600.0, 5600.0, 5600.0, 5600.0, 5600.0, 5600.0, 5900.0, 5900.0, 5900.0, 5900.0, 5900.0, 5900.0, 5900.0, 5900.0, 5900.0, 6200.0, 6200.0, 6200.0, 6200.0, 6200.0, 6200.0, 6200.0, 6200.0, 6200.0, 6200.0, 6500.0, 6500.0, 6500.0, 6500.0, 6500.0, 6500.0, 6500.0, 6500.0, 6500.0, 6800.0, 6800.0, 7100.0, 7100.0, 7100.0, 7100.0, 7100.0, 7400.0, 7400.0, 7400.0, 7400.0, 7400.0, 7400.0, 7400.0, 7400.0, 7400.0, 7700.0, 7700.0, 7700.0, 7700.0, 7700.0, 7700.0, 7700.0, 7700.0, 7700.0, 7700.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8000.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0, 8400.0]
    idx = get_unique_indices(array)
    print(idx)
    array_unique = []
    for ii in idx:
        array_unique.append(array[ii])
    print(array_unique)
        
        
