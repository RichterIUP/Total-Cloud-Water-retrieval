#!/usr/bin/python
'''@package docstring
Entrance point for L-IWP. Prepare the data and call the iteration
'''

import sys
import os
import datetime as dt
import threading as th

sys.path.append("./src")

import numpy as np
import pandas as pd

import inp
import aux2 as aux
import inversion
import read_input
import run_lbldis as rL
import get_atm
import guess_apr


def main(cl_param):

    read_input.read_input(cl_param[0], calc_cre=True)#, cl_param[1], cl_param[2])
    if cl_param[1] == "TIR":
        inp.WINDOWS = inp.TIR
    elif cl_param[1] == "FIR":
        inp.WINDOWS = inp.FIR_TIR
    elif cl_param[1] == "FIRo":
        inp.WINDOWS = inp.FIR
    elif cl_param[1] == "FIRs":
        inp.WINDOWS = inp.FIR_MCP
    elif cl_param[1] == "NIR":
        inp.WINDOWS = inp.NIR
    NOW = dt.datetime.now()
    directory = "{}_{}_{}_{}_{}_{}".format(NOW.month, NOW.day, \
      NOW.hour, NOW.minute, NOW.second, NOW.microsecond)
    aux.TIME_INDEX = directory

    '''
    Create all the necessary folders
    '''
    if not os.path.exists("{}".format(inp.PATH)):
        os.mkdir("{}".format(inp.PATH))
    if not os.path.exists("{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1])):
        os.mkdir("{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1]))
    inp.PATH = "{}/{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1], aux.TIME_INDEX)
    if not os.path.exists("{}".format(inp.PATH)):
        os.mkdir("{}".format(inp.PATH))
    inp.FORWARD = True


    wn_low = np.float64(cl_param[6])
    wn_high= np.float64(cl_param[7])
    
    wavenumber_clear = np.array([])
    radiance_clear = np.array([])
    
    wavenumber_cloudy = np.array([])
    radiance_cloudy = np.array([])
    
    while True:
        if wn_high - wn_low > 2000.0:
            wn_high = wn_low + 1500.0
        inp.WINDOWS = [0]
        aux.MICROWINDOWS = [[wn_low, wn_high]]
        
        '''
        Get clear sky radiances
        '''
        inversion.__set_up_retrieval()
        inp.MCP = [0., 0., 10., 30.]
        [rms, wavenumber, radiance] = inversion.retrieve()
        
        wavenumber_clear = np.concatenate((wavenumber_clear, wavenumber))
        radiance_clear = np.concatenate((radiance_clear, radiance))
        '''
        Get cloudy sky radiances
        '''
        inp.MCP = np.array([np.float64(cl_param[ii]) for ii in range(1, 5)])
        
        [rms, wavenumber, radiance] = inversion.retrieve()
        
        wavenumber_cloudy = np.concatenate((wavenumber_clear, wavenumber))
        radiance_cloudy = np.concatenate((radiance_clear, radiance))
        
        if wn_high != np.float64(cl_param[7]):
            wn_low = wn_high
            wn_high = np.float64(cl_param[7])
        else:
            break
    cre = pd.DataFrame({'wavenumber_clear' : wavenumber_clear, \
                        'radiance_clear': radiance_clear, \
                        'wavenumber_cloudy': wavenumber_cloudy, \
                        'radiance_cloudy': radiance_cloudy})
    cre.to_csv("{}/results_cre.csv".format(inp.PATH), index=False)
    return
    
if __name__ == '__main__':

    cl_param = [sys.argv[ii] for ii in range(1, 8)]
    
    main(cl_param)
