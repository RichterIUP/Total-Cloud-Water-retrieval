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


    wn_low = [200.0, 1500.0]#np.float64(cl_param[5])
    wn_high= [1500.0, 2600.0]#np.float64(cl_param[6])
    
    wavenumber_clear = [None for ii in range(2)]
    radiance_clear = [None for ii in range(2)]
    
    wavenumber_cloudy =  [None for ii in range(2)]
    radiance_cloudy = [None for ii in range(2)]
    
    inp.MCP = np.array([np.float64(cl_param[ii]) for ii in range(1, 5)])

    
    for jj in range(2):
        inp.WINDOWS = [0]
        aux.MICROWINDOWS = [[wn_low[jj], wn_high[jj]]]
        print(aux.MICROWINDOWS)
        
        '''
        Get clear sky and cloudy sky radiances
        '''
        [wavenumber_clear[jj], radiance_clear[jj]] = inversion.__set_up_retrieval()    
        [rms, wavenumber_cloudy[jj], radiance_cloudy[jj]] = inversion.retrieve()

    cre = pd.DataFrame({'wavenumber_clear' : np.concatenate((wavenumber_clear[0], wavenumber_clear[1])), \
                        'radiance_clear':    np.concatenate((radiance_clear[0], radiance_clear[1])), \
                        'wavenumber_cloudy': np.concatenate((wavenumber_cloudy[0], wavenumber_cloudy[1])), \
                        'radiance_cloudy':   np.concatenate((radiance_cloudy[0], wavenumber_cloudy[1]))})
    cre.to_csv("{}/results_cre.csv".format(inp.PATH), index=False)
    return
    
if __name__ == '__main__':

    cl_param = sys.argv[1:]# for ii in range(1, 8)]
    #print(cl_param)
    #exit(-1)
    main(cl_param)
