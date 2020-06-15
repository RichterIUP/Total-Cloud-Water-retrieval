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

    read_input.read_input_cre(cl_param[0])
    #print(aux.CLOUD_LAYERS)
    #print(aux.ATMOSPHERIC_GRID[1])
    #exit(-1)
    aux.TIME_INDEX = cl_param[-1]
    inp.RESOLUTION = 1.0
    
    #with open()

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


    wn_low = np.float64(cl_param[-3])
    wn_high= np.float64(cl_param[-2])
    
    #inp.MCP = np.array([np.float64(cl_param[ii]) for ii in range(1, 5)])

    
    inp.WINDOWS = [0]
    aux.MICROWINDOWS = [[wn_low, wn_high]]
    print(aux.MICROWINDOWS)
        
    '''
    Get clear sky and cloudy sky radiances
    '''
    [wavenumber_clear, radiance_clear] = inversion.__set_up_retrieval()    
    [rms, wavenumber_cloudy, radiance_cloudy] = inversion.retrieve()

    cre = pd.DataFrame({'wavenumber_clear' : wavenumber_clear, \
                        'radiance_clear':    radiance_clear, \
                        'wavenumber_cloudy': wavenumber_cloudy, \
                        'radiance_cloudy':   radiance_cloudy})
    outfile = "{}/results_cre_{}_{}.csv".format(inp.PATH, wn_low, wn_high)
    cre.to_csv(outfile, index=False)
    print(outfile)
    return
    
if __name__ == '__main__':

    cl_param = sys.argv[1:]# for ii in range(1, 8)]
    #print(cl_param)
    #exit(-1)
    main(cl_param)
