#!/usr/bin/python
'''@package docstring
Entrance point for L-IWP. Prepare the data and call the iteration
'''

import sys
import os
import datetime as dt

sys.path.append("./src")

import numpy as np

import inp
import aux2 as aux
import inversion
import read_input
import run_lbldis as rL
import get_atm


def main(cl_param):

    read_input.read_input(cl_param[0])
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

    inp.MCP_APRIORI = inp.MCP[:]

    aux.MICROWINDOWS.append([aux.MICROWINDOWS[inp.WINDOWS[0]][0], \
                aux.MICROWINDOWS[inp.WINDOWS[-1]][-1]])

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

    aux.SLOPE_RETR = False
    inversion.__set_up_retrieval()
    
    if not inp.FORWARD:
        '''
        Start the iteration using the chosen microwindows
        '''
        '''
        inp.FORWARD = True

        tt = [0.2, 1.0, 3.0, 4.0, 6.0]
        rt = [10, 15, 20, 25, 30, 35, 40, 45]
        fi = 0.5
        rad_lbldis = []
        rad_ftir = []
        tt_y = []
        for param_num in range(len(tt)):
            inp.MCP[0] = tt[param_num]*fi
            inp.MCP[1] = tt[param_num]*fi
            guess_apr = inversion.retrieve()
            rad_lbldis.append(guess_apr[0])
            rad_ftir.append(guess_apr[1])
            tt_y.append(tt[param_num])
        rad_ftir_av = np.mean(rad_ftir)
        tt_best = np.interp(rad_ftir_av, np.array(rad_lbldis), tt_y)

        inp.MCP[0] = tt_best * (1-fi)
        inp.MCP[1] = tt_best * fi        
            
        slope_lbldis = []
        slope_ftir = []
        rt_y = []
        fact = 5
        for param_num in range(len(rt)):
            inp.MCP[2] = rt[param_num] / ((fact-1)*fi+1)
            inp.MCP[3] = rt[param_num] * fact / ((fact-1)*fi+1)
            guess_apr = inversion.retrieve()
            slope_lbldis.append(guess_apr[2])
            slope_ftir.append(guess_apr[3])
            rt_y.append(rt[param_num])
        slope_ftir_av = np.mean(slope_ftir)
        rt_best = np.interp(slope_ftir_av, np.array(slope_lbldis), rt_y)
        inp.MCP[2] = rt_best / ((fact-1)*fi+1)
        inp.MCP[3] = rt_best * fact / ((fact-1)*fi+1)

        inp.FORWARD = False
        inp.MCP_APRIORI = inp.MCP[:]
        inversion.retrieve()
        '''
        
        inp.MCP = [ 0.056,  0.000, 14.990, 75.000]
        inp.FORWARD = False
        aux.SLOPE_RETR = True
        read_input.read_input(cl_param[0])
        inversion.__set_up_retrieval()
        inp.MCP_APRIORI = inp.MCP[:]
        inversion.retrieve()
        

    else:   
        '''
        Calculate the spectral radiance for the entire spectral range
        '''
        inp.WINDOWS = [-1]
        inversion.retrieve()
        
    return
    
if __name__ == '__main__':
    cl_param = []

    cl_param.append(sys.argv[1])#Name of the spectrum
    
    if len(sys.argv) > 2:
        cl_param.append(sys.argv[2])
    else:
        cl_param.append("TIR")

    main(cl_param)
