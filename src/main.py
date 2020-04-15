#!/usr/bin/python
'''@package docstring
Entrance point for L-IWP. Prepare the data and call the iteration
'''

import sys
import os

sys.path.append("./src")

import numpy as np

import inp
import aux2 as aux
import inversion
import read_input
import run_lbldis as rL


def main(cl_param, ONLY_OD=False, SEARCH_INIT=False, DIR="", ADJUST_RADII=False):
    '''
    Write the command line parameters to the corresponding variables
    @param cl_param All the command line parameters: 
    [PATH_TO_SPECTRUM (str), WINDOW_RANGE (str), MAX_ITER (int), FORWARD (bool), 
    RESOLUTION (float), MCP (list of float, len=4)
    @return The retrieved MCP: tt, fi, rl, ri
    '''

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
    aux.MAX_ITER = cl_param[2]
    inp.FORWARD = cl_param[3]
    inp.RESOLUTION = cl_param[4]
    if inp.FI_FROM_CLOUDS:
        fi = inp.MCP[1]
        inp.MCP = cl_param[5][:]
        inp.MCP[1] = fi
    else:
        inp.MCP = cl_param[5][:]

    inp.SEARCH_INIT = bool(int(SEARCH_INIT))
    if int(ONLY_OD) > 0:
        inp.ONLY_OD = True
        inp.NO_OD = False
    elif int(ONLY_OD) == 0:
        inp.ONLY_OD = False
        inp.NO_OD = False
    else:
        inp.ONLY_OD = False
        inp.NO_OD = False#True
        inp.FIXED_FI = True
    
    inp.ONLY_OD = bool(int(ONLY_OD))
    inp.ADJUST_RADII = bool(int(ADJUST_RADII))

    aux.TIME_INDEX = DIR
    
    if inp.SEARCH_INIT or inp.ONLY_OD:
        inp.LM_INIT = 0.0

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


    aux.LBLTP5 = "{}/tp5_{}".format(inp.PATH, aux.TIME_INDEX)
    aux.LBLTMP = '{}'.format(inp.PATH)
    aux.LBLLOG = '{}/lbllog.txt'.format(inp.PATH)
    aux.LBLDIR = "{}/lblout_{}".format(inp.PATH, aux.FTIR.split("/")[-1])

    '''
    Create the directory for the optical depths of LBLRTM
    '''
    if not os.path.exists("{}".format(aux.LBLDIR)):
        os.mkdir("{}".format(aux.LBLDIR))
    rL.forward_run([0.1, 0.1, 7.0, 7.0], [0, 1.0], True, 0)

    if not inp.FORWARD:
        '''
        Start the iteration using the chosen microwindows
        '''
    
        inp.FORWARD = True
        tt = np.array([0.2, 1.0, 3.0, 4.0, 6.0])
        rt = np.array([10, 15, 20, 25, 30, 35, 40, 45])
        fi = 0.5
        rad_lbldis = []
        rad_ftir = []
        tt_y = []
        for param_num in range(len(tt)):
            inp.MCP[0] = tt[param_num]/2.0
            inp.MCP[1] = tt[param_num]/2.0
            guess_apr = inversion.retrieve()
            rad_lbldis.append(guess_apr[0])
            rad_ftir.append(guess_apr[1])
            tt_y.append(tt[param_num])
        rad_ftir_av = np.mean(rad_ftir)
        tt_best = np.interp(rad_ftir_av, np.array(rad_lbldis), tt_y)
        inp.MCP[0] = tt_best * (1-fi)
        inp.MCP[1] = tt_best * fi
        f = open("MCP", "a")
        f.write("{}\n".format(inp.MCP))
        f.close()
        slope_lbldis = []
        slope_ftir = []
        rt_y
        for param_num in range(len(rt)):
            inp.MCP[2] = rt[param_num] / (2*fi+1)
            inp.MCP[3] = 3*inp.MCP[2]
            guess_apr = inversion.retrieve()
            slope_lbldis.append(guess_apr[2])
            slope_ftir.append(guess_apr[3])
            rt_y.append(tt[param_num])
        slope_ftir_av = np.mean(slope_ftir)
        rt_best = np.interp(slope_ftir_av, np.array(slope_lbldis), rt_y)
        inp.MCP[2] = rt_best / (2*fi+1)
        inp.MCP[3] = 3*inp.MCP[2]
        f = open("MCP", "a")
        f.write("{}\n".format(inp.MCP))
        f.close()
        exit(-1)
        
        inp.FORWARD = False
        inversion.retrieve()
        #with open("{}/results.dat".format(inp.PATH), "w") as f:
        #    f.write("{}\n".format(aux.TOTAL_OPTICAL_DEPTH[-1]))
        #    f.write("{}\n".format(aux.ICE_FRACTION[-1]))
        #    f.write("{}\n".format(aux.RADIUS_LIQUID[-1]))
        #    f.write("{}\n".format(aux.RADIUS_ICE[-1]))
        #    f.write("{}\n".format(aux.CHI_ADJ))
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
        
    if len(sys.argv) > 3:
        cl_param.append(int(sys.argv[3]))
    else:
        cl_param.append(aux.MAX_ITER)
        
    if len(sys.argv) > 4:
        cl_param.append(bool(int(sys.argv[4])))
    else:
        cl_param.append(inp.FORWARD)
        
    if len(sys.argv) > 5:
        if float(sys.argv[5]) < 0.0:
            cl_param.append(inp.RESOLUTION)
        else:
            cl_param.append(float(sys.argv[5]))
    else:
        cl_param.append(inp.RESOLUTION)

    if len(sys.argv) > 6:
        if float(sys.argv[6]) < 0.0:
            cl_param.append(inp.MCP)
        else:
            cl_param.append([float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9])])
    else:
        cl_param.append(inp.MCP)

    main(cl_param, ONLY_OD=sys.argv[10], SEARCH_INIT=sys.argv[11], DIR=sys.argv[13], ADJUST_RADII=sys.argv[12])
