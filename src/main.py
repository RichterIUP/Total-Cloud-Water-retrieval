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

import inp
import aux2 as aux
import inversion
import read_input
import run_lbldis as rL
import get_atm
import guess_apr


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
    print(inp.PATH)
    if not os.path.exists("{}".format(inp.PATH)):
        os.mkdir("{}".format(inp.PATH))
    print("{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1]))
    if not os.path.exists("{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1])):
        os.mkdir("{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1]))
    inp.PATH = "{}/{}/{}".format(inp.PATH, aux.FTIR.split("/")[-1], aux.TIME_INDEX)
    print(inp.PATH)
    if not os.path.exists("{}".format(inp.PATH)):
        os.mkdir("{}".format(inp.PATH))

    return 
    aux.SLOPE_RETR = False
    inversion.__set_up_retrieval()
    
    if not inp.FORWARD:
        '''
        Start the iteration using the chosen microwindows
        '''

        inp.FORWARD = True
        apr_list = []
        counter = 0 
        for fi in [0.05, 0.2, 0.4, 0.6, 0.8, 0.95]:
            
            apr_list.append(th.Thread(target=guess_apr.guess_apr, args=(fi, )))
            apr_list[-1].start()
            counter = counter + 1
            if (counter+1)%inp.NUM_OF_CPU == 0:
                for element in apr_list:
                    element.join()
                apr_list = []
                counter = 0
                
        for element in apr_list:
            element.join()
        inp.FORWARD = False


        idx = guess_apr.SEARCH_APR_RMS.indexpython3 TCWret.py batch_testcases/inputs_00
(min(guess_apr.SEARCH_APR_RMS))
        fi = guess_apr.SEARCH_APR_MCP[idx][0]
        tt = guess_apr.SEARCH_APR_MCP[idx][1]
        rl = guess_apr.SEARCH_APR_MCP[idx][2]
        ri = guess_apr.SEARCH_APR_MCP[idx][3]
        inp.MCP[0] = (1-fi)*tt
        inp.MCP[1] = fi*tt
        inp.MCP[2] = rl
        inp.MCP[3] = ri


        inp.FORWARD = False
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
    
    #main(["/home/phi.richter/TESTCASES_nc/simulation_2012021818_cld2_0001_ch2_singleLayer.nc", "TIR"])
    #exit(-1)
    with open(sys.argv[1], "r") as f:
        cont = f.readlines()
    for element in cont:
        cl_param = [element.rstrip(), "TIR"]
        #print(cl_param)
        main(cl_param)
