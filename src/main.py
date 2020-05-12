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

    read_input.read_input(cl_param[0], cl_param[1], cl_param[2])
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

    if inp.FORWARD:
        inp.WINDOWS = [0]
        aux.MICROWINDOWS = [[770., 1200.]]
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
    print(inp.FORWARD)
    if not inp.FORWARD:
        '''
        Start the iteration using the chosen microwindows
        '''

        inp.FORWARD = True
        apr_list = []
        counter = 0 
        fi = 0.5
        rl = inp.MCP[2]
        ri = inp.MCP[3]
        for tt in [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
            apr_list.append(th.Thread(target=guess_apr.guess_apr, args=(tt, fi, rl, ri, )))
            apr_list[-1].start()
            counter = counter + 1
            if (counter)%inp.NUM_OF_CPU == 0:
                for element in apr_list:
                    element.join()
                apr_list = []
                counter = 0
                
        for element in apr_list:
            element.join()
        inp.FORWARD = False


        idx = guess_apr.SEARCH_APR_RMS.index(min(guess_apr.SEARCH_APR_RMS))
        inp.MCP[0] = guess_apr.SEARCH_APR_MCP[idx][0]
        inp.MCP[1] = guess_apr.SEARCH_APR_MCP[idx][1]

        '''
        guess_apr.SEARCH_APR_MCP = []      
        guess_apr.SEARCH_APR_MCP = []
        inp.FORWARD = True
        apr_list = []
        counter = 0 
        tt = inp.MCP[0]+inp.MCP[1]
        for fi in [0.0, 0.25, 0.5, 0.75, 1.0]:
            for rl in [5.0, 10.0, 15.0]:
                for ri in [20.0, 30., 40.]:
                    apr_list.append(th.Thread(target=guess_apr.guess_apr, args=(tt, fi, rl, ri, )))
                    apr_list[-1].start()
                    counter = counter + 1
                    if (counter)%inp.NUM_OF_CPU == 0:
                        for element in apr_list:
                            element.join()
                        apr_list = []
                        counter = 0
                
        for element in apr_list:
            element.join()
            
        idx = guess_apr.SEARCH_APR_RMS.index(min(guess_apr.SEARCH_APR_RMS))
        inp.MCP[0] = guess_apr.SEARCH_APR_MCP[idx][0]
        inp.MCP[1] = guess_apr.SEARCH_APR_MCP[idx][1]
        inp.MCP[2] = guess_apr.SEARCH_APR_MCP[idx][2]
        inp.MCP[3] = guess_apr.SEARCH_APR_MCP[idx][3]
        '''
        
        inp.FORWARD = False
        inp.MCP_APRIORI = inp.MCP[:]
        inversion.retrieve()


    else:   
        '''
        Calculate the spectral radiance for the entire spectral range
        '''
        inversion.retrieve()
        
    return
    
if __name__ == '__main__':

    cl_param = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]]
    
    main(cl_param)
