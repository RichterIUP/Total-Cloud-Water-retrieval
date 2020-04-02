#!/usr/bin/python3
'''@package docstring
Write results to logfiles
'''

import sys
sys.path.append("./src/")
import aux2 as aux
import inp

'''
Number of digits in the formatstring
'''
DEC = 3


def log_prog_start():
    '''Initialise logfile
    '''
    with open("{}/retrieval_log.dat".format(inp.PATH), "a") as file_:
        file_.write("\n\n#########################################\n")
        file_.write("# L\\IWP\n")
        file_.write("#\n")
        file_.write("# Spec: {}\n".format(aux.FTIR))
        file_.write("# Started: {}\n".format(aux.TIME_INDEX))
        for i in inp.WINDOWS:
            file_.write("# Microwindow: {}\n".format(aux.MICROWINDOWS[i]))
        for element in aux.CLOUD_GRID:
            file_.write("# Cloud layer: {}\n".format(element))
        file_.write("# Mean temperature of cloud: {}\n".format(aux.CLOUD_TEMP))
        file_.write("#########################################\n\n")
    
    return

def log_pre_iter(variance_ra):
    '''Write PWV and noise to retrieval_log.dat
    
    @param variance_ra The calculated variance
    '''
    
    with open("{}/retrieval_log.dat".format(inp.PATH), "a") as file_:
        file_.write("#########################################\n")
        file_.write("# Noise: {:6.4f}\n".format(variance_ra))
        file_.write("# PWV (cm): {:5.3f}\n".format(aux.PRECIPITABLE_WATER_VAPOUR))
        file_.write("#########################################\n")


def write(text):
    '''Write arb. text to retrieval_log.dat
    
    @param text Text to be written to retrieval_log.dat
    '''
    with open("{}/retrieval_log.dat".format(inp.PATH), "a") as file_:
        file_.write("{}\n".format(text))

    return
