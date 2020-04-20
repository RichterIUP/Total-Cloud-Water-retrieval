#!/usr/bin/python
'''@package docstring
Entrance point for L-IWP. Prepare the data and call the iteration
'''

import sys
import subprocess

if __name__ == '__main__':

    with open(sys.argv[1], "r") as f:
        cont = f.readlines()
    #cont = ["/home/phi.richter/TESTCASES_nc/simulation_2012101705_cld1_0006_ch2_singleLayer.nc"]

    for element in cont:
        subprocess.call(["python3", "src/main.py", "{}".format(element.rstrip()), "TIR"])

