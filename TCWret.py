#!/usr/bin/python
'''@package docstring
Entrance point for L-IWP. Prepare the data and call the iteration
'''

import sys
import subprocess

if __name__ == '__main__':

    with open(sys.argv[1], "r") as f:
        cont = f.readlines()

    for element in cont:
        subprocess.call(["python3", "src/main.py", "/{}".format(element.rstrip()), "FIR"])
    
    #spec = "/home/phi.richter/Emission_Data/PS.20170616_211900.nc"
    #tl = '0.5'
    #ti = '0.5'
    #for rl in ['5', '10', '15', '20', '25', '30', '35', '40']:
    #    for ri in ['10', '20', '30', '40', '50', '60', '70']:
    #        with open("log", "a") as f:
    #            f.write("{} {}\n".format(rl, ri))
    #            subprocess.call(["python3", "src/main.py", "/{}".format(spec), "TIR", tl, ti, rl, ri])
