#!/usr/bin/python
'''@package docstring
Entrance point for TCWret. Prepare the data and call the iteration
'''

import sys
import subprocess

if __name__ == '__main__':

    #with open("../{}".format(sys.argv[1]), "r") as f:
    #    cont = f.readlines()

    #for element in cont:
    radiances_files = "/home/phi.richter/radiances/nyem.20200506_184033.nc"#element.split(" ")[0]
    atmospheric_prof = "/home/phi.richter/atm_prof/prof20200506_184033.nc"#element.split(" ")[1]
    cloud_file = "/home/phi.richter/cloud_files/CLOUDS.20200506.184033.nc"#element.split(" ")[2].rstrip()
    subprocess.call(["python3", "src/main.py", "{}".format(radiances_files), "{}".format(atmospheric_prof), "{}".format(cloud_file), "TIR"])
    
