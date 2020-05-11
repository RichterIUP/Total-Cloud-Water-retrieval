#!/usr/bin/python
'''@package docstring
Entrance point for TCWret. Prepare the data and call the iteration
'''

import sys
import subprocess

if __name__ == '__main__':

    with open("../{}".format(sys.argv[1]), "r") as f:
        cont = f.readlines()

    for element in cont:
        radiances_files = "/home/phi.richter/radiances/PS.20170611_1413.nc"#element.split(" ")[0]
        atmospheric_prof = "/home/phi.richter/atm_prof/prof20170611_1413.nc"#element.split(" ")[1]
        cloud_file = "/home/phi.richter/cloud_files/CLOUDS.20170611.141300.nc"#element.split(" ")[2].rstrip()
        subprocess.call(["python3", "src/main.py", "/{}".format(radiances_files, atmospheric_prof, cloud_file), "TIR"])
