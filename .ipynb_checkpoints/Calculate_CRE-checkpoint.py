#!/usr/bin/python
'''@package docstring
Entrance point for TCWret. Prepare the data and call the iteration
'''

import subprocess

if __name__ == '__main__':

    spec = "/home/phi.richter/Emission_Data/PS.20170611_141300.nc"
    subprocess.call(["python3", "src/main_cre.py", "{}".format(spec), '1.0', '1.0', '10.', '30.', '150.0', '2800.'])
