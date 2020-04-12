#!/usr/bin/python3
'''@package docstring
Interface for calling L-IWP. 
'''

import subprocess
import shutil
import datetime as dt
import sys
import os
import time

spectrum = sys.argv[1]
windows = sys.argv[2]

try:
    tt = float(sys.argv[3])
    fi = float(sys.argv[4])
    rl = float(sys.argv[5])
    ri = float(sys.argv[6])
except IndexError:
    tt = -1
    fi = -1
    rl = -1
    ri = -1

'''
If the MCP is out of bounds, rerun the retrieval using different MCP
'''
NOW = dt.datetime.now()
directory = "{}_{}_{}_{}_{}_{}".format(NOW.month, NOW.day, \
      NOW.hour, NOW.minute, NOW.second, NOW.microsecond)
path = "/home/phi.richter/TCWret/OUTFOLDER"

resolution_only_od   = 2.0
resolution_retrieval = -0.5

spec = spectrum.split("/")[-1]

for ii in range(5):
    subprocess.call(["python3", "src/main.py", spectrum, windows, "10", "0", str(resolution_only_od), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
    with open("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory), "r") as f:
        cont = f.readlines()
        tt = float(cont[0])
        fi = float(cont[1])
        rl = float(cont[2])
        ri = float(cont[3])
        chi = float(cont[4]
    shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory))
        
    f = open("test", "a")
    f.write("[{}, {}, {}, {}]\t{}\n".format(0.1*tt, 0.1*fi, rl, ri, chi))
    f.close()
    
'''
Start the retrieval with high resolution
'''    
subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_retrieval), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
#os.system("mv {}/{} /home/phi.richter/OUTFOLDER_slayer/{}".format(path, spectrum.split("/")[-1], spectrum.split("/")[-1]))
#shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory)) 
