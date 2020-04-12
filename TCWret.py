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

import numpy as np

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

tl = [0.1, 0.5, 1.5, 2.0]
ti = [0.1, 0.5, 1.5, 2.0]
rl = [6.67, 13.33, 20, 26.67]
ri = [13.33, 26.67, 40, 53.33]

tt = np.array([0.2, 1.0, 3.0, 4.0])
rt = np.array([10, 20, 30, 40])

rad_lbldis = [0, 0, 0, 0]
rad_ftir   = [0, 0, 0, 0]
slope_lbldis = [0, 0, 0, 0]
slope_ftir = [0, 0, 0, 0]

for ii in range(4):
    subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "1", str(resolution_only_od), str(tl[ii]), str(ti[ii]), str(rl[0]), str(ri[0]), "0", "0", "0", directory])
    with open("{}/{}/{}/lbldis.spec".format(path, spectrum.split("/")[-1], directory), "r") as f:
        cont = f.readlines()
        rad_lbldis[ii] = float(cont[-4])
        rad_ftir[ii] =   float(cont[-3])
    
rad_ftir_av = np.mean(rad_ftir)
tt_best = np.interp(rad_ftir_av, np.array(rad_lbldis), tt)
tl_best = tt_best / 2.0
ti_best = tl_best

for ii in range(4):
    rl_ii = rt[ii] / 1.5
    ri_ii = rt[ii] / 0.75
    subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "1", str(resolution_only_od), str(tl_best), str(ti_best), str(rl_ii), str(ri_ii), "0", "0", "0", directory])
    with open("{}/{}/{}/lbldis.spec".format(path, spectrum.split("/")[-1], directory), "r") as f:
        cont = f.readlines()
        slope_lbldis[ii] = float(cont[-2])
        slope_ftir[ii] =   float(cont[-1])
slope_ftir_av = np.mean(slope_ftir)
rt_best = np.interp(slope_ftir_av, np.array(slope_lbldis), rt)
rl_best = rt_best / 1.5
ri_best = rt_best / 0.75

tt = tl_best
fi = ti_best
rl = rl_best
ri = ri_best
print(tt, fi, ri, rl)
exit(-1)
'''
Start the retrieval with high resolution
'''    
#subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_retrieval), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
#os.system("mv {}/{} /home/phi.richter/OUTFOLDER_slayer/{}".format(path, spectrum.split("/")[-1], spectrum.split("/")[-1]))
#shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory)) 
