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
from scipy.interpolate import interp1d

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
fi = [0.2, 0.5, 0.8]
fi_best = 0.5
tt = np.array([0.2, 1.0, 3.0, 4.0])
rt = np.array([10, 20, 30, 40])

rt_y = []
tt_y = []

rad_lbldis = []
rad_ftir   = []
slope_lbldis = []
slope_ftir = []
'''
for ii in range(4):
    subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "1", str(resolution_only_od), str(tl[ii]), str(ti[ii]), str(rl[0]), str(ri[0]), "0", "0", "0", directory])
    if os.path.exists("{}/{}/{}/lbldis.spec".format(path, spectrum.split("/")[-1], directory)):
        with open("{}/{}/{}/lbldis.spec".format(path, spectrum.split("/")[-1], directory), "r") as f:
            cont = f.readlines()
            rad_lbldis.append(float(cont[-4]))
            rad_ftir.append(float(cont[-3]))
            tt_y.append(tl[ii]+ti[ii])
        shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory)) 

    
rad_ftir_av = np.mean(rad_ftir)
tt_best = np.interp(rad_ftir_av, np.array(rad_lbldis), tt)
tl_best = tt_best * (1-fi_best)
ti_best = tt_best * fi_best

for ii in range(4):
    rl_ii = rt[ii] / (2*fi_best+1)
    ri_ii = 3*rl_ii
    subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "1", str(resolution_only_od), str(tl_best), str(ti_best), str(rl_ii), str(ri_ii), "0", "0", "0", directory])
    if os.path.exists("{}/{}/{}/lbldis.spec".format(path, spectrum.split("/")[-1], directory)):
        with open("{}/{}/{}/lbldis.spec".format(path, spectrum.split("/")[-1], directory), "r") as f:
            cont = f.readlines()
            slope_lbldis.append(float(cont[-2]))
            slope_ftir.append(float(cont[-1]))
            rt_y.append(rt[ii])
        shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory)) 

slope_ftir_av = np.mean(slope_ftir)
rt_best = np.interp(slope_ftir_av, np.array(slope_lbldis), rt_y)
rl_best = rt_best / (2*fi_best+1)
ri_best = 3*rl_best

tt = tl_best
fi = ti_best#ti_best
rl = rl_best
ri = ri_best
print(tt, fi, ri, rl)
if tt < 0.0:
    tt = 1.0
if fi < 0.0:
    fi = 1.0
if rl < 0.0:
    rl = 10.0
if ri < 0.0:
    ri = 30.0
'''
tt = 0.182
fi = 0.383
rl = 3.520
ri= 57.423
subprocess.call(["python3", "src/main.py", spectrum, windows, "15", "0", str(resolution_od), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
#exit(-1)
'''
Start the retrieval with high resolution
'''    
#subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_retrieval), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
#os.system("mv {}/{} /home/phi.richter/OUTFOLDER_slayer/{}".format(path, spectrum.split("/")[-1], spectrum.split("/")[-1]))
#shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory)) 
