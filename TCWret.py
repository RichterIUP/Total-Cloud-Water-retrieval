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
chi = 1e20
chi_prev=1e30

NOW = dt.datetime.now()
directory = "{}_{}_{}_{}_{}_{}".format(NOW.month, NOW.day, \
      NOW.hour, NOW.minute, NOW.second, NOW.microsecond)
path = "/home/phi.richter/TCWret/OUTFOLDER"

resolution_only_od   = 2.0
resolution_retrieval = -0.5


'''
Perform a quick retrieval with low resolution to get a better estimation
of the correct result. Can also be done with the results of the
previous retrieval. Can be skipped by specifying inital values as
command line parameters or if an estimation from a previous retrieval
is available
'''

#if len(sys.argv) < 4:# or True:
#=======

#spectrum = "/home/phi.richter/Emission_Data/{}".format(spectrum)
spec = spectrum.split("/")[-1]
if os.path.exists("/home/phi.richter/TCWret/RESULTS/results_{}.nc".format(spec)):
    exit(-1)
if len(sys.argv) < 4:
    subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_only_od), str(tt), str(fi), str(rl), str(ri), "1", "0", "0", directory])
    with open("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory), "r") as f:
        cont = f.readlines()
        tt = float(cont[0])
        fi = float(cont[1])
        rl = float(cont[2])
        ri = float(cont[3])
        chi2 = float(cont[4])
    shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory))
    if tt > 8.0:
        exit(-1)


    while True:
        chi2_adj = [0.0 for ii in range(6)]
        tt_adj = [0.0 for ii in range(6)]
        fi_adj = [0.0 for ii in range(6)]
        rl_adj = [0.0 for ii in range(6)]
        ri_adj = [0.0 for ii in range(6)]
        '''
        Retrieve with low resolution
        '''
        subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_only_od), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
        with open("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory), "r") as f:
            cont = f.readlines()
            tt = float(cont[0])
            fi = float(cont[1])
            rl = float(cont[2])
            ri = float(cont[3])
        shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory))   
        best = [tt, fi, rl, ri, chi2]

        
        constant = fi * ri + (1-fi)*rl
        fi_ = [fi*0.9, fi, fi*1.1]
        if fi == 0.0:
            fi_ = [0.0, 0.1, 0.2]
        elif fi == 1.0:
            fi_ = [0.8, 0.9, 1.0]
        ri_ = [0.0 ,0.0 ,0.0]
        if fi_[0] != 0.0:
            ri_[0] = (constant - (1 - fi_[0]) * rl) / fi_[0]
        else:
            ri_[0] = rl
        ri_[1] = (constant - (1 - fi_[1]) * rl) / fi_[1]
        ri_[2]= (constant - (1 - fi_[2]) * rl) / fi_[2]
        directory_ = ["", "", ""]
        for ii in range(3):
            NOW = dt.datetime.now()
            directory_[ii] = "{}_{}_{}_{}_{}_{}".format(NOW.month, NOW.day, \
                                                        NOW.hour, NOW.minute, NOW.second, NOW.microsecond)
            subprocess.Popen(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_only_od), str(tt), str(fi_[ii]), str(rl), str(ri_[ii]), "1", "0", "0", directory_[ii]])
    
        for ii in range(3):
            counter = 0
            while not os.path.exists("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory_[ii])):
                time.sleep(5)
                counter += 1
                if counter == 120:
                    break
            try:
                with open("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory_[ii]), "r") as f:
                    cont = f.readlines()
                    tt_adj[ii] = float(cont[0])
                    fi_adj[ii] = float(cont[1])
                    rl_adj[ii] = float(cont[2])
                    ri_adj[ii] = float(cont[3])
                    chi2_adj[ii] = float(cont[4])
                shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory_[ii]))
            except FileNotFoundError:
                tt_adj[ii] = 10.0
                fi_adj[ii] = 1.0
                rl_adj[ii] = 10.0
                ri_adj[ii] = 10.0
                chi2_adj[ii] = 1e30

        rl_ = [0.0 ,0.0 ,0.0]
        rl_[0] = (constant - fi_[0] * ri) / (1 - fi_[0])
        rl_[1] = (constant - fi_[1] * ri) / (1 - fi_[1])
        if fi < 1.0:
            rl_[2] = (constant - fi_[2] * ri) / (1 - fi_[2])
        else:
            rl_[2] = ri
            fi_[2] = 1.0
        for ii in range(3):
            NOW = dt.datetime.now()
            directory_[ii] = "{}_{}_{}_{}_{}_{}".format(NOW.month, NOW.day, \
                                                        NOW.hour, NOW.minute, NOW.second, NOW.microsecond)
            subprocess.Popen(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_only_od), str(tt), str(fi_[ii]), str(rl_[ii]), str(ri), "1", "0", "0", directory_[ii]])
        
        for ii in range(3, 6):
            counter = 0
            while not os.path.exists("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory_[ii-3])):
                time.sleep(5)
                counter += 1
                if counter == 120:
                    break
            try:
                with open("{}/{}/{}/results.dat".format(path, spectrum.split("/")[-1], directory_[ii-3]), "r") as f:
                    cont = f.readlines()
                    tt_adj[ii] = float(cont[0])
                    fi_adj[ii] = float(cont[1])
                    rl_adj[ii] = float(cont[2])
                    ri_adj[ii] = float(cont[3])
                    chi2_adj[ii] = float(cont[4])
                shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory_[ii-3]))
            except FileNotFoundError:
                tt_adj[ii] = 10.0
                fi_adj[ii] = 1.0
                rl_adj[ii] = 10.0
                ri_adj[ii] = 10.0
                chi2_adj[ii] = 1e30

        idx = chi2_adj.index(min(chi2_adj))
        best = [tt_adj[idx], fi_adj[idx], rl_adj[idx], ri_adj[idx], chi2_adj[idx]]
        '''
        Break if less then 5% change in the cost function
        '''
        if 1-chi2_adj[idx]/best[-1] < 0.05:
            tt = best[0]
            fi = best[1]
            rl = best[2]
            ri = best[3]
            break


'''
Start the retrieval with high resolution
'''    
subprocess.call(["python3", "src/main.py", spectrum, windows, "20", "0", str(resolution_retrieval), str(tt), str(fi), str(rl), str(ri), "0", "0", "0", directory])
#os.system("mv {}/{} /home/phi.richter/OUTFOLDER_slayer/{}".format(path, spectrum.split("/")[-1], spectrum.split("/")[-1]))
shutil.rmtree("{}/{}/{}".format(path, spectrum.split("/")[-1], directory)) 
