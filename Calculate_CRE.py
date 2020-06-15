#!/usr/bin/python
'''@package docstring
Entrance point for TCWret. Prepare the data and call the iteration
'''

import sys
import os
import datetime as dt
import time
import subprocess
import shutil

if __name__ == '__main__':

    for idx in range(2,4,1):#int(sys.argv[1]), int(sys.argv[2]), 1):
        path = "/home/phi.richter/input_for_TCWret"
        files = sorted(os.listdir(path))
        spec = "{}/{}".format(path, files[idx])
        subprocess.call(["python3", "src/main_cre.py", "{}".format(spec), "200.", "1500.", "1"])
        subprocess.call(["python3", "src/main_cre.py", "{}".format(spec), '1500.', '2800.', "2"])

        path_out = "/home/phi.richter/LW_Downward_radiation"
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        os.chdir(path_out)
        if not os.path.exists(files[idx])):
            os.mkdir(files[idx])
        else:
            shutil.rmtree(files[idx])
            os.mkdir(files[idx])
        shutil.copyfile(src="/mnt/beegfs/user/phi.richter/OUTFOLDER_Calc_CRE/{}/1/results_cre_200.0_1500.0.csv".format(files[idx]), \
                        dst="/home/phi.richter/LW_Downward_radiation/{}/results_cre_200.0_1500.0.csv".format(files[idx)    
        shutil.copyfile(src="/mnt/beegfs/user/phi.richter/OUTFOLDER_Calc_CRE/{}/2/results_cre_1500.0_2800.0.csv".format(files[idx]), \
                        dst="/home/phi.richter/LW_Downward_radiation/{}/results_cre_1500.0_2800.0.csv".format(files[idx])
