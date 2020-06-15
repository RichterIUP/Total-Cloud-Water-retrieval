#!/usr/bin/python
'''@package docstring
Entrance point for TCWret. Prepare the data and call the iteration
'''
import os
import datetime as dt
import time
import subprocess
import shutil

if __name__ == '__main__':

    idx = 15
    path = "/home/phi.richter/input_for_TCWret"
    files = sorted(os.listdir(path))
    spec = "{}/{}".format(path, files[idx])
    subprocess.call(["python3", "src/main_cre.py", "{}".format(spec), "200.", "1500.", "1"])
    subprocess.call(["python3", "src/main_cre.py", "{}".format(spec), '1500.', '2800.', "2"])
    shutil.copyfile("/mnt/beegfs/user/phi.richter/OUTFOLDER_Calc_CRE/{}/1/results_cre_200.0_1500.0.csv".format(files[idx]), "/home/phi.richter/results_cre_200.0_1500.0.csv")    
    shutil.copyfile("/mnt/beegfs/user/phi.richter/OUTFOLDER_Calc_CRE/{}/2/results_cre_1500.0_2800.0.csv".format(files[idx]), "/home/phi.richter/results_cre_1500.0_2800.0.csv")