#!/usr/bin/python
'''@package docstring
Entrance point for TCWret. Prepare the data and call the iteration
'''

import sys
import os
import subprocess
import shutil

if __name__ == '__main__':

    start = int(sys.argv[1]) * 322
    stop  = 322 + int(sys.argv[1]) * 322
    if stop > 5781:
        stop = 5781
        
    path = "/home/phi.richter/input_for_TCWret"
    files = sorted(os.listdir(path))[start:stop]
    print(files)
    for file_ in files:
        out_dir = "/mnt/beegfs/user/phi.richter/OUTFOLDER_Calc_CRE/"
        working_dir = "/home/phi.richter/TCWret"
        os.chdir(working_dir)
        spec = "{}/{}".format(path, file_)
        if not "TCWret_inp" in file_:
            continue
        subprocess.Popen(["python3", "src/main_cre.py", "{}".format(spec), "200.", "1500.", "1"])
        subprocess.Popen(["python3", "src/main_cre.py", "{}".format(spec), '1500.', '2800.', "2"])
        while not (os.path.exists("{}/{}/1/results_cre_200.0_1500.0.csv".format(out_dir, file_)) and \
                   os.path.exists("{}/{}/2/results_cre_1500.0_2800.0.csv".format(out_dir, file_))):
            continue
       