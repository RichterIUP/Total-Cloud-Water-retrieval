#!/usr/bin/python3

import os

def create_slurm_batches(path, fname, spec_per_files=10):
    with open("{}/{}".format(path, fname), "r") as file_all_spec:
        spec_names = file_all_spec.readlines()
    os.system("rm {}/{}_*".format(path, fname))
    os.system("rm {}/master".format(path))
    counter = -1
    for ii in range(len(spec_names)):
        if ii%spec_per_files == 0:
            counter += 1
            with open("{}/{}_{}".format(path, fname, counter), "w") as outfile:
                outfile.write("#!/bin/bash\n")
                outfile.write("#SBATCH --job-name={}\n".format(fname))
                outfile.write("#SBATCH -n 9\n")
                outfile.write("#SBATCH -N 1\n")
                outfile.write("#SBATCH -t 4320\n")
                outfile.write("#SBATCH -p all\n")
                outfile.write("#SBATCH --mem=9000\n")
                outfile.write("#SBATCH --open-mode=append\n")
                outfile.write("#SBATCH --mail-type=END\n")
                outfile.write("#SBATCH --mail-user=phi.richter@iup.physik.uni-bremen.de\n")
                outfile.write("module unuse /home/eb/modules/all\n")
                outfile.write("module unuse /home/eb/modules/LAMOS\n")
                outfile.write("module use /home/eb/modules/_legacy/all\n")
                outfile.write("module load Python\n")
                outfile.write("module load matplotlib\n")
                outfile.write("module load netcdf4-python\n")

        with open("{}/{}_{}".format(path, fname, counter), "a") as outfile:
            #outfile.write("python3 L-IWP.py '/home/phi.richter/Emission_Data/{}' TIR\n".format(spec_names[ii].rstrip()))
            outfile.write("{}\n".format(spec_names[ii].rstrip()))

    with open("{}/master".format(path), "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH --job-name=master_{}\n".format(fname))
        outfile.write("#SBATCH -n 9\n")
        outfile.write("#SBATCH -N 1\n")
        outfile.write("#SBATCH -t 4320\n")
        outfile.write("#SBATCH -p all\n")
        outfile.write("#SBATCH --mem=9000\n")
        outfile.write("#SBATCH --open-mode=append\n")
        outfile.write("#SBATCH --mail-type=END\n")
        outfile.write("#SBATCH --mail-user=phi.richter@iup.physik.uni-bremen.de\n")
        outfile.write("module unuse /home/eb/modules/all\n")
        outfile.write("module unuse /home/eb/modules/LAMOS\n")
        outfile.write("module use /home/eb/modules/_legacy/all\n")
        outfile.write("module load Python\n")
        outfile.write("module load matplotlib\n")
        outfile.write("module load netcdf4-python\n")
        for ii in range(counter+1):
            outfile.write("sbatch {}/{}_{}\n".format(path, fname, ii))
            outfile.write("sleep 20\n")
        
if __name__ == '__main__':
    x = 15
    create_slurm_batches("batch_singlelayer", "LIWP_singlelayer", x)
