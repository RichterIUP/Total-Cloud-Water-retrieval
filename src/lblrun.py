#!/usr/bin/python

#
# Originally written by    Paul van Delst, CIMSS/SSEC, 20-Jun-1996
# 			   paul.vandelst@ssec.wisc.edu
#
# Translation to Python by Philipp Richter, IUP Bremen, 11-Jul-2019
#                          phi.richter@iup.physik.uni-bremen.de

import os
import sys
import datetime as dt
import shutil

#move_file from lblrun
def move_file(source, destination, lbllog):
    f = open("{}".format(lbllog), "a")
    #if(os.path.exists(source)):
    try:
        shutil.copy2("{}".format(source), "{}".format(destination))
    except FileNotFoundError:
        f.write("     Error copying '{}' to '{}'. Attempting link instead\n".format(source, destination))
        os.symlink("{}".format(source), "{}".format(destination))
        if(not os.path.isfile("{}".format(destination))):
            f.write("     Error linking file also! Sorry. :(\n")
    else:
        os.remove("{}".format(source))
    return
                         
def lblrun(lbltp5, lbldir, lbllog, lbl_home):
    
    LBL_HOME = lbl_home
    f = open("{}".format(lbllog), "w")
    # -------------
    # RCS ID string
    # -------------
    f.write("\n")
    f.write("Adapted from: lblrun,v 1.12 1999/04/15 13:53:45 paulv Exp $\n")
    f.write("\n")
    
    # ---------------------------------------------------
    # Set name of TAPE3 file,  use default
    # ---------------------------------------------------
    
    T3_FILE = "tape3.data"

    # ------------------------------------------------------
    # Create the data and time flagged LBLRTM work directory
    # ------------------------------------------------------
    
    # -- Get the start date and time

    now = dt.datetime.now()
    year = now.year
    month = now.month
    day = now.day
    hour = now.hour
    minute = now.minute
    second = now.second
    LBLRUN_DATE="{:04d}{:02d}{:02d}".format(year, month, day)
    LBLRUN_TIME="{:02d}{:02d}{:02d}CEST".format(hour, minute, second)
    LBL_RUN_ROOT=os.getenv("HOME")
    
    # -- Create a root definition for the LBLRUN_TAG variable

    ROOT_LBLRUN_TAG="_{}_{}_{}".format(os.uname()[1], LBLRUN_DATE, LBLRUN_TIME)
    LBLRUN_TAG=ROOT_LBLRUN_TAG

    # -- Create the work directory name
    
    LBL_WORK="{}/.lblrtm{}".format(LBL_RUN_ROOT, LBLRUN_TAG)

    # -- Make sure the directory doesn't already exist. If it does, suffix
    # -- the name with an _ and an integer identifier.

    LBLRUN_SUFFIX=0
    while(os.path.exists(LBL_WORK)):
        LBLRUN_SUFFIX=LBLRUN_SUFFIX+1
        LBLRUN_TAG="{}_{}".format(ROOT_LBLRUN_TAG, LBLRUN_SUFFIX)
        LBL_WORK="{}/.lblrtm{}".format(LBL_RUN_ROOT, LBLRUN_TAG)

    # -- Make the directory

    os.mkdir("{}".format(LBL_WORK))

    # --------------------
    # Link to HITRAN files
    # --------------------
    
    # -- Link to cross-section molecule description file FSCDXS if required
    
    HITRAN_DIR="{}/hitran".format(LBL_HOME)
    if(not os.path.isfile("{}/FSCDXS".format(LBL_WORK))):
       os.symlink("{}/FSCDXS".format(HITRAN_DIR), "{}/FSCDXS".format(LBL_WORK))

    # -- Link /x (cross-section data) directory to LBL_WORK if required
    if(not os.path.isdir("{}/x".format(LBL_WORK))):
       os.symlink("{}/x".format(HITRAN_DIR), "{}/x".format(LBL_WORK))

    # -- Link /xs (cross-section data) directory to LBL_WORK if required
    if(not os.path.isdir("{}/xs".format(LBL_WORK))):
       os.symlink("{}/xs".format(HITRAN_DIR), "{}/xs".format(LBL_WORK))

    # -- Link TAPE3 spectroscopic database to work directory

    T3_DIR="{}".format(HITRAN_DIR)
    if(not os.path.exists("{}/{}".format(T3_DIR, T3_FILE))):   
       f.write("             *** NO FILE {} ***\n".format(T3_FILE))
       f.write("\n")
       f.write("     Directory of {} : \n")
       for ii in os.listdir("{}".format(T3_DIR)):
           f.write("{}\n".format(ii))
       sys.exit(3)
    
    f.write("\n")
    f.write("     Using TAPE3 file {}\n".format(T3_FILE))
    try:
        os.unlink("{}/TAPE3".format(LBL_WORK))
    except FileNotFoundError:
        pass
    os.symlink("{}/{}".format(T3_DIR, T3_FILE), "{}/TAPE3".format(LBL_WORK))

    # ---------------------------------
    # Copy TAPE5 file to work directory
    # ---------------------------------
    T5_FILE=lbltp5
    CURRENT_DIR=os.getcwd()

    if(not os.path.isfile(T5_FILE)):
        f.write("             *** NO FILE {} ***\n".format(T5_FILE))
        f.write("\n")
        f.write("     Directory of {} : \n".format(CURRENT_DIR))
        for ii in os.listdir("{}".format(CURRENT_DIR)):
            f.write("{}\n".format(ii))
        sys.exit(4)

    f.write("     Using LBL_HOME: {}\n".format(LBL_HOME))
    f.write("     Using TAPE5 file {}\n".format(T5_FILE))
    os.system("rm -f {}/'TAPE5' 2>/dev/null".format(LBL_WORK))
    shutil.copy2("{}".format(T5_FILE), "{}/TAPE5".format(LBL_WORK))

    # --------------------------------------------------
    # Copy EMISSIVITY file to work directory, if found
    # --------------------------------------------------
    EMISSION_FILE="EMISSIVITY"
    if(os.path.isfile(EMISSION_FILE)):
        f.write("             *** Found {} and copied it to run directory ***\n".format(EMISSION_FILE))
        shutil.copy2("{}".format(EMISSION_FILE), "{}/{}".format(LBL_WORK, EMISSION_FILE))

    # --------------------------------------------------
    # Copy REFLECTIVITY file to work directory, if found
    # --------------------------------------------------
    REFLECTION_FILE="REFLECTIVITY"
    if(os.path.isfile(REFLECTION_FILE)):
        f.write("             *** Found {} and copied it to run directory ***\n".format(REFLECTION_FILE))
        shutil.copy2("{}".format(REFLECTION_FILE), "{}/{}".format(LBL_WORK, REFLECTION_FILE))

    # ----------------------------
    # Run LBLRTM in work directory
    # ----------------------------

    LBLRTM='lblrtm'
    f.write("\n")
    f.write(" {} running.\n".format(LBLRTM))
    f.write(" Begin date: {}\n".format(LBLRUN_DATE))
    f.write(" Begin time: {}\n".format(LBLRUN_TIME))
    f.write("\n")
    os.system("cd {}; nice -1 {}/bin/{}".format(LBL_WORK, LBL_HOME, LBLRTM))

    # ---------------------------------------------------
    # Set and, if necessary, create the results directory
    # Modification by Tim Shippert to allow absolute path
    # ---------------------------------------------------
    
    f.write("    Using Shippert mod\n")
    OUT_DIR=lbldir
    if(OUT_DIR[0] == "/"):
        RESULTS_DIR=OUT_DIR
    else:
        RESULTS_DIR="{}/{}".format(CURRENT_DIR, OUT_DIR)

    if(not os.path.isdir(RESULTS_DIR)):
        os.mkdir("{}".format(RESULTS_DIR))
    # ------------------
    # Copy back products
    # ------------------


    # -- Determine disk usage in KB
    
    DISK_USAGE=shutil.disk_usage("{}".format(RESULTS_DIR))
    
    
    # -- Output info message
    
    f.write("\n")
    f.write(" Saving LBLRTM products ({} kB total) in {}....\n".format(DISK_USAGE, RESULTS_DIR))
    f.write("\n")
    f.close()
    # -- Copy products and then delete them individually. This will work
    # -- across filesystems, mv will not.

    move_file("{}/TAPE6".format(LBL_WORK), "{}/TAPE6".format(RESULTS_DIR), lbllog)
    move_file("{}/TAPE7".format(LBL_WORK), "{}/TAPE7".format(RESULTS_DIR), lbllog)
    move_file("{}/TAPE10".format(LBL_WORK), "{}/TAPE10".format(RESULTS_DIR), lbllog)
    move_file("{}/TAPE12".format(LBL_WORK), "{}/TAPE12".format(RESULTS_DIR), lbllog)
    #move_file("{}/TAPE13".format(LBL_WORK), "{}/TAPE13".format(RESULTS_DIR), lbllog)
    #move_file("{}/TAPE14".format(LBL_WORK), "{}/TAPE14".format(RESULTS_DIR), lbllog)
    #move_file("{}/TAPE27".format(LBL_WORK), "{}/TAPE27".format(RESULTS_DIR), lbllog)
    #move_file("{}/TAPE28".format(LBL_WORK), "{}/TAPE28".format(RESULTS_DIR), lbllog)
    #move_file("{}/EMISSIVITY".format(LBL_WORK), "{}/EMISSIVITY".format(RESULTS_DIR), lbllog)
    #move_file("{}/REFLECTIVITY".format(LBL_WORK), "{}/REFLECTIVITY".format(RESULTS_DIR), lbllog)
    #exit(-1)
    for file_ in os.listdir("{}".format(LBL_WORK)):
        if("TMP" in file_ or "OD" in file_):
            shutil.copy2("{}/{}".format(LBL_WORK, file_), "{}".format(RESULTS_DIR))
    
    # --------------------------------------------
    # Delete everything else in the work directory
    # --------------------------------------------

    for file_ in os.listdir("{}".format(LBL_WORK)):
        os.unlink("{}/{}".format(LBL_WORK, file_))
    
    # -- Delete the work directory.
    os.rmdir("{}".format(LBL_WORK))
    
    return

if __name__ == '__main__':
    lbltp5="tp5"
    lbldir="lblout"
    lbllog="lbllog.txt"
    lbl_home="/home/phi.richter/radiative_transfer/lblrtm"
    lblrun(lbltp5, lbldir, lbllog, lbl_home)

  
