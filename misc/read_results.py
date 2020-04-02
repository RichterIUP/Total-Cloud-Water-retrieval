#!/usr/bin/python

import os
import sys
import subprocess
import datetime as dt
sys.path.append("/home/phi.richter/L-IWP/src/")

import numpy as np

def get_pos(spectrum, path="/home/phi.richter/dship"):#liwp.py
  
  files = ["{}/05_24.csv", "{}/06_07.csv", "{}/06_21.csv", \
           "{}/07_05.csv", "{}/07_19.csv", "{}/08_02.csv"]

  date = spectrum.split("/")[-1].split("_")[1]
  time = spectrum.split("/")[-1].split("_")[2].split(".")[0]
  year = int(date[0:4])
  month = int(date[4:6])
  day = int(date[6:])
  hour = int(time[0:2])
  minute = int(time[2:])
  datetime = str(dt.datetime(year, month, day, hour, minute))
  for element in files:
    f = open(element.format(path), "r")
    cont = f.readlines()
    f.close()
    for line in cont:
      if(line.split(",")[0] == datetime):
        lat = float(line.split(",")[4])
        lon = float(line.split(",")[5])
        return [lat, lon]
  return [-1.0, -1.0]

def read_results(spec_date, outfile="results.csv", modifier="a", max_chi2=2000.0, firstline=False, dship="/home/phi.richter/dship"):
    runs = sorted(os.listdir(spec_date))
    h = open(outfile, modifier)
    if(firstline):
        h.write("Date/Time,Latitude (degree),Longitude (degree),Optical depth liquid (1),Error optical depth liquid (1),Optical depth ice (1), Error optical depth ice (1),Effective droplet radius liquid (um),Error effective droplet radius liquid (um),Effective droplet radius ice (um),Error effective droplet radius (um),Liquid water path (g m^(-2)),Error liquid water path (g m^(-2)),Ice water path (g m^(-2)),Error ice water path (g m^(-2)),Total water path (g m^(-2)),Error total water path (g m^(-2)),Cloud base height(m),Phase,Averaging Kernel matrix [11],Averaging Kernel matrix [12],Averaging Kernel matrix [13],Averaging Kernel matrix [14],Averaging Kernel matrix [21],Averaging Kernel matrix [22],Averaging Kernel matrix [23],Averaging Kernel matrix [24],Averaging Kernel matrix[31],Averaging Kernel matrix [32],Averaging Kernel matrix [33],Averaging Kernel matrix [34],Averaging Kernel matrix [41],Averaging Kernel matrix [42],Averaging Kernel matrix [43],Averaging Kernel matrix [44]\n")
    for folder in runs:
        run = "./{}/{}".format(spec_date, folder)
        print(run)
        #exit(-1)
        if(os.path.isdir(run) and "result_iteration" in os.listdir(run) and os.path.isdir(run)):
            f = open("{}/retrieval_log.dat".format(run), "r")
            cont = f.readlines()
            f.close()
            for ii in range(len(cont)):
                if("Diverged" in cont[ii]):
                    return
            spec = cont[5].split("/")[-1].rstrip()
            date = spec.split("_")[1]
            time = spec.split("_")[2].split(".")[0]
            #cb = float(cont[31].split(":")[1].rstrip().lstrip())
            #print(run)
            #if("2" in cont[-18].split("\t")[-1]):#Letzte Zeile mit Iteration
            #    return 
            date_iso = "{}-{}-{}T{}:{}:00".format(date[0:4], date[4:6], date[6:], time[0:2], time[2:])
            #print(spec)
            [lat, lon] = get_pos(spec, dship)
            f = open("{}/result_iteration".format(run), "r")
            cont = f.readlines()
            #for ii in range(len(cont)):
            #    print(ii, cont[ii])
            #exit(-1)
            # 4 tau_liq
            # 5 tau_ice
            # 6 ref_liq
            # 7 ref_ice
            # 8 LWP
            # 9 IWP
            # 21 chi2
            f.close()
            try:
                avk_00 = float(cont[14].split("\t")[0])
                avk_01 = float(cont[14].split("\t")[1])
                avk_02 = float(cont[14].split("\t")[2])
                avk_03 = float(cont[14].split("\t")[3])
                avk_10 = float(cont[15].split("\t")[0])
                avk_11 = float(cont[15].split("\t")[1])
                avk_12 = float(cont[15].split("\t")[2])
                avk_13 = float(cont[15].split("\t")[3])
                avk_20 = float(cont[16].split("\t")[0])
                avk_21 = float(cont[16].split("\t")[1])
                avk_22 = float(cont[16].split("\t")[2])
                avk_23 = float(cont[16].split("\t")[3])
                avk_30 = float(cont[17].split("\t")[0])
                avk_31 = float(cont[17].split("\t")[1])
                avk_32 = float(cont[17].split("\t")[2])
                avk_33 = float(cont[17].split("\t")[3])
            except IndexError:
                print(cont)
                print(spec_date)
                continue
            t_l_retr = [0.0, 0.0]
            t_i_retr = [0.0, 0.0]
            r_l_retr = [0.0, 0.0]
            r_i_retr = [0.0, 0.0]
            LWP_retr = [0.0, 0.0]
            IWP_retr = [0.0, 0.0]
            TWP_retr = [0.0, 0.0]
            try:
                t_l_retr[0] = float(cont[4].split(" ")[-3][1:])
                t_l_retr[1] = float(cont[4].split(" ")[-1].rstrip()[:-1])
                t_i_retr[0] = float(cont[5].split(" ")[-3][1:])
                t_i_retr[1] = float(cont[5].split(" ")[-1].rstrip()[:-1])
                r_l_retr[0] = float(cont[6].split(" ")[-3][1:])
                r_l_retr[1] = float(cont[6].split(" ")[-1].rstrip()[:-1])
                r_i_retr[0] = float(cont[7].split(" ")[-3][1:])
                r_i_retr[1] = float(cont[7].split(" ")[-1].rstrip()[:-1])
                LWP_retr[0] = float(cont[8].split(" ")[-3][1:])
                LWP_retr[1] = float(cont[8].split(" ")[-1].rstrip()[:-1])
                IWP_retr[0] = float(cont[9].split(" ")[-3][1:])
                IWP_retr[1] = float(cont[9].split(" ")[-1].rstrip()[:-1])
                TWP_retr[0] = LWP_retr[0] + IWP_retr[0]
                TWP_retr[1] = LWP_retr[1] + IWP_retr[1]
            except ValueError:
                print(cont)
                print(spec_date)
                continue
            chi2 = float(cont[21].split("\t")[0].split("=")[-1].rstrip())
            if(chi2 <= max_chi2):
                for fl in os.listdir(run):
                  #print(date_iso)
                  if("iter_final_" in fl):
                    os.system("cp {}/{} plots/{}_{}".format(run, fl, spec_date, fl))                        
                    #if("_abs.svg" in fl and not "iter" in fl):
                    #print(date_iso)
                    #print("h.write")
                    #os.system("cp {}/{} plots/{}_{}".format(run, fl, spec_date, fl))
                    h.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(date_iso, lat, lon, \
                                                                                                                                          t_l_retr[0], t_l_retr[1], \
                                                                                                                                          t_i_retr[0], t_i_retr[1], \
                                                                                                                                          r_l_retr[0], r_l_retr[1], \
                                                                                                                                          r_i_retr[0], r_i_retr[1], \
                                                                                                                                          LWP_retr[0], LWP_retr[1], \
                                                                                                                                          IWP_retr[0], IWP_retr[1], \
                                                                                                                                          TWP_retr[0], TWP_retr[1], \
                                                                                                                                          avk_00, avk_01, avk_02, avk_03, \
                                                                                                                                          avk_10, avk_11, avk_12, avk_13, \
                                                                                                                                          avk_20, avk_21, avk_22, avk_23, \
                                                                                                                                          avk_30, avk_31, avk_32, avk_33))
            else:
              print(chi2, spec_date)        

    h.close()

if __name__ == '__main__':
    subprocess.call(["mkdir", "plots"])
    folders = sorted(os.listdir("."))
    x = True
    for folder in folders:
        y = "./{}".format(folder)
        read_results(spec_date=y, firstline=x, max_chi2=1750.0)
        if(x):
            x = not x
