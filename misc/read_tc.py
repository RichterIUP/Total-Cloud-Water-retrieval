#!/usr/bin/python

import sys
import os
import numpy as np
sys.path.append("/home/phi.richter/L-IWP/src/")
import read_database as mie

def find_path1(path_):
    folders = []
    for spectra_folder in sorted(os.listdir(path_)):
        inner_folder = "{}/{}".format(path_, spectra_folder)
        #print("inner_folder", inner_folder)
        #print(os.path.isdir(inner_folder))
        if(not os.path.isdir(inner_folder)):
            break
        for case in sorted(os.listdir(inner_folder)):
            case_folder = "{}/{}".format(inner_folder, case)
          #for case2 in os.listdir(case_folder):
          #  case2_folder = "{}/{}".format(case_folder, case2)
            #print(case_folder)
            if(os.path.exists("{}/result_iteration".format(case_folder))):
                folders.append(case_folder)
            else:
                break
    return folders
      
def read_log(path_, tc):
    for spectra in os.listdir(path_):
        spectra_folder = "{}/{}".format(path_, spectra)
        if(tc == spectra):
          f = open("{}/log".format(spectra_folder), "r")
          cont = f.readlines()
          f.close()
          tt = np.float64(cont[0].split(" ")[-1])
          fi = np.float64(cont[1].split(" ")[-1])
          rl = np.float64(cont[2].split(" ")[-1])
          ri = np.float64(cont[3].split(" ")[-1])
          pwv = np.float64(cont[4].split(" ")[-1])
          if(np.isnan(rl)):
            rl = 7.0
          if(np.isnan(ri)):
            ri = 21.0
          return [tt, fi, rl, ri, pwv]
    
def read_param(path_testcases, path_spectra, outfiles, outpaths, shape):

    for folder in sorted(find_path1(path_testcases)):
        print(folder.split("/")[-1])
        #Read theoretical results
        f = open("{}/retrieval_log.dat".format(folder), "r")
        cont = f.readlines()
        f.close()
        resolution = 0.15
        idx = 1
        num_wn = 0.0
        for ii in range(len(cont)):
            if "Microwindow" in cont[ii]:
                wn_down = float(cont[ii].split(":")[-1].split(",")[0].split("[")[1].rstrip())
                print(wn_down)
                if wn_down < 700.0:
                    outfile = outfiles[0]
                    outpath = outpaths[0]
                else:
                    outfile = outfiles[1]
                    outpath = outpaths[1]
                if not os.path.exists("{}".format(outfile)):
                    g = open("{}".format(outfile), "w")
                    g.write("spectrum, tt_sim, fi_sim, rl_sim, ri_sim, pwv, LWP_sim, IWP_sim, TWP_sim, tt, dtt, fi, dfi, rl, drl, ri, dri, LWP, dLWP, IWP, dIWP, TWP, dTWP, chi2, conv\n")
                else:
                    g = open("{}".format(outfile), "a")
                if not os.path.exists("{}/plots_{}".format(outpath, shape)):
                    os.system("mkdir {}/plots_{}".format(outpath, shape))
                break

        for ii in range(len(cont)):
            if "Not converged!" in cont[ii]:
                conv = 0
                break
            elif "Final" in cont[ii]:
                conv = 1
                break
            
        spectrum = "simulation_{}".format(cont[5].split("simulation_")[1].split("/")[0])
        [tt_sim, fi_sim, rl_sim, ri_sim, pwv] = read_log(path_spectra, spectrum)
        #Read retrieval results
        f = open("{}/result_iteration".format(folder), "r")
        cont = f.readlines()
        f.close()
        
        
        tl = np.float64(cont[7].split("(")[-1].split("+-")[0])
        dtl = np.float64(cont[7].split("(")[-1].split("+-")[1].split(")")[0])
        ti = np.float64(cont[8].split("(")[-1].split("+-")[0])
        dti = np.float64(cont[8].split("(")[-1].split("+-")[1].split(")")[0])
        rl = np.float64(cont[9].split("(")[-1].split("+-")[0])
        drl = np.float64(cont[9].split("(")[-1].split("+-")[1].split(")")[0])
        ri = np.float64(cont[10].split("(")[-1].split("+-")[0])
        dri = np.float64(cont[10].split("(")[-1].split("+-")[1].split(")")[0])
        LWP = np.float64(cont[11].split("(")[-1].split("+-")[0])
        dLWP = np.float64(cont[11].split("(")[-1].split("+-")[1].split(")")[0])
        IWP = np.float64(cont[12].split("(")[-1].split("+-")[0])
        dIWP = np.float64(cont[12].split("(")[-1].split("+-")[1].split(")")[0])
        chi2 = np.float64(cont[24].split("=")[-1])
        #noise = np.float64(cont[25].split("=")[-1])
        #meandiff = np.float64(cont[26].split("=")[-1])
        #meandiff_per_win = meandiff / num_wn
        #maxdiff = np.float64(cont[27].split("=")[-1].split("[")[2].split("]")[0])
        TWP = LWP + IWP
        dTWP = dLWP + dIWP
        tt = tl+ti
        dtt = dtl+dti
        fi = ti/tt
        dfi = np.abs(dti/tt) + np.abs(ti/(tt**2) *dtt)
        #if meandiff > 1.0 and chi2 > 2000.0:
        #    conv = 1
        #elif meandiff > 1.0 or chi2 > 2000.0:
        #    conv = 2
        #elif rl_sim > 30.0 or ri_sim > 30.0:
        #    conv = 3
        #else:
        #    conv = 0
        
        
        [LWP_sim, IWP_sim, TWP_sim] = calc_theor_wp(tt_sim, fi_sim, rl_sim, ri_sim, spectrum)
        os.system("cp {}/tl_* {}/plots_{}/{}.svg".format(folder, outpath, shape, spectrum))          
        g.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(spectrum, tt_sim, fi_sim, rl_sim, ri_sim, pwv, LWP_sim, IWP_sim, TWP_sim, \
                                                                                                               tt, dtt, fi, dfi, rl, drl, ri, dri, LWP, dLWP, IWP, \
                                                                                                               dIWP, TWP, dTWP, chi2, conv))
        g.close()
        #exit(-1)
    
def calc_theor_wp(tt_sim, fi_sim, rl_sim, ri_sim, spectrum):
    path = "/home/phi.richter/L-IWP/ssp"#"/home/philippr/Arbeit/RETRIEVAL/L-IWP/misc"
    if("R000HBR" in spectrum):
        db_ice = "{}/ssp_db.BulletRosette.gamma.0p100".format(path)
    elif("R000P" in spectrum or "R050P" in spectrum):
        db_ice = "{}/ssp_db.Plate.gamma.0p100".format(path)
    elif("R000SC" in spectrum or "R050SC" in spectrum):
        db_ice = "{}/ssp_db.SolidCol.gamma.0p100".format(path)
    else:
        db_ice = "{}/ssp_db.mie_ice.gamma_sigma_0p100".format(path)

    [liq, ice] = mie.read_databases("{}/ssp_db.mie_wat.gamma_sigma_0p100".format(path), db_ice)

    t_t = tt_sim
    f_i = fi_sim
    r_l = rl_sim
    r_i = ri_sim
    IWP = mie.calc_iwp(t_t*f_i, 0.0, r_i, 0.0, ice)[0]
    LWP = mie.calc_lwp(r_l, 0.0, t_t - t_t * f_i, 0.0)[0]
    return [LWP, IWP, LWP+IWP]
  
if __name__ == '__main__':
    path_total = "/mnt/beegfs/user/phi.richter/"
    path_outfolder = "{}/OUTFOLDER".format(path_total)#sys.argv[1]
    path_testcases = "/home/phi.richter/TESTCASES"#sys.argv[2]
    outfile = [None, None]
    outpath = [None, None]
    outfile[0] = "/home/phi.richter/RES_FIR/outfile.csv".format(path_total)#sys.argv[3]
    outpath[0] = "/home/phi.richter/RES_FIR".format(path_total)
    outfile[1] = "/home/phi.richter/RES_TIR/outfile.csv".format(path_total)#sys.argv[3]
    outpath[1] = "/home/phi.richter/RES_TIR".format(path_total)    
    shape = "COR"
    read_param(path_outfolder, path_testcases, outfile, outpath, shape=shape)

    #path_outfolder = "{}/OUTFOLDER_tc_TIR".format(path_total)#sys.argv[1]
    #path_testcases = "/home/phi.richter/TESTCASES"#sys.argv[2]
    #outfile = "{}/L-IWP/RES_TIR/outfile.csv".format(path_total)#sys.argv[3]
    #outpath = "{}/L-IWP/RES_TIR".format(path_total)
    #print(len(sys.argv))
    #shape = "COR"
    #read_param(path_outfolder, path_testcases, outfile, outpath, shape=shape)  


