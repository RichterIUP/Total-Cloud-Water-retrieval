import scipy.io as sio
import sys
import numpy as np
import os
import sys
sys.path.append("/home/philippr/Seafile/PhD/RETRIEVAL/L-IWP/src")
import read_database as mie

fname = sys.argv[1]

path_to_tc = "/home/philippr/Seafile/PhD/RETRIEVAL/TESTCASES_ASCII"

f = open("spectra.dat", "r")
content = f.readlines()
f.close()
out = open("results.csv", "w")
out.write("specname,tt,dtt,fi,dfi,rl,drl,ri,dri,lwp,dlwp,iwp,diwp,twp,dtwp,ctemp,rms\n")
g = open("tc.csv", "w")
g.write("specname,tt,fi,rl,ri,lwp,iwp,twp\n")
g.close()
for element in sorted(content):
    with sio.netcdf_file(element.rstrip(), "r") as f:
        tt = np.array([f.variables['tt'][:], f.variables['dtt'][:]])
        fi = np.array([f.variables['fi'][:], f.variables['dfi'][:]])
        rl = np.array([f.variables['rl'][:], f.variables['drl'][:]])
        ri = np.array([f.variables['ri'][:], f.variables['dri'][:]])
        lwp = np.array([f.variables['lwp'][:], f.variables['dlwp'][:]])
        iwp = np.array([f.variables['iwp'][:], f.variables['diwp'][:]])
        twp = np.array([f.variables['twp'][:], f.variables['dtwp'][:]])
        ctemp = np.array(f.variables["av_ctemp"][:])
        rms = np.sqrt(np.mean(np.array(f.variables['residuum'][:])**2))
        
    outstring = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n"
    out.write(outstring.format(element.rstrip(), np.float_(tt[0]), np.float_(tt[1]), np.float_(fi[0]), np.float_(fi[1]), \
                            np.float_(rl[0]), np.float_(rl[1]), np.float_(ri[0]), np.float_(ri[1]), np.float_(lwp[0]), np.float_(lwp[1]), \
                            np.float_(iwp[0]), np.float_(iwp[1]), np.float_(twp[0]), np.float_(twp[1]), np.float_(ctemp), rms))  
                            
    for tc in os.listdir(path_to_tc):
        
        if tc in element.rstrip():
            with open("{}/{}/log".format(path_to_tc, tc), "r") as g:
                tc_cont = g.readlines()
                tc_tt = float(tc_cont[0].split(" ")[-1])
                tc_fi = float(tc_cont[1].split(" ")[-1]) 
                tc_rl = float(tc_cont[2].split(" ")[-1])
                tc_ri = float(tc_cont[3].split(" ")[-1])
                [liq, ice] = mie.read_databases("/home/philippr/Seafile/PhD/RETRIEVAL/L-IWP//ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                            "/home/philippr/Seafile/PhD/RETRIEVAL/L-IWP//ssp/ssp_db.mie_ice.gamma_sigma_0p100")
                try:
                    if tc_fi == 1.0:
                        lwp = 0.0
                        tc_rl = -1
                    else:
                        lwp = np.float_(mie.calc_lwp(tc_rl, 0.0, tc_tt*(1-tc_fi), 0.0)[0])
                    if tc_fi == 0.0:
                        iwp = 0.0
                        tc_ri = -1
                    else:
                        iwp = np.float_(mie.calc_iwp(tc_tt*tc_fi, 0.0, tc_ri, 0.0, ice)[0])
                except IndexError:
                    lwp = -1
                    iwp = -1
            with open("tc.csv", "a") as g:
                g.write("{},{},{},{},{},{},{},{}\n".format(tc, tc_tt, tc_fi, tc_rl, tc_ri, lwp, iwp, lwp+iwp))
            break
out.close()
    
    
