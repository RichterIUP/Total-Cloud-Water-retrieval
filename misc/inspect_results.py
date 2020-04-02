#!/usr/bin/python3

import scipy.io as sio
import sys
import numpy as np
import os

#fname = sys.argv[1]

def read_csv(fname):
    with sio.netcdf_file(fname, "r") as f:
        tt = np.array([f.variables['tt'][:], f.variables['dtt'][:]])
        fi = np.array([f.variables['fi'][:], f.variables['dfi'][:]])
        rl = np.array([f.variables['rl'][:], f.variables['drl'][:]])
        ri = np.array([f.variables['ri'][:], f.variables['dri'][:]])
        lwp = np.array([f.variables['lwp'][:], f.variables['dlwp'][:]])
        iwp = np.array([f.variables['iwp'][:], f.variables['diwp'][:]])
        twp = np.array([f.variables['twp'][:], f.variables['dtwp'][:]])
        ctemp = np.array(f.variables["av_ctemp"][:])
        rms = np.sqrt(np.mean(np.array(f.variables['residuum'][:])**2))
        conv = f.variables['conv'][:].copy()
        pwv = f.variables['pwv'][:].copy()
        cloud = np.array([np.float_(f.variables['cloud_base'][:]), np.float_(f.variables['cloud_top'][:])])
        
        print("Results:")
        print("Filename: {}".format(fname))
        print("Average Cloud temperature (K): {}\n".format(ctemp[0]))
        print("Cloud base height (m): {}".format(cloud[0]))
        print("Cloud top height (m): {}\n".format(cloud[1]))
        print("Converged: {}\n".format(conv))
        print("DIRECT PRODUCTS")
        print("Total Optical Depth (1): ({} +- {})".format(np.float_(tt[0]), np.float_(tt[1])))
        print("Ice Fraction (1): ({} +- {})".format(np.float_(fi[0]), np.float_(fi[1])))
        print("Liquid Radius (um): ({} +- {})".format(np.float_(rl[0]), np.float_(rl[1])))
        print("Ice Radius (um): ({} +- {})".format(np.float_(ri[0]), np.float_(ri[1])))
        print("Root-Mean-Squared Error (mW/sr (cm-1) m2): {}\n".format(rms))
        print("DERIVED PRODUCTS")
        print("Liquid Water Path (g/m2): ({} +- {})".format(np.float_(lwp[0]), np.float_(lwp[1])))
        print("Ice Water Path (g/m2): ({} +- {})".format(np.float_(iwp[0]), np.float_(iwp[1])))
        print("Total Water Path (g/m2): ({} +- {})".format(np.float_(twp[0]), np.float_(twp[1])))

    return [tt, fi, rl, ri, lwp, iwp, twp, rms, ctemp, pwv, conv, cloud]

if __name__ == '__main__':
    fname = sys.argv[1]
    if not os.path.isdir(fname):
        read_csv(fname)
    else:
        files = sorted(os.listdir(fname))
        f = open("out.csv", "w")
        f.write("fname,date,tt,dtt,fi,dfi,rl,drl,ri,dri,lwp,dlwp,iwp,diwp,twp,dtwp,rms,ctemp,pwv,conv,cbase,ctop\n")
        for element in files:
            #try:
            if True:
                try:
                    date = "{}-{}-{} {}:{}".format(element.split("PS.")[1][0:4], element.split("PS.")[1][4:6], element.split("PS.")[1][6:8], element.split("PS.")[1][9:11], element.split("PS.")[1][11:13])
                except IndexError:
                    date = fname
                out = read_csv("{}/{}".format(fname, element))
                tt = np.float_(out[0][0])
                dtt = np.float_(out[0][1])
                fi = np.float_(out[1][0])
                dfi = np.float_(out[1][1])
                rl = np.float_(out[2][0])
                drl = np.float_(out[2][1])
                ri = np.float_(out[3][0])
                dri = np.float_(out[3][1])
                lwp = np.float_(out[4][0])
                dlwp = np.float_(out[4][1])
                iwp = np.float_(out[5][0])
                diwp = np.float_(out[5][1])
                twp = np.float_(out[6][0])
                dtwp = np.float_(out[6][1])
                rms = out[7]
                ctemp = out[8][0]
                pwv = out[9][0]
                conv = out[10][0]
                cbase = out[11][0]
                ctop = out[11][1]
                f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(element, date, tt, dtt, fi, dfi, rl, drl, ri, dri, lwp, dlwp, iwp, diwp, twp, dtwp, rms, ctemp, pwv, conv, cbase, ctop))
            #except Exception:
            #    pass
        f.close()
                
                
