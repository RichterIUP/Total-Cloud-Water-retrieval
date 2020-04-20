#!/usr/bin/python3

import scipy.io as sio
import sys
import numpy as np
import os

#fname = sys.argv[1]

def read_csv(fname):
    with sio.netcdf_file(fname, "r") as f:
        tt = np.array([f.variables['x_ret'][0], f.variables['x_ret_err'][0], f.variables['x_err_res'][0]])
        fi = np.array([f.variables['x_ret'][1], f.variables['x_ret_err'][1], f.variables['x_err_res'][1]])
        rl = np.array([f.variables['x_ret'][2], f.variables['x_ret_err'][2], f.variables['x_err_res'][2]])
        ri = np.array([f.variables['x_ret'][3], f.variables['x_ret_err'][3], f.variables['x_err_res'][3]])
        lwp = np.array([f.variables['wp_ret'][0], f.variables['wp_ret_err'][0], f.variables['wp_err_res'][0]])
        iwp = np.array([f.variables['wp_ret'][1], f.variables['wp_ret_err'][1], f.variables['wp_err_res'][1]])
        twp = np.array([f.variables['wp_ret'][2], f.variables['wp_ret_err'][2], f.variables['wp_err_res'][2]])
        ctemp = np.array(f.variables["av_ctemp"][:])
        res = f.variables['residuum'][:].copy()
        wn = f.variables['wavenumber'][:].copy()
        rms = np.sqrt(np.mean(np.array(res)**2))
        conv = f.variables['conv'][:].copy()
        pwv = f.variables['pwv'][:].copy()
        cloud = np.array([np.float_(f.variables['cloud_base'][:]), np.float_(f.variables['cloud_top'][:])])
        x_a = np.array(f.variables['x_a'][:].copy())
        x_a_err = np.array(f.variables['x_a_err'][:].copy())
        
        
        print("Results:")
        print("Filename: {}".format(fname))
        print("Average Cloud temperature (K): {}\n".format(ctemp[0]))
        print("Cloud base height (m): {}".format(cloud[0]))
        print("Cloud top height (m): {}\n".format(cloud[1]))
        print("Converged: {}\n".format(conv))
        print("DIRECT PRODUCTS")
        print("Liquid Optical Depth (1): ({} +- {})".format(np.float_(tt[0]), np.float_(tt[1])))
        print("Ice Optical Depth (1): ({} +- {})".format(np.float_(fi[0]), np.float_(fi[1])))
        print("Liquid Radius (um): ({} +- {})".format(np.float_(rl[0]), np.float_(rl[1])))
        print("Ice Radius (um): ({} +- {})".format(np.float_(ri[0]), np.float_(ri[1])))
        print("Root-Mean-Squared Error (mW/sr (cm-1) m2): {}\n".format(rms))
        print("DERIVED PRODUCTS")
        print("Liquid Water Path (g/m2): ({} +- {})".format(np.float_(lwp[0]), np.float_(lwp[1])))
        print("Ice Water Path (g/m2): ({} +- {})".format(np.float_(iwp[0]), np.float_(iwp[1])))
        print("Total Water Path (g/m2): ({} +- {})".format(np.float_(twp[0]), np.float_(twp[1])))

    return [tt, fi, rl, ri, lwp, iwp, twp, rms, ctemp, pwv, conv, cloud, x_a, x_a_err, wn, res]

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
                dtt_res = np.float_(out[0][2])
                fi = np.float_(out[1][0])
                dfi = np.float_(out[1][1])
                dfi_res = np.float_(out[1][2])
                rl = np.float_(out[2][0])
                drl = np.float_(out[2][1])
                drl_res = np.float_(out[2][2])
                ri = np.float_(out[3][0])
                dri = np.float_(out[3][1])
                dri_res = np.float_(out[3][2])
                lwp = np.float_(out[4][0])
                dlwp = np.float_(out[4][1])
                dlwp_res = np.float_(out[4][2])
                iwp = np.float_(out[5][0])
                diwp = np.float_(out[5][1])
                diwp_res = np.float_(out[5][2])
                twp = np.float_(out[6][0])
                dtwp = np.float_(out[6][1])
                dtwp_res = np.float(out[6][2])
                rms = out[7]
                ctemp = out[8][0]
                pwv = out[9][0]
                conv = out[10][0]
                cbase = out[11][0]
                ctop = out[11][1]
                f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(element, date, tt, dtt, dtt_res, fi, dfi, dfi_res, rl, drl, drl_res, ri, dri, dri_res, lwp, dlwp, dlwp_res, iwp, diwp, diwp_res, twp, dtwp, dtwp_res, rms, ctemp, pwv, conv, cbase, ctop))
            #except Exception:
            #    pass
        f.close()
                
                
