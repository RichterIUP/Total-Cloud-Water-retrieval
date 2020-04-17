#!/usr/bin/python3

import numpy as np
import inversion
import inp
import threading

SEARCH_APR = []

def guess_apr(fi):

    global SEARCH_APR
    
    tt = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
    rl = 10.
    ri = 30.

    rms = []
    rad_lbldis = []
    rad_ftir = []
    tt_y = []
    rt_y = []

    for param_num in range(len(tt)):
        tl = tt[param_num]*(1-fi)
        ti = tt[param_num]*fi
        rms.append(inversion.__only_fwd(tau_liq=tl, tau_ice=ti, reff_liq=rl, reff_ice=ri)[-1])
        tt_y.append(tt[param_num])
    
    idx = rms.index(min(rms))
    tt_best = tt_y[idx]

    rms = []
    tl_best = tt_best * (1-fi)
    ti_best = tt_best * fi        

    for rl in [5, 8, 11, 14, 17, 20]:
        for ri in [10, 20, 30, 40, 50]:
            rms.append(inversion.__only_fwd(tau_liq=tl_best, tau_ice=ti_best, reff_liq=rl, reff_ice=ri)[-1])
            rt_y.append([rl, ri])
            #with open("radii_{}".format(fi), "a") as f:
            #    f.write("{} {} {} {}\n".format(fi, tt_best, rms[-1], rt_y[-1]))
    idx = rms.index(min(rms))
    with open("radii_{}".format(fi), "a") as f:
        f.write("{} {}\n".format(fi, tt_best, rms[idx], rt_y[idx]))
        
    lock = threading.Lock()
    
    lock.acquire()
    SEARCH_APR.append([rms[idx], [fi, tt_best, rt_y[idx][0], rt_y[idx][1]])
    lock.release()
