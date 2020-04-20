#!/usr/bin/python3

import numpy as np
import inversion
import inp
import threading

SEARCH_APR_RMS = []
SEARCH_APR_MCP = []

def guess_apr(fi):

    global SEARCH_APR
    
    tt = [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
    rl = 10.
    ri = 30.

    rms = []
    rad_lbldis = []
    rad_ftir = []
    tt_y = []
    rt_y = []

    for tt in [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]:
        tl = tt*(1-fi)
        ti = tt*fi
        rms.append(inversion.__only_fwd(tau_liq=tl, tau_ice=ti, reff_liq=rl, reff_ice=ri, filenum=int(10*fi))[-2])
        tt_y.append(tt)
            
    idx = rms.index(min(rms))
    tt_best = tt_y[idx]

    slope = []
    tl_best = tt_best * (1-fi)
    ti_best = tt_best * fi        

    for rl in [5, 10, 15, 20]:
        for ri in [20, 25, 30, 35, 40, 45, 50]:
            slope.append(inversion.__only_fwd(tau_liq=tl_best, tau_ice=ti_best, reff_liq=rl, reff_ice=ri, filenum=int(10*fi))[-1])
            rt_y.append([rl, ri])

    idx = slope.index(min(slope))

        
    lock = threading.Lock()    
    lock.acquire()
    SEARCH_APR_RMS.append(rms[idx])
    SEARCH_APR_MCP.append([fi, tt_best, rt_y[idx][0], rt_y[idx][1]])
    lock.release()

    return
