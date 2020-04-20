#!/usr/bin/python3

import numpy as np
import inversion
import inp
import threading

import log

SEARCH_APR_RMS = []
SEARCH_APR_MCP = []

def guess_apr(tt):

    global SEARCH_APR

    ri = inp.MCP[2]
    rl = inp.MCP[3]
    fi = 0.5
    rms = []
    tt_y = []
    rt_y = []

    tl = tt*(1-fi)
    ti = tt*fi
    rms.append(inversion.__only_fwd(tau_liq=tl, tau_ice=ti, reff_liq=rl, reff_ice=ri, filenum=int(10*ri))[-2])
    tt_y.append(tt)
            
    idx = rms.index(min(rms))
    tt_best = tt_y[idx]

    #rms = []
    tl_best = tl_best = tt_best * (1-fi)
    ti_best = ti_best = tt_best * fi        
    
    #for rl in [5, 10, 15, 20]:
    #    rms.append(inversion.__only_fwd(tau_liq=tl_best, tau_ice=ti_best, reff_liq=rl, reff_ice=ri, filenum=int(10*fi))[-2])
    #    rt_y.append([rl, ri])

    #idx = rms.index(min(slope))
    #ri = rt_y[idx]
    
    lock = threading.Lock()    
    lock.acquire()
        
    log.write("{} {}".format(rms[idx], [tl_best, ti_best, rl, ri]))
    SEARCH_APR_RMS.append(rms[idx])
    SEARCH_APR_MCP.append([tl_best, ti_best])
    lock.release()

    return
