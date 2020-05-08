#!/usr/bin/python3

import numpy as np
import inversion
import inp
import threading

import log

SEARCH_APR_RMS = []
SEARCH_APR_MCP = []

def guess_apr(tt, fi, rl, ri):

    global SEARCH_APR
    #ri = inp.MCP[2]
    #rl = inp.MCP[3]
    #fi = 0.5
    #rms = []
    #tt_y = []
    #rt_y = []

    #tl = tt*(1-fi)
    #ti = tt*fi
    rms = inversion.__only_fwd(tau_liq=tl, tau_ice=ti, reff_liq=rl, reff_ice=ri, filenum=int(100*tt+100*fi+100*rl+100*ri))[-2]
    #tt_y.append(tt)
            
    #tt_best = tt

    #tl_best = tl_best = tt_best * (1-fi)
    #ti_best = ti_best = tt_best * fi        

    
    lock = threading.Lock()    
    lock.acquire()
        
    log.write("{} {}".format(rms, [tl, ti, rl, ri]))
    SEARCH_APR_RMS.append(rms)
    SEARCH_APR_MCP.append([tl, ti, rl, ri])
    lock.release()

    return
