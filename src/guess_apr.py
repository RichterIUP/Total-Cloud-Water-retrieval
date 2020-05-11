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

    tl=tt*(1-fi)
    ti=tt*fi
    rms = inversion.__only_fwd(tau_liq=tl, tau_ice=ti, reff_liq=rl, reff_ice=ri, filenum=int(tt+fi+rl+10.*ri))
       
    lock = threading.Lock()    
    lock.acquire()
        
    log.write("{} {}".format(rms, [tl, ti, rl, ri]))
    SEARCH_APR_RMS.append(rms)
    SEARCH_APR_MCP.append([tl, ti, rl, ri])
    lock.release()

    return
