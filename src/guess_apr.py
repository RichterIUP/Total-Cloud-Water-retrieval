#!/usr/bin/python3

import numpy as np
import inversion
import inp

def guess_apr(fi):
    inp.FORWARD = True

    tt = [0.2, 1.0, 3.0, 4.0, 6.0]
    rt = [10, 15, 20, 25, 30, 35, 40, 45]

    rms = []

    for param_num in range(len(tt)):
        tl = tt[param_num]*(1-fi)
        ti = tt[param_num]*fi
        rms.append(inversion.__only_fwd(tau_liq=tl, tau_ice=ti, reff_liq=10., reff_ice=30.)[4])
 
    tt_best = tt[rms.index(min(rms))]
    rms = []
    tl_best = tt_best * (1-fi)
    ti_best = tt_best * fi        

    for rl in [5, 8, 11, 14, 17, 20]:
        for ri in [10, 20, 30, 40, 50]:
            rms.append(inversion.retrieve(tau_liq=tl, tau_ice=ti, reff_liq=rl, reff_ice=ri)[4])
            rt_y.append([rl, ri])
            with open("radii_{}".format(fi), "a") as f:
                f.write("{} {} {} {}\n".format(fi, tt_best, rms[-1], [rl, ri]))
    idx = rms.index(min(rms))
    with open("radii_{}".format(fi), "a") as f:
        f.write("{} {}\n".format(fi, tt_best, rms[idx], rt_y[idx]))
        
    #inp.FOWARD = False
    exit(-1)

