#!/usr/bin/python3

import sys
import numpy as np

R_S = 461.5#J/kg*K

def help():
    print("Usage: python3 convert_humidity.py [-OPTIONS] [VALUE] -t TEMPERAUTRE -p PRESSURE")
    print("OPTIONS:")
    print("-a: Absolute humidity (g/m3)")
    print("-r: Relative humidity (%)")
    print("-s: Specifig humidity (g/kg)")
    print("-t: Temperature (K)")
    print("-p: Pressure (hPa)")
    
def convert(par_1, par_2, par_3, verbose=False):
    if par_1[0] == "-a" or par_1[0] == "-r" or par_1[0] == "-s":
        humd = par_1
        humd[1] = humd[1]
    elif par_2[0] == "-a" or par_2[0] == "-r" or par_2[0] == "-s":
        humd = par_2
    elif par_3[0] == "-a" or par_3[0] == "-r" or par_3[0] == "-s":
        humd = par_3
        
    if par_1[0] == "-tc" or par_1[0] == "-tk":
        temp = par_1
    elif par_2[0] == "-tk" or par_2[0] == "-tc":
        temp = par_2
    elif par_3[0] == "-tk" or par_3[0] == "-tc":
        temp = par_3 
              
    if temp[0] == '-tc': 
        temp[1] = temp[1] + 273.15
        
    if par_1[0] == "-p":
        pres = par_1
    elif par_2[0] == "-p":
        pres = par_2
    elif par_3[0] == "-p":
        pres = par_3
    pres[1] = pres[1] * 100
    
    max_vapour_pressure = np.exp(-6094.4642 * temp[1]**(-1) + \
                        21.1249952 - \
                        2.7245552e-2 * temp[1] + \
                        1.6853396e-5 * temp[1]**2 + \
                        2.4575506 * np.log(temp[1]))
                        
    if humd[0] == "-s":
        humd[1] = humd[1] / 1000.0
        ah = __spec_to_abs(humd[1], temp[1], pres[1], max_vapour_pressure)
        rh = __spec_to_rel(humd[1], temp[1], pres[1], max_vapour_pressure)
        if verbose:
            print("Temperature: {} K".format(temp[1]))
            print("Pressure: {} hPa".format(pres[1]*0.01))
            print("Specifig Humidity: {} g/kg\n".format(humd[1]*1000))
            print("Absolute Humidity: {} g/m3".format(ah*1000)) 
            print("Relative Humidity: {} %".format(rh))
        return [ah, rh]
    elif humd[0] == "-a":
        humd[1] = humd[1] / 1000.0
        rh = __abs_to_rel(humd[1], temp[1], pres[1], max_vapour_pressure)
        q = __abs_to_spec(humd[1], temp[1], pres[1], max_vapour_pressure)
        if verbose:
            print("Temperature: {} K".format(temp[1]))
            print("Pressure: {} hPa".format(pres[1]*0.01))
            print("Absolute Humidity: {} g/m3\n".format(humd[1]*1000)) 
            print("Specific Humidity: {} g/kg".format(q*1000))
            print("Relative Humidity: {} %".format(rh))
        return [rh, q]
    elif humd[0] == "-r":
        ah = __rel_to_abs(humd[1], temp[1], pres[1], max_vapour_pressure)
        q = __rel_to_spec(humd[1], temp[1], pres[1], max_vapour_pressure)
        if verbose:
            print("Temperature: {} K".format(temp[1]))
            print("Pressure: {} hPa".format(pres[1]*0.01))
            print("Relative Humidity: {} %\n".format(humd[1]))
            print("Absolute Humidity: {} g/m3".format(ah*1000)) 
            print("Specific Humidity: {} g/kg".format(q*1000))
        return [ah, q]
        
def __spec_to_abs(humd, temp, pres, max_vapour_pressure):
    global R_S
    mixing_ratio = ((1.0/humd) - 1)**(-1)
    vapour_pressure = mixing_ratio * pres / (mixing_ratio + 0.622)
    ah = vapour_pressure / (R_S * temp)
    return ah
    
def __spec_to_rel(humd, temp, pres, max_vapour_pressure):
    mixing_ratio = ((1.0/humd) - 1)**(-1)
    vapour_pressure = mixing_ratio * pres / (mixing_ratio + 0.622)
    rh = vapour_pressure / max_vapour_pressure * 100.0
    return rh 
    
def __abs_to_spec(humd, temp, pres, max_vapour_pressure):
    global R_S
    vapour_pressure = humd * R_S * temp 
    mixing_ratio = 0.622 * vapour_pressure/ (pres - vapour_pressure)
    q = mixing_ratio / (1 + mixing_ratio)
    return q
    
def __abs_to_rel(humd, temp, pres, max_vapour_pressure):
    
    vapour_pressure = humd * R_S * temp 
    rh = 100.0*vapour_pressure/max_vapour_pressure
    return rh
    
def __rel_to_spec(humd, temp, pres, max_vapour_pressure):
    vapour_pressure = humd * max_vapour_pressure / 100.0 
    mixing_ratio = 0.622 * vapour_pressure/ (pres - vapour_pressure)
    q = mixing_ratio / (1 + mixing_ratio)
    return q
    
def __rel_to_abs(humd, temp, pres, max_vapour_pressure):
    global R_S
    vapour_pressure = humd * max_vapour_pressure / 100.0 
    ah = vapour_pressure / (R_S * temp)
    return ah
    
if __name__ == '__main__':
    verbose = False
    if len(sys.argv) < 7:
        help()
        sys.exit(-1)
    if len(sys.argv) > 7:
        if sys.argv[7] == "-v":
            verbose = True
    #convert([sys.argv[1], float(sys.argv[2])], [sys.argv[3], float(sys.argv[4])], [sys.argv[5], float(sys.argv[6])])
    x = convert(["-s", np.array([3.6, 3.7])], ["-tk", np.array([273.0, 268.0])], ["-p", np.array([1013.0, 1000.0])], verbose)
    print(x)
