# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:23:11 2020

@author: Philipp
"""

import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

# CLARRA
Rm = np.loadtxt("C:/Users/Philipp/Rm.csv", delimiter=",")
Ra = np.loadtxt("C:/Users/Philipp/Rm.csv", delimiter=",")
Ra = np.loadtxt("C:/Users/Philipp/Ra.csv", delimiter=",")
nu = np.loadtxt("C:/Users/Philipp/nu.csv", delimiter=",")

# TCWret

path_tcwret = "C:/Users/Philipp/RESULTS_TCWret_28_02_2020"
tcwret = "{}/results_PS.20170701_144500.nc.nc".format(path_tcwret)
        

f = nc.Dataset("{}".format(tcwret))
wn_tcwret = f.variables['wavenumber'][:]
ra_tcwret = f.variables['lbldis radiance'][:]
ms_tcwret = f.variables['ftir radiance'][:]
mcp_tcwret = [np.round(np.float_(f.variables['tt'][:]), 2), \
              np.round(np.float_(f.variables['fi'][:]), 2), \
              np.round(np.float_(f.variables['rl'][:]), 2), \
              np.round(np.float_(f.variables['ri'][:]), 2)]

# Measurement
wn =  np.loadtxt("C:/Users/Philipp/wavenumber.csv", delimiter=",")
Rftir = np.loadtxt("C:/Users/Philipp/radiance_ftir.csv", delimiter=",")

plt.plot(wn_tcwret, ms_tcwret, ".", label="TCWret")
plt.plot(wn, Rftir, ".", label="FTIR")
plt.legend()
plt.xlim([770, 1200])
plt.ylim([0, 120])
plt.grid(True)
plt.show()