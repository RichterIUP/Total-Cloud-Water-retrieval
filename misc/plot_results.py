#!/usr/bin/python

import sys
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
import scipy.optimize as opt
import subprocess
import os
from cloud_alts import *
import average

id_ftir = {'lat': 1, 'lon': 2, 't_l': 3, 't_i': 5, 'rliq': 7, 'rice': 9, 'lwp': 11, 'iwp': 13, 'twp': 15, 'num' : -1}
id_cnet = {"lwp" : 1, "iwp": 3, "rliq": 5, "rice": 7, 'cb' : 9, 'ct': 10, 'num' : 11, 'twp': -2}
id_phase = {'WATER_CLOUD' : 1, 'MIXED_PHASE' : 0, 'ICE_CLOUD': 2} 

def func(x, a, b):
  return a * x + b

def read_csv(csv_file, delim=";", cnet=False, num_of_clouds=None):
  results_raw = np.loadtxt(csv_file, delimiter=delim, dtype=str, skiprows=1)
  f = open(csv_file, "r")
  firstline = f.readline()
  f.close()
  tags = firstline.split(delim)
  number = len(tags)

  if(num_of_clouds == None):
    cells = [[] for ii in range(number)]
  else:
    cells = [[] for ii in range(number+1)]
    
  date = []
  lwp  = []
  lwp_err = []
  iwp  = []
  iwp_err = []
  rliq = []
  rliq_err = []
  rice = []
  rice_err = []
  cb   = []
  ct   = []

  for element in results_raw:
    #print(element)
    #print(element)
    #if(not cnet):
    try:
      cells[0].append(dt.datetime.strptime(element[0], "%Y-%m-%dT%H:%M:%S"))
    except Exception:
      cells[0].append(dt.datetime.strptime(element[0], "%Y%m%d %H%M"))
    for ii in range(1,number):
      #print(ii)
      try:
        cells[ii].append(float(element[ii]))
      except ValueError:
        cells[ii].append(id_phase[element[ii].lstrip().rstrip()])
    #if(num_of_clouds != None):
    #  for kk in range(len(num_of_clouds['datetime'])):
    #    if(num_of_clouds['datetime'][kk] == cells[0][-1]):
    #      cells[number].append(num_of_clouds['num'][kk])
    #      break
    #    if(cells[0][-1] >= dt.datetime(2017, 7, 20)):
    #      cells[number].append(-1)
    #      break
    #  #exit(-1)
  for ii in range(len(cells)):
    cells[ii] = np.array(cells[ii])

  return [tags, cells]

def plot(dataset1, dataset2, quant1, date_limits, ylim, l1="FTIR", l2="Cnet", ftir=False, cnet=True, plot=False, scatter=True, \
         t_l=[0.0, 6.0], t_i=[0.0, 6.0], t_t=[0.0, 6.0], rl=[0.0, 50.0], ri=[0.0, 100.0], \
         twp=[0.0, 100.0], lwp=[0.0, 100.0], iwp=[0.0, 100.0], cb_min=0.0, cb_max=10000.0, \
         geom_thickness=10000.0, icefrac=[0.0, 1.0], factor=None, cnum=1, rr_ratio = 1.0, \
         plotname="plot.svg", dataset3=None, show=True, num_of_clouds=None):

  set1 = {'DATE': [], 'TWP': [], 'TWP_ERR': [], 'IWP': [], 'IWP_ERR': [], 'LWP': [], 'LWP_ERR': [], 'RLIQ': [], 'RLIQ_ERR': [], \
          'RICE': [], 'RICE_ERR': [], 'T_L': [], 'T_L_ERR': [], 'T_I': [], 'T_I_ERR' : [], 'T_T' : [], 'T_T_ERR' : [], \
          'CB': [], 'CT': [], 'TEMP': [], 'PHASE' : [], 'NUM': [], 'DATE_STR': []}
  set2 = {'DATE': [], 'TWP': [], 'TWP_ERR': [], 'IWP': [], 'IWP_ERR': [], 'LWP': [], 'LWP_ERR': [], 'RLIQ': [], 'RLIQ_ERR': [], \
          'RICE': [], 'RICE_ERR': [], 'T_L': [], 'T_L_ERR': [], 'T_I': [], 'T_I_ERR' : [], 'T_T' : [], 'T_T_ERR': [], \
          'CB': [], 'CT': [], 'TEMP': [], 'PHASE': [], 'NUM' : []}
  cnet_set3 = {'DATE': [], 'TWP': [], 'TWP_ERR': []}

  if(dataset2 == None):
    len_dataset2 = len(dataset1[0])
  else:
    len_dataset2 = len(dataset2[0])
  print(len(dataset1[0]))
  print(len_dataset2)
  #exit(-1)
  for ii in range(len(dataset1[0])):
    for jj in range(len_dataset2):
      cond = dataset1[0][ii] >= date_limits[0] and dataset1[0][ii] <= date_limits[1]
      if(dataset2 != None):
        cond = cond and dataset1[0][ii] == dataset2[0][jj]
      if(quant1 == 'twp' and cnet and dataset2 != None):
        cond = cond and dataset2[id_cnet['lwp']][jj] > -100 and dataset2[id_cnet['iwp']][jj] > -100 and dataset2[id_cnet['twp']][jj] >= twp[0]
      elif(quant1 == 'lwp' and cnet):
        cond = cond and dataset2[id_cnet['lwp']][jj] > -100
      elif(quant1 == 'iwp' and cnet):
        cond = cond and dataset2[id_cnet['iwp']][jj] > -100
      elif(quant1 == 'rliq' and cnet):
        cond = cond and dataset2[id_cnet['rliq']][jj] > 0
      elif(quant1 == 'rice' and cnet):
        cond = cond and dataset2[id_cnet['rice']][jj] > 0
      twp1 = dataset1[id_ftir['lwp']][ii]   + dataset1[id_ftir['iwp']][ii]
      dtwp1= dataset1[id_ftir['lwp']+1][ii] + dataset1[id_ftir['iwp']+1][ii]
      twp2 = dataset2[id_ftir['lwp']][jj]   + dataset2[id_ftir['iwp']][jj]
      dtwp2= dataset2[id_ftir['lwp']+1][jj] + dataset2[id_ftir['iwp']+1][jj]
      cond = cond and ((twp1 - dtwp1 <= twp2) and (twp1 + dtwp1 >= twp2) or (twp2 - dtwp2 <= twp1) and (twp2 + dtwp2 >= twp1))# or twp1 <= 5.0 or twp2 <= 5.0)
      #cond = cond and (dataset2[id_cnet['lwp']+1][ii]+dataset2[id_cnet['iwp']+1][ii] < 20.0)
      #cond = cond and (dataset1[0][ii] < dt.datetime(2017, 6, 30, 16) or dataset1[0][ii] >= dt.datetime(2017, 6, 30, 17))
      if(cond):
        set1['DATE'].append(dataset1[0][ii])
        set1['DATE_STR'].append(dt.datetime.strftime(dataset1[0][ii], "%m%d"))
        set1['TWP'].append(dataset1[id_ftir['lwp']][ii]+dataset1[id_ftir['iwp']][ii])
        set1['TWP_ERR'].append(dataset1[id_ftir['lwp']+1][ii]+dataset1[id_ftir['iwp']+1][ii])
        set1['IWP'].append(dataset1[id_ftir['iwp']][ii])
        set1['IWP_ERR'].append(dataset1[id_ftir['iwp']+1][ii])
        set1['LWP'].append(dataset1[id_ftir['lwp']][ii])
        set1['LWP_ERR'].append(dataset1[id_ftir['lwp']+1][ii])
        set1['RLIQ'].append(dataset1[id_ftir['rliq']][ii])
        set1['RLIQ_ERR'].append(dataset1[id_ftir['rliq']+1][ii])
        set1['RICE'].append(dataset1[id_ftir['rice']][ii])
        set1['RICE_ERR'].append(dataset1[id_ftir['rice']+1][ii])
        set1['T_L'].append(dataset1[id_ftir['t_l']][ii])
        set1['T_L_ERR'].append(dataset1[id_ftir['t_l']+1][ii])
        set1['T_I'].append(dataset1[id_ftir['t_i']][ii])
        set1['T_I_ERR'].append(dataset1[id_ftir['t_i']+1][ii])
        set1['T_T'].append(dataset1[id_ftir['t_l']][ii]+dataset1[id_ftir['t_i']][ii])
        set1['T_T_ERR'].append(dataset1[id_ftir['t_l']+1][ii]+dataset1[id_ftir['t_i']+1][ii])
        #set1['CB'].append(dataset1[id_ftir['cb']][ii])
        if(cnet and not ftir and dataset2 != None):
          set2['DATE'].append(dataset2[0][jj])
          set2['TWP'].append(dataset2[id_cnet['twp']][jj])
          set2['TWP_ERR'].append(dataset2[id_cnet['twp']+1][jj])
          set2['IWP'].append(dataset2[id_cnet['iwp']][jj])
          set2['IWP_ERR'].append(dataset2[id_cnet['iwp']+1][jj])
          set2['LWP'].append(dataset2[id_cnet['lwp']][jj])
          set2['LWP_ERR'].append(dataset2[id_cnet['lwp']+1][jj])
          set2['RLIQ'].append(dataset2[id_cnet['rliq']][jj])
          set2['RLIQ_ERR'].append(dataset2[id_cnet['rliq']+1][jj])
          set2['RICE'].append(dataset2[id_cnet['rice']][jj])
          set2['RICE_ERR'].append(dataset2[id_cnet['rice']+1][jj])
          set2['CB'].append(dataset2[id_cnet['cb']][jj])
          set2['CT'].append(dataset2[id_cnet['ct']][jj])
        elif(ftir and not cnet and dataset2 != None):
          set2['DATE'].append(dataset2[0][jj])
          set2['TWP'].append(dataset2[id_ftir['lwp']][jj]+dataset2[id_ftir['iwp']][jj])
          set2['TWP_ERR'].append(dataset2[id_ftir['lwp']+1][jj]+dataset2[id_ftir['iwp']+1][jj])
          set2['IWP'].append(dataset2[id_ftir['iwp']][jj])
          set2['IWP_ERR'].append(dataset2[id_ftir['iwp']+1][jj])
          set2['LWP'].append(dataset2[id_ftir['lwp']][jj])
          set2['LWP_ERR'].append(dataset2[id_ftir['lwp']+1][jj])
          set2['RLIQ'].append(dataset2[id_ftir['rliq']][jj])
          set2['RLIQ_ERR'].append(dataset2[id_ftir['rliq']+1][jj])
          set2['RICE'].append(dataset2[id_ftir['rice']][jj])
          set2['RICE_ERR'].append(dataset2[id_ftir['rice']+1][jj])
          set2['T_L'].append(dataset2[id_ftir['t_l']][jj])
          set2['T_L_ERR'].append(dataset2[id_ftir['t_l']+1][jj])
          set2['T_I'].append(dataset2[id_ftir['t_i']][jj])
          set2['T_I_ERR'].append(dataset2[id_ftir['t_i']+1][jj])
          set2['T_T'].append(dataset2[id_ftir['t_l']][jj]+dataset2[id_ftir['t_i']][jj])
          set2['T_T_ERR'].append(dataset2[id_ftir['t_l']+1][jj]+dataset2[id_ftir['t_i']+1][jj])
          #try:
          #  set2['NUM'].append(dataset2[id_ftir['num']][jj])
          #except IndexError:
          #  set2['NUM'].append(-1)
          break
      else:
        pass

  if(dataset3 != None):
    for ii in range(len(dataset3[0])):
      if(dataset3[id_cnet['lwp']][ii] >= 0.0 and dataset3[id_cnet['iwp']][ii] >= 0.0):
        cnet_set3['DATE'].append(dataset3[0][ii])
        cnet_set3['TWP'].append(dataset3[-2][ii])
        cnet_set3['TWP_ERR'].append(dataset3[-1][ii])
  print(len(set2['TWP']))
  print(len(set1['TWP']))

  q = quant1.upper()
  p = "{}_err".format(q)
  func1 = lambda x, b: x * b

  ftir_set = {'DATE': None, 'TWP' : None, 'IWP' : None, 'LWP': None, 'RLIQ' : None, 'RICE' : None, 'T_L' : None, 'T_I': None, \
              'TWP_err' : None, 'IWP_err' : None, 'LWP_err': None, 'RLIQ_err' : None, 'RICE_err' : None, 'T_L_err' : None, \
              'T_I_err': None, 'T_T': None, 'T_T_err': None, 'DATE_STR': None}

  if(ftir and not cnet):
    cnet_set = {'DATE': None, 'TWP' : None, 'IWP' : None, 'LWP': None, 'RLIQ' : None, 'RICE' : None, 'T_L' : None, 'T_I': None, \
                'TWP_err' : None, 'IWP_err' : None, 'LWP_err': None, 'RLIQ_err' : None, 'RICE_err' : None, 'T_L_err' : None, \
                'T_I_err': None, 'T_T': None, 'T_T_err': None}
  elif(cnet and not ftir):
    cnet_set = {'DATE': None, 'TWP' : None, 'IWP' : None, 'LWP': None, 'RLIQ' : None, 'RICE' : None,  \
                'TWP_err' : None, 'IWP_err' : None, 'LWP_err': None, 'RLIQ_err' : None, 'RICE_err' : None, \
                'CB': None, 'CT': None}
                                             
  for ii in ftir_set.keys():
    if(ftir and dataset2 != None):
      ftir_set[ii] = np.array(set1[ii.upper()])[(np.array(set1['T_L']) <= t_l[1]) & \
                                                (np.array(set1['T_I']) <= t_i[1]) & \
                                                (np.array(set1['T_L']) >= t_l[0]) & \
                                                (np.array(set1['T_I']) >= t_i[0]) & \
                                                (np.array(set1['T_T']) <= t_t[1]) & \
                                                (np.array(set1['T_T']) >= t_t[0]) & \
                                                (np.array(set1['RLIQ']) <= rl[1]) & \
                                                (np.array(set1['RICE']) <= ri[1]) & \
                                                (np.array(set1['RLIQ']) >= rl[0]) & \
                                                (np.array(set1['RICE']) >= ri[0]) & \
                                                (np.array(set1['IWP']) <= iwp[1]) & \
                                                (np.array(set1['LWP']) <= lwp[1]) & \
                                                (np.array(set1['TWP']) <= twp[1]) & \
                                                (np.array(set1['IWP']) >= iwp[0]) & \
                                                (np.array(set1['LWP']) >= lwp[0]) & \
                                                (np.array(set1['TWP']) >= twp[0])]
    elif(ftir and dataset2 == None):
      ftir_set[ii] = np.array(set1[ii.upper()])[(np.array(set1['T_L']) <= t_l[1]) & \
                                                (np.array(set1['T_I']) <= t_i[1]) & \
                                                (np.array(set1['T_L']) >= t_l[0]) & \
                                                (np.array(set1['T_I']) >= t_i[0]) & \
                                                (np.array(set1['T_T']) <= t_t[1]) & \
                                                (np.array(set1['T_T']) >= t_t[0]) & \
                                                (np.array(set1['RLIQ']) <= rl[1]) & \
                                                (np.array(set1['RICE']) <= ri[1]) & \
                                                (np.array(set1['RLIQ']) >= rl[0]) & \
                                                (np.array(set1['RICE']) >= ri[0]) & \
                                                (np.array(set1['IWP']) <= iwp[1]) & \
                                                (np.array(set1['LWP']) <= lwp[1]) & \
                                                (np.array(set1['TWP']) <= twp[1]) & \
                                                (np.array(set1['IWP']) >= iwp[0]) & \
                                                (np.array(set1['LWP']) >= lwp[0]) & \
                                                (np.array(set1['TWP']) >= twp[0])]
    
    else:
      ftir_set[ii] = np.array(set1[ii.upper()])[(np.array(set1['T_L']) <= t_l[1]) & \
                                                (np.array(set1['T_I']) <= t_i[1]) & \
                                                (np.array(set1['T_L']) >= t_l[0]) & \
                                                (np.array(set1['T_I']) >= t_i[0]) & \
                                                (np.array(set1['T_T']) <= t_t[1]) & \
                                                (np.array(set1['T_T']) >= t_t[0]) & \
                                                (np.array(set1['RLIQ']) <= rl[1]) & \
                                                (np.array(set1['RICE']) <= ri[1]) & \
                                                (np.array(set1['RLIQ']) >= rl[0]) & \
                                                (np.array(set1['RICE']) >= ri[0]) & \
                                                (np.array(set1['IWP']) <= iwp[1]) & \
                                                (np.array(set1['LWP']) <= lwp[1]) & \
                                                (np.array(set1['TWP']) <= twp[1]) & \
                                                (np.array(set1['IWP']) >= iwp[0]) & \
                                                (np.array(set1['LWP']) >= lwp[0]) & \
                                                (np.array(set1['TWP']) >= twp[0]) & \
                                                (np.array(set2['CT'])-np.array(set2['CB']) <= geom_thickness)]
  if(cnet or ftir and dataset2 != None):
    for ii in cnet_set.keys():
      if(ftir):
        cnet_set[ii] = np.array(set2[ii.upper()])[(np.array(set1['T_L']) <= t_l[1]) & \
                                                  (np.array(set1['T_I']) <= t_i[1]) & \
                                                  (np.array(set1['T_L']) >= t_l[0]) & \
                                                  (np.array(set1['T_I']) >= t_i[0]) & \
                                                  (np.array(set1['T_T']) <= t_t[1]) & \
                                                  (np.array(set1['T_T']) >= t_t[0]) & \
                                                  (np.array(set1['RLIQ']) <= rl[1]) & \
                                                  (np.array(set1['RICE']) <= ri[1]) & \
                                                  (np.array(set1['RLIQ']) >= rl[0]) & \
                                                  (np.array(set1['RICE']) >= ri[0]) & \
                                                  (np.array(set1['IWP']) <= iwp[1]) & \
                                                  (np.array(set1['LWP']) <= lwp[1]) & \
                                                  (np.array(set1['TWP']) <= twp[1]) & \
                                                  (np.array(set1['IWP']) >= iwp[0]) & \
                                                  (np.array(set1['LWP']) >= lwp[0]) & \
                                                  (np.array(set1['TWP']) >= twp[0])]
      else:
        cnet_set[ii] = np.array(set2[ii.upper()])[(np.array(set1['T_L']) <= t_l[1]) & \
                                                  (np.array(set1['T_I']) <= t_i[1]) & \
                                                  (np.array(set1['T_L']) >= t_l[0]) & \
                                                  (np.array(set1['T_I']) >= t_i[0]) & \
                                                  (np.array(set1['T_T']) <= t_t[1]) & \
                                                  (np.array(set1['T_T']) >= t_t[0]) & \
                                                  (np.array(set1['RLIQ']) <= rl[1]) & \
                                                  (np.array(set1['RICE']) <= ri[1]) & \
                                                  (np.array(set1['RLIQ']) >= rl[0]) & \
                                                  (np.array(set1['RICE']) >= ri[0]) & \
                                                  (np.array(set1['IWP']) <= iwp[1]) & \
                                                  (np.array(set1['LWP']) <= lwp[1]) & \
                                                  (np.array(set1['TWP']) <= twp[1]) & \
                                                  (np.array(set1['IWP']) >= iwp[0]) & \
                                                  (np.array(set1['LWP']) >= lwp[0]) & \
                                                  (np.array(set1['TWP']) >= twp[0]) & \
                                                  (np.array(set2['CT'])-np.array(set2['CB']) <= geom_thickness)]

    
    if(factor == None):
      #print(len(cnet_set[q]), len(ftir_set[q]))
      fit = opt.curve_fit(func1, cnet_set[q], ftir_set[q])
      fact = fit[0][0]**(-1)
    else:
      fact = factor
    #print(fact)
  else:
    fact = 1.0

  print(ftir_set['DATE'][0], cnet_set['DATE'][0], ftir_set['TWP'][0], cnet_set['TWP'][0])
  #exit(-1)
  #[(np.array(set1['IWP'])/np.array(set1['TWP']) >= 0.0) & (np.array(set1['IWP'])/np.array(set1['TWP']) <= 0.2)]

  #print(len(set1['TWP']))
  #print(len(yax))
  if(plot):
    plt.figure(figsize=(20, 10))
    labl1 =r"EM-FTIR"
    labl2 =r"Cloudnet"
    idx = id_cnet['twp']
    #, markeredgewidth=8, markersize=8, elinewidth=5
    plt.errorbar(ftir_set['DATE'], ftir_set[q], yerr=ftir_set[p], xerr=dt.timedelta(seconds=30), fmt="^", markeredgewidth=4, markersize=4, elinewidth=2, color='green', label=labl1)
    if(dataset2 != None):
      plt.errorbar(cnet_set['DATE'], cnet_set[q], yerr=cnet_set[p], xerr=dt.timedelta(seconds=30), fmt="^", markeredgewidth=4, markersize=4, elinewidth=2, color='black', label=labl1)
    #plt.plot(ftir_set['DATE'], ftir_set['NUM'], "x", color="blue")
    if(dataset3 != None):
      cnet_set3['DATE'] = np.array(cnet_set3['DATE'])
      cnet_set3['TWP']  = np.array(cnet_set3['TWP'])
      cnet_set3['TWP_ERR'] = np.array(cnet_set3['TWP_ERR'])
      plt.plot(cnet_set3['DATE'], cnet_set3['TWP'],  color='red', label=labl2)
      #plt.errorbar(cnet_set3['DATE'], cnet_set3['TWP'], yerr=cnet_set3['TWP_ERR'], fmt="v", markeredgewidth=4, markersize=4, elinewidth=2, color='red', label=labl2)
      plt.fill_between(cnet_set3['DATE'], cnet_set3['TWP']-cnet_set3['TWP_ERR'], cnet_set3['TWP']+cnet_set3['TWP_ERR'], color="#bebebe")

    #plt.legend(fontsize=25)
    plt.tick_params(labelsize=20)
    plt.ylabel(r"%s" % (l1), fontsize=30)
    plt.xlabel(r"Datetime", fontsize=30)
    plt.xlim(date_limits)
    plt.ylim(ylim)
    plt.grid(True)
    #plt.show()
    #plt.close()
  
  elif(scatter):
    plt.figure(figsize=(20,10))
    cor = stats.pearsonr(cnet_set[q], ftir_set[q])
    #print(fit[0][0]**(-1), np.float64(np.sqrt(fit[1][0])/fit[0][0]**2))
    print("FACTOR = CNET / FTIR -> {}".format(fit[0][0]**(-1)))
    print("|R| = {}".format(np.abs(cor[0])))
    print("p = {}".format(np.abs(cor[1])))
    print("Number of points = {}".format(len(ftir_set[q])))
    #plt.plot(ylim, func1(np.array(ylim),fit[0][0]), label=r"y = $%f \cdot x$" % fit[0][0])
    plt.plot(ylim, ylim, label=r"$y = x$")
    
    plt.errorbar(cnet_set[q], ftir_set[q], yerr=ftir_set[p], xerr=cnet_set[p], fmt=".", color='black', markeredgewidth=0, markersize=0, elinewidth=1, label=r"$|r| = %f$, $p = %e$" % (np.abs(cor[0]), np.abs(cor[1])))
    #plt.scatter(cnet_set[q], ftir_set[q], s=75, zorder=20)

    sc = plt.scatter(cnet_set[q], ftir_set[q])#, c=cnet_set['CT']-cnet_set['CB'],cmap=plt.get_cmap('viridis'), marker='s', s=75, zorder=20)
    
    #try:
    #  plt.colorbar(sc)
    #except NameError:
    #  pass
    plt.xlabel(l2, fontsize=30)
    plt.xlim(ylim)
    plt.ylim(ylim)
    plt.ylabel(l1, fontsize=30)
    plt.legend(loc=2, fontsize=15)
    plt.tick_params(labelsize=25)
    plt.grid(True)
  if(plotname != None):
    plt.savefig(plotname, dpi=300, bbox_inches='tight')
  else:
    if(show):
      plt.show()
  plt.close()
  f = open("/home/philippr/Arbeit/RETRIEVAL/playground/dates.txt", "w")
  for ii in range(len(ftir_set['DATE'])):
    if(ftir_set['DATE'][ii] == cnet_set['DATE'][ii]):
      f.write("{}\n".format(dt.datetime.strftime(ftir_set['DATE'][ii], "%Y-%m-%dT%H:%M:%S")))
  f.close()

if __name__ == '__main__':
  PS106_Cloudnet      = read_csv("cloudnet_data.csv", cnet=True, delim=",")[1]
  FTIR_conv2          = read_csv(csv_file="/home/philippr/Arbeit/RETRIEVAL/playground/results_DSHIP2.csv", delim=",")[1]
  FTIR_conv           = read_csv(csv_file="/home/philippr/Arbeit/RETRIEVAL/playground/results_DSHIP_cor.csv", delim=",")[1]
  FTIR_pangaea        = read_csv(csv_file="/home/philippr/Arbeit/RETRIEVAL/RESULTS/results_to_copy/microphysical_cloud_parameters_cloud_base_cl51.csv", delim=",")[1]
  
  PS106_Cloudnet.append(PS106_Cloudnet[id_cnet['lwp']]   + PS106_Cloudnet[id_cnet['iwp']])
  PS106_Cloudnet.append(PS106_Cloudnet[id_cnet['lwp']+1] + PS106_Cloudnet[id_cnet['iwp']+1])
  
  min_wp = 0.0
  max_wp =100.0
  iwp_ = [min_wp, max_wp]
  lwp_ = [min_wp, max_wp]
  twp_ = [min_wp, max_wp]
  plot(FTIR_conv, FTIR_pangaea, 'twp', [dt.datetime(2017, 6, 10), dt.datetime(2017, 8, 12)], [-10.0, 100.0], \
       ftir=True, cnet=False, geom_thickness=10000.0, cb_max=10000.0, cnum=10, rr_ratio=50.0, \
       t_l = [0.0, 4.0], t_i = [0.0, 4.0], t_t=[0.0, 4.0], rl = [0.0, 100.0], ri = [0.0, 100.0], \
       iwp = iwp_, lwp = lwp_, twp = twp_, plot=False, scatter=True, \
       l1 = r"FTIR TWP", l2=r"CNET TWP", \
       plotname=None, dataset3=None, show=True)  
