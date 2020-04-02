#!/usr/bin/python

'''
Read the retrieval results. Write them to file or plot them
'''

import os
import datetime as dt
import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np

RESULTS = {"DATETIME": [], "LAT": [], "LON": [], \
           "OD_LIQ": [], "OD_LIQ_ERR": [], "OD_ICE": [], "OD_ICE_ERR": [], \
           "OD_TOTAL": [], "OD_TOTAL_ERR": [], \
           "REFF_LIQ": [], "REFF_LIQ_ERR": [], "REFF_ICE": [], "REFF_ICE_ERR": [], \
           "LWP": [], "LWP_ERR": [], "IWP": [], "IWP_ERR": [], \
           "TWP": [], "TWP_ERR": [], \
           "AVK_[0,0]": [], "AVK_[0,1]" : [], "AVK_[0,2]": [], "AVK_[0,3]": [], \
           "AVK_[1,0]": [], "AVK_[1,1]" : [], "AVK_[1,2]": [], "AVK_[1,3]": [], \
           "AVK_[2,0]": [], "AVK_[2,1]" : [], "AVK_[2,2]": [], "AVK_[2,3]": [], \
           "AVK_[3,0]": [], "AVK_[3,1]" : [], "AVK_[3,2]": [], "AVK_[3,3]": [], \
           'CHI_2': [], "RES": [], "NOISE": [], "PWV": [], "Fname": [], "TotalPath": []}

KEYS = ["LAT", "LON", "OD_LIQ", "OD_LIQ_ERR", "OD_ICE", "OD_ICE_ERR", \
        "OD_TOTAL", "OD_TOTAL_ERR", "REFF_LIQ", "REFF_LIQ_ERR", "REFF_ICE", "REFF_ICE_ERR", \
        "LWP", "LWP_ERR", "IWP", "IWP_ERR", "TWP", "TWP_ERR", "AVK_[0,0]", "AVK_[0,1]", \
        "AVK_[0,2]", "AVK_[0,3]", "AVK_[1,0]", "AVK_[1,1]", "AVK_[1,2]", "AVK_[1,3]", \
        "AVK_[2,0]", "AVK_[2,1]", "AVK_[2,2]", "AVK_[2,3]", \
        "AVK_[3,0]", "AVK_[3,1]", "AVK_[3,2]", "AVK_[3,3]", "CHI_2", \
        "RES", "NOISE", "PWV", "Fname"]

def __go_through_subfolders(path):
    '''
    Recursively traverse through subfolders
    '''
    for element in sorted(os.listdir(path)):
        file_or_folder = "{}/{}".format(path, element)
        if os.path.isdir(file_or_folder):
            return __go_through_subfolders(file_or_folder)
    return [path, os.listdir(path)]

def __get_results(path):
    '''
    Read result_iteration
    Line 2: Datetime
    Line 3: Latitude
    Line 4: Longitude
    Line 7: OD Liquid
    Line 8: OD Ice
    Line 9: Reff Liquid
    Line 10: Reff Ice
    Line 11: LWP
    Line 12: IWP
    Line 17-20: Averaging kernel matrix
    Line 24: Chi_2
    '''
    
    fname = path.split("/")[-2]

    res_file = "{}/result_iteration".format(path)
    file_ = open(res_file, "r")
    cont = file_.readlines()
    file_.close()

    datestr = cont[2].split("Datetime:")[-1].lstrip().rstrip()

    datetime = dt.datetime.strptime(datestr, "%m/%d/%Y %H:%M")
    lat = float(cont[3].split(":")[-1].lstrip().rstrip())
    lon = float(cont[4].split(":")[-1].lstrip().rstrip())
    
    tau_liq = float(cont[7].split("(")[2].split("+")[0])
    tau_liq_err = float(cont[7].split(")")[1].split("-")[-1])
    od_liq = [tau_liq, tau_liq_err]
    
    tau_ice = float(cont[8].split("(")[2].split("+")[0])
    tau_ice_err = float(cont[8].split(")")[1].split("-")[-1])
    od_ice = [tau_ice, tau_ice_err]
    
    reff_liq = float(cont[9].split("(")[2].split("+")[0])
    reff_liq_err = float(cont[9].split(")")[1].split("-")[-1])
    radius_liq = [reff_liq, reff_liq_err]
    
    reff_ice= float(cont[10].split("(")[2].split("+")[0])
    reff_ice_err = float(cont[10].split(")")[1].split("-")[-1])
    radius_ice = [reff_ice, reff_ice_err]
    
    liquid_wp = float(cont[11].split("(")[2].split("+")[0])
    liquid_wp_err = float(cont[11].split(")")[1].split("-")[-1])
    lwp = [liquid_wp, liquid_wp_err]
        
    ice_wp = float(cont[12].split("(")[2].split("+")[0])
    ice_wp_err = float(cont[12].split(")")[1].split("-")[-1])
    iwp = [ice_wp, ice_wp_err]
    
    avk = []
    for in_line in range(4):
        for column in range(4):
            avk.append(float(cont[column+17].split("\t")[in_line]))
    
    chi2 = float(cont[24].split("=")[-1])
    try:
        noise = float(cont[25].split("=")[-1])
        res = float(cont[26].split("=")[-1])
        pwv = float(cont[27].split("=")[-1])
    except IndexError:
        ret_log_file = open("{}/retrieval_log.dat".format(path), "r")
        cont = ret_log_file.readlines()
        ret_log_file.close()
        for ii in range(len(cont)):
            if("Mean temperature" in cont[ii]):
                index_noise = ii + 4
                index_pwv = ii + 5
                index_iter = ii + 8
                break
        residuals = []
        for ii in range(index_iter, len(cont)):
            if("Finished" in cont[ii]):
                break
            residuals.append(float(cont[ii].split("\t")[-2]))
        res = min(residuals)
        noise = float(cont[index_noise].split(":")[-1])
        pwv = float(cont[index_pwv].split(":")[-1]) 
       
    return [datetime, lat, lon, od_liq, od_ice, radius_liq, radius_ice, \
            lwp, iwp, avk, chi2, noise, res, pwv, fname]

def __update_results(to_append, path):
    '''
    Append the new results to the old result dict
    '''
    global RESULTS
    RESULTS['DATETIME'].append(to_append[0])
    RESULTS['LAT'].append(to_append[1])
    RESULTS['LON'].append(to_append[2])
    RESULTS['OD_LIQ'].append(to_append[3][0])
    RESULTS['OD_LIQ_ERR'].append(to_append[3][1])
    RESULTS['OD_ICE'].append(to_append[4][0])
    RESULTS['OD_ICE_ERR'].append(to_append[4][1])
    RESULTS['OD_TOTAL'].append(to_append[3][0] + to_append[4][0])
    RESULTS['OD_TOTAL_ERR'].append(to_append[3][1] + to_append[4][1])
    RESULTS['REFF_LIQ'].append(to_append[5][0])
    RESULTS['REFF_LIQ_ERR'].append(to_append[5][1])
    RESULTS['REFF_ICE'].append(to_append[6][0])
    RESULTS['REFF_ICE_ERR'].append(to_append[6][1])
    RESULTS['LWP'].append(to_append[7][0])
    RESULTS['LWP_ERR'].append(to_append[7][1])
    RESULTS['IWP'].append(to_append[8][0])
    RESULTS['IWP_ERR'].append(to_append[8][1])
    RESULTS['TWP'].append(to_append[7][0] + to_append[8][0])
    RESULTS['TWP_ERR'].append(to_append[7][1] + to_append[8][1])
    for row in range(4):
        for column in range(4):
            RESULTS['AVK_[{},{}]'.format(row, column)].append(to_append[9][row * 4 + column])
    RESULTS['CHI_2'].append(to_append[10])
    RESULTS['NOISE'].append(to_append[11])
    RESULTS['RES'].append(to_append[12])
    RESULTS['PWV'].append(to_append[13])
    RESULTS['Fname'].append(to_append[14])
    RESULTS['TotalPath'].append(path)
    return RESULTS

def read_retrieval_results(path):
    '''
    Read the retrieval results
    '''
    
    all_folders = sorted(os.listdir(path))
    for folder in all_folders:
        res_path = "{}/{}".format(path, folder)
        if os.path.isdir(res_path):
            [path_to_res, files] = __go_through_subfolders(res_path)
            if "result_iteration" in files:
                res_ = __get_results(path_to_res)
                RESULTS = __update_results(res_, path_to_res)
    return RESULTS

def write_to_pangaea(path, fname="results.csv", maxchi=100000.0, maxres=100.0):
    pangaea_file = open(fname, "w")
    pangaea_file.write("Date/Time,Latitude (degree),Longitude (degree),")
    pangaea_file.write("Error optical depth liquid (1),Optical depth liquid (1),")
    pangaea_file.write("Optical depth ice (1), Error optical depth ice (1),")
    pangaea_file.write("Optical depth total (1), Error optical depth total (1),")
    pangaea_file.write("Effective droplet radius liquid (um),Error effective droplet radius liquid (um),")
    pangaea_file.write("Effective droplet radius ice (um),Error effective droplet radius (um),")
    pangaea_file.write("Liquid water path (g m^(-2)),Error liquid water path (g m^(-2)),")
    pangaea_file.write("Ice water path (g m^(-2)),Error ice water path (g m^(-2)),")
    pangaea_file.write("Total water path (g m^(-2)),Error total water path (g m^(-2)),")
    pangaea_file.write("Averaging Kernel matrix [11],Averaging Kernel matrix [12],")
    pangaea_file.write("Averaging Kernel matrix [13],Averaging Kernel matrix [14],")
    pangaea_file.write("Averaging Kernel matrix [21],Averaging Kernel matrix [22],")
    pangaea_file.write("Averaging Kernel matrix [23],Averaging Kernel matrix [24],")
    pangaea_file.write("Averaging Kernel matrix [31],Averaging Kernel matrix [32],")
    pangaea_file.write("Averaging Kernel matrix [33],Averaging Kernel matrix [34],")
    pangaea_file.write("Averaging Kernel matrix [41],Averaging Kernel matrix [42],")
    pangaea_file.write("Averaging Kernel matrix [43],Averaging Kernel matrix [44],")
    pangaea_file.write("Chi_2,Residual,Noise,PWV,Fname\n")
    num_of_datapoints = len(RESULTS['DATETIME'])
    for line in range(num_of_datapoints):
        if RESULTS['CHI_2'][line] < maxchi and RESULTS['RES'][line] < maxres:
            pangaea_file.write("{},".format(dt.datetime.strftime(RESULTS['DATETIME'][line], "%Y-%m-%dT%H:%M:00")))
            for key in KEYS:
                pangaea_file.write("{},".format(RESULTS[key][line]))
            pangaea_file.write("\n")
            plotname = RESULTS['TotalPath'][line].split("/")[-1]
            os.system("cp {}/tl_* {}/plots/full_{}.svg".format(RESULTS['TotalPath'][line], \
                      path, plotname))
            os.system("cp {}/iter_final_* {}/plots/mw_{}.svg".format(RESULTS['TotalPath'][line], \
                      path, plotname))
    pangaea_file.close()
    return 

def read_result_file(fname, delimiter=","):
    res_file = open(fname, "r")
    cont = res_file.readlines()
    res_file.close()
    
    tags = []
    for element in cont[0].split(delimiter):
        tags.append(element.rstrip())
    
    num_of_datapoints = len(cont[1:])
    results = dict()
    for tag in tags:
        results.update({tag.rstrip(): []})
        
    for dp in range(1, num_of_datapoints):
        for loop in range(len(tags)):
            element = cont[dp].split(delimiter)[loop].rstrip()

            if "Date" in tags[loop].rstrip():
                element = dt.datetime.strptime(cont[dp].split(delimiter)[loop].rstrip(), "%Y-%m-%dT%H:%M:%S")
            else:
                try:
                    element = float(cont[dp].split(delimiter)[loop].rstrip())
                except ValueError:
                    pass
            results[tags[loop].rstrip()].append(element)
            
    for loop in range(len(tags)):
        results[tags[loop].rstrip()] = np.array(results[tags[loop].rstrip()])
    return [tags, results]

def __average(value_axis, tags):
    dp = 0
    value_cache = dict()
    for element in tags:
        value_cache.update({element: []})
    #First: Average
    while dp < len(value_axis['Date/Time'])-2:
        cond_1 = value_axis['Date/Time'][dp] == value_axis['Date/Time'][dp+1]
        cond_2 = value_axis['Date/Time'][dp] == value_axis['Date/Time'][dp+2]
        if cond_1 and not cond_2:
            jump = 2
        elif cond_1 and cond_2:
            jump = 3
        else:
            jump = 1
            
        value_cache['Date/Time'].append(value_axis['Date/Time'][dp])
        for loop in range(1, len(value_axis)):
            try:
                value_cache[tags[loop]].append(np.mean(value_axis[tags[loop].rstrip()][dp:dp+jump]))
            except TypeError:
                value_cache[tags[loop]].append(value_axis[tags[loop]][dp])

        dp = dp + jump
        #date_cache.append(date_axis[dp])
    return value_cache

def find_corresponding(value_axis_1, value_axis_2, tags_1, tags_2):
    '''
    Find corresponding datapoints. Average over datapoints with the same timestamp
    '''

    value_corr_1 = dict()
    for element in tags_1:
        value_corr_1.update({element: []})
    value_corr_2 = dict()
    for element in tags_2:
        value_corr_2.update({element: []})

    value_axis_1 = __average(value_axis_1, tags_1)
    value_axis_2 = __average(value_axis_2, tags_2)
    #print(value_axis_1['Date/Time'])
    for d1 in range(len(value_axis_1['Date/Time'])):
        for d2 in range(len(value_axis_2['Date/Time'])):
            if value_axis_1['Date/Time'][d1] == value_axis_2['Date/Time'][d2]: 
                for key in tags_1:
                    value_corr_1[key].append(value_axis_1[key][d1])
                for key in tags_2:
                    value_corr_2[key].append(value_axis_2[key][d2])
                
    for key in tags_1:
        value_corr_1[key] = np.array(value_corr_1[key])
    for key in tags_2:
        value_corr_2[key] = np.array(value_corr_2[key])
    return [value_corr_1, value_corr_2]

def quickplot(xaxis, yaxis, xerr=False, yerr=False, xlim=False, ylim=False, save=False):
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    if xerr and yerr:
        ax.errorbar(xaxis, yaxis, xerr=xerr, yerr=yerr, fmt=".")
    elif xerr and not yerr:
        ax.errorbar(xaxis, yaxis, xerr=xerr, fmt=".")
    elif yerr and not xerr:
        ax.errorbar(xaxis, yaxis, yerr=yerr, fmt=".")
    else:
        ax.plot(xaxis, yaxis, ".")
    if xlim and ylim:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    elif xlim and not ylim:
        ax.set_xlim(xlim)
    elif ylim and not xlim:
        ax.set_ylim(ylim)
    ax.grid(True)
    if save:
        plt.savefig(save, dpi=fig.dpi, bbox_inches="tight")
    else:
        plt.show()
    return
    
if __name__ == '__main__':
    #FOLDER = "/home/philippr/res/OUTFOLDER"
    #RESULTS = read_retrieval_results(FOLDER)
    #write_to_pangaea(FOLDER, fname="/home/philippr/res_pangaea.csv", maxres=0.5, maxchi=2000.0)
    file_015 = "/media/philippr/PS1067Backup/RESULTS/PASCAL_SIPCA_FRAM/microphysical_cloud_parameters_cloud_base_cl51.csv"
    file_05 = "/media/philippr/PS1067Backup/RESULTS/PASCAL_SIPCA_FRAM/microphysical_cloud_parameters_res_05cm.csv"
    file_cnet = "/media/philippr/PS1067Backup/RESULTS/cloudnet_data.csv"
    [tags_015, res_015] = read_result_file(file_015)
    [tags_05, res_05] = read_result_file(file_05)
    [tags_cnet, res_cnet] = read_result_file(file_cnet)
    [res_015, res_05] = find_corresponding(res_015, res_05, tags_015, tags_05)
    idx_1 = 15
    idx_2 = 13
    fi_05 = res_05[tags_05[5]]/res_05[tags_05[7]]
    fi_015 = res_015[tags_015[5]]/(res_015[tags_015[3]]+res_015[tags_015[5]])
    lbl = st.pearsonr(fi_05, fi_015)
    fig = plt.figure(figsize=(20, 10))
    plt.plot(fi_05, fi_015, ".", label=lbl)
    plt.title("f_i (1)")
    plt.xlabel("res=0.5cm")
    plt.ylabel("res=0.15cm")
    plt.legend()
    plt.savefig("corr_fi.png", dpi=100.0, bbox_inches="tight")
    plt.close()

    #quickplot(np.array(RESULTS['LWP'])+np.array(RESULTS['IWP']), RESULTS['TWP'])