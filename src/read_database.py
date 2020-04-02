#!/usr/bin/python

import scipy.interpolate as scp
import numpy             as np
import matplotlib.pyplot as plt

def remove_whitespace(line):
    output = []
    counter = 0
    for i in line.split(" "):
        try:
            output.append(float(i))
            counter = counter + 1
        except ValueError:
            pass
        if counter > 12:
            return output
    return output

def read_databases(database_liq="mie_database/ssp_db.mie_wat.gamma_sigma_0p100", \
                    database_ice="mie_database/ssp_db.mie_ice.gamma_sigma_0p100"):

    db_ice_fd = open(database_ice)
    db_ice = db_ice_fd.readlines()
    db_ice_fd.close()

    db_liq_fd = open(database_liq)
    db_liq = db_liq_fd.readlines()
    db_liq_fd.close()

    entry_liq = []
    entry_ice = []
    set_radius = []
    set_wavenumber = []

    #print("Read database ice ...")
    for j in range(7, len(db_ice)):
        line = remove_whitespace(db_ice[j])
        if line[1] >= 500.0 and line[1] <= 1500.0:
            entry_ice.append({"wavelength (um)" : line[0], "wavenumber (cm-1)" : line[1], \
            "reff (um)" : line[2], "ext (um2)" : line[3], "scat (um2)" : line[4], \
            "abso (um2)": line[5], "ssa" : line[6], "g" : line[7], \
            "Qext" : line[8], "Qabs" : line[9], "Qsca" : line[10], \
            "Vol (um3)" : line[11], "Proj_area (um2)": line[12]})

    #print("Read database liquid ...")
    for j in range(7, len(db_liq)):
        line = remove_whitespace(db_liq[j])
        if line[1] >= 500.0 and line[1] <= 1500.0:
            entry_liq.append({"wavelength (um)" : line[0], "wavenumber (cm-1)" : line[1], \
            "reff (um)" : line[2], "ext (um2)" : line[3], "scat (um2)" : line[4], \
            "abso (um2)": line[5], "ssa" : line[6], "g" : line[7], \
            "Qext" : line[8], "Qabs" : line[9], "Qsca" : line[10], \
            "Vol (um3)" : line[11], "Proj_area (um2)": line[12]})

    for i in range(len(entry_liq)):
        set_wavenumber.append(entry_liq[i]['wavenumber (cm-1)'])
        set_radius.append(entry_liq[i]['reff (um)'])

    return [entry_liq, entry_ice]

def get_entry(database, key, wavenumber, radius):
    radii = []
    wn = []
    entries = []
    x_ax = []
    y_ax = []
    z_ax = []
    sw = []
    sr = []

    if wavenumber <= 0.0 or radius <= 0.0:
        return 0.0
    for i in range(len(database)):
        sw.append(database[i]['wavenumber (cm-1)'])
        sr.append(database[i]['reff (um)'])
    #Finden der Radien und Wellenzahlen, welche direkt neben den gewuenschten Werten liegen
    for i in range(len(sr)-1):
        if radius > sr[i] and radius < sr[i+1]:
            radii = [sr[i], sr[i+1]]
        elif radius == sr[i]:
            radii = [sr[i]]

    for i in range(len(sw)-1):
        if wavenumber > sw[i] and wavenumber < sw[i+1]:
            wn = [sw[i], sw[i+1]]
        elif wavenumber == sw[i]:
            wn = [sw[i]]

    for i in range(len(database)):
        if database[i]['reff (um)'] == radii[0] or database[i]['reff (um)'] == radii[-1]:
            if database[i]['wavenumber (cm-1)'] == wn[0] or \
            database[i]['wavenumber (cm-1)'] == wn[-1]:
                x_ax.append(database[i]['wavenumber (cm-1)'])
                y_ax.append(database[i]['reff (um)'])
                z_ax.append(database[i][key])
                entries.append([database[i]['reff (um)'], \
                                database[i]['wavenumber (cm-1)'], database[i][key]])

    if len(wn) == 2 and len(radii) == 2:
        a = scp.interp2d(x_ax, y_ax, z_ax)
        a = a(wavenumber, radius)
    elif len(wn) == 1 and len(radii) == 2:
        a = scp.interp1d(y_ax, z_ax)
        a = a(radius)
    elif len(wn) == 2 and len(radii) == 1:
        a = scp.interp1d(x_ax, z_ax)
        a = a(wavenumber)
    elif len(wn) == 1 and len(radii) == 1:
        a = z_ax[0]
    return a

def __calc_iwp(Vol, ext, tau):
    if(Vol <= 0.0 or ext <= 0.0):
        return np.array([0.0])
    ext = np.array(ext)
    Vol = np.array(Vol)
    N_0 = 1/Vol
    N = tau/ext
    rho = 916896
    ice_water_path = N / N_0 * rho
    return np.float_(ice_water_path)

def __calc_lwp(radius, tau_ir):
    rho = 1000000
    liquid_water_path = 2.0/3.0 * radius*1e-6 * tau_ir * rho
    return np.float_(liquid_water_path)

def calc_iwp(tau_ice, dtau_ice, reff_ice, dreff_ice, db_ice):
    ret_wn = 900.0
    delta = 1e-3

    Vol = get_entry(db_ice, "Vol (um3)", ret_wn, reff_ice)*1e-18
    delta_Vol_1_r = get_entry(db_ice, "Vol (um3)", ret_wn, reff_ice+delta)
    delta_Vol_2_r = get_entry(db_ice, "Vol (um3)", ret_wn, reff_ice-delta)
    delta_Vol_r = (delta_Vol_1_r - delta_Vol_2_r)/(2*delta) * 1e-18
    delta_Vol = np.abs(delta_Vol_r * dreff_ice) * 1e-6

    ext = get_entry(db_ice, "ext (um2)", ret_wn, reff_ice)*1e-12
    delta_ext_1_r = get_entry(db_ice, "ext (um2)", ret_wn, reff_ice+delta)
    delta_ext_2_r = get_entry(db_ice, "ext (um2)", ret_wn, reff_ice-delta)
    delta_ext_r = (delta_ext_1_r - delta_ext_2_r)/(2*delta) * 1e-12
    delta_ext = np.abs(delta_ext_r * dreff_ice) * 1e-6

    err_1 = np.abs(__calc_iwp(Vol, ext, dtau_ice))
    err_2 = np.abs(__calc_iwp(delta_Vol, ext, tau_ice))
    err_3 = np.abs(__calc_iwp(Vol, ext, tau_ice)/(-ext)*delta_ext)
    ice_water_path = np.float_(__calc_iwp(Vol, ext, tau_ice))
    dice_water_path = np.float_(err_1 + err_2 + err_3)
    return [ice_water_path, dice_water_path]

def calc_lwp(reff_liq, dreff_liq, tau_liq, dtau_liq):
    liquid_water_path = __calc_lwp(reff_liq, tau_liq)
    err_1 = np.abs(__calc_lwp(dreff_liq, tau_liq))
    err_2 = np.abs(__calc_lwp(reff_liq, dtau_liq))
    dliquid_water_path = np.float_(err_1 + err_2)
    return [liquid_water_path, dliquid_water_path]

if __name__ == '__main__':
    [liq, ice] = read_databases("../ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
    "../ssp/ssp_db.mie_ice.gamma_sigma_0p100")
    print(get_entry(db_liq, "Qext", 900.0, ))

    t_i = 1.0
    r_i_arr = np.array([ii*1.0+10.0 for ii in range(60)])
    dtl = 0.0
    dti = 0.0
    drl = 0.0
    dri = 0.0
    #t_l = 3.66 * (1-0.34)
    #t_i = 3.66 * 0.34
    t_t = [8.4] 
    f_i = [0.91]
    r_l = 2.678
    r_i = 17.25
    for jj in range(1):
        lwp = calc_lwp(r_l, drl, t_t[jj]*(1-f_i[jj]), dtl)
        iwp = calc_iwp(t_t[jj]*f_i[jj], dti, r_i, dri, ice)
        print(iwp[0], lwp[0])
    #01.521		 00.418		 11.959		 42.294
    #IWP = []
    #for r_i in np.array([(ii*1.0+10.0) for ii in range(60)]):
    #    IWP.append(calc_iwp(t_i, dt_i, r_i, dr_i, ice)[0])#

    #plt.plot(r_i_arr, IWP)
    #plt.xlabel("r_i")
    #plt.ylabel("IWP")
    #plt.grid(True)
    #plt.show()
