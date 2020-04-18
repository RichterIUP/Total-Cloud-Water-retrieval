#!/usr/bin/python3
'''@package docstring
Write the results to a netCDF file
'''

import os
import sys
sys.path.append("./src")

import scipy.io as sio
import numpy as np
import inp
import aux2 as aux


def create_nc(chi_2, index=-1, avk_matrix=None, errors=None, nc=1, covariance_matrix=None, transfer_matrix=None, , errors_res=None):
    '''
    Create the netCDF file
    
    @param chi2 Cost function
    @param avk_matrix Averaging kernel matrix
    @param errors Standard deviations of MCP
    '''
    if not os.path.exists(inp.RESULTS):
        os.mkdir(inp.RESULTS)
        
    nc_fname = "{}/results_{}.nc".format(inp.RESULTS, aux.FTIR.split("/")[-1])
    with sio.netcdf_file(nc_fname, "w") as outfile:
        outfile.createDimension("const", 1)
        outfile.createDimension("num_of_clouds", len(aux.CLOUD_BASE))
        outfile.createDimension("mcp", 4)
        outfile.createDimension("mcp_err", 2)
        outfile.createDimension("wp", 3)
        outfile.createDimension("level", len(aux.ATMOSPHERIC_GRID[0]))
        outfile.createDimension("wavenumber", len(aux.WAVENUMBER_FTIR))
        
        conv = outfile.createVariable("conv", "i", ("const", ))
        conv.units = "1"
        lat = outfile.createVariable("lat", "f8", ("const", ))
        lat.units = "deg"
        lon = outfile.createVariable("lon", "f8", ("const", ))
        lon.units = "deg"
        ctemp = outfile.createVariable("av_ctemp", "f8", ("const", ))
        ctemp.unit = "K"
        
        specname = outfile.createVariable("specname", "S1", ("const", ))
        specname.units = "Name"
        num = outfile.createVariable("numclouds", "f8", ("const", ))
        num.units = "1"
        cbh = outfile.createVariable("cloud_base", "f8", ("num_of_clouds", ))
        cbh.units = "m"
        cth = outfile.createVariable("cloud_top", "f8", ("num_of_clouds", ))
        cth.units = "m"
        clevel = outfile.createVariable("clevel", "f8", ("const", ))
        clevel.units = "1"
        sza = outfile.createVariable("sza", "f8", ("const", ))
        sza.units = "deg"
    
        pres = outfile.createVariable("P", "f8", ("level", ))
        pres.units = "hPa"
        temp = outfile.createVariable("T", "f8", ("level", ))
        temp.units = "K"
        humd = outfile.createVariable("humidity", "f8", ("level", ))
        humd.units = inp.HUMIDITY
        alt = outfile.createVariable("z", "f8", ("level", ))
        alt.units = "km"
        pwv = outfile.createVariable("pwv", "f8", ("level", ))
        pwv.units = "cm"
    
        wavenumber = outfile.createVariable("wavenumber", "f8", ("wavenumber",))
        wavenumber.units = "cm-1"
        ftir_radiance = outfile.createVariable("ftir radiance", "f8", ("wavenumber", ))
        ftir_radiance.units = "mW * (sr*cm-1*m2)**(-1)"
        lbldis_radiance = outfile.createVariable("lbldis radiance", "f8", ("wavenumber", ))
        lbldis_radiance.units = "mW * (sr*cm-1*m2)**(-1)"
        residuum = outfile.createVariable("residuum", "f8", ("wavenumber", ))
        residuum.units = "mW * (sr*cm-1*m2)**(-1)"
        rms = outfile.createVariable("Root-Mean-Square", "f8", ("const", ))
        rms.units = "mW * (sr*cm-1*m2)**(-1)"
        
        chi2 = outfile.createVariable("cost function", "f8", ("const", ))
        chi2.units = "1"
        
        avk = outfile.createVariable("averaging kernel matrix", "f8", ("mcp", "mcp"))
        avk.units = "1"

        cov_mat = outfile.createVariable("covariance matrix", "f8", ("mcp", "mcp"))
        cov_mat.units = "1"
                
        t_mat = outfile.createVariable("transfer matrix", "f8", ("mcp", "wavenumber"))
        t_mat.units = "1"    
        #tt = outfile.createVariable("tl", "f8", ("const", ))
        #tt.units = "1"
        #fi = outfile.createVariable("ti", "f8", ("const", ))
        #fi.units = "1"
        #rl = outfile.createVariable("rl", "f8", ("const", ))
        #rl.units = "um"
        #ri = outfile.createVariable("ri", "f8", ("const", ))
        #ri.units = "um"
        #iwp = outfile.createVariable("iwp", "f8", ("const", ))
        #iwp.units = "g/m2"
        #lwp = outfile.createVariable("lwp", "f8", ("const", ))
        #lwp.units = "g/m2"
        #twp = outfile.createVariable("twp", "f8", ("const", ))
        #twp.units = "g/m2"
        
        x_ret = outfile.createVariable('x_ret', 'f8', ('mcp', ))
        x_ret.units = '1'
        x_ret_err = outfile.createVariable('x_ret_err', 'f8', ('mcp', ))
        x_ret_err.units = '1'
        x_err_res = outfile.createVariable('x_err_res', 'f8', ('mcp', ))
        x_err_res.units = '1'
        
        wp_ret = outfile.createVariable('wp_ret', 'f8', ('wp', ))
        wp_ret.units = 'g/m2'
        wp_ret_err = outfile.createVariable('wp_ret_err', 'f8', ('wp', ))
        wp_ret_err.units = 'g/m2'
        wp_err_res = outfile.createVariable('wp_err_res', 'f8', ('wp', ))
        wp_err_res.units = 'g/m2'
                
        x_a = outfile.createVariable("x_a", "f8", ("mcp", ))
        x_a.units = "1"
        x_a_err = outfile.createVariable('x_a_err', 'f8', ('mcp', ))
        x_a_err.units = '1'
        
        #dtt = outfile.createVariable("dtl", "f8", ("const", ))
        #dtt.units = "1"
        #dfi = outfile.createVariable("dti", "f8", ("const", ))
        #dfi.units = "1"
        #drl = outfile.createVariable("drl", "f8", ("const", ))
        #drl.units = "um"
        #dri = outfile.createVariable("dri", "f8", ("const", ))
        #dri.units = "um"
        #diwp = outfile.createVariable("diwp", "f8", ("const", ))
        #diwp.units = "g/m2"
        #dlwp = outfile.createVariable("dlwp", "f8", ("const", ))
        #dlwp.units = "g/m2"
        #dtwp = outfile.createVariable("dtwp", "f8", ("const", ))
        #dtwp.units = "g/m2"
    
        '''
        Write data
        '''

        conv[:] = nc
        specname[:] = aux.FTIR.split("/")[-1]
        wavenumber[:] = aux.WAVENUMBER_FTIR[:]
        ftir_radiance[:] = aux.RADIANCE_FTIR[:]
        lbldis_radiance[:] = aux.RADIANCE_LBLDIS[0][-1][:]
        residuum[:] = list(aux.RESIDUUM[index]) 
        rms[:] = np.sqrt(np.mean(np.array(aux.RESIDUUM[index])**2))
        cbh[:] = aux.CLOUD_BASE[:]
        cth[:] = aux.CLOUD_TOP[:]
        num[:] = len(aux.CLOUD_BASE[:])
        clevel[:] = len(aux.CLOUD_LAYERS)
        lat[:] = aux.LAT
        lon[:] = aux.LON
        sza[:] = aux.SOLAR_ZENITH_ANGLE
        ctemp[:] = aux.CLOUD_TEMP
        pwv[:] = aux.PRECIPITABLE_WATER_VAPOUR
        pres[:] = aux.ATMOSPHERIC_GRID[0][:]
        alt[:] = aux.ATMOSPHERIC_GRID[1][:]
        temp[:] = aux.ATMOSPHERIC_GRID[2][:]
        humd[:] = aux.ATMOSPHERIC_GRID[3][:]
        
        chi2[:] = chi_2
        if type(avk_matrix) != type(None):
            avk[:] = avk_matrix[:]
        if type(covariance_matrix) != type(None):
            cov_mat[:] = covariance_matrix[:]
        if type(transfer_matrix) != type(None):
            t_mat[:] = transfer_matrix[:]
        
        #tt[:] = np.float_(aux.TOTAL_OPTICAL_DEPTH[-1])
        #fi[:] = np.float_(aux.ICE_FRACTION[-1])
        #rl[:] = np.float_(aux.RADIUS_LIQUID[-1])
        #ri[:] = np.float_(aux.RADIUS_ICE[-1])
        
        x_a[:] = np.array([np.float_(aux.TOTAL_OPTICAL_DEPTH[0]), np.float_(aux.ICE_FRACTION[0]), np.float_(aux.RADIUS_LIQUID[0]), np.float_(aux.RADIUS_ICE[0])])
        x_ret[:] = np.array([np.float_(aux.TOTAL_OPTICAL_DEPTH[-1]), np.float_(aux.ICE_FRACTION[-1]), np.float_(aux.RADIUS_LIQUID[-1]), np.float_(aux.RADIUS_ICE[-1])])
        
        #x_a_0 = 1.0/inp.VARIANCE_APRIORI[0]
        
        x_a_err[:] = inp.WEIGHT_APRIORI * np.sqrt(np.reciprocal(np.array(inp.VARIANCE_APRIORI)))
    
        if type(errors) != type(None):
            dtl = np.float_(errors[0])
            dti = np.float_(errors[1])
            drl = np.float_(errors[2])
            dri = np.float_(errors[3])
            x_ret_err[:] = np.array([dtl, dti, drl, dri])
            
            dtl_res = np.float_(errors_res[0])
            dti_res = np.float_(errors_res[1])
            drl_res = np.float_(errors_res[2])
            dri_res = np.float_(errors_res[3])
            x_err_res[:] = np.array([dtl_res, dti_res, drl_res, dri_res])
            
            lwp = np.float_(errors[4])
            iwp = np.float_(errors[6])
            twp = iwp+lwp            
            wp_ret[:] = np.array([lwp, iwp, twp])
            
            dlwp = np.float_(errors[5])
            diwp = np.float_(errors[7])
            dtwp = dlwp+diwp
            wp_ret_err[:] = np.array([dlwp, diwp, dtwp])
            
            dlwp_res = np.float(errors_res[5])
            diwp_res = np.float(errors_res[7])
            dtwp_res = dlwp_res+diwp_res
            wp_err_res[:] = np.array([dlwp_res, diwp_res, dtwp_res])
        
            #dtt[:] = np.float_(errors[0])

            #dfi[:] = np.float_(errors[1])

            #drl[:] = np.float_(errors[2])

            #dri[:] = np.float_(errors[3])

            #lwp[:] = np.float_(errors[4])
            #dlwp[:] = np.float_(errors[5])
            #iwp[:] = np.float_(errors[6])
            #diwp[:] = np.float_(errors[7])
            #twp[:] = np.float_(errors[4])+np.float_(errors[6])
            #dtwp[:] = np.float_(errors[5])+np.float_(errors[7])    
        
    return
