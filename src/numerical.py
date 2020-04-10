#!/usr/bin/python
'''@package Docstring
Some numerical functions
'''

#Standard-lib
import sys
import os

#Third-party modules
import numpy          as np

sys.path.append("./src/")

#Self-defined modules
import inp
import aux2 as aux
import read_database as mie

def residuum(idx=0):
    '''Calculation of the residuum
    
    @return The residuum y - F(x)
    '''
    
    res = np.array(aux.RADIANCE_FTIR[:]) - np.array(aux.RADIANCE_LBLDIS[idx][-1][:])
    return np.transpose(np.matrix(res))

def jacobian():
    '''Calculation of the jacobian_matrix
    
    @return The transposed jacobian matrix
    '''

    if inp.ONLY_OD:
        deriv_r_liq = np.array([0.0 for ii in range(len(aux.WAVENUMBER_FTIR))])
        deriv_r_ice = np.array([0.0 for ii in range(len(aux.WAVENUMBER_FTIR))])
        deriv_f_ice = np.array([0.0 for ii in range(len(aux.WAVENUMBER_FTIR))])
        deriv_tau_total = (np.array(aux.RADIANCE_LBLDIS[1][-1]) \
                     - np.array(aux.RADIANCE_LBLDIS[0][-1]))/(aux.STEPSIZE)
        
    else:
        deriv_r_liq = (np.array(aux.RADIANCE_LBLDIS[1][-1]) \
                       - np.array(aux.RADIANCE_LBLDIS[0][-1]))/(aux.STEPSIZE)
        deriv_r_ice = (np.array(aux.RADIANCE_LBLDIS[3][-1]) \
                       - np.array(aux.RADIANCE_LBLDIS[0][-1]))/(aux.STEPSIZE)
        deriv_tau_total = (np.array(aux.RADIANCE_LBLDIS[5][-1]) \
                        - np.array(aux.RADIANCE_LBLDIS[0][-1]))/(aux.STEPSIZE)
        deriv_f_ice = (np.array(aux.RADIANCE_LBLDIS[7][-1]) \
                       - np.array(aux.RADIANCE_LBLDIS[0][-1]))/(aux.STEPSIZE)
                       
        counter = 0
        while os.path.exists("deriv_r_liq_{}".format(counter)):
            counter += 1
        
        f = open("deriv_r_liq_{}".format(counter), "w")
        g = open("deriv_r_ice_{}".format(counter), "w")
        h = open("deriv_t_tot_{}".format(counter), "w")
        i = open("deriv_f_ice_{}".format(counter), "w")
        for ii in range(len(deriv_r_liq)):
            f.write("{},{},{},{}\n".format(aux.WAVENUMBER_FTIR[ii], aux.RADIANCE_LBLDIS[1][-1][ii], aux.RADIANCE_LBLDIS[0][-1][ii], aux.RADIANCE_LBLDIS[1][-1][ii]-aux.RADIANCE_LBLDIS[0][-1][ii]))
            g.write("{},{},{},{}\n".format(aux.WAVENUMBER_FTIR[ii], aux.RADIANCE_LBLDIS[3][-1][ii], aux.RADIANCE_LBLDIS[0][-1][ii], aux.RADIANCE_LBLDIS[3][-1][ii]-aux.RADIANCE_LBLDIS[0][-1][ii]))
            h.write("{},{},{},{}\n".format(aux.WAVENUMBER_FTIR[ii], aux.RADIANCE_LBLDIS[5][-1][ii], aux.RADIANCE_LBLDIS[0][-1][ii], aux.RADIANCE_LBLDIS[5][-1][ii]-aux.RADIANCE_LBLDIS[0][-1][ii]))
            i.write("{},{},{},{}\n".format(aux.WAVENUMBER_FTIR[ii], aux.RADIANCE_LBLDIS[7][-1][ii], aux.RADIANCE_LBLDIS[0][-1][ii], aux.RADIANCE_LBLDIS[7][-1][ii]-aux.RADIANCE_LBLDIS[0][-1][ii]))
        f.write("{}\n".format(np.sqrt(np.mean((aux.RADIANCE_LBLDIS[1][-1]-aux.RADIANCE_LBLDIS[0][-1])**2))))
        g.write("{}\n".format(np.sqrt(np.mean((aux.RADIANCE_LBLDIS[1][-1]-aux.RADIANCE_LBLDIS[0][-1])**2))))
        h.write("{}\n".format(np.sqrt(np.mean((aux.RADIANCE_LBLDIS[1][-1]-aux.RADIANCE_LBLDIS[0][-1])**2))))
        i.write("{}\n".format(np.sqrt(np.mean((aux.RADIANCE_LBLDIS[1][-1]-aux.RADIANCE_LBLDIS[0][-1])**2))))

        f.close()
        g.close()
        h.close()
        i.close()

    return [deriv_tau_total, deriv_f_ice, deriv_r_liq, deriv_r_ice]

def calc_avk(t_matrix):
    '''Calculation of the averaging kernel matrix
    
    @param t_matrix The transfer matrix (Ceccherini/Ridolfi 2010)
    @return The averaging kernel matrix
    '''
    jacobian_matrix_transposed = np.matrix(jacobian())
    jacobian_matrix = np.transpose(jacobian_matrix_transposed)
    averaging_kernel = np.matmul(t_matrix, jacobian_matrix)
    return averaging_kernel

def calc_vcm(t_matrix):
    '''Calculation of the variance-covariance matrix
    
    @param t_matrix The transfer matrix (Ceccherini/Ridolfi 2010)
    @return The variance-covariance matrix
    '''
    s_y_matrix = np.linalg.inv(aux.S_Y_INV_MATRIX)
    cov = np.matmul(t_matrix, s_y_matrix)
    cov = np.matmul(cov, np.transpose(t_matrix))
    return cov

def iteration(res, lm_param, t_matrix):
    '''Calculate the adjustment vector
    
    This function solves the equation
    [J^T S_y_1 J + S_a_1 + mu**2 D]s_n = J^T S_y_1 [y - F(x_n)] + S_a_1 (x_a - x_i)
    w.r.t. s_n, the so-called adjustment vector
    
    @param res The residuum y - F(x)
    @param lm_param The Levenberg-Marquardt parameter mu
    @param t_matrix The transfer matrix (Ceccherini/Ridolfi 2010)
    
    @return The adjustment of the MCP, the new transfer matrix and covariance matrix for the calculation of the convergence
    '''
    s_y_inv_matrix = aux.S_Y_INV_MATRIX[:]
    s_a_inv_matrix = aux.S_A_INV_MATRIX[:]
    apr_vec = inp.MCP_APRIORI[:]
    dim = 4#len(apr_vec)
    jac = jacobian()
    atm_param = np.array([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]])
    
    jacobian_matrix_transposed = np.matrix(jac)
    jacobian_matrix = np.transpose(jacobian_matrix_transposed)

    x_i = np.array(atm_param)

    '''
    Calculate J^T S_y_1
    '''
    JT_W = np.matmul(jacobian_matrix_transposed, s_y_inv_matrix)

    '''
    Calculate J^T S_y_1 J
    '''
    JT_W_J = np.matmul(JT_W, jacobian_matrix)

    '''
    Calculate J^T S_y_1 (y-F(x))
    '''
    JT_W_r = np.matmul(JT_W, res)

    '''
    Calculate D
    '''
    lm_matrix = np.identity(dim)

    '''
    Calculate S_a_1 (x_a - x_i)
    '''
    second_sum = np.matmul(s_a_inv_matrix, apr_vec - x_i)

    '''
    Calculate J^T S_y_1 J + S_a_1 + mu**2 D
    '''
    left_side = np.add(np.add(JT_W_J, s_a_inv_matrix), lm_param**2*lm_matrix)

    '''
    Calculate J^T S_y_1 (y-F(x)) + S_a_1 (x_a - x_i)
    '''
    right_side = \
    np.transpose(np.add(np.transpose(JT_W_r), second_sum))

    '''
    Solve the equation
    '''
    s_n = np.linalg.solve(left_side, right_side)

    '''
    Calculate the transfer matrix
    '''
    try:
        M = np.linalg.inv(left_side)
        G = np.matmul(M, JT_W)
        I_GJ_MR = np.identity(dim) - np.matmul(G, jacobian_matrix) \
                - np.matmul(M, s_a_inv_matrix)
        T_new = G + np.matmul(I_GJ_MR, t_matrix)
    except ValueError:
        T_new = t_matrix
    
    return [s_n, T_new, right_side]

def calc_error(atmospheric_param, t_matrix):
    '''Perform error propagation
    
    @param atmospheric_param The MCP
    @param t_matrix The transfer matrix
    
    @return The standard deviations of the MCP
    '''

    tau_total = atmospheric_param[0]
    f_ice = atmospheric_param[1]
    reff_liq = atmospheric_param[2]
    reff_ice = atmospheric_param[3]
    cov = calc_vcm(t_matrix)
    tau_t_var = np.sqrt(cov.item((0, 0)))
    f_i_var = np.sqrt(cov.item((1, 1)))
    ref_l_var = np.sqrt(cov.item((2, 2)))
    ref_i_var = np.sqrt(cov.item((3, 3)))
    
    '''
    Convert tau_total and f_ice to tau_liquid and tau_ice
    '''
    tau_liq = tau_total * (1 - f_ice)
    tau_ice = tau_total * f_ice
    tau_l_var = np.abs(tau_t_var * (1 - f_ice)) + np.abs(- tau_total * f_i_var)
    tau_i_var = np.abs(tau_t_var * f_ice) + np.abs(tau_total * f_i_var)

    ice_database_shapes = [None for count in range(8)]
    ice_water_path = [0.0 for count in range(8)]
    dice_water_path = [0.0 for count in range(8)]

    ice_database_shapes[0] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                             "./ssp/ssp_db.mie_ice.gamma_sigma_0p100")[1]
    ice_database_shapes[1] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                                "./ssp/ssp_db.Aggregate.gamma.0p100")[1]
    ice_database_shapes[2] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                                "./ssp/ssp_db.BulletRosette.gamma.0p100")[1]
    ice_database_shapes[3] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                              "./ssp/ssp_db.Droxtal.gamma.0p100")[1]
    ice_database_shapes[4] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                                "./ssp/ssp_db.HollowCol.gamma.0p100")[1]
    ice_database_shapes[5] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                            "./ssp/ssp_db.Plate.gamma.0p100")[1]
    ice_database_shapes[6] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                               "./ssp/ssp_db.SolidCol.gamma.0p100")[1]
    ice_database_shapes[7] = mie.read_databases("./ssp/ssp_db.mie_wat.gamma_sigma_0p100", \
                                               "./ssp/ssp_db.Spheroid.gamma.0p100")[1]

    '''
    Calculate the error for IWP
    '''
    if inp.BUILTIN:
        if reff_ice < 16.0:
            idx = 3
        else:
            idx = 6
        [ice_water_path[idx], dice_water_path[idx]] = mie.calc_iwp(tau_ice, tau_i_var, \
        reff_ice, ref_i_var, ice_database_shapes[idx])
    else:
        for idx in range(8):
            ice_water_path_and_error = mie.calc_iwp(inp.SHAPES[idx]*tau_ice, \
                                                    inp.SHAPES[idx]*tau_i_var, \
                                                    reff_ice, ref_i_var, \
                                                    ice_database_shapes[idx])
            ice_water_path[idx] = ice_water_path_and_error[0]
            dice_water_path[idx] = ice_water_path_and_error[1]

    total_ice_water_path = np.sum(ice_water_path)
    total_dice_water_path = np.sum(dice_water_path)

    '''
    Calculate the error for LWP
    '''
    liq_water_path_and_error = mie.calc_lwp(tau_liq, tau_l_var, reff_liq, ref_l_var)
    total_liq_water_path = liq_water_path_and_error[0]
    total_dliq_water_path = liq_water_path_and_error[1]

    return [tau_t_var, f_i_var, ref_l_var, ref_i_var, \
            total_liq_water_path, total_dliq_water_path, \
            total_ice_water_path, total_dice_water_path]
