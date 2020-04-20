#!/usr/bin/python
'''@package docstring
Perform the inversion
'''
import threading as th
import os
import sys
import datetime as dt
import shutil
sys.path.append("./src")
sys.path.append("./src/prowe/clarragroup-run_lblrtm_disort")
import numpy as np
import scipy.signal      as sig
import numerical
import inp
import aux2 as aux
import log
import run_lbldis as rL
#import run_lblrtm_disort_prowe as rD
import get_atm
import create_nc
import matplotlib.pyplot as plt

'''
Indices to avoid magic numbers
'''
TT = 0
FI = 1
RL = 2
RI = 3


def calculate_epsilon(chi, s_n):

    deriv_x_n =  np.transpose(np.matmul(np.array(np.transpose(np.matrix(numerical.jacobian(-2)))), s_n))
    res_linapprox = np.transpose(np.matrix(np.array(aux.RADIANCE_FTIR[:]) - (np.array(aux.RADIANCE_LBLDIS[0][-2][:]) + deriv_x_n)))
    linear_approx = np.float_(np.dot(np.matmul(np.transpose(res_linapprox), aux.S_Y_INV_MATRIX[:]), \
                        res_linapprox))
    change_of_costfunction = aux.CHI2[-1] - chi
    change_of_costfunction_for_linear_model = aux.CHI2[-1] - linear_approx
    eps = change_of_costfunction / change_of_costfunction_for_linear_model
    log.write("# epsilon = {}".format(eps)) 
    return eps

def __retrieve_step(lm_param, loop_count, s_n):#, chi2, residuum):
    '''
    Calculate the next parameters using Least Squares with Gauss-Newton/Levenberg-Marquardt.
    If the prior cost function is better than the current one, the current elements of the
    lists are deleted. Otherwise the new values are appended
    
    @param res The residuum y - F(x)
    @param chi_2 The cost function
    @param lm_param The Levenberg-Marquardt parameter nu
    @param loop_counter The number of the iteration
    @param t_matrix The transfer matrix
    
    @return The residuum, the cost function, the new Levenberg-Marquardt parameter and the covariance matrix
    '''

    
    '''
    Berechne die Kostenfunktion
    '''
    [chi2, residuum, _res, _apr] = __calc_chi_2_and_residuum()

    '''
    Ueberpruefe, ob die neue Kostenfunktion kleiner ist als die vorherige. Falls die neue 
    Kostenfunktion groesser ist, verwirf die aktuellen Parameter und berechne mit dem vorherigen Parameter
    und einem hoeheren Levenberg-Marquardt-Parameter den neuen Vektor s_n. Ansonsten berechne die Anpassung
    von den aktuellen Werten ausgehend und verringere den Levenberg-Marquardt-Parameter
    '''
    log.write("# Current X^2: {}".format(chi2))
    #if loop_count > 0:
    #    eps = calculate_epsilon(chi2, s_n)
    if loop_count > 0:
        log.write("# Prev X^2: {}".format(aux.CHI2[-1])) 
    if loop_count == 0 or chi2 <= aux.CHI2[-1]:
        eps = 0.5
        if loop_count >= 1:
            eps = calculate_epsilon(chi2, s_n)
        aux.CHI2.append(chi2)
        aux.RESIDUUM.append(residuum)
        if eps < 0.25:
            lm_param = lm_param * aux.INCREASE_LM
            if lm_param == 0.0:
                lm_param = inp.LM_INIT
        elif eps >= 0.25 and eps < 0.75:
            lm_param = lm_param
        elif eps >= 0.75:
            lm_param = lm_param / aux.DECREASE_LM
    elif chi2 > aux.CHI2[-1]:
        lm_param = lm_param * aux.INCREASE_LM
        if lm_param == 0.0:
            lm_param = inp.LM_INIT
        aux.CHI2.append(aux.CHI2[-1])
        aux.RESIDUUM.append(aux.RESIDUUM[-1])
        for num_iter in range(9):
            aux.RADIANCE_LBLDIS[num_iter][-1] = aux.RADIANCE_LBLDIS[num_iter][-2]
        aux.TOTAL_OPTICAL_DEPTH[-1] = np.float_(aux.TOTAL_OPTICAL_DEPTH[-2])
        aux.ICE_FRACTION[-1] = np.float_(aux.ICE_FRACTION[-2])
        aux.RADIUS_LIQUID[-1] = np.float_(aux.RADIUS_LIQUID[-2])
        aux.RADIUS_ICE[-1] = np.float_(aux.RADIUS_ICE[-2])
        aux.T_MATRIX[-1] = aux.T_MATRIX[-2]

        
    '''
    Calculate the adjustment vector
    '''
    delta = numerical.iteration(aux.RESIDUUM[-1], lm_param, aux.T_MATRIX[-1])
    
    s_n = delta[0]
    t_matrix_new = delta[1]
    cov_matrix = delta[2]
    jacobian_mat = delta[3]

    
    '''
    Berechne die neuen Parameter
    '''
    this_tt = np.float_(aux.TOTAL_OPTICAL_DEPTH[-1] + s_n[0])
    this_fi = np.float_(aux.ICE_FRACTION[-1]        + s_n[1])
    this_rl = np.float_(aux.RADIUS_LIQUID[-1]       + s_n[2])
    this_ri = np.float_(aux.RADIUS_ICE[-1]          + s_n[3])
    log.write("# s_n = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(np.float_(s_n[0])/inp.SCALE, np.float_(s_n[1])/inp.SCALE, np.float_(s_n[2]), inp.SCALE*np.float_(s_n[3])))

    if this_tt < 0.0:
        this_tt = 0.0
    if this_fi < 0.0:
        this_fi = 0.0
    
    '''
    Falls einer der Parameter kleiner als 0 ist, oder ice fraction groesser
    als 1 ist, verwerfe diese Parameter und erhoehe den Levenberg-Marquardt-
    Parameter und bestimme erneut das delta
    '''
    while this_rl < 1.0 or this_ri < 1.0:#this_tt < 0.0 or this_fi < 0.0 or this_fi > 100.0 or this_rl < 1.0 or this_ri < 1.0:
        lm_param = lm_param * aux.INCREASE_LM
        if lm_param == 0.0:
            lm_param = inp.LM_INIT
        delta = numerical.iteration(aux.RESIDUUM[-1], lm_param, aux.T_MATRIX[-1])
        s_n = delta[0]
        t_matrix_new = delta[1]
        cov_matrix = delta[2]
        jacobian_mat = delta[3]
        this_tt = np.float_(aux.TOTAL_OPTICAL_DEPTH[-1] + s_n[0])
        this_fi = np.float_(aux.ICE_FRACTION[-1]        + s_n[1])
        this_rl = np.float_(aux.RADIUS_LIQUID[-1]       + s_n[2])
        this_ri = np.float_(aux.RADIUS_ICE[-1]          + s_n[3])
        log.write("# x_n = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(this_tt, this_fi, this_rl, this_ri))
        
    '''
    Fuege die neuen Parameter in die Listen ein
    '''
    
    log.write("# x_prev = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(aux.TOTAL_OPTICAL_DEPTH[-1]/inp.SCALE, aux.ICE_FRACTION[-1]/inp.SCALE, aux.RADIUS_LIQUID[-1], inp.SCALE*aux.RADIUS_ICE[-1]))
    log.write("# s_n = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(np.float_(s_n[0])/inp.SCALE, np.float_(s_n[1])/inp.SCALE, np.float_(s_n[2]), inp.SCALE*np.float_(s_n[3])))
    aux.TOTAL_OPTICAL_DEPTH.append(this_tt)
    aux.ICE_FRACTION.append(this_fi)
    aux.RADIUS_LIQUID.append(this_rl)
    aux.RADIUS_ICE.append(this_ri)
    aux.T_MATRIX.append(t_matrix_new)
    
    rms = np.sqrt(np.mean(np.array(aux.RESIDUUM[-1])**2))
    log.write("# Root-Mean-Squared Error = {}".format(rms))
    plt.plot(aux.WAVENUMBER_FTIR, aux.RESIDUUM[-1], ".")
    plt.grid(True)
    plt.savefig("{}/residuum_{}.png".format(inp.PATH, loop_count))
    plt.close()
    plt.clf()
    plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR, ".")
    plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_LBLDIS[0][-1], ".")
    plt.grid(True)
    plt.savefig("{}/radiance_{}.png".format(inp.PATH, loop_count))
    plt.close()
    plt.clf()

    return [lm_param, cov_matrix, s_n, t_matrix_new]

################################################################################

def __convergence(lm_param, loop_count, conv_test):
    '''Test if convergence reached
    
    @param residuum The residuum y - F(x)
    @param lm_param The current Levenberg-Marquardt parameter
    @param lm_param_prev The previous Levenberg-Marquardt parameter
    @param loop_counter Number of the current iteration
    @param t_matrix Transfer matrix
    @param chi_2 The cost function
    @param conv_test Value of the convergence test 
    
    @returns True if converged, the final MCP and the final cost function
    '''
        
    global ALPHA
    if loop_count != 0 or aux.MAX_ITER == 1:     
        condition = conv_test < inp.CONVERGENCE and conv_test > 0.0# and lm_param < inp.LM_INIT/10.0                  

        if loop_count != 0 and condition or loop_count == aux.MAX_ITER-1:

            mcp = [np.float_(aux.TOTAL_OPTICAL_DEPTH[-1]), np.float_(aux.ICE_FRACTION[-1]), \
                    np.float_(aux.RADIUS_LIQUID[-1]), np.float_(aux.RADIUS_ICE[-1])]

            averaging_kernel = numerical.calc_avk(aux.T_MATRIX[-1], idx=-1)
            errors = numerical.calc_error(mcp, aux.T_MATRIX[-1], variance_matrix=aux.S_Y_INV_MATRIX)
    
            res_error = np.array([np.max(np.abs(aux.RESIDUUM[-1]))**2 for ii in range(len(aux.WAVENUMBER_FTIR))])
            s_res_inv_matrix = np.reciprocal(res_error) * np.identity(len(res_error))
            errors_res = numerical.calc_error(mcp, aux.T_MATRIX[-1], variance_matrix=s_res_inv_matrix)
            
            cov_mat = errors[-1]
            log.write("# Final change of parameters: {}\n".format(conv_test))
            if  loop_count == aux.MAX_ITER-1:
                log.write("Not converged!\n")
                idx = aux.CHI2.index(min(aux.CHI2))-1
                mcp = [np.float_(aux.TOTAL_OPTICAL_DEPTH[idx]), np.float_(aux.ICE_FRACTION[idx]), \
                    np.float_(aux.RADIUS_LIQUID[idx]), np.float_(aux.RADIUS_ICE[idx])]
                aux.TOTAL_OPTICAL_DEPTH[-1] = aux.TOTAL_OPTICAL_DEPTH[idx]
                aux.ICE_FRACTION[-1] = aux.ICE_FRACTION[idx]
                aux.RADIUS_LIQUID[-1] = aux.RADIUS_LIQUID[idx]
                aux.RADIUS_ICE[-1] = aux.RADIUS_ICE[idx]
                aux.RESIDUUM[-1] = aux.RESIDUUM[idx]
                averaging_kernel = numerical.calc_avk(aux.T_MATRIX[idx], idx=idx)
                errors = numerical.calc_error(mcp, aux.T_MATRIX[idx], variance_matrix=aux.S_Y_INV_MATRIX)
                
                res_error = np.array([np.max(np.abs(aux.RESIDUUM[-1]))**2 for ii in range(len(aux.WAVENUMBER_FTIR))])
                s_res_inv_matrix = np.reciprocal(res_error) * np.identity(len(res_error))
                errors_res = numerical.calc_error(mcp, aux.T_MATRIX[idx], variance_matrix=s_res_inv_matrix)

                cov_mat = errors[-1]
                log.write("Best estimation: x_{} = ({}, {}, {}, {})\n".format(loop_count, mcp[0], mcp[1], mcp[2], mcp[3])) 

                create_nc.create_nc(chi_2=aux.CHI2[idx], avk_matrix=averaging_kernel, errors=errors, index=-1, nc=0, covariance_matrix=cov_mat, transfer_matrix=aux.T_MATRIX[idx], errors_res=errors_res)               
            else:
                aux.CONVERGED = True
                log.write("Finished! Final Parameters: x_{} = ({}, {}, {}, {})\n".format(loop_count, mcp[0], mcp[1], mcp[2], mcp[3]))
                create_nc.create_nc(chi_2=aux.CHI2[-1], avk_matrix=averaging_kernel, errors=errors, covariance_matrix=cov_mat, transfer_matrix=aux.T_MATRIX[-1], errors_res=errors_res)

            inp.MCP = mcp
            return True
    return False

################################################################################

def __set_up_retrieval():
    '''Initialise the retrieval and the matrizes
    '''
    aux.RADIANCE_LBLDIS = [[],[],[],[],[],[],[],[],[]]

    aux.LBLTP5 = "{}/tp5_{}".format(inp.PATH, aux.TIME_INDEX)
    aux.LBLTMP = '{}'.format(inp.PATH)
    aux.LBLLOG = '{}/lbllog.txt'.format(inp.PATH)
    aux.LBLDIR = "{}/lblout_{}".format(inp.PATH, aux.FTIR.split("/")[-1])

    '''
    Create the directory for the optical depths of LBLRTM
    '''
    if not os.path.exists("{}".format(aux.LBLDIR)):
        os.mkdir("{}".format(aux.LBLDIR))

    '''
    Prepare the atmospheric data
    '''
    get_atm.get_atm()
       
    '''
    Set up the S_a matrix
    '''    
    num_of_params = len(inp.MCP)
    aux.S_A_INV_MATRIX = inp.WEIGHT_APRIORI * np.array(inp.VARIANCE_APRIORI) * \
    np.identity(num_of_params)
    log.log_prog_start()
    '''
    Change the resolution of the spectrum
    '''
    aux.change_resolution()


    '''
    If the current spectrum is a testcase, convolve is with a boxcar
    ''' 
    if inp.TESTCASE and inp.CONVOLVE:
        ft_boxcar = lambda xx: 2 * (inp.OPD) * np.sinc(2 * (inp.OPD) * xx)
        radiance = np.array(aux.RADIANCE_FTIR)
        wavenumber = np.array(aux.WAVENUMBER_FTIR)
        x_axis = np.array([loop_count for \
                           loop_count in range(len(wavenumber))])
        convolution = sig.convolve(radiance, \
                                   ft_boxcar(x_axis), mode='full')
        convolution = np.array(convolution[:len(radiance)])
        normalisation = max(radiance)/max(convolution)
        radiance = normalisation * convolution
        aux.RADIANCE_FTIR = radiance

    '''
    Calculate the clear sky spectrum
    '''
    #if not os.path.exists(aux.LBLDIR):
    if inp.MODELFRAMEWORK == "LBLDIS":
        rL.forward_run(inp.MCP, [0, 1.0], True, 0)
    
    '''
    If the current spectrum is a testcase, add some noise to the radiances
    '''

    if inp.TESTCASE:
        aux.RADIANCE_FTIR = aux.add_noise()
        
    '''
    Calculate the noise and the S_y matrix
    '''
    [aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR] = aux.average(aux.WAVENUMBER_FTIR[:], aux.RADIANCE_FTIR[:])
    #exit(-1)
    [variance_ra, aux.S_Y_INV_MATRIX] = aux.s_y_inv()
    
    log.log_pre_iter(variance_ra)
    return

################################################################################

def __initialise_variables():
    '''Initialise all variables used for the retrieval
    
    @return Wildcards for: [Levenberg-Marquardt-Parameter, previous Levenberg-Marquardt-Parameter, residuum y-F(x), cost function, transfer matrix]
    '''
    dim = len(inp.MCP)
    aux.T_MATRIX = [np.zeros((dim, len(aux.WAVENUMBER_FTIR)))]

    return [inp.LM_INIT, inp.LM_INIT * 2.0]

################################################################################

def __only_fwd(tau_liq=0.0, tau_ice=0.0, reff_liq=0.0, reff_ice=0.0, lblrtm=False, filenum=0):
    '''Execute only one forward run
    
    @param lblrtm If true, also execute lblrtm
    '''

    idx=7
    [wavenumber, radiance] = rL.forward_run([tau_liq, tau_ice, reff_liq, reff_ice], [filenum, 1.0], lblrtm, 0)
    rms = np.sqrt(np.mean((np.array(radiance) - np.array(aux.RADIANCE_FTIR))**2))
    slope_lbldis = (radiance[0] - radiance[idx])/(wavenumber[0] - wavenumber[idx])
    slope_ftir = (aux.RADIANCE_FTIR[0] - aux.RADIANCE_FTIR[idx])/(aux.WAVENUMBER_FTIR[0] - aux.WAVENUMBER_FTIR[idx])
    slope = np.abs(slope_lbldis - slope_ftir)
    
    return [np.sum(radiance[idx:-1]), np.sum(aux.RADIANCE_FTIR[idx:-1]), rms, slope]

################################################################################

def __run_lbldis_and_derivatives():
    '''Run lbldis for the main cloud parameters and the derivatives
    
    @param lm_param The current Levenberg-Marquardt parameter
    @param lm_param_prev The previous Levenberg-Marquardt parameter
    '''
    lbldis_run = []

    '''
    Set up the list with for the forward run and the derivatives:
    0 -> Forward run
    1 -> Derivative wrt reff liq
    2 -> Derivative wrt reff ice
    3 -> Derivative wrt tau total
    4 -> Derivative wrt f ice
    '''

    iter_list = [[0, 1.0], [1, 1.0], [1, -1.0], [2, 1.0], [2, -1.0], [3, 1.0], [3, -1.0], [4, 1.0], [4, -1.0]]
        
    '''
    Limit the number of the allocated CPUs
    '''
    if inp.NUM_OF_CPU > len(iter_list):
        inp.NUM_OF_CPU = len(iter_list)

    lblrtm = False
    fnum = 0
    if  True:
        for counter in range(len(iter_list)):
            fact_num = iter_list[counter]
            '''
            Run LBLDIS
            '''
            if inp.MODELFRAMEWORK == "LBLDIS":
                lbldis_run.append(th.Thread(target=rL.forward_run, \
                                            args=([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]], fact_num, \
                                                  lblrtm, fnum)))
            #elif inp.MODELFRAMEWORK == "CLARRA":
            #    lbldis_run.append(th.Thread(target=rD.forward_run, \
            #                            args=([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]], fact_num, \
            #                                  lblrtm, fnum)))
            lbldis_run[-1].start()
            fnum = fnum + 1
            if (counter+1)%inp.NUM_OF_CPU == 0:
                for element in lbldis_run:
                    element.join()
                lbldis_run = []

        if lbldis_run:
            for element in lbldis_run:
                element.join()

    return

################################################################################

def __calc_chi_2_and_residuum(idx=-1):
    '''Calculate the new cost function
    
    @param chi_2 The current cost function
    @param residuum The residuum y - F(x)
    
    @return The new cost function and the new residuum
    '''

    residuum = numerical.residuum(idx=-1)

    '''
    Calculate (y - F(x))^T S_y_1 (y - F(x))
    '''
    _res = np.float_(np.dot(np.matmul(np.transpose(residuum), aux.S_Y_INV_MATRIX[:]), \
                            residuum))
                            
    apr_vec = np.array(inp.MCP_APRIORI[:]) - np.array([aux.TOTAL_OPTICAL_DEPTH[idx], aux.ICE_FRACTION[idx], aux.RADIUS_LIQUID[idx], aux.RADIUS_ICE[idx]])

    '''
    Calculate (x_a - x_i)^T S_a_1 (x_a - x_i)
    '''
    _apr = np.float_(np.dot(np.matmul(np.transpose(apr_vec), aux.S_A_INV_MATRIX[:]), apr_vec))
    
    '''
    Calculate chi^2 = (y - F(x))^T S_y_1 (y - F(x)) + Calculate (x_a - x_i)^T S_a_1 (x_a - x_i)
    '''
    chi2 = np.float_(_res + _apr)
    
    return [chi2, residuum, _res, _apr]

################################################################################

def __conv_diagnostics(loop, cov_matrix):
    '''Calculate the convergence value
    
    @param cov_matrix The current covariance matrix
    
    @return The convergence value
    '''
    if loop == 0:
        return 1e10

    convergence_rodgers = np.abs((aux.CHI2[-1]-aux.CHI2[-2])/aux.CHI2[-1])
    return convergence_rodgers
    
################################################################################

def retrieve():
    '''Wrapper for inversion
    
    @return True if converged and the final cost function
    '''

    aux.RADIANCE_LBLDIS = [[],[],[],[],[],[],[],[],[]]
    aux.TOTAL_OPTICAL_DEPTH = [inp.MCP[0]]
    aux.ICE_FRACTION = [inp.MCP[1]]
    aux.RADIUS_LIQUID = [inp.MCP[2]]
    aux.RADIUS_ICE = [inp.MCP[3]]
    conv_test = 1000.0
    [lm_param, lm_param_prev] = __initialise_variables()
    cov_matrix = None
    s_n = [0.0, 0.0, 0.0, 0.0]
    for retr_loop in range(aux.MAX_ITER):

        if inp.FORWARD:
            return __only_fwd(tau_liq=aux.TOTAL_OPTICAL_DEPTH[-1], \
                                tau_ice=aux.ICE_FRACTION[-1], \
                                reff_liq=aux.RADIUS_LIQUID[-1], \
                                reff_ice=aux.RADIUS_ICE[-1])
            
            
        log.write("# Iteration: {}".format(retr_loop))
        log.write("# [{}]".format(dt.datetime.now()))
        log.write("# MCP of the current iteration: ")
        log.write("# MCP = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(aux.TOTAL_OPTICAL_DEPTH[-1]/inp.SCALE, aux.ICE_FRACTION[-1]/inp.SCALE, aux.RADIUS_LIQUID[-1], inp.SCALE*aux.RADIUS_ICE[-1]))
        log.write("# Levenberg-Marquardt parameter: {}".format(lm_param))
        
        __run_lbldis_and_derivatives()
    
        [lm_param, cov_matrix, s_n, t_matrix_new] = __retrieve_step(lm_param, retr_loop, s_n)

        conv_test = __conv_diagnostics(retr_loop, cov_matrix)
        converged = __convergence(lm_param*10, retr_loop, conv_test)
        if converged:
            return
        
        log.write("# Convergence criterion (must be < {} at lm < {}): {}\n".format(inp.CONVERGENCE, 0.001, conv_test))
        
    return

