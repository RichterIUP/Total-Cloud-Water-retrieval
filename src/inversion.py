#!/usr/bin/python
'''@package docstring
Perform the inversion
'''
import threading as th
import os
import sys
import datetime as dt
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

ALPHA = 1.0

def __retrieve_step(lm_param, loop_count, chi2, residuum):
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

    global INCREASE
    global ALPHA
    
    '''
    Calculate the adjustment vector
    '''
    delta = numerical.iteration(residuum, lm_param, aux.T_MATRIX[-1])
    
    s_n = delta[0]
    t_matrix_new = delta[1]
    cov_matrix = delta[2]

            
    #if s_n[1] + aux.ICE_FRACTION[-1] < 0.0 or s_n[3] + aux.RADIUS_ICE[-1] < 1.0:
    #    s_n[3] = aux.RADIUS_LIQUID[-1]-aux.RADIUS_ICE[-1]
    #    s_n[1] = -aux.ICE_FRACTION[-1]
    #if s_n[0] + aux.TOTAL_OPTICAL_DEPTH[-1] < 0.0:
    #    s_n[0] = -aux.TOTAL_OPTICAL_DEPTH
    #    s_n[2] = aux.RADIUS_ICE[-1]-aux.RADIUS_LIQUID[-1]
    #elif s_n[1] + aux.ICE_FRACTION[-1] > 1.0 or s_n[2] + aux.RADIUS_LIQUID[-1] < 1.0:
    #    s_n[2] = aux.RADIUS_ICE[-1]-aux.RADIUS_LIQUID[-1]
    #    s_n[1] = 1-aux.ICE_FRACTION[-1]
    '''
    Check if the new cost function is smaller then the previous one
    '''

    log.write("# Adjustment vector: ")
    log.write("# s_n = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(np.float_(s_n[0]), np.float_(s_n[1]), np.float_(s_n[2]), np.float_(s_n[3])))

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
        
    if loop_count != 0 or aux.MAX_ITER == 1:                            
        condition = conv_test < inp.CONVERGENCE and conv_test > 0.0 and loop_count > 15
        if loop_count != 0 and condition or loop_count == aux.MAX_ITER-1:

            mcp = [np.float_(aux.TOTAL_OPTICAL_DEPTH[-1]), np.float_(aux.ICE_FRACTION[-1]), \
                    np.float_(aux.RADIUS_LIQUID[-1]), np.float_(aux.RADIUS_ICE[-1])]

            averaging_kernel = numerical.calc_avk(aux.T_MATRIX[-1])
            errors = numerical.calc_error(mcp, aux.T_MATRIX[-1])
            log.write("# Final change of parameters: {}\n".format(conv_test))
            if  loop_count == aux.MAX_ITER-1:
                log.write("Not converged!\n")
                idx = aux.CHI2.index(min(aux.CHI2))
                mcp = [np.float_(aux.TOTAL_OPTICAL_DEPTH[idx]), np.float_(aux.ICE_FRACTION[idx]), \
                    np.float_(aux.RADIUS_LIQUID[idx]), np.float_(aux.RADIUS_ICE[idx])]
                aux.TOTAL_OPTICAL_DEPTH[-1] = aux.TOTAL_OPTICAL_DEPTH[idx]
                aux.ICE_FRACTION[-1] = aux.ICE_FRACTION[idx]
                aux.RADIUS_LIQUID[-1] = aux.RADIUS_LIQUID[idx]
                aux.RADIUS_ICE[-1] = aux.RADIUS_ICE[idx]
                aux.RESIDUUM[-1] = aux.RESIDUUM[idx]
                #averaging_kernel = numerical.calc_avk(aux.T_MATRIX[idx])
                errors = numerical.calc_error(mcp, aux.T_MATRIX[idx])
                log.write("Best estimation: x_{} = ({}, {}, {}, {})\n".format(loop_count, mcp[0], mcp[1], mcp[2], mcp[3])) 
                
                #Perform error estimation using the current step
                if not inp.ONLY_OD:
                    __run_lbldis_and_derivatives()
                    [chi2, residuum, _res, _apr] = __calc_chi_2_and_residuum()
                    delta = numerical.iteration(residuum, 0.0, np.zeros((4, len(aux.WAVENUMBER_FTIR))))
                    averaging_kernel = numerical.calc_avk(delta[1])
                    errors = numerical.calc_error(mcp, delta[1])
                    create_nc.create_nc(chi_2=chi2, avk_matrix=averaging_kernel, errors=errors, index=-1, nc=0)               
            else:
                aux.CONVERGED = True
                log.write("Finished! Final Parameters: x_{} = ({}, {}, {}, {})\n".format(loop_count, mcp[0], mcp[1], mcp[2], mcp[3]))
                if not inp.ONLY_OD:
                    create_nc.create_nc(chi_2=aux.CHI2[-1], avk_matrix=averaging_kernel, errors=errors)

            
            return [True, errors]
    return [False, False]

################################################################################

def __set_up_retrieval():
    '''Initialise the retrieval and the matrizes
    '''
    aux.RADIANCE_LBLDIS = [[],[],[],[],[],[],[],[],[]]
    aux.TOTAL_OPTICAL_DEPTH = []
    aux.ICE_FRACTION = []
    aux.RADIUS_LIQUID = []
    aux.RADIUS_ICE = []
    aux.TOTAL_OPTICAL_DEPTH.append(inp.MCP[0])
    aux.ICE_FRACTION.append(inp.MCP[1])
    aux.RADIUS_LIQUID.append(inp.MCP[2])
    aux.RADIUS_ICE.append(inp.MCP[3])
    
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
    #inp.RESOLUTION = np.mean(np.ediff1d(aux.WAVENUMBER_FTIR))


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
    if inp.MODELFRAMEWORK == "LBLDIS":
        rL.forward_run([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]], [0, 1.0], True, 0)
    
    '''
    If the current spectrum is a testcase, add some noise to the radiances
    '''
    if inp.TESTCASE:
        aux.RADIANCE_FTIR = aux.add_noise()
        
    '''
    Calculate the noise and the S_y matrix
    '''
    plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR)
    [aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR] = aux.average(aux.WAVENUMBER_FTIR[:], aux.RADIANCE_FTIR[:])
    plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR, ".")
    plt.grid(True)
    plt.savefig("averaged.png")
    #exit(-1)
    [variance_ra, aux.S_Y_INV_MATRIX] = aux.calc_noise()

    #aux.RADIANCE_LBLDIS = [[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR],[aux.RADIANCE_FTIR]]

    
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

def __only_fwd(lblrtm=False):
    '''Execute only one forward run
    
    @param lblrtm If true, also execute lblrtm
    '''

    print(inp.MODELFRAMEWORK)
    if inp.MODELFRAMEWORK == "LBLDIS":
        rL.forward_run([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]], [0, 1.0], lblrtm, 0)
    #elif inp.MODELFRAMEWORK == "CLARRA":
    #    rD.forward_run([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]], [0, 1.0], lblrtm, 0)
    f = open("{}/lbldis.spec".format(inp.PATH), "w")
    for ii in range(len(aux.WAVENUMBER_FTIR)):
        f.write("{},{}\n".format(aux.WAVENUMBER_FTIR[ii], aux.RADIANCE_LBLDIS[0][-1][ii]))
    f.close()
    return

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
    if inp.ONLY_OD:
        iter_list = [[0, 1.0], [3, 1.0]]
    else:
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

def __calc_chi_2_and_residuum(idx=0):
    '''Calculate the new cost function
    
    @param chi_2 The current cost function
    @param residuum The residuum y - F(x)
    
    @return The new cost function and the new residuum
    '''

    residuum = numerical.residuum(idx)

    '''
    Calculate (y - F(x))^T S_y_1 (y - F(x))
    '''
    _res = np.float_(np.dot(np.matmul(np.transpose(residuum), aux.S_Y_INV_MATRIX[:]), \
                            residuum))
                            
    apr_vec = np.array(inp.MCP_APRIORI[:]) - np.array([aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]])

    '''
    Calculate (x_a - x_i)^T S_a_1 (x_a - x_i)
    '''
    _apr = np.float_(np.dot(np.matmul(np.transpose(apr_vec), aux.S_A_INV_MATRIX[:]), apr_vec))
    
    '''
    Calculate chi^2 = (y - F(x))^T S_y_1 (y - F(x)) + Calculate (x_a - x_i)^T S_a_1 (x_a - x_i)
    '''
    chi2 = np.float_(_res + _apr)
    if inp.SEARCH_INIT:
        chi2 = np.float(_res)
    
    return [chi2, residuum, _res, _apr]

################################################################################

def __conv_diagnostics(cov_matrix):
    '''Calculate the convergence value
    
    @param cov_matrix The current covariance matrix
    
    @return The convergence value
    '''
    if type(cov_matrix) == type(None):
        return 1e10
    tt = aux.TOTAL_OPTICAL_DEPTH[-2:]
    fi = aux.ICE_FRACTION[-2:]
    rl = aux.RADIUS_LIQUID[-2:]
    ri = aux.RADIUS_ICE[-2:]
    parameters_n = np.array([np.float_(tt[1]), np.float_(fi[1]), np.float_(rl[1]), np.float_(ri[1])])
    parameters_n_1 = np.array([np.float_(tt[0]), np.float_(fi[0]), np.float(rl[0]), np.float(ri[0])])

    '''
    Calculate (x_n - x_n+1)^T
    '''
    x_n_x_n_1 = np.matrix(np.array(parameters_n) - np.array(parameters_n_1))

    '''
    Calculate (x_n - x_n+1)^T J^T S_y_1 (y-F(x)) + S_a_1 (x_a - x_i)
    '''
    convergence_rodgers = np.abs(np.float_(np.matmul(x_n_x_n_1, cov_matrix)))
    return convergence_rodgers
    
################################################################################

def retrieve():
    '''Wrapper for inversion
    
    @return True if converged and the final cost function
    '''

    __set_up_retrieval()
    conv_test = 1000.0
    [lm_param, lm_param_prev] = __initialise_variables()
    cov_matrix = None
    s_n = [0.0, 0.0, 0.0, 0.0]
    alpha = 1.0
    for retr_loop in range(aux.MAX_ITER):

        if inp.FORWARD:
            __only_fwd()
            return
    
        log.write("# Iteration: {}".format(retr_loop))
        log.write("# [{}]".format(dt.datetime.now()))
        log.write("# MCP of the current iteration: ")
        log.write("# MCP = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]))
        log.write("# Levenberg-Marquardt parameter: {}".format(lm_param))
        
        skipped = False
        if aux.RADIUS_ICE[-1] > 100.0 or aux.RADIUS_ICE[-1] < 0.0 or \
            aux.RADIUS_LIQUID[-1] > 100.0 or aux.RADIUS_LIQUID[-1] < 0.0 or \
            aux.ICE_FRACTION[-1] > 100.0 or aux.ICE_FRACTION[-1] < 0.0 or \
            aux.TOTAL_OPTICAL_DEPTH[-1] > 100.0 or aux.TOTAL_OPTICAL_DEPTH[-1] < 0.0:
                skipped = True
        else:
            __run_lbldis_and_derivatives()
        
        [chi2, residuum, _res, _apr] = __calc_chi_2_and_residuum()

        if retr_loop == 0:
            aux.CHI2.append(chi2)
            aux.RESIDUUM.append(residuum)
        else:
            log.write("# Current X^2: {} + {} = {}".format(_res, _apr, chi2))
            log.write("# Prev X^2: {}".format(aux.CHI2[-1])) 
            if chi2 <= aux.CHI2[-1] and not skipped:
                aux.CHI2.append(chi2)
                aux.RESIDUUM.append(residuum)
                if lm_param > inp.LM_MIN or conv_test < 1.0:
                    lm_param = lm_param / 10.0
            else:
                log.write("# X^2_n > X^2_n-1!")
                log.write("# Increase Levenberg-Marquardt parameter")
                log.write("# Discard s_n!")
                lm_param = lm_param * 20.0
                #if lm_param == 0.0 and aux.ENABLE_LM_DURING_ITER:
                #    lm_param = 100.0
                aux.CHI2.append(aux.CHI2[-1])
                aux.RESIDUUM.append(aux.RESIDUUM[-1])
                alpha = alpha / 2.0
                aux.TOTAL_OPTICAL_DEPTH[-1] = np.float_(aux.TOTAL_OPTICAL_DEPTH[-1] - alpha*np.float_(s_n[0]))
                aux.ICE_FRACTION[-1] = np.float_(aux.ICE_FRACTION[-1] - alpha*np.float_(s_n[1]))
                aux.RADIUS_LIQUID[-1] = np.float_(aux.RADIUS_LIQUID[-1] - alpha*np.float_(s_n[2]))
                aux.RADIUS_ICE[-1] = np.float_(aux.RADIUS_ICE[-1] - alpha*np.float_(s_n[3]))
                aux.T_MATRIX[-1] = aux.T_MATRIX[-2]
                log.write("# MCP = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(aux.TOTAL_OPTICAL_DEPTH[-1], aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]))


                nums = 9
                #if inp.ONLY_OD:
                #    nums = 2
                
                #for num_iter in range(nums):
                #    aux.RADIANCE_LBLDIS[num_iter][-1] = aux.RADIANCE_LBLDIS[num_iter][-2]
                continue
            conv_test = __conv_diagnostics(cov_matrix)
            converged = __convergence(lm_param*10, retr_loop, conv_test)
                          

        if retr_loop != 0 and converged[0]:
            rms = np.sqrt(np.mean(np.array(aux.RESIDUUM[-1])**2))
            log.write("# Root-Mean-Squared Error = {}".format(rms))
            log.write("# Costfunction X^2 = {}".format(aux.CHI2[-1]))
            log.write("# Convergence criterion (must be < {}): {}\n".format(inp.CONVERGENCE, conv_test))
            aux.CHI_ADJ = chi2
            return

        rms = np.sqrt(np.mean(np.array(aux.RESIDUUM[-1])**2))
        log.write("# Root-Mean-Squared Error = {}".format(rms))
        log.write("# Costfunction X^2 = {}".format(aux.CHI2[-1]))
        plt.figure()
        plt.plot(aux.WAVENUMBER_FTIR, aux.RESIDUUM[-1])
        plt.savefig("{}/residuum_{}.png".format(inp.PATH, retr_loop), dpi=300)
        plt.close()
        plt.clf()

        [lm_param, cov_matrix, s_n, t_matrix_new] = __retrieve_step(lm_param, retr_loop, aux.CHI2[-1], aux.RESIDUUM[-1])

        aux.TOTAL_OPTICAL_DEPTH.append(np.float_(aux.TOTAL_OPTICAL_DEPTH[-1] + alpha*s_n[0]))
        aux.ICE_FRACTION.append(np.float_(aux.ICE_FRACTION[-1] + alpha*s_n[1]))
        aux.RADIUS_LIQUID.append(np.float_(aux.RADIUS_LIQUID[-1] + alpha*s_n[2]))
        aux.RADIUS_ICE.append(np.float_(aux.RADIUS_ICE[-1] + alpha*s_n[3]))
        aux.T_MATRIX.append(t_matrix_new)
        log.write("# Convergence criterion (must be < {} at lm < {}): {}\n".format(inp.CONVERGENCE, 0.001, conv_test))
        
    return

