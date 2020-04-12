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
X_PREV = [0.0, 0.0, 0.0, 0.0]
LOOP_PREV = 0

def __retrieve_step(lm_param, loop_count):#, chi2, residuum):
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
    Berechne die Kostenfunktion
    '''
    [chi2, residuum, _res, _apr] = __calc_chi_2_and_residuum()

    '''
    Ueberpruefe, ob die neue Kostenfunktion kleiner ist als die vorherige. Falls die neue 
    Kostenfunktion groesser ist, verwirf die aktuellen Parameter und berechne mit dem vorherigen Parameter
    und einem hoeheren Levenberg-Marquardt-Parameter den neuen Vektor s_n. Ansonsten berechne die Anpassung
    von den aktuellen Werten ausgehend und verringere den Levenberg-Marquardt-Parameter
    '''
    ALPHA = 1.0
    log.write("# Current X^2: {} + {} = {}".format(_res, _apr, chi2))
    if loop_count > 0:
        log.write("# Prev X^2: {}".format(aux.CHI2[-1])) 
    if loop_count == 0 or chi2 <= aux.CHI2[-1]:
        aux.CHI2.append(chi2)
        aux.RESIDUUM.append(residuum)
        #ALPHA = 1.0
        if lm_param > inp.LM_MIN:# or conv_test < 1.0:
            lm_param = lm_param / 2.0
    elif chi2 > aux.CHI2[-1]:
        lm_param = lm_param*4.0
        aux.CHI2.append(aux.CHI2[-1])8Kapitel 1: Nichtlineare Ausgleichsprobleme1.3  LEVENBERG-MARQUARDT-VERFAHRENDas Problem des Gauß-Newton-Verfahrens ist, dass die Linearisierung nur f ̈ur”kleine“ Schritteskzul ̈assig ist, sofern die Iteriertexknoch weit vom Minimum entfernt ist. Die folgende Methodeversucht durch eine Steuerung der L ̈ange des Korrekturvektorsskeine Verbesserung zu erzielen.Im Levenberg-Marquardt-Verfahren wird das Ausgleichsproblem (1.5) zur Bestimmung des Kor-rekturvektorsskdurch ein anderes leicht abge ̈andertes Minimierungsproblem‖F′(xk)sk+F(xk)‖22+μ2‖sk‖22→min(1.9)ersetzt, wobeiμ>0ein zu w ̈ahlender Parameter ist. Als neue Ann ̈aherung wird dann wiederumxk+1=xk+skgesetzt. Aus der Gleichung∥∥∥∥(F′(xk)μI)sk+(F(xk)0)∥∥∥∥22=‖F′(xk)sk+F(xk)‖22+μ2‖sk‖22folgt, dass die Minimierungsaufgabe (1.9) die ̈aquivalente Formulierung∥∥∥∥(F′(xk)μI)sk+(F(xk)0)∥∥∥∥2→min(1.10)besitzt. Im Vergleich zum Ausgleichsproblem (1.5) im Gauß-Newton-Verfahren besitzt (1.10) im-mer eine eindeutige L ̈osungsk, da die Matrix(F′(xk)μI)vollen Rang hat.F ̈ur den Korrekturvektorsk-=0giltμ2‖sk‖22≤‖F′(xk)sk+F(xk)‖22+μ2‖sk‖22=mins∈Rn{‖F′(xk)s+F(xk)‖22+μ2‖s‖22}≤‖F(xk)‖22und daher‖sk‖2≤‖F(xk)‖2μ.Der Parameterμ>0kann folglich eine D ̈ampfung der Korrekturskbewirken und durch einegeeignete Wahl vonμkann man eine zu große Korrektur vermeiden. Man kann zeigen, dass unterbestimmten Voraussetzungen anFdas Levenberg-Marquardt-Verfahren f ̈ur”hinreichend großes“μkonvergiert. Um Konvergenz zu gew ̈ahrleisten darf manμalso nicht zu klein w ̈ahlen; auf deranderen Seite f ̈uhrt ein großesμaber nur zu einer kleinen Korrektur und somit erh ̈alt man nur sehrlangsame Konvergenz.Im Folgenden wird ein m ̈ogliches Verfahren zur Bestimmung des Parametersμvorgestellt, in demin jedem Schritt ̈uberpr ̈uft wird, obμzu klein oder zu groß gew ̈ahlt ist, undμgegebenenfallsangepasst wird.Sei dazuxk∈Rn,sk=sk(μ)die Korrektur aus (1.9) sowieεμ:=‖F(xk)‖22−‖F(xk+sk)‖2
        aux.RESIDUUM.append(aux.RESIDUUM[-1])
        for num_iter in range(9):
            aux.RADIANCE_LBLDIS[num_iter][-1] = aux.RADIANCE_LBLDIS[num_iter][-2]
        aux.TOTAL_OPTICAL_DEPTH[-1] = np.float_(aux.TOTAL_OPTICAL_DEPTH[-2])
        aux.ICE_FRACTION[-1] = np.float_(aux.ICE_FRACTION[-2])
        aux.RADIUS_LIQUID[-1] = np.float_(aux.RADIUS_LIQUID[-2])
        aux.RADIUS_ICE[-1] = np.float_(aux.RADIUS_ICE[-2])
        aux.T_MATRIX[-1] = aux.T_MATRIX[-2]
        #ALPHA = ALPHA / 2.0

        
    '''
    Calculate the adjustment vector
    '''
    delta = numerical.iteration(aux.RESIDUUM[-1], lm_param, aux.T_MATRIX[-1])
    
    s_n = delta[0]
    t_matrix_new = delta[1]
    cov_matrix = delta[2]
    jacobian_mat = delta[3]
    
    #skipped = False
    
    '''
    Berechne die neuen Parameter
    '''
    this_tt = np.float_(aux.TOTAL_OPTICAL_DEPTH[-1] + ALPHA*s_n[0])
    this_fi = np.float_(aux.ICE_FRACTION[-1] + ALPHA*s_n[1])
    this_rl = np.float_(aux.RADIUS_LIQUID[-1] + ALPHA*s_n[2])
    this_ri = np.float_(aux.RADIUS_ICE[-1] + ALPHA*s_n[3])
    
    #if this_tt < 0.0:
    #    this_tt = 0.0
    #    this_rl = aux.RADIUS_LIQUID[-1]
    #if this_fi < 0.0:
    #    this_fi = 0.0
    #    this_ri = aux.RADIUS_ICE[-1]
    
    '''
    Falls einer der Parameter kleiner als 0 ist, oder ice fraction groesser
    als 1 ist, verwerfe diese Parameter und erhoehe den Levenberg-Marquardt-
    Parameter und bestimme erneut das delta
    '''
    while this_tt < 0.0 or this_fi < 0.0 or this_fi > 100.0 or this_rl < 1.0 or this_ri < 1.0:
        #while this_rl < 1.0 or this_ri < 1.0:
        lm_param = lm_param * 2.0
        delta = numerical.iteration(aux.RESIDUUM[-1], lm_param, aux.T_MATRIX[-1])
        s_n = delta[0]
        t_matrix_new = delta[1]
        cov_matrix = delta[2]
        jacobian_mat = delta[3]
        #ALPHA /= 2.0
        this_tt = np.float_(aux.TOTAL_OPTICAL_DEPTH[-1] + ALPHA*s_n[0])
        this_fi = np.float_(aux.ICE_FRACTION[-1] + ALPHA*s_n[1])
        this_rl = np.float_(aux.RADIUS_LIQUID[-1] + ALPHA*s_n[2])
        this_ri = np.float_(aux.RADIUS_ICE[-1] + ALPHA*s_n[3])
        log.write("# x_n = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(this_tt, this_fi, this_rl, this_ri))
        
    '''
    Fuege die neuen Parameter in die Listen ein
    '''
    
    log.write("# x_prev = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(0.1*aux.TOTAL_OPTICAL_DEPTH[-1], 0.1*aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]))
    log.write("# s_n = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(0.1*np.float_(s_n[0]), 0.1*np.float_(s_n[1]), np.float_(s_n[2]), np.float_(s_n[3])))
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
    #exit(-1)
    if loop_count > 0:
        F_2_x_n = np.linalg.norm(aux.RADIANCE_LBLDIS[0][-2])**2
        F_2_x_n1 = np.linalg.norm(aux.RADIANCE_LBLDIS[0][-1])**2
        F_2_x_n1_series = np.linalg.norm(np.array(aux.RADIANCE_LBLDIS[0][-2]) + np.matmul(np.array(numerical.jacobian(-2)), s_n))**2
        eps = (F_2_x_n - F_2_x_n1) / (F_2_x_n - F_2_x_n1_series)
        log.write("# epsilon = {}\n".format(eps)) 
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
        condition = conv_test < inp.CONVERGENCE and ALPHA == 1.0 and lm_param <= 0.001
        #log.write("{} < {}? {} {}\n".format(conv_test, inp.CONVERGENCE, conv_test < inp.CONVERGENCE, condition))                       

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

            
            return True
    return False

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
    [aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR] = aux.average(aux.WAVENUMBER_FTIR[:], aux.RADIANCE_FTIR[:])
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
    plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_FTIR)
    plt.plot(aux.WAVENUMBER_FTIR, aux.RADIANCE_LBLDIS[0][-1])
    plt.grid(True)
    plt.savefig("radiance_fwd.png")
    plt.close()
    plt.clf()
    exit(-1)
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
    chi2 = np.float_(_res)# + _apr)
    if inp.SEARCH_INIT:
        chi2 = np.float(_res)
    
    return [chi2, residuum, _res, _apr]

################################################################################

def __conv_diagnostics(loop, cov_matrix):
    '''Calculate the convergence value
    
    @param cov_matrix The current covariance matrix
    
    @return The convergence value
    '''
    if loop == 0:#type(cov_matrix) == type(None):
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
    convergence_rodgers = np.abs((aux.CHI2[-1]-aux.CHI2[-2])/aux.CHI2[-1])
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
    #alpha = 1.0
    for retr_loop in range(aux.MAX_ITER):

        if inp.FORWARD:
            __only_fwd()
            return
            
            
        log.write("# Iteration: {}".format(retr_loop))
        log.write("# [{}]".format(dt.datetime.now()))
        log.write("# MCP of the current iteration: ")
        log.write("# MCP = [{:6.3f}, {:6.3f}, {:6.3f}, {:6.3f}]".format(0.1*aux.TOTAL_OPTICAL_DEPTH[-1], 0.1*aux.ICE_FRACTION[-1], aux.RADIUS_LIQUID[-1], aux.RADIUS_ICE[-1]))
        log.write("# Levenberg-Marquardt parameter: {}".format(lm_param))
        
        __run_lbldis_and_derivatives()
    
        [lm_param, cov_matrix, s_n, t_matrix_new] = __retrieve_step(lm_param, retr_loop)#, aux.CHI2[-1], aux.RESIDUUM[-1])

        conv_test = __conv_diagnostics(retr_loop, cov_matrix)
        converged = __convergence(lm_param*10, retr_loop, conv_test)
        if converged:
            return
        
        log.write("# Convergence criterion (must be < {} at lm < {}): {}\n".format(inp.CONVERGENCE, 0.001, conv_test))
        
    return

