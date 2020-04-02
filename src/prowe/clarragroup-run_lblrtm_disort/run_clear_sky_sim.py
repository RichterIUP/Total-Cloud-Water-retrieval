"""
run_clear_sky_sims.py

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose:
  Run LBLRTM to get clear sky radiances and optical depths

Inputs:

"""

# .. Built-in Modules
import numpy as np
from shutil import copyfile
import scipy.io as spio
import os.path
from subprocess import Popen
import glob


# .. My modules 
from reduce_resolution import reduce_resolution_pad as reduce_res
from lblrtm_utils import read_od_lblrtm
from lblrtm_utils import Prof
from lblrtm_utils import write_tape5_prof

class ClrSky:
    __slots__ = ('date',
                 'view_angle',
                 'od_file',
                 'nu',
                 'rads',
                 'tsc',
                 'rad_above',
                 'rad_clear',
                 'Bctc',
                 'dt_dtau',
                 )
                 


    # .. Initialize by getting clear sky simulation data for optical depths
    def __init__(self, date, view_angle, od_file):
        
        self.date = date
        self.view_angle = view_angle
        self.od_file = od_file
    
    
    def load_from_file(self):
        '''
        Load in the low resolution clear-sky data from Octave:
              nu: wavenumber vector for same range as mrad_nu0
              rads:   upwelling radiances from layer to TOA
              tsc:    layer-to-TOA low-res transmittance
              rad_above: radiance from tropopause to TOA x trans(tropopause)
        
        Obsolete
              Bctc:  Planck function times layer-to-TOA low-res trans
              odlyr:  layer low-res optical depths
              radClr: Clear sky radiance from TAPE12 file, zenith view
        
        '''

        self.date, self.view_angle, \
        self.nu, self.rads, \
        self.tsc, self.rad_above, \
        self.Bctc, self.dt_dtau = load_od_gas(self.od_file)
        
        # .. If nu is a matrix, convert to a vector
        if len(self.nu.shape)>1:
          self.nu.shape = np.size(self.nu)
          
        self.rad_clear = self.rads[:,-1] + self.rad_above

    def set_vals(self, view_angle, nu, rads, tsc, dt_dtau, Bc_tsc, rad_above):
        self.view_angle = view_angle
        self.nu = nu
        self.rads = rads
        self.tsc = tsc
        self.dt_dtau = dt_dtau  
        self.Bctc = Bc_tsc
        self.rad_above = rad_above
        self.rad_clear = self.rads[:,-1] + self.rad_above
        
    def restrict_nus(self, nu_new):

        mask_nu = (self.nu >= nu_new[0]-1e-5) & \
                  (self.nu <= nu_new[-1]+1e-5)
        self.nu = self.nu[mask_nu]
        self.rads = self.rads[mask_nu,:]
        self.rad_above = self.rad_above[mask_nu]
        self.tsc = self.tsc[mask_nu,:]
        self.Bctc = self.Bctc[mask_nu,:]
        self.dt_dtau = self.dt_dtau[mask_nu,:]
        self.rad_clear = self.rad_clear[mask_nu]
        


def load_od_gas(odfile):
    '''
    Purpose:
      Load in the matlab generated file because python is
      too slow for cubic interp and gives small differences
      But the python way is saved for reference in extras.py
    '''
    
    odinfo         = spio.loadmat(odfile)
    date           = odinfo['date']
    nu             = odinfo['nu'][0]
    rads           = odinfo['rads']
    rad_above      = odinfo['rad_above'][0]
    tsc            = odinfo['tsc']
    view_angle     = odinfo['view_angle']
    Bctc           = odinfo['Bc_tsc']
    dt_dtau        = odinfo['dt_dtau']
    
    return date, view_angle, nu, rads, tsc, rad_above, Bctc, dt_dtau


def run_clear_sky_sim(od_file, profDir, dir_lblrtm, lblrtm_command,
                      observer_alt_km, z_trop, bwn, ewn, 
                      NptsInLoResIntrfgrm, Npad,
                      date, dnu, prof_file, zangle, sim_toa,
                      rerun, save_results_to_disk, print_stdout):
    ''''
    Purpose: Run LBLRTM to compute gaseous (clear sky) optical depths
    Input variables:
      ip: input variables
      date, dnu, prof_file,
      zangle, sim_toa: dividing height between upper and lower atmosphere
      proc: process number for running LBLRTM (corresponds to directory)
      print_stdout: if set to true or 1, prints info to standard output
    Notes:
      The lower atmosphere, defined by z_trop, is run for the date and
      prof_file. Runs are done for all zenith viewing angle zangle.
         
    '''
                        
    # LBLRTM directory
    dir_od = dir_lblrtm    
    run_lblrtm = dir_lblrtm + lblrtm_command

    # .. Set parameters (don't change)
    date_fmt = "%Y-%m-%d %H:%M:%S"
    
        
    # .. Do the radiative transfer calculations
    #    to get the simulated radiances, optical depths, and transmittances
    
    od_file_exists = False 
    # .. Get the optical depth file names
    # .. Get the extention for the optical depth file name
    #    If viewing angle not zero, optical depth file name must
    #    include an underscore followed by the viewing angle
    if print_stdout: print('View angle:', zangle)
    if (zangle >= 90) or (zangle < 0):
        raise NameError('Viewing angle must be between 0 and 90,' + \
                        ' but is ' + str(zangle) )
        
    # .. Initialize clear sky class
    clrSky = ClrSky(date, zangle, od_file)
    
    # .. Check if it exists
    if os.path.isfile(od_file):
        od_file_exists = True
        
        
    # .. Set values in "prof" based on file prof_file to calculate 
    #    optical depths with LBLRTM from bwn to ewn up to height at zenith. 
    #    and get number of model layers to model TOA and tropopause
    prof = Prof(profDir + prof_file)
    if sim_toa < prof.zm[-1]:
        prof.crop_heights(0, sim_toa)
        
    # .. Round altitudes to 3 decimal places, since that is all LBLRTM
    #    accepts
    prof.zm = np.round(prof.zm,3)
    
    # .. Add observer altitude to profile, if necessary
    prof.include_observer_alt(observer_alt_km)

    if (not rerun) and od_file_exists:
        if print_stdout: print('Loading in existing data from ' + od_file )
        clrSky.load_from_file()
        return clrSky, prof
    if print_stdout: print('Doing LBLRTM run for', prof_file)


    
    
    # .. Do LBLRTM run for prof (zenith view) and copy results to dir_od
    #    Go one wavenumber past limits because of round-off differences
    do_lblrtm_run(dir_lblrtm, run_lblrtm, prof, bwn-1., ewn+1., dir_od)

    
    
    # .. Prepare to do radiative transfer for cosine of viewing angles
    #    and get number of layers up to the tropopause
    cos_angle = np.cos(np.deg2rad(zangle))
    nlyr_toa = len(prof.zm) - 1
    nlyr_trop = (np.abs(prof.zm - z_trop)).argmin()

    # .. Print stuff to screen, if directed
    if print_stdout:
        print('nlyr_toa = ', nlyr_toa)
        print('nlyr_trop = ', nlyr_trop)
    
    # .. Set boolean for whether we are going to tropopause or TOA
    if prof.zm[-1] < z_trop:
        raise NameError('Calculation must at least go up to z_trop, ' \
                        + str(z_trop) + ' km')
    elif nlyr_toa > nlyr_trop:
        get_upper_atm = True
    else:
        get_upper_atm = False


    # .. Load in the optical depths and calculate:
    #    rads: finite-resolution layer radiances up to tropopause
    #    tsc: finite-resolution surface-to-layer transmittances
    #    dt_dtau: fin-res change in transmittance with change in optical depth
    #    nu_trop_mono: monochromatic res wavenumbers for
    #    rad_hi_mono: monochromatic radiance for total troposphere
    #    Bc_tsc: Planck function of temperature x tsc
    #    Slow step (~46 - 90 s)
    nu, rads, tsc, dt_dtau, \
    nu_trop_mono, t_trop_mono, \
    nu_mono, rad_hi_mono, Bc_tsc = \
      radtran_multi_angle(dir_lblrtm, nlyr_trop, nlyr_toa,
                          prof.tm, 0, dnu, bwn, ewn, 
                          NptsInLoResIntrfgrm, Npad, prof.zm, 
                          cos_angle, 'd', get_upper_atm)
    
    # .. Save results to optical depth files
    if get_upper_atm:
        # .. Get values for this case
        rad_above_mono = rad_hi_mono[:,0] * np.interp(nu_mono,
                                                      nu_trop_mono, 
                                                      t_trop_mono[:,0])
        nua, rad_above = reduce_res(nu_mono, rad_above_mono, dnu,
                                    NptsInLoResIntrfgrm, Npad, 2**23)
        rad_above = np.real(rad_above[np.logical_and(nua>=bwn, 
                                                     nua<=ewn)])

        if save_results_to_disk:
            # .. Save these results
            mdict = {
                     "view_angle": zangle,
                     "nu": nu, 
                     "rads": rads[:,:,0],
                     "tsc": tsc[:,:,0], 
                     "dt_dtau": dt_dtau[:,:,0],
                     "Bc_tsc": Bc_tsc[:,:,0],
                     "date": date.strftime(date_fmt),
                     "rad_above": rad_above,
                    }
                    # "nu_trop_mono": nu_trop_mono, 
                    # "t_trop_mono": t_trop_mono[:,0],
            spio.savemat(od_file, mdict)
        
        # .. Put results in class and return it.
        clrSky.set_vals(zangle, nu, rads[:,:,0], tsc[:,:,0], dt_dtau[:,:,0], 
                        Bc_tsc[:,:,0], rad_above)
             
    else:
        # .. If we are only going up to the tropopause, 
        #    save everything, including tropopause transmittance,
        #    which will be used and deleted later
        #    WARNING: IN THIS CASE, tsc WILL ONLY GO UP TO
        #    TROPOPAUSE. UPPER LEVEL tsc NEEDS TO BE PATCHED IN!
        rad_above = 0 * nu
        if save_results_to_disk:
            mdict = {
                     "view_angle": zangle,
                     "nu": nu, 
                     "rads": rads[:,:,0],
                     "tsc": tsc[:,:,0], 
                     "dt_dtau": dt_dtau[:,:,0],
                     "Bc_tsc": Bc_tsc[:,:,0],
                     "date": date.strftime(date_fmt),
                     "rad_above": rad_above,
                    }
        
            spio.savemat(od_file, mdict)
            
       
        # .. Put results in class and return it.
        clrSky.set_vals(zangle, nu, rads[:,:,0], tsc[:,:,0], dt_dtau[:,:,0], 
                        Bc_tsc[:,:,0], rad_above)

    return clrSky, prof
   


def do_lblrtm_run(dir_lblrtm, run_lblrtm, prof, v1, v2, dir_od):
    
    # .. LBLRTM calculation (zenith view)
    #    Clean out the LBLRTM directory
    for filen in glob.glob(dir_lblrtm + "ODdeflt_*"):
        os.remove(filen)
    for filen in glob.glob(dir_lblrtm + "TAPE12"):
        os.remove(filen)
    for filen in glob.glob(dir_lblrtm + "TAPE5"):
        os.remove(filen)
    
    # .. Write the TAPE5 file
    write_tape5_prof(prof, v1, v2, 0, 'od', dir_lblrtm)
    
    
    # .. Run LBLRTM (slow!)
    p = Popen([run_lblrtm], cwd = dir_lblrtm)
    p.wait()

    # Move everything you could want to output directory
    if dir_lblrtm != dir_od:
        for filen in glob.glob(dir_lblrtm + "ODdeflt_*"):
            copyfile('ODdeflt*', dir_od + filen)
        copyfile('TAPE5', dir_od + 'TAPE5')
        # movefile('TAPE6', dir_od + 'TAPE6')
        # movefile('TAPE7', dir_od + 'TAPE7')



def radtran_multi_angle(od_dir, nlyr_trop, nlyr_toa, tzl, tave,
  dnu, bwn, ewn, N, Npad, zs, cos_angles, prec, get_upper_atm):
    """
    Copyright Nov. 30, 2014 by Dr. Penny M. Rowe and NorthWest Research
      Associates. All rights reserved.
    
    Purpose:
       Read in LBLRTM optical depths and get the following quantities:

    Outputs (all at instrument resolution unless specified):
      nu:           Final wavenumber
      rads:         Surface-to-layer radiances: (0=> surface-to-surface)
                    from surface-to-surface (R[:,0]) to 
                    surface-to-tropopause (R[:,nlyr_trop])
      Bc_tc:        Planck function of T times transmittance below
                    (as for rads)
      tc:           Surface-to-layer transmittance, as for rads
      dt_dtau:
      nu_mono:      Monochromatic wavenumber, from lblrtm
      t_trop_mono:  Monochromatic trans from surface to troposphere
      rad_hi_mono:  Monochromatic radiance from tropopshere to TOA
      
    Important note on layering:
       The variables rads, Bc_tc, and tc go from the surface to the
       top boundary of the layer, starting from surface to top of first
       layer. This is a change from previously, for which the bottom
       radiance was always zero, bottom transmittance one. This will
       require changes to CO2 slicing (e.g considering a cloud at the 
       surface). Thus there are nlyr_trop values in rads and tc 
       
       For optical depth, odlyr is defined to be from the bottom of
       the layer (not from the surface) to the top of the layer
    
    Inputs in function call:
      od_dir:    directory with LBLRTM optical depths
      nlyr_trop: number of layers up to troposphere
      nlyr_toa:   number of layers desired for optical depths and transmittances
      tzl:       temperature at layer bottom
      tave:      "ave" temperature, as defined by LBLRTM
      dnu:       wavenumber spacing (determined by resolution)
      bwn:       beginning wavenumber
      ewn:       ending wavenumber
      cos_angle: cosine of viewing angle, defined as 0 = view straight up
      prec:      precision for   read_od_lblrtm
      
    
    Othere Variables used in this code:
         od_mono: monochromatic layer optical depth, from lblrtm
         nu_L: layer wavenumber, at instrument resolution
         odtot: surface-to-layer optical depth, at instrument resolution
         trans: surface-to-layer transmittance, at instrument resolution
         transLm1: above for surface-to-layer-minus-1
         
    """
  
    # ... QC
    if dnu <= 0:
        raise NameError('Wavenumber spacing must be >= 0')
  
  
    if od_dir[-1] != '/':
        od_dir += '/'
    
    Nang = len(cos_angles)

    # .. Inputs for doing radiative transfer
    AA = 0.278                 # Pade coefficient
  
    if (type(tave)) == int or (len(tave) == 1):
        tave = (tzl[:len(tzl)-1] + tzl[1:len(tzl)]) / 2
  
    # Set number of layers and preallocate matrices
    if get_upper_atm:
        nlyr = nlyr_toa
        # Preallocate big arrays
        t_hi_mo = np.nan * np.ones((20000000, Nang))
        t_Lm1_hi_mono = np.ones((20000000, Nang))
        rad_hi_mono = np.zeros((20000000, Nang))
        t_Lm1_mono = np.ones((20000000, Nang))
        t_mono = np.ones((20000000, Nang))
        dt_dtau_mono = np.ones((20000000, Nang))
    else:
        nlyr = nlyr_trop
        t_Lm1_mono = np.ones((1000000, Nang))
        t_mono = np.ones((1000000, Nang))
        dt_dtau_mono = np.ones((1000000, Nang))
        nu_mono = []
        rad_hi_mono = []
   
    # .. Viewing angle can be downwelling (0) or upwelling (180)
    if np.all(cos_angles >= 0):
      # View from surface, so surface up
      ilayers = list(range(1, nlyr+1))
      upflag = 0
    elif np.all(cos_angles < 0):
      # View from above, so flip everything
      ilayers = range(nlyr, 1, -1)
      upflag = 1
      raise NameError('Write code for upwelling, based on radtran_belowCld.m')
    else:
      raise NameError('Bad value for cosine of viewing angle' +
                      'They must all be uplooking or downlooking.')
      
    # .. Before we do all this work, make sure the final optical depth
    #    file exists
    fileinp = 'ODdeflt_' + "%03d" % (ilayers[-1])
    if not os.path.isfile(od_dir + fileinp):
        # .. Change this to write an error to a log file and move on
        #    to next file?
        raise NameError('Final OD file: ', fileinp, ' does not exist!')

    # .. Get transmittance, rad., from the bottom layer up
    firstlayer = True
    initialize_low_res = True
    nu = 0.5                       # dummy value
    for ilayer in ilayers:
        fileinp = 'ODdeflt_' + "%03d" % (ilayer)
        
        # .. Number of points in long interferogram
        #    z > 50: 2**23  
        #    19 < z < 50:  2**22  
        #    11 < z < 19:  2**21  
        #    0 < z < 11:  2**20  
        if zs[ilayer] > 50:
            n_long = 2**23
        elif zs[ilayer] > 19:
            n_long = 2**22
        elif zs[ilayer] > 11:
            n_long = 2**21
        elif zs[ilayer] > 0:
            n_long = 2**20
       
        
        # .. Load LBLRTM file
        nu_mono, od_mono = read_od_lblrtm(od_dir + fileinp, prec)
        len_nu_mono = len(nu_mono)
        
        if firstlayer:
            nuLm1_mono = nu_mono
            len_nuLm1_mono = len(nuLm1_mono)
            firstlayer = False

 
        # .. Steps that do not depend on viewing angle
        #    Planck functions
        BB = plancknu(nu_mono, tave[ilayer-1])
        if upflag:
            # Note that tzl is from the bottom up here, because
            # we are going from the TOA down, and
            # tzl(itop) or the temp at TOA, is used first
            BBA = plancknu(nu_mono, tzl[ilayer-1])
            B_Lm1 = plancknu(nu, tave[ilayer-1])
        else:
            BBA = plancknu(nu_mono, tzl[ilayer-1])
            B_Lm1 = plancknu(nu, tave[ilayer-1])

        for iz, cos_angle in enumerate(cos_angles):

            if iz > 0:
                raise ValueError('Debug for multiple angles!')
                
            if (len_nu_mono != len(nuLm1_mono)):
                # interp old stuff up
                t_Lm1_mono[:len_nu_mono,iz] = np.interp(nu_mono,
                    nuLm1_mono, t_Lm1_mono[:len(nuLm1_mono),iz])
                # spl = interpolate.splrep(nuLm1_mono, t_Lm1_mono)
                # t_Lm1_mono = interpolate.splev(nu_mono, spl)
                if any(np.isnan(t_Lm1_mono[:len_nu_mono,iz])):
                    raise NameError("NaNs found in t_Lm1_mono")


            # .. The transmittance up to the bottom of this layer (t_Lm1_mono)
            #    times the trans of this layer, at monochromatic resolution
            try:
                t_mono[:len_nu_mono,iz] = \
                    t_Lm1_mono[:len_nu_mono,iz] * np.exp(-od_mono / cos_angle)
            except:
                raise ValueError('Problem computing transmittance, see code')
            
            dt_dtau_mono[:len_nu_mono,iz] =  \
               (t_mono[:len_nu_mono,iz] - t_Lm1_mono[:len_nu_mono,iz]) \
                 / (od_mono / cos_angle)
    
            # .. Reduce resolution after multiplying by B for the average
            #    temperature of the layer. Note that for ilayer = 1,
            #    tave[ilayer-1=0] corresponds to the average T for layer 1,
            #    extending from z[0] to z[1], while t_mono for ilayer=1
            #    corresponds to the transmittance through the first layer.
            #    Thus here we have tave[layer 1] * t[layer 1]
            nuL, BtransL0 = reduce_res(nu_mono,
                                      BB.T * t_mono[:len_nu_mono,iz],
                                      dnu, N, Npad, n_long)
            dum, Bdt_dtau0 = reduce_res(nu_mono, 
                                       BB.T * dt_dtau_mono[:len_nu_mono,iz],
                                       dnu, N, Npad, n_long)
           
        
            # .. Chop to bwn:ewn.  PMR, 2016/08/31
            nu_mask = np.logical_and(nuL>=bwn, nuL<=ewn)
            BtransL0 = np.real(BtransL0[nu_mask])
            Bdt_dtau0 = np.real(Bdt_dtau0[nu_mask])
            nuL = nuL[nu_mask]
        
            # .. Initialize some array sizes. This happens only once, for
            #    for the first zenith angle only. So values for all 
            #    zenith angles must be initialized here.
            if initialize_low_res:
                nu = nuL
                nu_len = len(nu)
                rads = np.zeros((nu_len, nlyr_trop, Nang))
                tc = np.ones((nu_len, nlyr, Nang))
                dt_dtau = np.ones((nu_len, nlyr, Nang))
                Bc_tc = np.zeros((nu_len, nlyr+1, Nang))
                Bc_tc[:,0,:] = (plancknu(nu, tzl[0])*np.ones((Nang,1))).T
                B_Lm1 = plancknu(nu, tave[ilayer-1])
                radL = np.zeros((nu_len, Nang))
                transL = np.zeros((nu_len, Nang))
                initialize_low_res = False

            transL0 = BtransL0 / B_Lm1
            dt_dtau0 = Bdt_dtau0 / B_Lm1
 
    
            # ... The low resolution wavenumber vector can lose or gain a
            #     single point at the end due to minor differences. If that
            #     happens, we need to interpolate to the standard grid
            if (len(nuL) != nu_len) or (np.sum(nuL - nu)>0) :
                if nuL[0] - nu[0] > 1e-13:
                    raise NameError('nu vector changed!')
                elif len(nuL)+1 == len(nu):
                    transL0 = transL0[:-1]
                else:
                    raise NameError('Not expecting this outcome, see code.')
            
            transL[:,iz] = transL0
            
            # .. Get Planck function of the layer in the same way as LBLRTM.
            #    BBA corresponds to the lower boundary temperature;  
            #    since the boundary indices start from zero while layers
            #    start from 1, the lower bound is ilayer - 1
            #    Similarly, BB corresponds to the average temperature, which
            #    is for index ilayer-1 for layer ilayer
            XX = AA * od_mono / cos_angle
            BL = (BB + XX*BBA) / (1 + XX)
                    
            # .. Surface-to-layer transmittance will go up to TOA (for DISORT)
            tc[:, ilayer-1, iz] = transL[:,iz]          # trans to layer base
            dt_dtau[:, ilayer-1, iz] = dt_dtau0    
            Bc_tc[:,ilayer,iz] = BtransL0
            

            if (ilayer <= nlyr_trop) or upflag:
                # .. Get Bc tsurface-to-cloud at monochromatic resolution,
                #    Then reduce the resolution (slower, more accurate way)
                #    Why transmittance to layer top, temperature from base?
                #    What if we used BL here?
                
                # ... The layer radiance that reaches the surface
                rad0 = -BL * (t_mono[:len_nu_mono,iz] \
                              - t_Lm1_mono[:len_nu_mono,iz])
                
                # ... Reduce resolution
                nuL_0, radL = reduce_res(nu_mono, rad0, dnu, N, Npad, n_long)
                radL = np.interp(nu, nuL_0, np.real(radL))
                
                 
                # .. Assign values for this layer
                #    For rads, subtract one from ilayer because python
                #    counts from zero, while layers start from 1.
                rads[:, ilayer-1, iz] = radL + rads[:,ilayer-2, iz]   
                
                if (ilayer == nlyr_trop):
                    if iz == 0:
                        nu_trop_mono = nu_mono
                        t_trop_mono = np.zeros((len_nu_mono, Nang))
                    t_trop_mono[:,iz] = t_mono[:len_nu_mono,iz]
                    
            elif (ilayer > nlyr_trop):
                if (ilayer > nlyr_trop+1) and \
                    (len_nu_mono != len(nuLm1_mono)):
                    # .. Interpolate up old monochromatic transmittance 
                    #    and rad (the old nu will be the same as for the
                    #    lower atmosphere)
                    t_Lm1_hi_mono[:len_nu_mono, iz] = np.interp(nu_mono,
                      nuLm1_mono, t_Lm1_hi_mono[:len_nuLm1_mono,iz])
                    rad_hi_mono[:len_nu_mono,iz] = np.interp(nu_mono,
                      nuLm1_mono, rad_hi_mono[:len_nuLm1_mono,iz])
    
                # Get the monochromatic radiance from tropopause up
                t_hi_mo[:len_nu_mono, iz] = t_Lm1_hi_mono[:len_nu_mono,iz] \
                  * np.exp(-od_mono / cos_angle)
    
                rad_hi_mono[:len_nu_mono,iz] += BL * \
                  (t_Lm1_hi_mono[:len_nu_mono,iz] - t_hi_mo[:len_nu_mono,iz])
                
                # Set for next time around
                t_Lm1_hi_mono[:len_nu_mono,iz] = t_hi_mo[:len_nu_mono,iz]
        
            else:
                raise NameError('Bad value for ilayer')
        
        
            # .. Set values for next loop through
            t_Lm1_mono[:len_nu_mono,iz] = t_mono[:len_nu_mono,iz]
            
            # Increment index to viewing angle
            iz += 1
            
        # .. new values become old values before next loop through  
        nuLm1_mono = nu_mono
        len_nuLm1_mono = len(nuLm1_mono)


    # .. Scale the results
    rads = 1e3 * rads
    Bc_tc = 1e3 * Bc_tc
  
    # .. Cut extra components of rad_hi_mono and scale
    if nlyr_toa > nlyr_trop:
        rad_hi_mono = 1e3 * rad_hi_mono[:len_nu_mono,:]
  
    # .. nu_trop_mono, t_trop_mono, nu_mono, and rad_hi_mono
    #    only need to go from x to y
    return(nu, rads, tc, dt_dtau, nu_trop_mono, t_trop_mono, 
           nu_mono, rad_hi_mono, Bc_tc)


def plancknu(nu_icm,T):

  # function f = plancknu(nu_icm,T);
  #
  # spectral Planck function as function of wavenumbers (cm-1)
  #
  # [h]    = J*s
  # [c]    = m/s
  # [cbar] = cm/s
  # [k]    = J*K-1
  #
  #
  #    Note: LBLRTM uses ...
  #c    Constants from NIST 01/11/2002
  #
  #      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
  #     *     CLIGHT / 2.99792458E+10 /, 
  # 
 
  import numpy as np
   
  h    = 6.62606896e-34				# J s;  CODATA 2006
  c    = 2.99792458e8				# m/s;  NIST
  k    = 1.3806504e-23				# J K-1; CODATA 2006
  cbar = 100*c			  		# cm/s
  
  nu = nu_icm+0.
  
  # Ignore the possibility of input nu == 0 for now
  #if np.isscalar(nu):
  #  if nu==0:
  #    f = 0
  #else
  #  i0 = np.where(nu==0)[0]
  #  if len(i0) > 0:
  #    nu[i0] = 1e-6


  top    = 2 * h * cbar**3 * nu**3
  bottom =  c**2 * ( np.exp(h*cbar*nu/(k*T))-1 )
  
  f = cbar*top/bottom
 

  # now set f to zero for wavenumbers less than or equal to 0
  #if len(i0)>0:
  #  f[i0] = 0


  #[B] = cm*s-1*J*s*cm3*s-3*cm-3*m-2*s2
  #[B] = W m-2 cm

  return f

	
