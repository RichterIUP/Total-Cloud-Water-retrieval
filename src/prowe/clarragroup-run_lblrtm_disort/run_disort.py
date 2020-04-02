"""
run_disort.py

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose:
  A code for calling disort_driver_py for running DISORT.
"""



# .. Built in modules
import numpy as np
from scipy.interpolate import interp1d

# .. My modules
from disort_driver_py import disort_driver


def run_disort(nu, dtau_gas, CLDLYR,
    reff_liq, total_od_vis_liq,
    reff_ice, total_od_vis_ice,
    iTemp_liq, od_vis_liq, iTemp_ice, od_vis_ice,
    iliq_layer, iice_layer,
    sspLiq, sspIce, TEMPER, BTEMP, NSTR_in, UMU0, UMU,
    albedo_nu, fbeam_nu, delv_nu, iobs):
    """
    Purporse: Run DISORT for multi-layer liquid and ice clouds.
    
    Inputs:
     nu:               wavenumber vector 
     dtau_gas:         layer optical depths of gas
     total_od_vis_liq: Total optical depth in visible for liquid
     total_od_vis_ice: Total optical depth in visible for liquid
     iTemp_liq:        Inds to ssp files that are used, for liquid
     iTemp_ice:        Inds to ssp files that are used, for liquid
     od_vis_liq:       od in visible for liquid, per layer, per ssp
     od_vis_ice:       od in visible for ice, per layer, per ssp set
     reff_liq:         effective radius of liquid (scalar)
     reff_ice:         effective radius of ice (scalar)
     iliq_layer: 
     iice_layer: 
     cldlyr_liq:       vector for number of layers
     cldlyr_ice:       vector, for number of layers
     iobs:             observer height index
    
    Input variable types and sizes:
     nu:               (vector)
     dtau_gas:         (nus x layers)
     total_od_vis_liq: (scalar) 
     total_od_vis_ice: (scalar) 
     iTemp_liq:        (x layers)
     iTemp_ice:        (x layers)
     od_vis_liq:       (layers x ssps)
     od_vis_ice:       (layers x ssps)
     reff_liq:         (scalar)
     reff_ice:         (scalar)
     iliq_layer: 
     iice_layer: 
     cldlyr_liq:       vector for number of layers
     cldlyr_ice:       vector, for number of layers
     iobs:             observer height index
   
    
    
    Disort parameters (asterisks indicate variables that change with nu):
    From DISORT.txt:
      ANGLE CONVENTION:
        Polar (zenith) angles are measured from the upward direction:
        straight up is zero degrees and straight down is 180 degrees.
        There is a small inconsistency in that, for historical reasons,
        the cosine of the incident beam angle (UMU0) is positive,
        whereas according to this convention it should be negative.
      NLYR                     Number layers, INGEGER*
      DTAUC                    Layer optical depths, REAL [1:MAXCLY]*
      SSALB                    Single scattering albedo, REAL [1:MAXCLY]*
      NMOM                     Number of moments, INTEGER*
      NCLDLYR                  Number of cloud layers, INTEGER
      CLDLYR+1                 cloud layers TOA to surf, INTEGER [1:MAXCLY]
      Pmom_cld                 Legendre Poly, REAL [0:MAXMOM,1:MAXCLY]*
      TEMPER = (input)         bndry Ts TOA-to-surface REAL [0:MAXCLY]
      WVNUMLO                  REAL*
      WVNUMHI,                 REAL*
      USRTAU = 1               1=>specify taus, LOGICAL
      NTAU   = 1               Number output levels?  INTEGER
      UTAU   =                 Output levels, as sum(DTAUC), REAL [MAXULV]*
      NSTR   = (input)         Number of streams*
      USRANG = 1               True => specify observation zenith angle
      NUMU   = 1               Number of observation zenith angles
      UMU    = (input)         Cosine of output angles in increasing order:
                               neg (downward) to pos (upward) values, [1:MAXUMU]
      NPHI   = 1               Number of observation azimuthal angles
      PHI    = [0.]            Observation azimuthal angle
      IBCND  = 0               0 => General boundary condition flag
      FBEAM  =                 solar, depends on fbeam_nu*
      UMU0   = (input)         Cosine of incident beam angle, REAL
      PHI0   = 0.              Azimuth angle of incident beam
      FISOT  = 0.              Isotropic incident flux
      LAMBER = 1               True => The surface is Lambertian
      ALBEDO =                 ALBEDO_long[nu], from input, REAL*
      BTEMP                    Bottom T (e.g. surface T), REAL
      TTEMP =                  TEMPER[i1], Top T (NA if TEMIS = 0), REAL*
      TEMIS  = 0.              Top emissivity, REAL
      PLANK  = 1               True => thermal emission included, INTEGER
      ONLYFL = 0               False => Return azimuthally averaged
                               intensities (not just fluxes), INTEGER
      HEADER = '...'           unused string, CHAR STRING
      MAXCLY =                 len(DTAUC), INTEGER
      MAXULV =                 len(UTAU), number output levels, INTEGER
      MAXUMU =                 len(UMU), number UMU?, INTEGER
      MAXPHI =                 len(PHI), number PHI, INTEGER
      
      So we should have:
      dtauc_nu
      ssalb_nu
      Pmom_cld_nu
      utau_nu
      albedo_nu (input)
      fbeam_nu (input)
    """
    
    # # # # # # # # #     Debugging     # # # # # # # # # # #
    # .. If the following flag is set to True, the disort inputs
    #    will be printed to a file for each wavenumber bin
    #    to allow debugging. 
    #
    #    Note: do not leave this set to True if not debugging, as it
    #    will slow the code down!
    debug_flag = False          # If True, will print DISORT inputs to log    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    
    
    # # # # # # # # #     Set parameters     # # # # # # # # # # #
    
    # .. For each nu and Temp, there can be a different number of Pmoms
    #    So we need to set the maximum allowed and also set the maximum
    #    ever used (starting with 0, here)
    disort_warning = False      # Will become true if DISORT returns a message
    maxPmom = 500               # Maximum number of Pmoms allowed
    
    # .. Number of layers and frequencies
    NLYR_long, nnus = dtau_gas.shape

    # .. DISORT inputs that are hard-wired 
    PLANK  = 1              # True => thermal emission included
    USRTAU = 1              # True => we want to specify taus
    NTAU   = 1              # There will be one tau
    USRANG = 1              # True => specify observation zenith angle
    NUMU   = 1              # Number of observation zenith angles
    NPHI   = 1              # Number of observation azimuthal angles
    PHI    = [0.]           # Observation azimuthal angle
    IBCND  = 0              # 0 => General boundary condition flag
    PHI0   = 0.             # Azimuth angle of incident beam
    LAMBER = 1              # True => The surface is Lambertian
    ONLYFL = 0              # False => Return azimuthally averaged intensities
                            # (not just fluxes)
    FISOT  = 0.             # Isotropic incident flux
    TEMIS  = 0.             # Top emissivity
    HEADER = 'Run DISORT from Python'  # char string, not used by disort.disort
    MAXPHI = 1
    MAXULV = 1
    MAXUMU = len(UMU)

    # BTEMP is now an input variable     PMR 2020/03/03
    # BTEMP  = TEMPER[NLYR_long]   # Bottom temperature. 
    
    
    # .. Number of temperatures to interpolate SSPs over
    iLiq = np.unique(iTemp_liq)
    iIce = np.unique(iTemp_ice)
    Nssp_sets_liq = len(iLiq)
    Nssp_sets_ice = len(iIce)
    
    # .. Remove zero-d out layers
    if total_od_vis_liq <= 0:
        if total_od_vis_ice <= 0:
            Nssp_sets_liq = 0
            Nssp_sets_ice = 0
            maxPmom = 5
            NSTR_in = 4
            # We have to have at least one cloud layer, put it at the surface
            CLDLYR = np.array([NLYR_long-1])
            NCLDLYR = 1
        else:
            Nssp_sets_liq = 0
            CLDLYR = np.array(CLDLYR)[iice_layer]
            od_vis_ice = od_vis_ice[iice_layer,:]
    elif total_od_vis_ice <= 0:
        Nssp_sets_ice = 0
        CLDLYR = np.array(CLDLYR)[iliq_layer]
        od_vis_liq = od_vis_liq[iliq_layer,:]
        
         
    # .. Preallocate
    #    The array cldyTau_* holds the cloud od (in visible region) for
    #    liquid and ice for each possible ssp index and each layer
    #    so it is large and mostly full of zeros
    Radiance = np.zeros(nnus)
    
    # .. Preallocate vectors and matrices that depend on number cloud layers
    NCLDLYR = len(CLDLYR)
    SSALB = np.zeros(NLYR_long)
    w0dtau = np.zeros((NCLDLYR, nnus))
    Pmom_w0dtau = np.zeros((NCLDLYR, nnus, maxPmom))
    Pmom = np.zeros((maxPmom, nnus, NCLDLYR))
    dtau = 1*dtau_gas
    
    
    # .. Loop over the sspfiles we are using (depends on temps and layers)
    #    Get liquid cloud component based on weights to ssp files
    #    There may be multiple liquid cloud layers, and for each
    #    there can be one or two ssp files; so interpolate
    #    Interpolate liquid optical properties for frequencies & reff
    # NPmom_max = NSTR  # Number of Pmoms must be at least NSTR
    NPmom_nu_max = 5*np.ones(nnus).astype('int32')   # Require at least 5 Pmoms
    for i in range(Nssp_sets_liq):
        j = iLiq[i]                                     # index to ssp set
        current_ssp = sspLiq[j]
        qext = current_ssp.fqext(reff_liq)
        w0_liq = current_ssp.fw0(reff_liq)
        NPmom_nu_i = (np.floor(current_ssp.fNPmom(reff_liq))).astype('int32')
        NPmom_nu_max = (np.max([NPmom_nu_max, NPmom_nu_i], axis=0)).astype('int32')
        # NPmom_max_i = np.max(NPmom_nu_i)                # scalar
        # NPmom_max = np.max([NPmom_max, NPmom_max_i])    # scalar
        
        # Loop over wavenumbers?
        for inu in range(nnus):
            f = interp1d(current_ssp.reff, current_ssp.Pmom[:,inu,:NPmom_nu_i[inu]].T)
            Pmom0 = f(reff_liq)   
            # Pmom0 = sspLiq[j].fPmom[inu](reff_liq)[:NPmom_nu_i[inu]]
            dtau_liq_i_inu = qext[inu] * od_vis_liq[:,j]/2
            w0dtau0_liq_i_inu = w0_liq[inu] * dtau_liq_i_inu
            dtau[CLDLYR,inu] += dtau_liq_i_inu
            w0dtau[:,inu] += w0dtau0_liq_i_inu
            
            for icld in range(len(CLDLYR)): 
                Pmom_w0dtau[icld, inu, :NPmom_nu_i[inu]] += \
                  w0dtau0_liq_i_inu[icld] * Pmom0
                #if np.any(np.isnan(Pmom_w0dtau)):
                #    raise ValueError('Nan in Pmom')
                
                # .. If no ice, then 
                #    on last ssp, to get Pmom, divide by (dtau*SSALB)
                if (Nssp_sets_ice==0) and (i == Nssp_sets_liq-1):
                    Pmom[:,inu,icld] = Pmom_w0dtau[icld,inu,:]/w0dtau[icld,inu]
                    #if np.any(np.isnan(Pmom)):
                    #    raise ValueError('Nan in Pmom')
        
    for i in range(Nssp_sets_ice):
        j = iIce[i]                                     # index to ssp set
        current_ssp = sspIce[j]
        qext = current_ssp.fqext(reff_ice)
        w0_ice = current_ssp.fw0(reff_ice)
        NPmom_nu_i = (np.floor(current_ssp.fNPmom(reff_ice))).astype('int32')
        NPmom_nu_max = (np.max([NPmom_nu_max, NPmom_nu_i], axis=0)).astype('int32')
        # NPmom_max_i = np.max(NPmom_nu_i)
        # NPmom_max = np.max([NPmom_max, NPmom_max_i])  # scalar
        
        # Loop over wavenumbers
        for inu in range(nnus):
            f = interp1d(current_ssp.reff, current_ssp.Pmom[:,inu,:NPmom_nu_i[inu]].T)
            Pmom0 = f(reff_ice)
            # Pmom0 = sspIce[j].fPmom[inu](reff_ice)[:NPmom_nu_i[inu]]
            dtau_ice_i_inu = qext[inu] * od_vis_ice[:,j]/2
            w0dtau0_ice_i_inu = w0_ice[inu] * dtau_ice_i_inu            
            dtau[CLDLYR,inu] += dtau_ice_i_inu
            w0dtau[:,inu] += w0dtau0_ice_i_inu
            
            
            for icld in range(len(CLDLYR)):
                Pmom_w0dtau[icld, inu, :NPmom_nu_i[inu]] += \
                  w0dtau0_ice_i_inu[icld] * Pmom0 #[:NPmom_nu_i[inu]]
                #if np.any(np.isnan(Pmom_w0dtau)):
                #    raise ValueError('Nan in Pmom')

                
                # On last ssp, to get Pmom, divide by (dtau*SSALB)
                if i == Nssp_sets_ice-1:
                    Pmom[:,inu,icld] = Pmom_w0dtau[icld,inu,:]/w0dtau[icld,inu]
                    #if np.any(np.isnan(Pmom)):
                    #    raise ValueError('Nan in Pmom')
    
    
    # .. To get SSALB, divide by dtau,
    w0 = w0dtau / dtau[CLDLYR,:]
    
    # .. Solar contribution only if umu0 > 0, otherwise
    #    sun is below the horizon, so turn it off
    if UMU0 <= 0:
        FBEAM_nu = 0*nu
        UMU0 = 0
    elif UMU0 > 0:
        FBEAM_nu = fbeam_nu
    else:
        raise NameError('Bad value for UMU0')
    
    
    # .. Wavenumbers
    WVNUMLOs = nu - delv_nu/2
    WVNUMHIs = nu + delv_nu/2
    
    
    # .. Make sure there are at least as many Pmoms as strings
    NSTR = 1 * NSTR_in
    #if np.any(NPmom_nu_max < NSTR+1):
    #    print('got one')
    NPmom_nu_max[NPmom_nu_max<NSTR+1] = NSTR+1
    
    # ..  Loop over the frequencies
    izm = np.zeros((nnus,1))
    for inu in range(nnus):
        uu = np.nan # This logic because disort doesn't always converge
    
        # Values for this wavenumber
        WVNUMLO = WVNUMLOs[inu]
        WVNUMHI = WVNUMHIs[inu]
        FBEAM = FBEAM_nu[inu]
        DTAUC = dtau[:,inu]                             # x NLYR
        SSALB[CLDLYR] = w0[:,inu]                       # x NCLDLYR
        PMOM_CLD = Pmom[:NPmom_nu_max[inu],inu,:]       # NMOM x CLDLYR
        
        #if np.any(np.isnan(PMOM_CLD)):
        #    raise ValueError('Nan in Pmom')
        
        # .. Cut out layers at top where od<5e-7 or 1e-5,
        #    otherwise round-off error can become significant
        i2 = len(DTAUC)
        
        # .. For an observer layer very close to the surface (within a few m),
        #    The optical depth in-between may be < 1e-5 because the layer
        #    is so thin. Furthermore, the layer between the surface and the
        #    observer plays a minior role for IR downwelling radiance,
        #    as it is only impacts radiation scattered of the surface and
        #    then scattered back down again to the observer. So, if the
        #    optical depth of this layer is too small for single precision
        #    accuracy, it's best to just remove it.
        #
        #    Here we check if:
        #    1) the observer layer is specified: iobs > 0, 
        #    2) the observer layer is above the: iobs <= len(DTAUC),
        #    3) The optical depths for layers below iobs to the surface layer
        #     are of the order of the precision: (np.sum(DTAUC[iobs:i2])<1e-5)
        # 
        #    If 1-3 all true, remove the surface layers by setting i2 = iobs
        if (iobs > 0) and (iobs < len(DTAUC)) \
          and (np.sum(DTAUC[iobs:i2])<1e-5):
            i2 = iobs

        # .. Some of the upper layers may also have optical depths that
        #    are on the order of the precision. Remove those by
        #    increasing i1 from 0.
        izeds = np.where(DTAUC[:i2]<1e-5)[0]
        if len(izeds)==0:
          i1 = 0
        elif (izeds[0] != 0) or (np.any(np.diff(izeds)!=1)):
          # print('Warning: Redesign your layers so od gets smaller going up.');
          i1 = np.max(izeds)+1
        else:
          i1 = np.max(izeds)+1
        
        if i1 == i2:
          raise NameError('i1=i2= ' + str(i1) + ', so no atmosphere. ' \
                          'This is likely because too many optical ' \
                          'depths are <=1e-6. \n Near-surface optical ' \
                          'depths are ' \
                           + str(DTAUC[-1]) + ' ' \
                           + str(DTAUC[-2]) + ' ' \
                           + str(DTAUC[-3]) + ' ' \
                           + str(DTAUC[-4]) + ' ' \
                           + '.\n nu = ' + str(nu[inu]))
        izm[inu] = i1
        NLYR = i2-i1
        try:
            CLDLYR_inu = np.array(CLDLYR) - i1 + 1
        except:
            raise ValueError('This should not happen')
            CLDLYR_inu = np.array(CLDLYR-i1) + 1
            
        # Height of observer set with index iobs (index to layer),
        # iobs = len(DTAUC) => surface, cumulative od of atmosphere
        # iobs = 0 => TOA, UTAU = 0
        if (iobs > 0) and (iobs <= len(DTAUC)):
            UTAU = [np.sum(DTAUC[i1:iobs])]
        elif iobs == 0.:
          UTAU =[0.]
        else:
          raise NameError('Option iobs = ' + str(iobs) + ' not allowed.')
    
        # Old maxes not used, b/c sizes must be exact
        # e.g. MAXCLY = 120; MAXPHI = 1; MAXULV = 2; MAXUMU = 10
        MAXCLY = NLYR
        MAXMOM = PMOM_CLD.shape[0]-1
        NMOM = MAXMOM
        #NSTR = 1 * NSTR_in
        #if NSTR > NMOM:
        #    if NMOM % 2 == 0:
        #        NSTR = NMOM
        #    else:
        #        NSTR = NMOM - 1
        #    if NSTR<=3:
        #        raise NameError('NSTR  is <=3')
        # .. We want at least 16 streams
        if NSTR < 16:
            msg = 'Wavenumber: ' + str(WVNUMLO) + ', Re, liq: ' \
                   + str(reff_liq)  + ', Re, ice: ' +str(reff_ice) \
                   + ', streams: ' + str(NSTR)
            raise ValueError(msg)
    
        if debug_flag: # and inu==0:
            dlog = open("disort_output.txt", "w")
            print("WVNUMLO, WVNUMHI", WVNUMLO,WVNUMHI, file = dlog)
            print("UMU0 = ", UMU0, file = dlog)
            print("MAXULV, MAXUMU, MAXPHI = ",MAXULV,MAXUMU,MAXPHI,file=dlog)
            print("PHIP, FISOT, LAMBER = ", PHI0, FISOT, LAMBER, file=dlog)
            print("TEMIS, PLANK, ONLYFL = ", TEMIS,PLANK,ONLYFL, file = dlog)
            print("HEADER = ", HEADER, file = dlog)
            print("ALBEDO = ", albedo_nu[inu], file = dlog)
            print("FBEAM = ", FBEAM, file = dlog)
            
            print("NLYR = ", NLYR, file = dlog)
            print("MAXCLY = ", MAXCLY, file = dlog)
            print("TEMPER.shape = ", TEMPER[i1:i2+1].shape, file = dlog)     
            print("TEMPER = ", TEMPER[i1:i2+1], file = dlog)   
            print("TTEMP = ", TEMPER[i1], file = dlog)
            print("BTEMP =", BTEMP, file = dlog)
            
            print("USRTAU,NTAU,UTAU = ", USRTAU,NTAU,UTAU, file = dlog)
            print("USRANG,NUMU,UMU = ", USRANG,NUMU,UMU, file = dlog)
            print("NPHI,PHI,IBCND = ", NPHI,PHI,IBCND, file = dlog)
            print("DTAU = ", DTAUC[i1:i2], file = dlog)
            print("DTAU.shape = ", DTAUC[i1:i2].shape, file = dlog)
            print("SSALB = ", SSALB[i1:i2], file = dlog)
            print("SSALB.shape = ", SSALB[i1:i2].shape, file = dlog)
            print("NCLDLYR = ", NCLDLYR, file = dlog)
            print("Pmom.shape[1] = ", PMOM_CLD.shape[1], file = dlog)
            print("CLDLYR.shape = ", CLDLYR_inu.shape, file = dlog)                
            print("CLDLYR = ", CLDLYR_inu, file = dlog)
            
            print("NSTR =", NSTR, file = dlog)
            print("NMOM = ", NMOM, file = dlog)
            print("MAXMOM = ", MAXMOM, file = dlog)
            print("Pmom.shape[0] = ", PMOM_CLD.shape[0], file = dlog)
            print("Pmom = ", PMOM_CLD, file = dlog) 
            dlog.close()
            
            # Quality control
            if np.any(PMOM_CLD==np.nan) or np.any(SSALB[i1:i2]==np.nan) \
              or np.any(DTAUC[i1:i2]==np.nan) or np.any(DTAUC[i1:i2]<=0):
                  raise ValueError('Bad inputs to DISORT!')

        '''
        # .. Catch for the extremely unlikely but happened occurrence
        #    that UMU0 and NSTR cause a singularity in DISORT
        if (round(UMU0,5)==0.23724) and (NSTR==16):
            if NMOM >= 19: 
                NSTR = 18
            else:
                format_str = "UMU0 = %.5f, NSTR = %i, NMOM = %i"
                print(format_str % (UMU0, NSTR, NMOM))
        elif (round(UMU0,5)==0.12923) and (NSTR==14):
            if (NMOM==14) or (NMOM==15):
                NSTR = 12
            else:
                format_str = "UMU0 = %.5f, NSTR = %i, NMOM = %i"
                print(format_str % (UMU0, NSTR, NMOM))
        elif (round(UMU0,5)==0.23078) and (NSTR==10):
            if (NMOM <= 11):
                NSTR = 12
                NMOM = 12
                MAXMOM = 12
                PMOM_CLD = Pmom[:NMOM+1,inu,:]
                NMOM = MAXMOM
            else:
                format_str = "UMU0 = %.5f, NSTR = %i, NMOM = %i"
                print(format_str % (UMU0, NSTR, NMOM))
        '''
            
        while np.isnan(uu):
            # Call disort
            # We must add 1 to CLDLYR below because python indexes from
            # zero but in the fortran code it is indexed from one
            # PMR, 2016/08/28
            rfldir, rfldn, flup, \
            dfdt, uavg, uu, \
            albmed, trnmed, \
            errmsg, errflag = disort_driver(NLYR, DTAUC[i1:i2], SSALB[i1:i2], \
                                            NMOM, NCLDLYR, CLDLYR_inu, \
                                            PMOM_CLD, TEMPER[i1:i2+1], \
                                            WVNUMLO, WVNUMHI, USRTAU, NTAU, \
                                            UTAU, NSTR, USRANG, NUMU, UMU, \
                                            NPHI, PHI, IBCND, FBEAM, UMU0, \
                                            PHI0, FISOT, LAMBER, \
                                            albedo_nu[inu], BTEMP, \
                                            TEMPER[i1], TEMIS, PLANK, \
                                            ONLYFL, HEADER, MAXCLY, \
                                            MAXULV, MAXUMU, MAXPHI, MAXMOM)
            
            # .. If DISORT returns a fatal error (indicated by errflag), quit
            if errflag:
                raise ValueError(errmsg)


        #    If DISORT returns a warning, print the wavenumber and warning
        #    to the log file and set output variable disort_warning to true            
        errmsg = (errmsg.decode("utf-8")).rstrip()
        if errmsg != '':   
            disort_warning = True
            log.warning(str(round(WVNUMLO,2)) + ' cm-1: ' + errmsg)
            
                 
        Radiance[inu] = uu  

    Radiance = Radiance*1e3 # Radiance units (milliwatts)
    

    return Radiance, izm, rfldn, flup, disort_warning

'''
if disort_warnings:
    log.error(thisdatetime.strftime('%Y%m%d.%H%M%S') + ': ' + msg )
    log.warning(()
'''
    
    
    
