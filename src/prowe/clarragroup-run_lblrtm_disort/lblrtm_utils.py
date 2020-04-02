#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 14:00:20 2019

Copyright 2018 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

@author: prowe

Purpose: Utilities for running LBLRTM using Python
"""


# .. Built-in Modules
import numpy as np
import array
from netCDF4 import Dataset
from scipy.interpolate import CubicSpline
from meteo import humidRS80
from meteo import hypsometric

# # # # # # #     TAPE5 variables     # # # # # # #
class Lbl():
    __slots__ = (
    'model', 
    'v1',
    'v2',
    'zangle',
    'ihirac', 
    'ilblf4', 
    'icntnm',
    'iaersl', 
    'iscan', 
    'ifiltr', 
    'iplt', 
    'itest', 
    'iatm', 
    'ilas', 
    'iod', 
    'mpts', 
    'sample',  
    'dvset', 
    'itype',
    'nozero', 
    'nprnt',  
    'ipunch', 
    're',  
    'hspace',  
    'vbar', 
    'iprfl',                                       
    'ixsbin',                                     
    'iprfl',                                      
    'izorp', 'zorp',
    'xsname','xtitle',
    'iemit','imrg', 'nmol', 'ixsect',
    'ixmols', 'immax', 'ibmax',
    'h1', 'h2', 'izbnd', 'layx', 'zbnd',
    'zm', 'pm', 'tm', 'traceGases',
    'jcharp', 'jchart', 'jchar', 'jcharXs',
    'cfc', 
    )
    
    
    def __init__(self, v1, v2, angle, desired_output):
        
        # .. Input variables
        self.v1 = v1
        self.v2 = v2
        self.zangle = angle
        if desired_output == 'radiance':
            self.iemit = 1
            self.imrg = 0
        elif (desired_output == 'optical depth') or (desired_output == 'od'):
            self.iemit = 0
            self.imrg = 1

        # .. Hard-wired values, integers
        self.model = 0
        self.ihirac = 1
        self.ilblf4 = 1
        self.icntnm = 6
        self.iaersl = 0
        self.iscan = 0
        self.ifiltr = 0
        self.iplt = 0
        self.itest = 0
        self.iatm = 1
        self.ilas = 0
        self.iod = 0
        self.mpts = 0
        self.sample = 0 
        self.dvset = 0
        self.itype = 2
        self.nozero = 0
        self.nprnt = 0 
        self.ipunch = 1
        self.re = 0 
        self.hspace = 0 
        self.vbar = 0
        self.iprfl = 0                                      
        self.ixsbin = 0                                    
        self.iprfl = 0                                     
        self.izorp = 0       # 0=> alt, 1=>P
        
        
        # .. Hard-wired strings. These aren't used anyway
        self.xsname = 'F11       F12       F113      '     # string
        self.xtitle = '  CFC amounts'                      # string
        
        
        # .. Default values for variables that may be changed based on nc file
        self.nmol = 0         # number of molecules that will be specified
        self.ixsect = 0       # flag to indicate inclusion of 'x' or CFC molecules
        self.ixmols = 0       # number of 'x' or CFC molecules
        self.immax = 0        # number of atmospheric profile boundaries
        self.ibmax = 0        # Layering. 0=>LBLRTM generates layers, <0=inputs are P,
                              # otherwise same as immax.
                              
        # .. Input variables, not needed
        # self.iemit = 0        # 1 to output radiance, 0 for optical depth
        # self.imrg = 1         # 0 for radiance, 1 for optical depth
                         
        # .. Optional variables that will be set below if needed
        # jcharp
        # jchart
        # jchar
        # jcharXs
        # h1 = zm[0]                                   # height 1
        # h2 = zm[-1]                                  # height 2
        # zbnd = zm                                    # height bounds
        # layx = immax                                 # number of layers for CFCs
        # zorp = zm                                    # bounds for CFCs
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def set_from_prof(self, prof):
        self.zm = prof.zm
        self.tm = prof.tm
        self.pm = prof.pm
        self.jcharp = prof.jcharp
        self.jchart = prof.jchart
        self.nmol = prof.nmol
        self.traceGases = prof.traceGases
        self.jchar = prof.jchar
        self.cfc = prof.cfc
        self.jcharXs = prof.jcharXs
        self.ixmols = prof.ixmols
        self.ixsect = prof.ixsect

        self.immax = len(prof.zm)
        self.ibmax = len(prof.zm)                   # number of levels
        self.h1 = prof.zm[0]                        # height 1
        self.h2 = prof.zm[-1]                       # height 2
        self.zbnd = prof.zm                         # height bounds
        self.layx = len(prof.zm)                    # number layers for CFCs
        self.zorp = prof.zm                         # bounds for CFCs
        self.zbnd = prof.zm
        

# # # # # # #     profile     # # # # # # #
class Prof():
    __slots__ = ('modelExtra',
                 'zm',
                 'pm',
                 'tm',
                 'jcharp',
                 'jchart',
                 'nmol',
                 'traceGases',
                 'jchar',
                 'cfc',
                 'jcharXs',
                 'ixmols',
                 'ixsect',
                 'profile_variables',
                 )
    
    def __init__(self, prof_file):    
        '''
        Load values from netcdf file prof_file
        '''
    
        # .. Parameters used here
        all_gases = ('h2o', 'co2', 'o3', 'n2o', 'co', 'ch4', 'o2', 'no', \
                     'so2', 'no2', 'nh3', 'hno3', 'oh', 'hf', 'hcl', 'hbr')
        all_cfcs = ('f11', 'f12', 'f113')
        ngases = len(all_gases)
        unit_dict = dict({('K', 'A'), ('mb','A'), ('ppmv','A'), ('gm_kg','C'),\
                          ('B', 'B')})
    
        # .. Profile variables that are functions of height
        self.profile_variables = ('zm', 'pm', 'tm',
                                  'f11', 'f12', 'f113', 
                                  'traceGases', 'cfc')

    
        with Dataset(prof_file, 'r', format='NETCDF4_CLASSIC') as nc:
            # .. Get values from nc file. 
            self.modelExtra = nc['model_extra'][0]
            self.zm = nc['z'][:]
            self.pm = nc['P'][:]
            self.tm = nc['T'][:]        
            self.jcharp = unit_dict[nc['P'].units]
            self.jchart = unit_dict[nc['T'].units]
            
            nlevels = len(self.zm)
            traceGases = np.zeros((nlevels, ngases))   
            s = [str(self.modelExtra) for i in range(ngases)]
            jchar = ''.join(s)
            for gas in all_gases:
                if gas in nc.variables:
                    i = all_gases.index(gas)
                    traceGases[:,i] = nc[gas][:]
                    unit = nc[gas].units
                    if i+1 <= ngases: 
                        jchar_fin = jchar[i+1:]
                    else: 
                        jchar_fin = ''
                    jchar = jchar[:i] + unit_dict[unit] + jchar_fin
                    
                    
            # .. set nmol based on gases present, or nmol if present,
            #    whichever is larger. Gases not specified will be set
            #    to jchar
            self.nmol = i+1
            if 'nmol' in nc.variables:
                self.nmol = max(self.nmol, int(nc['nmol'][:]))
            self.traceGases = traceGases[:,:self.nmol]
            self.jchar = jchar[:self.nmol]
        
            cfc = np.zeros((nlevels, len(all_cfcs)))
            jcharXs = ''.join('A' for i in range(len(all_cfcs)))
            ixmols = 0
            ixsect = 0
            for gas in all_cfcs:
                if gas in nc.variables:
                    i = all_cfcs.index(gas)
                    cfc[:,i] = nc[gas][:]
                    unit = nc[gas].units
                    jcharXs = jcharXs[:i] + unit_dict[unit]
                    ixmols += 1
                    ixsect = 1
            self.cfc = cfc[:,:ixmols]
            self.jcharXs = jcharXs
            self.ixmols = ixmols
            self.ixsect = ixsect
            

    def crop_heights(self, z1, z2):
        '''
        Set the atmospheric profile from zbot to ztop, inclusive,
          as the weighted mean of two profiles
        '''
        
        i = np.logical_and(self.zm>=z1, self.zm<=z2)

       
        # .. Crop any of the variables that exist in the class instance prof
        for att in self.profile_variables:
            if hasattr(self, att):
                v1 = getattr(self, att)
                setattr(self, att, v1[i])


    def include_observer_alt(self, alt_km):
        '''
        Make sure observer altitude is also included in profile
        '''
        
        # .. Quality control on observer altitude
        if alt_km < self.zm[0]:
            raise ValueError('The observer seems to be below the surface.')



        # .. Interpolate to the observer altitude, if needed
        if alt_km not in self.zm:
            old_zm = self.zm
            iz = np.where(self.zm > alt_km)[0][0]
            self.zm = np.hstack([self.zm[:iz], alt_km, self.zm[iz:]])
            

            # .. Interpolate any of the variables that exist in the class 
            #    instance prof
            for att in self.profile_variables:
                if att != 'zm' and hasattr(self, att):
                    old_val = getattr(self, att)
                    if len(np.shape(old_val)) == 1:
                        if att != 'pm':
                            # Linearly interpolate
                            val_at_alt = np.interp(alt_km, old_zm, old_val)
                        
                            new_val = np.hstack([old_val[:iz], 
                                                 val_at_alt, old_val[iz:]])
                            setattr(self, att, new_val)
                        #elif att == 'pm':
                        #    # pressure, use spline
                        #    Pfun = CubicSpline(old_zm, old_val)
                        #    Pspline = Pfun(alt_km)
                    elif len(np.shape(old_val)) == 2:
                        new_val = np.zeros((len(self.zm), np.shape(old_val)[1]))
                        for icol in range(np.shape(old_val)[1]):
                            val_at_alt = np.interp(alt_km, 
                                                   old_zm, old_val[:,icol])
                            new_val[:,icol] = np.hstack([old_val[:iz, icol], 
                                                         val_at_alt, 
                                                         old_val[iz:, icol]])
                        setattr(self, att, new_val)
                    
            # .. Set the pressure last
            rhi = humidRS80(self.traceGases[:2,0], self.tm[:2], 'ppmv', 'rh', self.pm[:2]); 
            Pi = hypsometric(self.zm[:2]*1000, self.tm[:2], rhi[:2], self.pm[0])[1]
            #rhi = humidRS80(self.traceGases[:2,0], self.tm[:2], 'ppmv', 'rh', self.pm[:2]); print('Pi=',Pi, 'RH=',rhi)
            #Pi = hypsometric(self.zm[:2]*1000, self.tm[:2], rhi[:2], self.pm[0])
            #Pi = Pi[1]; print('Pi=', Pi)                 
            self.pm = np.hstack([self.pm[:iz], Pi, self.pm[iz:]])
                    
            

def write_tape5_prof(prof, v1, v2, angle, desired_output, direc):
    '''    
    Write TAPE5 file directly from profile data
    '''
    lbl = Lbl(v1, v2, angle, desired_output)     # Initialize LBLRTM variables
    lbl.set_from_prof(prof)                      # Set vars from profile
    write_tape5(lbl, direc + 'TAPE5')            # Write the TAPE5 file





def write_tape5(lbl, fnametape5):
    '''
    Write a TAPE5 file for input to LBLRTM
    Warning: This has not been extensively tested
    Inputs:
        lbl: class containing variables needed to write TAPE5 file
             for the limited LBLRTM cases used here.

    '''
    f = open(fnametape5, 'w')
    
    
    # .. Header
    hr1 = '         1         2         3         4         5'
    hr2 = '         6         7         8'
    hr3 = '123456789012345678901234567890123456789012345678901234567890123456789'
    hr4 = '01234567890'
    f.write(hr1 + hr2 + '\n')
    f.write(hr3 + hr4 + '\n')
    f.write('$'+'\n')
      
    
    # .. Record 1.2
    s = ' HI=' + str(lbl.ihirac) + ' F4=' + str(lbl.ilblf4) + \
        ' CN='+str(lbl.icntnm) + ' AE=' + str(lbl.iaersl) + \
        ' EM=' + str(lbl.iemit) + ' SC='+str(lbl.iscan) + \
        ' FI=' + str(lbl.ifiltr) + ' PL=' +str(lbl.iplt) + \
        ' TS=' + str(lbl.itest) + ' AM=' + str(lbl.iatm) + \
        ' MG=' + str(lbl.imrg) + ' LA=' +str(lbl.ilas) + \
        ' MS=' + str(lbl.iod) + ' XS=' + str(lbl.ixsect) + \
        '   0' + str(lbl.mpts) +'   00'
    f.write(s + '\n')
    
    
    # .. Record 1.2a (required if icntnm==6, otherwise omit)
    f.write('  1.000  1.000  1.000  1.000  1.000  1.000  1.000\n')
    
    
    # .. Record 1.3 (required if IHIRAC>0; IAERSL>0; IEMIT=1; or ILAS>0)
    #    (otherwise omit)
    s = '%10.3f%10.3f%10.3f%10.3f' % (lbl.v1, lbl.v2, lbl.sample, lbl.dvset)
    f.write(s + '                                   ')
    f.write('                            \n')
    
    
    # .. Record 1.4 (required if IEMIT=1, or both IEMIT=2 and IOTFLG=2)
    #    (otherwise omit)
    if lbl.iemit == 1:
        f.write('0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00' \
                + '0.000e+00 0.000e+00\n')     
        # s = '%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e ' % \
        #      tbound, sremis(1),sremis(2),sremis(3), 
        #      srrefl(1), srrefl(2), srrefl(3)
        # s2 = 
        # f.write(s + '    ' + surf_refl + '\n')    
    elif lbl.iemit != 0:
        raise NameError('iemit must be zero or 1')
    
    
    
    # .. Record 3.1
    if lbl.iatm == 1:
        # ifxtype, munits left out
        s1 = '%5i%5i%5i%5i%5i%5i%5i' % \
            (lbl.model, lbl.itype, lbl.ibmax, lbl.nozero, lbl.nprnt, \
             lbl.nmol, lbl.ipunch)
        s2 = '%10.3f%10.3f%10.3f' % (lbl.re, lbl.hspace, lbl.vbar)
        f.write(s1 + '     ' + s2 + '\n')
    else:
        raise NameError('iatm must be 1')
        
        
    # .. Record 3.2 for non-horizontal path
    #	 required if itype = 2,3, slant path
    if lbl.model >= 0 and lbl.model <= 6:
        if lbl.itype == 2:
            f.write('%10.3f%10.3f%10.3f\n' % (lbl.h1, lbl.h2, lbl.zangle))
        else:
            raise NameError('Option not here yet, see code')
        
        if lbl.ibmax == 0:
            raise NameError('ibmax==0 not here yet, see code')
        elif lbl.ibmax > 0:
            count = 0
            for k in range(int(np.ceil(lbl.ibmax/8))):
                if count < lbl.ibmax-7:
                    s = ['%10.3f' % z for z in lbl.zbnd[count:count+8]]
                    f.write(''.join(s))
                else:
                    s = ['%10.3f' % z for z in lbl.zbnd[count:lbl.ibmax]]
                    f.write(''.join(s))
                f.write('\n')
                count += 8
                        
        else:
            raise NameError('Bad value for ibmax')
            
    
        # .. Records 3.4, 3.5, for model=0, user atmospheres ..
        #	 required for model==0, otherwise omit
        if lbl.model == 0:
            if lbl.immax >= 1:
                i1 = min(lbl.nmol, 8)
                i2 = min(lbl.nmol, 16)
                i3 = min(lbl.nmol, 24)
                f.write('%5i         user atmosphere\n' % lbl.immax)
                for k in range(lbl.immax):
                    f.write('%10.3f%10.3f%10.3f' % (lbl.zm[k], lbl.pm[k], lbl.tm[k]))
                    f.write('     ' + lbl.jcharp + lbl.jchart + '   ')
                    f.write(lbl.jchar + '\n')
                    s = ['%10.3e' % gas for gas in lbl.traceGases[k,:i1]]
                    f.write(''.join(s) + '\n')
                    if lbl.nmol > 8:
                        s = ['%10.3e' % gas for gas in lbl.traceGases[k,8:i2]]
                        f.write(''.join(s) + '\n')
                    if lbl.nmol > 16:
                        s = ['%10.3e' % gas for gas in lbl.traceGases[k,16:i3]]
                        f.write(''.join(s) + '\n')
            else:
                raise NameError('immax must be >= 1')
        
        
        else:
            raise NameError('model must be 0 - 6, only model==0 available')
            
        
        if lbl.iatm==1 and lbl.ixsect==1:
            # .. write record 3.7
            f.write('%5i%5i%5i\n' % (lbl.ixmols, lbl.iprfl, lbl.ixsbin)) 
    	
    		# write record 3.7.1
            f.write('%s\n' % lbl.xsname)
            
            if lbl.iprfl == 0:
    			# write record 3.8
                f.write('%5i%5i%s\n' % (lbl.layx, lbl.izorp, lbl.xtitle))
                for k in range(lbl.layx): #= 1:layx
    				# write record 3.8.1
                    f.write('%10.3f     %s\n' % (lbl.zorp[k], lbl.jcharXs))
    					
                    # write record 3.8.2
                    s = ''.join(['%10.3e' % gas for gas in lbl.cfc[k,:]])
                    f.write(s + '\n')    
    
        
    # .. Footer
    f.write('-1\n' + '%\n')
    
    # .. Done, close!
    f.close()
    

def read_od_lblrtm(filename = None, prec='s'):
    '''
    Purpose: Read in the wavenumber and optical depth from an LBLRTM optical
             depth file.
    
    Inputs: filename: Name of file output from LBLRTM
            prec: precision of optical depth, 'd' or 's'
    
           Based on rdOD by S. Neshyba,
           Modified for single or double precision, 2016/08/23 P. Rowe
    '''

    

    # .. Open binary file
    fp = open(filename, 'rb')


    # Set params for double or single precision
    if prec == 's':
        nheader = 266*4 + 4     # single
        floaty   = 'f' #'float32'     # 4 bytes
        inty     = 'i' #'int'
    elif prec == 'd':
        nheader = 356*4 + 4     # double
        floaty   = 'd' #float64'     # 8 bytes
        inty     = 'l' #'int64'
    else:
      raise NameError('Input variable "prec" must be "s" or "d".')
    

    #
    # .. Reads the header.  Note that the header is NOT decoded since this
    #    information is not used below.  It is 1068 bytes long.
    #
    header = array.array("B")  # B is a typecode having 1 byte
    header.fromfile(fp, nheader) #1068)

    #
    # .. Reads blocks of data.
    nlim = 0
    nrecords = 0
    spectrum = np.zeros(10000000)
    i1 = 0
    
    # .. Continues to read data until EOF is encountered.
    while True:
    
        # Read the first block header
        v1 = array.array('d'); v1.fromfile(fp, 1)
        v2 = array.array('d'); v2.fromfile(fp, 1)
        dv = array.array(floaty); dv.fromfile(fp, 1)
        nlima = array.array(inty); nlima.fromfile(fp, 1)
        nlim = nlima[0]
 
        if (nlim < 1):
            fp.close()
            break
        elif (nrecords == 0):
            v1keep = v1[0]
            v2keep = v2[0]   # in case its needed
            dvkeep = dv[0]
        else:
            v2keep = v2[0]
        eor1 = array.array('i'); eor1.fromfile(fp,2)

        # Build the spectrum vector block-by-block.
        sblock = array.array(floaty); sblock.fromfile(fp, nlim)
        i2 = i1 + len(sblock)
        spectrum[i1:i2] = sblock
        i1 = i2
        
        # Skip the post-data marker.
        try:
            eor2 = array.array('i'); eor2.fromfile(fp,2)
        except:
            raise NameError('Was not able to read eor2 at end: see code')
        
        # Increment our counters
        nrecords += 1


    # .. Go just above the upper limit to avoid dropping wavenumbers
    wn = np.arange(v1keep, v2keep+dvkeep/4, dvkeep)
    spectrum = spectrum[:len(wn)]
    
    return(wn, spectrum)
    
    


