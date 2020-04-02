"""
utils.py

Copyright 2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.

Purpose:
  Utilities for performing radiative transfer using LBLRTM and DISORT
"""


# .. Built in modules
import numpy as np
from scipy.interpolate import interp1d
from scipy.io.netcdf import NetCDFFile as DS
from scipy.interpolate import interp2d


class SSP:
    __slots__ = ('NPmom', 'Pmom', 'reff', 'w0', 'qext',
                 'fNPmom', 'fPmom', 'fqext',  'fw0')
    
    def __init__(self, pmomfile, nu):
        self.NPmom, self.Pmom, \
        self.reff, self.w0, self.qext = load_ssp_nu(pmomfile, nu)
        
        self.fNPmom = interp1d(self.reff, self.NPmom.T)
        self.fqext = interp1d(self.reff, self.qext.T)
        self.fw0 = interp1d(self.reff, self.w0.T)


def get_ssp(pmomfiles, nu):
    Nssp = len(pmomfiles)
    ssp = {}
    for i in range(Nssp):
        ssp[i] = SSP(pmomfiles[i], nu)
            
    return ssp


def load_ssp_nu(datafile, nu):
    
    # Load the moments from a netcdf file & return the unbundled arrays
    dnu = nu[1] - nu[0]
    with DS(datafile,'r') as nc:
    
        inu = np.where(np.logical_and(\
                              nc.variables['wnum_list'][:] > nu[0]-2*dnu, 
                              nc.variables['wnum_list'][:] < nu[-1]+2*dnu))[0]
        Npmomarray = nc.variables['Npmomarray'][:,inu].astype('int32')
        w0_mesh    = nc.variables['w0_mesh'][:,inu].astype('float64')
        qext_mesh  = nc.variables['qext_mesh'][:,inu].astype('float64')
        reff = nc.variables['reff_list'][:].astype('float64'); reff = reff[:,0]
        wnum_vec  = nc.variables['wnum_list'][inu].astype('float64')
        # wnum_mesh  = nc.variables['wnum_mesh'][:,inu].astype('float64')
        # reff_mesh  = nc.variables['reff_mesh'][:,inu].astype('float64')
        
        
        # Set up empty output arrays
        Nnu = nu.size
        Nreff = reff.size
        qext = np.zeros((Nreff, Nnu))
        w0 = np.zeros((Nreff, Nnu))
        NPmom_fp = np.zeros((Nreff, Nnu))
    
        # Interpolate qext, w0, get an interpolated number of moments!
        fq = interp2d(reff, wnum_vec, qext_mesh.T)
        fw = interp2d(reff, wnum_vec, w0_mesh.T)
        fNP = interp2d(reff, wnum_vec, Npmomarray.T)
    
        for i in range(Nreff):
            qext[i,:] = fq(reff[i], nu)[:,0]
            w0[i,:] = fw(reff[i], nu)[:,0]
            NPmom_fp[i,:] = fNP(reff[i], nu)[:,0]
        
        # Use floor so we never interpolate between a moment and 0.
        NPmom = np.floor(NPmom_fp).astype(int) 
        NPmom_max = np.max(NPmom)
        pmomarray  = nc.variables['pmomarray'][:,inu,:NPmom_max]
        pmomarray = pmomarray.astype('float64')
        
        # Loop over all the moments to do the same
        Pmom = np.zeros((Nreff, Nnu, NPmom_max));
        for j in range( NPmom_max):
            f = interp2d(reff, wnum_vec, pmomarray[:,:,j].T)
            for i in range(Nreff):
                Pmom[i,:,j] = f(reff[i], nu)[:,0]

    return (NPmom, Pmom, reff, w0, qext)



def interp1_weights2(y=None, yi_list=None, method=None):
    #    
    # Return indices and weights to pmom files
    #
    # Steven Neshyba and Penny Rowe
    #
    # Inputs:
    #   y = temperatures corresponding to CRI
    #   yi_list = ssp temperatures
    #   method = liqTdependence: interopolate, none_average
    # 
    # Updates:
    #   Modified to include average for ice ssp data (can only
    #   have two files! Further mods needed!!)
    #   PMR, 2018/07/04
    #
    #
    # Examples:
    # 
    # input: yi_list=np.array([240.1]); 
    #        y=np.array([250.,240.,260.,270.]); 
    #        method = 'interpolate'
    # output: I = [[1 0]] and w = [[ 0.99  0.01]]
    # 
    # input: yi_list=np.array([269.0]); 
    #        y=np.array([250.,240.,260.,270.]); 
    #        method = 'nearest'
    # output: [[3]] and w = [[ 1.]]
    # 
    # input: yi_list=np.array([280.]); 
    #        y=np.array([250.,240.,260.,270.]); 
    #        method = 'interpolate'
    # output: I = [[3 3]] and w =[[ 1.  0.]]
    # 
    # input: yi_list=np.array([[253.],[267.]]); 
    #        y=np.array([250.,240.,260.,270.]); 
    #        method = 'interpolate'
    # output: I = [[0 2]    and w = [[ 0.7  0.3]
    #              [2 3]]            [ 0.3  0.7]]
    #
    # Used by cloudModel.cloudModel

    # Default output values
    Iout = 0
    wout = 0

    # Check for length
    Ny = np.size(y) # Number of ys
    if Ny == 1:
        Iout = np.zeros((len(yi_list),1), dtype=np.int32)
        wout = np.ones((len(yi_list),1))
        return Iout, wout

    # Sort the y's in increasing order
    ysorted = np.sort(y) 
    iysorted = np.argsort(y)  
    N_yi_list = np.size(yi_list)
    
    for iyi in range(N_yi_list):

        # Extract one value from the desired list
        yi = yi_list[iyi]               # Current interpolated y (yi)

        # Take care of cases when yi is out of the given range
        if (yi >= ysorted[-1]):         # yi is bigger than biggest y
            I1 = Ny-1
            w1 = 1  
            I2 = Ny-1
            w2 = 0  
        elif (yi <= ysorted[0]):        # yi is smaller than smallest y
            I1 = 0
            w1 = 1  
            I2 = 0
            w2 = 0  
        else:
            dum = np.arange(0,Ny)
            Ifloat = np.interp(yi, ysorted, dum, yi)
            I1 = np.floor(Ifloat).astype(np.int32)
            I2 = np.ceil(Ifloat).astype(np.int32)
            if (I2 == I1):
                I2 = I1 + 1
            w1 = (I2 - Ifloat) / (I2 - I1)
            w2 = 1 - w1
            
        # Store pointers and weights of the closest two values in given range
        if iyi == 0:
            I = np.zeros((1,2),dtype=np.int32)
            I[0] = [iysorted[I1], iysorted[I2]]
            w = np.zeros((1,2)); w[0] = [w1, w2]
        else:
            Inext = np.zeros((1,2),dtype=np.int32)
            Inext[0] = [iysorted[I1], iysorted[I2]]
            wnext = np.zeros((1,2)); wnext[0] = [w1, w2]
            I = np.concatenate((I,Inext))  
            w = np.concatenate((w,wnext))  
    
    # .. Keep as-is if the algorithm is to interpolate, 
    #    evenly weight files if non_average,
    #    otherwise choose nearest
    if (method == 'interpolate'):
        Iout = I
        wout = w
    elif (method == 'nearest'):
        Iout = np.zeros((N_yi_list,1),dtype=np.int32)
        for iyi in range(N_yi_list):
            wnext = w[iyi] #rint wnext
            iwnearest = wnext.argmax() 
            Iout[iyi][0] = I[iyi][iwnearest]
        wout = np.ones(np.shape(Iout))
    elif (method=='none_average'):
      # For now, assume we only have two pmom files
      # Later we will allow more combinations
      # No interpolation to temperatures, just average them all
      # Initialize:
      Iout = np.zeros((1,2),dtype=np.int32)
      wout = np.zeros((1,2))
      # Set values and weights
      Iout[0] = [0, 1]
      wout[0] = [.5, .5]
    else:
        raise ValueError('Method is not recognized')
    
    # Exit
    return Iout, wout





