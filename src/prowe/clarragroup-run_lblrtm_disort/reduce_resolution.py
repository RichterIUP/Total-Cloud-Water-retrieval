# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 19:55:13 2015

@author: prowe

Copyright 2018-2020 by Penny M. Rowe and NorthWest Research Associates.
All rights reserved.
"""


import numpy as np
from scipy.fftpack import fft
from scipy.fftpack import ifft


def reduce_resolution_pad(nu_mono, calc_mono, dnu, 
                          n_interferogram, n_padded, 
                          n_hi_res):

  """
  Purpose: 
      Reduce resolution according to dnu, assuming n_interferogram pts 
      in interfogram and using n_hi_res points to zero-pad the spectrum
  
  By Penny M. Rowe, 2020/02/21, based on aeri.m
  
  Inputs:
      nu_mono
      calc_mono
      dnu: desired resolution
      n_interferogram: number of points in interferogram
      n_padded: number of points after padding
      n_hi_res: large number of points, for interpolation
       
  
  Notes:
      We want to put the calculated interferogram on the same spacing
      as the measured interferogram so we can chop it at the same
      place. This is defined by NptsInLoResInt
  
      xmax = n_interferogram/laserWn 
      dnu = 1/2/xmax
      numax = dnu*(n_interferogram-1)
      nu = (0:dnu:numax)
  
      Since we are given dnu, rather than laserWn, we have
      xmax = 1/2/dnu
      

      The final wavenumber vector needs to be:
      xmax = NptsInPaddedInt/laserWn ;
      dnu = 1/2/xmax;
      numax = dnu*(NptsInPaddedInt-1);
      wnf=(0:dnu:numax);
  
  """
 
  
  # .. The calculated spectrum needs to begin at zero wavenumbers
  #    (it typically has a spacing of 0.0003 cm-1)
  #    and it needs to end at numax to give dx properly
  bwn = nu_mono[0]
  ewn = nu_mono[-1]
  
  if n_hi_res==0 or np.isnan(n_hi_res):
      n_hi_res = 2**22                     # large number of points


  numax = dnu * n_interferogram            # dnu = lnu/N/2
                                           # so lnu = dnu*N*2
                                           # and numax=lnu/2 or dnu*N
  
  
  nu_mono_interp = np.linspace(0, numax, n_hi_res+1) # e.g. 0.030s
  nu_mono_interp = nu_mono_interp.T
  
  # .. Create new radiance vector with tapered ends. 
  #    Create np.zeros array (0.019 s)
  calc_mono_interp = np.zeros(n_hi_res + 1)
  ind = np.where(np.logical_and(nu_mono_interp >= bwn, nu_mono_interp <= ewn))[0]
  # i1 = int(ind[1])

  
  # Interpolate onto the new grid (0.014 s)
  # this should probably use the cubic interpolation, but there
  # is not a good one in Python and it would take even longer.
  # What are the resulting errors?
  # linear interpolation: 0.013 s for first layer
  # cubic spline: 0.07 to 0.12 s for first - 37th layer
  # Cubic spline adds about 5s
  # However, the improvement is not noticeable.
  calc_mono_interp[ind] = np.interp(nu_mono_interp[ind], nu_mono, calc_mono)
  #spl = interpolate.splrep(nu_mono, calc_mono)
  #calc_mono_interp[ind] = interpolate.splev(nu_mono_interp[ind], spl)
  #print('time to interpolate onto new grid: ' + str(time.time() - start))
	

  #start = time.time()
  # Low-wavenumber roll-off
  # set up index vectors (0.012 s)
  ind_lwn = np.arange(ind[1]); 
  #ind_lwn = np.arange(1e5); 
  #ind_hwn = np.arange(ind[-1], n_hi_res+1); 
  ind_hwn = np.arange(ind[-1], ind[-1]+3000); 
  

  # Compute and paste on the low-wavenumber roll-offs (0.13 s)
  hwn = ind_hwn - ind[-1] + 1
  
  # ind_lwn
  # i1-100000:i1
  calc_mono_interp[ind_lwn] = ( ( np.cos( ind_lwn*np.pi/ind[0] ) )[::-1] +1 ) \
                                  / 2 * calc_mono_interp[ind[0]]
  calc_mono_interp[ind_hwn] = (   np.cos( hwn * np.pi/(len(hwn)   ))  +1 ) \
                                  / 2 * calc_mono_interp[ind_hwn[0]]


  # .. Calculate the interferogram
  #    slow step (e.g. 1.1 s)
  x, ifg = ift(nu_mono_interp, calc_mono_interp, 1, 1)

  
  # .. Pad with zeros starting from the interferogram length+1 + 1
  #    In Python we count from zero, so this is index n_interferogram+1
  #    To the padded length (n_padded+1)
  if (not np.isnan(n_padded)) & (n_padded > n_interferogram):
      ifg[n_interferogram+1:n_padded+1] = 0
  
  # .. Truncate the interferogram to the padded length (n_padded+1)
  #    and take the Fourier Transform to get back the spectrum
  #    Note that in Python we count from zero but also omit the last value
  #    hence 0:n_interferogram+1 => 1st value to n_interferogram+1+1-1th value
  #    or 1st value to n_interferogram+1th value
  #    (0.002 s)
  nu, calc = ft(x[:n_padded+1], ifg[:n_padded+1], 1, 1)
  
  
  #
  # .. Extract the portion of the spectrum between the original
  #    wavenumber limits. (0.0002 s)
  ind = np.where(np.logical_and(nu >= bwn, nu <= ewn))[0]
  nu = nu[ind]
  calc = calc[ind]


  return(nu, calc)



def ift(x, y, sidesin, sidesout):
  """
  #function [a,b] = ift3(x,y,sidesin,sidesout)
  #
  # ifft y, where y=f(x)
  # get y = ifft(yft), and xft
  # where y = f(x)
  #
  #
  # sides = 1 for one-sided (will be made two-sided before FT)
  #			note: for one-sided can have any number of points
  #			assumes that the spectrum starts from 0 frequency
  #
  # Note: the format must be as follows:
  #
  # 	for single sided:
  # 		spectrum must start with zero frequency and end with
  #       nyquist critical frequency.  It can have even or odd
  #       number of points (will end up even)
  #
  #	for double sided:
  #		must be in following format already ...
  #		zero frequency up to nyquist frequency, then
  #		negative of one less than nyquist frequency down to
  #		negative of one more than zero frequency
  #
  #   to optimize we should use a factor of 2 number of points for
  #       double-sided but in this code we use factor of 2 number of
  #       points plus one point (zero nu)
  #       this should be fixed in future versions.
  #
  #
  # From "help fft" in matlab:
  #   The functions Y=fft(x) and y=ifft(X) implement the transform
  #   and inverse transform pair given for vectors of length N by:
  #
  #   For length N input vector x, the DFT is a length N vector X,
  #   with elements
  #                    N
  #      X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
  #                   n=1
  #   The inverse DFT (computed by IFFT) is given by
  #                    N
  #      x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
  #                   k=1
  #
  #   See also FFT2, FFTN, FFTSHIFT, FFTW, IFFT, IFFT2, IFFTN.
  #
  #
  # Inverse Fourier Transform
  # Penny Rowe
  # July 26, modified March 3, 2000 and March 26, 2000
  """
  
  # make y vector into row vectors
  #if dims(1)==1
  #  # already row vector
  #elif dims(2)==1
  #  y=y';
  #end
    
  
  # Note: assuming we start with an interferogram and fft to get a spectrum
  # or alternatively start with a spectrum and ifft to get an interferorgram
  # Then the following is true
  #
  # spectrum(1) = sum(interferogram)
  # interferogram(1) = sum(spectrum)/10
  #
  # So, in going from the spectrum to the interferogram (ift), we want
  # the integral under the spectrum to equal the integral under the
  # interferogram, so we make the first point be the integral under the
  # spectrum
  
  
  if sidesin == 1:
    
    # Flip to make the left hand side (0.032s)
    rhs = np.flipud(y[1:-1])    # : => end, :-1 => end-1
    
    # note: the spectrum must be anti-symmetric! (0.030s)
    #       2 x number of points - 2
    if all(np.isreal(y)):
      y = np.hstack([y, rhs])  #[y rhs];
    else:
      y = np.hstack([np.real(y) - 1j * np.imag(y), rhs])
    
  elif sidesin ==2:
    x = x[x>=0]
  else:
    raise NameError('sidesin =1 for one-sided, 2 for 2-sided')

  
  # I assume everything is perfect already
  # Now ifft using built in code (e.g. 0.92 s)
  #start = time.time()
  b = ifft(y)
  #end = time.time()
  #print('time for Pythons np.ifft: ' + str(end - start))

  
  if sidesout==2:
    # make a a vector of integers, starting from zero (0.18s)
    N = 2*len(x)-2
    n = np.array(range(len(x)))   #(0:length(x)-1);
    #n = [n -fliplr(n(2:end-1))]; #replaced length(n) with end
    #n = np.hstack([n,-np.flipud(n[1:-2])])    
    n = np.hstack([n,-np.flipud(n[1:-1])])      # -1 => end-1
    dx = x[1] - x[0]  #dx = x(2)-x(1);
  
    # now multiply by da to scale properly
    a = (1/(N*dx)) * n
    
  elif sidesout==1:
    # make a a vector of integers, starting from zero (0.18s)
    N = len(x)
    dx = x[1] - x[0]  #dx = x(2)-x(1);
    
    # now multiply by da to scale properly
    a = 1/((2*N-2)*dx) * np.arange(N)  #(1/((2*N-2)*dx)) * (0:N-1);
    b = b[:N]   #b(1:N);
    
    
  return a, b



def ft(x,y,sidesin,sidesout):
  """
  #
  # function [nu, yft] = ft(x,y,sidesin,sidesout)
  #
  # fft y, where y=f(x)
  # get yft = fft(y), and nu
  # where yft = f(nu)
  #
  # sides = 1 for one-sided (will be made two-sided)
  # sides = 2 for two-sided
  #
  # Fourier Transform
  # Penny Rowe
  # July 26, modified March 3, 2000
  #
  # For speed, x and y should be length N+1, where
  # N is a power of 2. This is because we take
  # only y(2:end-1) to form the negative half,
  # so that N-1 + N-1-2 = N.
  """
  
  # make vectors into row vectors
  #if dims(1)==1
  #  # already row vector
  #elseif dims(2)==1
  #  x=x';
  #end
  
  #if dims(1)==1
  #  # already row vector
  #elseif dims(2)==1
  #  y=y';
  #end
  

  if sidesin == 1:
    y2 = np.flipud(y[1:-1])       #fliplr(y(2:end-1));
    y = np.hstack([y, y2])        #y = [y y2]; # original
    
    # We have 2*length(y)-2 points, so if N is a power of 2,
    # we should have started with N+1 points, to give
    # 2*(N+1)-2 = 2N points
    
    # The following appears not to work:
    #y2 = fliplr(y(2:length(y)-1));
    #rhs = fliplr(y(2:length(y)-1));
    #y = [rhs real(y)-1i*imag(y) ];
    
  elif sidesin ==2:
    # get positive values of x only
    x = x[x>=0]  #x = x(x>=0);
  else:
    raise NameError('error: sidesin =1 for one-sided, 2 for 2-sided')

  
  yft = fft(y)
  
  
  if sidesout==2:
    # make a a vector of integers, starting from zero
    N = 2*len(x)-2
    n = np.arange(len(x))                #(0:length(x)-1);
    #n = np.hstack([n,-np.flipud(n[1:-2])])   #[n -fliplr(n(2:end-1))];
    n = np.hstack([n,-np.flipud(n[1:-1])])   #end-1 => -1
    dx = x[1] - x[0]                       #dx = x(2)-x(1);
    
    nu = (1/(N*dx)) * n                    #nu = (1/(N*dx)) * n;
  elif sidesout==1:
    # make a a vector of integers, starting from zero
    N = len(x)
    dx = x[1] - x[0]             #x(2)-x(1);
    nu = 1/((2*N-2)*dx) * np.arange(N) #nu  = (1/((2*N-2)*dx)) * (0:N-1);
    yft = yft[:N+1] #yft = yft(1:N);


  return nu, yft