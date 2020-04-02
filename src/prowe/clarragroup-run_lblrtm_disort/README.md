# README #

## run_LBLRTM_DISORT

### Overview ###

This set of Python codes is for simulating cloudy-sky downwelling infrared radiances in wavenumber microwindows, using LBLRTM and DISORT. The user specifies the cloud properties and atmospheric profile, and uses single scattering parameter files (URL will be provided soon). LBLRTM is used to get gaseous optical depths, which are then input to DISORT to get cloudy-sky radiances.

### Setting up ###

1. Install runDisort_py. This enables you to call DISORT (please see repository runDisort_py, including references for DISORT).
2. Install LBLRTM: http://rtweb.aer.com/lblrtm.html.
3. Install f90nml, through, e.g. $ pip install f90nml.
4. Change the directory in run_lblrtm_disort_micro.py line 81, site.addsitedir("/Users/prowe/Git_repos/rundisort_py/installation/")
to the directory where disort_driver_py.so is. I need to develop a more elegant solution for this.  

### Running run_lblrtm_disort_micro.py ###
1. run_lblrtm_disort_micro.py is the main function. Call it as follows:
    microwindows, nu_eff, Rcloudy, Rclear_disort, Rclear_lblrtm \  
    = run_lblrtm_disort_micro(inputfile)
where inputfile is the name of the input namelist.


2. You will need to specify all other inputs in the namelist file. An example namelist file, "inputs.nml" is provided in the folder sample_run.


3. The outputs are returned by the function (see run_lblrtm_disort_micro.py) and saved to a file, and are:
    * lower wavenumber bound (cm-1)
    * upper wavenumber bound (cm-1)
    * effective wavenumber (cm-1)
    * Cloudy radiance, from DISORT (mW/[m2 sr cm-1]) 
    * Clear radiance, from LBLRTM (mW/[m2 sr cm-1]) 
    * Clear radiance, from DISORT (mW/[m2 sr cm-1])

4. Single-scattering parameters (pmomfiles) are needed and will be provided through a website as they are too large to share here.

5. Atmospheric profiles are needed that specify height, temperature, pressure, and trace gas amounts with altitude. These are netcdf files (netcdf4 classic) that contain boundary layer values for altitude (z) in km, pressure (P) in mb, and temperature (T), in K. The units must also be specified, and currently only these units are accepted. Trace gases are specified by gas, in lower case (e.g. h2o, co2, etc). Gases specified in the lblrtm instruction file, through hbr, may be inputted. Units must be specified, and the only currently acceptable units are ppm, gm_kg, and "B" (referenced to value for jchar in LBLRTM). See also the class Prof in lblrtm_utils.py. CFCs f11, f12, and f113 may also be inputed and must have units "A" (referenced to value for jcharX in LBLRTM). The number of molecules ("nmol") may optionally be specified. The variable modelExtra must be set to the desired model for all unspecified trace gases. An example file is supplied in sample_run.


### Credits and Conditions of Use ###

This code builds on a lot of previous work. If this code is used in work leading to a publication, please reference the following:

* Overall reference for code: Rowe, P. M., Cox, C., Neshyba, S., & Walden, V. P. (2019). Toward autonomous surface-based infrared remote sensing of polar clouds: retrievals of cloud optical and microphysical properties. Atmospheric Measurement Techniques, 12(9), 5071–5086. http://doi.org/10.5194/amt-12-5071-2019.

* Single-scattering parameters (pmomfiles). In progress ... If you use these please acknowledge Rowe et al 2020 (in review) and references therein.
* DISORT: see references in runDisort_py for appropriate attribution.
* LBLRTM: see references at the website (http://rtweb.aer.com/lblrtm.html) for appropriate attribution.

### Contact ###

* If you have any questions, please contact Dr. Penny Rowe (penny@nwra.com).
