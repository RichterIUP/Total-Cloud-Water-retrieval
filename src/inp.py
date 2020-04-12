#!/usr/bin/python
'''@package docstring
Control variables for the retrieval
'''
import numpy as np

MODELFRAMEWORK = 'LBLDIS'#'CLARRA'

'''
Path to the output of TCWret
'''
PATH = '/home/phi.richter/TCWret/OUTFOLDER'
#PATH = '/home/philipp/TCWret/OUTFOLDER'


'''
Path to the binary of LBLRTM
'''
PATH_TO_LBLRTM = "/home/phi.richter/radiative_transfer/lblrtm"

'''
Path to the binary of LBLDIS
'''
PATH_TO_LBLDIS = "/home/phi.richter/radiative_transfer/lbldis"

'''
Standard atmosphere for LBLRTM for CO, CH4, N2O
1: tropical model
2: midlatitude summer model
3: midlatitude winter model
4: subarctic summer model
5: subarctic winter model
6: U.S. standard 1976
'''
PREDEF_ATM = 4
PROF_CO2 = False
PROF_O3  = False

'''
Several predefined spectral ranges. FIR_TIR incorporates all far-infrared and
thermal-infrared windows. TIR only uses the thermal-infrared windows and FIR only
uses the far-infrared windows. FIR_MCP incorporates the thermal-infrared and the 
nearest far-infrared window
'''
FIR_TIR = range(21)
TIR = range(4, 20)

'''
Microwindows which are used by L-IWP. This can be either one of the predefined
spectral ranges or one can pass a self-defined list of microwindows
'''
WINDOWS = FIR_TIR

'''
Set the maximum altitude for the atmospheric profile in km
'''
MAX_ALT = 30.0

'''
Deceide wether L-IWP shall convolve the calculated radiances with a boxcar or
not.
'''
CONVOLVE = False

'''
Resolution of the spectrum in cm-1. FTIR spectral radiances will be interpolated
to the chosen resolution. If RESOLUTION is a negative number, then no interpolation
is done. OPD is the Optical path difference of the spectrometer
'''
RESOLUTION = 0.3
OPD = 3.0

'''
Initial parameter mu for the Levenberg-Marquardt-Algorithm. The inversion follows the equation
[J^T S_y_1 J + S_a_1 + mu**2 D]s_n = J^T S_y_1 [y - F(x_n)] + S_a_1 (x_a - x_i)
The matrix D is D = np.diag(JT_W_J) * np.identity(len(MCP)
If LM_INIT equals 0.0, then the iteration changes over to Gauss-Newton
'''
LM_INIT = 1e3#1e-3

'''
Minimum value for mu. The parameter won't go below this value.
'''
LM_MIN = 0.0

'''
Number of CPU which shall be used for the calculation of the derivatives.
This number will be decreased automatically, if less than NUM_OF_CPU
are needed
'''
NUM_OF_CPU = 5

'''
If the present spectrum is a testcase of Cox et al. (2016), then this should
be true. This modifies the cloud boundaries and allows to add noise manually
'''
TESTCASE = True

'''
Convergence criterion due to Rodgers (2000):
[x_(n)-x_(n+1)]^T * S^-1 * [x_(n)-x_(n+1)] << len(x)
This is reached if CONVERGENCE is undershot.
'''
CONVERGENCE = 0.001#2.0

'''
Standard deviation of the measured spectral radiances. If STDDEV is below 
0, then the STDDEV is read from the spectrum file. If it equals 0, then
a non-weighted least squares algorithm is run. Positive values let L-IWP ignore
the value from the spectrum file. If a testcase is retrieved, then positive
values are the standard deviation of the added noise
'''

STDDEV = 0.2

'''
Possible values: %, ppmv, g/m3, g/kg
'''
HUMIDITY = 'g/kg'


'''
Disturb the temperature, absolute value (K)
'''
DISTURB_TEMPERATURE = 0.0

'''
Disturb the humidity relative to the value from the radiosonde
'''
DISTURB_HUMIDITY = 0.0

'''
Add an offset to the spectrum (mW/[sr * m2 * cm-1])
'''
OFFSET = 0.0

'''
Parameters for the first guess: tau_total, f_ice, reff_liq, reff_ice
'''
MCP = [10., 10., 15., 30.]#[ 0.277,  0.230,  9.775, 41.394]#[ 0.415,  0.072, 14.782, 30.017]#[1.,1.,15.,30.]

'''
If this is set to true, then L-IWP searches for a file containing cloud height
informations. These should start with CLOUDS. If no file could be found, then 
the data from the radiance file or from these config file will be used.
'''
USE_CLOUD_FILES = False

'''
Calculate the ice fraction from the cloud files. This should be only set
to true, if the cloud height files contain the lwc and iwc. 
'''
FI_FROM_CLOUDS =  False

'''
The a priori for the optimal estimation inversion. By default, the
a priori equals the first guess MCP
'''
MCP_APRIORI = np.array(MCP[:])

'''
The variance of the a priori and its weighting. This will be converted
to the S_A matrix (watch comment of LM_INIT)
'''
VARIANCE_APRIORI = [40.0**(-2), 40.0**(-2), 20.0**(-2), 30.0**(-2)]
WEIGHT_APRIORI = 1.0

'''
Only retrieve the total optical depth
'''
ONLY_OD = False

'''
Sophisticated adjustment
'''
SOPHISTICATED_ADJ = False

'''
Manual choise of the cloud thresholds. If L-IWP should use the cloud
thresholds from the file, these have to be set to -1
'''
CLOUD_BASE = [-1]#[197.075, 2972.02]
CLOUD_TOP = [-1]#[352.971, 8334.85]

'''
Composition of ice particles. The built in composition is:
IF ri > 16 THEN DROXTAL
ELSE SOLID COLUMN
If BUILTIN equals false, then the user can specify his/her own ice particle
shape composition
'''
BUILTIN = False

'''
Number of different ice shapes. The ice optical depth will be divided into
NUM_OF_SHAPES fractions
'''
NUM_OF_SHAPES = 1

'''
Weighting factor for each ice shape. The content of the list has to follow
[SPHERES, AGGREGATE, BULLET_ROSETTE, DROXTAL, HOLLOW_COLUMN, PLATE, SOLID_COLUMN, SPHEROID]
The sum of SHAPES has to be 1.0
'''
SHAPES = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

'''
Surface emissivity
'''
EMISSIVITY = 1.0#0.9

'''
If this is true, then the retrieval is enabled
'''
STANDARD = True

'''
If this is true, then a forward simulation (after a possible retrieval)
will be performed
'''
FORWARD = False

'''
If this is true, then the calculated radiances will be written to lbldis.spec
'''
WRITE_LBLDIS_TO_FILE = True
