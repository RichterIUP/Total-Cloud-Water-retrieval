# How to run L-IWP

This repository contains the Python code of the retrieval algorithm L-IWP.

-------
## Dependencies

To run this program, one needs to install [Python 3](https://www.python.org/). As forward models, [LBLRTM](http://rtweb.aer.com/lblrtm.html) and [DISORT (LBLDIS)](https://www.nssl.noaa.gov/users/dturner/public_html/lbldis/index.html) are used.
Packages required to execute the Python code are

```sh
> pip install numpy
> pip install scipy
> pip install matplotlib
> pip install netcdf4
```

The models LBLRTM and DISORT have to be set up in the folder ./radiative_transfer. One needs to set the environment variable LBL_HOME:

```sh
> export LBL_HOME=./radiative_transfer/lblrtm
```

-------
## Running the retrieval

To run the retrieval, type

```sh
> python liwp_main.py "PATH/TO/CASE/FILE"
```

This case file looks as follows:
```sh
SPEC:	./PATH/TO/SPECTRUM
RASO:	./PATH/TO/RADIOSONDE
SZA:	SOLAR_ZENITH_ANGLE
CBH:	CLOUD_BASE_HEIGHT(METER)
CTH:	CLOUD_TOP_HEIGHT(METER)
```

-------
## Preparing the data

The spectrum needs to be in the following shape:

```sh
Wavenumber (cm-1);Spectral radiance(mW/[m2*sr*cm-1])
```

The spectral resolution should be 0.15cm-1 and the spectral range should be at least from 700cm-1 to 1200cm-1.

The atmospheric state is read from the RASO-file, which needs to be in the following shape:

```sh
Height (km),Pressure (hPa),Temperature (K),absolute humidity (g/m3)
```

The data should be at least stated up to 30 km.

Solar zenith angle is needed for the solar influence on the spectral radiances. The cloud base height and the cloud top height can be set to the same value, which leads to a single-layer retrieval.

-------
## Running testcases
