# How to install L-IWP

## 1. Download source codes:

- [LBLRTM](http://rtweb.aer.com/lblrtm.html)
- [DISORT (LBLDIS)](https://www.nssl.noaa.gov/users/dturner/public_html/lbldis/index.html)
- L-IWP: git clone https://github.com/RichterIUP/L-IWP.git 

## 2. LBLRTM:

- Extract the source code: 
```sh
> tar -xzvf aer_lblrtm_v12.8.tar.gz
```

LBLRTM requires a spectroscopic database. There are two ways to deal with this:

2.1 Set up the database on your own. To do so, one needs to extract LNFL and the HITRAN database. Next, lnfl needs to be compiled and then the database can be created
```sh
> tar -xzvf lnfl_v3.1.tar.gz
> tar -xzvf aer_v_3.6.tar.gz 
> make -f make_lnfl linuxGNUsgl
> cd ../../aer_v_3.6/line_file
> ../../lnfl/lnfl_v3.1_linux_gnu_sgl
```

2.2 Use the file TAPE3.10-3500cm-1.first_7_molecules

- Create lblrtm/hitran and place the following files:
```sh
> ln -s path/to/TAPE3/file/TAPE3.10-3500.cm-1.first_7_molecules ./tape3.data
> ln -s ../run_examples/xs_files/xs ./xs
> ln -s ../run_examples/xs_files/xs ./x
> ln -s ../run_examples/xs_files/FSCDXS ./FSCDXS
```
- Compile LBLRTM
```sh
> cd lblrtm/build
> make -f make_lblrtm linuxGNUsgl
```
Maybe there could be some errors. Then you should fix the bugs. First bug seems to be in lblatm.f90, line 7967. Insert a whitespace between STOP and the string.
	
- Create the folder lblrtm/bin and link the binary
```sh
> mkdir lblrtm/bin
> ln -s ../lblrtm_v12.8_linux_gnu_sgl ./lblrtm
```

## 3. LBLDIS:

- Create folder lbldis and copy the archive into it. Extract the archive: 
```sh
> mkdir lbldis
> cd lbldis
> tar -xvf lbldis.Release_3_0.tar 
> chmod +w *
```

- By default, LBLDIS wants to be compiled using IFORT. If you want to use GFORTRAN, then change following lines in Makefile:
```sh
- Line 20: Change ifort to gfortran
- Line 26: Remove -nofor_main
- Line 120 and 129: Change -r8 to -fdefault-real-8 -fdefault-double-8
```

- Compile the source
```sh
> make
```

## 4. L-IWP:

- Set the paths in inp.py:
```python
path = path/to/retrieval/output
outfolder = path/to/results/output
exec_lbldis = '(path/to/lbldis/binary {}/.lbldis_{}.parm 0 {}) >& {}\n'
path_to_lblrtm = 'path/to/lblrtm/main/folder' (not to the bin-folder!)
```

## 5. Create virtual environment for Python 3 

If you don't want to 'contaminate' your local Python, you can create a virtual environment

```sh
> sudo apt-get install python3-venv
> python3 -m venv path/to/python/virtenv
> source path/to/python/virtenv/bin/activate
> path/to/python/virtenv/bin/pip3 install numpy
> path/to/python/virtenv/bin/pip3 install scipy
> path/to/python/virtenv/bin/pip3 install matplotlib
> path/to/python/virtenv/bin/pip3 install netcdf4
```
	
## 6. Test installation:

- Go to L-IWP and run liwp_main.py
```sh
path/to/python/virtenv/bin/python3 liwp_main.py 
```
