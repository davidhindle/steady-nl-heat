# steady-nl-heat
[![DOI](https://sandbox.zenodo.org/badge/633796407.svg)](https://sandbox.zenodo.org/badge/latestdoi/633796407)
## Introduction

ellip is a fortran 90 code, consisting of a main code, and 2 modules, which implements the half station discretisation of the 1 dimensional, non-linear (temperature dependent thermal conductivity) steady state heat flow equation with a variable coefficient. 
 

## Important

This is a "bare bones" repository to get the source code for the paper, Anees et al (2023) published and available for download. These instructions are therefore not failsafe. If you are not experienced in Fortran 90 usage and Linux shell environments, things may go wrong and you may not understand why. Moroever, it is important to note that this program has been developed, compiled and run using the GNU Fortran compiler gfortran on a Linux system. The code employs some rather specific tricks, for instance 128 bit quadruple precision variables, and it cannot be guaranteed to work with other Fortran compilers and on other systems. I would thus welcome reports (up to the point of saturation) from any users on other platforms who find problems, and especially from users who find AND SOLVE those problems.

## Requirements

Ideally <br/>
1) Linux or BSD or similar UNIX type environment <br/>
2) GNU fortran (gfortran) compiler <br/>
3) GMTv6.x installed (for visualisation) <br/><br/>
If you have these, you are pretty much okay. If you have a Windows machine, it may be possible to use the power shell and the g95 compiler to get things to work, but there are no guarantees and this has not been tested.

## Steps to running the code.

Download the main code and modules.<br/>
ellip.f90 (main code)<br/>
module_tdma.f90 (pentadiagonal matrix solver)<br/>
pyplot.py (python code for plotting needs matplotlib/numpy environment)<br/>
Makefile<br/>
place all 3 files in a directory of your choice, best to create fresh one for the purpose.<br/>
Open a terminal window and navigate to the directory.<br/><br/>
Type make<br/>
You should now have an exectable file ellip.exe. If you are an experienced Linux user, you will probably place the executable somewhere in your path, and make it globally available on your system. You need no help from me here.<br/><br/>
For less experienced users however, the simple way to run the program is to type
`./ellip.exe`

## required input files 

Compiling the code is not everything. The code is structured to read very specific input file with very specific formats, with very specific meanings. Some idea of these can be got from reading the comments within the main source code, and the input file itself.<br/><br/>
The input files for the code is<br/><br/>
`heat-in.dat` <br/>
the variables in `heat-in.dat` are, in order<br/><br/>
**radiative** - true/false; true = radiative component of thermal conductivity activated <br/>
**bckind**- true/false; true = basal boundary condition = Neumann (heat flow), false = Dirichlet (fixed temperature) <br/>
**tol** - iteration tolerance (sum of difference of successive iteration steps). < tol = iteration stops <br/>
**kr**- temperature cut off for radiative heat, °K; 1000 = standard value <br/>
**dx**- grid spacing, metres <br/>
**bc1** - surface boundary condition, temp °C <br/>
**bc2** - basal boundary condition, bckind = false, temp °C; bckind = true, heat flow, W/m^2, +ve = flow into model<br/>
**ko(n)** - integer, number of layers to be given in succeeding line<br/>
**layer properties** - layer k_o (W/m K), layer heat production (W/m³)<br/>
the layer properties are given as two values per line, with as many lines as the value ko(n) 0<br/>

`<br/><br/>


## output and plots
output is in text file format. The accompanying python script pyplot.py will produce 3 separate figures showing temperature, heat flow and thermal conductivity profiles from a model run


## Authors
**David Hindle**, University of Göttingen, <dhindle-at-gwdg.de>


## Reference

This notebook has been published at Zenodo. Please cite the following reference if you publish any work that uses this notebook:

[![DOI](https://sandbox.zenodo.org/badge/633796407.svg)](https://sandbox.zenodo.org/badge/latestdoi/633796407)


## License
This project is licensed under the GNU lesser general public license (LGPL v3). See the [LICENSE.txt](LICENSE.txt) file for details.


