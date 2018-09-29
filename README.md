# MibS (Mixed Integer Bilevel Solver) 1.1

## Build Status

[![Build Status](https://travis-ci.org/coin-or/MibS.svg?branch=master)](https://travis-ci.org/coin-or/MibS)

[![Build status](https://ci.appveyor.com/api/projects/status/aqxs9wcp2tjgpffd?svg=true)](https://ci.appveyor.com/project/tkralphs/mibs/branch/master)
## Download

[ ![Download](https://api.bintray.com/packages/coin-or/download/MibS/images/download.svg?version=1.0.0) ](https://bintray.com/coin-or/download/MibS/1.0.0/link)

Binary packages are available for some platforms from [Bintray](https://bintray.com/coin-or/download/MibS).

## Cite

[![DOI](https://zenodo.org/badge/39053653.svg)](https://zenodo.org/badge/latestdoi/39053653)

## DESCRIPTION

MibS is a solver for mixed integelibrr bilevel optimization problems. For an
introduction to bilevel optimization, see [this slide
deck](http://coral.ie.lehigh.edu/~ted/files/talks/BILEVEL-IWOBIP16.pdf). A
somewhat outdated but still useful introduction to MibS is
[here](http://coral.ie.lehigh.edu/~ted/files/talks/BILEVEL-INFORMS11.pdf).

## SUPPORTED PLATFORMS

MiBS should work on all major patforms (OS X, Linux, and Windows), though the
software is in active development and most recent testing has been on Linux.

## DEPENDENCIES

MibS depends on the [CHiPPS](https://projects.coin-or.org/CHiPPS),
[Cbc](https://projects.coin-or.org/Cbc), and [SYMPHONY]
(https://projects.coin-or.org/SYMPHONY) projects of COIN-OR. If you already
have these installed, you can build and link MibS against the installed
version. Otherwise, by following the instruction below, you'll be able to
download, build, and install all dependencies.

## BUILDING AND INSTALLING

### Building on Linux

Most Linux distributions come with all the required tools installed.
To obtain the source code, open a terminal and do

```
git clone https://www.github.com/tkralphs/MibS
```

To build from source, there is a script that fetches dependent projects
and builds automatically. To get the script and perform the build, do

```
cd MibS
git clone https://github.com/coin-or-tools/BuildTools/
```

and then execute

```
BuildTools/get.dependencies.sh fetch
BuildTools/get.dependencies.sh build --quiet
```

This will build all required dependencies and MibS itself. Afterwards, the
binaries will be installed in the directory `Mibs/build/bin` and the libraries
in the directory `MibS/build/lib`. If you wish to install in a different
directory, such as `/usr/local`, then run the command

```
BuildTools/get.dependencies.sh install --prefix=/path/to/install/dir
```

After installation, you will also need to add `/your/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/your/install/dir/lib`
to your `LD_LIBRARY_PATH` if you want to link to COIN libraries. 

### Building on Windows (MSys2/CYGWIN and MinGW/MSVC)

By far, the easiest way to build on Windows is with the GNU autotools and the
MinGW compilers.  
 1. The first step is to install either [Msys2](https://msys2.github.io/) or
 [CYGWIN](http://cygwin.org/). If you don't already
 have CYGWIN installed, it is recommended to use MSys2, since it provides a
 minimal toolset that is easy to install.    
 2. To get MSys2, either download the installer
 [here](https://msys2.github.io/) or download and unzip MSys2 base from
 [here](http://kent.dl.sourceforge.net/project/msys2/Base/x86_64/msys2-base-x86_64-20150512.tar.xz). 
 3. Either run `msys2_shell.bat` or manually add `msys64\usr\bin`,
 `msys64\mingw32\bin`, and `msys64\mingw64\bin` to your Windows path.   
 4. Open a Windows terminal and type

   ```
   bash
   pacman -S make wget tar patch dos2unix diffutils git svn
   ```

 5. Obtain the source code with 

   ```
   git clone https://www.github.com/tkralphs/MibS
   ```

 6. To build from source, there is a script that fetches dependent projects
and builds automatically. To get the script amd perform the build, do 

```
cd MibS
git clone https://github.com/coin-or-tools/BuildTools/
```

and then execute

```
BuildTools/get.dependencies.sh fetch
BuildTools/get.dependencies.sh build --quiet
```

This will build all required dependencies and MibS itself. Afterwards, the
binaries will be installed in the directory `Mibs/build/bin` and the libraries
in the directory `MibS/build/lib`. If you wish to install in a different
directory, such /c/Program\ Files\ \(x86\)/MibS, then run the command

```
BuildTools/get.dependencies.sh install --prefix=/path/to/install/dir
```

 7. To use the resulting binaries and/or libraries, you will need to add the
 full path of the directory `MibS\build\bin` to your Windows executable
 search `PATH`, or, alternatively, copy this directory to `C:\Program Files
 (x86)` and add the directory `C:\Program Files (x86)\MibS\bin` to your
 Windows executable search `PATH`. You may also consider copying the
 `build\lib` and `build\include` directories if you want to link to the
 COIN-OR libraries.

It is possible to use almost the exact same commands to build with the Visual
Studio compilers. Before doing any of the above commands in the Windows
terminla, first run the `vcvarsall.bat` script for your version of Visual
Studio. Note that you will also need a compatible Fortran compiler if you want
to build any projects requiring Fortran (`ifort` is recommended, but not
free). Then follow all the steps above, but replace the `build` command
with

```
BuildTools/get.dependencies.sh build --quiet --enable-msvc 
```

### Building on OS X

OS X is a Unix-based OS and ships with many of the basic components needed to
build COIN-OR, but it's missing some things. For examples, the latest versions
of OS X come with the `clang` compiler but no Fortran compiler. You may also
be missing the `wget` utility and `subversion` and `git` clients (needed for
obtaining source code). The easiest way to get these missing utilitites is to
install Homebrew (see http://brew.sh). After installation, open a terminal and
do

```
brew install gcc wget svn git
```

```
git clone https://www.github.com/tkralphs/MibS
```

To build from source, there is a script that fetches dependent projects
and builds automatically. To get the script amd perform the build, do

```
cd MibS
git clone https://github.com/coin-or-tools/BuildTools/
```

and then execute

```
BuildTools/get.dependencies.sh fetch
BuildTools/get.dependencies.sh build --quiet
```

This will build all required dependencies and MibS itself. Afterwards, the
binaries will be installed in the directory `Mibs/build/bin` and the libraries
in the directory `MibS/build/lib`.

With this setup, `clang` will be used for compiling C++ by default and
`gfortran` will be used for Fortran. Since `clang` uses the GNU standard
library, `gfortran` is compatible.

If you want to use the `gcc` compiler provided by Homebrew, then replace the
`build` command above with

```
BuildTools/get.dependencies.sh build --quiet CC=gcc-5 CXX=g++-5
```

If you wish to install in a different
directory, such as `/usr/local`, then run the command

```
BuildTools/get.dependencies.sh install --prefix=/path/to/install/dir
```

After installation, you will also need to add `/your/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/your/install/dir/lib`
to your `DYLD_LIBRARY_PATH` if you want to link to COIN libraries. 

## USING

To solve a bilevel program, you must provide both an MPS file and an auxiliary
information file that specifies which variables and constraints are associated
with the each level (see [here](http://coral.ise.lehigh.edu/wp-content/uploads/2016/02/MibS_inputFile.pdf)). Then call `mibs` like this:
```
<build_or_install_dir>/bin/mibs -Alps_instance file.mps -MibS_auxiliaryInfoFile aux_file.txt
```
It is also possible to specify additional settings in a parameter file with,
e.g., 
```
<build_or_install_dir>/bin/mibs -param <build_or_install_dir>/MibS/src/mibs.par
```
MibS has many parameters. See the example parameter file `mibs.par` and
the header file `MibParam.h` for explanations. You can also find a detailed
description of MibS
[here](http://www.optimization-online.org/DB_FILE/2017/04/5977.pdf).
Furthermore, to conduct the experiments illustrated in this report, see
the `README` file in the directory `scripts`.     

HELP

Please post questions and issues to the github project page for MibS.

ACKNOWLEDGEMENT

MibS was developed with support from

* National Science Foundation (Grants CMMI-1435453 and CMMI-0728011)
* Lehigh University
* Zuse Institute Berlin
* Research Campus Modal "Mathematical Optimization and Data Analysis 
Laboratories" funded by the German Federal Ministry of Education and Research
(BMBF Grant 05M14ZAM) and by the DFG SFB/Transregio 154

http://github.com/coin-or/MibS

Enjoy!
