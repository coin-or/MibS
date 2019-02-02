# MibS (Mixed Integer Bilevel Solver) 1.1

***IMPORTANT: The procedure for obtaining the source code for MibS and
dependencies has changed. Please read below and do not clone this project
directly.***

## Build Status

[![Build Status](https://travis-ci.org/coin-or/MibS.svg?branch=master)](https://travis-ci.org/coin-or/MibS)

[![Build status](https://ci.appveyor.com/api/projects/status/aqxs9wcp2tjgpffd?svg=true)](https://ci.appveyor.com/project/tkralphs/mibs-gkymh/branch/master)
## Download

[ ![Download](https://api.bintray.com/packages/coin-or/download/MibS/images/download.svg?version=1.1.2) ](https://bintray.com/coin-or/download/MibS/1.1.2/link)

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

Most Linux distributions come with all the required tools installed. To obtain
the source code, the first step is to get the installer that will then
fetch the source for MibS and all its dependencies. *You do not need to
clone MibS first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/COIN-OR-OptimizationSuite
```

Next, to check out source code for and build all the necessary projects (including 
dependencies), execute the script in the `COIN-OR-OptimizationSuite` subdirectory. 
To execute the script, do

```
cd COIN-OR-OptimizationSuite
chmod u+x coin.install.sh
./coin.install.sh
```

(Note: The `chmod` command is only needed if the execute permission is not automatically 
set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS
```

Options to the `configure` script can simply be added to the command-line. For
example, to build with debugging symbols, do

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS --enable-debug
```

To get help with additional options available in running the script, do

```
./coin/install.sh --help
```

The above procedures will build all required dependencies and MibS itself.
Afterwards, the binaries will be installed in the directory `Mibs/build/bin`
and the libraries in the directory `MibS/build/lib`. If you wish to install in
a different directory, such as `/usr/local`, then run the command

```
./coin.install.sh install --prefix=/path/to/install/dir
```

After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `LD_LIBRARY_PATH` if you want to link to COIN libraries. 

### Building on Windows (MSys2/CYGWIN and MinGW/MSVC)

By far, the easiest way to build on Windows is with the GNU autotools and the
GCC compilers. The first step is to install either
   * [Msys2](https://msys2.github.io/) or
   * [CYGWIN](http://cygwin.org/).
   * [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
If you don't already have CYGWIN installed and don't want to fool around with
WSL (which is a great option if you already know your way around Unix), it is
recommended to use MSys2, since it provides a minimal toolset that is easy to
install. To get MSys2, either download the installer
[here](https://msys2.github.io/) or download and unzip MSys2 base from
[here](http://kent.dl.sourceforge.net/project/msys2/Base/x86_64/msys2-base-x86_64-20150512.tar.xz) 
(this is an out-of-date version, there may be a better place to get an archive version). 

Following any of the above steps, you should have the `bash` command (with MSys2m be sure to run `msys2_shell.bat` 
or manually add `msys64\usr\bin`, `msys64\mingw32\bin`, and `msys64\mingw64\bin` to your Windows path).   

Once you have bash installed and in your `PATH`, open a Windows terminal and type 

```
bash
pacman -S make wget tar patch dos2unix diffutils git svn
```

To obtain the source code, the first step is to get the installer that will then
fetch the source for MibS and all its dependencies. *You do not need to
clone MibS first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/COIN-OR-OptimizationSuite
```

Next, to check out source code for and build all the necessary projects (including 
dependencies), execute the script in the `COIN-OR-OptimizationSuite` subdirectory. 
To execute the script, do

```
cd COIN-OR-OptimizationSuite
chmod u+x coi.install.sh
./coin.install.sh
```

(Note: The `chmod` command is only needed if the execute permission is not automatically 
set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. the
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS
```
Options to the `configure` script can simply be added to the command-line. For
example, to build with debugging symbols, do

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS --enable-debug
```

To get help with additional options available in running the script, do

```
./coin/install.sh --help
```

To use the resulting binaries and/or libraries, you will need to add the
full path of the directory `build\bin` to your Windows executable
search `PATH`, or, alternatively, copy the conents of the build directory to 
`C:\Program Files (x86)\MibS` and add the directory `C:\Program Files (x86)\MibS\bin` 
to your Windows executable search `PATH`. You may also consider adding
`C:\Program Files (x86)\MibS\lib` to the `LIB` path and 
`C:\Program Files (x86)\MibS\include` to the `INCLUDE` path. 

It is possible to use almost the exact same commands to build with the Visual
Studio compilers. Before doing any of the above commands in the Windows
terminal, first run the `vcvarsall.bat` script for your version of Visual
Studio. Note that you will also need a compatible Fortran compiler if you want
to build any projects requiring Fortran (`ifort` is recommended, but not
free). Then follow all the steps above, but replace the `build` command
with

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS --enable-msvc
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

To obtain the source code, the first step is to get the installer that will then
fetch the source for MibS and all its dependencies. *You do not need to
clone MibS first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/COIN-OR-OptimizationSuite
```

Next, to check out source code for and build all the necessary projects (including 
dependencies), execute the script in the `COIN-OR-OptimizationSuite` subdirectory. 
To execute the script, do

```
cd COIN-OR-OptimizationSuite
chmod u+x coi.install.sh
./coin.install.sh
```

(Note: The `chmod` command is only needed if the execute permission is not automatically 
set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. the
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS
```

With this setup, `clang` will be used for compiling C++ by default and
`gfortran` will be used for Fortran. Since `clang` uses the GNU standard
library, `gfortran` is compatible.

If you want to use the `gcc` compiler provided by Homebrew, then replace the
`build` command above with

```
./coin.install.sh build --no-prompt --main-proj=MibS CC=gcc-5 CXX=g++-5
```

Additional options to the `configure` script can simply be added to the command-line. For
example, to build with debugging symbols, do

```
./coin.install.sh fetch build --no-prompt --main-proj=MibS --enable-debug
```

To get help with additional options available in running the script, do

```
./coin/install.sh --help
```

If you wish to install in a different directory, such as `/usr/local`, then run the command

```
./coin.install.sh install --prefix=/path/to/install/dir
```

After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
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
