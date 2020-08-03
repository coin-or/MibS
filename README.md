# MibS (Mixed Integer Bilevel Stochastic Solver) 2.0

## Build Status

[![Build Status](https://travis-ci.org/coin-or/MibS.svg?branch=master)](https://travis-ci.org/coin-or/MibS)

[![Build status](https://ci.appveyor.com/api/projects/status/aqxs9wcp2tjgpffd?svg=true)](https://ci.appveyor.com/project/tkralphs/mibs-gkymh/branch/master)
## Download

[ ![Download](https://api.bintray.com/packages/coin-or/download/MibS/images/download.svg?version=1.1.2) ](https://bintray.com/coin-or/download/MibS/1.1.2/link)

Binary packages are available for some platforms from [Bintray](https://bintray.com/coin-or/download/MibS).

## Cite

[![DOI](https://zenodo.org/badge/39053653.svg)](https://zenodo.org/badge/latestdoi/39053653)

## DESCRIPTION

MibS is a solver for stochastic mixed integer bilevel linear optimization problems. For an
introduction to bilevel optimization, see [this slide
deck](http://coral.ie.lehigh.edu/~ted/files/talks/BILEVEL-IWOBIP16.pdf). A
somewhat outdated but still useful introduction to MibS is
[here](http://coral.ie.lehigh.edu/~ted/files/talks/BILEVEL-INFORMS11.pdf). 
A paper that contains a complete technical description of the algorithms in MibS is 
[here](http://coral.ie.lehigh.edu/~ted/files/papers/MIBLP16.pdf).

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

Most Linux distributions come with all the required tools installed. To obtain
the source code, the first step is to get the installer that will then
fetch the source for `MibS` and all its dependencies. *You do not need to
clone the repository first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/coinbrew
```

Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `coinbrew`
subdirectory. To execute the script, do

```
cd coinbrew
chmod u+x coinbrew
./coinbrew
```

(Note: The `chmod` command is only needed if the execute permission is not
automatically set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coinbrew fetch --no-prompt MibS:stable/x.y
./coinbrew build --no-prompt MibS --prefix=/path/to/install/dir
./coinbrew install MibS
```
Note that the prefix specified above is the directory where the packages will be
installed. If the specified prefix is writable, then all packages will be
automatically installed immediately after building. If no prefix is specified,
the package will be installed in the directory dist/. Options that would have
been passed to the `configure` script under the old build system can simply be
added to the command-line. For example, to build with debugging symbols, do

```
./coinbrew build --no-prompt MibS --prefix=/path/to/install/dir --enable-debug
```

To get help with additional options available in running the script, do

```
./coinbrew --help
```

After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `LD_LIBRARY_PATH` if you want to link to COIN libraries. 

### Building on Windows (MSys2/CYGWIN and MinGW/MSVC)

By far, the easiest way to build on Windows is with the GNU autotools and the
GCC compilers. The first step is to install either
   * [Msys2](https://msys2.github.io/)
   * [CYGWIN](http://cygwin.org/)
   * [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
Bash and the gcc compilers also come with the [Anaconda Python distribution](https://www.anaconda.com/distribution/)

If you don't already have CYGWIN installed and don't want to fool around with
WSL (which is a great option if you already know your way around Unix), it is
recommended to use MSys2, since it provides a minimal toolset that is easy to
install. To get MSys2, either download the installer
[here](https://msys2.github.io/) or download and unzip MSys2 base from
[here](http://kent.dl.sourceforge.net/project/msys2/Base/x86_64/msys2-base-x86_64-20190512.tar.xz) 
(this is an out-of-date version, there may be a better place to get an archive
version). 

Following any of the above steps, you should have the `bash` command
(with Msys2, be sure to run `msys2_shell.bat` 
or manually add `msys64\usr\bin`, `msys64\mingw32\bin`, and
`msys64\mingw64\bin` to your Windows path).   

Once you have bash installed and in your `PATH`, open a Windows terminal and
type 

```
bash
pacman -S make wget tar patch dos2unix diffutils git svn
git clone https://www.github.com/coin-or/coinbrew
```
Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `COIN-OR-OptimizationSuite`
subdirectory. To execute the script, do

```
cd coinbrew
chmod u+x coinbrew
./coinbrew
```
(Note: The `chmod` command is only needed if the execute permission is not
automatically set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coinbrew fetch --no-prompt MibS:stable/x.y
./coinbrew build --no-prompt MibS --prefix=C:\path\to\install\dir
./coinbrew install MibS
```
Note that the prefix specified above is the directory where the packages will be
installed. If the specified prefix is writable, then all packages will be
automatically installed immediately after building. If no prefix is specified,
the package will be installed in the directory dist/. Options that would have
been passed to the `configure` script under the old build system can simply be
added to the command-line. For example, to build with debugging symbols, do
```
./coinbrew build --no-prompt MibS --prefix=C:\path\to\install\dir --enable-debug
```

To get help with additional options available in running the script, do

```
./coinbrew --help
```

To use the resulting binaries and/or libraries, you will need to add the
full path of the directory `build\bin` to your Windows executable
search `PATH`, or, alternatively, copy the conents of the build directory to 
`C:\Program Files (x86)\MibS` and add the directory
`C:\Program Files (x86)\MibS\bin` 
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
./coinbrew build --no-prompt MibS --prefix=C:\path\to\install\dir --enable-msvc
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

To obtain
the source code, the first step is to get the installer that will then
fetch the source for MibS and all its dependencies. *You do not need to
clone MibS first, just do the following!* Open a terminal and execute

```
git clone https://www.github.com/coin-or/coinbrew
```

Next, to check out source code for and build all the necessary projects
(including dependencies), execute the script in the `coinbrew`
subdirectory. To execute the script, do

```
cd coinbrew
chmod u+x coinbrew
./coinbrew
```

(Note: The `chmod` command is only needed if the execute permission is not
automatically set by git on cloning). Once you run the script,
you will be prompted interactively to select a project to fetch and build. The
rest should happen automagically. Alternatively, the following command-line
incantation will execute the procedure non-interactively.

```
./coinbrew fetch --no-prompt MibS:stable/x.y
./coinbrew build --no-prompt MibS --prefix=/path/to/install/dir
./coinbrew install MibS
```
Note that the prefix specified above is the directory where the packages will be
installed. If the specified prefix is writable, then all packages will be
automatically installed immediately after building. If no prefix is specified,
the package will be installed in the directory dist/. Options that would have
been passed to the `configure` script under the old build system can simply be
added to the command-line. For example, to build with debugging symbols, do

```
./coinbrew build --no-prompt MibS --prefix=/path/to/install/dir --enable-debug
```

To get help with additional options available in running the script, do

```
./coinbrew --help
```
After installation, you will also need to add `/path/to/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/path/to/install/dir/lib`
to your `DYLD_LIBRARY_PATH` if you want to link to COIN libraries. 

## USING

To solve a deterministic mixed integer bilevel linear optimization problem, you must provide both an MPS file and an auxiliary
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
the header file `MibParams.hpp` for explanations. You can also find a detailed
description of MibS
[here](http://www.optimization-online.org/DB_FILE/2017/04/5977.pdf).
Furthermore, to conduct the experiments illustrated in this report, see
the `README` file in the directory `scripts`.

`MibS` is also capable of solving two-stage mixed integer stochastic bilevel linear optimization problems. To solve these problems, there are two ways:

* You must provide an MPS file, a time file and a stoch file in the same way as the SMPS format. The second-stage objective coefficients should be defined at the end of the time file (see [here](https://github.com/tkralphs/BilevelLib/blob/master/stochastic/sslp/bilevel_nonZeroSum_sslp_10_50_50.tim)). Then call `mibs` like this:
```
<build_or_install_dir>/bin/mibs -Alps_instance file.mps -MibS_auxiliaryTimFile file.tim -MibS_auxiliaryStoFile file.sto -MibS_stochasticityType stochasticWithoutSAA -MibS_isSMPSFormat 1 -MibS_isA2Random 1
```
The parameter `MibS_isA2Random` should be set to 0 in case the coefficients of the first-stage variables in the second-stage problem are not random variables.

* You must provide an MPS file and an auxiliary information file in the same way as deterministic bilevel problems. The probability distributions of the random variables also should be specified by setting the values of corresponding parameters (MibS currently supports only the discrete uniform distribution). For a sample parameter file, see `src\mibsStochastic.par.in`.

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
