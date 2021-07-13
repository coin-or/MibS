# MibS

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

[![Latest Release](https://img.shields.io/github/v/release/coin-or/MibS?sort=semver)](https://github.com/coin-or/MibS/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

## CITE

[![DOI](https://zenodo.org/badge/39053653.svg)](https://zenodo.org/badge/latestdoi/39053653)

MibS is a solver for stochastic mixed integer bilevel linear optimization problems. For an
introduction to bilevel optimization, see [this slide
deck](http://coral.ie.lehigh.edu/~ted/files/talks/BILEVEL-IWOBIP16.pdf). A
somewhat outdated but still useful introduction to MibS is
[here](http://coral.ie.lehigh.edu/~ted/files/talks/BILEVEL-INFORMS11.pdf). 
A paper that contains a complete technical description of the algorithms in MibS is 
[here](http://coral.ie.lehigh.edu/~ted/files/papers/MIBLP16.pdf).

MibS is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/eclipse-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org)

The MibS development site is https://github.com/coin-or/MibS.

## CITE

[![DOI](https://zenodo.org/badge/39053653.svg)](https://zenodo.org/badge/latestdoi/39053653)

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/MibS/actions/workflows/windows-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/MibS/actions/workflows/windows-ci.yml?query=branch%3Amaster)

[![Linux and MacOS Builds](https://github.com/coin-or/MibS/actions/workflows/linux-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/MibS/actions/workflows/linux-ci.yml?query=branch%3Amaster)

## DOWNLOAD

### Docker image

There is a Docker image that provides MibS, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/MibS/releases). Older binaries
are archived as part of MibS
[here](https://www.coin-or.org/download/binary/MibS).

Due to license incompatibilities, pre-compiled binaries may lack some
functionality. If binaries are not available for your platform for the latest
version and you would like to request them to be built and posted, feel free
to let us know in the discussion formum.

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of MibS from the
 [releases](https://github.com/coin-or/MibS/releases) page.
 * Cloning this repository from [Github](https://github.com/coin-or/MibS) or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

Below is a quick start guide for building on common platforms. More detailed
build instructions are
[here](https://coin-or.github.io/user_introduction.html).

## BUILDING from source

The quick start assumes you are in a bash shell. 

### Using `coinbrew`

To build MibS from source, obtain the `coinbrew` script, do
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch MibS@master
./coinbrew build MibS
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/MibS
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

## USING

To solve a deterministic mixed integer bilevel linear optimization problem,
you must provide both an MPS file and an auxiliary information file that
specifies which variables and constraints are associated with the each level
(see
[here](http://coral.ise.lehigh.edu/wp-content/uploads/2016/02/MibS_inputFile.pdf)).
Then call `mibs` like this: 
``` 
<build_or_install_dir>/bin/mibs -Alps_instance file.mps -MibS_auxiliaryInfoFile aux_file.txt 
``` 
It is also possible to
specify additional settings in a parameter file with, e.g., 
```
<build_or_install_dir>/bin/mibs -param <build_or_install_dir>/MibS/src/mibs.par 
``` 
MibS has many parameters. See the
example parameter file `mibs.par` and the header file `MibParams.hpp` for
explanations. You can also find a detailed description of MibS
[here](http://www.optimization-online.org/DB_FILE/2017/04/5977.pdf).
Furthermore, to conduct the experiments illustrated in this report, see the
`README` file in the directory `scripts`.

`MibS` is also capable of solving two-stage mixed integer stochastic bilevel
linear optimization problems. To solve these problems, there are two ways:

  * Provide an MPS file, a time file and a stoch file in the same way as the
    SMPS format. The second-stage objective coefficients should be defined at
    the end of the time file (see
    [here](https://github.com/tkralphs/BilevelLib/blob/master/stochastic/sslp/bilevel_nonZeroSum_sslp_10_50_50.tim)).
    Then call `mibs` like this: 
    ``` 
    <build_or_install_dir>/bin/mibs -Alps_instance file.mps -MibS_auxiliaryTimFile file.tim -MibS_auxiliaryStoFile file.sto -MibS_stochasticityType stochasticWithoutSAA -MibS_isSMPSFormat 1 -MibS_isA2Random 1 
    ``` 
    The parameter `MibS_isA2Random` should be set to 0 in case the coefficients of
    the first-stage variables in the second-stage problem are not random
    variables.

  * Provide an MPS file and an auxiliary information file in the same way as
    deterministic bilevel problems. The probability distributions of the
    random variables also should be specified by setting the values of
    corresponding parameters (MibS currently supports only the discrete
    uniform distribution). For a sample parameter file, see
    [src/mibsStochastic.par.in].


## Project Links

 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [Discussion formum](https://github.com/coin-or/MibS/discussions)
 * [Report a bug](https://github.com/coin-or/MibS/issues/new)
 * [Doxygen-generated html documentation](http://coin-or.github.io/MibS/Doxygen)

## ACKNOWLEDGEMENT

MibS was developed with support from

 * Office of Naval Research (Grant N000141912330)
 * National Science Foundation (Grants CMMI-1435453, CMMI-0728011, ACI-0102687)
 * Lehigh University
 * Zuse Institute Berlin
 * Research Campus Modal "Mathematical Optimization and Data Analysis 
   Laboratories" funded by the German Federal Ministry of Education and Research
   (BMBF Grant 05M14ZAM) and by the DFG SFB/Transregio 154

