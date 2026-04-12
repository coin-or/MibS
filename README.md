# MibS 1.2 - Improving Directions

This README explains how to download and build this specific fork of MibS in order to replicate the computational experiments illustrated in this [manuscript](https://arxiv.org/abs/2511.03566).

For information on how to replicate the experiments illustrated in the manuscript, see the
`README` file in the directory `scripts` ([here](https://github.com/febattista/MibS/tree/improving-directions/scripts)).

For general usage instructions and information about other releases of MibS, please refer to the official
[MibS repository](https://github.com/coin-or/MibS).

## BUILDING from source

The quick start assumes you are in a bash shell. 

### Using `coinbrew` (Recommended)

To download MibS from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch https://github.com/febattista/MibS@improving-directions
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

The recommended command to build MibS is the following
```
./coinbrew build MibS@improving-directions -b build-idBC-opt -p --tests none -j 12
```

This command builds MibS inside the directory `./build-idBC-opt` (referred to as `<build-dir>` in the
`scripts` README). The executable can be found at `<build-dir>/bin/mibs`.
The source for the manuscript was compiled using the version 12.2.0 of `g++`.

## USING

### Command Line

To solve a deterministic mixed integer bilevel linear optimization problem,
you must provide both an MPS file and an auxiliary information file that
specifies which variables and constraints are associated with each level
(see a description of the file format
[here](https://coin-or.github.io/MibS/input.html)).
Then call `mibs` like this: 
``` 
<build-dir>/bin/mibs -Alps_instance file.mps -MibS_auxiliaryInfoFile file.aux 
``` 
Specifying the path to the auxiliry file is unnecessary provided it has the same root name as the MPS file and is in the same location. 

It is also possible to
specify additional settings either on the command line or in a parameter file with, e.g.,  
```
<build-dir>/bin/mibs -Alps_instance file.mps -MibS_branchStrategy 1
```
or
```
<build-dir>/bin/mibs -param <build-dir>/MibS/src/mibs.par 
``` 
MibS has many parameters. See the
example parameter file `mibs.par` and the header file `MibParams.hpp` for
explanations. 

## Project Links

 * [Additional documentation](https://coin-or.github.io/MibS)
 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [MibS Discussion forum](https://github.com/coin-or/MibS/discussions)
 * [Report a bug](https://github.com/coin-or/MibS/issues/new)
 * [Doxygen generated documentation](http://coin-or.github.io/MibS/Doxygen)

## ACKNOWLEDGEMENT

MibS was developed with support from

 * Office of Naval Research (Grant N000141912330)
 * National Science Foundation (Grants CMMI-1435453, CMMI-0728011, ACI-0102687)
 * Lehigh University
 * Zuse Institute Berlin
 * Research Campus Modal "Mathematical Optimization and Data Analysis 
   Laboratories" funded by the German Federal Ministry of Education and Research
   (BMBF Grant 05M14ZAM) and by the DFG SFB/Transregio 154

