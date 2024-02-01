---
layout: default
title: Quick Start
navigation: 1
---

# MibS Quick Start

## Versions

`MibS` uses [semantic versioning](https://semver.org/) (major.minor.patch). There are three major
versions of `MibS` currently under development.

  * Versions with major number 1 are recommended for most users. This is the
original solver, which is focused on solution of _deterministic_ mixed integer
bilevel linear optimization problems (MIBLPs). The latest minor version in
this series is [1.2](https://github.com/coin-or/MibS/tree/stable/1.2) and
the most recent release is
[![Latest Release](https://img.shields.io/github/v/release/coin-or/MibS?sort=semver)](https://github.com/coin-or/MibS/releases). 
  * The version in the [master](https://github.com/coin-or/MibS/tree/master)
branch will become major version 2 and can additionally solve _stochastic_
MIBLPs. Once released, this version will subsume all capabilities of version 1
and version 1 will be deprecated. It has not yet had a release, but is well tested.
  * There is also a major verion 3 of `MibS` under development that incorporates the
ability to solve MIBLPs under the assumption of _bounded rationality_.

## Building and Installing

For information on building and installing `MibS` on different platforms, see
the README for the appropriate version, e.g.,
[here](https://github.com/coin-or/MibS/tree/stable/1.2#building-from-source)
(procedure is similar for all versions), or the more detailed [COIN-OR User
Guide](https://coin-or.github.io/user_introduction.html).

## Using

To understand how to use `MibS`, it helps to know that `MibS` is built on top of
the [BLIS](https://github.com/coin-or/CHiPPS-BLIS) MILP solver, which is in
turn implemented using the (COIN-OR High-Performance Parallel Search (CHiPPS)
Framework and more specifically with the [Abstract Library for Parallel Search
[ALPS](https://github.com/coin-or/CHiPPS-ALPS), a library of abstract base
classes for implementing parallel search algorithms. You can find a detailed
description of the algorithms and the important parameters of major version 1
of `MibS`
[here](https://arxiv.org/abs/2104.09010http://www.optimization-online.org/DB_FILE/2017/04/5977.pdf).

### Modelling Systems

MibS has interfaces to the following modelling systems that allow the user to
conveniently build the bilevel model in a high-level modelling language and
pass the model to MibS through the interface for solution.

- Pyomo (Python) through the
[PAO](https://pao.readthedocs.io/en/latest/index.html#) package,
see [here](https://pao.readthedocs.io/en/latest/solvers.html)
- JuMP (Julia) through the [BilevelJuMP](https://github.com/joaquimg/BilevelJuMP.jl) package, see
[here](https://github.com/joaquimg/BilevelJuMP.jl/blob/master/docs/src/examples/MibS_example1.jl)
for an example. 

### Basic Command-line Usage

The command-line interface on all platforms consists of the executable name `mibs`
followed by a list of pairs of parameters and values. Parameter names should always be preceded
by a `-`. Some of the parameters are `MibS` parameters, but parameters of ALPS and BLIS may also 
be necessary.

#### Version 1

To solve a deterministic mixed integer bilevel linear optimization problem
with `MibS` version 1, you must provide both an MPS, LP, or GMPL/AMPL file and
an auxiliary information file that specifies which variables and constraints
are associated with the each level (see [here](input.html)). Then call `mibs` like this:
```bash
mibs -instance file.mps -auxiliaryInfoFile file.aux
```
Note that specifying the name and location of the auxiliary is unnecessary if 
it has the same base name and is in the same folder.

`MibS` has many additional parameters, many of which are documented
[here](parameters.html). See the example parameter file 
[`mibs.par`](https://github.com/coin-or/MibS/blob/stable/1.2/examples/mibs.par.in) and the
header file [`MibParams.hpp`](https://github.com/coin-or/MibS/blob/stable/1.2/src/MibSParams.hpp) 
for additional information. It is also possible to specify additional settings in a parameter file with, e.g., 
```bash
mibs -param MibS/src/mibs.par 
```

#### Version 2

Major version 2 of`MibS` is additionally capable of solving two-stage mixed
integer stochastic bilevel linear optimization problems. There are two
different ways to specify such problems to `MibS`. 

  * Provide an MPS file, a time file and a stoch file in the same way as the
    [SMPS](http://www.maximalsoftware.com/resources/GassmannKristjansson_dpm007v1.pdf)
    format. The second-stage objective coefficients should be defined at the
    end of the time file (see
    [here](https://github.com/tkralphs/BilevelLib/blob/master/stochastic/sslp/bilevel_nonZeroSum_sslp_10_50_50.tim)).
    For this instance format, call `mibs` like this:
    ```bash
    dist/bin/mibs -Alps_instance file.mps -MibS_auxiliaryTimFile file.tim -MibS_auxiliaryStoFile file.sto -MibS_stochasticityType stochasticWithoutSAA -MibS_isSMPSFormat 1 -MibS_isA2Random 1 
    ``` 
    The parameter `MibS_isA2Random` should be set to 0 in case the
    coefficients of the first-stage variables in the second-stage problem are
    not random variables.

  * Provide an MPS file and an auxiliary information file in the same way as
    deterministic bilevel problems. In this case, the probability
    distributions of the random variables must also be specified by setting
    the values of corresponding parameters (`MibS` currently supports only the
    discrete uniform distribution). For a sample parameter file, see
    [https://github.com/coin-or/MibS/blob/master/src/mibsStochastic.par.in].

### The C++ API

Major versions 1 and 2 both have a C++ class library that can be used to link
to the MibS from another application. The main class user interact with is the
`MibSModel` class. Below is a short code snippet showing how to pass instance
data to MibS in memory from another application. As described above, MibS is
ultimately built on the abstract base classes of
[ALPS](https://github.com/coin-or/CHiPPS-ALPS) and hence the solve call is
actually invoked by calling the `search` method of the `AlpsKnowledgeBroker`
object. (In the future, we plan to wrap all of this in order to make
invocation simpler to understand).

```cpp
#include <iostream>

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

#include "MibSConfig.hpp"
#include "MibSModel.hpp"

#include "AlpsKnowledgeBrokerSerial.h"

int main(int argc, char* argv[])
{
   /** Set up lp solver **/
   OsiClpSolverInterface lpSolver;
   lpSolver.messageHandler()->setLogLevel(0);
      
   /** Create MibS model **/
   MibSModel model;
   model.setSolver(&lpSolver);

   /****************************************************************/
   /*                                                              */
   /* Add code here that fills in the data structures needed for   */
   /* the functions below.                                          */
   /*                                                              */
   /****************************************************************/

   /* The auxiliary data describes which rows and variables in the */
   /* upper level and which are in the lower level, as well as     */
   /* specifying the lower-level objective function                */
   
   model.loadAuxiliaryData(lowerColNum, lowerRowNum, lowerColInd,
                            lowerRowInd, 1.0, lObjCoeff,
                            upperColNum, upperRowNum, upperColInd,
                            upperRowInd, structRowNum, structRowInd,
                            0, NULL, lColLbInLProb, lColUbInLProb);

   /* The problem data itself is passed in as a CoinPackedMatrix,  */
   /* vectors of upper and lower bounds, etc., as usual for an     */
   /* MILP solver                                                  */
   
   model.loadProblemData(*newMatrix, varLB, varUB, objCoef, conLB,
                          conUB, colType, 1, mps->getInfinity(),
                          rowSense);

   int argc = 1;
   char** argv = new char* [1];
   argv[0] = "mibs";

   /* Create the knowledge broker to do the actual search. argc */
   /* argv are passed in to parse the arguments                 */
   
   AlpsKnowledgeBrokerSerial broker(argc, argv, model);

   /** This is the function that does the solve **/
   
   broker.search(model);

   /** Get the solution **/

   if (model->getNumSolutions() == 0){
      std::cout << "MibS did not find any bilevel feasible solutions. "
               << std::endl;
      abort();
   }

   MibSSolution *solution = dynamic_cast<MibSSolution* >
      (broker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
   
   double *y = new double[upperColNum];
   
   for (j = lowerColNum; j < numTotalCols; j++){
      y[j - lowerColNum] = floor(solution->getValues()[j] + 0.5);
   }
}
```

There are also function in `MibSModel` for directly parsing MPS/LP/AMPL/GMPL
files, as well as axiliary files if thoe are to be read directly. See
[MibSModel.hpp](https://github.com/coin-or/MibS/blob/master/src/MibSModel.hpp).
