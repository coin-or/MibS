---
layout: page
title: Parameters
navigation: 2
---

# Parameters

## ALPS Parameters

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| Alps_instance | Path to instance file |
| Alps_timeLimit | Overall wallclock time limit |
| Alps_nodeLimit | Limit on number of nodes explored |
| Alps_msgLevel | Verbosity level for ALPS message. Set > 5 for additional debugging output |
| Alps_logFileLevel | Whether to create log file and its verbosity. Currently, this option dumps parameter settings at level 1 and that's it. |
| Alps_searchStrategy |  `0`: Best first <br> `1`: Best estimate <br>  `2`: Breath first <br> `3`: Depth first <br> `4` Hybrid |
| Alps_nodeLogInterval | How frequently to print status lines |

## BLIS Parameters

### General

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| Blis_scaleConFactor | Level of dynamism above which cuts will be discarded. Note that setting this too low can result in incorrect results, since cuts necessary for validity may be thrown out |
| Blis_denseConFactor | Cut density above which cuts will be discarded. Note that setting this too low can result in incorrect results, since cuts necessary for validity may be thrown out |
| Blis_heurStrategy | When to call heuristics (default) <br> `0`: disable <br> `1`: root <br> `2`: auto <br> `3`: periodic |
| Blis_heurCallFrequency | How often to call heuristics |
| Blis_heurRoundStrategy | When to call rounding heuristic (default) <br> `0`: disable <br> `1`: root <br> `2`: auto <br> `3`: periodic |
| Blis_heurRoundFreq | How often to call heuristics |
| Blis_branchStrategy |  MibS uses pseudocost branching by default and the other strategies are untested. <br> `0`: max infeasibility <br> `1`: pseudocost <br> `2`: reliability <br> `3`: strong |

### Cuts

These are parameters for controlling generation of inequalites valid for
MILPs, which can be used to eliminate fractional solutions if desired. These
cuts are all off by default, as they're usually not effective. 

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| Blis_cutStrategy | Strategy for cut generation (default) <br> `0`: disable <br> `1`: root <br> `2`: auto <br> `3`: periodic |
| Blis_cutGenerationFrequency | How frequently to generate MILP cuts |
| Blis_cutCliqueStrategy <br> Blis_cutGomoryStrategy <br> Blis_cutFlowCoverStrategy <br> Blis_cutKnapsackStrategy <br> Blis_cutMirStrategy <br> Blis_cutOddHoleStrategy <br> Blis_cutProbingStrategy <br> Blis_cutTwoMirStrategy | Strategy for generating individual classes of inequalities. <br> `0`: disable <br> `1`: root <br> `2`: auto <br> `3`: periodic |
| Blis_cutCliqueFreq <br> Blis_cutGomoryFreq <br> Blis_cutFlowCoverFreq <br> Blis_cutKnapsackFreq <br> Blis_cutMirFreq <br> Blis_cutOddHoleFreq <br> Blis_cutProbingFreq <br> Blis_cutTwoMirFreq | Frequency for generating individual classes of inequalities. 

## MibS Parameters

### General

These are general MibS parameters

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| MibS_auxiliaryInfoFile | Path to auxiliary info file |
| MibS_usePreprocessor | `0`: off <br> `1`: on |
| MibS_bilevelProblemType |  `0`: general <br> `1`: interdiction |
| MibS_objBoundStrategy |  How to derive a bound on lower level objective for preprocessing <br> `0`: LL obj bound <br> `1`: interdiction bound |
| MibS_whichActiveConMethod |  How to determine which constraints are binding <br> `0`: simple <br> `1`: basis |
| MibS_upperFileFormat |  `0`: MPS <br> `1`: AMPL/GMPL |

### Subsolvers

These are parameters related to solving subproblems for either checking
feasibility or determining the best lower-level solution associated with a
given upper-level solution by solving an auxiliary MILP.

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| MibS_solveSecondLevelWhenXYVarsInt | `0`: false <br> `1`: true|
| MibS_solveSecondLevelWhenXVarsInt | `0`: false <br> `1`: true|
| MibS_solveSecondLevelWhenLVarsInt | `0`: false <br> `1`: true|
| MibS_solveSecondLevelWhenLVarsFixed | `0`: false <br> `1`: true|
| MibS_computeBestUBWhenXVarsInt | `0`: false <br> `1`: true|
| MibS_computeBestUBWhenLVarsInt | `0`: false <br> `1`: true|
| MibS_computeBestUBWhenLVarsFixed | `0`: false <br> `1`: true|
| MibS_useLinkingSolutionPool | `0`: false <br> `1`: true|
| MibS_doDualFixing | `0`: false <br> `1`: true|
| MibS_feasCheckSolver | Options are currently <br> -`Cbc` <br> -`SYMPHONY` <br> -`CPLEX` |
| MibS_warmStartLL | `0`: false <br> `1`: true|
| MibS_maxThreadsLL | Number of threads to use for parallel solve of lower level problem |
| MibS_whichCutsLL            |  `0`: no cuts <br> `1`: gomory only <br> `2`: all cuts |

### Heuristics

These are parameters for controlling heuristics.

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| MibS_useLowerObjHeuristic |  `-1`: auto <br> `0`: false <br> `1`: true |
| MibS_useObjCutHeuristic |  `-1`: auto <br> `0`: false <br> `1`: true |
| MibS_useWSHeuristic | `-1`: auto <br> `0`: false <br> `1`: true |
| MibS_useGreedyHeuristic |  `0`: false <br> `1`: true |

### Branching

These are parameters for controlling branching.

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| MibS_branchStrategy |  `0`: fractional <br> `1`: linking |

### Cuts

These are parameters for controlling generation of valid inequalities. For an explanation, please see this [technical report]( http://coral.ie.lehigh.edu/~ted/files/papers/MibSCuts20.pdf). 
Note that some parameter names related to cuts were changed in version 1.2 to corespond to the names in the report. The names used in earlier versions still work for backwards compatibility, 
but are not documented. 

| Parameter Name | Description                                     |
|----------------|-------------------------------------------------|
| MibS_cutStrategy |  `0`: branch only <br> `1`: cut only <br> `2`: use cut and branch |
| MibS_maxCutDepth | Deepest level of the tree at which cuts should be generated|
| MibS_turnOffDefaultCuts | Turn off all cuts not explicitly turned on by parameters <br> `0`: false <br> `1`: true |
| MibS_useFractionalCuts | Whether to generate cuts when solution is fractional (see tech report) <br> `0`: false <br> `1`: true |
| MibS_useIntegerNoGoodCut |  Whether to generate integer no good cuts <br> `0`: false <br> `1`: true |
| MibS_useValFuncCut | Whether to generate value function cuts <br> `0`: false <br> `1`: true |
| MibS_useNoGoodCut | Whether to generate no good cuts <br> `0`: false <br> `1`: true |
| MibS_useBendersBinaryCut | Whether to generate Benders binary cuts <br> `0`: false <br> `1`: true |
| MibS_useBendersInterdictionCut | Whether to generate Benders intersection cuts <br>`0`: false <br> `1`: true |
| MibS_bendersInterdictionCutType | Whether to generate a single or multiple cuts (from different solutions) in each iteration <br> `0`: justOne <br> `1`: multiple |
| MibS_useGeneralizedNoGoodCut | Whether to generate generalized no good cuts <br> `0`: false <br> `1`: true |
| MibS_useImprovingSolutionIC | Whether to generate improving solution intersection cuts <br> `0`: false <br> `1`: true |
| MibS_bilevelFreeSetTypeIC | What kind of bilevel free set type to use for improving solution ICs <br> `0`: Derive solution by solving lower level problem to optimality <br> `1`: Derive a solution by solving an auxiliary problem |
| MibS_useImprovingDirectionIC | Whether to generate improving direction intersection cuts <br> `0`: false <br> `1`: true |
| MibS_useHypecubeIC | Whether to generate hypercube intersection cuts <br> `0`: false <br> `1`: true |
| MibS_useBoundCut | Whether to generate this cut (see tech report) <br> `0`: false <br> `1`: true |
| MibS_boundCutOptimal | What kind of bound cut to generate (there is currently only one option) <br> `0`: false <br> `1`: true |


