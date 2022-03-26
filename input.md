---
layout: page
title: File Format
navigation: 3
---

# MibS Input File Format

## Version 1

Major version 1 of MibS requires two input files.
 1. An _instance file_ in MPS, LP, or GMPL/AMPL format that contains a
 description of the upper level objective, the variables (upper and lower
 level), and the constraints (upper and lower level); and 
 2. An _auxiliary (aux) file_ that contains information necessary to
 identify the upper and lower level constraints and variables.

### The Instance File

The instance can be in any of the standard formats mentioned above. These are
the same formats utilized by all MILP solvers. Information about these formats
can be easily found with an Internet search. 

### The Auxiliary File

Each line in the auxiliary file contains a pair of elements: a type indicator
and a data value. The type indicators are as follows.

| Type Indicator | Associated Data                                     |
|----|-----------------------------------------------------------------|
| N  | the number of lower level variables                             |
| M  | the number of lower level constraints                           |
| LC | name/index of one of the lower level variables                  |
| LR | name/index of one of the lower level constraint                 |
| LO | coefficient of a variable in the lower level objective function |
| OS | lower level objective sense (1=min, -1=max)                     |

Note that variables and constraints can be identified either by index or name.
When indentified by index, variables are assumed to be in the order they
appear in the input file (indices start at 0). Similarly, constraints are
assumed to be in the order they appear in the input file (indices start at 0).
For the lower level objective, coefficients are assumed to be in the order of
lower level variables. Currently, there is no option for specifying the
lower-level objective using names.

### Example

The model:
```
min  -x − 7y
s.t. −3x + 2y ≤ 12
     x + 2y ≤ 20
     x ≤ 10
     y ∈ arg min {z : 2x - z ≤ 7,
                      2x + 4z ≤ 16,
                      z ≤ 5}
```
The MPS file:

```
NAME generalExample
ROWS
 L  R0
 L  R1
 L  R2
 L  R3
 N  R4
COLUMNS
    INT1 'MARKER' 'INTORG'
    C0   R0       -3
    C0   R1       1
    C0   R2       2
    C0   R3       -2
    C0   R4       -1
    C1   R0       2
    C1   R1       2
    C1   R2       -1
    C1   R3       4
    C1   R4       -7
    INT1END 'MARKER' 'INTEND'
RHS
    B    R0       12
    B    R1       20
    B    R2       7
    B    R3       16
BOUNDS
 UP BOUND C0      10
 UP BOUND C1      5
ENDATA
```

The auxiliary file:

```
N 1
M 2
LC 1
LR 2
LR 3
LO 1
OS 1
```

Optimal solution: x = 6, y = 5

## Version 2

Major version 2 accepts the same file format as version 1 for deterministic
instances and additionally accepts stochastic instances in the following two
formats.

 1. Provide an MPS file, a time file, and a stoch file, as in the well-known
    [SMPS format](http://www.maximalsoftware.com/resources/GassmannKristjansson_dpm007v1.pdf). The second-stage objective coefficients should be defined at
    the end of the time file (see
    [here](https://github.com/tkralphs/BilevelLib/blob/master/stochastic/sslp/bilevel_nonZeroSum_sslp_10_50_50.tim)).

 2. Provide an MPS file and an auxiliary information file in the same way as
    deterministic bilevel problems. The probability distributions of the
    random variables also should be specified by setting the values of
    corresponding parameters (`MibS` currently supports only the discrete
    uniform distribution). For a sample parameter file, see
    [src/mibsStochastic.par.in].

