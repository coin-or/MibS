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

### The Auxiliary File (Name-based)

In the name-based format (recommended), keywords are proceeded by an `@`. The
following are the keywords.

| Keyword       | Meaning                                                         |
|---------------|-----------------------------------------------------------------|
| @NUMVARS      | Next line contains umber of lower-level variables               |
| @NUMCONSTR    | Next line contains number of lower-level constraints            |
| @OBJSENSE     | Next line contains objective sense (`MIN` or `MAX`)             |
| @VARSBEGIN    | Marks beginning of variables section                            |
| @VARSEND      | Marks end of variables section                                  |
| @CONSTRBEGIN  | Marks beginning of constraints section                          |
| @CONSTREND    | Marks end of constraints section                                |
| @NAME         | Next line contains the name of the instance (optional)          |
| @MPS          | Next line contains the name of the MPS file with which this instance is associated |
| @LP           | Next line contains the name of the LP file with which this instance is associated  |

In the variables section, each line consists of the name of one of the variables from 
the instance file that should be aken to be a lower-level variables, followed by its objective 
coefficient in the lower-level problem (separated by a space). In the constraints section, each 
row consists of the name of one of the constraints from the instance file that should be 
considered a lower-level constraint. 

Note that the bounds on the lower-level variables are always considered constraints of the
lower-leel problem, although they could technically also be taken as constraints at the 
upper-level. If the latter is desired, simply represent them explicitly as constraints rather
than listing them in the variable bounds section of the instance file. 

The instance file itself can be either in MPS or LP format, so exactly one of the keywords
`@MPS` and `@LP` can appear in the auxiliary file. Specifying a name is optional and the name 
can be different than the name of the underlying instace from the MPS/LP file. This is to
make it possible to associate multiple axiliary files with a single underlying instance file.

#### Example

The model:
```
min  -x − 7y
s.t. −3x + 2y ≤ 12
     x + 2y ≤ 20
     x ≤ 10
     y ∈ arg min {z : 2x - z ≤ 7,
                      -2x + 4z ≤ 16,
                      z ≤ 5}
```
The MPS file `genealExample.mps`:

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
The auxiliary file `generalExample.aux`:
```
@NUMVARS
1
@NUMCONSTRS
4
@OBJSENSE
MIN
@VARSBEGIN
LV 1.
@VARSEND
@CONSTRSBEGIN
R1
R2
R3
R4
@CONSTRSEND
@NAME
General Example
@MPS
generalExample.mps
```
To solve with defult parameters, invoke MibS with
```
mibs -Alps_instance generalExample.mps
```
Note that if the auiliary file has the same name as the MPS/LP file (but with the `aux` extension, 
it is not necessary to specify the name of the auxiliary file. 

For this example, the optimal solution should be: `x = 6, y = 5`

### The Auxiliary File (Legacy Index-based)

There is an older file format that is maintained for legacy purposes. In this format, each line in 
the auxiliary file contains a pair of elements: a type indicator
and a data value. The type indicators are as follows.

| Type Indicator | Associated Data                                     |
|----|-----------------------------------------------------------------|
| N  | the number of lower level variables                             |
| M  | the number of lower level constraints                           |
| LC | name/index of one of the lower level variables                  |
| LR | name/index of one of the lower level constraint                 |
| LO | coefficient of a variable in the lower level objective function |
| OS | lower level objective sense (1=min, -1=max)                     |

In this format, the variables and constraints are identified by index.
This is not recommended, as it is rather fragile. It depends on the underlying
MIP solver ordering the constraints and variables in the same way as the MPS/LP file,
as well as requiring the user to manually index the constraints from the file.
When indentified by index, variables and constraints are assumed to be in the order they
appear in the input file (indices start at 0). For the lower level objective, coefficients 
are assumed to be in the order of lower level variables. 

#### Example

Here is the auxiliary file corresponding to the above example:

```
N 1
M 2
LC 1
LR 2
LR 3
LO 1
OS 1
```

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

