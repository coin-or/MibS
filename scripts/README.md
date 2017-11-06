This script is a guidance for conducting the experiments illustrated
[here](http://www.optimization-online.org/DB_FILE/2017/04/5977.pdf).

# Building and Installing MibS

In order to install `MibS`, see
[here](https://github.com/tkralphs/MibS/blob/master/README.md).

The addresses of the directories in which `MibS` is placed and installed
are referred by `<mibs-dir>` and `<build-dir>` respectively in the next sections.

Furthermore, the directory `<mibs-dir>/scripts` includes all the parameter
files and the required scripts. In the next sections, we refer to the
address of this directory by `<scripts-dir>`.

# Test Sets

Three different data sets (171 instances in total) were employed in our
experiments as follows:

1. IBLP-DEN: This set contains 50 instances.
2. IBLP-FIS: This set contains 21 instances.
3. MIBLP-XU: This set contains 100 instances.

# Conducting the experiments

All described problems in previous section should be solved by all 11
methods used in our analyses.

If the directory `build` is placed inside of the directory `MibS`
(i.e., `<build-dir>=<mibs-dir>/build`), to solve all instances with all
methods and plot the figures presented in the paper, run the below commands

```
cd <scripts-dir>/analyze
./run.sh
```

Otherwise, run the below command

```
<scripts-dir>/analyze/run.sh --build-dir=<build-dir> --mibs-dir=<mibs-dir>
```
