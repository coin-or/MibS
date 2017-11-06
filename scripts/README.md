# Building and Installing MibS

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
BuildTools/get.dependencies install --prefix=/path/to/install/dir
```

After installation, you will also need to add `/your/install/dir/bin` to your
`PATH` variable in your `.bashrc` and also add `/your/install/dir/lib`
to your `LD_LIBRARY_PATH` if you want to link to COIN libraries.

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

All described problems in Section 2 should be solved by all 11
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
