This script is a guidance for conducting the experiments illustrated
[here](https://arxiv.org/abs/2511.03566).

# Building and Installing MibS

In order to download and install this fork of `MibS`, see
[here](https://github.com/febattista/MibS/blob/improving-directions/README.md).

The addresses of the directories in which `MibS` is placed and installed
are referred by `<mibs-dir>` and `<build-dir>` respectively in the next sections.

Furthermore, the directory `<mibs-dir>/scripts` includes all scripts (in Python)
required to run the experiments, parse the output and make the plots. 
In the next sections, we refer to the address of this directory by `<scripts-dir>`.

# Setting a Python environment
To run the scripts, it is strongly recommended using Python 3.10 under virtual environment,
e.g., `conda` or `venv`. The Python dependencies are contained in `<scripts-dir>/analyze/requirements.txt`.
The following commands show how to replicate an environment using `conda`.

```
# Create environment with specific Python version
conda create -n myenv python=3.10

# Activate
conda activate myenv

# Install dependencies using pip
pip install -r <scripts-dir>/analyze/requirements.txt
```

# Test Sets

The datasets for this experiments are contained in `<mibs-dir>/data/improvingDirectionDatasets`.

# Conducting the experiments

The file `<scripts-dir>/analyze/myrun.py` is used to configure and manage the computational experiments. In particular, it allows you to specify:

- The datasets to be used and their corresponding locations;
- The path to the MibS executable;
- The configurations to be tested (i.e., different parameter settings of MibS to compare),
- The output directory where all raw log files will be stored. 

By default, the file is preconfigured with relative paths consistent with the setup instructions provided above. 
It includes all datasets and parameter configurations required to reproduce the experiments reported in the manuscript. 
Please, note that:
1. The path of `mibsDir` must be changed to match `<mibs-dir>`;
2. If your directory structure differs from the one described in the "BUILDING from source" instructions, please update the paths accordingly.

All experiments were conducted on compute nodes running Linux (Debian 8.11) operating system with dual AMD Opteron 6128 processors and 32GB of RAM. All experiments were run sequentially with a time limit of 3600 seconds.
On this setup, all computations were completed in (approximately) 10 days.

To start the experiments, run the following command.

```
cd <scripts-dir>/analyze
python run.py
```

Once computations are done, run the below command to create all plots.

```
python process.py
```
All generated plots and CSV can be found in the current directory.

By default, `<scripts-dir>/analyze/process.py` is configured to reproduce all plots presented in the manuscript. It assumes that all datasets have been solved under every configuration specified in the experimental setup. If this is not the case, the script must be adjusted accordingly (see the inline comments for guidance).

# Notes
Before running the full experiments, you can verify that the workflow works correctly using the small example datasets: `iblpSmall` and `interSmall`. To do this:
1. Modify `myrun.py` and `process.py` accordingly to use only those datasets;
2. Run the scripts to ensure the workflow executes as expected.

This helps catch any issues early before processing the full datasets.

# Troubleshooting
If executing any of the installed binaries results in an error that shared libraries cannot be found, you may need to add
```
# Linux
export LD_LIBRARY_PATH=<build-dir>/lib 
```
or
```
# OS X
export DYLD_LIBRARY_PATH=<build-dir>/lib 
```
to your `~/.bashrc`.
