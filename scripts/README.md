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
e.g., `conda` or `venv`. The Python dependencies are contained in `<scripts-dir>/analyze/replication/requirements.txt`.
The following commands show how to replicate an environment using `conda`.

```
# Create environment with specific Python version
conda create -n myenv python=3.10.16

# Activate
conda activate myenv

# Install dependencies using pip
python -m pip install -r <scripts-dir>/analyze/replication/requirements.txt
```

# Test Sets

The datasets for these experiments are contained in `<mibs-dir>/data/improvingDirectionDatasets`.

# Conducting the experiments

The file `<scripts-dir>/analyze/replication/myrun.py` contains the parameter configurations for
the computational experiments. It defines all MibS parameter settings (scenarios) and common
parameters (e.g., time limit) required to reproduce the experiments reported in the manuscript.
This file should not need to be modified unless the user wants to test different configurations.

All experiments were conducted on compute nodes running Linux (Debian 8.11) operating system with dual AMD Opteron 6128 processors and 32GB of RAM. All experiments were run sequentially with a time limit of 3600 seconds.
On this setup, all computations were completed in (approximately) 10 days.

To start the experiments, run the following commands from `<mibs-dir>`.

```
cd scripts/analyze
python run.py \
    --binariesPath idbc ../../../build-idBC-opt/bin/mibs \
    --instanceDirs \
        iblpDen    ../../data/improvingDirectionDatasets/iblpDen \
        iblpDen2   ../../data/improvingDirectionDatasets/iblpDen2 \
        iblpZhang  ../../data/improvingDirectionDatasets/iblpZhang \
        iblpZhang2 ../../data/improvingDirectionDatasets/iblpZhang2 \
        iblpFis    ../../data/improvingDirectionDatasets/iblpFis \
        interdDen  ../../data/improvingDirectionDatasets/interdDen \
        interKpShi ../../data/improvingDirectionDatasets/interKpShi \
    --outputDir <results-dir>
```

The `--binariesPath` argument takes `VERSION PATH` pairs and `--instanceDirs` takes `DATASET DIRECTORY` pairs.
The `--outputDir` defaults to `./results` if omitted.
If an output file for an instance already exists and is complete, `run.py` will skip it automatically,
so the script can be safely restarted after an interruption by executing the same command again.

The paths above assume the standard setup described in the "Building and Installing MibS" section:
`coinbrew` run from the parent of `<mibs-dir>`, with build directory `build-idBC-opt`.
If your directory layout differs, adjust `--binariesPath` and `--instanceDirs` accordingly.

After `run.py` completes, `<results-dir>` will have the following structure.

```
<results-dir>/
└── idbc/
    ├── idBC-MILP_fracB/
    │   ├── iblpDen/
    │   │   ├── instance1.out
    │   │   ├── instance1.err
    │   │   └── ...
    │   ├── iblpDen2/
    │   ├── iblpFis/
    │   ├── iblpZhang/
    │   └── iblpZhang2/
    ├── idBC-LS-k_2-dBnd_Inf_fracB/
    │   └── ... (same dataset subdirectories)
    └── ... (one subdirectory per scenario)
```

Once computations are done, run the below command to create all plots presented in the manuscript.

```
python make_plots.py \
    --outputDir <results-dir> \
    --makePlotsImprovingDirectionsPaper
```

All generated plots are saved under `<results-dir>/figures/` and summary CSV files are written to `<results-dir>/`.

The `--makePlotsImprovingDirectionsPaper` flag activates the plotting configuration specific to this
paper. Without it, `make_plots.py` can also be used to generate plots for individual datasets by
passing `--dataSets <dataset-name> [...]` and optionally `--aggregateDatasets` to combine them into the same plots;
add `--help` for further guidance.

# Notes
Before running the full experiments, you can verify that the pipeline works correctly using the small example datasets: `iblpSmall` and `interSmall`. To do this:
1. Run `run.py` passing only those datasets via `--instanceDirs`;
2. Run `make_plots.py` with `--dataSets iblpSmall interSmall` (without `--makePlotsImprovingDirectionsPaper`) to ensure the workflow executes as expected.

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
