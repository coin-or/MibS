# Script to replicate experiments of paper 
# "Improving Directions in Mixed Integer Bilevel Linear Optimization"
# Battista F. & Ted K. Ralphs
# Last edited 2026

import sys, os, collections
import shutil, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import myrun
from myrun import instanceDirs, outputDir, mibsParamsInputs, pbsfile, testname, versions, commonParams, binariesPath, writeParams

def writeParamsToFile(outDir, params):
    """
        Write MibS parameter for each scenario into a file
    """
    parampath = os.path.join(outDir, 'parameters')
    if not os.path.exists(parampath):
        os.mkdir(parampath)
    
    for scenario in params:
        # for this scenario, make a parameter file and save to specified directory
        paramsubpath1 = os.path.join(parampath, scenario)
        if not os.path.exists(paramsubpath1):
            os.mkdir(paramsubpath1)
        os.chdir(paramsubpath1)
        file = open(scenario+'.par', 'w')
        for k, v in params[scenario].items():
            file.write(k + ' ' + v + '\n')
        file.close()


def runExperiments(instPaths, outDir, versions, params):
    """
        Use to run experiments on local machine. 
    """

    # set up output directories
    # use hierarchy:  outDir/version/param_scenario_name/testset_name/
    for v in versions:
        currpath = os.path.join(outDir, v)
        if not os.path.exists(currpath):
            os.mkdir(currpath)
        for scenario in params:
            print(scenario)
            currpath = os.path.join(outDir, v, scenario)
            if not os.path.exists(currpath):
                os.mkdir(currpath)
            for testset in instPaths:
                currsubpath1 = os.path.join(currpath, testset)
                if not os.path.exists(currsubpath1):
                    os.mkdir(currsubpath1)
   
    # if choose to write params into files
    if writeParams:
       writeParamsToFile(outDir, params)
    
    # Where parameter files are
    parampath = os.path.join(outDir, 'parameters')

    for v in versions:
        exe = binariesPath[v]
        for scenario in params:
            paramsubpath = os.path.join(parampath, scenario)     
            for testset in instPaths:
                paramfile = os.path.join(paramsubpath, scenario + '.par')
                outsubpath = os.path.join(outDir, v, scenario, testset)
                os.chdir(outsubpath)
                with os.scandir(instPaths[testset]) as inst_it: 
                    for instance in inst_it:
                        # Run experiments only on .mps files
                        if instance.name.endswith('.mps') or \
                            instance.name.endswith('.mps.gz'):

                            if instance.name.endswith('.mps'):
                                outname = instance.name[:-4] + '.out'
                                auxname = instance.path[:-4] + ".aux"
                            elif instance.name.endswith('.mps.gz'):
                                outname = instance.name[:-7] + '.out'
                                auxname = instance.path[:-7] + ".aux"
                            else:
                                # Something went wrong
                                print("Something went wrong with file extension!")
                                print(instance.name)
                                exit(1)

                            outfile = open(outname,'w')
                            argList = [exe,
                                        "-Alps_instance", instance.path,
                                        "-MibS_auxiliaryInfoFile", auxname]

                            # How do we pass parameters?      
                            if writeParams:
                                # Using file previously created
                                argList += ["-param", paramfile]
                            else:
                                # Using command line
                                paramcmd = ' -'.join(' '.join(_) for _ in params[scenario].items())
                                paramcmd = '' + paramcmd
                                argList += paramcmd.split()
                            
                            # Run command and redirect output to outfile
                            # subprocess.run(argList, stdout=outfile)
                            outfile.close()
                            print('Complete {}'.format(instance.name))                 

if __name__ == "__main__":

    for t in mibsParamsInputs:
        mibsParamsInputs[t].update(commonParams)

    ######################### Run Experimests #########################
    # local: provide paths in myrun.py
    runExperiments(instanceDirs, outputDir, versions, mibsParamsInputs)
    


