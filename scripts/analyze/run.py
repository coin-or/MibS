# Script to run MibS with different SL target gap.
# The script may produce auxiliary folder/files in run directory. 
# Last edited by yux616
# Apr 2020 
# Script path:  /MibS/scripts/analyze
# Some os function requires Python 3.5+

# add arg parser later

import sys, os, collections
import shutil, subprocess
# from types import FrameType
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib 
#import itertools

# from runparams import exe, instanceDirs, outputDir, mibsParamsInputs, pbsfile
# from runparams_cuts import exe, instanceDirs, outputDir, mibsParamsInputs, pbsfile
from myrun import instanceDirs, outputDir, mibsParamsInputs, pbsfile, testname, versions, commonParams

def writeParamsToFile(outDir, params, gaps=[]):
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
        if gaps:
            for g in gaps:
                src = scenario + '.par'
                dst = scenario + '_g'+ str(g) +'.par'
                shutil.copyfile(src, dst)
                file = open(dst, 'a')
                file.write('MibS_slTargetGap ' + str(g) + '\n')
                file.close()


def runExperiments(exe, instPaths, outDir, versions, params, gaps=[]):
    """
        Use to run experiments on local machine. 
    """

    writeParams = False
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
                if gaps: 
                    for g in gaps:
                        currsubpath2 = os.path.join(currsubpath1,'BR' + str(g) + 'Output')
                        if not os.path.exists(currsubpath2):
                            os.mkdir(currsubpath2)
   
    # if choose to write params into files
    if writeParams:
        cwd = os.getcwd()
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
            if gaps:
                for g in gaps:
                    src = scenario + '.par'
                    dst = scenario + '_g'+ str(g) +'.par'
                    shutil.copyfile(src, dst)
                    file = open(dst, 'a')
                    file.write('MibS_slTargetGap ' + str(g) + '\n')
                    file.close()
        
            # run experiments use modified parameter file      
            for testset in instPaths:
                if gaps: 
                    for g in gaps:
                        paramfile = os.path.join(paramsubpath1, scenario + '_g'+ str(g) +'.par')
                        outsubpath = os.path.join(outDir, scenario, testset,'BR' + str(g) + 'Output')
                        os.chdir(outsubpath)   
                        with os.scandir(instPaths[testset]) as inst_it: 
                            for instance in inst_it:
                                if instance.name.endswith('.mps'):
                                    outname = instance.name[:-4]+'.out'
                                    outfile = open(outname,'w')
                                    subprocess.run([exe, "-param", paramfile,
                                                    "-Alps_instance", instance.path,
                                                    "-MibS_auxiliaryInfoFile", instance.path[:-4]+".aux",
                                                    ],
                                                    stdout=outfile)
                                    outfile.close()
                                    print('Complete {} with gap {}'.format(instance.name, g))
                else:
                    paramfile = os.path.join(paramsubpath1, scenario + '.par')
                    # print(paramcmd)
                    v = versions[0]
                    outsubpath = os.path.join(outDir, v, scenario, testset)
                    # print(outsubpath)
                    os.chdir(outsubpath)
                    with os.scandir(instPaths[testset]) as inst_it: 
                        for instance in inst_it:
                            # print(instance.name)
                            if instance.name.endswith('.mps'):
                                outname = instance.name[:-4]+'.out'
                                outfile = open(outname,'w')
                                subprocess.run([exe, "-param", paramfile,
                                                "-Alps_instance", instance.path,
                                                "-MibS_auxiliaryInfoFile", instance.path[:-4]+".aux",
                                            ],
                                                stdout=outfile)
                                outfile.close()
                                print('Complete {}'.format(instance.name))
                            elif instance.name.endswith('.mps.gz'):
                                outname = instance.name[:-7]+'.out'
                                outfile = open(outname,'w')
                                subprocess.run([exe, "-param", paramfile,
                                                "-Alps_instance", instance.path,
                                                "-MibS_auxiliaryInfoFile", instance.path[:-7]+".aux ",
                                            ],
                                                stdout=outfile)
                                outfile.close()
                                print('Complete {}'.format(instance.name))
        # remove paramter files created earlier?   
    else:
        # run experiments use command line paramters
        for v in versions:
            for testset in instPaths:
                for scenario in params:
                    if gaps: 
                        for g in gaps:
                            paramcmd = ' -'.join(' '.join(_) for _ in params[scenario].items())
                            paramcmd = '-' + paramcmd + ' -MibS_slTargetGap ' + str(g)
                            # print(paramcmd)
                            outsubpath = os.path.join(outDir, v, scenario, testset,'BR' + str(g) + 'Output')
                            os.chdir(outsubpath)   
                            with os.scandir(instPaths[testset]) as inst_it: 
                                for instance in inst_it:
                                    if instance.name.endswith('.mps'):
                                        outname = instance.name[:-4]+'.out'
                                        outfile = open(outname,'w')
                                        subprocess.run([exe,
                                                        "-Alps_instance", instance.path,
                                                        "-MibS_auxiliaryInfoFile", instance.path[:-4]+".aux",
                                                        paramcmd,
                                                    ],
                                                       stdout=outfile)
                                        outfile.close()
                                        print('Complete {} with gap {}'.format(instance.name, g))
                    else:
                        if scenario == 'MIBS':
                            exe = '/home/federico/Scrivania/coin/improvingDir/build-1.2-opt/bin/mibs'
                        if scenario == 'idB&C-IDIC':
                            exe = '/home/federico/Scrivania/coin/improvingDir/build-ipco-opt/bin/mibs'
                        print(exe)
                        paramcmd = ' -'.join(' '.join(_) for _ in params[scenario].items())
                        paramcmd = '' + paramcmd
                        # print(paramcmd)
                        outsubpath = os.path.join(outDir, v, scenario, testset)
                        os.chdir(outsubpath)
                        with os.scandir(instPaths[testset]) as inst_it: 
                            for instance in inst_it:
                                # print(instance.name)
                                if instance.name.endswith('.mps'):
                                    outname = instance.name[:-4]+'.out'
                                    outfile = open(outname,'w')
                                    subprocess.run([exe,
                                                    "-Alps_instance", instance.path,
                                                    "-MibS_auxiliaryInfoFile", instance.path[:-4]+".aux",
                                                ]+paramcmd.split(),
                                                    stdout=outfile)
                                    outfile.close()
                                    print('Complete {}'.format(instance.name))
                                elif instance.name.endswith('.mps.gz'):
                                    outname = instance.name[:-7]+'.out'
                                    outfile = open(outname,'w')
                                    subprocess.run([exe,
                                                    "-Alps_instance", instance.path,
                                                    "-MibS_auxiliaryInfoFile", instance.path[:-7]+".aux ",
                                                ]+paramcmd.split(),
                                                    stdout=outfile)
                                    outfile.close()
                                    print('Complete {}'.format(instance.name))

    return

def runExperimentsPBS(instPaths, outDir, versions, params, pbsfile, gaps=[]):
    """
        Use to submit batch jobs via qsub. 
    """

    total_jobs = 0
    # set up output directories
    # use hierarchy:  outDir/version/param_scenario_name/testset_name
    for v in versions:
        for scenario in params:
            for testset in instPaths:
                if gaps: 
                    for g in gaps:
                        currsubpath = os.path.join(outDir, v, scenario, testset, 'BR' + str(g) + 'Output')
                        if not os.path.exists(currsubpath):
                            os.makedirs(currsubpath)
                        else:
                            shutil.rmtree(currsubpath, ignore_errors=True)
                            os.makedirs(currsubpath, exist_ok=True)      
                else:
                    currsubpath = os.path.join(outDir, v, scenario, testset)
                    if not os.path.exists(currsubpath):
                        os.makedirs(currsubpath)
                    else:
                        print('Directory already exists! Ignoring...')
                        #sys.exit()
                        #shutil.rmtree(currsubpath, ignore_errors=True)
                        #os.makedirs(currsubpath, exist_ok=True)  

    # submit experiments use command line paramters
    # see .pbs file anc comments for job submission arguments
    for v in versions:
        if v == 'ib':
            exe = '/home/ted/Projects/build-mibs-debug-no-cplex/bin/mibs'
        if v == '1.1':
            exe = '/home/ted/tmp/COIN/build-mibs-1.1-debug/bin/mibs'
        if v == '1.2-opt':
            #exe = '/home/ted/Projects/build-mibs-1.2/bin/mibs'
            exe = '/home/ted/Projects/build-mibs-cplex-1.2/bin/mibs'
        if v == 'pr-92':
            exe = '/home/ted/Projects/build-mibs-pr-92/bin/mibs'
        if v == 'improvingDir':
            exe = '/home/feb223/improvingDir/build-MibS-opt/bin/mibs'
        if v == 'idBC':
            exe = '/home/feb223/improvingDir/build-idBC-opt/bin/mibs'
        if 'idBC-debug' in v:
            exe = '/home/feb223/improvingDir/build-idBC-opt/bin/mibs'
        if v == '1.2':
            exe = '/home/feb223/improvingDir/build-1.2-opt/bin/mibs'
                
        for testset in instPaths:
            for scenario in params:
                paramcmd = '-' + ' -'.join(' '.join(_) for _ in params[scenario].items())
                # print(paramcmd)
                if gaps: 
                    for g in gaps:
                        paramcmd_g = '-' + paramcmd + ' -MibS_slTargetGap ' + str(g)
                        # print(paramcmd)
                        outsubpath = os.path.join(outDir, v, scenario, testset,'BR' + str(g) + 'Output')
                        os.chdir(outsubpath)   
                        with os.scandir(instPaths[testset]) as inst_it: 
                            for instance in inst_it:
                                if instance.name.endswith('.mps'):
                                    outfile = os.path.join(outsubpath, instance.name[:-4]+'.out')
                                    errfile = os.path.join(outsubpath, instance.name[:-4]+'.err')
                                    subprocess.run(["qsub", "-v", 
                                                    "EXECUTABLE="+exe+","
                                                    +"INSTANCENAME="+instance.path+","
                                                    +"AUXNAME="+instance.path[:-4]+".aux"+","
                                                    +"PARAMARG="+paramcmd_g,
                                                    "-o", outfile,
                                                    "-e", errfile,
                                                    "-N", testname,
                                                    pbsfile])
                else:
                    outsubpath = os.path.join(outDir, v, scenario, testset)
                    os.chdir(outsubpath)   
                    with os.scandir(instPaths[testset]) as inst_it: 
                        for instance in inst_it:
                            if instance.name.endswith('.mps'):
                                outfile = os.path.join(outsubpath, instance.name[:-4]+'.out')
                                errfile = os.path.join(outsubpath, instance.name[:-4]+'.err')
                                isComplete = False
                                if os.path.exists(outfile):
                                    with open(outfile, 'r') as f:
                                        isComplete = "Number of problems (VF) solved" in f.read()
                                isOutFileIncomplete = (os.path.isfile(errfile) and os.path.getsize(errfile) == 0 and not isComplete)
                                if not os.path.exists(outfile) or isOutFileIncomplete:
                                    subprocess.run(["qsub", "-v", 
                                                    "EXECUTABLE="+exe+","
                                                    +"INSTANCENAME="+instance.path+","
                                                    +"AUXNAME="+instance.path[:-4]+".aux"+","
                                                    +"PARAMARG="+paramcmd,
                                                    "-o", outfile,
                                                    "-e", errfile,
                                                    "-N", testname,
                                                    pbsfile])
                                    total_jobs += 1
                                    # exit(0)

                                else:
                                    print("File", outfile, "exists!")
                            elif instance.name.endswith('.mps.gz'):
                                outfile = os.path.join(outsubpath, instance.name[:-7]+'.out')
                                errfile = os.path.join(outsubpath, instance.name[:-7]+'.err')
                                isComplete = False
                                if os.path.exists(outfile):
                                    with open(outfile, 'r') as f:
                                        isComplete = "Number of problems (VF) solved" in f.read()
                                isOutFileIncomplete = (os.path.isfile(errfile) and os.path.getsize(errfile) == 0 and not isComplete)
                                if not os.path.exists(outfile) or isOutFileIncomplete:
                                    subprocess.run(["qsub", "-v", 
                                                    "EXECUTABLE="+exe+","
                                                    +"INSTANCENAME="+instance.path+","
                                                    +"AUXNAME="+instance.path[:-7]+".aux"+","
                                                    +"PARAMARG="+paramcmd,
                                                    "-o", outfile,
                                                    "-e", errfile,
                                                    "-N", testname,
                                                    pbsfile])
                                    total_jobs += 1
                                    # exit(0)
                                    
                                else:
                                    print("File", outfile, "exists!")
    print("")
    print("Launched ", total_jobs, "jobs!")
    print("")
    return                    

if __name__ == "__main__":

    for t in mibsParamsInputs:
        mibsParamsInputs[t].update(commonParams)

    ######################### Run Experimests #########################
    # local: provide paths in runparams.py
    # exe = '/home/federico/Scrivania/coin/improvingDir/build-ipco-opt/bin/mibs'
    # exe = '/home/federico/Scrivania/coin/improvingDir/build-1.2-opt/bin/mibs'
    # runExperiments(exe, instanceDirs, outputDir, versions, mibsParamsInputs)
    
    # using pbs file: provide paths in runparams.py
    writeParamsToFile(outputDir, mibsParamsInputs)
    runExperimentsPBS(instanceDirs, outputDir, versions, mibsParamsInputs, pbsfile)
    


