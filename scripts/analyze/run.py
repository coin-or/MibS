# Script to replicate experiments of paper 
# "Improving Directions in Mixed Integer Bilevel Linear Optimization"
# Battista F. & Ted K. Ralphs
# Last edited 2026

import os
import subprocess

from parameters import Parameters

# Import myrun
from replication.myrun import mibsParamsInputs, commonParams

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
        paramfile = os.path.join(paramsubpath1, scenario + '.par')
        with open(paramfile, 'w') as file:
            for k, v in params[scenario].items():
                file.write(k + ' ' + v + '\n')


def runExperiments(instPaths, outDir, versions, params, pbsfile=None, testname='mibs'):
    """
        Use to run experiments on local machine or using qsub. 
    """

    # TODO: integrate gaps call for bounded rationality

    # set up output directories
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    # use hierarchy:  outDir/version/param_scenario_name/testset_name/
    for v in versions:
        currpath = os.path.join(outDir, v)
        if not os.path.exists(currpath):
            os.mkdir(currpath)
        for scenario in params:
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
        print(f"Version: {v}")

        exe = binariesPath[v]
        for scenario in params:
            print(f"   Scenario: {scenario}")
            paramsubpath = os.path.join(parampath, scenario)     
            for testset in instPaths:
                paramfile = os.path.join(paramsubpath, scenario + '.par')
                outsubpath = os.path.join(outDir, v, scenario, testset)
                with os.scandir(instPaths[testset]) as inst_it:
                    instances = [
                        instance for instance in inst_it
                        if instance.name.endswith('.mps') or instance.name.endswith('.mps.gz')
                    ]

                total_instances = len(instances)
                print(f"      Dataset: {testset} ({total_instances} instances)")
                for index, instance in enumerate(instances, start=1):
                    remaining = total_instances - index

                    if instance.name.endswith('.mps'):
                        auxname = instance.path[:-4] + ".aux"
                        outname = os.path.join(outsubpath, instance.name[:-4] + '.out')
                        errname = os.path.join(outsubpath, instance.name[:-4] + ".err")
                    elif instance.name.endswith('.mps.gz'):
                        auxname = instance.path[:-7] + ".aux"
                        outname = os.path.join(outsubpath, instance.name[:-7] + '.out')
                        errname = os.path.join(outsubpath, instance.name[:-7] + ".err")
                    else:
                        # Something went wrong. Raise exception and exit.
                        raise ValueError("Unexpected file extension!")

                    print(
                        "         "
                        f"Instance: {instance.name} "
                        f"({index}/{total_instances}, {remaining} remaining)"
                    )

                    # Check if output file already exists and is complete
                    isComplete = False
                    if os.path.exists(outname):
                        with open(outname, 'r') as f:
                            isComplete = "Number of problems (VF) solved" in f.read()
                    
                    isComplete = isComplete and (os.path.isfile(errname) and \
                                           os.path.getsize(errname) == 0)
                    
                    if not isComplete:
                        # prepare command for execution
                        if pbsfile:
                            # qsub 
                            argList = ["qsub", "-v", 
                                            "EXECUTABLE="+exe+","
                                            +"INSTANCENAME="+instance.path+","
                                            +"AUXNAME="+auxname+".aux"]
                        else:
                            # Local
                            argList = [exe,
                                        "-Alps_instance", instance.path,
                                        "-MibS_auxiliaryInfoFile", auxname]

                        # How do we pass parameters?      
                        if writeParams:
                            # Using file previously created
                            if pbsfile:
                                argList[-1] += ",PARAMARG=-param " + paramfile
                            else:
                                argList += ["-param", paramfile]
                        else:
                            # Using command line
                            paramcmd = ' -'.join(' '.join(_) for _ in params[scenario].items())
                            paramcmd = '' + paramcmd
                            if pbsfile:
                                argList[-1] += ",PARAMARG=" + paramcmd
                            else:
                                argList += paramcmd.split()
                        
                        # Run command and redirect output to outfile and errfile
                        if pbsfile:
                            argList += ["-o", outname,
                                        "-e", errname,
                                        "-N", testname,
                                        pbsfile]
                            subprocess.run(argList)
                            print("            Submitted")

                        else: 
                            outfile = open(outname,'w')
                            errfile = open(errname, 'w')
                            subprocess.run(argList, stdout=outfile, stderr=errfile)
                            outfile.close()
                            errfile.close()
                            print("            Complete")
    
                    else:
                        print("            Skipping already complete output")

if __name__ == "__main__":

    # Update mibsParamsInputs with commonParams
    for t in mibsParamsInputs:
        mibsParamsInputs[t].update(commonParams)

    # Read arguments from CLI
    params = Parameters()
    params.parse()

    # Validate required arguments
    if  not params['binariesPath'] or \
        not params['instanceDirs']:
          raise ValueError("Missing required arguments: --binariesPath, --instanceDirs, are required.")
    
    binariesPath = params['binariesPath']
    instanceDirs = params['instanceDirs']
    versions = params['versions']
    outputDir = params['outputDir']
    writeParams = params['writeParams']
    testName = params['testName']
    pbsFile = params['pbsFile']

    ######################### Run Experimests #########################
    runExperiments(instanceDirs, 
                   outputDir, 
                   versions, 
                   mibsParamsInputs, 
                   pbsfile=pbsFile, 
                   testname=testName)
    
