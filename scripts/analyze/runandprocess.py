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
from runparams_cuts import exe, instanceDirs, outputDir, mibsParamsInputs, pbsfile

    
def runExperiments(exe, instPaths, outDir, params, gaps=[]):
    """
        Use to run experiments on local machine. 
    """
    writeParams = False

    # set up output directories
    # use hierarchy:  outDir/param_scenario_name/testset_name/BR_Output/
    for scenario in params:
        # print(scenario)
        currpath = os.path.join(outDir, scenario)
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
        parampath = os.path.join(cwd, '../parameters')
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
                    pass # add alternative call information
        # remove paramter files created earlier?   
    else:
        # run experiments use command line paramters
        for scenario in params:
            for testset in instPaths:
                if gaps: 
                    for g in gaps:
                        paramcmd = ' -'.join(' '.join(_) for _ in params[scenario].items())
                        paramcmd = '-' + paramcmd + ' -MibS_slTargetGap ' + str(g)
                        # print(paramcmd)
                        outsubpath = os.path.join(outDir, scenario, testset,'BR' + str(g) + 'Output')
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
                    pass # add alternative call information
    return

def runExperimentsPBS(exe, instPaths, outDir, params, pbsfile, gaps=[]):
    """
        Use to submit batch jobs via qsub. 
    """
    # set up output directories
    # use hierarchy:  outDir/param_scenario_name/testset_name/BR_Output/
    for scenario in params:
        for testset in instPaths:
            if gaps: 
                for g in gaps:
                    currsubpath = os.path.join(outDir, scenario, testset, 'BR' + str(g) + 'Output')
                    if not os.path.exists(currsubpath):
                        os.makedirs(currsubpath)
                    else:
                        shutil.rmtree(currsubpath, ignore_errors=True)
                        os.makedirs(currsubpath, exist_ok=True)      
            else:
                currsubpath = os.path.join(outDir, scenario, testset)
                if not os.path.exists(currsubpath):
                    os.makedirs(currsubpath)
                else:
                    shutil.rmtree(currsubpath, ignore_errors=True)
                    os.makedirs(currsubpath, exist_ok=True)  

    # submit experiments use command line paramters
    # see .pbs file anc comments for job submission arguments
    for scenario in params:
        for testset in instPaths:
            paramcmd = ' -'.join(' '.join(_) for _ in params[scenario].items())
            if gaps: 
                for g in gaps:
                    paramcmd_g = '-' + paramcmd + ' -MibS_slTargetGap ' + str(g)
                    # print(paramcmd)
                    outsubpath = os.path.join(outDir, scenario, testset,'BR' + str(g) + 'Output')
                    os.chdir(outsubpath)   
                    with os.scandir(instPaths[testset]) as inst_it: 
                        for instance in inst_it:
                            if instance.name.endswith('.mps'):
                                outfile = os.path.join(outsubpath, instance.name[:-4]+'.out')
                                errfile = os.path.join(outsubpath, instance.name[:-4]+'.err')
                                # pbs args: ARGS="-Alps_instance $INSTANCENAME -MibS_auxiliaryInfoFile $AUXNAME $PARAMARG"; $EXECUTABLE $ARGS;
                                subprocess.run(["qsub", "-v", 
                                                "EXECUTABLE="+exe+","
                                                +"INSTANCENAME="+instance.path+","
                                                +"AUXNAME="+instance.path[:-4]+".aux"+","
                                                +"PARAMARG="+paramcmd_g,
                                                "-o", outfile,
                                                "-e", errfile,
                                                pbsfile])
            else:
                outsubpath = os.path.join(outDir, scenario, testset)
                os.chdir(outsubpath)   
                with os.scandir(instPaths[testset]) as inst_it: 
                    for instance in inst_it:
                        if instance.name.endswith('.mps'):
                            outfile = os.path.join(outsubpath, instance.name[:-4]+'.out')
                            errfile = os.path.join(outsubpath, instance.name[:-4]+'.err')
                            # pbs args: ARGS="-Alps_instance $INSTANCENAME -MibS_auxiliaryInfoFile $AUXNAME $PARAMARG"; $EXECUTABLE $ARGS;
                            subprocess.run(["qsub", "-v", 
                                            "EXECUTABLE="+exe+","
                                            +"INSTANCENAME="+instance.path+","
                                            +"AUXNAME="+instance.path[:-4]+".aux"+","
                                            +"PARAMARG="+paramcmd,
                                            "-o", outfile,
                                            "-e", errfile,
                                            pbsfile])
    return                    

def parseOutput(outputDir, params, gaps=[], writeCSV=True, filename='summary.csv'):
    """
        The function parse the output file in the given directory.
        Assume the subfolders hierarchy: outputDir/param_scenario_name/testset_name/BR_Output/file.out.
        The result will also be written to a .csv file if not specified.
        Note: currently not able to read incomplete output due to external interruption.
        Input: 
            outputDir: string, a path to the parent output directory
            writeCSV: boolean, whether to save the results in structured format
        Return:
            a pandas dataframe containing results
    """
    # may move those up as input in the future...
    # then need to match keywords and fields
    keywords = [
        'No solution found',
        'of nodes processed',
        # 'branched',
        'Search CPU time',
        'VF) solved',
        'UB) solved',
        'solving problem (VF)',
        'solving problem (UB)',
        'Cost',
        'fully processed',
        'gap is',
        'integer UL Variables',
        'integer LL Variables',
        ]

    # fields = [
    #     'given_SL_gap',
    #     'dataset',
    #     'instance',
    #     'nodes_f_proc',
    #     'cpu_time',
    #     'vf_solved',
    #     'ub_solved',
    #     'vf_time',
    #     'ub_time',
    #     'objval_cost',
    #     'prob_gap'
    #     ]

    results = collections.defaultdict(list)
    etol = np.finfo(float).eps 

    # iterate over scenarios, datasets, gaps, and files in each folder to read results
    # for g in gaps:
    for scenario in params:
        # resultDir = os.path.join(outputDir, "BR"+str(g)+"Output")
        resultDir = os.path.join(outputDir, scenario)

        # iterate over different datasets available
        with os.scandir(resultDir) as dataset_it:
            for d_entry in dataset_it:
                if gaps:
                    for g in gaps:
                        subDir =  os.path.join(d_entry.path, 'BR'+str(g)+'Output')
                        # iterate over files in the folder
                        with os.scandir(subDir) as output_it: 
                            for o_entry in output_it:
                                if o_entry.name.endswith('.out'):
                                    # start to write result to the dictionary
                                    results['given_gap'].append(g)
                                    results['dataset'].append(d_entry.name)
                                    results['scenario'].append(scenario)
                                    results['instance'].append(o_entry.name.split('.')[0])
                                    
                                    incomplete = True # mark incomplete output file
                                    nosoln = False # mark no soluntion found by Alps at termination
                                    # read value for each field from file
                                    with open(o_entry.path,'r') as file:
                                        for line in file:
                                            if keywords[0] in line:
                                                nosoln = True
                                                print("No solution found instance:", o_entry.name)
                                            
                                            elif keywords[1] in line or keywords[8] in line:
                                                results['nodes_proc'].append(int(line.split(':')[1]))
                                            
                                            elif keywords[2] in line:
                                                results['cpu_time'].append(float( (line.split(':')[1]).split()[0] ))
                                            
                                            elif keywords[3] in line:
                                                results['vf_solved'].append(int(line.split('=')[1]))  
                                                            
                                            elif keywords[4] in line:
                                                results['ub_solved'].append(int(line.split('=')[1]))
                                            
                                            elif keywords[5] in line:
                                                results['vf_time'].append(float(line.split('=')[1]))
                                            
                                            elif keywords[6] in line:
                                                results['ub_time'].append(float(line.split('=')[1]))
                                            
                                            elif keywords[7] in line:
                                                results['objval_cost'].append(int(line.split('=')[1]))
                                            
                                            elif keywords[9] in line:
                                                incomplete = False
                                                if 'infinity' in line:
                                                    results['prob_gap'].append(1000000) # no soln found
                                                    results['solved'].append(False)
                                                else:
                                                    solgap = float(line.split(' ')[5].strip('%\n'))
                                                    results['prob_gap'].append(solgap)
                                                    # mark unsolved instances in given time limit
                                                    if (solgap - 0.0 < etol):
                                                        results['solved'].append(True)
                                                    else:
                                                        results['solved'].append(False)    
                                            elif keywords[10] in line:
                                                results['ul_int_var'].append(int(line.split(':')[1]))
                                            
                                            elif keywords[11] in line:
                                                results['ll_int_var'].append(int(line.split(':')[1]))
                                            else:
                                                pass
                                    if incomplete:
                                        results['nodes_proc'].append(-1)
                                        results['cpu_time'].append(-1)
                                        results['vf_solved'].append(-1)
                                        results['ub_solved'].append(-1)
                                        results['vf_time'].append(-1)
                                        results['ub_time'].append(-1)
                                        results['objval_cost'].append(-1000000)
                                        results['prob_gap'].append(-1)
                                        results['solved'].append(False)
                                    elif nosoln:
                                        results['vf_solved'].append(-1)
                                        results['ub_solved'].append(-1)
                                        results['vf_time'].append(-1)
                                        results['ub_time'].append(-1)
                                        results['objval_cost'].append(-1000000)    
                else:
                    subDir = d_entry.path
                    # iterate over files in the folder
                    with os.scandir(subDir) as output_it: 
                        for o_entry in output_it:
                            if o_entry.name.endswith('.out'):
                                # start to write result to the dictionary
                                results['dataset'].append(d_entry.name)
                                results['scenario'].append(scenario)
                                results['instance'].append(o_entry.name.split('.')[0])
                                
                                incomplete = True # mark incomplete output file
                                nosoln = False # mark no soluntion found by Alps at termination
                                # read value for each field from file
                                with open(o_entry.path,'r') as file:
                                    for line in file:
                                        if keywords[0] in line:
                                            nosoln = True
                                            print("No solution found instance:", o_entry.name)
                                        
                                        elif keywords[1] in line or keywords[8] in line:
                                            results['nodes_proc'].append(int(line.split(':')[1]))
                                        
                                        elif keywords[2] in line:
                                            results['cpu_time'].append(float( (line.split(':')[1]).split()[0] ))
                                        
                                        elif keywords[3] in line:
                                            results['vf_solved'].append(int(line.split('=')[1]))  
                                                        
                                        elif keywords[4] in line:
                                            results['ub_solved'].append(int(line.split('=')[1]))
                                        
                                        elif keywords[5] in line:
                                            results['vf_time'].append(float(line.split('=')[1]))
                                        
                                        elif keywords[6] in line:
                                            results['ub_time'].append(float(line.split('=')[1]))
                                        
                                        elif keywords[7] in line:
                                            results['objval_cost'].append(int(line.split('=')[1]))
                                        
                                        elif keywords[9] in line:
                                            incomplete = False
                                            if 'infinity' in line:
                                                results['prob_gap'].append(1000000) # no soln found
                                                results['solved'].append(False)
                                            else:
                                                solgap = float(line.split(' ')[5].strip('%\n'))
                                                results['prob_gap'].append(solgap)
                                                # mark unsolved instances in given time limit
                                                if (solgap - 0.0 < etol):
                                                    results['solved'].append(True)
                                                else:
                                                    results['solved'].append(False)    
                                        elif keywords[10] in line:
                                            results['ul_int_var'].append(int(line.split(':')[1]))
                                        
                                        elif keywords[11] in line:
                                            results['ll_int_var'].append(int(line.split(':')[1]))
                                        else:
                                            pass
                                if incomplete:
                                    results['nodes_proc'].append(-1)
                                    results['cpu_time'].append(-1)
                                    results['vf_solved'].append(-1)
                                    results['ub_solved'].append(-1)
                                    results['vf_time'].append(-1)
                                    results['ub_time'].append(-1)
                                    results['objval_cost'].append(-1000000)
                                    results['prob_gap'].append(-1)
                                    results['solved'].append(False)
                                elif nosoln:
                                    results['vf_solved'].append(-1)
                                    results['ub_solved'].append(-1)
                                    results['vf_time'].append(-1)
                                    results['ub_time'].append(-1)
                                    results['objval_cost'].append(-1000000)             
    
    # for k in results:
    #     print(len(results[k]))
    df_result = pd.DataFrame(results)
    
    # make some adjustment to formats
    # display check feasibility time as % of search time? 
    # sum vf+ub time -> feasibility time (or read from output directly?)
    df_result['chk_feas_time'] = df_result['ub_time'] + df_result['vf_time']
    df_result['chk_feas_time'] = df_result['chk_feas_time'].astype(float).round(2)
    df_result['cpu_time'] = df_result['cpu_time'].astype(float).round(2)
     
    # write results to .csv file
    if writeCSV: 
        #df_result.to_csv(filename, mode='a', header=False, index=False) # append results only
        df_result.to_csv(filename, index=False)
    
    return df_result

def processTable(df, displayCols, gaps=[], writeLTX=False, filename='ltx_tb.txt'):
    """
        Print a summary table for required columns and gaps.
        Input:
            df: a dataframe with all info from parseOutput
            gaps: a list of int
            displayCol: columns to print
    """

    # Draft:
    # separate instance to different tables by given_gap
    # convert each instance related data into a dictionary 
    # then in each dictionaries, there is another dict of data for different gaps
    # each data field can print to a table; 
    # or print a summary table where instance by row, gaps by column

    # obtain the list of instances
    instList = list(df.instance.unique())
    scnList = list(df.scenario.unique())
    print(instList)

    # collect required info into dict
    rsltDict = {}
    if gaps: # Todo: separate datasets
        for inst in instList:
            rsltDict[inst] = {}
            for g in gaps:
                for scn in scnList:
                    cond = (df['scenario'] == scn) & (df['given_gap'] == g) & (df['instance'] == inst) 
                    df_temp = df[cond]
                    if df_temp['solved'].values[0]:
                        # using (outer key, inner key) for table column formatting
                        # rsltDict[inst].update({(g, col):df_temp[col].values[0] for col in displayCols})
                        rsltDict[inst].update({(scn, g, col):df_temp[col].values[0] for col in displayCols})
                    else:
                        # special color for unsolved cases? but still need to process after print_to_latex
                        rsltDict[inst].update({(scn, g, col):'textcolor{lightgray}{'+str(df_temp[col].values[0])+'}' for col in displayCols})
            
            # for debug print
            # inst == 'miblp_20_20_50_0110_10_1' and print(rsltDict[inst])
    else:
        for inst in instList:
            rsltDict[inst] = {}
            for scn in scnList:
                cond = (df['scenario'] == scn) & (df['instance'] == inst) 
                df_temp = df[cond]
                ds = df_temp['dataset'].values[0]
                if df_temp['solved'].values[0]:
                    # using (outer key, inner key) for table column formatting
                    rsltDict[inst].update({(scn, ds, col):df_temp[col].values[0] for col in displayCols})
                else:
                    # special color for unsolved cases? but still need to process after print_to_latex
                    rsltDict[inst].update({(scn, ds, col):'textcolor{lightgray}{'+str(df_temp[col].values[0])+'}' for col in displayCols})
        
    
    # convert dict to structured df: change to formal column names?
    df_forprint = pd.DataFrame.from_dict(rsltDict, orient='index')
    if gaps:
        df_forprint.columns.names = ['scn', 'gaps', 'fields']
    else:
        df_forprint.columns.names = ['scn', 'datasets', 'fields']
    df_forprint = df_forprint.sort_index()

    # OPTION 1: print results to a single table: suggest to use when display col number < 2 
    # with open('ltx_tb1.txt', 'w') as file:
    #     file.write(df_forprint.to_latex())
    
    # OPTION 2: for each displayCol, print a table; using slicer indexing
    if writeLTX:
        with open(filename, 'w') as file:
            for col in displayCols:
                for scn in scnList:
                    if gaps:
                        file.write(df_forprint.loc[:, (scn, slice(None), col)].to_latex()) # indices may change after add a column
                    else:
                        file.write(df_forprint.loc[:, (scn, slice(None), col)].to_latex())

    # OPTION 3: just process table, do not print latex table to file
    # pass

    return df_forprint

def plotBarChart(df, gaps, plotCols, plotScns):
    """
        Make some plots for given gaps and columns for quick observation.
        Use the measures for the smallest gap as 1.
        #plot = len(plotCol)* len(plotScns)
        Input:
            df: pandas dataframe output from processTable
            gaps: a list of int
            plotCol: columns to make single plots
            plotScns: scenarios 
    """
    # markers = itertools.cycle(('X', '+', '.', 'o', '*')) 
    # colors = itertools.cycle(('dodgerblue', 'slateblue', 'blueviolet', 'palevioletred','lightcoral',
        # 'sandybrown', 'gold', 'yellowgreen', 'darkturquoise')) #'orange',
    colors = ['dodgerblue', 'slateblue', 'blueviolet', 'palevioletred','lightcoral',
        'sandybrown', 'gold', 'yellowgreen', 'darkturquoise']

    for col in plotCols:
        for scn in plotScns:
            # get columns to print using cross-section
            df_sub = df.xs((scn, col), level=[0, 'fields'], axis=1, drop_level=True)
            # df_sub.columns = df_sub.columns.to_flat_index()

            plot_series = []
            df_sub['base'] = pd.to_numeric(df_sub[gaps[0]], errors='coerce').replace(np.nan, -0.01)

            for g in gaps:
                # handle the unsolved cases
                df_sub[g] = df_sub[g].apply(lambda x: 100 if isinstance(x, str) else x)
                plot_series.append((df_sub[g]+0.0001)/(df_sub['base']+0.0001))

            # begin plotting
            fig, ax = plt.subplots(1,1, figsize=(30,5)) #figsize=(7.5,5)
            # obtain xlocs and shifts
            xloc = np.arange(len(df_sub.index))
            mid = int(len(gaps) / 2)
            bar_width = 1/(len(gaps)+1)
            for i,g in enumerate(gaps):
                ax.bar(xloc+(i-mid)*(bar_width*0.6), plot_series[i], 
                    label='TargetGap ='+str(g), color=colors[i], width=bar_width, alpha=0.5)
            
            # set tickers and lims
            ax.set_xticks(xloc)
            ax.set_xticklabels(df_sub.index, rotation=90, fontsize='x-small')
            l_xlim = xloc[0] - (mid+1)*bar_width
            r_xlim = xloc[-1] + (mid+1)*bar_width
            ax.set_xlim(l_xlim, r_xlim)
            # set special ylims
            'time' in col and ax.set_ylim(bottom=-1,top=5) # only for cpu_time/check_feas_time
            'solved' in col and ax.set_ylim(bottom=-1,top=3) # only for vf/ub solved
            
            # set other figure elements
            ax.set_title("Performance for "+ plotCols[col])
            ax.set_ylabel("Percentage of Base Case (gap=10)")
            ax.set_xlabel("Instances")
            ax.legend(loc='upper right', ncol=len(gaps), fontsize='x-small')
            ax.grid(axis='y',alpha=0.2)

            fig.tight_layout()
            fig.savefig("./performance/barchart/"+scn+col, dpi=fig.dpi)
            # Note: transparency is not supported in .eps format
            # fig.savefig("./performance/barchart/"+scn+col+'.eps', format='eps', dpi=600)

def perfProf(df, plotname=None, fixmin=None, xmin=1, xmax=None, legendnames={}):
    """
        Generate a performance profile plot for the given dataframe.
        Assume data given are in number types.
        x-axis label: multiple of virtual best;
        y-axis label: franction of instances.
        Input:
            df: instances as index, field-to-plot as columns
            plotname: name of the plot
            fixmin: the base value used to compute ratio; using df min if not given
            xmin: the smallest x-ticker to display; set by xlim
            xmax: the largest x-ticker to display; set by xlim
            displaynames: a dictionary contains legend name; using df col name if not given
    """

    colors = ['dodgerblue', 'slateblue', 'blueviolet', 'palevioletred','lightcoral',
        'sandybrown', 'gold', 'yellowgreen', 'darkturquoise']
    # 'chartreuse'
    # colors = ['red', 'lime', 'blue', 'fuchsia', 'aqua', 'yellow', 'black', 'gold']
    
    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}
        
    # find min value in the dataframe
    col_list = df.columns.values.tolist()
    df['min_val'] = df[col_list].min(axis=1)
    print(df['min_val'])

    # start ploting
    fig, ax = plt.subplots(1,1) # figsize = (6,5)
    
    for i, col in enumerate(col_list):
        # for each col, compute ratio
        ratios = df[col] / df['min_val']
        # i == 1 and print(ratios)
        uniq_ratios = ratios.unique()
        uniq_ratios.sort() # sort in place
        cum_cnt =  np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_prob = cum_cnt / len(ratios)
        
        # form x-tickers: if xmax is not given, use current max and round up
        if xmax == None:
            xmax = np.ceil(uniq_ratios[-1])
        elif uniq_ratios[-1] < xmax:
            np.append(uniq_ratios, xmax) # append array at the boundary point
            np.append(cum_prob, cum_prob[-1])
                
        # add turning points and form series to plot
        x_val = []
        y_val = []
        x_val.append(uniq_ratios[0])
        y_val.append(cum_prob[0])
        for j, r in enumerate(uniq_ratios[1:]):
            x_val.extend([r, r])
            y_val.extend([cum_prob[j], cum_prob[j+1]])

        if legendnames:
            ax.plot(x_val, y_val, label=legendnames[col], color=colors[i])
        else:
            ax.plot(x_val, y_val, label=col, color=colors[i])
        
    # set plot properties
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', direction='in', right=True)

    # set other figure elements
    ax.set_title("Performance profile: " + plotname[9:])
    ax.set_xlabel("Multiple of virtual best")
    ax.set_ylabel("Fraction of instances")
    ax.legend(loc='lower right', markerscale=1.25, frameon=False, labelspacing=0.35)

    fig.tight_layout()
    if plotname == None:
        fig.savefig("./performance/perfprof/"+"test", dpi=fig.dpi)
    else:
        fig.savefig("./performance/perfprof/"+plotname, dpi=fig.dpi)
        # fig.savefig("./performance/barchart/"+plotname+'.eps', format='eps', dpi=600)

def plotPerf(df, plotCols, plotScns, plotDSs, gaps=[]):
    """
        Prepare data for plotting performance profile; running time only.
        Input:
            df: pandas dataframe output from processTable
            plotCol: columns to make single plots
            plotScns: scenarios on one plot
            gaps: a list of given integer gaps; seen as s type of scenarios 
    """
    df = df[list(plotScns.keys())]
    # replace unsolved cases by a large number
    for scn in df.columns:
            df[scn] = pd.to_numeric(df[scn], errors='coerce').replace(np.nan, 1e+11)
    if gaps:
        for g in gaps:
            for col in plotCols:
                # prepare data and replace unsolved cases by a large number
                df_sub = df.xs((g, col), level=['gaps', 'fields'], axis=1, drop_level=True).copy()
                for scn in df_sub.columns:
                    df_sub[scn] = pd.to_numeric(df_sub[scn], errors='coerce').replace(np.nan, 36000)
                # call plot function
                perfProf(df_sub, plotname='perfprof_'+col+'_g'+str(g), xmax=400, legendnames=plotScns)
    else:
        for ds in plotDSs:
            # apply index filter on solution time
            df_time = df.xs((ds, 'cpu_time'), level=['datasets', 'fields'], axis=1, drop_level=True).copy()
            # df_time = pd.to_numeric(df_time, errors='coerce').replace(np.nan, 36000)
            # filter out cases where time is < 5'' or > 3600'' for all methods
            col_list = df_time.columns.values.tolist()
            drop_lessthan = df_time[(df_time[col_list] < 5).all(axis=1)].index.tolist()
            drop_morethan = df_time[(df_time[col_list] > 3600).all(axis=1)].index.tolist()
            drop_list = list(set(drop_lessthan) | set(drop_morethan))
            # print(drop_list)
            # df = df.drop(drop_list)
            for col in plotCols:
                df_sub = df.xs((ds, col), level=['datasets', 'fields'], axis=1, drop_level=True).copy()
                df_sub = df_sub.drop(drop_list)
                print(df_sub)
                # call plot function
                perfProf(df_sub, plotname='perfprof_'+col+'_'+ds, xmax=plotCols[col][1], legendnames=plotScns)

def main():

    # choose based on experiment results
    # gaps = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    # gaps = [10, 20]
    # gaps = sorted(gaps)
    gaps = []
    
    ######################### Run Experimests #########################
    # local: provide paths in runparams.py
    # runExperiments(exe, instanceDirs, outputDir, mibsParamsInputs, gaps)
    
    # using pbs file: provide paths in runparams.py
    # runExperimentsPBS(exe, instanceDirs, outputDir, mibsParamsInputs, pbsfile, gaps)
    
    ################# Process & Save | Load from CSV ###################
    # specify summary file name
    file_csv = "summary_cuts.csv"

    args = sys.argv[1:]  # change to arg parser later...
    if len(args) == 0: 
        df_r = parseOutput(outputDir, mibsParamsInputs, gaps, writeCSV=True, filename=file_csv)
    else:
        try:
            df_r = pd.read_csv(file_csv)
            set_cond = (df_r['scenario'].isin(list(mibsParamsInputs.keys()))) | (df_r['dataset'].isin(list(instanceDirs.keys())))
            df_r = df_r[set_cond]

        except FileNotFoundError:
            print('{} does not exist in current directory.'.format(file_csv))
        else:
            print('Reading from', file_csv)
    
    ################### Format Data & Print Table ####################
    # specify txt file name to print tables in LATEX
    file_txt = "ltx_tb_cut.txt"

    # columns to process and print
    displayCols = {
        'cpu_time': 'CPU Search Time',
        'nodes_proc': 'Number of Processed Nodes'
        # 'chk_feas_time': 'Check Feasibility Time',
        # 'vf_solved': 'Number of VF problem solved',
        # 'ub_solved': 'Number of UB problem solved',
        # 'objval_cost': 'Object Value'
    }

    df_proc = processTable(df_r, displayCols, writeLTX=False, filename=file_txt)
    
    ################### Make Bar Charts ####################
    '''
    # columns to compare in the plot
    plotCols = {
        'cpu_time': 'CPU Search Time',
        # 'chk_feas_time': 'Check Feasibility Time',
        # 'vf_solved': 'Number of VF problem solved',
        # 'ub_solved': 'Number of UB problem solved',
        # 'objval_cost': 'Object Value'
    } # this is the same as displayCols or a subset of it

    plotScns = mibsParamsInputs.keys()

    if set(plotCols).issubset(set(displayCols)):
        plotBarChart(df_proc, gaps, plotCols, plotScns)
    else:
        print('plotCols is not a subset of available columns. Try available ones.')
        plotCols_sub = list(set(displayCols)&set(plotCols))
        plotBarChart(df_proc, gaps, plotCols_sub, plotScns)
    '''
    ################### Make Performance Profile ####################
    # columns to compare in the plot
    plotCols = {
        'cpu_time': ['CPU Search Time', 25],
        'nodes_proc': ['Number of Processed Nodes', 75]
    }
    # plotGaps = [10, 20]

    # scenarios name dict used for legend; if read from .csv check name match
    plotScns = {}
    # manual input example:
    # for k in mibsParamsInputs.keys():
    #     if '01' in k:
    #         plotScns[k] = 'linkingBranchStrategy'
    #     else:
    #         plotScns[k] = 'fractionalBranchStrategy'    
    
    for k in mibsParamsInputs.keys():
        plotScns[k] = k

    # plotPerf(df_proc, plotCols, plotScns, plotGaps)
    plotPerf(df_proc, plotCols, plotScns, instanceDirs)


if __name__ == "__main__":
    main()