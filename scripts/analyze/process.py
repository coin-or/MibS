# Script to parse the running result from MibS.
# See run2020.sh for output name and path.
# Last edited by yux616
# Dec 2020 
# Some os function requires Python 3.5+

import sys, os, collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import itertools

def parseOutput(outputDir, gaps, writeCSV = True):
    """
        The function parse the output file in the given directory.
        Assume the result subfolders are named as BR$GAP$Output.
        The result will also be written to a .csv file if not specified.
        Input: 
            outputDir: string, a path to the parent output directory
            writeCSV: boolean, whether to save the results in structured format
        Return:
            a pandas dataframe containing results
    """
    # may move those up as input in the future...
    # then need to match keywords and fields
    keywords = [
        'ALPS did not find a solution',
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
        'integer LL Variables'
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

    # iterate over gap size, dataset folders, and files in each folder to check for results
    for g in gaps:
        resultDir = os.path.join(outputDir, "BR"+str(g)+"Output")

        # iterate over different datasets
        with os.scandir(resultDir) as dataset_it:
            for f_entry in dataset_it:
                
                # iterate over files in the dolfer
                with os.scandir(f_entry.path) as output_it: 
                    for o_entry in output_it:
                        # start to write result to the dictionary
                        results['given_gap'].append(g)
                        results['dataset'].append(f_entry.name)
                        results['instance'].append(o_entry.name.split('.')[0])
                        
                        # read value for each field from file
                        with open(o_entry.path,'r') as file:
                            for line in file:
                                if keywords[0] in line: # !!! Need to append "" to pad  
                                    print("infeasible instance:", o_entry.name)
                                
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

    df_result = pd.DataFrame(results)
    
    # make some adjustment to formats
    # display check feasibility time as % of search time? 
    # sum vf+ub time -> feasibility time (or read from output directly?)
    df_result['chk_feas_time'] = df_result['ub_time'] + df_result['vf_time']
    df_result['chk_feas_time'] = df_result['chk_feas_time'].astype(float).round(2)
    df_result['cpu_time'] = df_result['cpu_time'].astype(float).round(2)
     
    # write results to .csv file
    if writeCSV: 
        # df_result.to_csv('test.csv', mode='a')
        df_result.to_csv('summary.csv')
    
    return df_result

def processTable(df, gaps, displayCols):
    """
        Print a summary table for required columns and gaps.
        Input:
            df_r: a dataframe with all info from parseOutput
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
    
    # collect required info into dict
    rsltDict = {}
    for inst in instList:
        rsltDict[inst] = {}
        for g in gaps:
            cond = (df['given_gap'] == g) & (df['instance'] == inst)
            df_temp = df[cond]
            if df_temp['solved'].values[0]:
                # using (outer key, inner key) for table column formatting
                rsltDict[inst].update({(g, col):df_temp[col].values[0] for col in displayCols})
            else:
                # special color for unsolved cases? but still need to process after print_to_latex
                rsltDict[inst].update({(g, col):'textcolor{lightgray}{'+str(df_temp[col].values[0])+'}' for col in displayCols})
        
        # for debug print
        # inst == 'miblp_20_20_50_0110_10_1' and print(rsltDict[inst])
    
    # convert dict to structured df
    #  change to formal column names?
    df_forprint = pd.DataFrame.from_dict(rsltDict, orient='index')
    df_forprint.columns.names = ["gaps","fields"]

    # OPTION 1: print results to a single table: suggest to use when display col number < 2 
    # with open('ltx_tb1.txt', 'w') as file:
    #     file.write(df_forprint.to_latex())
    
    # OPTION 2: for each displayCol, print a table; using slicer indexing
    with open('ltx_tb2.txt', 'a') as file:
        for col in displayCols:
            file.write(df_forprint.loc[:, (slice(None), col)].to_latex())

    # OPTION 3: just process table, do not print latex table to file
    # -- add a parameter for this

    return df_forprint

def plotSelected(df, gaps, plotCols):
    """
        Make some plots for given gaps and columns.
        Use the measures for the smallest gap as 1.
        Input:
            df: (processed!) pandas dataframe from processTable
            gaps: a list of int
            plotCol: columns to make single plots
    """
    # markers = itertools.cycle(('X', '+', '.', 'o', '*')) 
    # colors = itertools.cycle(('dodgerblue', 'slateblue', 'blueviolet', 'palevioletred','lightcoral',
        # 'sandybrown', 'gold', 'yellowgreen', 'darkturquoise')) #'orange',
    colors = ['dodgerblue', 'slateblue', 'blueviolet', 'palevioletred','lightcoral',
        'sandybrown', 'gold', 'yellowgreen', 'darkturquoise']

    # lambda function used to create new columns
    def computeRatio(r):
        # check unsolved cases
        if isinstance(r[g], str) or isinstance(r[gaps[0]], str): # and 'text' in r[g]: 
            return -1 
        # check rounded-to-zero cases
        elif r[gaps[0]] == 0 and r[g] == 0:
            return 1
        elif r[gaps[0]] == 0 and r[g] > 0:
            return r[g]/0.001
        elif r[g] == 0 and r[gaps[0]] > 0:
            return 0.001/r[gaps[0]]
        # otherwise, normal division
        else:
            return r[g]/r[gaps[0]]

    # plot some bar chart for summarized result
    for col in plotCols:
        # get columns to print using cross-section then flatten column index
        df_sub = df.xs(col, level='fields', axis=1, drop_level=True)
        df_sub.columns = df_sub.columns.to_flat_index()
 
        for g in gaps: # revesed order to drop cols? gaps[::-1]
            # handle the unsolved cases
            df_sub[str(g)+'_r'] = df_sub.apply(computeRatio, axis=1)

        df_plot = df_sub.loc[df_sub[str(gaps[0])+'_r'] >= 0, [str(g)+'_r' for g in gaps]]
        
        # begin plotting
        fig, ax = plt.subplots(1,1, figsize=(30,5)) #figsize=(7.5,5)
        # obtain xlocs and shifts
        xloc = np.arange(len(df_plot.index))
        mid = int(len(gaps) / 2)
        bar_width = 1/(len(gaps)+1)
        for i,g in enumerate(gaps):
            ax.bar(xloc+(i-mid)*(bar_width*0.6), df_plot[str(g)+'_r'], 
                label='TargetGap ='+str(g), color=colors[i], width=bar_width, alpha=0.5)
        # set other elements
        ax.set_title("Performance Plot for "+ plotCols[col])
        ax.set_ylabel("Percentage of (base?) Case")
        ax.set_xlabel("Instances")
        ax.legend(loc='upper right', ncol=len(gaps), fontsize='x-small')
        ax.grid(axis='y',alpha=0.2)
        ax.set_xticks(xloc)
        ax.set_xticklabels(df_plot.index, rotation=90, fontsize='x-small')
        l_xlim = xloc[0] - (mid+1)*bar_width
        r_xlim = xloc[-1] + (mid+1)*bar_width
        ax.set_xlim(l_xlim, r_xlim)
        # set special ylims
        'time' in col and ax.set_ylim(bottom=-1,top=5) # only for cpu_time/check_feas_time
        'solved' in col and ax.set_ylim(bottom=-1,top=3) # only for vf/ub solved

        fig.tight_layout()
        fig.savefig("./performance/test_"+col, dpi=fig.dpi)

def main():

    # choose based on experiment results
    gaps = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    # gaps = [25, 50]
    gaps = sorted(gaps)

    # read from output file directly or read from cache
    # ? may change to a more elegant way later...
    args = sys.argv[1:] 
    if len(args) == 0: 
        cwd = os.getcwd() #os.path.realpath('..')
        outputDir = os.path.join(cwd, '../../output/')
        df_r = parseOutput(outputDir, gaps, True)
    else:
        try:
            df_r = pd.read_csv("summary.csv")
        except FileNotFoundError:
            print('summary.csv not exist in current directory')
        else:
            print('read from summary file')
    
    # the columns that we want to process and print
    displayCols = {
        'cpu_time': 'CPU Search Time',
        'chk_feas_time': 'Check Feasibility Time',
        'vf_solved': 'Number of VF problem solved',
        'ub_solved': 'Number of UB problem solved',
        'objval_cost': 'Object Value'
    }

    # print tables
    df_proc = processTable(df_r, gaps, displayCols)

    # make some plots
    plotCols = {
        'cpu_time': 'CPU Search Time',
        'chk_feas_time': 'Check Feasibility Time',
        'vf_solved': 'Number of VF problem solved',
        'ub_solved': 'Number of UB problem solved',
        'objval_cost': 'Object Value'
    } # this is the same as the above list or a subset

    if set(plotCols).issubset(set(displayCols)):
        plotSelected(df_proc, gaps, plotCols)
    else:
        print('plotCols is not a subset of available columns. Try available ones.')
        plotCols_sub = list(set(displayCols)&set(plotCols))
        plotSelected(df_proc, gaps, plotCols_sub)

if __name__ == "__main__":
    main()