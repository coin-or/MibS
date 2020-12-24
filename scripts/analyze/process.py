# Script to parse the running result from MibS.
# See run2020.sh for output name and path.
# Last edited by yux616
# Dec 2020 
# Some os function requires Python 3.5+

import sys, os, collections
import pandas as pd
import numpy as np 

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
    
    # OPTION 2: for each displayCol, print a table
    with open('ltx_tb2.txt', 'a') as file:
        for col in displayCols:
            file.write(df_forprint.loc[:, (slice(None), col)].to_latex())

    # print(df_forprint.columns)
    return

def main():

    # choose based on experiments
    gaps = [25, 50]
    # gaps = [10, 20, 30, 40, 50, 60, 70, 80, 90]

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
    
    # the columns that we want to print
    displayCols = [
        'cpu_time',
        'chk_feas_time',
        'vf_solved',
        'ub_solved',
        'objval_cost'
    ]

    # print tables
    processTable(df_r, gaps, displayCols)

if __name__ == "__main__":
    main()