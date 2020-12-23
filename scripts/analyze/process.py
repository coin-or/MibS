# Script to parse the running result from MibS.
# See run2020.sh for output name and path.
# Last edited by yux616
# Dec 2020 
# Some os function requires Python 3.5+

import sys, os, collections
import pandas as pd 

def parseOutput(outputDir, writeCSV = True):
    """
        The function parse the output file in the given directory.
        Assume the result subfolders are named as BR$GAP$Output.
        The result will also be written to a .csv file if not specified.
        Input: 
            outputDir: string, a path to the parent output directory
            writeCSV: boolean, whether to save the results in structured format
        Return:
            a pd dataframe containing results
    """

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

    # idxflag = True
    results = collections.defaultdict(list) 

    gaps = [10, 25, 50]
    # gaps = [25]

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
                                    results['cpu_time'].append(( (line.split(':')[1]).split()[0] ))
                                
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
                                    results['prob_gap'].append(float(line.split(' ')[5].strip('%\n')))
                                
                                elif keywords[10] in line:
                                    results['ul_int_var'].append(int(line.split(':')[1]))
                                
                                elif keywords[11] in line:
                                    results['ll_int_var'].append(int(line.split(':')[1]))
                                else:
                                    pass

    df_result = pd.DataFrame(results)

    # write results to .csv file
    if writeCSV: 
        # df_result.to_csv('test.csv', mode='a')
        df_result.to_csv('summary.csv')
    
    return df_result

def printTable():
    """
        Print a summary table as page 34 in MibS paper
    """
    pass

def compareObj():
    """
        Print/plot some results: obj vs. SL_gap
    """
    pass

def main():
    args = sys.argv[1:]
    if len(args) == 0: 
        cwd = os.getcwd() #os.path.realpath('..')
        outputDir = os.path.join(cwd, '../../output/')
        df_r = parseOutput(outputDir, True)
    else:
        try:
            df_r = pd.read_csv("summary.csv")
        except FileNotFoundError:
            print('summary.csv not exist in current directory')
        else:
            print('read from summary')

if __name__ == "__main__":
    main()