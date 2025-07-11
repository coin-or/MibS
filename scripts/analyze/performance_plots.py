# Script to run MibS with different SL target gap.
# The script may produce auxiliary folder/files in run directory.
# Last edited by yux616
# Apr 2020
# Script path:  /MibS/scripts/analyze
# Some os function requires Python 3.5+

# add arg parser later

import sys
import os
import collections
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def parseInstanceOutput(o_entry, d_entry, results, keywords, opt_values, etol):

    incomplete = True  # mark incomplete output file
    nosoln = False  # mark no soluntion found
    results['num_cuts'].append(0)
    results['cut_time'].append(0)
    results['cg_called'].append(0)
    results['time_per_cg_call'].append(0)
    results['root_bound'].append(10000000)
    results['100_bound'].append(10000000)
    results['num_idic'].append(0)
    results['num_idic_fail'].append(0)
    results['idic_fail_rate'].append("")
    results['num_full_int_idic'].append(0)
    results['num_link_int_idic'].append(0)
    results['num_lower_int_idic'].append(0)
    results['num_frac_idic'].append(0)
    results['num_full_int_idic_fail'].append(0)
    results['num_link_int_idic_fail'].append(0)
    results['num_lower_int_idic_fail'].append(0)
    results['num_frac_idic_fail'].append(0)
    results['full_int_idic_fail_rate'].append("")
    results['link_int_idic_fail_rate'].append("")
    results['lower_int_idic_fail_rate'].append("")
    results['frac_idic_fail_rate'].append("")
    results['num_isic'].append(0)
    results['num_isic_fail'].append(0)
    results['isic_fail_rate'].append("")
    results['num_full_int_isic'].append(0)
    results['num_link_int_isic'].append(0)
    results['num_lower_int_isic'].append(0)
    results['num_frac_isic'].append(0)
    results['num_full_int_isic_fail'].append(0)
    results['num_link_int_isic_fail'].append(0)
    results['num_lower_int_isic_fail'].append(0)
    results['num_frac_isic_fail'].append(0)
    results['full_int_isic_fail_rate'].append("")
    results['link_int_isic_fail_rate'].append("")
    results['lower_int_isic_fail_rate'].append("")
    results['frac_isic_fail_rate'].append("")
    results['num_hypercube'].append(0)
    results['num_gen_no_good'].append(0)
    results['num_ben_binary'].append(0)
    results["vf_solved"].append("")
    results["ub_solved"].append("")
    #if ('1.0.0-opt' not in versions and
    #    'filmosi' not in versions):
    if 1:
        results["vf_time"].append(0)
        results["ub_time"].append(0)
    results["objval"].append(-1000000)
    results["gap"].append(1000000)

    # read value for each field from file
    with open(o_entry.path, "r") as file:
        for line in file.read().splitlines():
            if (len(line.split()) > 1 and
                line.split()[0] == '0' and
                results["version"][-1] != 'filmosi'):
                if len(line.split()) == 4:
                    results["root_bound"][-1] = float(line.split()[1])
                else:
                    results["root_bound"][-1] = float(line.split()[2])

            if (len(line.split()) > 0 and line.split()[0] == '100'):
                results["100_bound"][-1] = float(line.split()[2])

            if keywords["solved"] in line:
                nosoln = True
                results["solved"].append(False)
                #print("No solution found instance:", o_entry.name, s)

            elif keywords["nodes"] in line or \
                 keywords["nodes_full_proc"] in line:
                results["nodes"].append(int(line.split(":")[1]))

            elif keywords["cpu"] in line:
                results["cpu"].append( float((line.split(":")[1]).split()[0]))

            elif keywords["vf_solved"] in line:
                results["vf_solved"][-1] = int(line.split("=")[1])

            elif keywords["ub_solved"] in line:
                results["ub_solved"][-1] = int(line.split("=")[1])

            elif keywords["vf_time"] in line: #and
                  #'1.0.0-opt' not in versions and
                  #'filmosi' not in versions):
                results["vf_time"][-1] = float(line.split("=")[1])

            elif keywords["ub_time"] in line: #and
                  #'1.0.0-opt' not in versions and
                  #'filmosi' not in versions):
                results["ub_time"][-1] = float(line.split("=")[1])
                
            elif keywords["objval"] in line:
                results["objval"][-1] = int(float(line.split()[6]))
                
            elif keywords["gap"] in line:
                incomplete = False
                if "infinity" in line:
                    results["gap"][-1] = 1000000
                    # no soln found
                    # results['solved'].append(False)
                else:
                    if nosoln:
                        print("Gap is incorrectly reported")
                    solgap = float(line.split(" ")[5].strip("%\n"))
                    results["gap"][-1] = solgap
                    # mark unsolved instances in given time limit
                    if nosoln == False:
                        if solgap - 0.0 < etol:
                            results["solved"].append(True)
                        else:
                            results["solved"].append(False)
            elif keywords["num_full_int_idic_fail"] in line:
                results["num_full_int_idic_fail"][-1] = int(line.split(':')[1])
                results["num_idic_fail"][-1] += results["num_full_int_idic_fail"][-1]
                n = (results["num_full_int_idic"][-1]+
                     results["num_full_int_idic_fail"][-1])
                if n > 0:
                    results["full_int_idic_fail_rate"][-1] = \
                        round(float(results["num_full_int_idic_fail"][-1])/n, 2)
            elif keywords["num_link_int_idic_fail"] in line:
                results["num_link_int_idic_fail"][-1] = int(line.split(':')[1])
                results["num_idic_fail"][-1] += results["num_link_int_idic_fail"][-1]
                n = (results["num_link_int_idic"][-1] +
                     results["num_link_int_idic_fail"][-1])
                if n > 0:
                    results["link_int_idic_fail_rate"][-1] = \
                        round(float(results["num_link_int_idic_fail"][-1])/n, 2)
            elif keywords["num_lower_int_idic_fail"] in line:
                results["num_lower_int_idic_fail"][-1] = int(line.split(':')[1])
                results["num_idic_fail"][-1] += results["num_lower_int_idic_fail"][-1]
                n = (results["num_lower_int_idic"][-1]+
                     results["num_lower_int_idic_fail"][-1])
                if n > 0:
                    results["lower_int_idic_fail_rate"][-1] = \
                        round(float(results["num_lower_int_idic_fail"][-1])/n,  2)
            elif keywords["num_frac_idic_fail"] in line:
                results["num_frac_idic_fail"][-1] = int(line.split(':')[1])
                results["num_idic_fail"][-1] += results["num_frac_idic_fail"][-1]
                n = (results["num_frac_idic"][-1] +
                     results["num_frac_idic_fail"][-1])
                if n > 0:
                    results["frac_idic_fail_rate"][-1] = \
                        round(float(results["num_frac_idic_fail"][-1])/n, 2)
            elif keywords["num_full_int_idic"] in line:
                results["num_full_int_idic"][-1] = int(line.split(':')[1])
                results["num_idic"][-1] += results["num_full_int_idic"][-1]
            elif keywords["num_link_int_idic"] in line:
                results["num_link_int_idic"][-1] = int(line.split(':')[1])
                results["num_idic"][-1] += results["num_link_int_idic"][-1]
            elif keywords["num_lower_int_idic"] in line:
                results["num_lower_int_idic"][-1] = int(line.split(':')[1])
                results["num_idic"][-1] += results["num_lower_int_idic"][-1]
            elif keywords["num_frac_idic"] in line:
                results["num_frac_idic"][-1] = int(line.split(':')[1])
                results["num_idic"][-1] += results["num_frac_idic"][-1]
            elif keywords["num_full_int_isic_fail"] in line:
                results["num_full_int_isic_fail"][-1] = int(line.split(':')[1])
                results["num_isic_fail"][-1] += results["num_full_int_isic_fail"][-1]
                n = (results["num_full_int_isic"][-1]+
                     results["num_full_int_isic_fail"][-1])
                if n > 0:
                    results["full_int_isic_fail_rate"][-1] = \
                        round(float(results["num_full_int_isic_fail"][-1])/n, 2)
            elif keywords["num_link_int_isic_fail"] in line:
                results["num_link_int_isic_fail"][-1] = int(line.split(':')[1])
                results["num_isic_fail"][-1] += results["num_link_int_isic_fail"][-1]
                n = (results["num_link_int_isic"][-1] +
                     results["num_link_int_isic_fail"][-1])
                if n > 0:
                    results["link_int_isic_fail_rate"][-1] = \
                        round(float(results["num_link_int_isic_fail"][-1])/n, 2)
            elif keywords["num_lower_int_isic_fail"] in line:
                results["num_lower_int_isic_fail"][-1] = int(line.split(':')[1])
                results["num_isic_fail"][-1] += results["num_lower_int_isic_fail"][-1]
                n = (results["num_lower_int_isic"][-1]+
                     results["num_lower_int_isic_fail"][-1])
                if n > 0:
                    results["lower_int_isic_fail_rate"][-1] = \
                        round(float(results["num_lower_int_isic_fail"][-1])/n, 2)
            elif keywords["num_frac_isic_fail"] in line:
                results["num_frac_isic_fail"][-1] = int(line.split(':')[1])
                results["num_isic_fail"][-1] += results["num_frac_isic_fail"][-1]
                n = (results["num_frac_isic"][-1] +
                     results["num_frac_isic_fail"][-1])
                if n > 0:
                    results["frac_isic_fail_rate"][-1] = \
                        round(float(results["num_frac_isic_fail"][-1])/n, 2)
            elif keywords["num_full_int_isic"] in line:
                results["num_full_int_isic"][-1] = int(line.split(':')[1])
                results["num_isic"][-1] += results["num_full_int_isic"][-1]
            elif keywords["num_link_int_isic"] in line:
                results["num_link_int_isic"][-1] = int(line.split(':')[1])
                results["num_isic"][-1] += results["num_link_int_isic"][-1]
            elif keywords["num_lower_int_isic"] in line:
                results["num_lower_int_isic"][-1] = int(line.split(':')[1])
                results["num_isic"][-1] += results["num_lower_int_isic"][-1]
            elif (keywords["num_frac_isic"] in line and
                  'Fractional solutions will be separated' not in line and
                  'Fractional cuts will be generated.' not in line and
                  'useFractionalCuts' not in line and
                  'MibSBranchingStrategyFractional' not in line):
                results["num_frac_isic"][-1] = int(line.split(':')[1])
                results["num_isic"][-1] += results["num_frac_isic"][-1]
            elif keywords["num_hypercube"] in line:
                results["num_hypercube"][-1] = int(line.split(':')[1])
            elif keywords["num_gen_no_good"] in line:
                results["num_gen_no_good"][-1] = int(line.split(':')[1])
            elif keywords["num_ben_binary"] in line:
                results["num_ben_binary"][-1] = int(line.split(':')[1])
            elif keywords['ul_int_var'] in line:
               results['ul_int_var'].append(int(line.split(':')[1]))
            elif keywords['ll_int_var'] in line:
               results['ll_int_var'].append(int(line.split(':')[1]))
            elif keywords["num_cuts"] in line:
                found = True
                results["num_cuts"][-1] = int(line.split(" ")[8])
                results["cut_time"][-1] = float(line.split(" ")[12])
                results["cg_called"][-1] = float(line.split(" ")[5])
                if results["num_cuts"][-1] > 0:
                    results["time_per_cg_call"][-1] = results["cut_time"][-1]/results["num_cuts"][-1]
            elif keywords["infeasible"] in line and results["version"][-1] != 'filmosi':
                print("Infeasible instance!")
            elif 'STAT;' in line and len(line.split(';')) > 6:
                incomplete=False
                results["objval"][-1] = float(line.split(";")[2])
                results["root_bound"][-1] = float(line.split(";")[4])
                results["cpu"].append(float(line.split(";")[5]))
                results["nodes"].append(float(line.split(";")[8]))
                results["gap"][-1] = float(line.split(";")[10])
                results["vf_solved"][-1] = 0
                results["ub_solved"][-1] = 0
                #results["vf_time"][-1] = 0
                #results["ub_time"][-1] = 0
                if results["gap"][-1] - 0.0 < etol:
                    results["solved"].append(True)
                else:
                    results["solved"].append(False)
            elif 'User cuts applied:' in line:
                results["num_cuts"][-1] = int(line.split()[3])
            else:
                pass

    if results["num_idic"][-1]+results["num_idic_fail"][-1] > 0:
        results["idic_fail_rate"][-1] = (
            results["num_idic_fail"][-1]/(
                results["num_idic"][-1]+results["num_idic_fail"][-1]))
    if results["num_isic"][-1]+results["num_isic_fail"][-1] > 0:
        results["isic_fail_rate"][-1] = (
            results["num_isic_fail"][-1]/(
                results["num_isic"][-1]+results["num_isic_fail"][-1]))
    if incomplete:
        results["nodes"].append(1000000000)
        results["cpu"].append(3600)
        results["solved"].append(False)
        #print("Incomplete instance:", o_entry.name, s)

    if results["solved"][-1] == False:
        results["cpu"][-1] = 3600
    else:
        if (opt_values[results["instance"][-1]] == np.inf):
            opt_values[results["instance"][-1]] = results["objval"][-1]
        elif opt_values[results["instance"][-1]] != results["objval"][-1]:
            print("************ Warning: objective values don't agree!")
            print("************ ", results["instance"][-1],
                  opt_values[results["instance"][-1]], ' != ',
                  results["objval"][-1])
                                        
    #if results["cpu"][-1] < 0.01:
    #   print ("Small value! ", results["cpu"][-1], results["instance"][-1])
    #   results["cpu"][-1] = .01
    instance = results["instance"][-1]
    if (opt_values[instance] not in [0, np.inf] and
        results['root_bound'][-1] != 10000000):
        results["root_gap"].append(round(
            100*abs((opt_values[instance] -
                     results["root_bound"][-1])/opt_values[instance]),2))
    else:
        results["root_gap"].append(100000)
    #print(instance, results["scenario"][-1],
    #      results["root_bound"][i], opt_values[instance])
    if (opt_values[instance] not in [0, np.inf] and
        results['100_bound'][-1] != 10000000):
        results["100_gap"].append(round(
            100*abs((opt_values[instance] -
                     results["100_bound"][-1])/opt_values[instance]),2))
    else:
        results["100_gap"].append(100000)
    #print(instance, results["scenario"][i],
    #      results["100_bound"][-1], opt_values[instance])

    return

def parseOutput(outputDir, versions, scenarios, keywords, dataSets,
                name='', debug=False):
    """
    This function parse the output files in the given directory.
    It assumes the following folder hierarch: 
       outputDir
          version1
             param_scenario1
                instance_set1 
                   instance1.out
                   instance2.out
                   ...
                instance_set2 
                   instance1.out
                   instance2.out
                   ...
                ...
             param_scenario2
                instance_set1 
                   instance1.out
                   instance2.out
                   ...
                instance_set2 
                   instance1.out
                   instance2.out
                   ...
                ...
             ...
          version2
             ...         
    The results will be written to a .csv file if not specified.
    Input:
        outputDir: string, a path to the parent output directory
        versions:  dictionary of versions, where the keys are the directory 
                   names and the value is a name string to be used in the figures
        scenarios: dictionary of scenarios, where the keys are the directory 
                   names and the value is a name string to be used in the figures
        keywords:  dictionary of keywords, where the keys are the statistics
                   to be stored and the values are unique strings found on lines
                   of the file containing those data elements.
        writeCSV:  boolean, whether to save the results in structured format
        filename:  string, name of .csv file
        name:      a name string to be used to identify the figures
    Return:
        a pandas dataframe containing results
    """

    results = collections.defaultdict(list)
    opt_values = {}
    etol = np.finfo(float).eps

    # iterate over versions, scenarios, datasets, and files in each folder
    for v in versions:
        for s in scenarios:
            found = False
            for dir in outputDir:
                resultDir = os.path.join(dir, v, s)
                if os.path.isdir(resultDir) == True:
                    found = True
                    break
            if not found:
                continue
            print("Compiling results for version", v+",", "scenario", s)
            # iterate over different datasets available
            with os.scandir(resultDir) as dataset_it:
                for d_entry in dataset_it:
                    if d_entry.name not in dataSets:
                        continue
                    # iterate over files in the folder
                    with os.scandir(d_entry.path) as output_it:
                        print("   Parsing files in dataset ", d_entry.name)
                        for o_entry in output_it:
                            if o_entry.name.endswith(".out"):
                                if name == '':
                                    results["dataset"].append(d_entry.name)
                                else:
                                    results["dataset"].append(name)
                                results["scenario"].append(scenarios[s])
                                results["version"].append(v)
                                results["instance"].append(
                                    os.path.splitext(o_entry.name)[0]
                                )
                                if results["instance"][-1] not in opt_values:
                                    opt_values[results["instance"][-1]] = np.inf
        
                                parseInstanceOutput(o_entry, d_entry, results,
                                                    keywords, opt_values, etol)
                                a = [len(results[k]) for k in results]
                                b = [k for k in results]
                                for k in range(len(a)):
                                    if a[k] != a[0]:
                                        print(k)
                                        print(results[b[k]])
                                        print(o_entry.path+"!!!!!!!!!!!!!", b[k])
                                    
    if debug:
        for k in results:
            print (k)
            print(len(results[k]))

    df = pd.DataFrame(results)

    # sum vf+ub time -> feasibility time
    if '1.0.0-opt' not in versions and 'filmosi' not in versions:
        df["chk_feas_time"] = df["ub_time"] + df["vf_time"]
        df["chk_feas_time"] = df["chk_feas_time"].astype(float).round(2)
    #df_result["cpu"] = df_result["cpu"].astype(float).round(2)

    return df

def export(df, columns=None, filename="summary.csv"):

    df_csv = df.copy()

    # Select columns with numeric values and empty strings
    def is_numeric_and_empty_string(column):
        return column.apply(lambda x: isinstance(x, (int, float)) or x == '').all()

    selected_columns = [col for col in df_csv.columns if is_numeric_and_empty_string(df_csv[col])]
    for s in selected_columns:
        df_csv[s] = pd.to_numeric(df_csv[s])

    if 0:
        average_values = df_csv.select_dtypes(include = ['number']).mean(skipna=True)
        print()
        df_csv.loc['Average'] = average_values
    else:
        # Means by (scenario, dataset)
        group_means = df_csv.groupby(['scenario','dataset'])[selected_columns].mean()
        group_means = group_means.reset_index()
        #group_means['dataset'] = group_means['dataset'].astype(str) + '_mean'
        #df_csv = pd.concat([group_means, df_csv], ignore_index=True)

        # Means by scenario
        group_means_s = df_csv.groupby('dataset')[selected_columns].mean()
        group_means_s = group_means_s.reset_index()
        #group_means['dataset'] = group_means['dataset'].astype(str) + '_mean'
        group_means = pd.concat([group_means, group_means_s], ignore_index=True)

        # Means by dataset
        group_means_d = df_csv.groupby('scenario')[selected_columns].mean()
        group_means_d = group_means_d.reset_index()
        #group_means['dataset'] = group_means['dataset'].astype(str) + '_mean'
        group_means = pd.concat([group_means, group_means_d], ignore_index=True)

    # write results to .csv file
    # df.to_csv(filename, mode='a', header=False, index=False) # append results only
    df_csv.to_csv(filename+".csv", columns=columns, index=False)
    group_means.to_csv(filename+"_means.csv", columns=columns[2:], index=False)

    return

def processTable(df, displayCols, writeLTX=False, filename="ltx_tb.txt"):
    """
    Print a summary table for required columns.
    Input:
        df: a dataframe with all info from parseOutput
        displayCol: columns to print
    """

    # separate instance to different tables
    # convert each instance related data into a dictionary
    # each data field can print to a table
    # or print a summary table where instance by row

    # obtain the list of instances
    instList = list(df.instance.unique())
    scnList = list(df.scenario.unique())
    versionList = list(df.version.unique())
    # print(instList)

    # collect required info into dict
    rsltDict = {}
    for inst in instList:
        rsltDict[inst] = {}
        if 'nw04' in inst:
            continue
        for scn in scnList:
            for v in versionList:
                cond = (
                    (df["scenario"] == scn)
                    & (df["instance"] == inst)
                    & (df["version"] == v)
                )
                df_temp = df[cond]
                if len(df_temp["dataset"].values) > 0:
                    ds = df_temp["dataset"].values[0]
                    rsltDict[inst].update(
                        {(scn, v, ds, col): df_temp[col].values[0] for col in displayCols}
                )

    # convert dict to structured df: change to formal column names?
    df_forprint = pd.DataFrame.from_dict(rsltDict, orient="index")
    #print (df_forprint)
    df_forprint.columns.names = ["scn", "v", "datasets", "fields"]
    df_forprint = df_forprint.sort_index()

    # OPTION 1: print results to a single table: suggest to use when display col number < 2
    # with open('ltx_tb1.txt', 'w') as file:
    #     file.write(df_forprint.to_latex())

    # OPTION 2: for each displayCol, print a table; using slicer indexing
    if writeLTX:
        with open(filename, "w") as file:
            for col in displayCols:
                for scn in scnList:
                    file.write(df_forprint.loc[:, (scn, slice(None), col)].to_latex())

    # OPTION 3: just process table, do not print latex table to file
    # pass

    return df_forprint

def dropFilter(df, scenarios, ds):
    """
    Prepare data for plotting performance profile; running time only.
    Input:
        df: pandas dataframe output from processTable
        plotCol: columns to make single plots
        scenarios: scenarios on one plot
    """
    df = df[scenarios.values()]
    # replace unsolved cases by a large number
    for scn in df.columns:
        df[scn] = pd.to_numeric(df[scn], errors="coerce").replace(np.nan, 1e11)
    # apply index filter on solution time
    df_time = df.xs(
        (ds, "cpu"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    df_solved = df.xs(
        (ds, "solved"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    df_gap = df.xs(
        (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    # df_time = pd.to_numeric(df_time, errors='coerce').replace(np.nan, 36000)
    # filter out cases where time is < 5'' or > 3600'' for all methods
    col_list = df_time.columns.values.tolist()
    
    drop_all_unsolved = df_solved[(df_solved[col_list] != True).all(axis=1)].index.tolist()
    drop_any_unsolved = df_solved[(df_solved[col_list] != True).any(axis=1)].index.tolist()
    drop_no_gap = df_gap[(df_gap[col_list] >= 1000000).all(axis=1)].index.tolist()
    if 1: #ds == 'INTERD-DEN':
        drop_easy = df_time[(df_time[col_list] < 1).all(axis=1)].index.tolist()
        drop_small_time = df_time[(df_time[col_list] <= 0.01).any(axis=1)].index.tolist()
        drop_list_time = list(set(drop_easy) | set(drop_all_unsolved) | set(drop_small_time))
    else:
        drop_easy = df_time[(df_time[col_list] < 5).all(axis=1)].index.tolist()
        drop_list_time = list(set(drop_easy) | set(drop_all_unsolved))
    ##drop_list_time.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    #drop_list_gap.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])

    df_solved = df.drop(drop_list_time)
    df_all_solved = df.drop(list(set(drop_easy) | set(drop_any_unsolved) |
                                 set(drop_small_time)))
    df_has_soln = df.drop(list(set(drop_easy) | set(drop_small_time) |
                               set(drop_no_gap)))

    # with pd.option_context('display.max_rows', None,
    #                        'display.max_columns', None,
    #                        'display.precision', 3,
    #                        'display.float_format', lambda x: '%.5f' % x,
    #                        ):
    #     print(df_solved)
    #     print(col_list)
    #     print(drop_all_unsolved)
    #     print(drop_any_unsolved)

    return df_all_solved, df_solved, df_has_soln

def plotPerfProf(
        df, versions, plotname="perf_profile", plottitle="Performance Profile",
        xmin=0.0, xmax=None, legendnames={}, versionlegend=False
):
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

    fig, ax = plt.subplots(1, 1)

    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}

    # find min value in the dataframe
    col_list = df.columns.values.tolist()
    df["virtual_best"] = df[col_list].min(axis=1)

    #print(df["virtual_best"])
    
    for col in col_list:
        #print(col)
        # for each col, compute ratio
        ratios = df[col] / df["virtual_best"]
        uniq_ratios = ratios.unique()
        uniq_ratios.sort()  # sort in place
        # with pd.option_context('display.max_rows', None,
        #                       'display.max_columns', None,
        #                       'display.precision', 3,
        #                       'display.float_format', lambda x: '%.5f' % x,
        # ):
        #    print(df[col])
        #    print(ratios.sort_values())
        #    print(ratios)
        #    print(uniq_ratios)
        cum_cnt = np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_frac = cum_cnt / len(ratios)

        # form x-tickers: if xmax is not given, use current max and round up
        if xmax == None:
            xmax = np.ceil(uniq_ratios[-1])
        elif uniq_ratios[-1] < xmax:
            np.append(uniq_ratios, xmax)  # append array at the boundary point
            np.append(cum_frac, cum_frac[-1])

        # add turning points and form series to plot
        x_val = []
        y_val = []
        x_val.append(1.0)
        y_val.append(0.0)
        if uniq_ratios[0] > 1:
            x_val.append(uniq_ratios[0])
            y_val.append(0)
        x_val.append(uniq_ratios[0])
        y_val.append(cum_frac[0])
        for j, r in enumerate(uniq_ratios[1:]):
            x_val.extend([r, r])
            y_val.extend([cum_frac[j], cum_frac[j + 1]])
            #print(r, cum_frac[j])
            #print(r, cum_frac[j+1])
        if cum_frac[-1] == 1.0:
            x_val.append(xmax)
            y_val.append(1.0)

        if legendnames:
            # , color=colors[i])
            plt.plot(x_val, y_val, label=legendnames[col])
        elif versionlegend:
            plt.plot(x_val, y_val, label=col[0]+':'+versions[col[1]])  # , color=colors[i])
        else:
            plt.plot(x_val, y_val, label=col[0])  # , color=colors[i])

    # set plot properties
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.02, 1.05)
    ax.tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax.set_title(plottitle)
    ax.set_xlabel("Multiple of virtual best")
    ax.set_ylabel("Fraction of instances")
    ax.legend(
        loc="lower right",
        #bbox_to_anchor=(0.9, 0.9),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)


def plotCumProf(df_time, df_gap, versions, plotname="cum_profile",
                plottitle = "Cumulative Profile",
                legendnames={}, versionlegend=False):

    fig = plt.figure()
    gs = fig.add_gridspec(1, 2, wspace=0)
    ax = gs.subplots(sharey=True)

    col_list = df_time.columns.values.tolist()
    time_buckets = range(0, 3600)

    # with pd.option_context('display.max_rows', None,
    #                        'display.max_columns', None,
    #                        'display.precision', 3,
    #                        'display.float_format', lambda x: '%.5f' % x,
    #                        ):
    #     print(df_time)
    #     print(df_gap)

    for col in col_list:
        #print(col)
        times = df_time[col]
        #print(times)
        cum_cnt = np.sum(np.array([times <= t for t in time_buckets]), axis=1)
        cum_frac = cum_cnt / len(df_time)
        #print(cum_frac)
        x_val = []
        if legendnames:
            ax[0].plot(time_buckets, cum_frac, label=legendnames[col])
        elif versionlegend:
            ax[0].plot(time_buckets, cum_frac, label=col[0]+':'+versions[col[1]])  # , color=colors[i])
        else:
            ax[0].plot(time_buckets, cum_frac, label=col[0])

    ax[0].set_xlim(0, 3599)
    ax[0].set_ylim(0.0, 1)
    ax[0].tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Fraction of instances")

    gap_buckets = np.linspace(0, 100, 1000)

    for col in col_list:
        #print(col)
        gaps = df_gap[col]
        #print(gaps)
        cum_cnt = np.sum(np.array([gaps <= g for g in gap_buckets]), axis=1)
        cum_frac = cum_cnt / len(df_gap)
        #print(cum_frac)
        x_val = []
        if legendnames:
            ax[1].plot(gap_buckets, cum_frac, label=legendnames[col])
        elif versionlegend:
            ax[1].plot(gap_buckets, cum_frac, label=col[0]+':'+versions[col[1]])  # , color=colors[i])
        else:
            ax[1].plot(gap_buckets, cum_frac, label=col[0])

    ax[1].set_xlim(0.0, 100)
    ax[1].tick_params(axis="both", direction="in", right=True)
    ax[1].legend(
        loc="lower right",
        #bbox_to_anchor=(0.9, 0.95),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    # set other figure elements
    ax[1].set_xlabel("Gap")
    ax[1].label_outer()

    fig.suptitle(plottitle)
    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)
    # fig.savefig("./performance/barchart/"+plotname+'.eps', format='eps', dpi=600)

def plotBaselineProf(
        df, versions, baseline, plotname="base_profile",
        plottitle="Baseline Profile",
        xmin=0.0, xmax=None, legendnames={}, versionlegend=False
):
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

    fig = plt.figure()
    gs = fig.add_gridspec(1, 2, wspace=0)
    ax = gs.subplots(sharey=True)

    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}

    # find min value in the dataframe
    col_list = df.columns.values.tolist()

    for col in col_list:
        if col == baseline or col[0] == "virtual_best":
            continue
        #print(col)
        # for each col, compute ratio
        ratios = df[col] / df[baseline]
        uniq_ratios = ratios.unique()
        uniq_ratios.sort()  # sort in place
        # with pd.option_context('display.max_rows', None,
        #                        'display.max_columns', None,
        #                        'display.precision', 2,
        #                        'display.float_format', lambda x: '%.5f' % x,
        # ):
        # #     print(df[col])
        #      print(ratios)
        # #     print(uniq_ratios)
        cum_cnt = np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_frac = cum_cnt / len(ratios)

        # form x-tickers: if xmax is not given, use current max and round up
        if xmax == None:
            xmax = np.ceil(uniq_ratios[-1])
        elif uniq_ratios[-1] < xmax:
            uniq_ratios = np.append(uniq_ratios, xmax)  # append array at the boundary point
            cum_frac = np.append(cum_frac, cum_frac[-1])

        #print(uniq_ratios)
        #print(cum_frac)

        # Values less than one are scaled differently
        if uniq_ratios[0] < 1:
            x_val = []
            y_val = []
            x_val.append(0.0)
            y_val.append(0.0)
            x_val.append(uniq_ratios[0])
            y_val.append(0)
            x_val.append(uniq_ratios[0])
            y_val.append(cum_frac[0])
            for j, r in enumerate(uniq_ratios[1:]):
                if r > 1:
                    x_val.append(r)
                    y_val.append(cum_frac[j])
                    break
                x_val.extend([r, r])
                # j is indexed starting at zero, not one!
                y_val.extend([cum_frac[j], cum_frac[j + 1]])

            if legendnames:
                # , color=colors[i])
                ax[0].plot(x_val, y_val, label=legendnames[col])
            elif versionlegend:
                ax[0].plot(x_val, y_val, label=col[0]+':'+versions[col[1]])  # , color=colors[i])
            else:
                ax[0].plot(x_val, y_val, label=col[0])  # , color=colors[i])

        # add turning points and form series to plot
        x_val = []
        y_val = []
        if uniq_ratios[0] >= 1:
            x_val.append(1.0)
            y_val.append(0.0)
            j = 0
        if uniq_ratios[0] > 1:
            x_val.append(uniq_ratios[0])
            y_val.append(0)
        x_val.append(uniq_ratios[j])
        y_val.append(cum_frac[j])
        
        for k, r in enumerate(uniq_ratios[j+1:]):
            x_val.extend([r, r])
            y_val.extend([cum_frac[k+j], cum_frac[k+j+1]])

        if legendnames:
            # , color=colors[i])
            ax[1].plot(x_val, y_val, label=legendnames[col])
        elif versionlegend:
            ax[1].plot(x_val, y_val, label=col[0]+':'+versions[col[1]])  # , color=colors[i])
        else:
            ax[1].plot(x_val, y_val, label=col[0])  # , color=colors[i])

    # set plot properties
    ax[0].set_xlim(0, 1)
    ax[0].set_ylim(-0.02, 1.05)
    ax[0].tick_params(axis="both", direction="in", right=True)

    ax[1].set_xlim(1, xmax)
    ax[1].label_outer()
    ax[1].tick_params(axis="both", direction="in", right=True)
    ax[1].legend(
        loc="lower right",
        #bbox_to_anchor=(0.9, 0.05),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    if versionlegend:
        fig.supxlabel("Ratio of baseline: v"+baseline[1]+", "+baseline[0])
    else:
        fig.supxlabel("Ratio of baseline: "+baseline[0])
    fig.supylabel("Fraction of instances")
    fig.suptitle(plottitle)
    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)

def plotBaselineProfSingle(
        df, versions, baseline, plotname="base_profile",
        plottitle="Baseline Profile",
        xmin=0.0, xmax=None, legendnames={}, versionlegend=False
):
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

    fig, ax = plt.subplots(1, 1)

    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}

    # find min value in the dataframe
    col_list = df.columns.values.tolist()

    for col in col_list:
        if col == baseline or col[0] == "virtual_best":
            continue
        #print(col)
        # for each col, compute ratio
        ratios = df[col] / df[baseline]
        #print(df[col])
        uniq_ratios = ratios.unique()
        uniq_ratios.sort()  # sort in place
        #print(uniq_ratios)
        cum_cnt = np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_frac = cum_cnt / len(ratios)

        # form x-tickers: if xmax is not given, use current max and round up
        if xmax == None:
            xmax = np.ceil(uniq_ratios[-1])
        elif uniq_ratios[-1] < xmax:
            uniq_ratios = np.append(uniq_ratios, xmax)  # append array at the boundary point
            cum_frac = np.append(cum_frac, cum_frac[-1])

        #print(uniq_ratios)
        #print(cum_frac)

        x_val = []
        y_val = []
        x_val.append(0.0)
        y_val.append(0.0)
        x_val.append(uniq_ratios[0])
        y_val.append(0)
        x_val.append(uniq_ratios[0])
        y_val.append(cum_frac[0])
        for j, r in enumerate(uniq_ratios[1:]):
            if r > 1:
                x_val.append(r)
                y_val.append(cum_frac[j])
                break
            x_val.extend([r, r])
            # j is indexed starting at zero, not one!
            y_val.extend([cum_frac[j], cum_frac[j + 1]])

        if legendnames:
            # , color=colors[i])
            plt.plot(x_val, y_val, label=legendnames[col])
        elif versionlegend:
            plt.plot(x_val, y_val, label=col[0]+':'+versions[col[1]])  # , color=colors[i])
        else:
            plt.plot(x_val, y_val, label=col[0])  # , color=colors[i])

    # set plot properties
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.05)
    ax.tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax.set_title(plottitle)
    if versionlegend:
        ax.set_xlabel("Ratio of baseline: v"+baseline[1]+", "+baseline[0])
    else:
        ax.set_xlabel("Ratio of baseline: "+baseline[0])
    ax.set_ylabel("Fraction of instances")
    ax.legend(
        loc="upper left",
        #bbox_to_anchor=(0.9, 0.05),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    fig.tight_layout()
    fig.savefig(plotname, dpi=fig.dpi)

