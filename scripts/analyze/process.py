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
import shutil, time, datetime
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def parseOutput(outputDir, versions, scenarios, writeCSV=True, filename="summary.csv"):
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
    keywords = {
        "solved": "No solution found",
        "nodes": "of nodes processed",
        "nodes_full_proc": "fully processed",
        "cpu": "Search CPU time",
        "vf_solved": "VF) solved",
        "ub_solved": "UB) solved",
        "vf_time": "solving problem (VF)",
        "ub_time": "solving problem (UB)",
        "objval": "Best solution found had quality",
        "gap": "optimality gap",
        "ul_int_var": "UL Variables (integer)",
        "ll_int_var": "LL Variables (integer)",
        "ul_constr" : "Number of UL Rows",
        "ll_constr" : "Number of LL Rows",
        "align" : "Degree of objective alignment",
        "cg_calls": "Called MIBS cut generator",
        "infeasible": "infeasible",
        "chk_feas_time" : "Checking feasibility",
        "int_idic" : "Improving direction integer calls",
        "frac_idic": "Improving direction fractional calls",
        "ls_calls": "Improving direction Local Search",
        "milp_calls": "Improving direction MILP",
        "ls_cpu": "Local Search CPU time spent",
        "milp_cpu": "MILP CPU time spent",
        "idic_intersect": "Interesection found"
    }

    results = collections.defaultdict(list)
    opt_values = {}
    etol = np.finfo(float).eps

    cpu_time_burnt = 0

    # iterate over versions, scenarios, datasets, and files in each folder to read results
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
            # iterate over different datasets available
            with os.scandir(resultDir) as dataset_it:
                for d_entry in dataset_it:
                    if d_entry.name not in dataSets:
                        # Not a dataset directory
                        continue
                    # iterate over files in the folder
                    with os.scandir(d_entry.path) as output_it:
                        for o_entry in output_it:
                            if o_entry.name.endswith(".out"):
                                print(o_entry.name)
                                # start to write result to the dictionary
                                results["dataset"].append(d_entry.name)
                                results["scenario"].append(scenarios[s])
                                results["version"].append(v)
                                results["instance"].append(
                                    os.path.splitext(o_entry.name)[0]
                                )
                                if results["instance"][-1] not in opt_values:
                                    opt_values[results["instance"][-1]] = np.inf

                                incomplete = True  # mark incomplete output file
                                nosoln = False  # mark no soluntion found
                                # Instance data
                                results['ul_int_var'].append(0)
                                results['ll_int_var'].append(0)
                                results['ul_constr'].append(0)
                                results['ll_constr'].append(0)
                                results['align'].append(0)

                                # Basic data
                                results['incomplete'].append(True)
                                results['solved'].append(False)
                                results['root_bound'].append(10000000)
                                results['100_bound'].append(10000000)
                                results["objval"].append(-1000000)
                                results["gap"].append(1000000)
                                results["cpu"].append(3600)
                                results['nodes'].append(10000000)
                                
                                # CG data
                                results['num_cuts'].append(0)                   # Tot. num cuts
                                results['cg_called'].append(0)                  # CG calls
                                results['cg_failed'].append(0)                  # CG failed: cg_called - num_cuts
                                results['cg_fail_rate'].append(0)               # cg_failed / cg_called
                                results['cg_time'].append(0)                    # CG CPU time
                                # IDICs Data
                                results['num_idic'].append(0)                   # num_int_idic + num_frac_idic
                                results['num_int_idic'].append(0)               # Num Int IDIC
                                results['num_frac_idic'].append(0)              # Num frac IDIC
                                results['idic_called'].append(0)
                                results['idic_failed'].append(0) 
                                results['idic_fail_rate'].append(0) 
                                results['int_idic_called'].append(0)            # Int IDIC CG
                                results['int_idic_failed'].append(0)            # Num CG Int IDIC failed: int_idic_called - num_int_idic
                                results['int_idic_fail_rate'].append(0)
                                results['frac_idic_called'].append(0)           # frac IDIC CG
                                results['frac_idic_failed'].append(0)           # Num CG frac IDIC failed: frac_idic_called - num_frac_idic
                                results['frac_idic_fail_rate'].append(0)
                                # IFDs
                                results['idic_milp_called'].append(0) 
                                results['idic_milp_found'].append(0) 
                                results['idic_milp_failed'].append(0)
                                results['idic_milp_fail_rate'].append(0) 
                                results['idic_ls_called'].append(0)
                                results['idic_ls_found'].append(0)
                                results['idic_ls_failed'].append(0)
                                results['idic_ls_fail_rate'].append(0)
                                results['idic_cg_time'].append(0)               # Tot CPU Time finding IFD
                                results['idic_ls_time'].append(0)               # CPU Time finding IFD with LS
                                results['idic_milp_time'].append(0)             # CPU Time finding IFD with MILP
                                results['idic_intersection_called'].append(0)
                                results['idic_intersection_found'].append(0)
                                results['idic_intersection_fail_rate'].append(0)
        
                                # Feasibility Check
                                results["vf_solved"].append(-1)             
                                results["ub_solved"].append(-1)
                                results["vf_time"].append(0)
                                results["ub_time"].append(0)
                                results['chk_feas_time'].append(0)          # Tot Check Feas time: vf_time + ub_time
                                

                                # read value for each field from file
                                with open(o_entry.path, "r") as file:
                                    for line in file.read().splitlines():
                                        
                                        if (len(line.split()) > 1 and line.split()[0] == '0'):
                                            if len(line.split()) == 4:
                                                results["root_bound"][-1] = float(line.split()[1])
                                            else:
                                                results["root_bound"][-1] = float(line.split()[2])
                                        
                                        elif (len(line.split()) > 0 and line.split()[0] == '100'):
                                            results["100_bound"][-1] = float(line.split()[2])
                                        
                                        elif keywords["solved"] in line:
                                            nosoln = True
                                            results["solved"][-1] = False
                                            print(
                                                "No solution found instance:",
                                                o_entry.name, s
                                            )

                                        elif (keywords["nodes"] in line
                                            or keywords["nodes_full_proc"] in line):
                                            results["nodes"][-1] = int(line.split(":")[1])

                                        elif keywords["cpu"] in line and "Local Search" not in line:
                                            incomplete = False
                                            results["incomplete"][-1] = False
                                            results["cpu"][-1] = float((line.split(":")[1]).split()[0])
                                            cpu_time_burnt += results["cpu"][-1]

                                        elif keywords["vf_solved"] in line:
                                            results["vf_solved"][-1] = int(line.split("=")[1])

                                        elif keywords["ub_solved"] in line:
                                            results["ub_solved"][-1] = int(line.split("=")[1])

                                        elif (keywords["vf_time"] in line and
                                              '1.0-opt' not in versions):
                                            results["vf_time"][-1] = float(line.split("=")[1])

                                        elif (keywords["ub_time"] in line and
                                              '1.0-opt' not in versions):
                                            results["ub_time"][-1] = float(line.split("=")[1])

                                        elif keywords["objval"] in line:
                                            results["objval"][-1] = int(float(line.split()[6]))

                                        elif keywords["gap"] in line:
                                            incomplete = False
                                            if "infinity" in line:
                                                results["gap"][-1] = 1000000
                                                # no soln found
                                                results['solved'][-1] = False
                                            else:
                                                solgap = float(line.split(" ")[5].strip("%\n"))
                                                results["gap"][-1] = solgap
                                                # mark unsolved instances in given time limit
                                                if nosoln == False:
                                                    if solgap - 0.0 < etol:
                                                        results['solved'][-1] = True
                                                    else:
                                                        results['solved'][-1] = False

                                        # elif keywords["chk_feas_time"] in line:
                                        #    results['chk_feas_time'][-1] = float(line.split(' ')[-2])

                                        elif keywords["cg_calls"] in line:
                                            results["num_cuts"][-1] = int(line.split(" ")[8])
                                            results["cg_time"][-1] = float(line.split(" ")[12])
                                            results["cg_called"][-1] = float(line.split(" ")[5])
					
                                        elif keywords["int_idic"] in line:
                                            results["int_idic_called"][-1] = (int(line.split(' ')[-1]))
                                            results["num_int_idic"][-1] = (int(line.split(' ')[-5]))
                                            results["int_idic_failed"][-1] = results["int_idic_called"][-1] - results["num_int_idic"][-1]

                                        elif keywords["frac_idic"] in line:
                                            results["frac_idic_called"][-1] = (int(line.split(' ')[-1]))
                                            results["num_frac_idic"][-1] = (int(line.split(' ')[-5]))
                                            results["frac_idic_failed"][-1] = results["frac_idic_called"][-1] - results["num_frac_idic"][-1]

                                        elif keywords["ls_calls"] in line and 'successful' in line:
                                            results['idic_ls_called'][-1] = int(line.split()[-1])
                                            results['idic_ls_found'][-1] = int(line.split()[-5])
                                            results['idic_ls_failed'][-1] = results['idic_ls_called'][-1] - results['idic_ls_found'][-1]

                                        elif keywords["milp_calls"] in line:
                                            results['idic_milp_called'][-1] = int(line.split()[-1])
                                            results['idic_milp_found'][-1] = int(line.split()[-5])
                                            results['idic_milp_failed'][-1] = results['idic_milp_called'][-1] - results['idic_milp_found'][-1]

                                        elif keywords["ls_cpu"] in line:
                                            results['idic_ls_time'][-1] = float(line.split(':')[-1])
                                        
                                        elif keywords["milp_cpu"] in line:
                                            results['idic_milp_time'][-1] = float(line.split(':')[-1])

                                        elif keywords["idic_intersect"] in line:
                                            results['idic_intersection_called'][-1] = int(line.split()[-1])
                                            results['idic_intersection_found'][-1] = int(line.split()[-5])
                                        
                                        elif keywords["infeasible"] in line:
                                            print("Infeasible instance!")

                                        elif keywords["ul_int_var"] in line:
                                            results['ul_int_var'][-1] = int(line.split(':')[-1])
                                        
                                        elif keywords["ll_int_var"] in line:
                                            results['ll_int_var'][-1] = int(line.split(':')[-1])
                                        
                                        elif keywords["ul_constr"] in line:
                                            results['ul_constr'][-1] = int(line.split(':')[-1])
                                        
                                        elif keywords["ll_constr"] in line:
                                            results['ll_constr'][-1] = int(line.split(':')[-1])

                                        elif keywords["align"] in line:
                                            results['align'][-1] = float(line.split(':')[-1])

                                        # filmosi
                                        elif 'STAT;' in line and len(line.split(';')) > 6:
                                            incomplete=False
                                            results["objval"].append(
                                                float(line.split(";")[2])
                                            )
                                            results["root_bound"][-1] = float(line.split(";")[4])
                                            results["cpu"].append(
                                                float(line.split(";")[5])
                                            )
                                            results["nodes"].append(
                                                float(line.split(";")[8])
                                            )
                                            results["gap"].append(
                                                float(line.split(";")[10])
                                            )
                                            results["vf_solved"].append(0)
                                            results["ub_solved"].append(0)
                                            results["vf_time"].append(0)
                                            results["ub_time"].append(0)
                                            if results["gap"][-1] - 0.0 < etol:
                                                results["solved"].append(True)
                                            else:
                                                results["solved"].append(False)
                                        elif 'User cuts applied:' in line:
                                            results["num_cuts"][-1] = int(line.split()[3])
                                        else:
                                            pass
                                if incomplete:
                                    print(
                                        "Incomplete instance:", o_entry.name, s
                                    )

                                if results["solved"][-1] == False:
                                    results["cpu"][-1] = 3600
                                else:
                                    if (opt_values[results["instance"][-1]] == np.inf):
                                        opt_values[results["instance"][-1]] = results["objval"][-1]
                                    elif opt_values[results["instance"][-1]] != results["objval"][-1]:
                                        print("[%s] Warning: objective values don't agree!" % (results["instance"][-1]))

                                if results["cpu"][-1] < 0.01:
                                    print ("Bad value! ", results["cpu"][-1])
                                    results["cpu"][-1] = .01

    for k in results:
      print(k)
      print(len(results[k]))

    print("")
    print("Congrats! You have burnt ", str(datetime.timedelta(seconds=cpu_time_burnt)), " CPU time!")
    print("")
    df_result = pd.DataFrame(results)

    # print("UL integer vars:", df_result["ul_int_var"].min(), df_result["ul_int_var"].max())
    # print("LL integer vars:", df_result["ll_int_var"].min(), df_result["ll_int_var"].max())
    # print("UL Constraints:", df_result["ul_constr"].min(), df_result["ul_constr"].max())
    # print("LL Constraints:", df_result["ll_constr"].min(), df_result["ll_constr"].max())
    # print("Align:", df_result["align"].min(), df_result["align"].max())


    # make some adjustment to formats
    # display check feasibility time as % of search time?
    # sum vf+ub time -> feasibility time (or read from output directly?)
    if '1.0-opt' not in versions:
        df_result["chk_feas_time"] = df_result["ub_time"] + df_result["vf_time"]
        df_result["chk_feas_time"] = df_result["chk_feas_time"].astype(float).round(2)
    #df_result["cpu"] = df_result["cpu"].astype(float).round(2)

    # Tot CPU time finding IFDs
    df_result["idic_called"] = df_result["int_idic_called"] + df_result["frac_idic_called"]
    df_result["idic_cg_time"] = df_result["idic_ls_time"] + df_result["idic_milp_time"]
    df_result["idic_cg_avg_time"] = df_result["idic_cg_time"] / df_result["idic_called"]
    df_result["idic_cg_avg_time"].fillna(0, inplace=True)
    # with pd.option_context('display.max_rows', None,
    #                           'display.max_columns', None,
    #                           'display.precision', 3,
    #     ):
        # print(df_result["idic_ls_time"])
        # print(df_result["idic_milp_time"])
        # print(df_result["idic_called"])
        # print(df_result["idic_cg_time"])
        # print(df_result["idic_cg_avg_time"])
        # exit(0)

    # compute percentages
    df_result["cg_fail_rate"] = df_result["cg_failed"] / df_result["cg_called"][df_result["cg_called"] > etol]
    df_result["cg_fail_rate"].fillna(0, inplace=True)
    df_result["idic_failed"] = df_result["int_idic_failed"] + df_result["idic_ls_failed"]
    df_result["idic_fail_rate"] = df_result["idic_failed"] / df_result["idic_called"][df_result["idic_called"] > etol]
    df_result["idic_fail_rate"].fillna(0, inplace=True)
    df_result["int_idic_fail_rate"] = df_result["int_idic_failed"] / df_result["int_idic_called"][df_result["int_idic_called"] > etol]
    df_result["int_idic_fail_rate"].fillna(0, inplace=True)
    df_result["frac_idic_fail_rate"] = df_result["frac_idic_failed"] / df_result["frac_idic_called"][df_result["frac_idic_called"] > etol]
    df_result["frac_idic_fail_rate"].fillna(0, inplace=True)
    df_result["idic_milp_fail_rate"] = df_result["idic_milp_failed"] / df_result["idic_milp_called"][df_result["idic_milp_called"] > etol]
    df_result["idic_milp_fail_rate"].fillna(0, inplace=True)
    df_result["idic_ls_fail_rate"] = df_result["idic_ls_failed"] / df_result["idic_ls_called"][df_result["idic_ls_called"] > etol]
    df_result["idic_ls_fail_rate"].fillna(0, inplace=True)
    df_result["idic_intersection_fail_rate"] = (df_result["idic_intersection_called"] - df_result["idic_intersection_found"] ) \
                                                / df_result["idic_intersection_called"][df_result["idic_intersection_called"] > etol]
    df_result["idic_intersection_fail_rate"].fillna(0, inplace=True)

    # write results to .csv file
    if writeCSV:
        # df_result.to_csv(filename, mode='a', header=False, index=False) # append results only
        df_result.to_csv(filename, index=False)

    return df_result


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
                    rsltDict[inst].update(
                        {(scn, v, 'all', col): df_temp[col].values[0] for col in displayCols}
                )

    # convert dict to structured df: change to formal column names?
    df_forprint = pd.DataFrame.from_dict(rsltDict, orient="index")
    df_forprint.columns.names = ["scn", "v", "datasets", "fields"]
    # df_forprint = df_forprint.sort_index()

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
    # df_time = pd.to_numeric(df_time, errors='coerce').replace(np.nan, 36000)
    # filter out cases where time is < 5'' or > 3600'' for all methods
    col_list = df_time.columns.values.tolist()

    drop_easy = df_time[(df_time[col_list] < 5).all(axis=1)].index.tolist()
    drop_small_time = df_time[(df_time[col_list] <= 0.01).any(axis=1)].index.tolist()
    drop_unsolved = df_solved[(df_solved[col_list] != True).all(axis=1)].index.tolist()
    drop_list_time = list(set(drop_easy) | set(drop_unsolved) | set(drop_small_time))
    #drop_list_time.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    #print(drop_easy)
    #print(drop_small_time)
    #print(drop_unsolved)
    df_solved = df.drop(drop_list_time)

    df_gap = df.xs(
        (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    drop_no_gap = df_gap[(df_gap[col_list] >= 1000000).all(axis=1)].index.tolist()
    drop_list_gap = list(drop_no_gap)
    #drop_list_gap.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    #print(drop_list_gap)
    df_has_soln = df.drop(drop_list_gap)

    return df_solved, df_has_soln


def plotPerfProf(
        df, plotname="perf_profile", plottitle="Performance Profile", plotformat='png',
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

    num_lines = len(scenarios)
    cmap = cm.get_cmap(colorPalette, num_lines)
    i = 0

    fig, ax = plt.subplots(1, 1)

    # if given legend name len != col #, use defualt column name
    if legendnames and (len(legendnames) != len(df.columns)):
        legendnames = {}

    # find min value in the dataframe
    col_list = df.columns.values.tolist()
    df["virtual_best"] = df[col_list].min(axis=1)

    #print(df["virtual_best"])
    
    for col in col_list:
        # print(col_list)
        # print(col)
        # for each col, compute ratio
        ratios = df[col] / df["virtual_best"]
        with pd.option_context('display.max_rows', None,
                              'display.max_columns', None,
                              'display.precision', 3,
        ):
            print(df[col])
        #    print(ratios.sort_values())
        #    print(ratios)
        uniq_ratios = ratios.unique()
        uniq_ratios.sort()  # sort in place

        #print(uniq_ratios)
        cum_cnt = np.sum(np.array([ratios <= ur for ur in uniq_ratios]), axis=1)
        cum_frac = cum_cnt / len(ratios)
        #print(cum_frac)

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
            plt.plot(x_val, y_val, label=legendnames[col], color=cmap(i))
        elif versionlegend:
            plt.plot(x_val, y_val, label=col, color=cmap(i))  # , color=colors[i])
        else:
            plt.plot(x_val, y_val, label=col[0], color=cmap(i))  # , color=colors[i])
        i += 1

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
    fig.savefig(plotname + '.' + plotformat, dpi=fig.dpi)


def plotCumProf(df, plotname="cum_profile", plottitle = "Cumulative Profile", plotformat='png',
                legendnames={}, versionlegend=False):

    fig = plt.figure()
    gs = fig.add_gridspec(1, 2, wspace=0)
    ax = gs.subplots(sharey=True)

    df_time = df.xs(
        (ds, "cpu"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    df_gap = df.xs(
        (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()

    num_lines = len(scenarios)
    cmap = cm.get_cmap(colorPalette, num_lines)
    i = 0

    col_list = df_time.columns.values.tolist()
    time_buckets = range(0, 3600)

    for col in col_list:
        # print(col)
        times = df_time[col]
        # print(times)

        cum_cnt = np.sum(np.array([times <= t for t in time_buckets]), axis=1)
        cum_frac = cum_cnt / len(df)
        #print(cum_frac)
        x_val = []
        if legendnames:
            ax[0].plot(time_buckets, cum_frac, label=legendnames[col], color=cmap(i))
        elif versionlegend:
            ax[0].plot(time_buckets, cum_frac, label=col, color=cmap(i))  # , color=colors[i])
        else:
            ax[0].plot(time_buckets, cum_frac, label=col[0], color=cmap(i))
        i += 1

    ax[0].set_xlim(0, 3599)
    ax[0].set_ylim(0.0, 1)
    ax[0].tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax[0].set_xlabel("Time")
    ax[0].set_ylabel("Fraction of instances")

    gap_buckets = np.linspace(0, 100, 1000)

    i = 0
    for col in col_list:

        # print(col)
        gaps = df_gap[col]
        # print(gaps)

        cum_cnt = np.sum(np.array([gaps <= g for g in gap_buckets]), axis=1)
        cum_frac = cum_cnt / len(df_gap)
        #print(cum_frac)
        x_val = []
        if legendnames:
            ax[1].plot(gap_buckets, cum_frac, label=legendnames[col], color=cmap(i))
        elif versionlegend:
            ax[1].plot(gap_buckets, cum_frac, label=col, color=cmap(i))  # , color=colors[i])
        else:
            ax[1].plot(gap_buckets, cum_frac, label=col[0], color=cmap(i))
        i += 1

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
    fig.savefig(plotname + '.' + plotformat, dpi=fig.dpi)
    # fig.savefig("./performance/barchart/"+plotname+'.eps', format='eps', dpi=600)

def plotBaselineProf(
        df, baseline, plotname="base_profile", plottitle="Baseline Profile", plotformat='png',
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

    num_lines = len(scenarios) - 1
    cmap = cm.get_cmap(colorPalette, num_lines)
    i = 0

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
                ax[0].plot(x_val, y_val, label=legendnames[col], color=cmap(i))
            elif versionlegend:
                ax[0].plot(x_val, y_val, label=col, color=cmap(i))  # , color=colors[i])
            else:
                ax[0].plot(x_val, y_val, label=col[0], color=cmap(i))  # , color=colors[i])

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
            ax[1].plot(x_val, y_val, label=legendnames[col], color=cmap(i))
        elif versionlegend:
            ax[1].plot(x_val, y_val, label=col, color=cmap(i))  # , color=colors[i])
        else:
            ax[1].plot(x_val, y_val, label=col[0], color=cmap(i))  # , color=colors[i])
        i += 1

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

    fig.supxlabel("Ratio of baseline (%s)" % baseline[0])
    fig.supylabel("Fraction of instances")
    fig.suptitle(plottitle)
    fig.tight_layout()
    fig.savefig(plotname + '.' + plotformat, dpi=fig.dpi)

def plotBaselineProfSingle(
        df, baseline, plotname="base_profile", plottitle="Baseline Profile", plotformat='png',
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
            plt.plot(x_val, y_val, label=col)  # , color=colors[i])
        else:
            plt.plot(x_val, y_val, label=col[0])  # , color=colors[i])

    # set plot properties
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.05)
    ax.tick_params(axis="both", direction="in", right=True)

    # set other figure elements
    ax.set_title(plottitle)
    ax.set_xlabel("Ratio of baseline")
    ax.set_ylabel("Fraction of instances")
    ax.legend(
        loc="upper right",
        #bbox_to_anchor=(0.9, 0.05),
        markerscale=1.25,
        frameon=True,
        labelspacing=0.35,
        fontsize="x-small",
    )

    fig.tight_layout()
    fig.savefig(plotname + '.' + plotformat, dpi=fig.dpi)

if __name__ == "__main__":

    # Datasets: names correspond to the key in myrun.py instanceDirs
    dataSets = [
            'all',
            # "interKpShi",
            # "interDen"
            "iblpDen",
            "iblpDen2",
            "iblpZhang",
            "iblpZhang2",
            'iblpFis'
        ]

    # Version: names correspond to versions in myrun.py
    versions = ['idBC']

    # Where the results are
    outputDir = ["/home/feb223/results/MibS"]

    # Scenarios
    #     keys: correspond to the keys of mibsParamsInputs in myrun.py
    #     values: is the name to use in plots
    scenarios = {
        'idBC-LS-k_2-dBnd_Inf_fracB'   : "idB&C-LS-k_2 (frac)",
        "idBC-LS-k_3-dBnd_Inf_fracB"   : "idB&C-LS-k_3 (frac)",
        # 'idBC-LS-k_2-dBnd_0_10_fracB'  : 'idB&C-LS-k_2-dBnd_0_10 (frac)',
        'idBC-LS-k_2-dBnd_10_Inf_fracB': 'idB&C-LS-k_2-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_0_10_fracB'  : 'idB&C-LS-k_3-dBnd_0_10 (frac)',
        'idBC-LS-k_3-dBnd_10_Inf_fracB': 'idB&C-LS-k_3-dBnd_10_Inf (frac)',
        # 'idBC-LS-k_3-dBnd_0_8_fracB'   : 'idB&C-LS-k_3-dBnd_0_8 (frac)',
        # 'idBC-LS-k_3-dBnd_8_Inf_fracB' : 'idB&C-LS-k_3-dBnd_8_Inf (frac)',
        # 'idBC-LS-k_3-dBnd_0_12_fracB'  : 'idB&C-LS-k_3-dBnd_0_12 (frac)',
        'idBC-LS-k_3-dBnd_12_Inf_fracB': 'idB&C-LS-k_3-dBnd_12_Inf (frac)',
        # 'idBC-LS-k_4-dBnd_Inf_fracB'   : 'idB&C-LS-k_4 (frac)',
        # 'idBC-LS-k_4-dBnd_0_10_fracB'  : 'idB&C-LS-k_4-dBnd_0_10 (frac)',
        'idBC-LS-k_4-dBnd_10_Inf_fracB': 'idB&C-LS-k_4-dBnd_10_Inf (frac)',
        # 'idBC-LS-k_5-dBnd_Inf_fracB'   : 'idB&C-LS-k_5 (frac)',
        # 'idBC-LS-k_5-dBnd_0_10_fracB'  : 'idB&C-LS-k_5-dBnd_0_10 (frac)',
        'idBC-LS-k_5-dBnd_10_Inf_fracB': 'idB&C-LS-k_5-dBnd_10_Inf (frac)',
        "idBC-MILP_fracB"              : "idB&C-MILP (frac)",
        'idBC-k_2-MILP_fracB'          : 'idB&C-MILP-k_2 (frac)',
        # 'idBC-k_3-MILP_fracB'          : 'idB&C-MILP-k_3 (frac)',
        'idBC-k_4-MILP_fracB'          : 'idB&C-MILP-k_4 (frac)',
        # 'idBC-k_5-MILP_fracB'          : 'idB&C-MILP-k_5 (frac)',
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
        "MibS_1_2_defaultB"            : "MibS (default)",
        "MibS_1_2_IDIC_defaultB"       : "MibS IDIC-MILP (default)",
        # "MibS_1_2-LS-k_3_defaultB"     : "MibS IDIC-LS-k_3 (default)",
        "MibS_1_2-LS-k_2_defaultB"     : "MibS IDIC-LS-k_2 (default)"
        }

    ################# Process & Save | Load from CSV ###################
    # Datasets name
    # ds_name = "INT-DEN-SHI"
    ds_name = "F-D-D2-Z-Z2"

    # specify summary file name
    file_csv_out = "summary_" + ds_name + ".csv"
    #file_csv_in = "summary-1.2.1.csv"
    file_csv_in = "summary_idBC_all_iblp.csv"
    file_csv_in = "summary_" + ds_name + ".csv"
    
    # if len(args) == 0:
    if 0:
        df_r = parseOutput(
            outputDir, versions, scenarios, writeCSV=True, filename=file_csv_out
        )
    else:
        try:
            df_r = pd.read_csv(file_csv_in)
            set_cond = (df_r["scenario"].isin(scenarios.values())) | (
                df_r["dataset"].isin(dataSets)
            )
            df_r = df_r[set_cond]
        except FileNotFoundError:
            print("{} does not exist in current directory.".format(file_csv_in))
        else:
            print("Reading from", file_csv_in)

    ################### Format Data & Print Table ####################
    # specify txt file name to print tables in LATEX
    file_txt = "ltx_tb_cut.txt"

    # columns to process and print
    # Columns to process and print
    displayCols = {
        # Basic data
        'incomplete': "Incomplete",
        'solved': "Solved",
        'root_bound': 'Bound at root',
        '100_bound': 'Bound after 100 nodes',
        "objval": 'Object Value',
        "gap": "Final Gap",
        "cpu": "CPU Search Time",
        'nodes': "Number of Processed Nodes",
        
        # CG data
        'num_cuts': "Total number of cuts",
        # 'cg_called': "Total number of CG calls",
        # 'cg_failed': "CG calls failed",
        'cg_fail_rate': "CG calls fail rate (%)",
        'cg_time': "CG CPU time",
        # IDICs data
        'num_idic': "Total number of IDICs",
        'num_int_idic' : "Total number of int IDICs",
        'num_frac_idic': "Total number of frac IDICs",
        'idic_fail_rate': "IDIC CG calls fail rate (%)",
        'int_idic_fail_rate': "Int IDIC CG calls fail rate (%)",
        'frac_idic_fail_rate': "Frac IDIC CG calls fail rate (%)",
        # IFDs data
        'idic_ls_fail_rate': "Local Search fail rate (%)",
        'idic_cg_time': "IFDs CG CPU Time",
        'idic_cg_avg_time': "Per-call IFDs CG CPU Time"

        # Feasibility Check
        # "vf_solved": 'Number of VF problem solved',      
        # "ub_solved": 'Number of UB problem solved',
        # "vf_time": "VF problem CPU time",
        # "ub_time": "UB problem CPU time",
        # 'chk_feas_time': 'Check Feasibility Time'
    }

    df_proc = processTable(df_r, displayCols, writeLTX=False, filename=file_txt)

    ################### Make Performance Profile ####################
    # columns to compare in the plot
    plotCols = {
        "cpu": ["CPU Time", 40],
        # "cpu": ["CPU Time", 20],
        "nodes": ["Nodes Processed", 15],
        'cg_time': ["Finding IFDs total CPU Time", 50],
        # 'cg_time': ["Cut Generation CPU Time", 30],
        # 'idic_fail_rate': ["IFDs CG calls fail rate (%)", 50],
        'idic_cg_avg_time': ["Finding IFDs average CPU Time", 60]
    }

    colorPalette = 'tab20'
    # plotCols = {}

    # manual input example:
    # for k in scenarios:
    #     if '01' in k:
    #         scenarios[k] = 'linkingBranchStrategy'
    #     else:
    #         scenarios[k] = 'fractionalBranchStrategy'


    baseline = None 
    # baseline = ("MibS only IDICs (frac)", "idBC")
    baseline = ("MibS (default)", "idBC")
    # baseline = ("MibS IDIC-MILP (default)", "idBC")

    # baseline = ("MibS_IDIC", 'ipco')
    #baseline = ("Type1IC", "1.2-opt")
    #baseline = ('GenNoGood+Type1+IntNoGood (link)', '1.2-opt')
    #baseline = ('Watermelon (frac+LV)', '1.2-opt')
    #baseline = ('FracWatermelon (frac)', '1.2-opt')
    #baseline = ('Benders Interdict (link)', '1.2.1-final')
    if len(versions) > 1:
        versionlegend = True
    else:
        versionlegend = False
    
    # plotformat = 'pdf'
    plotformat = 'png'
    
    dataSets = ['all']

    for ds in dataSets:
        df_solved, df_has_soln = dropFilter(df_proc, scenarios, ds)
        # print(df_solved)

        plotCumProf(df_has_soln, plotname=("cum_" + ds_name).replace(' ', '_'),
                    plottitle="Cumulative Profile: Time-Gap ("+ds_name+")",
                    versionlegend = versionlegend, plotformat=plotformat
        )

        if baseline is not None: 
            print("")
            print("Creating baseline profile for gap")
            print("")
            df_gap = df_has_soln.xs(
                (dataSets[0], "gap"), level=["datasets", "fields"], axis=1, drop_level=True
            ).copy()
            print(len(df_gap))
            df_baseline_has_gap = df_gap.drop(df_gap[df_gap[baseline] == 0].index.to_list())
            plotBaselineProf(
                df_baseline_has_gap, baseline = baseline,
                plotname=("base_" + baseline[0] + "_" + "gap_" + ds_name).replace(' ', '_'),
                plottitle = "Baseline Profile: Gap ("+ds_name+")",
                xmax=25, versionlegend = versionlegend, plotformat=plotformat
            )

        for col in plotCols:
            if col != "root_gap":
                df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
            else:
                df_sub = df_has_soln.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
            print("")
            print("Creating performance profile for " + col , ", num instances: ", len(df_sub))
            print("")
            plotPerfProf(
                df_sub, plotname=("perf_" + col + "_" + ds_name).replace(' ', '_'),
                plottitle = "Performance Profile: "+plotCols[col][0]+" ("+ds_name+")",
                xmin = 0.0, xmax=plotCols[col][1],
                versionlegend = versionlegend, plotformat=plotformat
            )
            if baseline is not None: 
                print("")
                print("Creating baseline profile for "+col)
                print("")
                # plotBaselineProfSingle(
                #     df_sub, baseline = baseline,
                #     plotname="base_"+baseline[0]+"_"+col+"_"+ds,
                #     plottitle = "Baseline Profile: "+plotCols[col][0]+" ("+ds+")",
                #     xmax=plotCols[col][1],
                #     versionlegend = versionlegend
                # )
                plotBaselineProf(
                    df_sub, baseline = baseline,
                    plotname=("base_"+baseline[0]+"_"+col+"_"+ds_name).replace(' ', '_'),
                    plottitle = "Baseline Profile: "+plotCols[col][0]+" ("+ds_name+")",
                    xmax=plotCols[col][1],
                    versionlegend = versionlegend, plotformat=plotformat
                )
