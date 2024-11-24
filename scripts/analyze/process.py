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


def parseOutput(outputDir, versions, scenarios, writeCSV=True,
                filename="summary.csv", name=''):
    """
    The function parse the output file in the given directory.
    Assume the subfolders hierarchy: 
    outputDir/param_scenario_name/testset_name/BR_Output/file.out.
    The result will also be written to a .csv file if not specified.
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
        "num_cuts": "Called MIBS cut generator",
        "infeasible": "infeasible",
        "int_idic" : "integer IDICs",
        "int_lv_idic" : "integer (LV) IDICs",
        "frac_idic" : "fractional IDICs",
        "int_idic_failed": "IDIC cut generation failed:",
        "int_lv_idic_failed": "cut generation failed (LV)",
        "frac_idic_failed": "Fractional IDIC cut generation failed"
    }

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
            # iterate over different datasets available
            with os.scandir(resultDir) as dataset_it:
                for d_entry in dataset_it:
                    if d_entry.name not in dataSets:
                        continue
                    # iterate over files in the folder
                    with os.scandir(d_entry.path) as output_it:
                        for o_entry in output_it:
                            if o_entry.name.endswith(".out"):
                                #if (o_entry.name.startswith("K") or
                                #    o_entry.name.startswith("bmilp")):
                                #    print("Skipping!!!!")
                                #    continue
                                # start to write result to the dictionary
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

                                incomplete = True  # mark incomplete output file
                                nosoln = False  # mark no soluntion found
                                results['num_cuts'].append(0)
                                results['cut_time'].append(0)
                                results['root_bound'].append(10000000)
                                results['100_bound'].append(10000000)
                                results['num_int_idic'].append(0)
                                results['num_int_lv_idic'].append(0)
                                results['num_frac_idic'].append(0)
                                results['num_int_idic_fail'].append(0)
                                results['num_int_lv_idic_fail'].append(0)
                                results['num_frac_idic_fail'].append(0)
                                results['cg_called'].append(0)
                                results['cg_failed'].append(0)
                                results['cg_fail_rate'].append(0)
                                results['num_idic'].append(0)
                                results['depth_idic'].append(-1)
                                results["vf_solved"].append(-1)
                                results["ub_solved"].append(-1)
                                if ('1.0.0-opt' not in versions and
                                    'filmosi' not in versions): 
                                    results["vf_time"].append(-1)
                                    results["ub_time"].append(-1)
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
                                            print(
                                                "No solution found instance:",
                                                o_entry.name, s
                                            )

                                        elif (
                                            keywords["nodes"] in line
                                            or keywords["nodes_full_proc"] in line
                                        ):
                                            results["nodes"].append(
                                                int(line.split(":")[1])
                                            )

                                        elif keywords["cpu"] in line:
                                            results["cpu"].append(
                                                float((line.split(":")[1]).split()[0])
                                            )

                                        elif keywords["vf_solved"] in line:
                                            results["vf_solved"][-1] = int(line.split("=")[1])

                                        elif keywords["ub_solved"] in line:
                                            results["ub_solved"][-1] = int(line.split("=")[1])

                                        elif (keywords["vf_time"] in line and
                                              '1.0.0-opt' not in versions and
                                              'filmosi' not in versions):
                                            results["vf_time"][-1] = float(line.split("=")[1])

                                        elif (keywords["ub_time"] in line and
                                              '1.0.0-opt' not in versions and
                                              'filmosi' not in versions):
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
                                                solgap = float(
                                                    line.split(" ")[5].strip("%\n")
                                                )
                                                results["gap"][-1] = solgap
                                                # mark unsolved instances in given time limit
                                                if nosoln == False:
                                                    if solgap - 0.0 < etol:
                                                        results["solved"].append(True)
                                                    else:
                                                        results["solved"].append(False)
                                        elif keywords['frac_idic'] in line:
                                            results['num_frac_idic'][-1] += 1
                                            results['num_idic'][-1] += 1
                                            depth = int(line.split()[4])
                                            results['depth_idic'][-1] = ((
                                                results['depth_idic'][-1]*(
                                                    results['num_idic'][-1] - 1
                                                    ) + depth)/results['num_idic'][-1])
                                        elif keywords['int_lv_idic'] in line:
                                            results['num_int_lv_idic'][-1] += 1
                                            results['num_idic'][-1] += 1
                                            depth = int(line.split()[5])
                                            results['depth_idic'][-1] = ((
                                                results['depth_idic'][-1]*(
                                                    results['num_idic'][-1] - 1
                                                    ) + depth)/results['num_idic'][-1])
                                        elif keywords['int_idic'] in line:
                                            results['num_int_idic'][-1] += 1
                                            results['num_idic'][-1] += 1
                                            depth = int(line.split()[4])
                                            results['depth_idic'][-1] = ((
                                                results['depth_idic'][-1]*(
                                                    results['num_idic'][-1] - 1
                                                    ) + depth)/results['num_idic'][-1])
                                        elif keywords['int_idic_failed'] in line:
                                            results['num_int_idic_fail'][-1]+=1
                                            results['cg_failed'][-1] += 1
                                        elif keywords['int_lv_idic_failed'] in line:
                                            results['num_int_lv_idic_fail'][-1]+=1
                                            results['cg_failed'][-1] += 1
                                        elif keywords['frac_idic_failed'] in line:
                                            results['num_frac_idic_fail'][-1]+=1
                                            results['cg_failed'][-1] += 1
                                        #elif keywords['ul_int_var'] in line:
                                        #    results['ul_int_var'].append(int(line.split(':')[1]))
                                        #elif keywords['ll_int_var'] in line:
                                        #    results['ll_int_var'].append(int(line.split(':')[1]))
                                        elif keywords["num_cuts"] in line:
                                            found = True
                                            results["num_cuts"][-1] = int(line.split(" ")[8])
                                            results["cut_time"][-1] = float(line.split(" ")[12])
                                            results["cg_called"][-1] = float(line.split(" ")[5])
                                        elif (keywords["infeasible"] in line and
                                              results["version"][-1] != 'filmosi'):
                                            print("Infeasible instance!")
                                        elif 'STAT;' in line and len(line.split(';')) > 6:
                                            incomplete=False
                                            results["objval"][-1] = \
                                                float(line.split(";")[2])
                                            results["root_bound"][-1] = float(line.split(";")[4])
                                            results["cpu"].append(
                                                float(line.split(";")[5])
                                            )
                                            results["nodes"].append(
                                                float(line.split(";")[8])
                                            )
                                            results["gap"][-1] = \
                                                float(line.split(";")[10])
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
                                if incomplete:
                                    results["nodes"].append(1000000000)
                                    results["cpu"].append(3600)
                                    results["solved"].append(False)
                                    print(
                                        "Incomplete instance:", o_entry.name, s
                                    )

                                if results["solved"][-1] == False:
                                    results["cpu"][-1] = 3600
                                else:
                                    if (opt_values[results["instance"][-1]] == np.inf):
                                        opt_values[results["instance"][-1]] = results["objval"][-1]
                                    elif opt_values[results["instance"][-1]] != results["objval"][-1]:
                                        print("************ Warning: objective values don't agree!")
                                        print("************ ",
                                              results["instance"][-1],
                                              opt_values[results["instance"][-1]],
                                              ' != ',
                                              results["objval"][-1])
                                        
                                if results["cpu"][-1] < 0.01:
                                    print ("Small value! ", results["cpu"][-1], results["instance"][-1])
                                    #results["cpu"][-1] = .01
                                if results['depth_idic'][-1] == 0:
                                    results['depth_idic'][-1] = 0.1
                                a = [len(results[k]) for k in results]
                                b = [k for k in results]
                                for i in range(len(a)):
                                    if a[i] != a[0]:
                                        print(o_entry.path+"!!!!!!!!!!!!!", b[i])
                                    
                                    
    #for k in results:
    #   print (k)
    #   print(len(results[k]))
    #print (results['num_frac_idic'])
    #print (results['num_int_idic'])
    #print (results['num_idic'])
    #print (results['depth_idic'])
    #print(results['root_bound'])
    #print(results['instance'])
    for i in range(len(results["instance"])):
        if results['cg_called'][i] > 0:
            results['cg_fail_rate'][i] = results['cg_failed'][i]/results['cg_called'][i]
        if (opt_values[results["instance"][i]] not in [0, np.inf] and
            results['root_bound'][i] != 10000000):
            results["root_gap"].append(round(100*abs((opt_values[results["instance"][i]] -
                                                      results["root_bound"][i])/opt_values[results["instance"][i]]),
                                                     2))
        else:
            results["root_gap"].append(100000)
        #print(results["instance"][i], results["scenario"][i],
        #      results["root_bound"][i], opt_values[results["instance"][i]])
        if (opt_values[results["instance"][i]] not in [0, np.inf] and
            results['100_bound'][i] != 10000000):
            results["100_gap"].append(round(100*abs((opt_values[results["instance"][i]] -
                                                     results["100_bound"][i])/opt_values[results["instance"][i]]),
                                                     2))
        else:
            results["100_gap"].append(100000)
        #print(results["instance"][i], results["scenario"][i],
        #      results["100_bound"][i], opt_values[results["instance"][i]])
    df_result = pd.DataFrame(results)

    # make some adjustment to formats
    # display check feasibility time as % of search time?
    # sum vf+ub time -> feasibility time (or read from output directly?)
    if '1.0.0-opt' not in versions and 'filmosi' not in versions:
        df_result["chk_feas_time"] = df_result["ub_time"] + df_result["vf_time"]
        df_result["chk_feas_time"] = df_result["chk_feas_time"].astype(float).round(2)
    #df_result["cpu"] = df_result["cpu"].astype(float).round(2)

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

    # convert dict to structured df: change to formal column names?
    df_forprint = pd.DataFrame.from_dict(rsltDict, orient="index")
    print (df_forprint)
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
    # df_time = pd.to_numeric(df_time, errors='coerce').replace(np.nan, 36000)
    # filter out cases where time is < 5'' or > 3600'' for all methods
    col_list = df_time.columns.values.tolist()

    if 1: #ds == 'INTERD-DEN':
        drop_easy = df_time[(df_time[col_list] < 1).all(axis=1)].index.tolist()
        drop_small_time = df_time[(df_time[col_list] <= 0.01).any(axis=1)].index.tolist()
        drop_unsolved = df_solved[(df_solved[col_list] != True).all(axis=1)].index.tolist()
        drop_list_time = list(set(drop_easy) | set(drop_unsolved) | set(drop_small_time))
    else:
        drop_easy = df_time[(df_time[col_list] < 5).all(axis=1)].index.tolist()
        drop_unsolved = df_solved[(df_solved[col_list] != True).all(axis=1)].index.tolist()
        drop_list_time = list(set(drop_easy) | set(drop_unsolved))

    ##drop_list_time.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    ##print(drop_easy)
    ##print(drop_small_time)
    ##print(drop_unsolved)
    df_solved = df.drop(drop_list_time)
    #print(df_solved)

    df_gap = df.xs(
        (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
    ).copy()
    drop_no_gap = df_gap[(df_gap[col_list] >= 1000000).all(axis=1)].index.tolist()
    drop_list_gap = list(drop_no_gap)
    #drop_list_gap.extend(["cap6000-0.100000","cap6000-0.500000","cap6000-0.900000"])
    #print(drop_list_gap)
    df_has_soln = df.drop(drop_list_gap)
    #print(df_has_soln)

    return df_solved, df_has_soln


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
        print(col)
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


def plotCumProf(df, versions, plotname="cum_profile",
                plottitle = "Cumulative Profile",
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

    col_list = df_time.columns.values.tolist()
    time_buckets = range(0, 3600)

    for col in col_list:
        #print(col)
        times = df_time[col]
        #print(times)
        cum_cnt = np.sum(np.array([times <= t for t in time_buckets]), axis=1)
        cum_frac = cum_cnt / len(df)
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

if __name__ == "__main__":

    dataSets = [
        #'MIBLP-XU',
        "IBLP-FIS",
        'INTERD-DEN',
        'IBLP-DEN',
        'IBLP-DEN2',
        'IBLP-ZHANG',
        'IBLP-ZHANG2',
        #'BENCHMARK'
        #'all'
    ]
    aggregate = True
#    aggregate = False
    name = ''

    if aggregate:
        for i in dataSets:
            if i == 'MIBLP-XU':
                name = name + 'X-'
            elif i == "IBLP-FIS":
                name = name + 'F-'
            elif i == 'INTERD-DEN':
                name = name + 'INT-'
            elif i == 'IBLP-DEN':
                name = name + 'D-'
            elif i == 'IBLP-DEN2':
                name = name + 'D2-'
            elif i == 'IBLP-ZHANG':
                name = name + 'Z-'
            elif i == 'IBLP-ZHANG2':
                name = name + 'Z2-'
            elif i == 'BENCHMARK':
                name = name + 'BENCH-'
        name = name[:len(name)-1]
    

    versions = {
        #'filmosi':'filmosi',
        #'1.0.0-opt':'1.0.0',
        #'1.1.3-opt':'1.1.3',
        #'1.2.0-opt':'1.2.0',
        #'1.2.1-opt':'1.2.1',
        #'1.2.1-cplex-opt':'1.2.1-cplex',
        '1.2.2-opt':'1.2.2',
        #'1.2.2-opt-cg-fail':'1.2.2'
    }
    
    # Output parent path
    outputDir = ["/home/ted/Projects/MibS/output"]
    #outputDir = ["/home/ted/Projects/MibS/output-Mac"]

    scenarios = {
        #### Interdiction
        #'noCut' : "No Cuts (link)", 
        #'fracISICType2-link': 'Frac ISIC Type 2 (link)',
        #'ISICType2-link': 'ISIC Type 2 (link)',
        #'fracIDIC-link': 'Frac IDIC (link)',
        #'IDIC-link': 'IDIC (link)',
        #'bendersInterdiction-frac': 'Benders Interdict (frac)',
        #'bendersInterdiction-link': 'Benders Interdict (link)',

        #### Pure Integer
        #'noCut' : "No Cuts (link)",
        #'fracIDIC+ISICType1-frac' : 'Frac IDIC + ISIC Type 1 (frac)',
        #'ISICType1-frac': 'ISIC Type 1 (frac)', 
        #'fracISICType2-frac': 'Frac ISIC Type 2 (frac)',
        #'ISICType2-frac': 'ISIC Type 2 (frac)',
        #'fracIDIC-frac': 'Frac IDIC (frac)',    
        #'IDIC-frac': 'IDIC (frac)',
        #'hyper-frac': 'Hypercube IC (frac)',    
        #'intNoGood-frac': 'Integer No Good (frac)',

        #### Pure Binary
        # 'noCut' : "No Cuts (link)", 
        # 'ISICType1-frac': 'ISIC Type 1 (frac)', 
        # 'fracISICType2-frac': 'Frac ISIC Type 2 (frac)',
        # 'ISICType2-frac': 'ISIC Type 2 (frac)',
        # 'fracIDIC-frac': 'Frac IDIC (frac)',    
        # 'IDIC-frac': 'IDIC (frac)',
        # 'hyper-frac': 'Hypercube IC (frac)',    
        # 'intNoGood-frac': 'Integer No Good (frac)',
        # 'bendersBinary-frac': 'Benders Binary (frac)', 
        # 'genNoGood-frac': 'Generalized No Good (frac)', 

        ###### filmosi

        #'default' : 'Default',
        #'default2' : 'Default 2',
        #'default2' : 'Default',
        #'default-MIP' : 'Default (no MIP cuts)',
        #'default+parallel' : 'Default (parallel)'

        ###### 1.2.2-opt

        #'LIntIDIC-frac' : 'LIntIDIC-frac',
        #'LIntIDIC-link' : 'LIntIDIC-link',
        #'XYIntIDIC-frac' : 'XYIntIDIC-frac',
        #'XYIntIDIC-link' : 'XYIntIDIC-link',
        #'YIntIDIC-frac' : 'YIntIDIC-frac',
        #'YIntIDIC-link' : 'YIntIDIC-link',
        'AlwaysIDIC-frac' : 'AlwaysIDIC-frac',
        #'AlwaysIDIC-link' : 'AlwaysIDIC-link',
        'AlwaysISICType1-frac' : 'AlwaysType1-frac',
        'XYIntISICType1-frac' : 'XYIntType1-frac',
        'LIntISICType1-frac' : 'LIntType1-frac',
        'YIntISICType1-frac' : 'YIntType1-frac',
        
        ###### 1.2.1-final

        #'default' : 'Default (frac)',
        #'default' : 'Default',
        #'default-MIP' : 'Default (no MIP cuts)',
        #'default+linking' : 'Default (link)',
        #'default+ll' : 'Default (lower)',
        #'default+ll' : 'Branch on second level variables',
        #'noCut' : "No Cuts (link)", 
        #'fracISICType1-frac': 'Frac ISIC Type 1 (frac)',
        #'fracISICType1-frac-every': 'Frac ISIC Type 1 (frac)',
        #'ISICType1-frac': 'ISIC Type 1 (frac)', ##########
        #'ISICType1-frac-lv': 'ISIC Type 1 (frac-lv)',
        #'fracISICType1-link': 'Frac ISIC Type 1 (link)',
        #'ISICType1-link': 'ISIC Type 1 (link)',
        #'ISICType1-link-lv': 'ISIC Type 1 (link-lv)',
        #'fracISICType2-frac': 'Frac ISIC Type 2 (frac)',
        #'ISICType2-frac': 'ISIC Type 2 (frac)',
        #'fracISICType2-link': 'Frac ISIC Type 2 (link)',
        #'ISICType2-link': 'ISIC Type 2 (link)',  ##########
        #'fracIDIC-frac': 'Frac IDIC (frac)',     ##########
        #'fracIDIC+MIP-frac': 'Frac IDIC + MIP (frac)',
        #'IDIC-frac': 'IDIC (frac)',
        #'IDIC+MIP-frac': 'IDIC + MIP (frac)',
        #'fracIDIC-link': 'Frac IDIC (link)',
        #'IDIC-link': 'IDIC (link)',
        #'fracIDIC-ll': 'Frac IDIC (ll)',
        #'hyper-link': 'Hypercube IC (link)',
        #'hyper-frac': 'Hypercube IC (frac)',    ##########
        #'hyper-link-lv': 'Hypercube IC (link-lv)',
        #'hyper-frac-lv': 'Hypercube IC (frac-lv)',
        #'bendersInterdiction-frac': 'Benders Interdict (frac)',
        #'bendersInterdiction-frac': 'Branch on all variables',
        #'bendersInterdiction+MIP-frac': 'Benders Interdict + MIP (frac)',
        #'bendersInterdiction-link': 'Benders Interdict (link)',
        #'bendersInterdiction-link': 'Branch on linking variables',
        #'bendersInterdiction+MIP-link': 'Benders Interdict + MIP (link)',
        #'bendersBinary-frac': 'Benders Binary (frac)', ##########
        #'bendersBinary-link': 'Benders Binary (link)',
        #'bendersBinary-frac-lv': 'Benders Binary (frac-lv)',
        #'bendersBinary-link-lv': 'Benders Binary (link-lv)',
        #'genNoGood-frac': 'Generalized No Good (frac)', ##########
        #'genNoGood-link': 'Generalized No Good (link)',
        #'intNoGood-frac': 'Integer No Good (frac)',
        #'intNoGood+MIP-frac': 'Integer No Good + MIP (frac)',
        #'intNoGood-link': 'Integer No Good (link)', ##########
        #'fracIDIC+ISICType1-frac': 'fracIDIC + ISIC Type 1 (frac)',
        #'fracIDIC+ISICType1-frac-lv': 'fracIDIC + ISIC Type 1 (frac-lv)',
        #'fracIDIC+ISICType1-bendersBinary-frac': 'fracIDIC + ISIC Type 1 + bendersBin (frac)',
        #'fracIDIC+ISICType1-bendersBinary-frac-lv': 'fracIDIC + ISIC Type 1 +bendersBin (frac-lv)',

        ###### 1.2.1-opt

        # 'IDIC' : 'IDIC (link)',
        # 'ISICType1' : 'ISIC Type1 (link)',
        # 'ISICType2' : 'ISIC Type2 (link)',
        # 'ISICType2-frac' : '',
        # 'ISICType2-frac-old' : '',
        # 'ISICType2-old' : '',
        # 'benders+fracidic' : 'Benders + FracIDIC (link)',
        # 'benders+fracidic-frac' : 'Benders + Frac IDIC (frac)',
        # 'benders+idic' : 'Benders + IDIC (link)',
        # 'benders+idic-frac' : 'Benders + IDIC (frac)',
        # 'default+NoFractionalCutsAndFrac' : 'No Frac Cuts (Frac)',
        # 'default+NoFractionalCutsAndLinking' : 'No Frac Cuts (Linking)',
        # 'default+frac' : 'Default (frac)',
        # 'default+linking': 'Default (Linking)',
        # 'fracIDIC' : 'Frac IDIC (link)',
        # 'fracISICType2' : 'Frac ISIC Type2 (link)',
        # 'fracISICType2-frac' : 'Frac ISIC Type2 (frac)',
        # 'fracISICType2-frac-old' : 'Frac ISIC Type2 (frac)',
        # 'genNoGood' : 'GenNoGood (link)',
        # 'hyper' : 'Hypercube IC (link)',
        # 'intNoGood' : 'IntNoGood (link)',

        ###### 1.2.0-opt

        # 'IDIC-frac' : 'IDIC (frac)',    
        # 'IDIC+MIP-frac' : 'IDIC+MIP (frac)',
        # 'IDIC+MIP2-frac' : 'IDIC+MIP2 (frac)',
        # 'default+FracRoot' : 'Default+FracRoot',
        # 'default+MIP' : "Default + MIP",
        # 'default+NoFractionalCutsAndFrac' : 'No Frac Cuts (Frac)',
        # 'default+NoFractionalCutsAndLinking' : 'No Frac Cuts (Linking)',
        #'defaultWithExtraOutput' : 'Default',
        #'default+NoFractionalCutswithExtraOutput' : 'No Frac Cuts',
        # 'default+Parallel' : 'Default Parallel',
        # 'default+frac' : 'Default (frac)',
        # 'default+linking': 'Default (Linking)',
        # "fracIDIC-frac" : 'Frac IDIC (frac)',
        # "fracIDIC+MIP-frac" : 'FracIDIC+MIP (frac)',
        # "fracIDIC+MIP2-frac" : 'FracIDIC+MIP2 (frac)',

        ###### 1.2-opt

        #'benders' : 'Benders Interdict (link)',
        #'benders-cplex' : 'Benders Interdict CPLEX (link)',
        #'benders-frac' : 'Benders Interdict (frac)',
        #'benders-frac-cplex' : 'Benders Interdict CPLEX (frac)',
        #'bound+hyper-cplex' : 'Bound Cut + Hypercube IC',
        #'boundCut-cplex': 'Bound Cut + Hypercube IC CPLEX',
        #'fracWatermelon' : 'Frac1 IDIC (link)',
        #"fracWatermelon-frac" : 'Frac IDIC (frac)',
        #"fracWatermelon-frac-cplex" : 'Frac IDIC CPLEX (frac)',
        #"fracWatermelon-frac-cplex-d10" : 'Frac IDIC CPLEX d10 (frac)',
        #"fracWatermelon-frac-cplex-d20" : 'Frac IDIC CPLEX d20 (frac)',
        #"fracType1-frac-cplex" : 'Frac Type 1 CPLEX (frac)',
        #"fracType1-frac-cplex-d10" : 'Frac Type 1 CPLEX d10 (frac)',
        #"fracType1-frac-cplex-d20" : 'Frac Type 1 CPLEX d20 (frac)',
        #'genNoGood': 'Generalized No Good (link)',
        #'genNoGood-frac': 'Generalized No Good (frac)',
        #'hyper': 'Hypercube IC (link)',
        #'hyper-frac': 'Hypercube IC (frac)',
        #"incObj" : "BendersBinary (link)",
        #'incObj-frac' : "BendersBinary (frac)",
        #'noCut' : "No Cuts (link)", 
        #'noCut-WS' : "No Cuts (link+WS)", 
        #'noCut-cplex' : "No Cuts CPLEX (link)", 
        #'noGood+incObj' : 'GenNoGood+BendBin (link)',
        #'noGood+type1+pureInteger' : 'GenNoGood+Type1+IntNoGood (link)',
        #'pureInteger' : 'IntNoGood (link)',
        #'pureInteger-frac' : 'IntNoGood (frac)',
        #'strengthenedIntegerNoGoodCut-cplex': 'Strengthened Int No Good CPLEX',
        #"type1" : 'ISIC1 Type1 (link)',
        #"type1-cplex" : 'Type1 (link)',
        #"type1-frac" : 'Type1 (frac)',
        #'type1-WS' : 'Type1 (link+WS)',
        #"type2" : "ISIC1 Type2 (link)",
        #"type2-frac" : "ISIC1 Type2 (frac)",
        #'watermelon' : 'Watermelon (link)',
        #'watermelon-frac' : 'Watermelon (frac)',    
        #"watermelon-frac-LV" : 'Watermelon (frac+LV)',
        #'watermelon+type1-frac' : 'Watermelon+Type1 (frac)',
        #'watermelon+type1-frac-LV' : 'Watermelon+Type1 (frac+LV)',
        #'watermelon+type1+incObj-frac' : "Watermelon+Type1+BendBin (frac)",
        #'watermelon+type1+incObj-frac-LV' : 'Watermelon+Type1+BendBin (frac+LV)',

        ###### 1.2-cplex-opt

        #'default' : 'Default',
        #'default+ParallelCplex': 'Default Parallel CPLEX',

        ###### tailoff

        #'default-new-tailoff-05' : 'Default',
        #'default+frac-new-tailoff' : 'Default (frac+new-tailoff-01)',
        #'default+link-new-tailoff' : 'Default (link+new-tailoff-01)',
        #'default+frac-new-tailoff-1' : 'Default (frac+new-tailoff-1)',
        #'default+link-new-tailoff-1' : 'Default (link+new-tailoff-1)',
        #'default+frac-new-tailoff-05' : 'Default (frac+new-tailoff-05)',
        #'default+link-new-tailoff-05' : 'Default (link+new-tailoff-05)',
        #'default+frac-new-tailoff-005' : 'Default (frac+new-tailoff-005)',
        #'default+link-new-tailoff-005' : 'Default (link+new-tailoff-005)',
        #'default+link-new-tailoff-001' : 'Default (link+new-tailoff-001)',

        ###### 1.0-opt

        #'default' : 'Default',
        
        ###### 1.1-opt

        #'default' : 'Default',
        
        ###### DenRal09

        #'default' : 'Default',

        ###### rev1

        #'fracWatermelon' : 'Frac1 IDIC (link)',
        #"fracWatermelon-frac" : 'Frac IDIC (frac)',
        #'genNoGood': 'Generalized No Good (link)',
        #'genNoGood-frac': 'Generalized No Good (frac)',
        #'hyper': 'Hypercube IC (link)',
        #'hyper-frac': 'Hypercube IC (frac)',
        #"incObj" : "BendersBinary (link)",
        #'incObj-frac' : "BendersBinary (frac)",
        #'noCut' : "No Cuts (link)", 
        #'pureInteger' : 'IntNoGood (link)',
        #'pureInteger-frac' : 'IntNoGood (frac)',
        #"type1" : 'ISIC1 Type1 (link)',
        #"type1-frac" : 'Type1 (frac)',
        #"type2" : "ISIC1 Type2 (link)",
        #"type2-frac" : "ISIC1 Type2 (frac)",
        #'watermelon' : 'Watermelon (link)',
        #'watermelon+type1-frac' : 'Watermelon+Type1 (frac)',
        #'watermelon+type1-frac-LV' : 'Watermelon+Type1 (frac+LV)',
        #'watermelon-frac' : 'Watermelon (frac)',    
        #"watermelon-frac-LV" : 'Watermelon (frac+LV)',
       
        #'fracISICType1-frac-everyIteration': 'Frac ISIC Type 1 (frac-every)',
        #'fracWatermelon-frac': 'Default(frac)',
        #'default+frac-.1' : 'Default (frac, 0.1)',
        #'default+frac-.5' : 'Default (frac, 0.5)',
        #'default+frac-1' : 'Default (frac, 1)',
        #'default+NoFractionalCuts' : 'No Frac Cuts (frac)',
        #'default+NoFractionalCuts+Frac' : 'No Frac Cuts (frac)',
        #'default+NoFractionalCuts+Frac-.1' : 'No Frac Cuts (frac, 0.1)',
        #'default+Linking': 'Default (linking)',
        #'default+Linking-.1': 'Default (linking, 0.1)',
        #'default+Linking-.5': 'Default (linking, 0.5)',
        #'default+Linking-1': 'Default (linking, 1)',
        #'default+NoFractionalCuts+Linking' : 'No Frac Cuts (linking)',
        #'default+WS' : 'Default w/ Warm Start',
        # 'default-frac',
        #'benders-frac' : 'BendersInterdict (frac)',
        #'genNoGood-frac' : 'GenNoGood (frac)',
        #'hyper-frac' : 'Hypercube (frac)',
        #'IDIC' : 'Watermelon (link)',
        #'fracidic+ll' : 'Frac IDIC (ll)',
        #"ISICType2" : "ISIC Type2 (link)",
        #"ISICType2-frac" : "ISIC Type2 (frac)",
        # 'interdiction',
    }
    ################# Process & Save | Load from CSV ###################
    # specify summary file name
    file_csv_out = "summary_"+name+".csv"
    #file_csv_in = "summary-1.2.1.csv"
    file_csv_in = "summary_branching.csv"
    
    # if len(args) == 0:
    if 1:
        df_r = parseOutput(
            outputDir, versions, scenarios, writeCSV=True, filename=file_csv_out,
            name=name
        )
    else:
        try:
            df_r = pd.read_csv(file_csv_in)
            set_cond = (df_r["scenario"].isin(scenarios.values())) | (
                df_r["dataset"].isin(dataSets)
            )
            df_r = df_r[set_cond]
        except FileNotFoundError:
            print("{} does not exist in current directory.".format(file_csv))
        else:
            print("Reading from", file_csv_in)

    ################### Format Data & Print Table ####################
    # specify txt file name to print tables in LATEX
    file_txt = "ltx_tb_cut.txt"

    # columns to process and print
    displayCols = {
        "cpu": "CPU Search Time",
        "nodes": "Number of Processed Nodes",
        "gap": "Final Gap",
        "root_gap": "Root Gap",
        "100_gap": "Gap After 100 Nodes",
        "solved": "Solved",
        "num_int_idic": "Int IDIC Success",
        "num_int_lv_idic": "Int IDIC (LV) Success",
        "num_frac_idic": "Frac IDIC Success",
        "num_int_lv_idic_fail": "Int IDIC (LV) Fail",
        "num_frac_idic_fail": "Frac IDIC Fail",
        "cg_called": "CG Calls",
        "cg_failed": "CG Failures",
        "cg_fail_rate": "CG Failure Rate",
        "num_cuts": "Number of Cuts",
        "depth_idic" : "Average Depth of IDIC",
        "num_idic" : "Number of IDIC",
        #'chk_feas_time': 'Check Feasibility Time',
        #'vf_solved': 'Number of VF problem solved',
        #'ub_solved': 'Number of UB problem solved',
        'objval': 'Object Value'
    }

    df_proc = processTable(df_r, displayCols, writeLTX=False, filename=file_txt)

    ################### Make Performance Profile ####################
    # columns to compare in the plot
    plotCols = {
        "cpu": ["CPU Time", 25],
        "nodes": ["Nodes Processed", 50],
        "root_gap": ["Root Gap", 10],
        #"num_cuts": ["Number of Cuts", 50],
        #"depth_idic" : ["Depth of IDICs", 10],
        #"num_idic" : ["Number of IDICs", 10],
        #'cg_called': ["CG Calls", 10],
        #'cg_failed': ["CG Failures", 10],
        #'cg_fail_rate': ['CG Failre Rate', 10],
        #"100_gap": ["Gap After 100 Nodes", 10]
    }
    # plotCols = {}

    # manual input example:
    # for k in scenarios:
    #     if '01' in k:
    #         scenarios[k] = 'linkingBranchStrategy'
    #     else:
    #         scenarios[k] = 'fractionalBranchStrategy'

    baseline=None
    #baseline = ('IDIC-frac', '1.2.1-opt')
    #baseline = ('default', '1.2.1-opt')
    #baseline = ("Type1IC", "1.2-opt")
    #baseline = ('GenNoGood+Type1+IntNoGood (link)', '1.2-opt')
    #baseline = ('Watermelon (frac+LV)', '1.2-opt')
    #baseline = ('FracWatermelon (frac)', '1.2-opt')
    #baseline = ('Benders Interdict (link)', '1.2.1-final')
    if len(versions) > 1:
        versionlegend = True
    else:
        versionlegend = False

    print(df_proc)

    if name != '':
        dataSets = [name]
        
    for ds in dataSets:
        df_solved, df_has_soln = dropFilter(df_proc, scenarios, ds)
        for col in plotCols:
            if col != "root_gap":
                df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
            else:
                df_sub = df_has_soln.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
            if versionlegend:
                plottitle = ": "+plotCols[col][0]+" ("+ds+")"
            else:
                plottitle = " (v"+list(versions.values())[0]+"): "+plotCols[col][0]+" ("+ds+")"

            print("")
            print("Creating performance profile for "+col)
            print("")
            plotPerfProf(
                df_sub, versions, plotname="perf_" + col + "_" + ds,
                plottitle = "Performance Profile"+plottitle,
                xmin = 0.0, xmax=plotCols[col][1],
                versionlegend = versionlegend
            )
            if baseline is not None: 
                print("")
                print("Creating baseline profile for "+col)
                print("")
                try:
                    if col != 'root_gap':
                        plotBaselineProf(
                            df_sub, versions,
                            baseline = (scenarios[baseline[0]],baseline[1]),
                            plotname="base_"+baseline[0]+"_"+col+"_"+ds,
                            plottitle = "Baseline Profile"+plottitle,
                            xmax=plotCols[col][1],
                            versionlegend = versionlegend
                        )
                    else:
                        plotBaselineProfSingle(
                            df_sub, versions,
                            baseline = (scenarios[baseline[0]],baseline[1]),
                            plotname="base_"+baseline[0]+"_"+col+"_"+ds,
                            plottitle = "Baseline Profile"+plottitle,
                            xmax=plotCols[col][1],
                            versionlegend = versionlegend
                        )
                except KeyError:
                    pass
            if col == 'cpu':
                plotCumProf(df_has_soln, versions, plotname="cum_" + col + "_" + ds,
                            plottitle="Cumulative Profile"+plottitle,
                            versionlegend = versionlegend
)
        # if baseline is not None: 
        #     print("")
        #     print("Creating baseline profile for gap")
        #     print("")
        #     df_gap = df_has_soln.xs(
        #         (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
        #     ).copy()
        #     df_baseline_has_gap = df_gap.drop(df_gap[df_gap[baseline] == 0].index.to_list())
            # plotBaselineProf(
            #     df_baseline_has_gap, baseline = baseline,
            #     plotname="base_" + baseline[0] + "_" + "gap_" + ds,
            #     plottitle = "Baseline Profile: Gap ("+ds+")",
            #     xmax=25, versionlegend = versionlegend
            # )
