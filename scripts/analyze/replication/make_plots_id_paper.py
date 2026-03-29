# Script to replicate experiments of paper 
# "Improving Directions in Mixed Integer Bilevel Linear Optimization"
# Battista F. & Ted K. Ralphs
# Last edited 2026

import os
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from performance_plots import *

def make_plots_improving_directions(params, keywords):
    """
        Make plots for improving directions paper BatRalImproving25.
    """

    versions = {
        'idbc':'idBC',
    }

    scenarios = {
        #### BatRalImproving25 (ID paper)
        'idBC-LS-k_2-dBnd_Inf_fracB'   : "idB&C-LS-k_2 (frac)",
        'idBC-LS-k_3-dBnd_Inf_fracB'   : "idB&C-LS-k_3 (frac)",
        # 'idBC-LS-k_2-dBnd_0_10_fracB'  : 'idB&C-LS-k_2-dBnd_0_10 (frac)',
        'idBC-LS-k_2-dBnd_10_Inf_fracB': 'idB&C-LS-k_2-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_0_10_fracB'  : 'idB&C-LS-k_3-dBnd_0_10 (frac)',
        'idBC-LS-k_3-dBnd_10_Inf_fracB': 'idB&C-LS-k_3-dBnd_10_Inf (frac)',
        # 'idBC-LS-k_3-dBnd_0_8_fracB'   : 'idB&C-LS-k_3-dBnd_0_8 (frac)',
        'idBC-LS-k_3-dBnd_8_Inf_fracB' : 'idB&C-LS-k_3-dBnd_8_Inf (frac)',
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
        "MibS_1_2-LS-k_2_defaultB"     : "MibS IDIC-LS-k_2 (default)",
    }

    displayCols = {
        'cpu': 'CPU Search Time',
        'nodes': 'Number of Processed Nodes',
        'gap': 'Final Gap',
        'root_gap': 'Root Gap',
        '100_gap': 'Gap After 100 Nodes',
        'solved': 'Solved',
        'vf_solved': 'Number of VF problem solved',
        'ub_solved': 'Number of UB problem solved',
        'objval': 'Object Value',
        'cg_called': 'CG Calls',
        'num_cuts': 'Number of Cuts',
        'cut_time': 'Cut Generation Time',
        'time_per_cg_call': 'Time Per CG Call',
        'idic_ifd_time': "Finding IFDs CPU Time",
        'idic_ifd_avg_time': "Finding IFDs Average CPU Time",
    }

    versionlegend = False

    outputDir = [params['outputDir']] 

    plotDir = os.path.join(params['outputDir'], "figures")

    os.makedirs(plotDir, exist_ok=True)

    columns = ['instance','version','scenario','dataset']
    columns.extend(displayCols.keys())

    ################################################################################################
    #            IBLP DATASET                                                                      #
    ################################################################################################
    aggrDataSName = "F-D-D2-Z-Z2"

    dataSets = [
        "iblpDen",
        "iblpDen2",
        "iblpZhang",
        "iblpZhang2",
        'iblpFis'
    ]

    print("Dataset: "+aggrDataSName)

    file_csv_in = os.path.join(params['outputDir'], "summary_" + aggrDataSName + ".csv") 
    file_csv_out = os.path.join(params['outputDir'], "summary_" + aggrDataSName)

    try:
        df_r = pd.read_csv(file_csv_in)
    except FileNotFoundError:
        print("{} does not exist in current directory.".format(file_csv_in))
        df_r = parseOutput(outputDir, versions, scenarios, keywords, dataSets, aggrDataSName)
        export(df_r, columns, file_csv_out)
    else:
        print("Reading from", file_csv_in)
        set_cond = df_r["scenario"].isin(scenarios.values())
        try: 
            set_cond |= df_r["dataset"].isin(dataSets)
        except:
            pass
        df_r = df_r[set_cond]

    df_proc = processTable(df_r, displayCols, writeLTX=False)

    # =============================================================================================
    #   Figures 8 and 9
    # =============================================================================================

    plotCols = {
        "cpu": ["CPU Time", 40],
        "nodes": ["Nodes Processed", 15],
        'idic_ifd_time': ["Finding IFDs total CPU Time", 50],
        'idic_ifd_avg_time': ["Finding IFDs average CPU Time", 60]
    }

    scenariosToPlot = {
        'idBC-LS-k_2-dBnd_Inf_fracB'   : "idB&C-LS-k_2 (frac)",
        "idBC-LS-k_3-dBnd_Inf_fracB"   : "idB&C-LS-k_3 (frac)",
        'idBC-LS-k_2-dBnd_10_Inf_fracB': 'idB&C-LS-k_2-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_10_Inf_fracB': 'idB&C-LS-k_3-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_8_Inf_fracB' : 'idB&C-LS-k_3-dBnd_8_Inf (frac)',
        'idBC-LS-k_3-dBnd_12_Inf_fracB': 'idB&C-LS-k_3-dBnd_12_Inf (frac)',
        'idBC-LS-k_4-dBnd_10_Inf_fracB': 'idB&C-LS-k_4-dBnd_10_Inf (frac)',
        "idBC-MILP_fracB"              : "idB&C-MILP (frac)",
        'idBC-k_2-MILP_fracB'          : 'idB&C-MILP-k_2 (frac)',
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
    }

    baseline = ("MibS_onlyIDIC_fracB", "idbc")

    figuresDir = os.path.join(plotDir, "figures_8")
    os.makedirs(figuresDir, exist_ok=True)

    print("   Figure 8")

    dataSets = [aggrDataSName]
    ds = dataSets[0]

    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenariosToPlot, ds)

    # CPU Time
    # =============================================================================================
    col = 'cpu'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_8_a_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_8_b_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # Nodes
    # =============================================================================================
    col = 'nodes'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()

    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_8_c_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_8_d_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    figuresDir = os.path.join(plotDir, "figures_9")
    os.makedirs(figuresDir, exist_ok=True)

    print("   Figure 9")

    # IFD Time
    # =============================================================================================
    col = 'idic_ifd_time'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()

    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_9_a_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_9_b_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # Avg IFD Time
    # =============================================================================================
    col = 'idic_ifd_avg_time'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()

    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_9_c_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_9_d_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # =============================================================================================
    #   Figures 12 (a)---(d)
    # =============================================================================================
    plotCols = {
        "cpu": ["CPU Time", 40]
    }
    
    scenariosToPlot = {
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
        "MibS_1_2_defaultB"            : "MibS (default)",
        "MibS_1_2_IDIC_defaultB"       : "MibS IDIC-MILP (default)",
        "MibS_1_2-LS-k_2_defaultB"     : "MibS IDIC-LS-k_2 (default)"
    }

    baseline = ("MibS_1_2_defaultB", "idbc")

    figuresDir = os.path.join(plotDir, "figures_12")
    os.makedirs(figuresDir, exist_ok=True)

    print("   Figure 12")

    dataSets = [aggrDataSName]
    ds = dataSets[0]

    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenariosToPlot, ds)

    # CPU Time
    # =============================================================================================
    col = 'cpu'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    df_time = df_has_soln.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    df_gap = df_has_soln.xs(
                    (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_12_a_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_12_b_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_12_c_base_gap_"+ds+".pdf")

    plotBaselineProf(
        df_gap, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile: Gap ("+ ds +")",
        xmax=25,
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_12_d_cum_"+col+"_"+ds+".pdf")

    plotCumProf(
        df_time, df_gap, versions, 
        plotname=plotfilepath,
        plottitle="Cumulative Profile: Time-Gap (" + ds + ")",
        versionlegend = versionlegend
    )

    # =============================================================================================
    #   Figures 13 (a) and (b)
    # =============================================================================================
    plotCols = {
        "cpu": ["CPU Time", 40]
    }
    
    scenariosToPlot = {
        'idBC-LS-k_2-dBnd_Inf_fracB'   : "idB&C-LS-k_2 (frac)",
        "idBC-LS-k_3-dBnd_Inf_fracB"   : "idB&C-LS-k_3 (frac)",
        'idBC-LS-k_2-dBnd_10_Inf_fracB': 'idB&C-LS-k_2-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_0_10_fracB'  : 'idB&C-LS-k_3-dBnd_0_10 (frac)',
        'idBC-LS-k_3-dBnd_10_Inf_fracB': 'idB&C-LS-k_3-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_12_Inf_fracB': 'idB&C-LS-k_3-dBnd_12_Inf (frac)',
        'idBC-LS-k_4-dBnd_10_Inf_fracB': 'idB&C-LS-k_4-dBnd_10_Inf (frac)',
        'idBC-LS-k_5-dBnd_10_Inf_fracB': 'idB&C-LS-k_5-dBnd_10_Inf (frac)',
        "idBC-MILP_fracB"              : "idB&C-MILP (frac)",
        'idBC-k_2-MILP_fracB'          : 'idB&C-MILP-k_2 (frac)',
        'idBC-k_4-MILP_fracB'          : 'idB&C-MILP-k_4 (frac)',
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
        "MibS_1_2_defaultB"            : "MibS (default)",
        "MibS_1_2_IDIC_defaultB"       : "MibS IDIC-MILP (default)",
        "MibS_1_2-LS-k_2_defaultB"     : "MibS IDIC-LS-k_2 (default)"
    }

    baseline = ("MibS_1_2_defaultB", "idbc")

    figuresDir = os.path.join(plotDir, "figures_13")
    os.makedirs(figuresDir, exist_ok=True)

    print("   Figure 13")

    dataSets = [aggrDataSName]
    ds = dataSets[0]

    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenariosToPlot, ds)

    # CPU Time
    # =============================================================================================
    col = 'cpu'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_13_a_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_13_b_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    ################################################################################################
    #            INTERDICTION DATASET                                                              #
    ################################################################################################
    aggrDataSName = "INT-DEN-SHI"

    dataSets = [
        "interKpShi",
        "interDen"
    ]

    print("Dataset: "+aggrDataSName)

    file_csv_in = os.path.join(params['outputDir'], "summary_" + aggrDataSName + ".csv") 
    file_csv_out = os.path.join(params['outputDir'], "summary_" + aggrDataSName)

    try:
        df_r = pd.read_csv(file_csv_in)
    except FileNotFoundError:
        print("{} does not exist in current directory.".format(file_csv_in))
        df_r = parseOutput(outputDir, versions, scenarios, keywords, dataSets, aggrDataSName)
        export(df_r, columns, file_csv_out)
    else:
        print("Reading from", file_csv_in)
        set_cond = df_r["scenario"].isin(scenarios.values())
        try: 
            set_cond |= df_r["dataset"].isin(dataSets)
        except:
            pass
        df_r = df_r[set_cond]

    df_proc = processTable(df_r, displayCols, writeLTX=False)

    # ===============================
    #   Figures 10 and 11
    # ===============================
    plotCols = {
        "cpu": ["CPU Time", 20],
        "nodes": ["Nodes Processed", 15],
        'idic_ifd_time': ["Finding IFDs total CPU Time", 50],
        'idic_ifd_avg_time': ["Finding IFDs average CPU Time", 60]
    }

    scenariosToPlot = {
        'idBC-LS-k_2-dBnd_Inf_fracB'   : "idB&C-LS-k_2 (frac)",
        "idBC-LS-k_3-dBnd_Inf_fracB"   : "idB&C-LS-k_3 (frac)",
        'idBC-LS-k_2-dBnd_10_Inf_fracB': 'idB&C-LS-k_2-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_10_Inf_fracB': 'idB&C-LS-k_3-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_8_Inf_fracB' : 'idB&C-LS-k_3-dBnd_8_Inf (frac)',
        'idBC-LS-k_3-dBnd_12_Inf_fracB': 'idB&C-LS-k_3-dBnd_12_Inf (frac)',
        'idBC-LS-k_4-dBnd_10_Inf_fracB': 'idB&C-LS-k_4-dBnd_10_Inf (frac)',
        "idBC-MILP_fracB"              : "idB&C-MILP (frac)",
        'idBC-k_2-MILP_fracB'          : 'idB&C-MILP-k_2 (frac)',
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
    }

    baseline = ("MibS_onlyIDIC_fracB", "idbc")

    figuresDir = os.path.join(plotDir, "figures_10")
    os.makedirs(figuresDir, exist_ok=True)

    print("   Figure 10")

    dataSets = [aggrDataSName]
    ds = dataSets[0]

    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenariosToPlot, ds)

    # CPU Time
    # =============================================================================================
    col = 'cpu'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_10_a_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_10_b_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # Nodes
    # =============================================================================================
    col = 'nodes'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()

    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_10_c_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_10_d_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    figuresDir = os.path.join(plotDir, "figures_11")
    os.makedirs(figuresDir, exist_ok=True)

    print("   Figure 11")

    # IFD Time
    # =============================================================================================
    col = 'idic_ifd_time'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()

    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_11_a_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_11_b_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # Avg IFD Time
    # =============================================================================================
    col = 'idic_ifd_avg_time'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()

    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_11_c_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_11_d_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # ===============================
    #   Figures 12 (e) and (f)
    # ===============================
    plotCols = {
        "cpu": ["CPU Time", 40]
    }
    
    scenariosToPlot = {
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
        "MibS_1_2_defaultB"            : "MibS (default)",
        "MibS_1_2_IDIC_defaultB"       : "MibS IDIC-MILP (default)",
        "MibS_1_2-LS-k_2_defaultB"     : "MibS IDIC-LS-k_2 (default)"
    }

    baseline = ("MibS_1_2_defaultB", "idbc")

    figuresDir = os.path.join(plotDir, "figures_12")

    print("   Figure 12")

    dataSets = [aggrDataSName]
    ds = dataSets[0]

    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenariosToPlot, ds)

    # CPU Time
    # =============================================================================================
    col = 'cpu'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_12_e_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_12_f_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    # ===============================
    #   Figures 13 (c) and (d)
    # ===============================
    plotCols = {
        "cpu": ["CPU Time", 40],
    }
    
    scenariosToPlot = {
        'idBC-LS-k_2-dBnd_Inf_fracB'   : "idB&C-LS-k_2 (frac)",
        "idBC-LS-k_3-dBnd_Inf_fracB"   : "idB&C-LS-k_3 (frac)",
        'idBC-LS-k_2-dBnd_10_Inf_fracB': 'idB&C-LS-k_2-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_0_10_fracB'  : 'idB&C-LS-k_3-dBnd_0_10 (frac)',
        'idBC-LS-k_3-dBnd_10_Inf_fracB': 'idB&C-LS-k_3-dBnd_10_Inf (frac)',
        'idBC-LS-k_3-dBnd_12_Inf_fracB': 'idB&C-LS-k_3-dBnd_12_Inf (frac)',
        'idBC-LS-k_4-dBnd_10_Inf_fracB': 'idB&C-LS-k_4-dBnd_10_Inf (frac)',
        'idBC-LS-k_5-dBnd_10_Inf_fracB': 'idB&C-LS-k_5-dBnd_10_Inf (frac)',
        "idBC-MILP_fracB"              : "idB&C-MILP (frac)",
        'idBC-k_2-MILP_fracB'          : 'idB&C-MILP-k_2 (frac)',
        'idBC-k_4-MILP_fracB'          : 'idB&C-MILP-k_4 (frac)',
        "MibS_onlyIDIC_fracB"          : "MibS only IDICs (frac)",
        "MibS_1_2_defaultB"            : "MibS (default)",
        "MibS_1_2_IDIC_defaultB"       : "MibS IDIC-MILP (default)",
        "MibS_1_2-LS-k_2_defaultB"     : "MibS IDIC-LS-k_2 (default)"
    }

    baseline = ("MibS_1_2_defaultB", "idbc")

    figuresDir = os.path.join(plotDir, "figures_13")

    print("   Figure 13")

    dataSets = [aggrDataSName]
    ds = dataSets[0]

    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenariosToPlot, ds)

    # CPU Time
    # =============================================================================================
    col = 'cpu'

    print("      Creating profiles for "+col)

    df_sub = df_solved.xs(
                    (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
                ).copy()
    
    plottitle = ": "+plotCols[col][0]+" ("+ds+")"
    plotfilepath = os.path.join(figuresDir, "figure_13_c_base_"+col+"_"+ds+".pdf")
    
    plotBaselineProf(
        df_sub, versions,
        baseline = (scenarios[baseline[0]],baseline[1]),
        plotname=plotfilepath,
        plottitle = "Baseline Profile"+plottitle,
        xmax=plotCols[col][1],
        versionlegend = versionlegend
    )

    plotfilepath = os.path.join(figuresDir, "figure_13_d_perf_"+col+"_"+ds+".pdf")

    plotPerfProf(
        df_sub, versions, 
        plotname=plotfilepath,
        plottitle = "Performance Profile"+plottitle,
        xmin = 0.0, xmax=plotCols[col][1],
        versionlegend = versionlegend
    )


    

    
