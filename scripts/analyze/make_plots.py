# Script to parse MibS results files.
# Script path:  /MibS/scripts/analyze
# Some os function requires Python 3.5+

import sys
import os
import collections
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from performance_plots import *

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
    "num_full_int_idic" : "Full Int IDIC",
    "num_link_int_idic" : "Linking Int IDIC",
    "num_lower_int_idic" : "Second-level Int IDIC",
    "num_frac_idic" : "Fractional IDIC",
    "num_full_int_idic_fail" : "Full Int IDIC (Fail)",
    "num_link_int_idic_fail" : "Linking Int IDIC (Fail)",
    "num_lower_int_idic_fail" : "Second-level Int IDIC (Fail)",
    "num_frac_idic_fail" : "Fractional IDIC (Fail)",
    "num_full_int_isic" : "Full Int ISIC",
    "num_link_int_isic" : "Linking Int ISIC",
    "num_lower_int_isic" : "Second-level Int ISIC",
    "num_frac_isic" : "Fractional",
    "num_full_int_isic_fail" : "Full Int ISIC (Fail)",
    "num_link_int_isic_fail" : "Linking Int ISIC (Fail)",
    "num_lower_int_isic_fail" : "Second-level Int ISIC (Fail)",
    "num_frac_isic_fail" : "Fractional ISIC (Fail)",
}

dataSets = [
    #'MIBLP-XU',
    'IBLP-FIS',
    #'INTERD-DEN',
    'IBLP-DEN',
    'IBLP-DEN2',
    'IBLP-ZHANG',
    'IBLP-ZHANG2',
    #'BENCHMARK'
    #'all'
]

aggregate = True
#aggregate = False
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
    'filmosi':'filmosi',
    #'1.0.0-opt':'1.0.0',
    #'1.1.3-opt':'1.1.3',
    #'1.2.0-opt':'1.2.0',
    #'1.2.1-opt':'1.2.1',
    #'1.2.1-cplex-opt':'1.2.1-cplex',
    '1.2.2-opt':'1.2.2',
    #'1.2.2-opt-5':'1.2.2-5',
    '1.2.2-cplex-opt':'1.2.2-cplex',
    #'1.2.2-opt-cg-fail':'1.2.2',
    #'old':'old'
}
    
# Output parent path
outputDir = ["/home/ted/Projects/MibS/output"]
#outputDir = ["/home/ted/Projects/MibS/output-Mac"]

scenarios = {
    
    #### Interdiction
    # 'default' : 'Default',
    # 'noCut' : "No Cuts (link)",
    # 'AlwaysIDIC-link': 'AlwaysIDIC (frac)',    
    # 'LIntISICType1-link': 'LIntType 1 (frac)', 
    # 'AlwaysISICType2-link': 'ISIC Type 2 (link)',
    # 'bendersInterdiction-frac': 'Benders Interdict (frac)',
    # 'bendersInterdiction-link': 'Benders Interdict (link)',
    
    #### Pure Integer
    # 'default' : 'Default',
    # 'noCut' : "No Cuts (link)",
    # 'AlwaysIDIC-frac': 'AlwaysIDIC (frac)',    
    # 'LIntISICType1-frac': 'LIntType 1 (frac)', 
    # 'XYIntISICType2-frac': 'XYIntType 2 (frac)',
    # 'hyper-frac': 'Hypercube IC (frac)',    
    # 'intNoGood-frac': 'Integer No Good (frac)',
    
    #### Binary First-Level
    # 'default' : 'Default',
    # 'noCut' : "No Cuts (link)", 
    # 'AlwaysIDIC-frac': 'AlwaysIDIC (frac)',
    # 'YLIntISICType1-frac': 'YLInt Type 1 (frac)', 
    # 'XYIntISICType2-frac': 'XYInt Type 2 (frac)',
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
    
    'default' : 'Default',
    'sep++' : 'SEP++',
    #'default-MIP' : 'Default (no MIP cuts)',
    # 'default+linking' : 'Default (Link)',
    # 'default+ll': 'Default (Lower Level)',
    #'AlwaysIDIC-frac' : 'AlwaysIDIC-frac',
    # 'XYIntIDIC-frac' : 'XYIntIDIC-frac',
    # 'LIntIDIC-frac' : 'LIntIDIC-frac',
    # 'YIntIDIC-frac' : 'YIntIDIC-frac',
    # 'YLIntIDIC-frac' : 'YLIntIDIC-frac',
    # 'AlwaysIDIC-link' : 'AlwaysIDIC-link',
    # 'XYIntIDIC-link' : 'XYIntIDIC-link',
    # 'LIntIDIC-link' : 'LIntIDIC-link',
    # 'YIntIDIC-link' : 'YIntIDIC-link',
    # 'AlwaysISICType1-frac' : 'AlwaysType1-frac',
    # 'XYIntISICType1-frac' : 'XYIntType1-frac',
    # 'LIntISICType1-frac' : 'LIntType1-frac', 
    # 'YIntISICType1-frac' : 'YIntType1-frac',
    # 'YLIntISICType1-frac' : 'YLIntType1-frac',
    # 'AlwaysISICType1-link' : 'AlwaysType1-link',
    # 'XYIntISICType1-link' : 'XYIntType1-link',
    # 'LIntISICType1-link' : 'LIntType1-link',
    # 'YIntISICType1-link' : 'YIntType1-link',
    # 'YLIntISICType1-link' : 'YLIntType1-link',
    # 'AlwaysISICType2-frac' : 'AlwaysType2-frac', ########
    # 'XYIntISICType2-frac' : 'XYIntType2-frac',
    # 'LIntISICType2-frac' : 'LIntType2-frac',
    # 'YIntISICType2-frac' : 'YIntType2-frac',
    # 'YLIntISICType2-frac' : 'YLIntType2-frac',
    # 'AlwaysISICType2-link' : 'AlwaysType2-link',
    # 'XYIntISICType2-link' : 'XYIntType2-link',
    # 'LIntISICType2-link' : 'LIntType2-link',
    # 'YIntISICType2-link' : 'YIntType2-link',
    # 'YLIntISICType2-link' : 'YLIntType2-link',
    # 'hyper-link': 'Hypercube IC (link)',
    # 'hyper-frac': 'Hypercube IC (frac)',    ##########
    # 'bendersBinary-frac': 'Benders Binary (frac)', ##########
    # 'bendersBinary-link': 'Benders Binary (link)',
    # 'genNoGood-frac': 'Generalized No Good (frac)', ##########
    # 'genNoGood-link': 'Generalized No Good (link)',
    #'intNoGood-frac': 'Integer No Good (frac)',
    #'intNoGood-link': 'Integer No Good (link)', ##########
    # 'bendersInterdiction-frac': 'Benders Interdict (frac)',
    # 'bendersInterdiction-link': 'Benders Interdict (link)',
    # 'bendersInterdiction+AlwaysIDIC-link': 'Benders Interdict + Always IDIC (link)',
    # 'bendersInterdiction+AlwaysIDIC+genNoGood-link': 'Benders Interdict + Always IDIC + genNoGood (link)',

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

################### Format Data & Print Table ####################
# specify txt file name to print tables in LATEX
file_txt = "ltx_tb_cut.txt"

# columns to process and print
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

    ######## IDIC ###################

    'num_idic': 'Number of IDICs',
    'num_full_int_idic': 'Number of Full Int IDICs',
    'num_link_int_idic': 'Number of Link Int IDICs',
    'num_lower_int_idic': 'Number of Lower Int IDICs',
    'num_frac_idic': 'Number of Fractional IDICs',
    'num_idic_fail': 'Number of IDIC Fails',
    'num_full_int_idic_fail': 'Number of Full Int IDIC Fails',
    'num_link_int_idic_fail': 'Number of Link Int IDIC Fails',
    'num_lower_int_idic_fail': 'Number of Lower Int IDIC Fails',
    'num_frac_idic_fail': 'Number of Fractional IDIC Fails',
    'idic_fail_rate': 'IDIC Failure Rate',
    'full_int_idic_fail_rate': 'Full Int IDIC Fail Rate',
    'link_int_idic_fail_rate': 'Link Int IDIC Fail Rate',
    'lower_int_idic_fail_rate': 'Lower Int IDIC Fail Rate',
    'frac_idic_fail_rate': 'Fractional IDIC Fail Rate',

    ######## ISIC ###################

    # 'num_isic': 'Number of ISICs',
    # 'num_full_int_isic': 'Number of Full Int ISICs',
    # 'num_link_int_isic': 'Number of Link Int ISICs',
    # 'num_lower_int_isic': 'Number of Lower Int ISICs',
    # 'num_frac_isic': 'Number of Fractional ISICs',
    # 'num_isic_fail': 'Number of ISIC Fails',
    # 'num_full_int_isic_fail': 'Number of Full Int ISIC Fails',
    # 'num_link_int_isic_fail': 'Number of Link Int ISIC Fails',
    # 'num_lower_int_isic_fail': 'Number of Lower Int ISIC Fails',
    # 'num_frac_isic_fail': 'Number of Fractional ISIC Fails',
    # 'isic_fail_rate': 'ISIC Failure Rate',
    # 'full_int_isic_fail_rate': 'Full Int ISIC Fail Rate',
    # 'link_int_isic_fail_rate': 'Link Int ISIC Fail Rate',
    # 'lower_int_isic_fail_rate': 'Lower Int ISIC Fail Rate',
    # 'frac_isic_fail_rate': 'Fractional ISIC Fail Rate',
}

if '1.0.0-opt' not in versions and 'filmosi' not in versions:
    displayCols['chk_feas_time'] = 'Check Feasibility Time'


columns = ['instance','scenario']
if len(dataSets) > 1 and name == '':
    columns.extend(['dataset'])
if len(versions) > 1:
    columns.extend(['version'])
columns.extend(displayCols.keys())

if 1:
    df_r = parseOutput(outputDir, versions, scenarios, keywords, dataSets,
                       writeCSV=True, columns=columns,
                       filename=file_csv_out, name=name)
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

baseline=None
#baseline = ('IDIC-frac', '1.2.1-opt')
#baseline = ('default', '1.2.1-opt')
#baseline = ("Type1IC", "1.2-opt")
#baseline = ('GenNoGood+Type1+IntNoGood (link)', '1.2-opt')
#baseline = ('Watermelon (frac+LV)', '1.2-opt')
#baseline = ('FracWatermelon (frac)', '1.2-opt')
#baseline = ('Benders Interdict (link)', '1.2.1-opt')
#baseline = ('LIntISICType1-frac','1.2.2-opt-4')
#baseline = ('noCut', '1.2.2-opt')
#baseline = ('XYIntISICType1-frac', '1.2.2-opt')
#baseline = ('Branch on linking variables', '1.2.1-opt')
if len(versions) > 1:
    versionlegend = True
else:
    versionlegend = False
    
if name != '':
    dataSets = [name]
        
for ds in dataSets:
    df_all_solved, df_solved, df_has_soln = dropFilter(df_proc, scenarios, ds)
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
            df_sub = df_all_solved.xs(
                (ds, col), level=["datasets", "fields"], axis=1, drop_level=True
            ).copy()
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
            
            df_time = df_has_soln.xs(
                (ds, "cpu"), level=["datasets", "fields"], axis=1, drop_level=True
            ).copy()
            df_gap = df_has_soln.xs(
                (ds, "gap"), level=["datasets", "fields"], axis=1, drop_level=True
            ).copy()

            print("")
            print("Creating cumutative profile for "+col)
            print("")
            plotCumProf(df_time, df_gap, versions, plotname="cum_" + col + "_" + ds,
                        plottitle="Cumulative Profile"+plottitle,
                        versionlegend = versionlegend)
