# Script to set up parameters for MibS instance path.
# Used to run experiments for diffrent cuts 
# Last edited by yux616
# Apr 2020

# Executable path and name
pbsfile = '/home/ted/Projects/MibS/scripts/analyze/mibs_batch.pbs'

# Instance path
# Directory name and path containing test instances in .mps format
# Keys are used to name subdirs in output dir
instanceDirs = {
    # 'MIBLP-XU': '/home/usr/coinbrew/MibS/testSets/BilevelLib/general/MIBLP-XU',
    #'IBLP-FIS': '/home/sat214/cutPaperMibS/dataPaper/iblpFis',
    #'INTERD-DEN': '/home/sat214/cutPaperMibS/dataPaper/interdDen',
    #'IBLP-DEN': '/home/sat214/cutPaperMibS/dataPaper/iblpDen',
    'IBLP-ZHANG': '/home/sat214/MIBS/ozaltinData/convertedData/Testbed1/'
}

#versions = ['1.1', 'ib']
#versions = ['1.2-opt']
versions = ['pr-92']

# Output parent path
outputDir = '/home/ted/Projects/MibS/output'

# Name
testname = 'mibs'

# Set up senarios
mibsParamsInputs = {}

commonParams = {
    'Alps_timeLimit': '3600',
    'Blis_scaleConFactor': '10000000000',
    'Blis_heurStrategy': '0',             # -2: root, -1: auto, 0: disable, any positive integer
    'Blis_heurRound': '0',                # -2: root, -1: auto, 0: disable, any positive integer
    'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
    'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
    'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
    'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
    'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true
    'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
    'MibS_blisCutStrategy': '-1',         # -2: root, -1: auto, 0: disable, any positive integer
#    'MibS_warmStartLL': '0',
    'MibS_maxThreadsLL': '1',
    'MibS_allowRemoveCut': '0',           # 0: false, 1: true
    'MibS_whichCutsLL': '2',              # 0: no cuts, 1: gomory only, 2: all cuts
    'MibS_doDualFixing': '0',
#    'MibS_feasCheckSolver': 'SYMPHONY',
    'MibS_solveSecondLevelWhenXYVarsInt': '1',
    'MibS_solveSecondLevelWhenXVarsInt': '0',
    'MibS_solveSecondLevelWhenLVarsFixed': '1',
    'MibS_computeBestUBWhenXVarsInt': '0',
    'MibS_computeBestUBWhenLVarsInt': '0',
    'MibS_computeBestUBWhenLVarsFixed': '1',
    'MibS_useLinkingSolutionPool': '1',
}    

# mibsParamsInputs['default'] = {
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['default-frac'] = {
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['noGood+type1+pureInteger'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useGeneralNoGoodCut':'1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_usePureIntegerCut':'1',
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',
# }

# mibsParamsInputs['watermelon+type1+incObj-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useIncObjCut':'1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',
# }

# mibsParamsInputs['watermelon+type1+incObj-frac-LV'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useIncObjCut':'1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',
# }

# mibsParamsInputs['watermelon+type1-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',
# }

# mibsParamsInputs['watermelon+type1-frac-LV'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',
# }

# mibsParamsInputs['benders'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useBendersCut':'1',
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',
# }

# mibsParamsInputs['benders-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useBendersCut':'1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['benders-cplex'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useBendersCut':'1',
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',
#     'MibS_feasCheckSolver' : 'CPLEX'
# }

# mibsParamsInputs['benders-frac-cplex'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useBendersCut':'1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX'
# }

mibsParamsInputs['incObj'] = {
    'MibS_turnOffOtherCuts': '1',
    'MibS_useIncObjCut':'1',
    'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
    'MibS_solveSecondLevelWhenLVarsInt': '1',
}

# mibsParamsInputs['incObj-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useIncObjCut':'1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',
# }

mibsParamsInputs['genNoGood'] = {
    'MibS_turnOffOtherCuts': '1',
    'MibS_useGeneralNoGoodCut':'1',
    'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
    'MibS_solveSecondLevelWhenLVarsInt': '1',           # 0: fractional, 1: linking
}

# mibsParamsInputs['genNoGood-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useGeneralNoGoodCut':'1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['watermelon'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['watermelon-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['watermelon-frac-LV'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',           # 0: fractional, 1: linking
# }

mibsParamsInputs['fracWatermelon'] = {
    'MibS_turnOffOtherCuts': '1',
    'MibS_useTypeWatermelon': '1',
    'MibS_useFractionalCuts': '1',
    'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
    'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
}

# mibsParamsInputs['fracWatermelon-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useFractionalCuts': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['fracWatermelon-frac-cplex'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useFractionalCuts': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX'
# }

# mibsParamsInputs['fracWatermelon-frac-cplex-d10'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useFractionalCuts': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX',
#     'MibS_maxCutDepth' : '10'
# }

# mibsParamsInputs['fracWatermelon-frac-cplex-d20'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeWatermelon': '1',
#     'MibS_useFractionalCuts': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX',
#     'MibS_maxCutDepth' : '20'
# }

# mibsParamsInputs['type1'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['type1-WS'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',           # 0: fractional, 1: linking
#     'MibS_warmStartLL': '1'
# }

# mibsParamsInputs['type1-cplex'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX'
# }

# mibsParamsInputs['fracType1-frac-cplex'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelEveryIteration': '1',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX'
# }

# mibsParamsInputs['fracType1-frac-cplex-d10'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelEveryIteration': '1',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX',
#     'MibS_maxCutDepth' : '10'
# }

# mibsParamsInputs['fracType1-frac-cplex-d20'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelEveryIteration': '1',           # 0: fractional, 1: linking
#     'MibS_feasCheckSolver' : 'CPLEX',
#     'MibS_maxCutDepth' : '20'
# }

# mibsParamsInputs['type1-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '0',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['type2'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '1',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['type2-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeIC': '1',
#     'MibS_bilevelFreeSetTypeIC': '1',     # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['hyper'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeHypercubeIC': '1',
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '1',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['hyper-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_useTypeHypercubeIC': '1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['pureInteger'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_usePureIntegerCut':'1',
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['pureInteger-frac'] = {
#     'MibS_turnOffOtherCuts': '1',
#     'MibS_usePureIntegerCut':'1',
#     'MibS_branchStrategy': '0',           # 0: fractional, 1: linking
#     'MibS_solveSecondLevelWhenLVarsInt': '0',           # 0: fractional, 1: linking
# }

# mibsParamsInputs['noCut'] = {
#      'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#      'MibS_cutStrategy': '0'
# }

# mibsParamsInputs['noCut-cplex'] = {
#      'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#      'MibS_cutStrategy': '0',
#      'MibS_feasCheckSolver' : 'CPLEX',
# }

# mibsParamsInputs['noCut-WS'] = {
#     'MibS_branchStrategy': '1',           # 0: fractional, 1: linking
#     'MibS_cutStrategy': '0',
#     'MibS_warmStartLL': '1'
# }

