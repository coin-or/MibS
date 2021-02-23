# Script to set up parameters for MibS instance path.
# Used to run experiments for diffrent cuts 
# Last edited by yux616
# Feb 2020: params for Cuts

# Executable path and name
exe = '/home/usr/coinbrew/build-mibs/bin/mibs'
pbsfile = 'mibs_batch.pbs'

# Instance path
# Directory name and path containing test instances in .mps format
# Keys are used to name subdirs in output dir
instanceDirs = {
    # 'MIBLP-XU': '/home/usr/coinbrew/MibS/testSets/BilevelLib/general/MIBLP-XU',
    # 'dataIBLP-FIS': '/home/usr/coinbrew/MibS/testSets/BilevelLib/general/IBLP-FIS',
    'dataIBLP-DEN': '/home/usr/coinbrew/MibS/testSets/BilevelLib/general/RANDOM/RAND_BILEVEL'
    }

# Output parent path
outputDir = '/home/usr/coinbrew/MibS/output'

# Set up senarios
mibsParamsInputs = {}

mibsParamsInputs['intNoGood'] = {
        'Alps_timeLimit': '3600',
        'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
        'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

        'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
        'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
        'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
        'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
        'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true

        'MibS_bilevelProblemType': '0',       # 0: general, 1: interdict
        'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
        'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
        'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
        'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
        'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
        'MibS_blisBranchStrategy': '1',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
        'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
        'MibS_branchStrategy': '0',       # 0: fractional, 1: linking
        'MibS_warmStartLL': '0',
        'MibS_maxThreadsLL': '1',
        
        'MibS_usePureIntegerCut': '1',        # 0: false, 1: true
        'MibS_turnOffOtherCuts': '1',
        'MibS_useGeneralNoGoodCut': '0',
        'MibS_useNoGoodCut': '0', 
        
        'MibS_solveSecondLevelWhenXYVarsInt': '0',
        'MibS_solveSecondLevelWhenXVarsInt': '0',
        'MibS_solveSecondLevelWhenLVarsInt': '1',
        'MibS_solveSecondLevelWhenLVarsFixed': '0',

        'MibS_computeBestUBWhenXVarsInt': '0',
        'MibS_computeBestUBWhenLVarsInt': '1',
        'MibS_computeBestUBWhenLVarsFixed': '0',

        'MibS_useLinkingSolutionPool': '1'
    }

mibsParamsInputs['type1IC'] = {
        'Alps_timeLimit': '3600',
        'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
        'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

        'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
        'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
        'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
        'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
        'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true

        'MibS_bilevelProblemType': '0',       # 0: general, 1: interdict
        'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
        'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
        'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
        'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
        'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
        'MibS_blisBranchStrategy': '1',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
        'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
        'MibS_branchStrategy': '0',       # 0: fractional, 1: linking
        'MibS_warmStartLL': '0',
        'MibS_maxThreadsLL': '1',

        'MibS_useTypeIC': '1',
        "MibS_bilevelFreeSetTypeIC": '0',
        'MibS_turnOffOtherCuts': '1',
        'MibS_useGeneralNoGoodCut': '0',
        'MibS_useNoGoodCut': '0', 

        'MibS_solveSecondLevelWhenXYVarsInt': '0',
        'MibS_solveSecondLevelWhenXVarsInt': '0',
        'MibS_solveSecondLevelWhenLVarsInt': '1',
        'MibS_solveSecondLevelWhenLVarsFixed': '0',

        'MibS_computeBestUBWhenXVarsInt': '0',
        'MibS_computeBestUBWhenLVarsInt': '1',
        'MibS_computeBestUBWhenLVarsFixed': '0',

        'MibS_useLinkingSolutionPool': '1'
    }

mibsParamsInputs['type2IC'] = {
        'Alps_timeLimit': '3600',
        'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
        'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

        'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
        'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
        'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
        'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
        'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true

        'MibS_bilevelProblemType': '0',       # 0: general, 1: interdict
        'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
        'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
        'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
        'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
        'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
        'MibS_blisBranchStrategy': '1',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
        'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
        'MibS_branchStrategy': '0',       # 0: fractional, 1: linking
        'MibS_warmStartLL': '0',
        'MibS_maxThreadsLL': '1',

        'MibS_useTypeIC': '1',
        "MibS_bilevelFreeSetTypeIC": '1',
        'MibS_turnOffOtherCuts': '1',
        'MibS_useGeneralNoGoodCut': '0',
        'MibS_useNoGoodCut': '0', 

        'MibS_solveSecondLevelWhenXYVarsInt': '0',
        'MibS_solveSecondLevelWhenXVarsInt': '0',
        'MibS_solveSecondLevelWhenLVarsInt': '1',
        'MibS_solveSecondLevelWhenLVarsFixed': '0',

        'MibS_computeBestUBWhenXVarsInt': '0',
        'MibS_computeBestUBWhenLVarsInt': '1',
        'MibS_computeBestUBWhenLVarsFixed': '0',

        'MibS_useLinkingSolutionPool': '1'
    }

mibsParamsInputs['watermelonIC'] = {
        'Alps_timeLimit': '3600',
        'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
        'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

        'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
        'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
        'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
        'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
        'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true

        'MibS_bilevelProblemType': '0',       # 0: general, 1: interdict
        'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
        'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
        'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
        'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
        'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
        'MibS_blisBranchStrategy': '1',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
        'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
        'MibS_branchStrategy': '0',       # 0: fractional, 1: linking
        'MibS_warmStartLL': '0',
        'MibS_maxThreadsLL': '1',

        "MibS_useTypeWatermelon": '1',
        'MibS_turnOffOtherCuts': '1',
        'MibS_useGeneralNoGoodCut': '0',
        'MibS_useNoGoodCut': '0', 

        'MibS_solveSecondLevelWhenXYVarsInt': '0',
        'MibS_solveSecondLevelWhenXVarsInt': '0',
        'MibS_solveSecondLevelWhenLVarsInt': '1',
        'MibS_solveSecondLevelWhenLVarsFixed': '0',

        'MibS_computeBestUBWhenXVarsInt': '0',
        'MibS_computeBestUBWhenLVarsInt': '1',
        'MibS_computeBestUBWhenLVarsFixed': '0',

        'MibS_useLinkingSolutionPool': '1'
    }

mibsParamsInputs['hyperIC'] = {
        'Alps_timeLimit': '3600',
        'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
        'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

        'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
        'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
        'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
        'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
        'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true

        'MibS_bilevelProblemType': '0',       # 0: general, 1: interdict
        'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
        'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
        'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
        'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
        'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
        'MibS_blisBranchStrategy': '1',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
        'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
        'MibS_branchStrategy': '0',       # 0: fractional, 1: linking
        'MibS_warmStartLL': '0',
        'MibS_maxThreadsLL': '1',

        "MibS_useTypeHypercubeIC": '1',
        'MibS_turnOffOtherCuts': '1',
        'MibS_useGeneralNoGoodCut': '0',
        'MibS_useNoGoodCut': '0', 

        'MibS_solveSecondLevelWhenXYVarsInt': '0',
        'MibS_solveSecondLevelWhenXVarsInt': '0',
        'MibS_solveSecondLevelWhenLVarsInt': '1',
        'MibS_solveSecondLevelWhenLVarsFixed': '0',

        'MibS_computeBestUBWhenXVarsInt': '0',
        'MibS_computeBestUBWhenLVarsInt': '1',
        'MibS_computeBestUBWhenLVarsFixed': '0',

        'MibS_useLinkingSolutionPool': '1'
    }
# only to 'dataIBLP-FIS'
mibsParamsInputs['genNoGood'] = {
        'Alps_timeLimit': '3600',
        'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
        'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

        'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
        'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
        'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
        'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
        'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true

        'MibS_bilevelProblemType': '0',       # 0: general, 1: interdict
        'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
        'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
        'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
        'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
        'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
        'MibS_blisBranchStrategy': '1',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
        'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
        'MibS_branchStrategy': '0',       # 0: fractional, 1: linking
        'MibS_warmStartLL': '0',
        'MibS_maxThreadsLL': '1',

        'MibS_useGeneralNoGoodCut': '1',
        'MibS_turnOffOtherCuts': '1',
        'MibS_useNoGoodCut': '0', 

        'MibS_solveSecondLevelWhenXYVarsInt': '0',
        'MibS_solveSecondLevelWhenXVarsInt': '0',
        'MibS_solveSecondLevelWhenLVarsInt': '1',
        'MibS_solveSecondLevelWhenLVarsFixed': '0',

        'MibS_computeBestUBWhenXVarsInt': '0',
        'MibS_computeBestUBWhenLVarsInt': '1',
        'MibS_computeBestUBWhenLVarsFixed': '0',

        'MibS_useLinkingSolutionPool': '1'
    }


# For reference
# All MibS parameters; see comments for specific params 
mibsParamsLib = {
    'Alps_timeLimit': '3600',
    'Alps_nodeLimit': '4000000',
    'Alps_nodeLimit': '20000',

    'Alps_msgLevel': '1000',
    'Alps_hubMsgLevel': '0',
    'Alps_workerMsgLevel': '0',
    'Alps_logFileLevel': '3',

    'Alps_hubNum': '2',

    'Alps_interClusterBalance': '1',      # 1: balancing load, 0: don't.
    'Alps_intraClusterBalance': '1',      # 1: balancing load, 0: don't.

    'Alps_searchStrategy': '4',   # 0: Best, 1: Best-est, 2: Breath, 3: Depth, 4 hybrid
    'Alps_searchStrategyRampUp': '0',
    'Alps_nodeLogInterval': '1',
    'Alps_hubWorkClusterSizeLimit': '2',
    'Alps_masterInitNodeNum': '8',
    'Alps_hubInitNodeNum': '16',
    'Alps_unitWorkNodes': '100',        # or unit time
    'Alps_unitWorkTime': '0.03',
    'Alps_needWorkThreshold': '0.5',
    'Alps_changeWorkThreshold': '0.10',
    'Alps_donorThreshold': '0.10',
    'Alps_receiverThreshold': '0.10',
    'Alps_masterBalancePeriod': '0.3',
    'Alps_hubReportPeriod': '0.5',

    'Blis_cutRampUp': '1',      # true(1) or false(0)
    'Blis_cutStrategy': '3',
    'Blis_cutGenerationFrequency': '1',

    'Blis_cutCliqueStrategy': '0',
    'Blis_cutGomoryStrategy': '3',
    'Blis_cutFlowCoverStrategy': '0',
    'Blis_cutKnapsackStrategy': '0',

    'Blis_cutMirStrategy': '0',

    'Blis_cutOddHoleStrategy': '0',
    'Blis_cutProbingStrategy': '0',

    'Blis_cutTwoMirStrategy': '0',

    'Blis_heurStrategy': '0',   # -2: root, -1: auto, 0: disable, any positive integer
    'Blis_heurRound': '0',   # -2: root, -1: auto, 0: disable, any positive integer

    'Blis_branchStrategy': '0',   # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
    'Blis_branchStrategyRampUp': '0',

    'Blis_pseudoWeight': '0.8',  # [0.0, 1.0]
    'Blis_pseudoRelibility': '8',
    'Blis_lookAhead': '4',
    'Blis_denseConFactor': '5.',
    'Blis_scaleConFactor': '100000000.0',
    'Blis_difference': '1',

    'Blis_sharePseudocostRampUp': '1',
    'Blis_sharePseudocostSearch': '1',

    'Blis_checkMemory': '1',

    'MibS_usePreprocessor': '0',          # -1: auto, 0: false, 1: true
    'MibS_useLowerObjHeuristic': '0',     # -1: auto, 0: false, 1: true
    'MibS_useObjCutHeuristic': '0',       # -1: auto, 0: false, 1: true
    'MibS_useWSHeuristic': '0',           # -1: auto, 0: false, 1: true
    'MibS_useGreedyHeuristic': '0',       # 0: false, 1: true
    'MibS_usePureIntegerCut': '1',        # 0: false, 1: true
    'MibS_useValFuncCut': '0',            # 0: false, 1: true
    'MibS_useNoGoodCut': '0',             # 0: false, 1: true
    'MibS_useIncObjCut': '0',             # 0: false, 1: true
    'MibS_bilevelProblemType': '1',       # 0: general, 1: interdict
    'MibS_bilevelCutTypes': '0',          # 0: general, 1: general/interdict, 2: #general/binary UL, 3: binary UL
    'MibS_whichActiveConMethod': '1',    # 0: simple, 1: basis
    'MibS_cutStrategy': '2',              # 0: branch only, 1: cut only, 2: use cut and branch
    'MibS_objBoundStrategy': '0',         # 0: LL obj bound, 1: interdiction bound
    'MibS_blisCutStrategy': '-1',                 # -2: root, -1: auto, 0: disable, any positive integer
    'MibS_blisBranchStrategy': '0',      # 0: max inf, 1: pseudocost, 2: relibility, 3: strong
    'MibS_upperFileFormat': '0',         # 0: mps, 1: AMPL/GMPL (not working, yet)
    'MibS_branchStrategy': '1',       # 0: fractional, 1: linking
    'MibS_warmStartLL': '0',
    'MibS_maxThreadsLL': '1',

    'MibS_useBendersCut': '0',
    'MibS_useGeneralNoGoodCut': '0',
    'MibS_useIntersectionCut': '0',     # Not used; commented off in code
    'MibS_intersectionCutType': '1',    # not used aswell?
    
    "MibS_useTypeIC": '0',              # PARAM OFF?
    "MibS_bilevelFreeSetTypeIC": '0',  # 0: Intersection Cut Type I; 1: Intersection Cut Type II; only valid when use intersectionCutTypeIC
    "MibS_useTypeWatermelon": '0',
    "MibS_useTypeHypercubeIC": '0',
    "MibS_useTypeTenderIC": '0',
    "MibS_useTypeHybridIC": '0',

    'MibS_useBoundCut': '0',
    'MibS_boundCutOptimal': '1',
    'MibS_boundCutRelaxUpper': '0',    # This option doesn't work for now
    'MibS_turnOffOtherCuts': '0',      # MibS cuts not set will be turned off: not affect GeneralNoGoodCut(maybe bug ard L6816)

    'MibS_whichCutsLL': '2',           # 0: no cuts, 1: gomory only, 2: all cuts
    'MibS_doDualFixing': '0',
    'MibS_feasCheckSolver': 'SYMPHONY',

    'MibS_solveSecondLevelWhenXYVarsInt': '0',
    'MibS_solveSecondLevelWhenXVarsInt': '0',
    'MibS_solveSecondLevelWhenLVarsInt': '1',
    'MibS_solveSecondLevelWhenLVarsFixed': '0',

    'MibS_computeBestUBWhenXVarsInt': '0',
    'MibS_computeBestUBWhenLVarsInt': '1',
    'MibS_computeBestUBWhenLVarsFixed': '0',

    'MibS_useLinkingSolutionPool': '1',
    'MibS_slTargetGap': '-1' # -1: OFF, 0-100: set SL target gap to given
}