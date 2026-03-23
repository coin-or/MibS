# Script to replicate experiments of paper 
# "Improving Directions in Mixed Integer Bilevel Linear Optimization"
# Battista F. & Ted K. Ralphs
# Last edited 2026

# Set up senarios
mibsParamsInputs = {}

commonParams = {
    'Alps_timeLimit': '3600',
    'Alps_msgLevel': '1',                   # any positive integer, increase for verbosity
    "MibS_useFractionalCuts" : "1",         # 0: false, 1: true
    'MibS_useImprovingDirectionPool' : '0', # 0: false, 1: true (EXPERIMENTAL!)
    'MibS_printParameters': '1'
} 

#--------------------------------------------------
# idB&C MILP
#--------------------------------------------------

mibsParamsInputs['idBC-MILP_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "0",
    "MibS_solveSecondLevelWhenXYVarsInt" : "0",
    "MibS_solveSecondLevelWhenXVarsInt" : "0",
    "MibS_solveSecondLevelWhenLVarsFixed" : "0",
    "MibS_solveSecondLevelWhenLVarsInt" : "0",
    "MibS_solveSecondLevelEveryIteration" : "0",
    "MibS_computeBestUBWhenXVarsInt" : "0",
    "MibS_computeBestUBWhenLVarsInt" : "0",
    "MibS_computeBestUBWhenLVarsFixed" : "0"
}

#--------------------------------------------------
# idB&C k = 2 dBnd = 10 Configurations
#--------------------------------------------------

mibsParamsInputs['idBC-LS-k_2-dBnd_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "2",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "0",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

# mibsParamsInputs['idBC-LS-k_2-dBnd_0_10_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "2",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "10",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

mibsParamsInputs['idBC-LS-k_2-dBnd_10_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "2",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "10",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

mibsParamsInputs['idBC-k_2-MILP_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "0",
    "MibS_maxNeighborhoodSize" : "2",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

#--------------------------------------------------
# idB&C k = 3 dBnd = 10 Configurations
#--------------------------------------------------

mibsParamsInputs['idBC-LS-k_3-dBnd_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "3",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "0",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

mibsParamsInputs['idBC-LS-k_3-dBnd_0_10_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "3",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "0",
    "MibS_useLocalSearchDepthUb" : "10",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

mibsParamsInputs['idBC-LS-k_3-dBnd_10_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "3",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "10",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

# mibsParamsInputs['idBC-k_3-MILP_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "0",
#     "MibS_maxNeighborhoodSize" : "3",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

#--------------------------------------------------
# idB&C k = 3 dBnd = 8 Configurations
#--------------------------------------------------

# mibsParamsInputs['idBC-LS-k_3-dBnd_0_8_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "3",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "8",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

mibsParamsInputs['idBC-LS-k_3-dBnd_8_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "3",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "8",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

# #--------------------------------------------------
# # idB&C k = 3 dBnd = 12 Configurations
# #--------------------------------------------------

# mibsParamsInputs['idBC-LS-k_3-dBnd_0_12_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "3",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "12",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

mibsParamsInputs['idBC-LS-k_3-dBnd_12_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "3",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "12",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}


#--------------------------------------------------
# idB&C k = 4 dBnd = 10 Configurations
#--------------------------------------------------

# mibsParamsInputs['idBC-LS-k_4-dBnd_Inf_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "4",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "10000",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

# mibsParamsInputs['idBC-LS-k_4-dBnd_0_10_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "4",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "10",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

mibsParamsInputs['idBC-LS-k_4-dBnd_10_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "4",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "10",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

mibsParamsInputs['idBC-k_4-MILP_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "0",
    "MibS_maxNeighborhoodSize" : "4",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}


#--------------------------------------------------
# idB&C k = 5 dBnd = 10 Configurations
#--------------------------------------------------

# mibsParamsInputs['idBC-LS-k_5-dBnd_Inf_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "5",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "10000",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

# mibsParamsInputs['idBC-LS-k_5-dBnd_0_10_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "5",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "0",
#     "MibS_useLocalSearchDepthUb" : "10",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }

mibsParamsInputs['idBC-LS-k_5-dBnd_10_Inf_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "1",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "5",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "10",
    "MibS_useLocalSearchDepthUb" : "10000",
    "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
    "MibS_solveSecondLevelWhenXVarsInt"  : "0",
    "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
    "MibS_solveSecondLevelWhenLVarsInt"  : "0",
    "MibS_solveSecondLevelEveryIteration" :  "0",
    "MibS_computeBestUBWhenXVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsInt"  : "0",
    "MibS_computeBestUBWhenLVarsFixed"  : "0"
}

# mibsParamsInputs['idBC-k_5-MILP_fracB'] = {
#     "MibS_branchStrategy" : "0",
#     "MibS_useImprovingDirectionOracle" : "1",
#     "MibS_turnOffDefaultCuts" : "1",
#     "MibS_useIntersectionCut" : "1",
#     "MibS_useImprovingSolutionIC" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "0",
#     "MibS_maxNeighborhoodSize" : "5",
#     "MibS_solveSecondLevelWhenXYVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenXVarsInt"  : "0",
#     "MibS_solveSecondLevelWhenLVarsFixed"  : "0",
#     "MibS_solveSecondLevelWhenLVarsInt"  : "0",
#     "MibS_solveSecondLevelEveryIteration" :  "0",
#     "MibS_computeBestUBWhenXVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsInt"  : "0",
#     "MibS_computeBestUBWhenLVarsFixed"  : "0"
# }



#--------------------------------------------------
# MibS 1.2 Configurations
#--------------------------------------------------

# Mibs 1.2 default (IDICs may be turned off) 
mibsParamsInputs['MibS_1_2_defaultB'] = {
    "MibS_useImprovingDirectionOracle" : "0"
}

# Mibs 1.2 default (IDICs)
mibsParamsInputs['MibS_1_2_IDIC_defaultB'] = {
    "MibS_useImprovingDirectionOracle" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "0"
}

# # Mibs 1.2 default (IDICs on with LS)
# mibsParamsInputs['MibS_1_2-LS-k_3_defaultB'] = {
#     "MibS_useImprovingDirectionOracle" : "0",
#     "MibS_useImprovingDirectionIC" : "1",
#     "MibS_improvingDirectionType" : "1",
#     "MibS_maxNeighborhoodSize" : "3",
#     "MibS_maxFeasImprovingDirections" : "10",
#     "MibS_useLocalSearchDepthLb" : "10",
#     "MibS_useLocalSearchDepthUb" : "10000"
# }

# Mibs 1.2 default (IDICs on with LS)
mibsParamsInputs['MibS_1_2-LS-k_2_defaultB'] = {
    "MibS_useImprovingDirectionOracle" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "1",
    "MibS_maxNeighborhoodSize" : "2",
    "MibS_maxFeasImprovingDirections" : "10",
    "MibS_useLocalSearchDepthLb" : "10",
    "MibS_useLocalSearchDepthUb" : "10000"
}

# # Mibs 1.2 (only IDICs)
mibsParamsInputs['MibS_onlyIDIC_fracB'] = {
    "MibS_branchStrategy" : "0",
    "MibS_useImprovingDirectionOracle" : "0",
    "MibS_turnOffDefaultCuts" : "1",
    "MibS_useIntersectionCut" : "1",
    "MibS_useImprovingSolutionIC" : "0",
    "MibS_useImprovingDirectionIC" : "1",
    "MibS_improvingDirectionType" : "0"
}
