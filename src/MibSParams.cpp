/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*          Sahar Tahernajad, Lehigh University                              */
/*                                                                           */
/* Copyright (C) 2007-2019 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#include "MibSParams.hpp"
#include "MibSConstants.hpp"
#include "MibSModel.hpp"

using std::make_pair;

void 
MibSParams::createKeywordList() {
   
   //--------------------------------------------------------
   // Create the list of keywords for parameter file reading
   //--------------------------------------------------------
   
   //--------------------------------------------------------
   // BoolPar
   //--------------------------------------------------------
  
   keys_.push_back(make_pair(std::string("MibS_useValFuncCut"),
			     AlpsParameter(AlpsBoolPar, useValFuncCut)));

   keys_.push_back(make_pair(std::string("MibS_useBoundCut"),
			     AlpsParameter(AlpsBoolPar, useBoundCut)));
   
   keys_.push_back(make_pair(std::string("MibS_boundCutOptimal"),
			     AlpsParameter(AlpsBoolPar, boundCutOptimal)));
   
   keys_.push_back(make_pair(std::string("MibS_boundCutRelaxUpper"),
			     AlpsParameter(AlpsBoolPar, boundCutRelaxUpper)));

   keys_.push_back(make_pair(std::string("MibS_useIpBound"),
			     AlpsParameter(AlpsBoolPar, useIpBound)));
   
   keys_.push_back(make_pair(std::string("MibS_isBilevelBranchProb"),
			     AlpsParameter(AlpsBoolPar, isBilevelBranchProb)));

   keys_.push_back(make_pair(std::string("MibS_warmStartLL"),
			     AlpsParameter(AlpsBoolPar, warmStartLL)));

   keys_.push_back(make_pair(std::string("MibS_doDualFixing"),
			     AlpsParameter(AlpsBoolPar, doDualFixing)));

   keys_.push_back(make_pair(std::string("MibS_useUBDecompose"),
			     AlpsParameter(AlpsBoolPar, useUBDecompose)));
   
   keys_.push_back(make_pair(std::string("MibS_turnOffOtherCuts"),
			     AlpsParameter(AlpsBoolPar, turnOffOtherCuts)));

   keys_.push_back(make_pair(std::string("MibS_printProblemInfo"),
			     AlpsParameter(AlpsBoolPar, printProblemInfo)));

   keys_.push_back(make_pair(std::string("MibS_allowRemoveCut"),
			     AlpsParameter(AlpsBoolPar, allowRemoveCut)));

   keys_.push_back(make_pair(std::string("MibS_useNewPureIntCut"),
			     AlpsParameter(AlpsBoolPar, useNewPureIntCut)));

   //the parameter for using progressive hedging
   //note that this heuristic can be used only when
   //the parameter "stochasticityType" is set to "stochasticWithoutSAA"
   keys_.push_back(make_pair(std::string("MibS_useProgresHedg"),
			     AlpsParameter(AlpsBoolPar, useProgresHedg)));
   

   //--------------------------------------------------------
   // BoolArrayPar
   //--------------------------------------------------------
   
   //--------------------------------------------------------
   // Int Parameters
   //--------------------------------------------------------
   
   //keys_.push_back(make_pair(std::string("MibS_verbosity"),
   //		    AlpsParameter(AlpsIntPar, verbosity)));
   
   keys_.push_back(make_pair(std::string("MibS_whichActiveConMethod"),
			     AlpsParameter(AlpsIntPar, 
					   whichActiveConMethod)));
   
   keys_.push_back(make_pair(std::string("MibS_maxNumActiveCons"),
			     AlpsParameter(AlpsIntPar, maxNumActiveCons)));
   
   //Might want to make this a string param
   keys_.push_back(make_pair(std::string("MibS_bilevelProblemType"),
			     AlpsParameter(AlpsIntPar, bilevelProblemType)));

   keys_.push_back(make_pair(std::string("MibS_bilevelCutTypes"),
			     AlpsParameter(AlpsIntPar, bilevelCutTypes)));
   
   keys_.push_back(make_pair(std::string("MibS_cutStrategy"),
			     AlpsParameter(AlpsIntPar, cutStrategy)));

   keys_.push_back(make_pair(std::string("MibS_objBoundStrategy"),
			     AlpsParameter(AlpsIntPar, objBoundStrategy)));

   keys_.push_back(make_pair(std::string("MibS_blisCutStrategy"),
			     AlpsParameter(AlpsIntPar, blisCutStrategy)));

   keys_.push_back(make_pair(std::string("MibS_blisBranchStrategy"),
			     AlpsParameter(AlpsIntPar, blisBranchStrategy)));
   
   keys_.push_back(make_pair(std::string("MibS_branchStrategy"),
			     AlpsParameter(AlpsIntPar, branchStrategy)));
 
   keys_.push_back(make_pair(std::string("MibS_upperFileFormat"),
			     AlpsParameter(AlpsIntPar, upperFileFormat)));

   keys_.push_back(make_pair(std::string("MibS_maxThreadsLL"),
			     AlpsParameter(AlpsIntPar, maxThreadsLL)));

   keys_.push_back(make_pair(std::string("MibS_whichCutsLL"),
			     AlpsParameter(AlpsIntPar, whichCutsLL)));

   keys_.push_back(make_pair(std::string("MibS_usePreprocessor"),
			     AlpsParameter(AlpsIntPar, usePreprocessor)));
   
   keys_.push_back(make_pair(std::string("MibS_useLowerObjHeuristic"),
			     AlpsParameter(AlpsIntPar, useLowerObjHeuristic)));

   keys_.push_back(make_pair(std::string("MibS_useObjCutHeuristic"),
			     AlpsParameter(AlpsIntPar, useObjCutHeuristic)));

   keys_.push_back(make_pair(std::string("MibS_useWSHeuristic"),
			     AlpsParameter(AlpsIntPar, useWSHeuristic)));

   keys_.push_back(make_pair(std::string("MibS_useGreedyHeuristic"),
			     AlpsParameter(AlpsIntPar, useGreedyHeuristic)));

   keys_.push_back(make_pair(std::string("MibS_usePureIntegerCut"),
			     AlpsParameter(AlpsIntPar, usePureIntegerCut)));

   keys_.push_back(make_pair(std::string("MibS_useNoGoodCut"),
			     AlpsParameter(AlpsIntPar, useNoGoodCut)));

   keys_.push_back(make_pair(std::string("MibS_useGeneralNoGoodCut"),
			     AlpsParameter(AlpsIntPar, useGeneralNoGoodCut)));
   
   keys_.push_back(make_pair(std::string("MibS_useIncObjCut"),
			     AlpsParameter(AlpsIntPar, useIncObjCut)));

   keys_.push_back(make_pair(std::string("MibS_useBendersCut"),
			     AlpsParameter(AlpsIntPar, useBendersCut)));

   keys_.push_back(make_pair(std::string("MibS_bendersCutType"),
			     AlpsParameter(AlpsIntPar, bendersCutType)));

   keys_.push_back(make_pair(std::string("MibS_useIntersectionCut"),
                             AlpsParameter(AlpsIntPar, useIntersectionCut)));

   //keys_.push_back(make_pair(std::string("MibS_intersectionCutType"),
   //			     AlpsParameter(AlpsIntPar, intersectionCutType)));

   keys_.push_back(make_pair(std::string("MibS_useTypeIC"),
			     AlpsParameter(AlpsIntPar, useTypeIC)));

   keys_.push_back(make_pair(std::string("MibS_useTypeWatermelon"),
			     AlpsParameter(AlpsIntPar, useTypeWatermelon)));

   keys_.push_back(make_pair(std::string("MibS_useTypeHypercubeIC"),
			     AlpsParameter(AlpsIntPar, useTypeHypercubeIC)));

   keys_.push_back(make_pair(std::string("MibS_useTypeTenderIC"),
			     AlpsParameter(AlpsIntPar, useTypeTenderIC)));

   keys_.push_back(make_pair(std::string("MibS_useTypeHybridIC"),
			     AlpsParameter(AlpsIntPar, useTypeHybridIC)));

   keys_.push_back(make_pair(std::string("MibS_bilevelFreeSetTypeIC"),
			     AlpsParameter(AlpsIntPar, bilevelFreeSetTypeIC)));

   //solve lower-level Parameters
   keys_.push_back(make_pair(std::string("MibS_solveSecondLevelWhenXYVarsInt"),
			     AlpsParameter(AlpsIntPar, solveSecondLevelWhenXYVarsInt)));

   keys_.push_back(make_pair(std::string("MibS_solveSecondLevelWhenXVarsInt"),
			     AlpsParameter(AlpsIntPar, solveSecondLevelWhenXVarsInt)));

   keys_.push_back(make_pair(std::string("MibS_solveSecondLevelWhenLVarsInt"),
			     AlpsParameter(AlpsIntPar, solveSecondLevelWhenLVarsInt)));

   keys_.push_back(make_pair(std::string("MibS_solveSecondLevelWhenLVarsFixed"),
			     AlpsParameter(AlpsIntPar, solveSecondLevelWhenLVarsFixed)));

   //solve problem UB
   keys_.push_back(make_pair(std::string("MibS_computeBestUBWhenXVarsInt"),
			     AlpsParameter(AlpsIntPar, computeBestUBWhenXVarsInt)));

   keys_.push_back(make_pair(std::string("MibS_computeBestUBWhenLVarsInt"),
			     AlpsParameter(AlpsIntPar, computeBestUBWhenLVarsInt)));

   keys_.push_back(make_pair(std::string("MibS_computeBestUBWhenLVarsFixed"),
			     AlpsParameter(AlpsIntPar, computeBestUBWhenLVarsFixed)));

   keys_.push_back(make_pair(std::string("MibS_useLinkingSolutionPool"),
			     AlpsParameter(AlpsIntPar, useLinkingSolutionPool)));

   keys_.push_back(make_pair(std::string("MibS_newPureIntCutDepthLb"),
			     AlpsParameter(AlpsIntPar, newPureIntCutDepthLb)));

   keys_.push_back(make_pair(std::string("MibS_newPureIntCutDepthUb"),
			     AlpsParameter(AlpsIntPar, newPureIntCutDepthUb)));

   keys_.push_back(make_pair(std::string("MibS_boundCutOptimalType"),
			     AlpsParameter(AlpsIntPar, boundCutOptimalType)));

   keys_.push_back(make_pair(std::string("MibS_boundCutDepthLb"),
			     AlpsParameter(AlpsIntPar, boundCutDepthLb)));

   keys_.push_back(make_pair(std::string("MibS_boundCutDepthUb"),
			     AlpsParameter(AlpsIntPar, boundCutDepthUb)));

   keys_.push_back(make_pair(std::string("MibS_boundCutFreq"),
			     AlpsParameter(AlpsIntPar, boundCutFreq)));

   keys_.push_back(make_pair(std::string("MibS_boundCutNodeLim"),
			     AlpsParameter(AlpsIntPar, boundCutNodeLim)));

   keys_.push_back(make_pair(std::string("MibS_relaxTypeParamBoundCut"),
   			     AlpsParameter(AlpsIntPar, relaxTypeParamBoundCut)));

   keys_.push_back(make_pair(std::string("MibS_maxActiveNodes"),
			     AlpsParameter(AlpsIntPar, maxActiveNodes)));

   //parameters for stochastic and SAA
   //this parameter should be set to false, when the problem is
   //stochastic and A2 is not random.
   keys_.push_back(make_pair(std::string("MibS_isA2Random"),
			     AlpsParameter(AlpsIntPar, isA2Random)));

   //this parameter should be set to true when the smps format
   //is used for the stochastic case. 
   keys_.push_back(make_pair(std::string("MibS_isSMPSFormat"),
			     AlpsParameter(AlpsIntPar, isSMPSFormat)));
   //N
   keys_.push_back(make_pair(std::string("MibS_sampSizeSAA"),
			     AlpsParameter(AlpsIntPar, sampSizeSAA)));
   //N'
   keys_.push_back(make_pair(std::string("MibS_evalSampSizeSAA"),
			     AlpsParameter(AlpsIntPar, evalSampSizeSAA)));
   //M
   keys_.push_back(make_pair(std::string("MibS_replNumSAA"),
			     AlpsParameter(AlpsIntPar, replNumSAA)));

   keys_.push_back(make_pair(std::string("MibS_lbDistB2SAA"),
			     AlpsParameter(AlpsIntPar, lbDistB2SAA)));
   
   keys_.push_back(make_pair(std::string("MibS_ubDistB2SAA"),
			     AlpsParameter(AlpsIntPar, ubDistB2SAA)));

   keys_.push_back(make_pair(std::string("MibS_lbDistA2SAA"),
			     AlpsParameter(AlpsIntPar, lbDistA2SAA)));

   keys_.push_back(make_pair(std::string("MibS_ubDistA2SAA"),
			     AlpsParameter(AlpsIntPar, ubDistA2SAA)));

   //it is assumed that incDistB2NumerSAA <= incDistB2DenumSAA
   //and they are relatively prime (the same for A2)
   //it is assumed that numerators are 1  
   keys_.push_back(make_pair(std::string("MibS_incDistB2NumerSAA"),
			     AlpsParameter(AlpsIntPar, incDistB2NumerSAA)));

   keys_.push_back(make_pair(std::string("MibS_incDistB2DenumSAA"),
			     AlpsParameter(AlpsIntPar, incDistB2DenumSAA)));

   keys_.push_back(make_pair(std::string("MibS_incDistA2NumerSAA"),
			     AlpsParameter(AlpsIntPar, incDistA2NumerSAA)));

   keys_.push_back(make_pair(std::string("MibS_incDistA2DenumSAA"),
			     AlpsParameter(AlpsIntPar, incDistA2DenumSAA)));

   //parameters for progressive hedging
   keys_.push_back(make_pair(std::string("MibS_iterationLimitPH"),
			     AlpsParameter(AlpsIntPar, iterationLimitPH)));

   keys_.push_back(make_pair(std::string("MibS_nodeLimitPHSubprob"),
			     AlpsParameter(AlpsIntPar, nodeLimitPHSubprob)));

   //--------------------------------------------------------
   // String Parameters.
   //--------------------------------------------------------
   
   keys_.push_back(make_pair(std::string("MibS_auxiliaryInfoFile"),
			     AlpsParameter(AlpsStringPar, auxiliaryInfoFile)));

   keys_.push_back(make_pair(std::string("MibS_auxiliaryTimFile"),
			     AlpsParameter(AlpsStringPar, auxiliaryTimFile)));

   keys_.push_back(make_pair(std::string("MibS_auxiliaryStoFile"),
			     AlpsParameter(AlpsStringPar, auxiliaryStoFile)));

   keys_.push_back(make_pair(std::string("MibS_feasCheckSolver"),
			     AlpsParameter(AlpsStringPar, feasCheckSolver)));

   keys_.push_back(make_pair(std::string("MibS_inputFormat"),
			     AlpsParameter(AlpsStringPar, inputFormat)));

   //saharStoc:this parameter should be set to
   //"deterministic" : deterministic problems
   //"stochasticWithSAA": stochastic problem and is solved by SAA method 
   //"stochasticWithoutSAA": stochastic problem, but SAA method is not used.
   //When the SAA method is not used, the format should be SMPS necessarily.
   //So the parameter "isSMPSFormat" should be set to true.
   keys_.push_back(make_pair(std::string("MibS_stochasticityType"),
			     AlpsParameter(AlpsStringPar, stochasticityType)));

   //--------------------------------------------------------
   // Double Parameters.
   //--------------------------------------------------------

   keys_.push_back(make_pair(std::string("MibS_boundCutTimeLim"),
			     AlpsParameter(AlpsDoublePar, boundCutTimeLim)));

   //For progressive hedging
   keys_.push_back(make_pair(std::string("MibS_optimalRelGapLimitPHSubprob"),
			     AlpsParameter(AlpsDoublePar, optimalRelGapLimitPHSubprob)));

}

//#############################################################################

void 
MibSParams::setDefaultEntries() {
   
   //-------------------------------------------------------------
   // Bool Parameters.
   //-------------------------------------------------------------

   setEntry(useValFuncCut, false);
     
   setEntry(useBoundCut, false);
   
   setEntry(boundCutOptimal, true);
   
   setEntry(boundCutRelaxUpper, false);

   setEntry(useIpBound, false);

   setEntry(isBilevelBranchProb, false);

   setEntry(warmStartLL, false);

   setEntry(doDualFixing, false);

   setEntry(useUBDecompose, false);

   setEntry(turnOffOtherCuts, false);

   setEntry(printProblemInfo, true);

   setEntry(allowRemoveCut, false);

   setEntry(useNewPureIntCut, false);

   setEntry(useProgresHedg, false);

   //-------------------------------------------------------------
   // Int Parameters.
   //-------------------------------------------------------------
   
   //   setEntry(verbosity, 0);
   
   setEntry(whichActiveConMethod, SIMPLE);

   setEntry(maxNumActiveCons, BIGCONSTANT);

   setEntry(bilevelProblemType, PARAM_NOTSET);

   setEntry(bilevelCutTypes, GENERALONLY);

   setEntry(cutStrategy, CUTONLY);

   setEntry(objBoundStrategy, LPBOUND);

   setEntry(blisCutStrategy, 0);

   setEntry(blisBranchStrategy, 0);

   setEntry(branchStrategy, MibSBranchingStrategyNotSet);

   setEntry(upperFileFormat, 0);

   setEntry(maxThreadsLL, 1);

   setEntry(whichCutsLL, 2);

   setEntry(usePreprocessor, PARAM_NOTSET);

   setEntry(useLowerObjHeuristic, PARAM_OFF);

   setEntry(useObjCutHeuristic, PARAM_OFF);

   setEntry(useWSHeuristic, PARAM_OFF);

   setEntry(useGreedyHeuristic, PARAM_NOTSET);

   setEntry(usePureIntegerCut, PARAM_NOTSET);

   setEntry(useNoGoodCut, PARAM_NOTSET);

   setEntry(useGeneralNoGoodCut, PARAM_NOTSET);

   setEntry(useIncObjCut, PARAM_NOTSET);

   setEntry(useBendersCut, PARAM_NOTSET);

   setEntry(bendersCutType, MibSBendersCutTypeJustOneCut);

   setEntry(useIntersectionCut, PARAM_NOTSET);

   //setEntry(intersectionCutType, MibSIntersectionCutTypeNotSet);

   setEntry(useTypeIC, PARAM_NOTSET);

   setEntry(useTypeWatermelon, PARAM_NOTSET);

   setEntry(useTypeHypercubeIC, PARAM_NOTSET);

   setEntry(useTypeTenderIC, PARAM_NOTSET);

   setEntry(useTypeHybridIC, PARAM_NOTSET);

   setEntry(bilevelFreeSetTypeIC, MibSBilevelFreeSetTypeICNotSet);

   setEntry(solveSecondLevelWhenXYVarsInt, PARAM_NOTSET);

   setEntry(solveSecondLevelWhenXVarsInt, PARAM_NOTSET);

   setEntry(solveSecondLevelWhenLVarsInt, PARAM_NOTSET);

   setEntry(solveSecondLevelWhenLVarsFixed, PARAM_NOTSET);

   setEntry(computeBestUBWhenXVarsInt, PARAM_NOTSET);

   setEntry(computeBestUBWhenLVarsInt, PARAM_NOTSET);

   setEntry(computeBestUBWhenLVarsFixed, PARAM_NOTSET);

   setEntry(useLinkingSolutionPool, PARAM_NOTSET);

   setEntry(newPureIntCutDepthLb, -1);

   setEntry(newPureIntCutDepthUb, -1);

   setEntry(boundCutOptimalType, MibSBoundCutOptimalTypeParametric);

   setEntry(boundCutDepthLb, -1);

   setEntry(boundCutDepthUb, -1);

   setEntry(boundCutFreq, 1);

   setEntry(boundCutNodeLim, ALPS_INT_MAX);

   setEntry(relaxTypeParamBoundCut, MibSRelaxTypeParamBoundCutLP);

   setEntry(maxActiveNodes, 1);

   setEntry(isA2Random, PARAM_NOTSET);

   setEntry(isSMPSFormat, PARAM_NOTSET);

   setEntry(sampSizeSAA, PARAM_NOTSET);

   setEntry(evalSampSizeSAA, PARAM_NOTSET);

   setEntry(replNumSAA, PARAM_NOTSET);

   setEntry(lbDistB2SAA, PARAM_NOTSET);

   setEntry(ubDistB2SAA, PARAM_NOTSET);

   setEntry(lbDistA2SAA, PARAM_NOTSET);

   setEntry(ubDistA2SAA, PARAM_NOTSET);

   setEntry(incDistB2NumerSAA, PARAM_NOTSET);

   setEntry(incDistB2DenumSAA, PARAM_NOTSET);

   setEntry(incDistA2NumerSAA, PARAM_NOTSET);

   setEntry(incDistA2DenumSAA, PARAM_NOTSET);

   setEntry(iterationLimitPH, ALPS_INT_MAX);

   setEntry(nodeLimitPHSubprob, ALPS_INT_MAX);

   //-------------------------------------------------------------
   // Double Parameters
   //-------------------------------------------------------------

   setEntry(boundCutTimeLim, 3600);

   setEntry(optimalRelGapLimitPHSubprob, 1.0e-4);
   
   //-------------------------------------------------------------
   // String Parameters
   //-------------------------------------------------------------
  
   setEntry(auxiliaryInfoFile, "");

   setEntry(auxiliaryTimFile, "");

   setEntry(auxiliaryStoFile, "");

   setEntry(feasCheckSolver, "SYMPHONY");

   setEntry(inputFormat, "indexBased");

   setEntry(stochasticityType, "deterministic");
}

//#############################################################################
