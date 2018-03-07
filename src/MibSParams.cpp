/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*          Sahar Tahernajad, Lehigh University                              */
/*                                                                           */
/* Copyright (C) 2007-2017 Lehigh University, Scott DeNegre, and Ted Ralphs. */
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
   
   keys_.push_back(make_pair(std::string("MibS_turnOffOtherCuts"),
			     AlpsParameter(AlpsBoolPar, turnOffOtherCuts)));

   keys_.push_back(make_pair(std::string("MibS_printProblemInfo"),
			     AlpsParameter(AlpsBoolPar, printProblemInfo)));

   keys_.push_back(make_pair(std::string("MibS_allowRemoveCut"),
			     AlpsParameter(AlpsBoolPar, allowRemoveCut)));

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

   keys_.push_back(make_pair(std::string("MibS_intersectionCutType"),
			     AlpsParameter(AlpsIntPar, intersectionCutType)));

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

   //--------------------------------------------------------
   // String Parameters.
   //--------------------------------------------------------
   
   keys_.push_back(make_pair(std::string("MibS_auxiliaryInfoFile"),
			     AlpsParameter(AlpsStringPar, auxiliaryInfoFile)));

   keys_.push_back(make_pair(std::string("MibS_feasCheckSolver"),
			     AlpsParameter(AlpsStringPar, feasCheckSolver)));

   keys_.push_back(make_pair(std::string("MibS_inputFormat"),
			     AlpsParameter(AlpsStringPar, inputFormat)));

}

//#############################################################################

void 
MibSParams::setDefaultEntries() {
   
   //-------------------------------------------------------------
   // Bool Parameters.
   //-------------------------------------------------------------

   setEntry(useValFuncCut, false);
     
   setEntry(useBoundCut, false);
   
   setEntry(boundCutOptimal, false);
   
   setEntry(boundCutRelaxUpper, true);

   setEntry(useIpBound, false);

   setEntry(isBilevelBranchProb, false);

   setEntry(warmStartLL, false);

   setEntry(doDualFixing, false);

   setEntry(turnOffOtherCuts, false);

   setEntry(printProblemInfo, true);

   setEntry(allowRemoveCut, false);

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

   setEntry(intersectionCutType, MibSIntersectionCutTypeNotSet);

   setEntry(bilevelFreeSetTypeIC, MibSBilevelFreeSetTypeICNotSet);

   setEntry(solveSecondLevelWhenXYVarsInt, PARAM_NOTSET);

   setEntry(solveSecondLevelWhenXVarsInt, PARAM_NOTSET);

   setEntry(solveSecondLevelWhenLVarsInt, PARAM_NOTSET);

   setEntry(solveSecondLevelWhenLVarsFixed, PARAM_NOTSET);

   setEntry(computeBestUBWhenXVarsInt, PARAM_NOTSET);

   setEntry(computeBestUBWhenLVarsInt, PARAM_NOTSET);

   setEntry(computeBestUBWhenLVarsFixed, PARAM_NOTSET);

   setEntry(useLinkingSolutionPool, PARAM_NOTSET);

   //-------------------------------------------------------------
   // Double Parameters
   //-------------------------------------------------------------
   
   //-------------------------------------------------------------
   // String Parameters
   //-------------------------------------------------------------
  
   setEntry(auxiliaryInfoFile, "");

   setEntry(feasCheckSolver, "SYMPHONY");

   setEntry(inputFormat, "indexBased");
}

//#############################################################################
