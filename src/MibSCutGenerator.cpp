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

#include "OsiCbcSolverInterface.hpp"

#include "MibSCutGenerator.hpp"
#include "MibSParams.hpp"
#include "MibSTreeNode.hpp"
#include "MibSSolution.hpp"
#include "MibSConstants.hpp"
#include "MibSConfig.hpp"

#include "BlisConGenerator.h"
#include "BlisConstraint.h"
#include "BlisHelp.h"
#include "BlisVariable.h"

#ifdef COIN_HAS_SYMPHONY
#include "OsiSymSolverInterface.hpp"
#include "symphony.h"
#endif

#ifdef COIN_HAS_CPLEX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

//#ifdef COIN_HAS_SYMPHONY
//#include "OsiSymSolverInterface.hpp"
//#endif
//to run in parallel, uncomment this (should be modified later)
//#define _OPENMPMIBS 1
#ifdef _OPENMPMIBS
#include "omp.h"
#endif 

//#############################################################################
MibSCutGenerator::MibSCutGenerator(MibSModel *mibs)
{
  localModel_ = mibs;
  setName("MIBS");
  setModel(localModel_);
  auxCount_ = 0;
  upper_ = 0.0;
  maximalCutCount_ = 0;
  numCalledBoundCut_ = -1;
  isBigMIncObjSet_ = false;
  bigMIncObj_ = 0.0;
  watermelonICSolver_ = 0;
}

//#############################################################################
MibSCutGenerator::~MibSCutGenerator()
{
    if(watermelonICSolver_) delete watermelonICSolver_;

}

//#############################################################################
int
MibSCutGenerator::bilevelFeasCut1(BcpsConstraintPool &conPool)
{

    /** Add bilevel feasibility cut **/

    bool allowRemoveCut(localModel_->MibSPar_->entry
			(MibSParams::allowRemoveCut));

    bool useNewPureIntCut(localModel_->MibSPar_->entry
			  (MibSParams::useNewPureIntCut));

    MibSTreeNode * node = static_cast<MibSTreeNode *>(localModel_->activeNode_);

    int newPureIntCutDepthLb(localModel_->MibSPar_->entry
			     (MibSParams::newPureIntCutDepthLb));

    int newPureIntCutDepthUb(localModel_->MibSPar_->entry
			     (MibSParams::newPureIntCutDepthUb));

    OsiSolverInterface * solver = localModel_->solver();

    const int numCols = solver->getNumCols();
    const int numRows = solver->getNumRows();
    const CoinPackedMatrix * matrix = solver->getMatrixByCol();
    MibSBilevel *bS = localModel_->bS_;
    int i(0), j(0);
    int numCuts(0);
    double etol(localModel_->etol_);
    const char * rowsense = solver->getRowSense();

    if(!bS->isIntegral_) return 0;

    bS->isIntegral_ = false;

    /** Find binding constraints at current solution **/
    int * binding = getBindingCons();

    if(!binding){
	std::cout << "No binding constraints set." << std::endl;
	abort();
    }

    double upper(getCutUpperBound());
    //double upper(- getCutUpperBound());
    double cutlb(- solver->getInfinity());
    double cutub(upper - 1.0);

    //double EPS(0.01);
    //double cutub(upper - EPS);

    std::vector<int> indexList;
    std::vector<double> valsList;
    double * tempVals = new double[numCols];
    CoinFillN(tempVals, numCols, 0.0);
    double value(0.0);
    int mult(0);
    double *leftHandSide = NULL;
    if(useNewPureIntCut == true){
	leftHandSide = new double[numCols];
	CoinFillN(leftHandSide, numCols, 0.0);
    }

    for(i = 0; i < numCols; i++){
	for(j = 0; j < numRows; j++){
	    value = matrix->getCoefficient(j, i);
	    switch(rowsense[j]){
	    case 'L':
		mult = 1;
		break;
	    case 'G':
		mult = -1;
		break;
	    case 'E':
		std::cout
		    << "MibS cannot currently handle equality constraints." << std::endl;
		abort();
		break;
	    case 'R':
		std::cout
		    << "MibS cannot currently handle range constraints." << std::endl;
		abort();
		break;
	    }
	    tempVals[i] += binding[j] * value * mult;
	}

	tempVals[i] += binding[numRows + i]; // variable upper bound
	tempVals[i] += - binding[numRows + numCols + i]; // variable lower bound

	if((tempVals[i] > etol) || (tempVals[i] < - etol)){
	    if(leftHandSide != NULL){
		leftHandSide[i] = -1 * tempVals[i];
	    }
	    indexList.push_back(i);
	    valsList.push_back(tempVals[i]);
	}
    }

    bool shouldUseNewPureIntCut(false);
    int depth(0);
    bool shouldPrune(false);

    if(useNewPureIntCut == true){
	depth = node->getDepth();
	if((newPureIntCutDepthLb >= 0) && (newPureIntCutDepthUb >= 0)){
	    if((depth >= newPureIntCutDepthLb) && (depth <= newPureIntCutDepthUb)){
		shouldUseNewPureIntCut = true;
	    }
	}
	else if(newPureIntCutDepthLb >= 0){
	    if(depth >= newPureIntCutDepthLb){
		shouldUseNewPureIntCut = true;
	    }
	}
	else if(newPureIntCutDepthUb >= 0){
	    if(depth <= newPureIntCutDepthUb){
		shouldUseNewPureIntCut = true;
	    }
	}
	else{
	    shouldUseNewPureIntCut = true;
	}
    }

    if(shouldUseNewPureIntCut == true){
	double tmpCutUb(0.0);
	boundCuts(conPool, leftHandSide, tmpCutUb, shouldPrune);

	if(shouldPrune == true){
	    //localModel_->counterGood_ ++;
	    localModel_->bS_->shouldPrune_ = true;
	}
	else{
	    tmpCutUb = floor(tmpCutUb + etol);

	    if(tmpCutUb < cutub){
		//localModel_->counterGood_ ++;
		cutub = tmpCutUb;
	    }
	    else{
		//localModel_->counterBad_ ++;
	    }
	}
    }

    assert(indexList.size() == valsList.size());

    if(shouldPrune == false){
	numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, allowRemoveCut);
    }

    delete [] binding;
    delete [] tempVals;
    if(leftHandSide != NULL){
	delete [] leftHandSide;
    }

    return numCuts;

}

//#############################################################################
int
MibSCutGenerator::bilevelFeasCut2(BcpsConstraintPool &conPool)
{

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  int i(0);
  int numCuts(0);
  double etol(localModel_->etol_);



  double cutlb(- solver->getInfinity());
  double cutub(solver->getInfinity());

  std::vector<int> indexList;
  std::vector<double> valsList;

  /* 
     returns a double with the values [alpha | beta | gamma]
     indexed 0..n1-1, n1..n1+n2-1, n1+n2
  */

  double * cutVals = findDeepestLandPCut_ValFunc();  
  //double * cutVals = findDeepestLandPCut1();  
  double val(0.0);

  /* find the nonzero entries and convert to <= constraint */

  if(0)
     std::cout << "cut coeffs" << std::endl;

  for(i = 0; i < numCols; i++){
    val = cutVals[i];
    if(0)
      std::cout << i << ": " << -val << std::endl;
    if((val > etol) || (val < - etol)){
      // entry is nonzero
      indexList.push_back(i);
      valsList.push_back(- val);
    }
  }

  if(0)
    std::cout << "cutub: " << -cutVals[numCols] << std::endl;

  cutub = - cutVals[numCols]; //gamma

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);

  delete [] cutVals;
  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::feasibilityCuts(BcpsConstraintPool &conPool)
{

  /** Add bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIBS Cuts." << std::endl;

  //bool simpleCutOnly 
  //= localModel_->MibSPar_->entry(MibSParams::simpleCutOnly);
    
   bool usePureIntegerCut 
	= localModel_->MibSPar_->entry(MibSParams::usePureIntegerCut);

   bool useValFuncCut 
    = localModel_->MibSPar_->entry(MibSParams::useValFuncCut);

  if(usePureIntegerCut && !useValFuncCut){
    return bilevelFeasCut1(conPool) ? true : false;
  }
  else if(!usePureIntegerCut && useValFuncCut){
    return bilevelFeasCut2(conPool) ? true : false;
  }
  else if(usePureIntegerCut && useValFuncCut){
     return ((bilevelFeasCut1(conPool) ? true : false) || 
	     (bilevelFeasCut2(conPool) ? true : false));
  }
  else{
    //std::cout << "No MIBS Cuts generated" << std::endl;
    return 0;
  }

}

//#############################################################################
int
MibSCutGenerator::intersectionCuts(BcpsConstraintPool &conPool,
				   double *optLowerSolution, MibSIntersectionCutType ICType)
{
    
    //MibSIntersectionCutType ICType = static_cast<MibSIntersectionCutType>
    //	    (localModel_->MibSPar_->entry(MibSParams::intersectionCutType));

    //std::cout << "ICType = " << ICType << std::endl;

	MibSBilevelFreeSetTypeIC bilevelFreeSetType =
	    static_cast<MibSBilevelFreeSetTypeIC>
	    (localModel_->MibSPar_->entry(MibSParams::bilevelFreeSetTypeIC));
	
	int useLinkingSolutionPool(localModel_->MibSPar_->entry
		       (MibSParams::useLinkingSolutionPool));

	bool allowRemoveCut(localModel_->MibSPar_->entry(MibSParams::allowRemoveCut));

	OsiSolverInterface * solver = localModel_->solver();

	const double * sol = solver->getColSolution();

	const CoinPackedMatrix * matrix = solver->getMatrixByRow();

	CoinWarmStartBasis * ws =
	    dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());

	MibSBilevel *bS = localModel_->bS_;
	    
	int numStruct, numArtf;
	int i,j,index;
	int numCuts(0);
	bool intersectionFound(true);
	
	numStruct = ws->getNumStructural();
	
	numArtf = ws->getNumArtificial();
	
	double etol(localModel_->etol_);
	double value(0.0);
	int numBasicOrig(numArtf);
	int numNonBasicOrig(numStruct);
	int numBasic(0);
	int numNonBasic(0);
	int newBasic(0);
	int numCols(numStruct+numArtf);
	int numRows(numArtf);
	int lCols(localModel_->getLowerDim());
	int lRows(localModel_->getLowerRowNum());

	bool isTimeLimReached(false);

	const double * colLb = solver->getColLower();
	const double * colUb = solver->getColUpper();

	double *lpSol = new double[numStruct];
	memcpy(lpSol, sol, sizeof(double) * numStruct);

	for(i = 0; i < numStruct; i++){
	    value = lpSol[i];
	    if(fabs(floor(value + 0.5) - value) <= etol){
		lpSol[i] = floor(value + 0.5);
	    }
	}
		
	int * rstat = new int[numArtf]();
	int * cstat = new int[numStruct]();

	for(i = 0; i < numArtf; i++){
	    switch (ws->getArtifStatus(i)){
	    case CoinWarmStartBasis::basic:
		rstat[i] = VAR_BASIC;
		break;
	    case CoinWarmStartBasis::atLowerBound:
		rstat[i] = VAR_AT_LB;
		break;
	    case CoinWarmStartBasis::atUpperBound:
		rstat[i] = VAR_AT_UB;
		break;
	    case CoinWarmStartBasis::isFree:
		rstat[i] = VAR_FREE;
		break;
	    }
	}

	for(i = 0; i < numStruct; i++){
	    switch (ws->getStructStatus(i)){
	    case CoinWarmStartBasis::basic:
		cstat[i] = VAR_BASIC;
		break;
	    case CoinWarmStartBasis::atLowerBound:
		cstat[i] = VAR_AT_LB;
		if((colLb[i] > etol) || (fabs(colLb[i] - colUb[i]) < etol)){
		    newBasic ++;
		}
		break;
	    case CoinWarmStartBasis::atUpperBound:
		cstat[i] = VAR_AT_UB;
		newBasic ++;
		break;
	    case CoinWarmStartBasis::isFree:
		cstat[i] = VAR_FREE;
		break;
	    }
	}
	
	numNonBasic = numNonBasicOrig;
	numBasic = numBasicOrig + newBasic;
	
	solver->enableSimplexInterface(1);
	
	int * basicIndexTmp = new int[numBasicOrig];
	
	solver->getBasics(basicIndexTmp);
	
	int * basicIndex = new int[numBasic];
	
	for(i = 0 ; i < numBasicOrig; i++){
	    basicIndex[i] = basicIndexTmp[i];
	}

	int cnt(numBasicOrig);
	const char * rowsense = solver->getRowSense();
	
	for(i = 0; i < numStruct; i++){
	    if((cstat[i] == VAR_AT_LB) &&
	       ((colLb[i] > etol) || (fabs(colLb[i] - colUb[i]) < etol))){
		basicIndex[cnt] = i;
		cnt ++;
	    }
	    else if(cstat[i] == VAR_AT_UB){
		basicIndex[cnt] = i;
		cnt ++;
	    }
	}
	
	int * isBasic = new int[numCols]();
	
	for(i = 0; i < numBasic; i++){
	    index = basicIndex[i];
	    isBasic[index] = 1;
	}
	
	cnt = 0;

	//get Optimal tableau
	double ** extRay = new double*[numNonBasic];
	for(i = 0; i < numNonBasic; ++i){
	    extRay[i] = new double[numStruct];
	}

	double * coef = new double[numRows];

	int mult(1);

	for(i = 0; i < numCols; i++){
	    if(isBasic[i] < 1){
		mult = 1;
		for(j = 0; j < numStruct; j++){
		    extRay[cnt][j] = 0;
		}
		solver->getBInvACol(i,coef);
		if((i >= numStruct) && (i-numStruct < numArtf)){
		    if(rowsense[i-numStruct] == 'G'){
			mult = -1;
		    }
		}
		for(j = 0; j < numBasicOrig; j++){
		    if(basicIndex[j] < numStruct){
			extRay[cnt][basicIndex[j]] = -1 * mult * coef[j];
		    }
		}
		if(i < numStruct){
		    extRay[cnt][i] = 1;
		}
		cnt++;
	    }
	}

	for(i = 0; i < numStruct; i++){
	    if(cstat[i] == VAR_AT_LB){
		if(fabs(colLb[i] - colUb[i]) < etol){
		    for(j = 0; j < numStruct; j++){
			extRay[cnt][j] = 0;
		    }
		    solver->getBInvACol(i,coef);
		    if(colLb[i] > etol){
			for(j = 0; j < numBasicOrig; j++){
			    if(basicIndex[j] < numStruct){
				extRay[cnt][basicIndex[j]] = -1 * coef[j];
			    }
			}
			extRay[cnt][i] = 1;
			cnt++;
		    }
		    else{
			for(j = 0; j < numBasicOrig; j++){
			    if(basicIndex[j] < numStruct){
				extRay[cnt][basicIndex[j]] = coef[j];
			    }
			}
			extRay[cnt][i] = -1;
			cnt++;
		    }
		}
		else if(colLb[i] > etol){
		    for(j = 0; j < numStruct; j++){
			extRay[cnt][j] = 0;
		    }
		    solver->getBInvACol(i,coef);
		    for(j = 0; j < numBasicOrig; j++){
			if(basicIndex[j] < numStruct){
			    extRay[cnt][basicIndex[j]] = -1 * coef[j];
			}
		    }
		    extRay[cnt][i] = 1;
		    cnt++;
		}
	    }
	    else if(cstat[i] == VAR_AT_UB){
		for(j = 0; j < numStruct; j++){
		    extRay[cnt][j] = 0;
		}
		solver->getBInvACol(i,coef);
		for(j = 0; j < numBasicOrig; j++){
		    if(basicIndex[j] < numStruct){
			extRay[cnt][basicIndex[j]] = coef[j];
		    }
		}
		extRay[cnt][i] = -1;
		cnt++;
	    }
	}
	
	std::vector<double> alpha(numNonBasic);

	bool shouldFindBestSol(true);
	if(((useLinkingSolutionPool != PARAM_ON) &&
	    ((bS->isUBSolved_ == true) || ((bS->isLowerSolved_ == true) &&
					   (bS->isProvenOptimal_ == false)))) ||
	   ((useLinkingSolutionPool == PARAM_ON) &&
	    ((bS->tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsInfeasible) ||
	     (bS->tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved)))){
	    shouldFindBestSol = false;
	}
	switch(ICType){
	case MibSIntersectionCutTypeIC:
	    {
		if(bilevelFreeSetType == MibSBilevelFreeSetTypeICWithNewLLSol){
		    //it shows which lower level rows are not
		    //necessary to be included in the convex set
		    double *uselessIneqs = new double[lRows];
		    double *lowerLevelSol = new double[lCols];
		    CoinZeroN(uselessIneqs, lRows);
		    CoinZeroN(lowerLevelSol, lCols);
		    isTimeLimReached = false;
		    findLowerLevelSol(uselessIneqs, lowerLevelSol, sol, isTimeLimReached);
		    if(isTimeLimReached == true){
			goto TREM_INTERSECTIONCUT;
		    }
		    intersectionFound = getAlphaIC(extRay, uselessIneqs, lowerLevelSol,
						   numStruct, numNonBasic, sol, alpha);

		    if(intersectionFound == false){
			localModel_->bS_->shouldPrune_ = true;
		    }
		    
		    delete [] uselessIneqs;
		    delete [] lowerLevelSol;
		}
		else{
		    intersectionFound = getAlphaIC(extRay, NULL, optLowerSolution,
						   numStruct, numNonBasic, sol, alpha);
		    if(intersectionFound == false){
			localModel_->bS_->shouldPrune_ = true;
		    }
		}
		break;
	    }
        case MibSIntersectionCutTypeWatermelon:
	    {
		double *uselessIneqs = new double[lRows + 2 * lCols];
	        double *lowerLevelSol = new double[lCols];
	        CoinZeroN(uselessIneqs, lRows);
	        CoinZeroN(lowerLevelSol, lCols);
		isTimeLimReached = false;
	        findLowerLevelSolWatermelonIC(uselessIneqs, lowerLevelSol, lpSol, isTimeLimReached);
		if(isTimeLimReached == true){
		    goto TREM_INTERSECTIONCUT;
		}
		intersectionFound = getAlphaWatermelonIC(extRay, uselessIneqs, lowerLevelSol,
							 numStruct, numNonBasic, lpSol, alpha);
	        delete [] uselessIneqs;
	        delete [] lowerLevelSol;
	        break;
	    }
	case MibSIntersectionCutTypeHypercubeIC:
	    if(shouldFindBestSol == true){
		isTimeLimReached = false;
		storeBestSolHypercubeIC(sol, bS->objValVec_, isTimeLimReached);
		if(isTimeLimReached == true){
		    goto TREM_INTERSECTIONCUT;
		}
	    }
	    getAlphaHypercubeIC(extRay, numStruct, numNonBasic, alpha);
	    break;
	case MibSIntersectionCutTypeTenderIC:
	    if(shouldFindBestSol == true){
		storeBestSolTenderIC(sol, bS->objValVec_[0]);
	    }
	    getAlphaTenderIC(extRay, numNonBasic, alpha);
	    break;
	case MibSIntersectionCutTypeHybridIC:
	    if(shouldFindBestSol == true){
		storeBestSolTenderIC(sol, bS->objValVec_[0]);
	    }
	    getAlphaHypercubeIC(extRay, numStruct, numNonBasic, alpha);
	    break;
	}
	
	if(intersectionFound){
	    cnt = 0;
	    mult = 0;

	    double tmp(0.0);
	    double cutUb(-1.0);
	    double cutLb(-1 * solver->getInfinity());
	    double rowRhs;
	    int rowIndex;
	
	    double * tmpValsList = new double[numStruct];
	    CoinFillN(tmpValsList, numStruct, 0.0);
	
	    std::vector<int> indexList;
	    std::vector<double> valsList;
		    
	    for(i = 0; i < numCols; i++){
		if(isBasic[i] < 1){
		    if(alpha[cnt] >= 0){
			if (i < numStruct){
			    if((solver->isInteger(i)) && ((1/alpha[cnt]) - 1 > etol)){
				alpha[cnt] = 1.0;
			    }
			    tmpValsList[i] += (1/alpha[cnt]);
			}
			else{
			    rowIndex = i - numStruct;
			    switch(rowsense[rowIndex]){
			    case 'L':
				mult = 1;
				break;
			    case 'G':
				mult = -1;
				break;
			    case 'E':
				std::cout
				    << "MibS cannot currently handle equality constraints." << std::endl;
				abort();
				break;
			    case 'R':
				std::cout
				    << "MibS cannot currently handle range constraints." << std::endl;
				abort();
				break;
			    }
			    rowRhs = mult * solver->getRightHandSide()[rowIndex];
			    cutUb += rowRhs/alpha[cnt];
			    for(j = 0; j < numStruct; j++){
				tmp = mult * matrix->getCoefficient(rowIndex,j);
			        tmpValsList[j] += -1 * (tmp/alpha[cnt]);
			    }
			}
		    }
		    cnt ++;
		}
	    }

	    for (i = 0; i < numStruct; i++){
		if(cstat[i] == VAR_AT_LB){
		    if(fabs(colLb[i] - colUb[i]) < etol){
			if(colLb[i] > etol){
			    if(alpha[cnt] >= 0){
				cutUb += -1 * (colLb[i]/alpha[cnt]);
				tmpValsList[i] += 1/alpha[cnt];
			    }
			}
			else{
			    if(alpha[cnt] >= 0){
				tmpValsList[i] += -1 * (1/alpha[cnt]);
			    }
			}
			cnt ++;
		    }
		    else if(colLb[i] > etol){
			if(alpha[cnt] >= 0){
			    cutUb += -1 * (colLb[i]/alpha[cnt]);
			    tmpValsList[i] += 1/alpha[cnt];
			}
			cnt ++;
		    }
		}
		else if(cstat[i] == VAR_AT_UB){
		    if(alpha[cnt] >= 0){
			cutUb += colUb[i]/alpha[cnt];
			tmpValsList[i] += -1 * (1/alpha[cnt]);
		    }
		    cnt ++;
		}
	    }

	    for(i = 0; i < numStruct; i++){
		if(fabs(tmpValsList[i]) >= etol){
		    indexList.push_back(i);
		    valsList.push_back(-1*tmpValsList[i]);
		}
	    }

	    if(valsList.size() == 0){
		throw CoinError("An invalid intersection cut is generated (length = 0)",
				"intersectionCuts", "MibSCutGenerator");
	    }

	    numCuts += addCut(conPool, cutLb, cutUb, indexList, valsList,
			      allowRemoveCut);

	    if(numCuts == 0){
		MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
		    (localModel_->MibSPar_->entry(MibSParams::branchStrategy));
		if(branchPar == MibSBranchingStrategyFractional){
		    double scaleConFactor = localModel_->BlisPar()->entry
			(BlisParams::scaleConFactor);
		    double maxElem = 0.0;
		    double minElem = ALPS_DBL_MAX;
		    double scaleFactor;
		    for(i = 0; i < valsList.size(); i++){
			if(fabs(valsList[i]) > maxElem){
			    maxElem = fabs(valsList[i]);
			}
			if(fabs(valsList[i]) < minElem){
			    minElem = fabs(valsList[i]);
			}
			if(minElem != 0.0){
			    scaleFactor = maxElem/minElem;
			}
			if(scaleFactor > scaleConFactor){
			    std::cout << "scaleFactor = " << scaleFactor << std::endl;
			    throw CoinError("Bad scaling", "intersectionCuts",
					    "MibSCutGenerator");
			}
		    }
		}
	    }

	    
	    delete [] tmpValsList;
	    indexList.clear();
	    valsList.clear();
	}

 TREM_INTERSECTIONCUT:
	solver->disableSimplexInterface();

	delete ws;
	delete [] lpSol;
	delete[] rstat;
	delete[] cstat;
	delete[] basicIndexTmp;
	delete[] basicIndex;
	delete[] isBasic;
	delete[] coef;
	for(i = 0; i < numNonBasic; ++i){
	    delete[] extRay[i];
	}
	delete[] extRay;
	
	return numCuts;
}

//#############################################################################
void
MibSCutGenerator::findLowerLevelSol(double *uselessIneqs, double *lowerLevelSol,
				    const double *sol, bool &isTimeLimReached)
{
    
    std::string feasCheckSolver(localModel_->MibSPar_->entry
				(MibSParams::feasCheckSolver));
    int maxThreadsLL(localModel_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(localModel_->MibSPar_->entry
		    (MibSParams::whichCutsLL));

    double timeLimit(localModel_->AlpsPar()->entry(AlpsParams::timeLimit));
    double remainingTime(0.0);

    bool getA2Matrix(false), getG2Matrix(false);
    OsiSolverInterface * oSolver = localModel_->solver();
    double infinity(oSolver->getInfinity());
    int i, j;
    int index(0), cntA2(0), cntG2(0), cntInt(0);
    double coef(0.0), lObjVal(0.0), value(0.0);
    int numCols(localModel_->getNumCols());
    int uCols(localModel_->getUpperDim());
    int lCols(localModel_->getLowerDim());
    int lRows(localModel_->getLowerRowNum());
    int numBinCols(lRows);
    int newNumCols(lCols + numBinCols);
    int newNumRows(lRows + 1);
    double lObjSense(localModel_->getLowerObjSense());
    double *lObjCoeff(localModel_->getLowerObjCoeffs());
    int *lRowInd(localModel_->getLowerRowInd());
    int *uColInd(localModel_->getUpperColInd());
    int *lColInd(localModel_->getLowerColInd());
    double *origRowLb(localModel_->getOrigRowLb());
    double *origRowUb(localModel_->getOrigRowUb());
    double *origColLb(localModel_->getOrigColLb());
    double *origColUb(localModel_->getOrigColUb());
    const double *colLb(oSolver->getColLower());
    const double *colUb(oSolver->getColUpper());
    char *rowSense = localModel_->getOrigRowSense();

    CoinPackedVector row;
    CoinPackedVector addedRow;

    //stores A^2*x (x is the optimal first-level solution of relaxation)
    double *multA2XOpt = new double[lRows];
    //stores A^2*lb(x)
    double *multA2XLb = new double[lRows];
    //stores A^2*ub(x)
    double *multA2XUb =new double[lRows];
    //store max(A^2*lb(x), A^2*ub(x))
    double *multA2XMax = new double[lRows];

    double *optUpperSol = new double[uCols];
    double *upperColLb = new double[uCols];
    double *upperColUb = new double[uCols];
    
    CoinPackedMatrix * matrixA2 = 0;
    CoinPackedMatrix * matrixG2 = 0;
    CoinPackedMatrix * newMatrix = new CoinPackedMatrix(false, 0, 0);
    newMatrix->setDimensions(0, newNumCols);

    double *newRowLb = new double[newNumRows];
    double *newRowUb = new double[newNumRows];
    CoinZeroN(newRowLb, newNumRows);
    CoinZeroN(newRowUb, newNumRows);
    double *newColLb = new double[newNumCols];
    double *newColUb = new double[newNumCols];
    CoinZeroN(newColLb, newNumCols);
    CoinZeroN(newColUb, newNumCols);
    double *newObjCoeff = new double[newNumCols];
    CoinZeroN(newObjCoeff, newNumCols);
    int *integerVars = new int[newNumCols];
    CoinZeroN(integerVars, newNumCols);

    if(!localModel_->getA2Matrix()){
	getA2Matrix = true;
    }

    if(!localModel_->getG2Matrix()){
	getG2Matrix = true;
    }

    if((getA2Matrix) || (getG2Matrix)){
	getLowerMatrices(false, getA2Matrix, getG2Matrix);
    }

    //extracting optimal first-level solution of the relaxation problem and
    // the original bounds
    for(i = 0; i < uCols; i++){
	index = uColInd[i];
	optUpperSol[i] = sol[index];
	upperColLb[i] = colLb[index];
	upperColUb[i] = colUb[index];
    }

    matrixA2 = localModel_->getA2Matrix();
    matrixG2 = localModel_->getG2Matrix();

    matrixA2->times(optUpperSol, multA2XOpt);
    matrixA2->times(upperColLb, multA2XLb);
    matrixA2->times(upperColUb, multA2XUb);
    for(i = 0; i < lRows; i++){
	multA2XMax[i] = CoinMax(multA2XLb[i], multA2XUb[i]);
    }

    //filling matrix of coefficients
    for(i = 0; i < lCols; i++){
	addedRow.insert(i, lObjCoeff[i] * lObjSense);
    }

    newMatrix->appendRow(addedRow);
    addedRow.clear();

    for(i = 0; i < lRows; i++){
	row = matrixG2->getVector(i);
	addedRow.append(row);
	//since linking vars are int, we are sure that multA2XOpt
	//is integer.
	value = floor(multA2XOpt[i] + 0.5);
	addedRow.insert(i + lCols, value - multA2XMax[i]);
	newMatrix->appendRow(addedRow);
	addedRow.clear();
    }

    //filling row bounds
    CoinFillN(newRowLb, newNumRows, -1 * infinity);
    for(i = 0; i < lCols; i++){
	index = lColInd[i];
	lObjVal += lObjCoeff[i] * sol[index];
    }
    newRowUb[0] = lObjVal * lObjSense - 1;

    for(i = 0; i < lRows; i++){
	index = lRowInd[i];
	if(rowSense[index] == 'L'){
	    newRowUb[i + 1] = origRowUb[index] - multA2XMax[i];
	}
	else if(rowSense[index] == 'G'){
	    newRowUb[i + 1] = -1 * origRowLb[index] - multA2XMax[i];
	}
	else{
	    assert(0);
	}
    }

    //filling col bounds
    for(i = 0; i < lCols; i++){
	newColLb[i] = origColLb[lColInd[i]];
	newColUb[i] = origColUb[lColInd[i]];
    }

    CoinFillN(newColLb + lCols, numBinCols, 0.0);
    CoinFillN(newColUb + lCols, numBinCols, 1.0);

    //filling objective coefficients
    CoinFillN(newObjCoeff + lCols, numBinCols, 1.0);

    //filling integer vars
    for(i = 0; i < lCols; i++){
	index = lColInd[i];
	if(oSolver->isInteger(index)){
	    integerVars[cntInt] = i;
	    cntInt ++;
	}
    }
    CoinIotaN(&integerVars[cntInt], numBinCols, lCols);
    cntInt += numBinCols;

    OsiSolverInterface * nSolver;

    if(feasCheckSolver == "Cbc"){
        nSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"findLowerLevelSol", "MibSCutGenerator");
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"findLowerLevelSol", "MibsBilevel");
#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"findLowerLevelSol", "MibsBilevel");
    }

    nSolver->loadProblem(*newMatrix, newColLb, newColUb,
			 newObjCoeff, newRowLb, newRowUb);

    nSolver->setInteger(integerVars, cntInt);
    nSolver->setObjSense(1); //1 min; -1 max
    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getTime();

    if(remainingTime <= localModel_->etol_){
	isTimeLimReached = true;
	localModel_->bS_->shouldPrune_ = true;
    }
    else{
	remainingTime = CoinMax(remainingTime, 0.00);

	if(feasCheckSolver == "Cbc"){
	    dynamic_cast<OsiCbcSolverInterface *>
		(nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
	}else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	    sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		(nSolver)->getSymphonyEnvironment();
	    sym_set_dbl_param(env, "time_limit", remainingTime);
	    sym_set_int_param(env, "do_primal_heuristic", FALSE);
            sym_set_int_param(env, "verbosity", -2);
            sym_set_int_param(env, "prep_level", -1);
            sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
            sym_set_int_param(env, "tighten_root_bounds", FALSE);
            sym_set_int_param(env, "max_sp_size", 100);
            sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	    if (whichCutsLL == 0){
		sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	    }else{
		sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	    }
	    if(whichCutsLL == 1){
		sym_set_int_param(env, "generate_cgl_knapsack_cuts",
				  DO_NOT_GENERATE);
	        sym_set_int_param(env, "generate_cgl_probing_cuts",
				  DO_NOT_GENERATE);
	        sym_set_int_param(env, "generate_cgl_clique_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_twomir_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_flowcover_cuts",
				  DO_NOT_GENERATE);
	    }
#endif
	}else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    nSolver->setHintParam(OsiDoReducePrint);
	    nSolver->messageHandler()->setLogLevel(0);
	    CPXENVptr cpxEnv =
		dynamic_cast<OsiCpxSolverInterface*>(nSolver)->getEnvironmentPtr();
	    assert(cpxEnv);
	    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
	}

	nSolver->branchAndBound();

        if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					       (dynamic_cast<OsiSymSolverInterface *>
						(nSolver)->getSymphonyEnvironment()))){
	    isTimeLimReached = true;
	    localModel_->bS_->shouldPrune_ = true;
	}
	else if(nSolver->isProvenOptimal()){
	    const double *optSol = nSolver->getColSolution();
	    CoinDisjointCopyN(optSol, lCols, lowerLevelSol);
	    CoinDisjointCopyN(optSol + lCols, numBinCols, uselessIneqs);
	}
	else{//this case should not happen when we use intersection cut for removing
	    //the optimal solution of relaxation which satisfies integrality requirements
	    throw CoinError("The MIP which gives the best lower-level sol, cannot be infeasible!",
			    "findLowerLevelSol", "MibSCutGenerator");
	}
    }

    delete [] multA2XOpt;
    delete [] multA2XLb;
    delete [] multA2XUb;
    delete [] multA2XMax;
    delete [] optUpperSol;
    delete [] upperColLb;
    delete [] upperColUb;
    delete newMatrix;
    delete [] newRowLb;
    delete [] newRowUb;
    delete [] newColLb;
    delete [] newColUb;
    delete [] newObjCoeff;
    delete [] integerVars;
    delete nSolver;
    
}
    
//#############################################################################
bool
MibSCutGenerator::getAlphaIC(double** extRay, double* uselessIneqs,
			     double* lowerSolution, int numStruct, int numNonBasic,
			     const double* lpSol, std::vector<double> &alphaVec)
{
    bool intersectionFound(false);
    int i;
    int colIndex(0), rowIndex(0);
    double alpha(0.0), value(0.0);
    double etol(localModel_->etol_);
    int lRows(localModel_->getLowerRowNum());
    int numRows(lRows + 1);
    int numCols(localModel_->getNumCols());
    int uCols(localModel_->getUpperDim());
    int lCols(localModel_->getLowerDim());
    int *uColInd(localModel_->getUpperColInd());
    int *lColInd(localModel_->getLowerColInd());
    int *lRowInd(localModel_->getLowerRowInd());
    double *origRowLb(localModel_->getOrigRowLb());
    double *origRowUb(localModel_->getOrigRowUb());
    char *rowSense = localModel_->getOrigRowSense();
    double *lObjCoeffs(localModel_->getLowerObjCoeffs());
    double objSense(localModel_->getLowerObjSense());
    bool getA2Matrix(false), getG2Matrix(false);
    
    if(localModel_->getA2Matrix() == NULL){
	getA2Matrix = true;
    }

    if(localModel_->getG2Matrix() == NULL){
	getG2Matrix = true;
    }

    if((getA2Matrix) || (getG2Matrix)){
	getLowerMatrices(false, getA2Matrix, getG2Matrix);
    }
    
    CoinPackedMatrix *matrixA2 = localModel_->getA2Matrix();
    CoinPackedMatrix *matrixG2 = localModel_->getG2Matrix();

    double *ray= new double[numCols];
    CoinZeroN(ray, numCols);
    double *rhs = new double[numRows];
    CoinZeroN(rhs, numRows);
    double *optUpperSol = new double[uCols];
    CoinZeroN(optUpperSol, uCols);
    //stores A^2*x (x is the optimal first-level solution of relaxation)
    double *multA2XOpt = new double[lRows];
    //Stores G^2*lowerSolution
    double *multG2LSol = new double[lRows];

    //extracting optimal first-level solution of the relaxation problem
    for(i = 0; i < uCols; i++){
	colIndex = uColInd[i];
	optUpperSol[i] = lpSol[colIndex];
    }

    matrixA2->times(optUpperSol, multA2XOpt);
    matrixG2->times(lowerSolution, multG2LSol);

    for(i = 0; i < lRows; i++){
	rowIndex = lRowInd[i];
	if(rowSense[rowIndex] == 'L'){
	    value = origRowUb[rowIndex];
	}
	else{
	    value = -1 * origRowLb[rowIndex];
	}
	rhs[i] = value + 1 - multG2LSol[i] - multA2XOpt[i];
    }

    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	rhs[lRows] += objSense * lObjCoeffs[i] * (lpSol[colIndex] - lowerSolution[i]);
    }

    for (i = 0; i < numNonBasic; i++){
	memcpy(ray, extRay[i], sizeof(double) * numCols);
	alpha = solveModelIC(uselessIneqs, ray, rhs, numNonBasic);
	alphaVec[i] = alpha;
	if(alpha > -1 * etol){
	    intersectionFound = true;
	}
    }

    delete [] ray;
    delete [] rhs;
    delete [] optUpperSol;
    delete [] multA2XOpt;
    delete [] multG2LSol;

    return intersectionFound;

}

//#############################################################################
double
MibSCutGenerator::solveModelIC(double *uselessIneqs, double *ray, double *rhs,
			       int numNonBasic)
{
    bool isUnbounded(true);
    int i;
    double value(0.0);
    double alphaUb(localModel_->solver()->getInfinity()), alpha(0.0);
    double etol(localModel_->etol_);
    int uCols(localModel_->getUpperDim());
    int lCols(localModel_->getLowerDim());
    int lRows(localModel_->getLowerRowNum());
    int numRows(lRows + 1);
    int *uColInd(localModel_->getUpperColInd());
    int *lColInd(localModel_->getLowerColInd());
    double *lObjCoeffs(localModel_->getLowerObjCoeffs());
    double objSense(localModel_->getLowerObjSense());
    CoinPackedMatrix *matrixA2 = localModel_->getA2Matrix();
    
    double * coeff = new double[numRows];
    CoinFillN(coeff, numRows, 0.0);
    double *upperRay = new double[uCols];

    //extract the upper component of ray
    for(i = 0; i < uCols; i++){
	upperRay[i] = ray[uColInd[i]];
    }

    matrixA2->times(upperRay, coeff);

    for(i = 0; i < lCols; i++){
	coeff[lRows] += -1 * objSense * lObjCoeffs[i] * ray[lColInd[i]];
    }

    for(i = 0; i < numRows-1; i++){
	if(((uselessIneqs) && (uselessIneqs[i] > 0.5)) ||
	   (!uselessIneqs)){
	    if(coeff[i] > etol){
		value = rhs[i]/coeff[i];
		if(alphaUb - value > etol){
		    alphaUb = value;
		    isUnbounded = false;
		}
	    }
	}
    }

    if(coeff[lRows] > etol){
	value = rhs[lRows]/coeff[lRows];
	if(alphaUb - value > etol){
	    alphaUb = value;
	    isUnbounded = false;
	}
    }

    assert(alphaUb >= 0);

    if(isUnbounded == false){
	alpha = alphaUb;
    }
    else{
	alpha = -1;
    }

    delete [] coeff;
    delete [] upperRay;

    return alpha;
}

//#############################################################################
void
MibSCutGenerator::findLowerLevelSolWatermelonIC(double *uselessIneqs, double *lowerLevelSol,
						double* lpSol, bool &isTimeLimReached)
{
    std::string feasCheckSolver(localModel_->MibSPar_->entry
				(MibSParams::feasCheckSolver));
    int maxThreadsLL(localModel_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(localModel_->MibSPar_->entry
		    (MibSParams::whichCutsLL));
    double timeLimit(localModel_->AlpsPar()->entry(AlpsParams::timeLimit));
    double remainingTime(0.0);
    
    OsiSolverInterface *oSolver = localModel_->solver();
    double infinity(oSolver->getInfinity());
    bool getA2G2Matrix(false), getG2Matrix(false);
    int i, j;
    int rowIndex(0), colIndex(0), numElements(0), pos(0), cntInt(0);
    double rhs(0.0), value(0.0);
    int numCols(localModel_->getNumCols());
    int lCols(localModel_->getLowerDim());
    int lRows(localModel_->getLowerRowNum());
    int numContCols(lRows + 2 * lCols);
    int newNumCols(lCols + numContCols);
    int newNumRows(2 * lRows + 2 * lCols + 1);
    double lObjSense(localModel_->getLowerObjSense());
    double *lObjCoeff(localModel_->getLowerObjCoeffs());
    int *lColInd(localModel_->getLowerColInd());
    int *lRowInd(localModel_->getLowerRowInd());
    double *origRowLb(localModel_->getOrigRowLb());
    double *origRowUb(localModel_->getOrigRowUb());
    double *origColLb(localModel_->getOrigColLb());
    double *origColUb(localModel_->getOrigColUb());
    char *origRowSense(localModel_->getOrigRowSense());
    CoinPackedMatrix origMatrix = *localModel_->getOrigConstCoefMatrix();

    CoinShallowPackedVector origRow;
    CoinPackedVector row;
    CoinPackedVector addedRow;

    //store A^2*x + G^2*y((x,y) is the optimal solution of relaxation)
    double *lCoeffsTimesLpSol = new double[lRows];

    origMatrix.reverseOrdering();
    
    if(!localModel_->getLowerConstCoefMatrix()){
	getA2G2Matrix = true;
    }

    if(!localModel_->getG2Matrix()){
	getG2Matrix = true;
    }

    if((getA2G2Matrix) || (getG2Matrix)){
	getLowerMatrices(getA2G2Matrix, false, getG2Matrix);
    }

    CoinPackedMatrix *matrixA2G2 = localModel_->getLowerConstCoefMatrix();
    
    //filling watermelonICSolver_
    if(!watermelonICSolver_){
	CoinPackedMatrix * newMatrix = new CoinPackedMatrix(false, 0, 0);
	newMatrix->setDimensions(0, newNumCols);
	CoinPackedMatrix *matrixG2 = localModel_->getG2Matrix();
	double *newRowLb = new double[newNumRows];
	double *newRowUb = new double[newNumRows];
	CoinZeroN(newRowLb, newNumRows);
	CoinZeroN(newRowUb, newNumRows);
	double *newColLb = new double[newNumCols];
	double *newColUb = new double[newNumCols];
	CoinZeroN(newColLb, newNumCols);
	CoinZeroN(newColUb, newNumCols);
	double *newObjCoeff = new double[newNumCols];
	CoinZeroN(newObjCoeff, newNumCols);
	int *integerVars = new int[lCols];
	CoinZeroN(integerVars, lCols);
	
	for(i = 0; i < lCols; i++){
	    addedRow.insert(i, lObjCoeff[i] * lObjSense);
	}
	newMatrix->appendRow(addedRow);
        addedRow.clear();

	for(i = 0; i < lRows; i++){
	origRow = matrixG2->getVector(i);
	addedRow = origRow;
	newMatrix->appendRow(addedRow);
	addedRow.insert(i + lCols, -1.00);
	newMatrix->appendRow(addedRow);
	addedRow.clear();
	}

	colIndex = lCols + lRows;
	for(i = 0; i < lCols; i++){
	    addedRow.insert(i, 1.00);
	    addedRow.insert(colIndex + i, -1.00);
	    newMatrix->appendRow(addedRow);
	    addedRow.clear();
	}

	colIndex = 2 * lCols + lRows;
	for(i = 0; i < lCols; i++){
	    addedRow.insert(i, -1.00);
	    addedRow.insert(colIndex + i, -1.00);
	    newMatrix->appendRow(addedRow);
	    addedRow.clear();
	}

	//filling row bounds
	CoinFillN(newRowLb, newNumRows, -1 * infinity);

	//filling col bounds
	CoinFillN(newColUb, newNumCols, infinity);
        for(i = 0; i < lCols; i++){
	    colIndex = lColInd[i];
	    newColLb[i] = origColLb[colIndex] - lpSol[colIndex];
	    newColUb[i] = origColUb[colIndex] - lpSol[colIndex];
	    if(origColUb[colIndex] >  (infinity/10)){
	      newColUb[i] = infinity;
	    }
	}

	//filling objective coeffs
	CoinFillN(newObjCoeff + lCols, numContCols, 1.0);

	//filling integer vars
	for(i = 0; i < lCols; i++){
	    colIndex = lColInd[i];
	    if(oSolver->isInteger(colIndex)){
		integerVars[cntInt] = i;
		cntInt ++;
	    }
	}

	if(feasCheckSolver == "Cbc"){
	    watermelonICSolver_ = new OsiCbcSolverInterface();
	}else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	    watermelonICSolver_ = new OsiSymSolverInterface();
#else
	    throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			    "findLowerLevelSolWatermelonIC", "MibSCutGenerator");
#endif
	}else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    watermelonICSolver_ = new OsiCpxSolverInterface();
#else
	    throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			    "findLowerLevelSolWatermelonIC", "MibSCutGenerator");
#endif
	}else{
	    throw CoinError("Unknown solver chosen",
			    "findLowerLevelSolWatermelonIC", "MibSCutGenerator");
	}

	watermelonICSolver_->loadProblem(*newMatrix, newColLb, newColUb,
					 newObjCoeff, newRowLb, newRowUb);

	watermelonICSolver_->setInteger(integerVars, cntInt);
        watermelonICSolver_->setObjSense(1); //1 min; -1 max
        watermelonICSolver_->setHintParam(OsiDoReducePrint, true, OsiHintDo);

	delete newMatrix;
	delete [] newRowLb;
	delete [] newRowUb;
	delete [] newColLb;
	delete [] newColUb;
	delete [] newObjCoeff;
	delete [] integerVars;
    }

    OsiSolverInterface *nSolver = watermelonICSolver_;

    matrixA2G2->times(lpSol, lCoeffsTimesLpSol);

    //modifying row bounds
    nSolver->setRowUpper(0, -1);
    for(i = 0; i < lRows; i++){
	rowIndex = lRowInd[i];
	if(origRowSense[rowIndex] == 'G'){
	    rhs = -1 * origRowLb[rowIndex];
	}
	else{
	    rhs = origRowUb[rowIndex];
	}
	nSolver->setRowUpper(2 * i + 1, rhs - lCoeffsTimesLpSol[i]);
    }

    //modifying col bounds
    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	value = lpSol[colIndex];
	nSolver->setColLower(i, origColLb[colIndex] - value);
	nSolver->setColUpper(i, origColUb[colIndex] - value);
    }

    remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getTime();
    remainingTime = CoinMax(remainingTime, 5.00);

    if(feasCheckSolver == "Cbc"){
	        dynamic_cast<OsiCbcSolverInterface *>
		    (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	        sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		    (nSolver)->getSymphonyEnvironment();
		sym_set_int_param(env, "use_hot_starts", FALSE);
		sym_set_dbl_param(env, "time_limit", remainingTime);
		sym_set_int_param(env, "do_primal_heuristic", FALSE);
		sym_set_int_param(env, "verbosity", -2);
		sym_set_int_param(env, "prep_level", -1);
		sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
		sym_set_int_param(env, "tighten_root_bounds", FALSE);
		sym_set_int_param(env, "max_sp_size", 100);
		sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
		if (whichCutsLL == 0){
		    sym_set_int_param(env, "generate_cgl_cuts", FALSE);
		}else{
		    sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
		    sym_set_int_param(env, "generate_cgl_twomir_cuts",
				      DO_NOT_GENERATE);
		}
		if(whichCutsLL == 1){
		    sym_set_int_param(env, "generate_cgl_knapsack_cuts",
				      DO_NOT_GENERATE);
		    sym_set_int_param(env, "generate_cgl_probing_cuts",
				      DO_NOT_GENERATE);
		    sym_set_int_param(env, "generate_cgl_clique_cuts",
				      DO_NOT_GENERATE);
		    sym_set_int_param(env, "generate_cgl_twomir_cuts",
				      DO_NOT_GENERATE);
		    sym_set_int_param(env, "generate_cgl_flowcover_cuts",
				      DO_NOT_GENERATE);
		}
#endif

    }else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver->setHintParam(OsiDoReducePrint);
	nSolver->messageHandler()->setLogLevel(0);
	        CPXENVptr cpxEnv =
		    dynamic_cast<OsiCpxSolverInterface*>(nSolver)->getEnvironmentPtr();
		assert(cpxEnv);
		CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
		CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
    }

    nSolver->branchAndBound();

    if(((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					   (dynamic_cast<OsiSymSolverInterface *>
					    (nSolver)->getSymphonyEnvironment())))
       || (timeLimit - localModel_->broker_->subTreeTimer().getTime() < 3)){
	isTimeLimReached = true;
	localModel_->bS_->shouldPrune_ = true;
    }
    else if(nSolver->isProvenOptimal()){
	const double *optSol = nSolver->getColSolution();
	CoinDisjointCopyN(optSol, lCols, lowerLevelSol);
	CoinDisjointCopyN(optSol + lCols, numContCols, uselessIneqs);
    }
    else{//this case should not happen when we use intersection cut for removing
	//the optimal solution of relaxation which satisfies integrality requirements
	std::cout << "iteration = " << localModel_->countIteration_ << std::endl;
	std::cout << "remaining time = " << remainingTime << std::endl;
	std::cout << "current time = " << timeLimit - localModel_->broker_->subTreeTimer().getTime() << std::endl;
	throw CoinError("The MIP which is solved for watermelon IC, cannot be infeasible!",
			"findLowerLevelSolWatermelonIC", "MibSCutGenerator");
    }
    delete [] lCoeffsTimesLpSol;

}

//#############################################################################
bool
MibSCutGenerator::getAlphaWatermelonIC(double** extRay, double *uselessIneqs,
				       double* lowerSolution, int numStruct,
				       int numNonBasic, double* lpSol,
				       std::vector<double> &alphaVec)
{
    int i, j;
    int cntRows(0), rowIndex(0), colIndex(0);
    double value(0.0);
    double etol(localModel_->etol_);
    double infinity(localModel_->solver()->getInfinity());
    double alphaUb(infinity);
    bool isUnbounded(true), intersectionFound(false);
    int numCols(localModel_->getNumCols());
    int lCols(localModel_->getLowerDim());
    int lRows(localModel_->getLowerRowNum());
    int sizeRhs(lRows + 2 * lCols);
    int *lColInd(localModel_->getLowerColInd());
    int *lRowInd(localModel_->getLowerRowInd());
    double *origRowLb(localModel_->getOrigRowLb());
    double *origRowUb(localModel_->getOrigRowUb());
    double *origColLb(localModel_->getOrigColLb());
    double *origColUb(localModel_->getOrigColUb());
    char *origRowSense(localModel_->getOrigRowSense());
    double *rhs = new double[sizeRhs];
    CoinZeroN(rhs, sizeRhs);

    double *ray= new double[numCols];
    CoinZeroN(ray, numCols);

    double * coeff = new double[sizeRhs];
    CoinFillN(coeff, sizeRhs, 0.0);
    
    //store A^2*x + G^2*y((x,y) is the optimal solution of relaxation)
    double *lCoeffsTimesLpSol = new double[lRows];

    double *lCoeffsTimesRay = new double[lRows];
    
    //store G^2*deltaY
    double *G2TimeslSol = new double[lRows];

    CoinPackedMatrix *matrixA2G2 = localModel_->getLowerConstCoefMatrix();
    CoinPackedMatrix *matrixG2 = localModel_->getG2Matrix();

    matrixA2G2->times(lpSol, lCoeffsTimesLpSol);
    matrixG2->times(lowerSolution, G2TimeslSol);

    for(i = 0; i < lRows; i++){
	rowIndex = lRowInd[i];
	if(origRowSense[i] == 'L'){
	    rhs[i] = origRowUb[rowIndex] + 1 - G2TimeslSol[i] -
		lCoeffsTimesLpSol[i];
	}
	else if(origRowSense[i] == 'G'){
	    rhs[i] = -1 * origRowLb[rowIndex] + 1 - G2TimeslSol[i] -
		lCoeffsTimesLpSol[i];
	}
	else{
	    assert(0);
	}
    }

    cntRows = lRows;
    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	value = lpSol[colIndex];
	rhs[cntRows + i] = origColUb[colIndex] + 1 - lowerSolution[i] - value;
	if(origColUb[colIndex] > (infinity/10)){
	  rhs[cntRows + i] = infinity;
	}
	rhs[cntRows + lCols + i] = -1 * origColLb[colIndex] + 1 +
	    lowerSolution[i] + value;
    }

    for(i = 0; i < sizeRhs - lCols; i++){
	assert(rhs[i] > 0);
    }
	
    for(i = 0; i < numNonBasic; i++){
	memcpy(ray, extRay[i], sizeof(double) * numCols);
	matrixA2G2->times(ray, lCoeffsTimesRay);
	CoinDisjointCopyN(lCoeffsTimesRay, lRows, coeff);
	for(j = 0; j < lCols; j++){
	    value = ray[lColInd[j]];
	    coeff[j + lRows] = value;
	    coeff[j + lRows + lCols] = -1 * value;
	}
	alphaUb = infinity;
	isUnbounded = true;
	for(j = 0; j < sizeRhs; j++){
	    //uselessIneqs[j] = 3;
	    if((uselessIneqs[j] > etol) && (coeff[j] > etol)){
		value = rhs[j]/coeff[j];
		if((alphaUb - value > etol) && (rhs[j] < (infinity/10))){
		    alphaUb = value;
		    isUnbounded = false;
		}
	    }
	}
	assert(alphaUb >= 0);
	if(isUnbounded == false){
	    alphaVec[i] = alphaUb;
	    intersectionFound = true;
	}
	else{
	    alphaVec[i] = -1;
	}
    }

    if(intersectionFound == false){
	localModel_->bS_->shouldPrune_ = true;
    }
	

    delete [] rhs;
    delete [] ray;
    delete [] coeff;
    delete [] lCoeffsTimesLpSol;
    delete [] lCoeffsTimesRay;
    delete [] G2TimeslSol;

    return intersectionFound;

}

//#############################################################################
void
MibSCutGenerator::storeBestSolHypercubeIC(const double* lpSol,
					  std::vector<double> &optLowerObjVec,
					  bool &isTimeLimReached)
{

#ifdef _OPENMPMIBS
  storeBestSolHypercubeICParallel(lpSol, optLowerObjVec, isTimeLimReached);
  return;
#endif

    int maxThreadsLL(localModel_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(localModel_->MibSPar_->entry
		    (MibSParams::whichCutsLL));
    std::string feasCheckSolver(localModel_->MibSPar_->entry
				(MibSParams::feasCheckSolver));
    bool useUBDecompose(localModel_->MibSPar_->entry
			(MibSParams::useUBDecompose));
    int numScenarios(localModel_->getNumScenarios()); 
    double timeLimit(localModel_->AlpsPar()->entry(AlpsParams::timeLimit));
    double remainingTime(0.0);
    
    OsiSolverInterface * oSolver = localModel_->solver();
    MibSBilevel *bS = localModel_->bS_;
    int i(0), j(0), k(0);
    int lpStat;
    bool isUBProvenOptimal(true);
    int numCols(oSolver->getNumCols());
    int uN(localModel_->getUpperDim());
    int lN(localModel_->getLowerDim());
    int truncLN(localModel_->getTruncLowerDim());
    double objVal(0.0);
    double startTimeUB(0.0);
    int * fixedInd = localModel_->getFixedInd();
    
    int useLinkingSolutionPool(localModel_->MibSPar_->entry
		   (MibSParams::useLinkingSolutionPool));

    double *valuesUB = new double[uN + lN];

    std::vector<double> linkSol;
    for(i = 0; i < uN + lN; i++){
	if(fixedInd[i] == 1){
	    linkSol.push_back(lpSol[i]);
	}
    }

    CoinPackedMatrix *truncMatrixG2 = NULL;
    double *multA2XOpt = NULL; 
    OsiSolverInterface *UBSolver = 0;

    //std::vector<double> optLowerObjVec;
    //optLowerObjVec.push_back(optLowerObj);
    
    /*if(bS->UBSolver_){
	bS->UBSolver_ = bS->setUpUBModel(localModel_->getSolver(), optLowerObjVec, false);
    }
    else{
	bS->UBSolver_ = bS->setUpUBModel(localModel_->getSolver(), optLowerObjVec, true);
	}

	UBSolver = bS->UBSolver_;*/

    int numDecomposedProbs(1);

    if((numScenarios > 1) && (useUBDecompose == true) && (localModel_->sizeFixedInd_ == uN)){
      numDecomposedProbs = numScenarios;
      const double *relaxSol = oSolver->getColSolution();
      const double * uObjCoeffs(oSolver->getObjCoefficients());
      double uObjSense = oSolver->getObjSense();
      CoinDisjointCopyN(relaxSol, uN, valuesUB);
      for(i = 0; i < uN; i++){
	objVal += relaxSol[i] * uObjCoeffs[i] * uObjSense;
      } 
    }

    for(i = 0; i < numDecomposedProbs; i++){
      if(numDecomposedProbs == 1){
	UBSolver = bS->setUpUBModel(localModel_->getSolver(), optLowerObjVec, true);
      }
      else{
	if(i == 0){
	  CoinPackedMatrix coefMatrix = *localModel_->origConstCoefMatrix_;
	  int truncNumCols = uN + truncLN;
	  int uRowNum = localModel_->getOrigUpperRowNum();
	  CoinPackedVector col1;
	  CoinPackedVector col2;
	  int *colInd = NULL;
	  double *colElem = NULL;
	  int colNumElem(0);
	  truncMatrixG2 = new CoinPackedMatrix(true, 0, 0);
	  truncMatrixG2->setDimensions(localModel_->getTruncLowerRowNum(), 0);
	  for(j = uN; j < truncNumCols; j++){
	    col1 = coefMatrix.getVector(j);
	    colInd = col1.getIndices();
	    colElem = col1.getElements();
	    colNumElem = col1.getNumElements();
	    for(k = 0; k < colNumElem; k++){
	      col2.insert(colInd[k] - uRowNum, colElem[k]);
	    }
	    truncMatrixG2->appendCol(col2);
	    col1.clear();
	    col2.clear();
	  }
	  int isA2Random(localModel_->MibSPar_->entry(MibSParams::isA2Random));
	  if(isA2Random == PARAM_OFF){
	      multA2XOpt = new double[localModel_->getTruncLowerRowNum()];
	  }
	  else{
	      multA2XOpt = new double[localModel_->getLowerRowNum()];
	  }
	  const double *lpSol = localModel_->getSolver()->getColSolution();
	  double *optUpperSol = new double[uN];
	  CoinDisjointCopyN(lpSol, uN, optUpperSol);
	  localModel_->getStocA2Matrix()->times(optUpperSol, multA2XOpt);
	  delete [] optUpperSol;
	}
	UBSolver = bS->setUpDecomposedUBModel(localModel_->getSolver(),
					      optLowerObjVec, i, truncMatrixG2,
					      multA2XOpt);
      }

      remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getTime();

      if(remainingTime <= localModel_->etol_){
	isTimeLimReached = true;
	bS->shouldPrune_ = true;
	break;
      }
      else{
	remainingTime = CoinMax(remainingTime, 0.00);

	if (feasCheckSolver == "Cbc"){
	    dynamic_cast<OsiCbcSolverInterface *>
		(UBSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
	}else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	    //dynamic_cast<OsiSymSolverInterface *>
	    // (lSolver)->setSymParam("prep_level", -1);

	    sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		(UBSolver)->getSymphonyEnvironment();
	    //Always uncomment for debugging!!
	    sym_set_dbl_param(env, "time_limit", remainingTime);
	    sym_set_int_param(env, "do_primal_heuristic", FALSE);
	    sym_set_int_param(env, "verbosity", -2);
	    //sym_set_int_param(env, "prep_level", -1);
	    sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
	    sym_set_int_param(env, "tighten_root_bounds", FALSE);
	    sym_set_int_param(env, "max_sp_size", 100);
	    sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	    if (whichCutsLL == 0){
		sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	    }else{
		sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	    }
	    if (whichCutsLL == 1){
		sym_set_int_param(env, "generate_cgl_knapsack_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_probing_cuts",
				  DO_NOT_GENERATE);
	        sym_set_int_param(env, "generate_cgl_clique_cuts",
			          DO_NOT_GENERATE);
	        sym_set_int_param(env, "generate_cgl_twomir_cuts",
			          DO_NOT_GENERATE);
	        sym_set_int_param(env, "generate_cgl_flowcover_cuts",
			          DO_NOT_GENERATE);
	    }
#endif
	}else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    UBSolver->setHintParam(OsiDoReducePrint);
	    UBSolver->messageHandler()->setLogLevel(0);
	    CPXENVptr cpxEnv =
		dynamic_cast<OsiCpxSolverInterface*>(UBSolver)->getEnvironmentPtr();
	    assert(cpxEnv);
	    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
	    CPXsetintparam(cpxEnv, CPX_PARAM_CLOCKTYPE, 1);
	    CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, remainingTime);
#endif
	}

	startTimeUB = localModel_->broker_->subTreeTimer().getTime(); 
        UBSolver->branchAndBound();
        localModel_->timerUB_ += localModel_->broker_->subTreeTimer().getTime() - startTimeUB;
        localModel_->counterUB_ ++;

	if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	    if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface *>
					 (UBSolver)->getSymphonyEnvironment())){
		isTimeLimReached = true;
		bS->shouldPrune_ = true;
		break;
	    }
#endif
	}
	else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    lpStat = CPXgetstat(dynamic_cast<OsiCpxSolverInterface*>
				(UBSolver)->getEnvironmentPtr(),
				dynamic_cast<OsiCpxSolverInterface*>
				(UBSolver)->getLpPtr());
	    if((lpStat == CPXMIP_TIME_LIM_FEAS) ||
	       (lpStat == CPXMIP_TIME_LIM_INFEAS)){
		isTimeLimReached = true;
		bS->shouldPrune_ = true;
	        break;
	    }
#endif
	}

	if(UBSolver->isProvenOptimal()){
	  isUBProvenOptimal = true;
	  const double *partialValuesUB = UBSolver->getColSolution();
	  if(numDecomposedProbs == 1){
	    CoinDisjointCopyN(partialValuesUB, lN + uN, valuesUB);
	  }
	  else{
	    int begPos = uN + i * truncLN;
	    CoinDisjointCopyN(partialValuesUB, truncLN, valuesUB + begPos);
	  }
	  objVal += UBSolver->getObjValue() * localModel_->solver()->getObjSense();
	}
	else{
	  isUBProvenOptimal = false;
	  //when all first-level variables are linking,
	  //UBSolver cannot be infeasible 
	  if(numDecomposedProbs > 1){
	    throw CoinError("When all first-level variables are linking, problem UB cannot be infeasible",
			    "storeBestSolHypercubeIC", "MibSCutGenerator");
	  }
	}
      }
      delete UBSolver;
    }

    if(isTimeLimReached == false){
      if(isUBProvenOptimal == true){
	MibSSolution *mibsSol = new MibSSolution(numCols,
						 valuesUB,
						 objVal,
						 localModel_);
	localModel_->storeSolution(BlisSolutionTypeHeuristic, mibsSol);
      }
      else{
	objVal = 10000000;
      }

      if(useLinkingSolutionPool){
	//Add to linking solution pool
	//localModel_->it = localModel_->seenLinkingSolutions.find(linkSol);
	//localModel_->it->second.tag = MibSSetETagUBIsSolved;
	//localModel_->it->second.UBObjVal1 = objVal;
	localModel_->bS_->tagInSeenLinkingPool_ = MibSLinkingPoolTagUBIsSolved;
	localModel_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagUBIsSolved;
	localModel_->seenLinkingSolutions[linkSol].UBObjValue = objVal;   
	if(isUBProvenOptimal == true){
	  localModel_->seenLinkingSolutions[linkSol].UBSolution.clear();
	  //localModel_->it->second.UBSol1.clear();
	  if(numScenarios == 1){
	    for(i = 0; i < uN + lN; i++){
	      //localModel_->it->second.UBSol1.push_back(valuesUB[i]);
	      localModel_->seenLinkingSolutions[linkSol].UBSolution.push_back(valuesUB[i]);
	    }
	  }
	}
      }
    }
    delete truncMatrixG2;
    delete [] multA2XOpt;
    delete [] valuesUB;
}

//#############################################################################
void
MibSCutGenerator::storeBestSolHypercubeICParallel(const double* lpSol,
					  std::vector<double> &optLowerObjVec,
					  bool &isTimeLimReached)
{

  int maxThreadsLL(localModel_->MibSPar_->entry
		   (MibSParams::maxThreadsLL));
  int whichCutsLL(localModel_->MibSPar_->entry
		  (MibSParams::whichCutsLL));
  int maxActiveNodes(localModel_->MibSPar_->entry
		     (MibSParams::maxActiveNodes)); 
  std::string feasCheckSolver(localModel_->MibSPar_->entry
			      (MibSParams::feasCheckSolver));
  bool useUBDecompose(localModel_->MibSPar_->entry
		      (MibSParams::useUBDecompose));
  int numScenarios(localModel_->getNumScenarios());
  double timeLimit(localModel_->AlpsPar()->entry(AlpsParams::timeLimit));
  double remainingTime(0.0);

  OsiSolverInterface * oSolver = localModel_->solver();
  MibSBilevel *bS = localModel_->bS_;
  int i(0), j(0), k(0);
  int lpStat;
  bool isUBProvenOptimal(true);
  int numCols(oSolver->getNumCols());
  int uN(localModel_->getUpperDim());
  int lN(localModel_->getLowerDim());
  int truncLN(localModel_->getTruncLowerDim());
  double objVal(0.0);
  double startTimeUB(0.0);
  int * fixedInd = localModel_->getFixedInd();

  int useLinkingSolutionPool(localModel_->MibSPar_->entry
			     (MibSParams::useLinkingSolutionPool));

  double *valuesUB = new double[uN + lN];

  std::vector<double> linkSol;
  for(i = 0; i < uN + lN; i++){
    if(fixedInd[i] == 1){
      linkSol.push_back(lpSol[i]);
    }
  }

  CoinPackedMatrix *truncMatrixG2 = NULL;
  double *multA2XOpt = NULL;
  //OsiSolverInterface *UBSolver = 0;

  //std::vector<double> optLowerObjVec;
  //optLowerObjVec.push_back(optLowerObj);

  /*if(bS->UBSolver_){
    bS->UBSolver_ = bS->setUpUBModel(localModel_->getSolver(), optLowerObjVec, false);
    }
    else{
    bS->UBSolver_ = bS->setUpUBModel(localModel_->getSolver(), optLowerObjVec, true);
    }
    UBSolver = bS->UBSolver_;*/

  int numDecomposedProbs(1);

  if((numScenarios > 1) && (useUBDecompose == true) && (localModel_->sizeFixedInd_ == uN)){
    numDecomposedProbs = numScenarios;
    const double *relaxSol = oSolver->getColSolution();
    const double * uObjCoeffs(oSolver->getObjCoefficients());
    double uObjSense = oSolver->getObjSense();
    CoinDisjointCopyN(relaxSol, uN, valuesUB);
    for(i = 0; i < uN; i++){
      objVal += relaxSol[i] * uObjCoeffs[i] * uObjSense;
    }
  }

  std::vector<OsiSolverInterface *>UBSolverVec(numDecomposedProbs); 
  for(i = 0; i < numDecomposedProbs; i++){
    remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getWallClock();

    if(remainingTime <= localModel_->etol_){
      isTimeLimReached = true;
      bS->shouldPrune_ = true;
      goto TERM_STOREBESTSOLHypercubeICPARAL;
    }
    if(numDecomposedProbs == 1){
      UBSolverVec[i] = bS->setUpUBModel(localModel_->getSolver(), optLowerObjVec, true);
    }
    else{
      if(i == 0){
	CoinPackedMatrix coefMatrix = *localModel_->origConstCoefMatrix_;
	int truncNumCols = uN + truncLN;
	int uRowNum = localModel_->getOrigUpperRowNum();
	CoinPackedVector col1;
	CoinPackedVector col2;
	int *colInd = NULL;
	double *colElem = NULL;
	int colNumElem(0);
	truncMatrixG2 = new CoinPackedMatrix(true, 0, 0);
	truncMatrixG2->setDimensions(localModel_->getTruncLowerRowNum(), 0);
	for(j = uN; j < truncNumCols; j++){
	  col1 = coefMatrix.getVector(j);
	  colInd = col1.getIndices();
	  colElem = col1.getElements();
	  colNumElem = col1.getNumElements();
	  for(k = 0; k < colNumElem; k++){
	    col2.insert(colInd[k] - uRowNum, colElem[k]);
	  }
	  truncMatrixG2->appendCol(col2);
	  col1.clear();
	  col2.clear();
	}
	int isA2Random(localModel_->MibSPar_->entry(MibSParams::isA2Random));
	if(isA2Random == PARAM_OFF){
	  multA2XOpt = new double[localModel_->getTruncLowerRowNum()];
	}
	else{
	  multA2XOpt = new double[localModel_->getLowerRowNum()];
	}
	const double *lpSol = localModel_->getSolver()->getColSolution();
	double *optUpperSol = new double[uN];
	CoinDisjointCopyN(lpSol, uN, optUpperSol);
	localModel_->getStocA2Matrix()->times(optUpperSol, multA2XOpt);
	delete [] optUpperSol;
      }
      UBSolverVec[i] = bS->setUpDecomposedUBModel(localModel_->getSolver(),
						  optLowerObjVec, i, truncMatrixG2,
						  multA2XOpt);
    }
  }

  for(i = 0; i < numDecomposedProbs; i++){ 
    remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getWallClock();

    if(remainingTime <= localModel_->etol_){
      isTimeLimReached = true;
      bS->shouldPrune_ = true;
      goto TERM_STOREBESTSOLHypercubeICPARAL;
    }
    
      remainingTime = CoinMax(remainingTime, 0.00);

      if (feasCheckSolver == "Cbc"){
	    dynamic_cast<OsiCbcSolverInterface *>
	      (UBSolverVec[i])->getModelPtr()->messageHandler()->setLogLevel(0);
      }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	//dynamic_cast<OsiSymSolverInterface *>
	// (lSolver)->setSymParam("prep_level", -1);

	    sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	      (UBSolverVec[i])->getSymphonyEnvironment();
	    //Always uncomment for debugging!!
	    sym_set_dbl_param(env, "time_limit", remainingTime);
	    sym_set_int_param(env, "do_primal_heuristic", FALSE);
	    sym_set_int_param(env, "verbosity", -2);
	    //sym_set_int_param(env, "prep_level", -1);
	    sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
	    sym_set_int_param(env, "tighten_root_bounds", FALSE);
	    sym_set_int_param(env, "max_sp_size", 100);
	    sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	    if (whichCutsLL == 0){
	      sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	    }else{
	      sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	    }
	    if (whichCutsLL == 1){
	      sym_set_int_param(env, "generate_cgl_knapsack_cuts",
				DO_NOT_GENERATE);
	      sym_set_int_param(env, "generate_cgl_probing_cuts",
				DO_NOT_GENERATE);
	      sym_set_int_param(env, "generate_cgl_clique_cuts",
				DO_NOT_GENERATE);
	      sym_set_int_param(env, "generate_cgl_twomir_cuts",
				DO_NOT_GENERATE);
	      sym_set_int_param(env, "generate_cgl_flowcover_cuts",
				DO_NOT_GENERATE);
	    }
#endif
      }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	UBSolverVec[i]->setHintParam(OsiDoReducePrint);
	UBSolverVec[i]->messageHandler()->setLogLevel(0);
	    CPXENVptr cpxEnv =
	      dynamic_cast<OsiCpxSolverInterface*>(UBSolverVec[i])->getEnvironmentPtr();
	    assert(cpxEnv);
	    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
	    CPXsetintparam(cpxEnv, CPX_PARAM_CLOCKTYPE, 2);
	    CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, remainingTime);
#endif
      }
  }

  startTimeUB = localModel_->broker_->subTreeTimer().getWallClock();
#ifdef _OPENMPMIBS
  omp_set_num_threads(maxActiveNodes);
#pragma omp parallel for
#endif
  for(i = 0; i < numDecomposedProbs; i++){
    UBSolverVec[i]->branchAndBound();
  }
  localModel_->timerUB_ += localModel_->broker_->subTreeTimer().getWallClock() - startTimeUB;
  localModel_->counterUB_ += numDecomposedProbs;

  for(i = 0; i < numDecomposedProbs; i++){
    if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
      if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface *>
				   (UBSolverVec[i])->getSymphonyEnvironment())){
	isTimeLimReached = true;
        bS->shouldPrune_ = true;
        goto TERM_STOREBESTSOLHypercubeICPARAL; 
      }
#endif
    }
    else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
      lpStat = CPXgetstat(dynamic_cast<OsiCpxSolverInterface*>
			  (UBSolverVec[i])->getEnvironmentPtr(),
			  dynamic_cast<OsiCpxSolverInterface*>
			  (UBSolverVec[i])->getLpPtr());
      if((lpStat == CPXMIP_TIME_LIM_FEAS) ||
	 (lpStat == CPXMIP_TIME_LIM_INFEAS)){
	isTimeLimReached = true;
	bS->shouldPrune_ = true;
	goto TERM_STOREBESTSOLHypercubeICPARAL; 
      }
#endif
    }

    if(UBSolverVec[i]->isProvenOptimal()){
      isUBProvenOptimal = true;
      const double *partialValuesUB = UBSolverVec[i]->getColSolution();
      if(numDecomposedProbs == 1){
	CoinDisjointCopyN(partialValuesUB, lN + uN, valuesUB);
      }
      else{
	int begPos = uN + i * truncLN;
	CoinDisjointCopyN(partialValuesUB, truncLN, valuesUB + begPos);
      }
      objVal += UBSolverVec[i]->getObjValue() * localModel_->solver()->getObjSense();
    }
    else{
      isUBProvenOptimal = false;
      //when all first-level variables are linking,
      //UBSolver cannot be infeasible
      if(numDecomposedProbs > 1){
	throw CoinError("When all first-level variables are linking, problem UB cannot be infeasible",
			"storeBestSolHypercubeICParallel", "MibSCutGenerator");
      }
    }
  }
  for(i = 0; i < numDecomposedProbs; i++){
    delete UBSolverVec[i];
  }

  //if(isTimeLimReached == false){
  if(isUBProvenOptimal == true){
    MibSSolution *mibsSol = new MibSSolution(numCols,
					     valuesUB,
					     objVal,
					     localModel_);
    localModel_->storeSolution(BlisSolutionTypeHeuristic, mibsSol);
  }
  else{
    objVal = 10000000;
  }

  if(useLinkingSolutionPool){
    localModel_->bS_->tagInSeenLinkingPool_ = MibSLinkingPoolTagUBIsSolved;
    localModel_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagUBIsSolved;
    localModel_->seenLinkingSolutions[linkSol].UBObjValue = objVal;
    if(isUBProvenOptimal == true){
      localModel_->seenLinkingSolutions[linkSol].UBSolution.clear();
      //localModel_->it->second.UBSol1.clear();
      if(numScenarios == 1){
	for(i = 0; i < uN + lN; i++){
	  //localModel_->it->second.UBSol1.push_back(valuesUB[i]);
	  localModel_->seenLinkingSolutions[linkSol].UBSolution.push_back(valuesUB[i]);
	}
      }
    }
  }
  //}

 TERM_STOREBESTSOLHypercubeICPARAL:
  if(truncMatrixG2){
    delete truncMatrixG2;
  }
  if(multA2XOpt){
    delete [] multA2XOpt;
  }
  if(valuesUB){
    delete [] valuesUB;
  }
}

//#############################################################################
void
MibSCutGenerator::getAlphaHypercubeIC(double** extRay, int numStruct, int numNonBasic,
				      std::vector<double> &alphaVec)
{
    int * fixedInd = localModel_->getFixedInd();
    double etol(localModel_->etol_);

    int i(0), j(0);
    int mult(0);
    double coeff(0.0), tmp(0.0);

    for(i = 0; i < numNonBasic; i++){
	alphaVec[i] = -1;
	for(j = 0; j < numStruct; j++){
	    if(fixedInd[j] == 1){
		coeff = extRay[i][j];
	        if(coeff > etol){
		    mult = 1;
		}
	        else if(coeff < -1*etol){
		    mult = -1;
		}
	        if(fabs(coeff) > etol){
		    tmp = (mult/coeff);
		    if(alphaVec[i] < 0){
			alphaVec[i] = tmp;
		    }
		    else if(alphaVec[i] > tmp){
			alphaVec[i] = tmp;
		    }
		}
	    }
	}
    }
}

//#############################################################################
void
MibSCutGenerator::storeBestSolTenderIC(const double* lpSol, double optLowerObj)
{

    std::string feasCheckSolver =
	localModel_->MibSPar_->entry(MibSParams::feasCheckSolver);
    
    OsiSolverInterface * oSolver = localModel_->solver();
    const CoinPackedMatrix * matrix = oSolver->getMatrixByRow();
    MibSBilevel *bS = localModel_->bS_;

    int numCols(oSolver->getNumCols());
    int uCols(localModel_->getUpperDim());
    int lCols(localModel_->getLowerDim());
    int lRows(localModel_->getLowerRowNum());
    int origNumRows(oSolver->getNumRows());
    double uObjSense(1);
    double lObjSense(localModel_->getLowerObjSense());
    double etol(localModel_->etol_);
    int * uColIndices = localModel_->getUpperColInd();
    int * lColIndices = localModel_->getLowerColInd();
    int * lRowIndices = localModel_->getLowerRowInd();
    int * fixedInd = localModel_->getFixedInd();
    const double * uObjCoeffs = oSolver->getObjCoefficients();
    double * lObjCoeffs = localModel_->getLowerObjCoeffs();
    double * colUb = new double[numCols];
    double * colLb = new double[numCols];

    int * withULVar = new int[lRows]();
    double tmp(0.0);

    int i(0), j(0), indexRow(0), indexCol(0);
    int withULVarSize(0), intCnt(0), pos(0), cnt(0);

    //Get LL constraints in which UL variables participate
    for(i = 0; i < lRows; i++){
	indexRow = lRowIndices[i];
	for(j = 0; j < uCols; j++){
	    indexCol = uColIndices[j];
	    tmp =  matrix->getCoefficient(indexRow, indexCol);
	    if(fabs(tmp) > etol){
		withULVar[i] = 1;
		withULVarSize ++;
		break;
	    }
	}
    }

    int numRows(origNumRows+withULVarSize+1);
    double * rowUb = new double[numRows];
    double * rowLb = new double[numRows];

    CoinFillN(rowLb, numRows, 0.0);
    CoinFillN(rowUb, numRows, 0.0);

    CoinFillN(colLb, numCols, 0.0);
    CoinFillN(colUb, numCols, 0.0);

    /** Set row bounds **/
    for(i = 0; i < origNumRows; i++){
	rowLb[i] = oSolver->getRowLower()[i];
	rowUb[i] = oSolver->getRowUpper()[i];
    }

    rowLb[origNumRows] = -1 * (oSolver->getInfinity());
    rowUb[origNumRows] = optLowerObj * lObjSense;

    /** Set col bounds **/
    for(i = 0; i < numCols; i++){
	colLb[i] = oSolver->getColLower()[i];
	colUb[i] = oSolver->getColUpper()[i];
    }

    OsiSolverInterface * nSolver;

#ifndef COIN_HAS_SYMPHONY
    nSolver = new OsiCbcSolverInterface();
#else
    nSolver = new OsiSymSolverInterface();
#endif

    int * integerVars = new int[numCols];
    double * objCoeffs = new double[numCols];

    CoinFillN(integerVars, numCols, 0);
    CoinFillN(objCoeffs, numCols, 0.0);

    /** Fill in array of integer variables **/
    for(i = 0; i < numCols; i++){
	if(oSolver->isInteger(i)){
	    integerVars[intCnt] = i;
	    intCnt++;
	}
    }

    CoinDisjointCopyN(uObjCoeffs, numCols, objCoeffs);

    CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
    newMat->setDimensions(0, numCols);

    CoinPackedVector row;
    for(i = 0; i < origNumRows; i++){
	for(j = 0; j < numCols; j++){
	    tmp = matrix->getCoefficient(i, j);
	    row.insert(j, tmp);
	}
	newMat->appendRow(row);
	row.clear();
    }

    for(i = 0; i < numCols; i++){
	pos = bS->binarySearch(0, lCols - 1, i, lColIndices);
	if(pos >= 0){
	    tmp = lObjCoeffs[pos] * lObjSense;
	}
	else{
	    tmp = 0.0;
	}
	row.insert(i, tmp);
    }
    newMat->appendRow(row);
    row.clear();

    for(i = 0; i < lRows; i++){
	if(withULVar[i] == 1){
	    cnt++;
	    indexRow = lRowIndices[i];
	    for(j = 0; j < numCols; j++){
		if(fixedInd[j] == 1){
		    tmp =  matrix->getCoefficient(indexRow, j);
		    row.insert(j, tmp);
		    rowLb[origNumRows+cnt] += tmp * lpSol[j];
		    rowUb[origNumRows+cnt] += tmp * lpSol[j];
		}
		else{
		    row.insert(j, 0.0);
		}
	    }
	    newMat->appendRow(row);
	    row.clear();
	}
    }

    nSolver->loadProblem(*newMat, colLb, colUb,
			 objCoeffs, rowLb, rowUb);

    for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
    }

    nSolver->setObjSense(uObjSense); //1 min; -1 max

    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    delete [] integerVars;
    delete [] withULVar;

    //To Do: sahar: write it more efficient
    OsiSolverInterface *nSolver2 = nSolver;

#ifndef COIN_HAS_SYMPHONY
    dynamic_cast<OsiCbcSolverInterface *>
	(nSolver2)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
    dynamic_cast<OsiSymSolverInterface *>
	(nSolver2)->setSymParam("prep_level", -1);

    dynamic_cast<OsiSymSolverInterface *>
	(nSolver2)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *>
	(nSolver)->setSymParam("max_active_nodes", 1);
#endif
    nSolver2->branchAndBound();

    if(nSolver2->isProvenOptimal()){
	const double * newSolution = nSolver2->getColSolution();
	MibSSolution *mibsSol = new MibSSolution(numCols, newSolution,
						 nSolver2->getObjValue(),
						 localModel_);
	
	localModel_->storeSolution(BlisSolutionTypeHeuristic, mibsSol);
    }
}

//#############################################################################
void
MibSCutGenerator::getAlphaTenderIC(double** extRay, int numNonBasic,
				   std::vector<double> &alphaVec)
{
    OsiSolverInterface * hSolver = localModel_->solver();
    const CoinPackedMatrix * matrix = hSolver->getMatrixByRow();

    int numCols(hSolver->getNumCols());
    int uCols(localModel_->getUpperDim());
    int lRows(localModel_->getLowerRowNum());
    double etol(localModel_->etol_);
    int * uColIndices = localModel_->getUpperColInd();
    int * lRowIndices = localModel_->getLowerRowInd();
    int * fixedInd = localModel_->getFixedInd();

    int withULVarSize(0);
    int * withULVar = new int[lRows]();

    int i(0), j(0), k(0);
    int indexRow(0), indexCol(0), mult(0);
    double tmp(0.0), coeff(0.0);

    //Get LL constraints in which UL variables participate
    for(i = 0; i < lRows; i++){
	indexRow = lRowIndices[i];
	for(j = 0; j < uCols; j++){
	    indexCol = uColIndices[j];
	    tmp =  matrix->getCoefficient(indexRow, indexCol);
	    if(fabs(tmp) > etol){
		withULVar[i] = 1;
		withULVarSize ++;
		break;
	    }
	}
    }

    for(i = 0; i < numNonBasic; i++){
	alphaVec[i] = -1;
	for(j = 0; j < lRows; j++){
	    if(withULVar[j] == 1){
		indexRow = lRowIndices[j];
		coeff = 0;
		for(k = 0; k < numCols; k++){
		    if(fixedInd[k] == 1){
			tmp =  matrix->getCoefficient(indexRow, k);
			coeff += tmp * extRay[i][k];
		    }
		}
		if(coeff > etol){
		    mult = 1;
		}
		else if(coeff < -1*etol){
		    mult = -1;
		}
		if(fabs(coeff) > etol){
		    tmp = (mult/coeff);
		    if(alphaVec[i] < 0){
			alphaVec[i] = tmp;
		    }
		    else if(alphaVec[i] > tmp){
			alphaVec[i] = tmp;
		    }
		}
	    }
	}
    }

    delete [] withULVar;
}

//############################################################################# 
int
MibSCutGenerator::boundCuts(BcpsConstraintPool &conPool, double *passedObjCoeffs, double &passedRhs,
			    bool &isInfeasible)
{
    /* Derive a bound on the lower level objective function value */
    bool boundCutOptimal
      = localModel_->MibSPar_->entry(MibSParams::boundCutOptimal);

    bool boundCutRelaxUpper
      = localModel_->MibSPar_->entry(MibSParams::boundCutRelaxUpper);

      //parametric or non-parametric
    int boundCutOptimalType =
      localModel_->MibSPar_->entry(MibSParams::boundCutOptimalType);

    std::string feasCheckSolver
	= localModel_->MibSPar_->entry(MibSParams::feasCheckSolver);

    int problemType
	= localModel_->MibSPar_->entry(MibSParams::bilevelProblemType);

    bool allowRemoveCut(localModel_->MibSPar_->entry
			(MibSParams::allowRemoveCut));

    double boundCutTimeLim
	= localModel_->MibSPar_->entry(MibSParams::boundCutTimeLim);

    int boundCutNodeLim
	= localModel_->MibSPar_->entry(MibSParams::boundCutNodeLim);

    MibSTreeNode * node =
      dynamic_cast<MibSTreeNode *>(localModel_->activeNode_);

    if(passedObjCoeffs != NULL){
      boundCutOptimal = true;
      boundCutOptimalType = MibSBoundCutOptimalTypeNonParametric;
    }

    if ((localModel_->boundingPass_ > 1) && (passedObjCoeffs == NULL)){
	return 0;
    }

    /** Set up lp solver **/
    OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);
    double lObjSense = localModel_->getLowerObjSense();
    double * lObjCoeffs = localModel_->getLowerObjCoeffs();
    int lCols(localModel_->getLowerDim());
    int * lColIndices = localModel_->getLowerColInd();
    int uCols(localModel_->getUpperDim());
    int * uColIndices = localModel_->getUpperColInd();
    double * origColLb = localModel_->getOrigColLb();
    double * origColUb = localModel_->getOrigColUb();
    int i(0), index(0);
    double etol(localModel_->etol_);

    OsiSolverInterface * oSolver = localModel_->getSolver();
    double infinity(oSolver->getInfinity());
    bool isTimeLimReached(false);

    CoinPackedMatrix matrix = *oSolver->getMatrixByCol();
    CoinPackedMatrix rowMatrix = *oSolver->getMatrixByRow();
    int numCols = localModel_->getLowerDim() + localModel_->getUpperDim();
    int auxCols = oSolver->getNumCols() - numCols;
    int *indDel = new int[auxCols];
    CoinIotaN(indDel, auxCols,numCols);
    matrix.deleteCols(auxCols, indDel);

    int tCols(numCols);
    double * nObjCoeffs = new double[tCols];

    CoinZeroN(nObjCoeffs, tCols);

    if(passedObjCoeffs != NULL){
	memcpy(nObjCoeffs, passedObjCoeffs, sizeof(double) * tCols);
    }
    else{
	for(i = 0; i < lCols; i++){
	    index = lColIndices[i];
	    nObjCoeffs[index] = -lObjSense * lObjCoeffs[i];
	}
    }

    int numCuts(0);
    double objval = -1 * infinity;
    double lower_objval = -1 * infinity;
    char * colType = new char[tCols];

    if(boundCutOptimal){
      if(((boundCutOptimalType == MibSBoundCutOptimalTypeParametric) &&
	  (node->getIndex() == 0)) ||
	 (boundCutOptimalType == MibSBoundCutOptimalTypeNonParametric) ||
	 (passedObjCoeffs != NULL)){
	/** Create new MibS model to solve bilevel **/
	MibSModel *boundModel = new MibSModel();
	boundModel->setSolver(&lpSolver);
	boundModel->AlpsPar()->setEntry(AlpsParams::msgLevel, -1);
	//boundModel->AlpsPar()->setEntry(AlpsParams::msgLevel, 1000);
	boundModel->AlpsPar()->setEntry(AlpsParams::nodeLimit, boundCutNodeLim);
	boundModel->AlpsPar()->setEntry(AlpsParams::timeLimit, boundCutTimeLim);
	boundModel->BlisPar()->setEntry(BlisParams::heurStrategy, 0);
	boundModel->MibSPar()->setEntry(MibSParams::feasCheckSolver, feasCheckSolver.c_str());
	boundModel->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
	boundModel->MibSPar()->setEntry(MibSParams::printProblemInfo, false);
	boundModel->MibSPar()->setEntry(MibSParams::useBoundCut, false);
	if(passedObjCoeffs != NULL){
	    boundModel->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyLinking);
	    boundModel->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeWatermelon, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeTenderIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeHybridIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::bilevelCutTypes, -1);
	    boundModel->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useNewPureIntCut, false);
	}
	else{
	  //sahar: these ones should be modified
	    boundModel->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyLinking);
	    boundModel->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_ON);
	    boundModel->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeWatermelon, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeTenderIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useTypeHybridIC, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
	    boundModel->MibSPar()->setEntry(MibSParams::useNewPureIntCut, false);
	    if(boundCutOptimalType == MibSBoundCutOptimalTypeParametric){
	      boundModel->AlpsPar()->setEntry(AlpsParams::deletePrunedNodes, true);
	      boundModel->AlpsPar()->setEntry(AlpsParams::deleteDeadNode, false);
	      boundModel->BlisPar()->setEntry(BlisParams::difference, -2);  
	    }
	}
	    
	double *colUpper = new double[tCols];
	double *colLower = new double[tCols];
	memcpy(colLower, oSolver->getColLower(), sizeof(double) * tCols);
	memcpy(colUpper, oSolver->getColUpper(), sizeof(double) * tCols);
	memcpy(colType, localModel_->colType_, tCols);

	double *lColLbInLProb = new double[lCols];
	double *lColUbInLProb = new double[lCols];
	for(i = 0; i < lCols; i++){
	    index = lColIndices[i];
	    lColLbInLProb[i] = origColLb[index];
	    lColUbInLProb[i] = origColUb[index];
	}

	//interdiction problems
	if(problemType == 1){
	    boundModel->isInterdict_ = true;
	}

	boundModel->loadAuxiliaryData(localModel_->getLowerDim(),
				      localModel_->getLowerRowNum(),
				      localModel_->getLowerColInd(),
				      localModel_->getLowerRowInd(),
				      localModel_->getLowerObjSense(),
				      localModel_->getLowerObjCoeffs(),
				      localModel_->getUpperDim(),
				      localModel_->getUpperRowNum(),
				      localModel_->getUpperColInd(),
				      localModel_->getUpperRowInd(),
				      localModel_->structRowNum_,
				      localModel_->structRowInd_,
				      0, NULL, lColLbInLProb, lColUbInLProb);

	boundModel->loadProblemData(matrix, rowMatrix, colLower, colUpper, nObjCoeffs,
				    oSolver->getRowLower(), oSolver->getRowUpper(),
				    colType, 1.0, oSolver->getInfinity(),
				    oSolver->getRowSense());

	int argc = 1;
	char** argv = new char* [1];
	argv[0] = (char *) "mibs";

#ifdef  COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, *boundModel);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, *boundModel);
#endif

	broker.search(boundModel);

	if (boundModel->getNumSolutions() > 0){
	    double *solution = boundModel->incumbent();
	}

	broker.printBestSolution();
	if(broker.getSolStatus() == AlpsExitStatusInfeasible){
	    isInfeasible = true;
	}
	else{
	    if(broker.getBestNode()){
		lower_objval = broker.getBestNode()->getQuality();
	    }
	    else{
		if(boundModel->getNumSolutions() > 0){
		    lower_objval = boundModel->getKnowledgeBroker()->getBestQuality();
		}
		else{
		    lower_objval = -1 * infinity;
		}
	    }
	}

	if(boundCutOptimalType == MibSBoundCutOptimalTypeParametric){
	  AlpsSubTree *boundProbTree = broker.getWorkingSubTree();
	  localModel_->setBoundProbRoot(boundProbTree->getRoot());
	  localModel_->setBoundProbLinkingPool(boundModel->seenLinkingSolutions);
	  int numStoredCuts(0);
	  findLeafNodes(localModel_->getBoundProbRoot(), &numStoredCuts,
			&localModel_->boundProbLeafNum_,
			localModel_->boundProbCutPoolStarts_,
			localModel_->boundProbCutPoolIndices_,
			localModel_->boundProbCutPoolValues_,
			localModel_->boundProbCutPoolBounds_,
			localModel_->boundProbCutPoolSource_,
			localModel_->boundProbLeafNodeCutInf_,
			localModel_->boundProbLeafNodeCutStarts_,
			localModel_->boundProbLeafLb_,
			localModel_->boundProbLeafUb_);

	  if(isTimeLimReached == true){
	    goto TERM_BOUNDCUT;
	  }
	}
	
	delete [] colLower;
        delete [] colUpper;
	delete [] lColLbInLProb;
	delete [] lColUbInLProb;
	delete boundModel;
      }
    }else if (localModel_->getNumSolutions() > 0){
	memcpy(colType, localModel_->colType_, tCols);
	for (i = 0; i < uCols; i++){
	    colType[uColIndices[i]] = 'C';
	}

	//Create new upperbound for lower level variables (to fix them)
	BlisSolution* blisSol = dynamic_cast<BlisSolution*>(
	   localModel_->getKnowledgeBroker()->getBestKnowledge(
			       AlpsKnowledgeTypeSolution).first);
	   const double * sol;
	   sol = blisSol->getValues();
	   int lN(localModel_->getLowerDim());
	   int * lowerColInd = localModel_->getLowerColInd();
	   int LowZero(0);
	   std::vector<int> zeroList;

	   for (i=0; i<lN; i++){
	       if ((sol[lowerColInd[i]]>-etol)&&(sol[lowerColInd[i]]<etol)){
		   zeroList.push_back(lowerColInd[i]);
		   LowZero ++;
	       }
	   }
	   double *NewColUpper = new double[numCols];
	   memcpy(NewColUpper, oSolver->getColUpper(), sizeof(double) * numCols);
	   for(i=0; i<LowZero/2; i++){
	       index = zeroList[i];
	       NewColUpper[index] = 0;
	   }
	   /** Create new MibS model to solve bilevel(with new upperbound) **/
	   MibSModel NewboundModel;
	   NewboundModel.setSolver(&lpSolver);
	   NewboundModel.AlpsPar()->setEntry(AlpsParams::msgLevel, -1);

	   NewboundModel.loadAuxiliaryData(localModel_->getLowerDim(),
					   localModel_->getLowerRowNum(),
					   localModel_->getLowerColInd(),
					   localModel_->getLowerRowInd(),
					   localModel_->getLowerObjSense(),
					   localModel_->getLowerObjCoeffs(),
					   localModel_->getUpperDim(),
					   localModel_->getUpperRowNum(),
					   localModel_->getUpperColInd(),
					   localModel_->getUpperRowInd(),
					   localModel_->structRowNum_,
					   localModel_->structRowInd_,
					   0, NULL, NULL, NULL);
	   NewboundModel.loadProblemData(matrix, rowMatrix,
					 oSolver->getColLower(), NewColUpper,
					 nObjCoeffs,
					 oSolver->getRowLower(),
					 oSolver->getRowUpper(),
					 localModel_->colType_, 1.0,
					 oSolver->getInfinity(), oSolver->getRowSense());

	   int argc1 = 1;
	   char** argv1 = new char* [1];
	   argv1[0] = (char *) "mibs";

#ifdef  COIN_HAS_MPI
	   AlpsKnowledgeBrokerMPI Newbroker(argc1, argv1, NewboundModel);
#else
	   AlpsKnowledgeBrokerSerial Newbroker(argc1, argv1, NewboundModel);
#endif

	   NewboundModel.MibSPar()->setEntry(MibSParams::bilevelCutTypes, 1);
	   NewboundModel.MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_ON);

	   NewboundModel.MibSPar()->setEntry(MibSParams::useLowerObjHeuristic, false);
	   NewboundModel.MibSPar()->setEntry(MibSParams::useObjCutHeuristic, false);
	   NewboundModel.MibSPar()->setEntry(MibSParams::useWSHeuristic, false);
	   NewboundModel.MibSPar()->setEntry(MibSParams::useGreedyHeuristic, false);

	   Newbroker.search(&NewboundModel);

	   Newbroker.printBestSolution();
	   if (NewboundModel.getNumSolutions() > 0){
	       objval = NewboundModel.getKnowledgeBroker()->getBestQuality();
#if 0
	       double objval1(boundModel.getKnowledgeBroker()->getBestQuality());
	       std::cout<<"obj="<<objval1<<std::endl;
	       std::cout<<"objWithRest="<<objval<<std::endl;
#endif
	   }
	   zeroList.clear();
    }

    if((boundCutOptimalType == MibSBoundCutOptimalTypeParametric) &&
       (node->getIndex() > 0)){
      double cutLb(0.0);
      cutLb = getRhsParamBoundCut(&isTimeLimReached);
      if(isTimeLimReached == true){
	goto TERM_BOUNDCUT;
      }
      if(cutLb + etol >= infinity){
	isInfeasible = true;
      }
      lower_objval = cutLb;
    }

    if(passedObjCoeffs != NULL){
	passedRhs = -1 * lower_objval;
    }
    else if(isInfeasible == true){
	localModel_->bS_->shouldPrune_ = true;
    }
    else{
	if (lower_objval > -1 * infinity){
	    double cutub(infinity);
	    std::vector<int> indexList;
	    std::vector<double> valsList;
	    for(i = 0; i < lCols; i++){
		if (lObjCoeffs[i] == 0.0){
		    continue;
		}
		index = lColIndices[i];
		indexList.push_back(index);
		valsList.push_back(-lObjSense *lObjCoeffs[i]);
	    }
	    numCuts += addCut(conPool, lower_objval, cutub, indexList, valsList,
			      allowRemoveCut);
	    indexList.clear();
	    valsList.clear();
	}
    }

 TERM_BOUNDCUT:
    delete[] indDel;
    delete [] colType;
    delete [] nObjCoeffs;

    return numCuts;
}

//############################################################################# 
void
MibSCutGenerator::findLeafNodes(AlpsTreeNode *node, int *numStoredCuts,
				int *numLeafNodes,
				std::vector<int> &cutStarts,
				std::vector<int> &cutIndices,
				std::vector<double> &cutValues,
				std::vector<double> &cutBounds,
				std::vector<int> &sourceNode,
				std::vector<int> &leafNodeCutInf,
				std::vector<int> &leafNodeCutStatrs,
				std::vector<double> &leafNodeLBs,
				std::vector<double> &leafNodeUBs)
{

  assert(node);

  int i(0), j(0), k(0);
  int tmp(0);
  double etol(localModel_->etol_);
  int numChildren = node->getNumChildren();


  if(node->getIndex() == 0){
    cutStarts.push_back(0);
    leafNodeCutStatrs.push_back(0);
  }

  if(numChildren == 0){
    *numLeafNodes = *numLeafNodes + 1;
    int numModify(0), index(0), nodeIndex(0), tmpInt(0);
    int numAllGenCut(0), pos(0), tmpStart(0), constSize(0);
    int tmpNumCuts(0);   
    double value(0.0);
    int numCoreVars(localModel_->getUpperDim() +
		    localModel_->getLowerDim());
    BlisNodeDesc* pathDesc = NULL;
    BlisConstraint *aCon = NULL;
    double *varHardLB = NULL;
    double *varHardUB = NULL;
    MibSTreeNode *tmpMibsNode = NULL;
    std::vector<AlpsTreeNode *> leafToRootPath;
    double *startColLB = new double[numCoreVars];
    double *startColUB = new double[numCoreVars];
    CoinFillN(startColLB, numCoreVars, -ALPS_DBL_MAX);
    CoinFillN(startColUB, numCoreVars, ALPS_DBL_MAX);
    AlpsTreeNode *parent = node->getParent();
    leafToRootPath.push_back(node);
    //note: in bounding problem the only explicit node is the root node.
    while(parent){
      leafToRootPath.push_back(parent);
      if(parent->getIndex() > 0){
	parent = parent->getParent();
      }
      else{
	break;
      }
    }
    for(i = static_cast<int> (leafToRootPath.size() - 1); i > -1; --i){
      pathDesc = dynamic_cast<BlisNodeDesc*>((leafToRootPath.at(i))->getDesc());

      varHardLB = pathDesc->getVars()->lbHard.entries;
      varHardUB = pathDesc->getVars()->ubHard.entries;
      numModify = pathDesc->getVars()->lbHard.numModify;
      for(k = 0; k < numModify; ++k){
	index = pathDesc->getVars()->lbHard.posModify[k];
	value = pathDesc->getVars()->lbHard.entries[k];
	startColLB[index] = CoinMax(startColLB[index], value);
      }
      
      numModify = pathDesc->getVars()->ubHard.numModify;
      for(k = 0; k < numModify; ++k){
	index = pathDesc->getVars()->ubHard.posModify[k];
	value = pathDesc->getVars()->ubHard.entries[k];
	startColUB[index] = CoinMin(startColUB[index], value);
      }
      
      numModify = pathDesc->getVars()->lbSoft.numModify;
      for(k = 0; k < numModify; ++k){
	index = pathDesc->getVars()->lbSoft.posModify[k];
	value = pathDesc->getVars()->lbSoft.entries[k];
	startColLB[index] = CoinMax(startColLB[index], value);
      }
      
      numModify = pathDesc->getVars()->ubSoft.numModify;
      for(k = 0; k < numModify; ++k){
	index = pathDesc->getVars()->ubSoft.posModify[k];
	value = pathDesc->getVars()->ubSoft.entries[k];
	startColUB[index] = CoinMin(startColUB[index], value);
      }

      tmpInt = pathDesc->getCons()->numAdd;
      tmpMibsNode = dynamic_cast<MibSTreeNode*>(leafToRootPath.at(i));
      nodeIndex = tmpMibsNode->getIndex();
      if(tmpMibsNode->getIsCutStored() == true){ 
	if(tmpInt > 0){
	  numAllGenCut = sourceNode.size();
	  for(k = 0; k < numAllGenCut; k++){
	    if(sourceNode[k] == nodeIndex){
	      pos = k;
	      break;
	    }
	  }
	  for(k = pos; k < pos + tmpInt; k++){
	    leafNodeCutInf.push_back(k);
	    tmpNumCuts ++;
	  }
	}
      }
      else{
	tmpMibsNode->setIsCutStored(true);
	for(k = 0; k < tmpInt; ++k){
	  leafNodeCutInf.push_back(*numStoredCuts);
	  tmpNumCuts ++;    
	  *numStoredCuts = *numStoredCuts + 1;
	  aCon = dynamic_cast<BlisConstraint *>
	    (pathDesc->getCons()->objects[k]);
	  cutBounds.push_back(aCon->getLbSoft());
	  cutBounds.push_back(aCon->getUbSoft());
	  constSize = aCon->getSize();
	  cutIndices.insert(cutIndices.end(), aCon->getIndices(),
			    aCon->getIndices() + constSize);
	  cutValues.insert(cutValues.end(), aCon->getValues(),
			   aCon->getValues() + constSize);
	  tmpStart = cutStarts.at(cutStarts.size() - 1);
	  cutStarts.push_back(tmpStart + constSize);
	  sourceNode.push_back(nodeIndex);
	}
      }
    }
    tmpStart = leafNodeCutStatrs.at(leafNodeCutStatrs.size() - 1);
    leafNodeCutStatrs.push_back(tmpStart + tmpNumCuts);
    leafNodeLBs.insert(leafNodeLBs.end(), startColLB,
		       startColLB + numCoreVars);
    leafNodeUBs.insert(leafNodeUBs.end(), startColUB,
		       startColUB + numCoreVars);

    leafToRootPath.clear();
    delete [] startColLB;
    delete [] startColUB;
  }
  else{
    for(j = 0; j < numChildren; j++){
      findLeafNodes(node->getChild(j), numStoredCuts, numLeafNodes,
		    cutStarts, cutIndices, cutValues, cutBounds, sourceNode,
		    leafNodeCutInf, leafNodeCutStatrs, leafNodeLBs, leafNodeUBs);
    }
  }
}

//#############################################################################
double
MibSCutGenerator::getRhsParamBoundCut(bool *isTimeLimReached)
{
  int i(0);
  double etol(localModel_->etol_);
  double cutLb(localModel_->getSolver()->getInfinity());
  int leafNum(localModel_->boundProbLeafNum_);
  double leafNodeLB(0.0);

  for(i = 0; i < leafNum; i++){
    leafNodeLB = solveLeafNode(i, isTimeLimReached);
    if(leafNodeLB < cutLb - etol){
      cutLb = leafNodeLB;
    } 
  }

  return cutLb;
} 
//#############################################################################
double
MibSCutGenerator::solveLeafNode(int leafNodeIndex, bool *isTimeLimReached)
{

  std::string feasCheckSolver =
    localModel_->MibSPar_->entry(MibSParams::feasCheckSolver);
  
  OsiSolverInterface * hSolver = localModel_->getSolver();

  int i;
  int index(0);
  double bound(0.0), value(0.0), lowerObjValue(0.0);
  double etol(localModel_->etol_);
  double infinity(hSolver->getInfinity());
  int uCols(localModel_->getUpperDim());
  int lCols(localModel_->getLowerDim());
  int numCols(uCols + lCols);
  int numOrigRows(localModel_->getNumOrigCons());
  double *origRowLb(localModel_->getOrigRowLb());
  double *origRowUb(localModel_->getOrigRowUb());

  int numRows(0);
  double lObjSense(localModel_->getLowerObjSense());

  bool isLinkFixed(true), isBoundChanged(false);
  std::vector<double> linkSol;

  double *lObjCoeffs = localModel_->getLowerObjCoeffs();
  int *uColInd = localModel_->getUpperColInd();
  int *lColInd = localModel_->getLowerColInd();
  int *fixedInd = localModel_->getFixedInd();

  double *leafColLb = new double[numCols];
  double *leafColUb = new double[numCols];
  CoinZeroN(leafColLb, numCols);
  CoinZeroN(leafColUb, numCols);
  std::copy(localModel_->boundProbLeafLb_.begin() + leafNodeIndex * numCols,
	    localModel_->boundProbLeafLb_.begin() + (leafNodeIndex + 1) * numCols,
	    leafColLb);

  std::copy(localModel_->boundProbLeafUb_.begin() + leafNodeIndex * numCols,
	    localModel_->boundProbLeafUb_.begin() + (leafNodeIndex + 1) * numCols,
	    leafColUb);
  
  const double *colLb = hSolver->getColLower();
  const double *colUb = hSolver->getColUpper();

  double *newColLb = new double[numCols];
  double *newColUb = new double[numCols];
  //CoinZeroN(newColLb, numCols);
  //CoinZeroN(newColUb, numCols);
  //CoinCopyN(leafColLb, numCols, newColLb);
  //CoinCopyN(leafColUb, numCols, newColUb);
  memcpy(newColLb, leafColLb, sizeof(double) * numCols);
  memcpy(newColUb, leafColUb, sizeof(double) * numCols);

  double *upperSol = new double[numCols];
  //CoinZeroN(upperSol, uCols);

  /*std::map<std::vector<double>, LINKING_SOLUTION> linkingPool
    = localModel_->getBoundProbLinkingPool();*/

  for(i = 0; i < numCols; i++){
    if(colLb[i] - leafColLb[i] > etol){
      newColLb[i] = colLb[i];
      isBoundChanged = true;
    }
    if(leafColUb[i] - colUb[i] > etol){
      newColUb[i] = colUb[i];
      isBoundChanged = true;
    }
    if(newColLb[i] - newColUb[i] > etol){
      bound = infinity;
      goto TERM_SOLVELEAFNODE;
    }
  }

  memcpy(upperSol, newColLb, sizeof(double) * numCols);
  for(i = 0; i < uCols; i++){
      index = uColInd[i];
      if(fixedInd[index] == 1){
	  if(fabs(newColLb[index] - newColUb[index]) <= etol){
	      value = floor(newColLb[index] + 0.5);
	      linkSol.push_back(floor(value));
	      upperSol[i] = value;
	  }
	  else{
	      isLinkFixed = false;
	      break;
	  }
      }
  }

  if(isLinkFixed == true){
    //should solve UB with new bounds
    //first, if the bounds are not changed, we check the linking
    //pool of bounding problem to see if the corresponding UB problem
    //is solved or not.
    //If not, we check to see if there is linkSol in the linking pool
    //of either bounding problem or main problem. If so, we can get
    //the optimal objective of second-level problem from the pool.
    /*if(linkingPool.find(linkSol) != linkingPool.end()){
      if(linkingPool[linkSol].tag == MibSLinkingPoolTagUBIsSolved){
	if(linkingPool[linkSol].UBSolution.size() <= 1){
	  bound = infinity;
	  goto TERM_SOLVELEAFNODE; 
	}
	else{
	  if(isBoundChanged == false){
	    bound = linkingPool[linkSol].UBObjValue;
	    goto TERM_SOLVELEAFNODE; 
	  }
	  //lowerObjValue = linkingPool[linkSol].lowerObjValue;
	}
      }  
      else if(linkingPool[linkSol].tag == MibSLinkingPoolTagLowerIsFeasible){
	lowerObjValue = linkingPool[linkSol].lowerObjValue;
      }
      else{
	bound = infinity;
	goto TERM_SOLVELEAFNODE;
      }
    }
    else if(localModel_->seenLinkingSolutions.find(linkSol) !=
	    localModel_->seenLinkingSolutions.end()){
      if(localModel_->seenLinkingSolutions[linkSol].tag ==
	 MibSLinkingPoolTagUBIsSolved){
	if(linkingPool[linkSol].UBSolution.size() <= 1){
	  bound = infinity;
	  goto TERM_SOLVELEAFNODE;
	}
	else{
	  lowerObjValue = localModel_->seenLinkingSolutions[linkSol].lowerObjValue;
	}
      }
      else if(localModel_->seenLinkingSolutions[linkSol].tag ==
	      MibSLinkingPoolTagLowerIsFeasible){
	lowerObjValue = localModel_->seenLinkingSolutions[linkSol].lowerObjValue;
      }
      else{
	bound = infinity;
	goto TERM_SOLVELEAFNODE;
      }
    }
    else{*/
      //we should solve second-level problem
      //sahar: we can add the result of solving second-level problem to
      //the linking pool of main problem to use it later, but the current
      //format of addSolutionToSeenLinkingSolutionPool does not allow it.
      //the current implementation does not get the linkSol as an argument.
      //change the format later.
      MibSBilevel *bS = localModel_->getMibSBilevel();
      if(bS->lSolver_){
	  bS->lSolver_ = bS->setUpModel(hSolver, false, 0, upperSol);
      }
      else{
	  bS->lSolver_ = bS->setUpModel(hSolver, true, 0, upperSol);
      }
      OsiSolverInterface *lSolver = bS->lSolver_;
      solveMips(lSolver);

      if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					     (dynamic_cast<OsiSymSolverInterface *>
					      (lSolver)->getSymphonyEnvironment()))){
	*isTimeLimReached = true;
	goto TERM_SOLVELEAFNODE;
      }

      if(lSolver->isProvenOptimal()){
	lowerObjValue = lSolver->getObjValue() * lObjSense;
	
      }
      else{
	bound = infinity;
	goto TERM_SOLVELEAFNODE;
      }
      //}
    //setting up UB
    CoinPackedMatrix *newMatrix = new CoinPackedMatrix(false, 0, 0);
    newMatrix->setDimensions(0, numCols);
    numRows = numOrigRows;
    CoinPackedMatrix origMatrix = *localModel_->origConstCoefMatrix_;
    origMatrix.reverseOrdering();
    for(i = 0; i < numOrigRows; i++){
	newMatrix->appendRow(origMatrix.getVector(i));
    }
    
    int numRowsUB(numRows + 1);
    double *rowLbUB = new double[numRowsUB];
    double *rowUbUB = new double[numRowsUB];
    CoinDisjointCopyN(origRowLb, numRows, rowLbUB);
    CoinDisjointCopyN(origRowUb, numRows, rowUbUB);
    
    OsiSolverInterface *UBSolver;
    if (feasCheckSolver == "Cbc"){
      UBSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
      UBSolver = new OsiSymSolverInterface();
#else
      throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
		      "solveLeafNode", "MibSCutGenerator");
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
      UBSolver = new OsiCpxSolverInterface();
#else
      throw CoinError("CPLEX chosen as solver, but it has not been enabled",
		      "solveLeafNode", "MibSCutGenerator");
#endif
    }else{
      throw CoinError("Unknown solver chosen",
		      "solveLeafNode", "MibSCutGenerator");
    }
    
    CoinPackedVector row;
    double *objCoeffs = new double[numCols];
    CoinZeroN(objCoeffs, numCols);
    for(i = 0; i < lCols; i++){
      value = lObjCoeffs[i] * lObjSense;
      row.insert(lColInd[i], value);
      objCoeffs[lColInd[i]] = -1 * value;
    }

    newMatrix->appendRow(row);

    row.clear();
    
    rowLbUB[numRows] = -1 * infinity;
    rowUbUB[numRows] = lowerObjValue;

    UBSolver->loadProblem(*newMatrix, newColLb, newColUb,
			  objCoeffs, rowLbUB, rowUbUB);

    for(i = 0; i < numCols; i++){
      if(hSolver->isInteger(i)){
	UBSolver->setInteger(i);
      }
    }

    UBSolver->setObjSense(1);

    UBSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    solveMips(UBSolver);

    if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					   (dynamic_cast<OsiSymSolverInterface *>
					    (UBSolver)->getSymphonyEnvironment()))){
      *isTimeLimReached = true;
    }

    if(UBSolver->isProvenOptimal()){
      bound = UBSolver->getObjValue() * UBSolver->getObjSense();

    }
    else{
      bound = infinity;
    }

    delete [] rowLbUB;
    delete [] rowUbUB;
    delete [] objCoeffs;
    delete newMatrix;
    delete UBSolver;
  }
  else{
    std::vector<double> rowLbVec;
    std::vector<double> rowUbVec;
    CoinPackedMatrix *newMatrix = getLeafConst(leafNodeIndex, numRows, rowLbVec, rowUbVec);
    double *rowLb = new double[numRows];
    double *rowUb = new double[numRows];
    CoinDisjointCopyN(origRowLb, numOrigRows, rowLb);
    CoinDisjointCopyN(origRowUb, numOrigRows, rowUb);
    std::copy(rowLbVec.begin(), rowLbVec.end(), rowLb + numOrigRows);
    std::copy(rowUbVec.begin(), rowUbVec.end(), rowUb + numOrigRows);
    double *objCoeffs = new double[numCols];
    CoinZeroN(objCoeffs, numCols);
    for(i = 0; i < lCols; i++){
      value = lObjCoeffs[i] * lObjSense;
      objCoeffs[lColInd[i]] = -1 * value;
    }
    
    //based on value of relaxTypeParamBoundCut, we solve
    //problem S or P of leaf node with new bounds
    int relaxType(localModel_->MibSPar_->entry
		  (MibSParams::relaxTypeParamBoundCut));
    
    if(relaxType == MibSRelaxTypeParamBoundCutMIP){

      OsiSolverInterface *relaxSolver;
      if (feasCheckSolver == "Cbc"){
	relaxSolver = new OsiCbcSolverInterface();
      }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	relaxSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"solveLeafNode", "MibSCutGenerator");
#endif
      }else if (feasCheckSolver == "CPLEX"){
	#ifdef COIN_HAS_CPLEX
	relaxSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"solveLeafNode", "MibSCutGenerator");
#endif
      }else{
	throw CoinError("Unknown solver chosen",
			"solveLeafNode", "MibSCutGenerator");
      }

      relaxSolver->loadProblem(*newMatrix, newColLb, newColUb,
			       objCoeffs, rowLb, rowUb);

      for(i = 0; i < numCols; i++){
	if(hSolver->isInteger(i)){
	  relaxSolver->setInteger(i);
	}
      }

      relaxSolver->setObjSense(1);

      relaxSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

      solveMips(relaxSolver);

      if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					     (dynamic_cast<OsiSymSolverInterface *>
					      (relaxSolver)->getSymphonyEnvironment()))){
	*isTimeLimReached = true;
      }
      
      if(relaxSolver->isProvenOptimal()){
	bound = relaxSolver->getObjValue() * relaxSolver->getObjSense();
      }
      else{
	bound = infinity;
      }

      delete relaxSolver;
    }
    else if(relaxType == MibSRelaxTypeParamBoundCutLP){
      double timeLimit(localModel_->AlpsPar()->entry(AlpsParams::timeLimit));

      double remainingTime(0.0);

      remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getTime();

      //sahar:Fix this by setting the time limit for relax solver
      if(remainingTime < etol){
	*isTimeLimReached = true; 
	goto TERM_SOLVELEAFNODE;
      }
      
      OsiSolverInterface *relaxSolver = new OsiClpSolverInterface();
      relaxSolver->loadProblem(*newMatrix, newColLb, newColUb,
			       objCoeffs, rowLb, rowUb);

      relaxSolver->setObjSense(1);

      relaxSolver->messageHandler()->setLogLevel(0);

      relaxSolver->initialSolve();

      remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getTime();

      //sahar:Fix this by setting the time limit for relax solver
      if(remainingTime < etol){
	*isTimeLimReached = true;
      }
      else{
	if(relaxSolver->isProvenOptimal()){
	  bound = relaxSolver->getObjValue() * relaxSolver->getObjSense();
	}
	else{
	  bound = infinity;
	}
      }
      delete relaxSolver;
    }
    delete newMatrix;
    delete [] rowLb;
    delete [] rowUb;
    delete [] objCoeffs;
  }

 TERM_SOLVELEAFNODE:
  delete [] leafColLb;
  delete [] leafColUb;
  delete [] newColLb;
  delete [] newColUb;
  delete [] upperSol;
  linkSol.clear();

  return bound;
}

//#############################################################################
CoinPackedMatrix*
MibSCutGenerator::getLeafConst(int nodeIndex,int &numRows, std::vector<double> &rowLb,
			       std::vector<double> &rowUb)
{
  
  OsiSolverInterface *hSolver = localModel_->getSolver();    

  int i(0), j(0);
  int rowIndex(0), rowStart(0), rowEnd(0);
  int numCols(localModel_->getUpperDim() + localModel_->getLowerDim());
  int numOrigRows(localModel_->getNumOrigCons());
  int start(localModel_->boundProbLeafNodeCutStarts_[nodeIndex]);
  int numAddedRows(localModel_->boundProbLeafNodeCutStarts_[nodeIndex + 1] - start);
  CoinPackedMatrix origMatrix = *localModel_->origConstCoefMatrix_;
  origMatrix.reverseOrdering(); 

  CoinPackedMatrix *leafConst = new CoinPackedMatrix(false, 0, 0);
  leafConst->setDimensions(0, numCols);
  
  CoinPackedVector row;

  numRows = numOrigRows + numAddedRows;
  for(i = 0; i < numOrigRows; i++){
    leafConst->appendRow(origMatrix.getVector(i));  
  }
  
  for(i = start; i < start + numAddedRows; i++){
    rowIndex = localModel_->boundProbLeafNodeCutInf_[i];
    rowStart = localModel_->boundProbCutPoolStarts_[rowIndex];
    rowEnd = localModel_->boundProbCutPoolStarts_[rowIndex + 1];
    for(j = rowStart; j < rowEnd; j++){
      row.insert(localModel_->boundProbCutPoolIndices_[j], localModel_->boundProbCutPoolValues_[j]);
    }
    leafConst->appendRow(row);
    row.clear();
    rowLb.push_back(localModel_->boundProbCutPoolBounds_[2 * rowIndex]);
    rowUb.push_back(localModel_->boundProbCutPoolBounds_[2 * rowIndex + 1]);   
  }

  return leafConst;

}

//#############################################################################
void
MibSCutGenerator::solveMips(OsiSolverInterface * mipSolver)
{

  std::string feasCheckSolver =
    localModel_->MibSPar_->entry(MibSParams::feasCheckSolver);

  int maxThreadsLL(localModel_->MibSPar_->entry
		   (MibSParams::maxThreadsLL));

  int whichCutsLL(localModel_->MibSPar_->entry
		  (MibSParams::whichCutsLL));

  double timeLimit(localModel_->AlpsPar()->entry(AlpsParams::timeLimit));

  double remainingTime(0.0);

  remainingTime = timeLimit - localModel_->broker_->subTreeTimer().getTime();
  remainingTime = CoinMax(remainingTime, 0.00);

  if (feasCheckSolver == "Cbc"){
        dynamic_cast<OsiCbcSolverInterface *>
	  (mipSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
        sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (mipSolver)->getSymphonyEnvironment();
	sym_set_dbl_param(env, "time_limit", remainingTime);
	sym_set_int_param(env, "do_primal_heuristic", FALSE);
	sym_set_int_param(env, "verbosity", -2);
	sym_set_int_param(env, "prep_level", -1);
	sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
	sym_set_int_param(env, "tighten_root_bounds", FALSE);
	sym_set_int_param(env, "max_sp_size", 100);
	sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	if (whichCutsLL == 0){
	  sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	}else{
	  sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	}
	if (whichCutsLL == 1){
	  sym_set_int_param(env, "generate_cgl_knapsack_cuts",
			    DO_NOT_GENERATE);
	  sym_set_int_param(env, "generate_cgl_probing_cuts",
			    DO_NOT_GENERATE);
	  sym_set_int_param(env, "generate_cgl_clique_cuts",
			    DO_NOT_GENERATE);
	  sym_set_int_param(env, "generate_cgl_twomir_cuts",
			    DO_NOT_GENERATE);
	  sym_set_int_param(env, "generate_cgl_flowcover_cuts",
			    DO_NOT_GENERATE);
	}
#endif
  }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
    mipSolver->setHintParam(OsiDoReducePrint);
    mipSolver->messageHandler()->setLogLevel(0);
        CPXENVptr cpxEnv =
	  dynamic_cast<OsiCpxSolverInterface*>(mipSolver)->getEnvironmentPtr();
	assert(cpxEnv);
	CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
  }

  mipSolver->branchAndBound();
}
  
//#############################################################################
int
MibSCutGenerator::generalNoGoodCut(BcpsConstraintPool &conPool)
{

    /** Add specialized bilevel feasibility cuts, as appropriate **/

    //std::cout << "Generating No-Good Cuts." << std::endl;
    int useLinkingSolutionPool(localModel_->MibSPar_->entry
		   (MibSParams::useLinkingSolutionPool));

    bool allowRemoveCut(localModel_->MibSPar_->entry
			(MibSParams::allowRemoveCut));
    
    OsiSolverInterface * solver = localModel_->solver();

    bool isTimeLimReached(false);
    int numCuts(0);
    double etol(localModel_->etol_);
    const double * sol = solver->getColSolution();
    int uN(localModel_->upperDim_);
    int * upperColInd = localModel_->getUpperColInd();
    int * fixedInd = localModel_->getFixedInd();

    int i(0), index(0);
    double cutub(- 1.0);
    double cutlb(- solver->getInfinity());
    std::vector<int> indexList;
    std::vector<double> valsList;

    MibSBilevel *bS = localModel_->bS_;

    bool shouldFindBestSol(true);

    if(((useLinkingSolutionPool != PARAM_ON) && ((bS->isUBSolved_ == true) ||
				     ((bS->isLowerSolved_ == true) &&
				      (bS->isProvenOptimal_ == false)))) ||
       ((useLinkingSolutionPool == PARAM_ON) &&
	((bS->tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsInfeasible) ||
				     (bS->tagInSeenLinkingPool_ ==
				      MibSLinkingPoolTagUBIsSolved)))){
	shouldFindBestSol = false;
    }

    
    if(shouldFindBestSol == true){
	storeBestSolHypercubeIC(sol, bS->objValVec_, isTimeLimReached);
    }

    if(isTimeLimReached == true){
	goto TERM_GENERALNOGOOD;
    }	

    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    indexList.push_back(index);
	    if(sol[index] > etol){
		valsList.push_back(1.0);
		cutub += 1.0;
	    }
	    else{
		valsList.push_back(- 1.0);
	    }
	}
    }

    assert(indexList.size() == valsList.size());

    numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, allowRemoveCut);

    TERM_GENERALNOGOOD:

    return numCuts;

}
    
//########################################################################### 
double *
MibSCutGenerator::findDeepestLandPCut_ValFunc()
{

  // VERIFY THAT THIS WORKS WITH > 1 UL CONSTRAINT

  /* Construct the CGLP to find the deepest Lift and Project Cut */

  OsiSolverInterface * cglpSolver = new OsiClpSolverInterface();
  //OsiSolverInterface * cglpSolver = new OsiSymSolverInterface();
  OsiSolverInterface * oSolver = localModel_->solver();

  int numOrigRows = localModel_->numOrigCons_;
  int numLowerRows = localModel_->lowerRowNum_;
  int numUpperRows = numOrigRows - numLowerRows;
  const int numCols = oSolver->getNumCols();
  const int numRows = oSolver->getNumRows();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByCol();
  const double * rightHandSide = oSolver->getRightHandSide();
  MibSBilevel *bS = localModel_->bS_;
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  int * lowerRowInd = localModel_->getLowerRowInd();
  double * lowerObjCoeffs = localModel_->getLowerObjCoeffs();

  int numEntries(numCols + 1);
  double * cutVals = new double[numEntries];
  CoinZeroN(cutVals, numEntries);

  int numULCols = localModel_->getUpperDim();
  int numLLCols = localModel_->getLowerDim();
  double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  double * optLowerSol = bS->optLowerSolutionOrd_; // optimal LL solution 

  /* the negative and positive slopes of the LL value function */
  double leftSlope = localModel_->leftSlope_;
  double rightSlope = localModel_->rightSlope_;

  /* 
     For each polyhedron:
     Column for each variable in the LR + gamma 
     + 2 columns for each constraint,including the disjunction;
     Note: vars associated with = constraints are split into 2 vars
     Row for each variable in the LR + 2 gamma constraints + normalization
  */

  int numAuxRows(numRows - numOrigRows);
  //int numAuxRows(0);
  int numCGCols(0);
  numCGCols = numCols + 1+ 4 * numOrigRows + 4 + 2 * numAuxRows;
  int numCGRows(2 * numCols + 2 + 1);
  
  double * colUb = new double[numCGCols];
  double * colLb = new double[numCGCols];
  double * rowUb = new double[numCGRows];
  double * rowLb = new double[numCGRows];
  double * objCoeffs = new double[numCGCols];
  
  /* 
     variables are ordered: alpha, beta, gamma, u,v,w,z,etc 
     alpha, beta, and gamma are free, rest nonnegative
  */

  CoinFillN(colLb, numCols + 1, - cglpSolver->getInfinity());
  CoinFillN(colLb + numCols + 1, (numCGCols - (numCols + 1)), 0.0);

  CoinFillN(colUb, numCGCols, cglpSolver->getInfinity());

  CoinFillN(rowLb, numCGRows - 1, - cglpSolver->getInfinity());
  CoinFillN(rowUb, numCGRows - 1, 0.0);

  /* normalization constraint */
  rowLb[numCGRows - 1] = 1.0;
  rowUb[numCGRows - 1] = 1.0;

  int i(0), j(0);

  CoinZeroN(objCoeffs, numCGCols);

  /* alpha */
  for(i = 0; i < numULCols; i++)
    objCoeffs[i] = upperSol[i];

  /* beta */
  for(i = numULCols; i < numULCols + numLLCols; i++)
    objCoeffs[i] = lowerSol[i - numULCols];

  /* gamma */
  objCoeffs[numULCols + numLLCols] = - 1.0;

  CoinPackedMatrix * cglpMat = new CoinPackedMatrix(false, 0, 0);
  cglpMat->setDimensions(0, numCGCols);

  /* Set up the matrix for the CGLP */

  /* first, the rows constraining alpha */

  int rowIndex(0);
  int colIndex(0), insertIndex1(0), insertIndex2(0);
  int cnt(0);
  double value(0.0);  

  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    row1.insert(i, - 1.0); // alpha coefficient
    row2.insert(i, - 1.0); // alpha coefficient

    for(j = 0; j < numOrigRows; j++){

      rowIndex = j;
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex + cnt;
        
      row1.insert(insertIndex1, value);
      row1.insert(insertIndex1 + 1, - value);

      //changed this 8/27 (sd)
      //insertIndex2 = numULCols + numLLCols + 1 + rowIndex + cnt 
      //+ 2 * numUpperRows + 4 * numLowerRows;
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex + cnt 
	+ 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows;
      
      row2.insert(insertIndex2, value);
      row2.insert(insertIndex2 + 1, - value);
      
      cnt++;
    }

    cnt = 0;

    for(j = 0; j < numLowerRows; j++){
      
      //matIndex = lowerRowInd[j];
      rowIndex = lowerRowInd[j];
      
      value = matrix->getCoefficient(rowIndex, colIndex);

      //changed this 8/27 (sd)
      //insertIndex1 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex - numUpperRows;
      insertIndex1 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + j;
      row1.insert(insertIndex1, value);
      
      insertIndex1 += numLowerRows;
      row1.insert(insertIndex1, - leftSlope * value);
      //row1.insert(insertIndex1, leftSlope * value);
     
      //changed this 8/27 
      //insertIndex2 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex
      //+ 2 * numUpperRows + 4 * numLowerRows - numUpperRows;
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + j
	+ 2 * numUpperRows + 4 * numLowerRows + numAuxRows;
      row2.insert(insertIndex2, - value);
      
      insertIndex2 += numLowerRows;
      row2.insert(insertIndex2, - rightSlope * value); 
      
    }
 
    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      //changed this 8/27 (sd)      
      //insertIndex1 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex 
      //- numUpperRows + numLowerRows;
      insertIndex1 = numULCols + numLLCols + 1 + 
      	2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + rowIndex 
	- numOrigRows;
      
      row1.insert(insertIndex1, - value);

      //changed this 8/27 (sd)      
      //insertIndex2 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex
      //+ 2 * numUpperRows + 4 * numLowerRows 
      //- numUpperRows + numLowerRows;
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + numLowerRows + 1 
	+ 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows 
	+ rowIndex - numOrigRows; 
      
      row2.insert(insertIndex2, - value);
      	
    }
 
    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }

  /* then, the rows constraining beta */

  for(i = 0; i < numLLCols; i++){
    
    colIndex = lowerColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    insertIndex1 = i + numULCols;
    row1.insert(insertIndex1, - 1.0); // beta coefficient

    insertIndex2 = i + numULCols;
    row2.insert(insertIndex1, - 1.0); // beta coefficient

    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
    
      //changed this 8/27 (sd)
      //insertIndex1 = numULCols + numLLCols + 1 
      //+ 2 * numUpperRows + rowIndex - numUpperRows;
      insertIndex1 = numULCols + numLLCols + 1 
	 + 2 * numUpperRows + j;
      
      row1.insert(insertIndex1, value);
      row1.insert(insertIndex1 + 1, - value);
      
      //changed this 8/27 (sd)
      //insertIndex2 = numULCols + numLLCols + 1 + rowIndex - numUpperRows 
      //+ 2 * numUpperRows + 2 * numUpperRows + 4 * numLowerRows;
      insertIndex2 = numULCols + numLLCols + 1 + j
	+ 4 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows;
      
      row2.insert(insertIndex2, value);
      row2.insert(insertIndex2 + 1, - value);
    }
    
    
    value = lowerObjCoeffs[i];
    
    insertIndex1 = numULCols + numLLCols + 1 + 
      2 * (numUpperRows + numLowerRows) + numLowerRows;
    
    row1.insert(insertIndex1, - value);
    
    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + 
    //4 * (numUpperRows + numLowerRows) + 3 * numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + 
      4 * (numUpperRows + numLowerRows) + 2 * numLowerRows + 1 + numAuxRows;

    row2.insert(insertIndex2, - value);

    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);
    	
      //changed this 8/27 (sd)
      //insertIndex1 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + rowIndex 
      //- numUpperRows;
      insertIndex1 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + rowIndex 
	- numOrigRows;
      
      row1.insert(insertIndex1, - value);
      
      //changed this 8/27 (sd)
      //insertIndex2 = numULCols + numLLCols + 1 + 
      //4 * (numUpperRows + numLowerRows) + 
      //3 * numLowerRows + 1 + rowIndex - numUpperRows;
      insertIndex2 = numULCols + numLLCols + 1 + 
	4 * (numUpperRows + numLowerRows) + 
	2 * numLowerRows + 2 + rowIndex - numOrigRows + numAuxRows;
            
      row2.insert(insertIndex2, - value);
      
    }

    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);
    
  }
  
  
  /* then, the rows constraining gamma */

  CoinPackedVector row1;
  CoinPackedVector row2;

  insertIndex1 = numULCols + numLLCols;   
  row1.insert(insertIndex1, 1.0); // gamma coefficient
  
  insertIndex2 = numULCols + numLLCols;
  row2.insert(insertIndex2, 1.0); // gammaa coefficient
  
  for(j = 0; j < numUpperRows; j++){
    
    value = rightHandSide[j];

    //changed this 8/27 (sd)
    //insertIndex1 = numULCols + numLLCols + 1 + j;
    insertIndex1 = numULCols + numLLCols + 1 + j + cnt;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + j
    //+ 2 * numUpperRows + 4 * numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + j
      + 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows + cnt;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    
    cnt++;
  }

  cnt = 0;
  
  for(j = 0; j < numLowerRows; j++){
    
    value = rightHandSide[numUpperRows + j];
   
    //changed this 8/27 (sd)
    //insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j;
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j + cnt;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
    //+ 2 * numUpperRows + 4 * numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows + cnt;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    cnt++;    
  }

  cnt = 0;

  double mult1(0.0), mult2(0.0);
  
  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];
    
    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
      mult1 += value * upperSol[i];
    }
  }
  
  for(j = 0; j < numLLCols; j++){

    value = lowerObjCoeffs[j];
    mult2 += value * optLowerSol[j];

  }
    
  insertIndex1 = numULCols + numLLCols + 1 + 2 * (numUpperRows + numLowerRows); 
  
  row1.insert(insertIndex1, - mult1);
  
  insertIndex1 += numLowerRows;
  
  //row1.insert(insertIndex1, leftSlope * mult1 - mult2); 
  row1.insert(insertIndex1, leftSlope * mult1 + mult2); 
  
  //changed this 8/27 (sd)  
  //insertIndex2 = numULCols + numLLCols + 1 
  //+ 4 * (numUpperRows + numLowerRows) + 2 * numLowerRows; 
  insertIndex2 = numULCols + numLLCols + 1 
    + 4 * (numUpperRows + numLowerRows) + numLowerRows + 1 + numAuxRows; 
  
  row2.insert(insertIndex2, mult1);
  
  insertIndex2 += numLowerRows;
  
  //row2.insert(insertIndex2, rightSlope * mult1 - mult2); 
  row2.insert(insertIndex2, rightSlope * mult1 + mult2); 
  
  for(j = numOrigRows; j < numRows; j++){
    
    //value = rightHandSide[j];
    value = - rightHandSide[j];
  
    //changed this 8/27 (sd)
    //insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
    //+ 2 * numLowerRows + numUpperRows + numLowerRows; 
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows 
      + 2 * numLowerRows + numUpperRows + 1 + j - numOrigRows; 
    
    row1.insert(insertIndex1, - value);

    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
    //+ 2 * numUpperRows + 4 * numLowerRows
    //+ 2 * numLowerRows + numUpperRows + numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + 4 * numUpperRows
      + 6 * numLowerRows + 2 + j - numOrigRows + numAuxRows;

    row2.insert(insertIndex2, - value);

  }
  
  cglpMat->appendRow(row1);
  cglpMat->appendRow(row2);

  /* then a normalization constraint */

  CoinPackedVector row;

  for(i = numULCols + numLLCols + 1; i < numCGCols; i++)
    row.insert(i, 1.0);

  cglpMat->appendRow(row);

  cglpSolver->assignProblem(cglpMat, colLb, colUb, objCoeffs, rowLb, rowUb);

  if(0)
    cglpSolver->writeLp("cglp");

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiClpSolverInterface *> 
     (cglpSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("max_active_nodes", 1);
#endif
  cglpSolver->initialSolve();
  
  if(cglpSolver->isProvenOptimal()){

    const double * cglpSol = cglpSolver->getColSolution();
    double cutViolation = cglpSolver->getObjValue();

    for(i = 0; i < numULCols; i++){
      cutVals[upperColInd[i]] = cglpSol[i];
      if(0)
	std::cout << upperColInd[i] << ": " << cglpSol[i] << std::endl;
    }

    for(i = 0; i < numLLCols; i++){
      cutVals[lowerColInd[i]] = cglpSol[numULCols + i];
      if(0){
	std::cout << lowerColInd[i] << ": " 
		  << cglpSol[numULCols + i] << std::endl;
      }
    }
    
    cutVals[numULCols + numLLCols] = cglpSol[numULCols + numLLCols];

    std::cout << "Violated inequality found, with violation of " 
	      << cutViolation << "." << std::endl;

  }
  else{

    std::cout << "CGLP failed to find a violated inequality." << std::endl;
    
  }

  delete cglpSolver;

  /* Return the entries of the deepest cut - alpha * x + beta * y >= gamma */
  return cutVals;

}

//#############################################################################
double *
MibSCutGenerator::findDeepestLandPCut_IncObj(double * upperSol, 
					     double * lowerSol,
					     double * optLowerSol)
{

  //NOTE: THIS IS ONLY CASE CASE OF THE CGLP (I.E. A^2 POS AS IN MIPINT)
  //      THIS ASSUMES INEQUALITY CONSTRAINTS
  //      SHOULD ADD A CHECK FOR THAT AND INTRODUCE +,- VARS FOR =-CONSTRAINTS

  //PROBABLY WANT TO ADD AN ARGUMENT OF UPPERSOL
  //THAT WAY CAN CALL THE SAME FUNCTION FOR BOTH THE CURRENT SOL
  //AND THE MAXIMAL VECTORS

  /* Construct the CGLP to find the deepest Lift and Project Cut */

  OsiSolverInterface * cglpSolver = new OsiClpSolverInterface();
  OsiSolverInterface * oSolver = localModel_->solver();

  int numOrigRows = localModel_->numOrigCons_;
  int numLowerRows = localModel_->lowerRowNum_;
  int numUpperRows = numOrigRows - numLowerRows;
  const int numCols = oSolver->getNumCols();
  const int numRows = oSolver->getNumRows();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByCol();
  //const CoinPackedMatrix * matrix = oSolver->getMatrixByRow();
  const double * rightHandSide = oSolver->getRightHandSide();
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  int * lowerRowInd = localModel_->getLowerRowInd();
  int * upperRowInd = localModel_->getUpperRowInd();
  double * lowerObjCoeffs = localModel_->getLowerObjCoeffs();
  double lowerObjSense = localModel_->getLowerObjSense();
  double etol = localModel_->etol_;

  int numEntries(numCols + 1);
  double * cutVals = new double[numEntries];
  CoinZeroN(cutVals, numEntries);

  int numULCols = localModel_->getUpperDim();
  int numLLCols = localModel_->getLowerDim();

  //PROBABLY WANT TO CHANGE THIS TO BE AN ARGUMENT
  //double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  //double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  //double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

  /* 
     Column for each variable in the LR + gamma 
     + 1 column for each inequality constraint,including the disjunction
     + 2 columns for each equality constraint;
     Note: vars associated with = constraints are split into 2 vars
     Row for each variable in the LR + 2 gamma constraints + normalization
  */

  int numCGCols(0);
  int numAuxRows(numRows - numOrigRows);
  numCGCols = numCols + 1 + 2 * numOrigRows + 4 + 2 * numAuxRows;
  int numCGRows(2 * numCols + 2 + 1);
  

  double * colUb = new double[numCGCols];
  double * colLb = new double[numCGCols];
  double * rowUb = new double[numCGRows];
  double * rowLb = new double[numCGRows];
  double * objCoeffs = new double[numCGCols];
  
  /* 
     variables are ordered: alpha, beta, gamma, u,v,w,z,etc 
     alpha, beta, and gamma are free, rest nonnegative
  */

  CoinFillN(colLb, numCols + 1, - cglpSolver->getInfinity());
  CoinFillN(colLb + numCols + 1, (numCGCols - (numCols + 1)), 0.0);

  CoinFillN(colUb, numCGCols, cglpSolver->getInfinity());

  CoinFillN(rowLb, numCGRows - 1, - cglpSolver->getInfinity());
  CoinFillN(rowUb, numCGRows - 1, 0.0);

  /* normalization constraint */
  rowLb[numCGRows - 1] = 1.0;
  rowUb[numCGRows - 1] = 1.0;

  int i(0), j(0);

  CoinZeroN(objCoeffs, numCGCols);

  /* alpha */
  for(i = 0; i < numULCols; i++)
    objCoeffs[i] = upperSol[i];

  /* beta */
  for(i = numULCols; i < numULCols + numLLCols; i++)
    objCoeffs[i] = lowerSol[i - numULCols];

  /* gamma */
  objCoeffs[numULCols + numLLCols] = - 1.0;

  CoinPackedMatrix * cglpMat = new CoinPackedMatrix(false, 0, 0);
  cglpMat->setDimensions(0, numCGCols);

  /* Set up the matrix for the CGLP */

  /* first, the rows constraining alpha */

  int rowIndex(0);
  int colIndex(0), insertIndex1(0), insertIndex2(0);
  double value(0.0);  

  //matrix->dumpMatrix();

  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    row1.insert(i, 1.0); // alpha coefficient
    row2.insert(i, 1.0); // alpha coefficient
    /*
    for(j = 0; j < numOrigRows; j++){

      rowIndex = j;
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex;
      
      row1.insert(insertIndex1, - value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex 
	+ numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
      
    }
    */
    for(j = 0; j < numUpperRows; j++){

      rowIndex = upperRowInd[j];
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex;
      
      row1.insert(insertIndex1, - value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex 
	+ numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
      
    }

    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex;
      
      row1.insert(insertIndex1, -value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex 
	+ numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, -value);
      //row2.insert(insertIndex2 + 1, value);
      
    }

    //cnt = 0;

    //FIXME: SHOULD BE 1 IF < ETOL!
    //w variable
    if(upperSol[i] > etol)
      value = 0.0;
    else 
      value = 1.0;

    insertIndex1 = numULCols + numLLCols + 1 + 
      numUpperRows + numLowerRows;
    row1.insert(insertIndex1, - value);

    insertIndex2 = numULCols + numLLCols + 1 
      + 2 * (numLowerRows + numUpperRows) + 2 + numAuxRows;

    row2.insert(insertIndex2, value);

    //multipliers associated with previously-added cuts
    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);

      insertIndex1 = numULCols + numLLCols + 1 + 
	numUpperRows + numLowerRows + 2 + (rowIndex - numOrigRows);
      row1.insert(insertIndex1, -value);
      
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows + 2) 
	+ (rowIndex - numOrigRows) + numAuxRows; 
      row2.insert(insertIndex2, -value);
      	
    }
        
    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }


  /* then, the rows constraining beta */

  for(i = 0; i < numLLCols; i++){
    
    colIndex = lowerColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    insertIndex1 = i + numULCols;
    row1.insert(insertIndex1, 1.0); // beta coefficient

    insertIndex2 = i + numULCols;
    row2.insert(insertIndex1, 1.0); // beta coefficient

    for(j = 0; j < numUpperRows; j++){

      rowIndex = upperRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
    
      insertIndex1 = numULCols + numLLCols + 1 
	+ rowIndex;
      
      row1.insert(insertIndex1, - value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 
	+ rowIndex + numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
    }

    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
    
      insertIndex1 = numULCols + numLLCols + 1 
	+ rowIndex;
      
      row1.insert(insertIndex1, -value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 
	+ rowIndex + numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
    }
    
    value = lowerObjCoeffs[i] * lowerObjSense;
    
    insertIndex1 = numULCols + numLLCols + 1 + 
      numUpperRows + numLowerRows + 1;

    row1.insert(insertIndex1, -value);

    //insertIndex2 = numULCols + numLLCols + 1 + 
    //4 * (numUpperRows + numLowerRows) + 3 * numLowerRows;

    //row2.insert(insertIndex2, - value);

    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);
	
      insertIndex1 = numULCols + numLLCols + 1 + 
	numUpperRows + numLowerRows + 2 + (rowIndex - numOrigRows);
      
      row1.insert(insertIndex1, - value);
      
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows + 2)
	+ (rowIndex - numOrigRows) + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      
    }
    
    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);
    
  }

  /* then, the rows constraining gamma */

  CoinPackedVector row1;
  CoinPackedVector row2;

  insertIndex1 = numULCols + numLLCols;   
  row1.insert(insertIndex1, - 1.0); // gamma coefficient
  
  insertIndex2 = numULCols + numLLCols;
  row2.insert(insertIndex2, - 1.0); // gammaa coefficient
  
  for(j = 0; j < numUpperRows; j++){
    
    //value = rightHandSide[j];
    value = rightHandSide[upperRowInd[j]];

    insertIndex1 = numULCols + numLLCols + 1 + j;

    row1.insert(insertIndex1, value);
    //row1.insert(insertIndex1 + 1, - value);

    insertIndex2 = numULCols + numLLCols + 1 + j
      + numUpperRows + numLowerRows + 2 + numAuxRows;

    row2.insert(insertIndex2, value);
    //row2.insert(insertIndex2 + 1, - value);
    
  }
  
  //cnt = 0;

  for(j = 0; j < numLowerRows; j++){
    
    //value = rightHandSide[numUpperRows + j];
    value = rightHandSide[lowerRowInd[j]];
   
    insertIndex1 = numULCols + numLLCols + 1 + numUpperRows + j;

    row1.insert(insertIndex1, value);
    //row1.insert(insertIndex1 + 1, - value);

    insertIndex2 = numULCols + numLLCols + 1 
      + numUpperRows + numLowerRows + 2 + numUpperRows
      + numAuxRows + j;

    row2.insert(insertIndex2, value);
    //row2.insert(insertIndex2 + 1, - value);
    
  }

  double mult2(0.0);

  for(j = 0; j < numLLCols; j++){

    value = lowerObjCoeffs[j] * lowerObjSense;
    mult2 += value * optLowerSol[j];

  }

  insertIndex1 = numULCols + numLLCols + 1 + numUpperRows + numLowerRows + 1; 
  
  row1.insert(insertIndex1, mult2); 
  
  insertIndex2 = numULCols + numLLCols + 1 
    + 2 * (numUpperRows + numLowerRows) + 2 + numAuxRows;

  //row2.insert(insertIndex2, value);
  row2.insert(insertIndex2, - 1);

  for(j = numOrigRows; j < numRows; j++){
    
    value = rightHandSide[j];
  
    insertIndex1 = numULCols + numLLCols + 1 
      + numUpperRows + numLowerRows + 2 + (j - numOrigRows);
    
    row1.insert(insertIndex1, value);

    insertIndex2 = numULCols + numLLCols + 1 
      + 2 * (numUpperRows + numLowerRows + 2) 
      + (j - numOrigRows) + numAuxRows;

    row2.insert(insertIndex2, value);

  }

  cglpMat->appendRow(row1);
  cglpMat->appendRow(row2);

  /* then a normalization constraint */

  CoinPackedVector row;

  for(i = numULCols + numLLCols + 1; i < numCGCols; i++)
    row.insert(i, 1.0);

  cglpMat->appendRow(row);

  cglpSolver->assignProblem(cglpMat, colLb, colUb, objCoeffs, rowLb, rowUb);

  cglpSolver->setObjSense(-1); //maximization

  if(0)
    cglpSolver->writeLp("cglp");
  //cglpSolver->writeMps("cglp");

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiClpSolverInterface *> 
     (cglpSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("max_active_nodes", 1);
#endif

  cglpSolver->initialSolve();

  if(cglpSolver->isProvenOptimal()){
    
    const double * cglpSol = cglpSolver->getColSolution();
    double cutViolation = cglpSolver->getObjValue();

    for(i = 0; i < numULCols; i++){
      cutVals[upperColInd[i]] = cglpSol[i];
      if(0)
	std::cout << upperColInd[i] << ": " << cglpSol[i] << std::endl;
    }
      
    for(i = 0; i < numLLCols; i++){
      cutVals[lowerColInd[i]] = cglpSol[numULCols + i];
      if(0){
	std::cout << lowerColInd[i] << ": " 
		  << cglpSol[numULCols + i] << std::endl;
      }
    }

    cutVals[numULCols + numLLCols] = cglpSol[numULCols + numLLCols];

    std::cout << "Violated inequality found, with violation of " 
	      << cutViolation << "." << std::endl;

  }
  else{

    std::cout << "CGLP failed to find a violated inequality." << std::endl;
    
  }

  delete cglpSolver;

  /* Return the entries of the deepest cut - alpha * x + beta * y >= gamma */
  return cutVals;

}
//#############################################################################
double *
MibSCutGenerator::findDeepestLandPCut1()
{
  // OLD FUNCTION.   
 
  /* Construct the CGLP to find the deepest Lift and Project Cut */

  OsiSolverInterface * cglpSolver = new OsiClpSolverInterface();
  //OsiSolverInterface * cglpSolver = new OsiCbcSolverInterface();
  //OsiSolverInterface * cglpSolver = new OsiSymSolverInterface();
  OsiSolverInterface * oSolver = localModel_->solver();

  int numOrigRows = localModel_->numOrigCons_;
  int numLowerRows = localModel_->lowerRowNum_;
  int numUpperRows = numOrigRows - numLowerRows;
  const int numCols = oSolver->getNumCols();
  const int numRows = oSolver->getNumRows();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByCol();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();
  const double * rightHandSide = oSolver->getRightHandSide();
  MibSBilevel *bS = localModel_->bS_;
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  double * lowerObjCoeffs = localModel_->getLowerObjCoeffs();

  int numEntries(numCols + 1);
  double * cutVals = new double[numEntries];
  CoinZeroN(cutVals, numEntries);

  int numULCols = localModel_->getUpperDim();
  int numLLCols = localModel_->getLowerDim();
  double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  double * optLowerSol = bS->optLowerSolutionOrd_; // optimal LL solution 

  /* the negative and positive slopes of the LL value function */
  double leftSlope = localModel_->leftSlope_;
  double rightSlope = localModel_->rightSlope_;

  /* 
     Column for each variable in the LR + gamma 
     + 2 columns for each constraint,including the disjunction;
     Note: vars associated with = constraints are split into 2 vars
     Row for each variable in the LR + 2 gamma constraints + normalization
  */

  int numCGCols(0);
  numCGCols = numCols + 1 + 4 * numOrigRows + 2 * (numRows - numOrigRows + 2);
  int numCGRows(2 * numCols + 2 + 1);
  
  int numAuxRows(numRows - numOrigRows);

  double * colUb = new double[numCGCols];
  double * colLb = new double[numCGCols];
  double * rowUb = new double[numCGRows];
  double * rowLb = new double[numCGRows];
  double * objCoeffs = new double[numCGCols];
  
  /* 
     variables are ordered: alpha, beta, gamma, u,v,w,z,etc 
     alpha, beta, and gamma are free, rest nonnegative
  */

  CoinFillN(colLb, numCols + 1, - cglpSolver->getInfinity());
  CoinFillN(colLb + numCols + 1, (numCGCols - (numCols + 1)), 0.0);

  CoinFillN(colUb, numCGCols, cglpSolver->getInfinity());

  CoinFillN(rowLb, numCGRows - 1, - cglpSolver->getInfinity());
  CoinFillN(rowUb, numCGRows - 1, 0.0);

  /* normalization constraint */
  rowLb[numCGRows - 1] = 1.0;
  rowUb[numCGRows - 1] = 1.0;

  int i(0), j(0), start(0), end(0);

  CoinZeroN(objCoeffs, numCGCols);

  /* alpha */
  for(i = 0; i < numULCols; i++)
    objCoeffs[i] = upperSol[i];

  /* beta */
  for(i = numULCols; i < numULCols + numLLCols; i++)
    objCoeffs[i] = lowerSol[i - numULCols];

  /* gamma */
  objCoeffs[numULCols + numLLCols] = - 1.0;

  CoinPackedMatrix * cglpMat = new CoinPackedMatrix(false, 0, 0);
  cglpMat->setDimensions(0, numCGCols);

  /* Set up the matrix for the CGLP */

  /* first, the rows constraining alpha */

  int colIndex(0), matIndex(0), insertIndex1(0), insertIndex2(0);
  int size(0), cnt(0);
  double value(0.0);  

  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    row1.insert(i, - 1.0); // alpha coefficient
    row2.insert(i, - 1.0); // alpha coefficient

    size = matrix->getVectorSize(colIndex);

    start = matStarts[colIndex];

    end = start + size;

    for(j = start; j < end; j++){

      matIndex = matIndices[j];
      if(matIndex < numOrigRows){
	
	value = matElements[matIndex + start];
	
	insertIndex1 = numULCols + numLLCols + 1 + matIndex + cnt;

	row1.insert(insertIndex1, value);
	row1.insert(insertIndex1 + 1, - value);

	insertIndex2 = numULCols + numLLCols + 1 + matIndex + cnt 
	  + 2 * numUpperRows + 4 * numLowerRows;

	row2.insert(insertIndex2, value);
	row2.insert(insertIndex2 + 1, - value);

	cnt++;
      }
    }
    cnt = 0;

    //FIXME:  THERE IS A BUG HERE
    // NEED TO ADD IN THE NUMBER OF CUTS THAT HAVE BEEN ADDED
    // TO THE INDICES BELOW
    // RIGHT NOW JUST USING NUM ORIGINAL ROWS
    // INDEX IS OFF
    // 8/26: IS THIS STILL TRUE?  SEEMS TO BE WORKING...

    for(j = start; j < end; j++){

       matIndex = matIndices[j];
       if((matIndex > numUpperRows - 1) && (matIndex < numOrigRows)){

	  value = matElements[matIndex + start];

	  insertIndex1 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex - numUpperRows;
	  row1.insert(insertIndex1, value);
	  
	  insertIndex1 += numLowerRows;
	  row1.insert(insertIndex1, - leftSlope * value);

	  insertIndex2 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex
	    + 2 * numUpperRows + 4 * numLowerRows - numUpperRows;
	  row2.insert(insertIndex2, - value);

	  insertIndex2 += numLowerRows;
	  row2.insert(insertIndex2, - rightSlope * value); 
       }
    }

    for(j = start; j < end; j++){

       matIndex = matIndices[j];
       if(matIndex > numOrigRows - 1){
	  value = matElements[matIndex + start];
	  
	  insertIndex1 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex 
	    - numUpperRows + numLowerRows;
	  row1.insert(insertIndex1, - value);

	  // row1.insert(numULCols + numLLCols + 1 +  
	  //     2 * (numUpperRows + numLowerRows) 
	  //      + numLowerRows + index2 - numUpperRows + numLowerRows,
	  //     - leftSlope * value);
 
	  insertIndex2 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex
	    + 2 * numUpperRows + 4 * numLowerRows 
	    - numUpperRows + numLowerRows;
	  row2.insert(insertIndex2, - value);

	  //row2.insert(numULCols + numLLCols + 1 +  
	  //     2 * (numUpperRows + numLowerRows) + numLowerRows 
	  //      + index2 + 2 * numUpperRows + 4 * numLowerRows 
	  //      - numUpperRows + numLowerRows, - rightSlope * value); 
       }
    }


    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }

  /* then, the rows constraining beta */

  for(i = 0; i < numLLCols; i++){
    
    colIndex = lowerColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    insertIndex1 = i + numULCols;
    row1.insert(insertIndex1, - 1.0); // beta coefficient

    insertIndex2 = i + numULCols;
    row2.insert(insertIndex1, - 1.0); // beta coefficient

    size = matrix->getVectorSize(colIndex);

    start = matStarts[colIndex];

    end = start + size;

    for(j = start; j < end; j++){

      matIndex = matIndices[j];

      //should this be searching lowerRowIndices?
      if(matIndex < numOrigRows){
	
	value = matElements[matIndex + start - numUpperRows];
	
	insertIndex1 = numULCols + numLLCols + 1 
	  + 2 * numUpperRows + matIndex - numUpperRows;

	row1.insert(insertIndex1, value);
	row1.insert(insertIndex1 + 1, - value);
	
	insertIndex2 = numULCols + numLLCols + 1 + matIndex - numUpperRows 
	  + 2 * numUpperRows + 2 * numUpperRows + 4 * numLowerRows;

	row2.insert(insertIndex2, value);
	row2.insert(insertIndex2 + 1, - value);
      }
    }

    
    value = lowerObjCoeffs[i];
    
    insertIndex1 = numULCols + numLLCols + 1 + 
      2 * (numUpperRows + numLowerRows) + numLowerRows;

    row1.insert(insertIndex1, - value);

    insertIndex2 = numULCols + numLLCols + 1 + 
      4 * (numUpperRows + numLowerRows) + 3 * numLowerRows;

    row2.insert(insertIndex2, - value);

    for(j = start; j < end; j++){

      matIndex = matIndices[j];

      if(matIndex > numOrigRows - 1){
	
	value = matElements[matIndex + start - numUpperRows];
	//row1.insert(numULCols + numLLCols + 1 + 2 * numUpperRows 
	//	    + index2 - numUpperRows + 1, value);
	
	insertIndex1 = numULCols + numLLCols + 1 + 
	  2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + matIndex 
	  - numUpperRows;
	
	row1.insert(insertIndex1, - value);

	//row1.insert(numULCols + numLLCols + 1 + index2 + 2 * numUpperRows + 1
	//    - numUpperRows + 1, - value);
	//	row2.insert(numULCols + numLLCols + 1 + index2 - numUpperRows 
	//    + 2 * numUpperRows + 2 * numUpperRows 
	//    + 4 * numLowerRows + 1, value);
	
	insertIndex2 = numULCols + numLLCols + 1 + 
	  4 * (numUpperRows + numLowerRows) + 
	  3 * numLowerRows + 1 + matIndex - numUpperRows;

	row2.insert(insertIndex2, - value);
    
	//row2.insert(numULCols + numLLCols + 1 + index2  - numUpperRows 
	//    + 2 * numUpperRows + 1 + 2 * numUpperRows 
	//    + 4 * numLowerRows + 1, - value);
      }
    }

    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }

  
  /* then, the rows constraining gamma */

  CoinPackedVector row1;
  CoinPackedVector row2;

  insertIndex1 = numULCols + numLLCols;   
  row1.insert(insertIndex1, 1.0); // gamma coefficient
  
  insertIndex2 = numULCols + numLLCols;
  row2.insert(insertIndex2, 1.0); // gammaa coefficient
  
  for(j = 0; j < numUpperRows; j++){
    
    value = rightHandSide[j];

    insertIndex1 = numULCols + numLLCols + 1 + j;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    insertIndex2 = numULCols + numLLCols + 1 + j
      + 2 * numUpperRows + 4 * numLowerRows;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    
  }
  
  for(j = 0; j < numLowerRows; j++){
    
    value = rightHandSide[numUpperRows + j];
   
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numUpperRows + 4 * numLowerRows;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    
  }

  double mult1(0.0), mult2(0.0);
  
  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];
    
    size = matrix->getVectorSize(colIndex);
    
    start = matStarts[colIndex];
    
    end = start + size;
    
    for(j = start; j < end; j++){
      
      matIndex = matIndices[j];
      if((matIndex > numUpperRows - 1) && (matIndex < numOrigRows)){
	value = matElements[matIndex + start];
	mult1 += value * upperSol[i];
	
      }
    }
  }
  
  for(j = 0; j < numLLCols; j++){
      
    //index2 = lowerColInd[j];
    value = lowerObjCoeffs[j];
    mult2 += value * optLowerSol[j];
  }
    

  insertIndex1 = numULCols + numLLCols + 1 + 2 * (numUpperRows + numLowerRows); 
  
  row1.insert(insertIndex1, - mult1);
  
  insertIndex1 += numLowerRows;
  
  row1.insert(insertIndex1, leftSlope * mult1 - mult2); 
  
  insertIndex2 = numULCols + numLLCols + 1 
    + 4 * (numUpperRows + numLowerRows) + 2 * numLowerRows; 

  row2.insert(insertIndex2, mult1);

  insertIndex2 += numLowerRows;

  row2.insert(insertIndex2, rightSlope * mult1 - mult2); 
  
  for(j = 0; j < numAuxRows; j++){
    
    value = rightHandSide[numOrigRows + j];
  
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numLowerRows + numUpperRows + numLowerRows; 
    
    row1.insert(insertIndex1, - value);

    //row1.insert(numULCols + numLLCols + 1 + 2 * numUpperRows + j + 1, value);
    
    insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numUpperRows + 4 * numLowerRows
      + 2 * numLowerRows + numUpperRows + numLowerRows;

    row2.insert(insertIndex2, - value);

    //row2.insert(numULCols + numLLCols + 1 + 2 * numUpperRows + j + 1
    //	+ 2 * numUpperRows + 4 * numLowerRows, value);
    
  }


  
  cglpMat->appendRow(row1);
  cglpMat->appendRow(row2);

  /* then a normalization constraint */

  CoinPackedVector row;

  for(i = numULCols + numLLCols + 1; i < numCGCols; i++)
    row.insert(i, 1.0);

  cglpMat->appendRow(row);

  cglpSolver->assignProblem(cglpMat, colLb, colUb, objCoeffs, rowLb, rowUb);

  if(0)
    cglpSolver->writeLp("cglp");
  //cglpSolver->writeMps("cglp");

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiClpSolverInterface *> 
     (cglpSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (cglpSolver)->setSymParam("max_active_nodes", 1);
#endif

  cglpSolver->initialSolve();
  
  if(cglpSolver->isProvenOptimal()){

    const double * cglpSol = cglpSolver->getColSolution();
    double cutViolation = cglpSolver->getObjValue();
    /*
    for(i = 0; i < numEntries; i++)
      cutVals[i] = cglpSol[i];
    */

    for(i = 0; i < numULCols; i++)
      cutVals[upperColInd[i]] = cglpSol[i];

    for(i = 0; i < numLLCols; i++)
      cutVals[lowerColInd[i]] = cglpSol[numULCols + i];

    cutVals[numULCols + numLLCols] = cglpSol[numULCols + numLLCols];

    std::cout << "Violated inequality found, with violation of " 
	      << cutViolation << "." << std::endl;

  }
  else{

    std::cout << "CGLP failed to find a violated inequality." << std::endl;
    
  }

  delete cglpSolver;

  /* Return the entries of the deepest cut - alpha * x + beta * y >= gamma */
  return cutVals;

}

//#############################################################################
int
MibSCutGenerator::interdictionCuts(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIPINT Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();
  MibSBilevel *bS = localModel_->bS_;

  int useBendersCut = 
     localModel_->MibSPar_->entry(MibSParams::useBendersCut);

  int numCuts(0);
  int lN(localModel_->getLowerDim());
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  int uN(localModel_->upperDim_);
  std::vector<int> indexList;
  std::vector<double> valsList;
  double cutub(- 1.0);
  double cutlb(- solver->getInfinity());

  if (bS->isUpperIntegral_){
     MibSTreeNode * node = 
	dynamic_cast<MibSTreeNode *>(localModel_->activeNode_); 
     double etol(localModel_->etol_);
     const double * sol = solver->getColSolution();
     int * upperColInd = localModel_->getUpperColInd();
     
     int i(0), index(0);
     
     for(i = 0; i < uN; i++){
	index = upperColInd[i];
	indexList.push_back(index);
	if(sol[index] > etol){
	   valsList.push_back(1.0);
	   cutub += 1.0;
	}
	else{
	   valsList.push_back(- 1.0);
	}
     }

#if 0
     indexList.push_back(uN + lN);
     valsList.push_back(- 1.0);
#endif
     
     assert(indexList.size() == valsList.size());
     
     numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, false);
     
     indexList.clear();
     valsList.clear();
  }

#if 0
  double maxLowerObj(node->getLowerUB());
  double m(solver->getObjValue() - maxLowerObj);
  //double m(maxLowerObj - 
  cutub = m + maxLowerObj;

  for(i = 0; i < lN; i++){
    index = lowerColInd[i];
    if(fabs(lObjCoeffs[i]) != 0.0){
      indexList.push_back(index);
      valsList.push_back(lObjCoeffs[i]);
    }
  }
  
  if(fabs(m) != 0.0){
    indexList.push_back(uN + lN);
    valsList.push_back(m);
  }

  //for(i = 0; i < valsList.size(); i++){
  //  std::cout << "index: " << indexList.at(i) << std::endl;
  //  std::cout << "value: " << valsList.at(i);
  //}

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, false);
#endif

  if (useBendersCut){
      numCuts += bendersInterdictionMultipleCuts(conPool);
  }

  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::generalWeakIncObjCutCurrent(BcpsConstraintPool &conPool)
{
    bool allowRemoveCut(localModel_->MibSPar_->entry
			(MibSParams::allowRemoveCut));
    
    OsiSolverInterface *solver = localModel_->solver();
    int i;
    int index(0);
    int numCuts(0);
    double cutub(0.0), cutlb(-solver->getInfinity()), bigM(0.0);
    double lObjVal(localModel_->bS_->objValVec_[0]);
    //double bigM(10000);
    double etol(localModel_->etol_);
    int uN(localModel_->getUpperDim());
    int lN(localModel_->getLowerDim());
    int *fixedInd(localModel_->getFixedInd());
    int *upperColInd(localModel_->getUpperColInd());
    int *lowerColInd(localModel_->getLowerColInd());
    double *lObjCoeffs(localModel_->getLowerObjCoeffs());
    double lowerObjSense(localModel_->getLowerObjSense());
    const double *sol(solver->getColSolution());

    std::vector<int> indexList;
    std::vector<double> valsList;

    if(!isBigMIncObjSet_){
    bigMIncObj_ = findBigMIncObjCut();
    isBigMIncObjSet_ = true;
    }

    bigM = bigMIncObj_ - lObjVal + 1;

    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if((fixedInd[index] == 1) && (sol[index] < etol)){
	    indexList.push_back(index);
	    valsList.push_back(-bigM);
	}
    }

    for(i = 0; i < lN; i++){
	index = lowerColInd[i];
	if(fabs(lObjCoeffs[i]) > etol){
	    indexList.push_back(index);
	    valsList.push_back(lowerObjSense * lObjCoeffs[i]);
	}
    }

    cutub = lObjVal;

    assert(indexList.size() == valsList.size());
    numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, allowRemoveCut);

    return numCuts;

}

//#############################################################################
double
MibSCutGenerator::findBigMIncObjCut()
{

    std::string feasCheckSolver =
	localModel_->MibSPar_->entry(MibSParams::feasCheckSolver);

    int maxThreadsLL = localModel_->MibSPar_->entry
	(MibSParams::maxThreadsLL);

    int whichCutsLL = localModel_->MibSPar_->entry
		    (MibSParams::whichCutsLL);
    
    OsiSolverInterface *oSolver = localModel_->solver();

    int i(0);
    int intCnt(0), colIndex(0);
    double bigM(0.0);
    int colNum(localModel_->getNumOrigVars());
    int lCols(localModel_->getLowerDim());
    double lObjSense(localModel_->getLowerObjSense());
    int *lowerColInd(localModel_->getLowerColInd());
    double *lObjCoeffs(localModel_->getLowerObjCoeffs());
    
    OsiSolverInterface * nSolver;
    double *objCoeffs = new double[colNum];
    int *integerVars = new int[colNum];
    CoinFillN(objCoeffs, colNum, 0.0);
    CoinFillN(integerVars, colNum, 0);

    if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"findBigMIncObjCut", "MibSCutGenerator");
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"findBigMIncObjCut", "MibSCutGenerator");
#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"findBigMIncObjCut", "MibSCutGenerator");
    }

    CoinPackedMatrix matrix = *localModel_->origConstCoefMatrix_;

    for(i = 0; i < colNum; i++){
	if(oSolver->isInteger(i)){
	    integerVars[intCnt] = i;
	    intCnt ++;
	}
    }

    for(i = 0; i < lCols; i++){
	colIndex = lowerColInd[i];
	objCoeffs[colIndex] = lObjCoeffs[i] * lObjSense;
    }
	

    nSolver->loadProblem(matrix, localModel_->getOrigColLb(), localModel_->getOrigColUb(),
			 objCoeffs, localModel_->getOrigRowLb(), localModel_->getOrigRowUb());

    for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
    }

    nSolver->setObjSense(-1); //1 min; -1 max

    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    if (feasCheckSolver == "Cbc"){
	    dynamic_cast<OsiCbcSolverInterface *>
		(nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	    sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		(nSolver)->getSymphonyEnvironment();
	    //Always uncomment for debugging!!
	    sym_set_int_param(env, "do_primal_heuristic", FALSE);
	    sym_set_int_param(env, "verbosity", -2);
	    sym_set_int_param(env, "prep_level", -1);
	    sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
	    sym_set_int_param(env, "tighten_root_bounds", FALSE);
	    sym_set_int_param(env, "max_sp_size", 100);
	    sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	    if (whichCutsLL == 0){
		sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	    }else{
		sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	    }
	    if (whichCutsLL == 1){
		sym_set_int_param(env, "generate_cgl_knapsack_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_probing_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_clique_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_twomir_cuts",
				  DO_NOT_GENERATE);
		sym_set_int_param(env, "generate_cgl_flowcover_cuts",
				  DO_NOT_GENERATE);
	    }
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver->setHintParam(OsiDoReducePrint);
	nSolver->messageHandler()->setLogLevel(0);
	    CPXENVptr cpxEnv =
		dynamic_cast<OsiCpxSolverInterface*>(nSolver)->getEnvironmentPtr();
	    assert(cpxEnv);
	    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
    }

    nSolver->branchAndBound();
    
    if (nSolver->isProvenOptimal()){
	bigM = nSolver->getObjValue();
    }
    else{
	assert(0);
    }

    delete nSolver;
    delete [] objCoeffs;
    delete [] integerVars;

    return bigM;
}

//#############################################################################
int
MibSCutGenerator::weakIncObjCutCurrent(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIPINT Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();

  int numCuts(0);
  MibSBilevel * bS = localModel_->bS_;

  //double maxLowerObj(node->getLowerUB());
  double minLowerObj(solver->getInfinity());
  //double maxLowerObj = localModel_->lowerObjectiveBound(true);
  double etol(localModel_->etol_);
  int lN(localModel_->getLowerDim());
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
  double lowerObjSense = localModel_->getLowerObjSense(); 

  /*double * tmpsol = new double[lN + uN];
  CoinZeroN(tmpsol, lN + uN);
  OsiSolverInterface * lSolver = bS->setUpModel(solver, 0, false, true, tmpsol);

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiCbcSolverInterface *> 
     (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (lSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (lSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (lSolver)->setSymParam("max_active_nodes", 1);
#endif

  if(0)
    lSolver->writeLp("boundfinder");

  lSolver->branchAndBound();
  minLowerObj = lSolver->getObjValue();

  //if it's exactly zero, just add a little bit to comply with blis
  if(minLowerObj == 0.0)
  minLowerObj = minLowerObj + etol;*/

  if(0)
    std::cout << "minLowerObj " << minLowerObj << std::endl;

  double bigM(10000);
  int i(0), index(0);
  double cutub(0.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  for(i = 0; i < uN; i++){
    index = upperColInd[i];
    if(sol[index] < etol){
      indexList.push_back(index);
      //valsList.push_back(minLowerObj);
      valsList.push_back(-bigM);
    }
  }

  //double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 
  
  for(i = 0; i < lN; i++){
    index = lowerColInd[i];
    //cutub += lObjCoeffs[i] * sol[index]; 
    if(fabs(lObjCoeffs[i]) > etol){
	//cutub += lObjCoeffs[i] * optLowerSol[i]; //should this be position? 
      indexList.push_back(index);
      valsList.push_back(lowerObjSense * lObjCoeffs[i]);
    }
  }

  cutub = bS->objValVec_[0];


  assert(indexList.size() == valsList.size());

  //numCuts += addCut(conPool, cutlb, 
  //	    lowerObjSense * cutub, indexList, valsList, false);
  numCuts += addCut(conPool, cutlb, 
		    lowerObjSense * cutub, indexList, valsList, true);

  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::weakIncObjCutMaximal(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIPINT Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();

  int numCuts(0);
  MibSBilevel * bS = localModel_->bS_;
  //MibSTreeNode * node = 
  //dynamic_cast<MibSTreeNode *>(localModel_->activeNode_); 
  //double maxLowerObj(node->getLowerUB());
  double maxLowerObj(solver->getInfinity());
  //double maxLowerObj = localModel_->lowerObjectiveBound(true);
  double etol(localModel_->etol_);
  int lN(localModel_->getLowerDim());
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  //const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
  double lowerObjSense = localModel_->getLowerObjSense(); 

  //if we keep this as our big M, should move (same every time)
  double * tmpsol = new double[lN + uN];
  CoinZeroN(tmpsol, lN + uN);
  OsiSolverInterface * lSolver = bS->setUpModel(solver, tmpsol);
  delete [] tmpsol;

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiCbcSolverInterface *> 
     (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (lSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (lSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (lSolver)->setSymParam("max_active_nodes", 1);
#endif
  
  //lSolver->setObjSense(lSolver->getObjSense());
  lSolver->branchAndBound();
  maxLowerObj = lSolver->getObjValue();

  int i(0), index(0);
  double cutub(0.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  if(0)
    std::cout << "maxLowerObj " << maxLowerObj << std::endl;

  const double * sol = findMaximalUpperSol(solver);

  //const double * maximalupper = findMaximalUpperSol(solver);

  if(!sol){

    return 0;
  }
  else{
  
      OsiSolverInterface * lSolver = bS->setUpModel(solver, sol);  

#ifndef COIN_HAS_SYMPHONY
    dynamic_cast<OsiCbcSolverInterface *> 
       (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
    dynamic_cast<OsiSymSolverInterface *> 
       (lSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
       (lSolver)->setSymParam("verbosity", -2);
    
    dynamic_cast<OsiSymSolverInterface *> 
       (lSolver)->setSymParam("max_active_nodes", 1);
#endif

    lSolver->branchAndBound();

    const double * optLowerSol = lSolver->getColSolution();
    
    for(i = 0; i < uN; i++){
      index = upperColInd[i];
      if(sol[index] < etol){
	indexList.push_back(index);
	valsList.push_back(maxLowerObj);
      }
    }

    for(i = 0; i < lN; i++){
      index = lowerColInd[i];
      //cutub += lObjCoeffs[i] * sol[index]; 
      if(fabs(lObjCoeffs[i]) > etol){
	cutub += lObjCoeffs[i] * optLowerSol[i]; 
	indexList.push_back(index);
	valsList.push_back(lowerObjSense * lObjCoeffs[i]);
      }
    }
    
    assert(indexList.size() == valsList.size());
    
    //numCuts += addCut(conPool, cutlb, 
    //    lowerObjSense * cutub, indexList, valsList, false);
    numCuts += addCut(conPool, cutlb, 
		      lowerObjSense * cutub, indexList, valsList, true);

    maximalCutCount_++;

    delete [] sol;
    return numCuts;
  }

}

//#############################################################################
int
MibSCutGenerator::binaryCuts(BcpsConstraintPool &conPool)
{

  //FIXME: NEED TO CHECK FOR ROW TYPES FOR CGLP.  
  //COMING FROM KNAP SOLVER THEY ARE RANGED.

  int useNoGoodCut 
    = localModel_->MibSPar_->entry(MibSParams::useNoGoodCut);

   bool useIncObjCut 
    = localModel_->MibSPar_->entry(MibSParams::useIncObjCut);

  if(useNoGoodCut && !useIncObjCut){
    return noGoodCut(conPool) ? true : false;
  }
  else if(!useNoGoodCut && useIncObjCut){
    return incObjCut(conPool) ? true : false;
  }
  else if(useNoGoodCut && useIncObjCut){
    return (noGoodCut(conPool) && 
	    incObjCut(conPool)) ? true : false;
  }
  else{
    //std::cout << "No BINARY Cuts generated" << std::endl;
    return 0;
  }
}

//#############################################################################
int
MibSCutGenerator::incObjCut(BcpsConstraintPool &conPool)
{

  int numCuts(0);

  //numCuts += incObjCutCurrent(conPool);

  //numCuts += incObjCutMaximal(conPool);

  numCuts += weakIncObjCutCurrent(conPool);
  
  //if(!maximalCutCount_)
  //numCuts += weakIncObjCutMaximal(conPool);

  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::incObjCutMaximal(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  MibSBilevel * bS = localModel_->bS_;
  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  int i(0);
  int numCuts(0);
  double etol(localModel_->etol_);
  int uCols = localModel_->getUpperDim();
  int lCols = localModel_->getLowerDim();
 
  //ACTUALLY, WE WANT THE FULL SOLUTION FROM FIND MAXIMAL
  //THIS GIVES US A FULL CURRENT SOLUTION
  //THEN, WANT TO TAKE UPPER PORTION AND FIND OPTIMAL LL SOL
  //SHOULD BE USING setUpModel(solver, sol) (from MibSBilevel?)

  //figure out how these need to be ordered, cglp is expecting ordered sols

  const double * maximalupper = findMaximalUpperSol(solver);

  double * upperSol = new double[uCols];
  double * lowerSol = new double[lCols];
  double * optLowerSol = new double[lCols];
  CoinZeroN(upperSol, uCols);
  CoinZeroN(lowerSol, lCols);
  CoinZeroN(optLowerSol, lCols);

  /*  
  std::cout << "The maximal solution is: " << std::endl;

  for(i = 0; i < uCols; i++){
    upperSol[i] = maximalupper[uColInd[i]];
    std::cout << uColInd[i] << ": " << upperSol[i] << endl;
  }

  for(i = 0; i < uCols; i++){
    lowerSol[i] = maximalupper[lColInd[i]];
    std::cout << lColInd[i] << ": " << lowerSol[i] << endl;
  }
  */
    //double * lSol = new double[localModel_->getLowerDim()];
    //CoinZeroN(lSol, localModel_->getLowerDim());

  if(!maximalupper){

    return 0;
  }
  else{
  
      OsiSolverInterface * lSolver = bS->setUpModel(solver, maximalupper);  

#ifndef COIN_HAS_SYMPHONY
    dynamic_cast<OsiCbcSolverInterface *> 
       (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
    dynamic_cast<OsiSymSolverInterface *> 
       (lSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
       (lSolver)->setSymParam("verbosity", -2);
    
    dynamic_cast<OsiSymSolverInterface *> 
       (lSolver)->setSymParam("max_active_nodes", 1);
#endif

    lSolver->branchAndBound();

    int lCols = localModel_->getLowerDim();
    
    //CoinCopyN(lSolver->getColSolution(), localModel_->getLowerDim(), lSol);
    for(i = 0; i < lCols; i++)
      optLowerSol[i] = lSolver->getColSolution()[i];

    double cutlb(- solver->getInfinity());
    double cutub(solver->getInfinity());
    
    std::vector<int> indexList;
    std::vector<double> valsList;
    
    //replace these with results from above
    //NOT CORRECT NOW
    //MibSBilevel * bS = localModel_->bS_; //not needed once fixed below
    //double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
    //double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
    //double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

    /* 
       returns a double with the values [alpha | beta | gamma]
       indexed 0..n1-1, n1..n1+n2-1, n1+n2
    */
    
    double * cutVals = findDeepestLandPCut_IncObj(upperSol, 
						  lowerSol, 
						  optLowerSol);  
    double val(0.0);
    
    /* find the nonzero entries */
    if(0)
      std::cout << "maximal solution." << std::endl;

    for(i = 0; i < numCols; i++){
      val = cutVals[i];
      if((val > etol) || (val < - etol)){
	// entry is nonzero
	indexList.push_back(i);
	if(0)
	  std::cout << "i: " << i << " -- "  << val << std::endl;
	valsList.push_back(val);
      }
    }
    
    cutub = cutVals[numCols]; //gamma

    if(0)
      std::cout << "cutub: " << cutub << std::endl;
    
    assert(indexList.size() == valsList.size());
    
    numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);
    
    delete [] cutVals;
    return numCuts;
  }

}

//#############################################################################
const double *
MibSCutGenerator::findMaximalUpperSol(OsiSolverInterface * si)
{

  OsiSolverInterface * maxSolver = new OsiCbcSolverInterface(si);

  int numCols = maxSolver->getNumCols();
  int numULCols = localModel_->getUpperDim();
  int * upperColInd = localModel_->getUpperColInd();

  //  double * maximalsol = new double[numCols];
  //CoinZeroN(maximalsol, numCols);

  int i(0);
  double * objectives = new double[numCols];
  CoinZeroN(objectives, numCols);

  //maximizing the sum of the upper-level variables 
  for(i = 0; i < numULCols; i++)
    objectives[upperColInd[i]] = 1.0;
 
  maxSolver->setObjective(objectives);
  maxSolver->setObjSense(-1); //maximization

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiCbcSolverInterface *> 
     (maxSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (maxSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (maxSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (maxSolver)->setSymParam("max_active_nodes", 1);
#endif
    
  maxSolver->branchAndBound();
  //maxSolver->initialSolve();

  if(0)
    maxSolver->writeLp("maxsolver");

  if(maxSolver->isProvenOptimal()){

    double * maxsol = new double[maxSolver->getNumCols()];
    CoinCopyN(maxSolver->getColSolution(), maxSolver->getNumCols(), maxsol);
    
    if(1){
      std::cout << "Maximal UL Solution found; cardinality of " 
		<< maxSolver->getObjValue() << std::endl;
    }

    delete [] objectives;
    //delete maxSolver;

    return (maxsol);

  }
  else{

    std::cout << "Unable to find maximal solution." << std::endl;
    delete [] objectives;
    //delete maxSolver;
    return NULL;

  }


}

//#############################################################################
int
MibSCutGenerator::incObjCutCurrent(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  OsiSolverInterface * solver = localModel_->solver();
  MibSBilevel * bS = localModel_->bS_;

  const int numCols = solver->getNumCols();
  int i(0);
  int numCuts(0);
  double etol(localModel_->etol_);
 
  double cutlb(- solver->getInfinity());
  double cutub(solver->getInfinity());

  std::vector<int> indexList;
  std::vector<double> valsList;

  double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  double * optLowerSol = bS->optLowerSolutionOrd_; // optimal LL solution 

  /* 
     returns a double with the values [alpha | beta | gamma]
     indexed 0..n1-1, n1..n1+n2-1, n1+n2
  */

  double * cutVals = findDeepestLandPCut_IncObj(upperSol, 
						lowerSol, 
						optLowerSol);  
  double val(0.0);

  /* find the nonzero entries */

  for(i = 0; i < numCols; i++){
    val = cutVals[i];
    if((val > etol) || (val < - etol)){
      // entry is nonzero
      indexList.push_back(i);
      valsList.push_back(val);
    }
  }

  cutub = cutVals[numCols]; //gamma

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);

  delete [] cutVals;
  return numCuts;

}

//###########################################################################
int
MibSCutGenerator::noGoodCut(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating No-Good Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();

  int numCuts(0);
  double etol(localModel_->etol_);
  const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
 
  int i(0), index(0);
  double cutub(- 1.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  for(i = 0; i < uN; i++){
    index = upperColInd[i];
    indexList.push_back(index);
    if(sol[index] > etol){
      valsList.push_back(1.0);
      cutub += 1.0;
    }
    else{
      valsList.push_back(- 1.0);
    }
  }

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, false);
  
  return numCuts;

}

//###########################################################################
int
MibSCutGenerator::bendersInterdictionOneCut(BcpsConstraintPool &conPool, double *lSolution)
{

  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
      (localModel_->MibSPar_->entry(MibSParams::branchStrategy));

  bool allowRemoveCut(localModel_->MibSPar_->entry(MibSParams::allowRemoveCut));

  //when the branching strategy is fractional and the optimal
  //solution of relaxation is integer, we are forced to generate cut.
  /*if(branchPar != MibSBranchingStrategyFractional){
      if (localModel_->boundingPass_ > 1){
	  return 0;
      }
  }*/

  int i;
  int indexU(0), indexL(0);
  int numCuts(0);
  double etol(localModel_->etol_);
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  double cutub(localModel_->solver()->getInfinity());
  double cutlb(0.0);
  std::vector<int> indexList;
  std::vector<double> valsList;

  for(i = 0; i < uN; i++){
      indexU = upperColInd[i];
      indexL = lowerColInd[i];
      indexList.push_back(indexL);
      valsList.push_back(-lObjCoeffs[i]);
      cutlb += -1 * lObjCoeffs[i] * lSolution[i];
      if(lSolution[i] > etol){
	  indexList.push_back(indexU);
	  valsList.push_back(-lObjCoeffs[i]*lSolution[i]);
      }
  }
  assert(indexList.size() == valsList.size());
  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, allowRemoveCut);

  indexList.clear();
  valsList.clear();

  return numCuts;
}

//###########################################################################
int
MibSCutGenerator::bendersInterdictionMultipleCuts(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/
  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
      (localModel_->MibSPar_->entry(MibSParams::branchStrategy));

  bool allowRemoveCut(localModel_->MibSPar_->entry(MibSParams::allowRemoveCut));

  //when the branching strategy is fractional and the optimal
  //solution of relaxation is integer, we are forced to generate cut.
  if(branchPar != MibSBranchingStrategyFractional){
      if (localModel_->boundingPass_ > 1){
	  return 0;
      }
  }
   
  OsiSolverInterface * solver = localModel_->solver();

  const double *lpSol = solver->getColSolution();
  
  int numCuts(0);
  double etol(localModel_->etol_);
  int uN(localModel_->upperDim_);
  int lN(localModel_->getLowerDim());
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
 
  int i(0), j(0);
  double cutub(solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  int indexU(0), indexL(0);

#if 0
  OsiSolverInterface * lSolver = localModel_->getMibSBilevel()->solver_;
  sym_environment *env = 
     dynamic_cast<OsiSymSolverInterface *>(lSolver)->getSymphonyEnvironment();  
  int sp_size(0);
  sym_get_sp_size(env, &sp_size);
  
  double ll_objval;
  double * ll_sol = new double[lN];

  for (i = 0; i < sp_size; i++){
     sym_get_sp_solution(env, i, ll_sol, &ll_objval);
     for(j = 0; j < uN; j++){
	indexU = upperColInd[j];
	indexL = lowerColInd[j];
	indexList.push_back(indexL);
	valsList.push_back(-lObjCoeffs[j]);
	if (ll_sol[j] > etol){
	   indexList.push_back(indexU);
	   valsList.push_back(-lObjCoeffs[j]);
	}
     }
     assert(indexList.size() == valsList.size());
     numCuts += addCut(conPool, -ll_objval, cutub, indexList, valsList, 
		       false);
     indexList.clear();
     valsList.clear();
  }

  delete[] ll_sol;
#endif

  std::vector<std::pair<AlpsKnowledge*, double> > solutionPool;
  localModel_->getKnowledgeBroker()->
     getAllKnowledges(AlpsKnowledgeTypeSolution, solutionPool);
  
  const double * sol; 
  double objval, lhs(0);
  BlisSolution* blisSol;
  std::vector<std::pair<AlpsKnowledge*, double> >::const_iterator si;
  for (si = solutionPool.begin(); si != solutionPool.end(); ++si){
     lhs = 0;
     blisSol = dynamic_cast<BlisSolution*>(si->first);
     sol = blisSol->getValues();
     objval = blisSol->getQuality();
     for(j = 0; j < uN; j++){
	indexU = upperColInd[j];
	indexL = lowerColInd[j];
	indexList.push_back(indexL);
	valsList.push_back(-lObjCoeffs[j]);
	lhs -= lObjCoeffs[j]*lpSol[indexL];
	if (sol[indexL] > etol){
	   indexList.push_back(indexU);
	   valsList.push_back(-lObjCoeffs[j]*sol[indexL]);
	   //lhs -= lObjCoeffs[j]*lpSol[indexU];
	   lhs -= lObjCoeffs[j]*sol[indexL]*lpSol[indexU];
	}
     }
     assert(indexList.size() == valsList.size());
     if (lhs - objval < -etol){
	 numCuts += addCut(conPool, objval, cutub, indexList, valsList,
			   allowRemoveCut);
     }
     indexList.clear();
     valsList.clear();
  }
  
  return numCuts;

}

//###########################################################################
// WARNING!! These cuts are not valid at the moment. Experiment gone wrong!
//###########################################################################
int
MibSCutGenerator::bendersZeroSumCuts(BcpsConstraintPool &conPool)
{

  if (localModel_->boundingPass_ > 1){
     return 0;
  }
   
  OsiSolverInterface * solver = localModel_->solver();

  const double *lpSol = solver->getColSolution();
  
  int numCuts(0);
  double etol(localModel_->etol_);
  int uN(localModel_->upperDim_);
  int lN(localModel_->getLowerDim());
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  const double * colUpper = solver->getColUpper();
  const double * colLower = solver->getColLower();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
 
  int i(0), j(0);
  double cutub(solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  int indexU(0), indexL(0);

  std::vector<std::pair<AlpsKnowledge*, double> > solutionPool;
  localModel_->getKnowledgeBroker()->
     getAllKnowledges(AlpsKnowledgeTypeSolution, solutionPool);
  
  const double * sol; 
  double objval, lhs(0);
  BlisSolution* blisSol;
  std::vector<std::pair<AlpsKnowledge*, double> >::const_iterator si;
  for (si = solutionPool.begin(); si != solutionPool.end(); ++si){
     lhs = 0;
     blisSol = dynamic_cast<BlisSolution*>(si->first);
     sol = blisSol->getValues();
     objval = blisSol->getQuality();
     //Check to see if solution is feasible for current node 
     for (j = 0; j < uN+lN; j++){
	indexU = upperColInd[j];
	if (sol[j] > colUpper[j] || sol[j] < colLower[j]){
	   break;
	}
     }
     //Don't add the cut if we found infeasibility
     if (j < lN+uN){
	continue;
     }
     
     for(j = 0; j < lN; j++){
	indexL = lowerColInd[j];
	indexList.push_back(indexL);
	valsList.push_back(-lObjCoeffs[j]);
	lhs -= lObjCoeffs[j]*lpSol[indexL];
     }
     assert(indexList.size() == valsList.size());
     if (lhs < objval){
	numCuts += addCut(conPool, objval, cutub, indexList, valsList, 
			  false);
     }
     indexList.clear();
     valsList.clear();
  }

  return numCuts;

}

//###########################################################################
bool 
MibSCutGenerator::generateConstraints(BcpsConstraintPool &conPool)
{

  //FIXME: MAKE THIS MORE SIMPLE
  int numCuts(0);

  //OsiSolverInterface * solver = localModel_->solver();

  int cutStrategy =
    localModel_->MibSPar_->entry(MibSParams::cutStrategy);

  MibSTreeNode * node =
      dynamic_cast<MibSTreeNode *>(localModel_->activeNode_);

  if(cutStrategy > 0){
     //use bilevel feasibility cut
    MibSBilevel *bS = localModel_->bS_;
    //int problemType = 
    //  localModel_->MibSPar_->entry(MibSParams::bilevelProblemType);
    int cutTypes = 
      localModel_->MibSPar_->entry(MibSParams::bilevelCutTypes);

    bool useBoundCut = 
       localModel_->MibSPar_->entry(MibSParams::useBoundCut);
    
    int useBendersCut = 
       localModel_->MibSPar_->entry(MibSParams::useBendersCut);

    //int useIntersectionCut =
    //localModel_->MibSPar_->entry(MibSParams::useIntersectionCut);

    MibSIntersectionCutType cutType(MibSIntersectionCutTypeNotSet);
    
    int useIntersectionCutTypeIC =
	localModel_->MibSPar_->entry(MibSParams::useTypeIC);

    int useIntersectionCutTypeWatermelon =
	localModel_->MibSPar_->entry(MibSParams::useTypeWatermelon);

    int useIntersectionCutTypeHypercubeIC =
	localModel_->MibSPar_->entry(MibSParams::useTypeHypercubeIC);

    int useIntersectionCutTypeTenderIC =
	localModel_->MibSPar_->entry(MibSParams::useTypeTenderIC);

    int useIntersectionCutTypeHybridIC =
	localModel_->MibSPar_->entry(MibSParams::useTypeHybridIC);

    int useGeneralNoGoodCut = 
	localModel_->MibSPar_->entry(MibSParams::useGeneralNoGoodCut);

    bool useIncObjCut
	= localModel_->MibSPar_->entry(MibSParams::useIncObjCut);

    CoinPackedVector *sol = localModel_->getSolution();

    if(cutTypes == 0){
      //general type of problem, no specialized cuts
      delete sol;

      if ((useBoundCut) && (localModel_->boundingPass_ <= 1)){

	  int boundCutFreq(localModel_->MibSPar_->entry(MibSParams::boundCutFreq));
	  
	  int boundCutDepthLb(localModel_->MibSPar_->entry
			      (MibSParams::boundCutDepthLb));

	  int boundCutDepthUb(localModel_->MibSPar_->entry
			      (MibSParams::boundCutDepthUb));

	  int depth = static_cast<MibSTreeNode *>(localModel_->activeNode_)->getDepth();

	  bool boundCutOptimal
	      = localModel_->MibSPar_->entry(MibSParams::boundCutOptimal);

	  int boundCutOptimalType =
	      localModel_->MibSPar_->entry(MibSParams::boundCutOptimalType);

	  useBoundCut = false;

	  if((boundCutDepthLb >= 0) && (boundCutDepthUb >= 0)){
	      if((depth >= boundCutDepthLb) && (depth <= boundCutDepthUb)){
		  useBoundCut = true;
	      }
	  }
	  else if(boundCutDepthLb >= 0){
	      if(depth >= boundCutDepthLb){
		  useBoundCut = true;
	      }
	  }
	  else if(boundCutDepthUb >= 0){
	      if(depth <= boundCutDepthUb){
		  useBoundCut = true;
	      }
	  }
	  else{
	      useBoundCut = true;
	  }

	  if(node->getIndex() == 0){
	      if((boundCutOptimal == true) &&
		 (boundCutOptimalType == MibSBoundCutOptimalTypeParametric)){
		  useBoundCut = true;
	      }
	  }

	  if(useBoundCut == true){
	      numCalledBoundCut_ ++;
	  }

	  if(((numCalledBoundCut_%boundCutFreq) == 0) && (useBoundCut == true)){
	      double tmpArg1 = 0;
	      bool tmpArg2 = false;
	      boundCuts(conPool, NULL, tmpArg1, tmpArg2);
	  }
      }
      
      if (bS->isIntegral_){
	  //if (useIntersectionCut == PARAM_ON){
	  //  intersectionCuts(conPool, bS->optLowerSolutionOrd_);
	  //}
	  if (useIntersectionCutTypeIC == PARAM_ON){
	      cutType = MibSIntersectionCutTypeIC;
	      intersectionCuts(conPool, bS->optLowerSolutionOrd_, cutType);
	  }

	  if (useIntersectionCutTypeWatermelon == PARAM_ON){
	      cutType = MibSIntersectionCutTypeWatermelon;
	      intersectionCuts(conPool, bS->optLowerSolutionOrd_, cutType);
	  }

	  if (useIntersectionCutTypeHypercubeIC == PARAM_ON){
	      cutType = MibSIntersectionCutTypeHypercubeIC;
	      intersectionCuts(conPool, bS->optLowerSolutionOrd_, cutType);
	  }

	  if (useIntersectionCutTypeTenderIC == PARAM_ON){
	      cutType = MibSIntersectionCutTypeTenderIC;
	      intersectionCuts(conPool, bS->optLowerSolutionOrd_, cutType);
	  }

	  if (useIntersectionCutTypeHybridIC == PARAM_ON){
	      cutType = MibSIntersectionCutTypeHybridIC;
	      intersectionCuts(conPool, bS->optLowerSolutionOrd_, cutType);
	  }

	  if (useGeneralNoGoodCut == PARAM_ON){
	      generalNoGoodCut(conPool);
	  }
	  
	  if (useBendersCut == PARAM_ON){
	      int bendersCutType =
		  localModel_->MibSPar_->entry(MibSParams::bendersCutType);
	      if(bendersCutType == MibSBendersCutTypeJustOneCut){
		  numCuts += bendersInterdictionOneCut(conPool,
						       bS->optLowerSolutionOrd_);
	      }
	      else{
		  numCuts += bendersInterdictionMultipleCuts(conPool);
	      }		  
	  }
	  
	  if (useIncObjCut == true){
	      numCuts += generalWeakIncObjCutCurrent(conPool);
	  }
	  
	 numCuts += feasibilityCuts(conPool) ? true : false;
	 
         /*if (useBoundCut){
	     double tmpArg1 = 0;
	     bool tmpArg2 = false;
	     boundCuts(conPool, NULL, tmpArg1, tmpArg2);
	     }*/
      }
      return (numCuts ? true : false);
    }
    else if(bS->isUpperIntegral_ && cutTypes == 1){
      //interdiction problem
      delete sol;
      int status = false;
      if (bS->isIntegral_){
	 numCuts += feasibilityCuts(conPool);
      }
      numCuts += interdictionCuts(conPool);
      return (numCuts ? true : false);
    }
    else if(bS->isUpperIntegral_ && cutTypes == 2){
      //problem with binary UL variables and integer LL variables
      delete sol;
      int status = false;
      if (bS->isIntegral_){
	 status = feasibilityCuts(conPool) ? true : false;
      }
      numCuts += binaryCuts(conPool);
      return (numCuts ? true : false);
    }
    else if(bS->isUpperIntegral_ && cutTypes == 3){
      //problem with binary UL variables and general LL variables
      delete sol;
      return binaryCuts(conPool) ? true : false;
    }
  
    /*
      if(bS->isIntegral_){
      delete sol;
      return feasibilityCuts(conPool) ? true : false;
      }
    */
    
    /*
      May want to add our own cuts later to separate fractional solutions
    */
    
    delete sol;
  }

  return numCuts ? true : false;

}

//#############################################################################
int 
MibSCutGenerator::addCut(BcpsConstraintPool &conPool,
			 double cutlb, double cutub,
			 std::vector<int> & indexList, 
			 std::vector<double> &elementList,
			 bool removable)
{

  int size(indexList.size());
  int * indices = new int[size];
  double * elements = new double[size];
  int capacity(localModel_->solver()->getNumCols());

  int i(0);

  for(i = 0; i < size; i++){
    indices[i] = indexList.at(i);
    elements[i] = elementList.at(i);
  }

  OsiRowCut *cut = new OsiRowCut(cutlb, cutub, capacity, 
				 size, indices, elements);

  BlisConstraint *blisCon = BlisOsiCutToConstraint(cut);

  if(!removable)
    blisCon->setStatus(BCPS_NONREMOVALBE);
  
  conPool.addConstraint(blisCon);

  delete cut;
  return 1;

}

//#############################################################################
int *
MibSCutGenerator::getBindingCons()
{

   //std::string method = model_->bindingMethod_; 
   int method = localModel_->MibSPar_->entry(MibSParams::whichActiveConMethod);
   int * binding = NULL;

  if(method == 0){
     binding = getBindingConsSimple();
  }
  else if(method == 1){
     binding = getBindingConsBasis();
  }
  else{
     std::cout << "No method for binding constraint set." << std::endl;
     return NULL;
  }
  
  return binding;
  
}

//#############################################################################
int *
MibSCutGenerator::getBindingConsSimple()
{

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const char * rowsense = solver->getRowSense();
  const double * rowrange = solver->getRowRange();
  const double * sol = solver->getColSolution();
  const double * rhs = solver->getRightHandSide();
  int i(0), index(0);
  double etol(localModel_->etol_);

  //FIXME: CHECK THIS
  const double * collb;
  const double * colub;

  if(localModel_->isRoot_){
    collb = localModel_->getSolver()->getColLower();
    colub = localModel_->getSolver()->getColUpper();
  }
  else{
    collb = solver->getColLower();
    colub = solver->getColUpper();
  }

  int *colStatus = new int[numCols];
  int *rowStatus = new int[numRows];
  solver->getBasisStatus(colStatus, rowStatus);

  /** Find all binding constraints at current solution **/

  int bindingCons(0);
  double slack(0.0);
  double upper(0.0);
  
  int * binding = new int[numRows + 2 * numCols];
  double * slackVal = new double[numRows + 2 * numCols];
  const double * rowActivity = solver->getRowActivity();
  int MAXCONS(numCols);
  //int MAXCONS(100);  

  CoinZeroN(binding, numRows + 2 * numCols);
  CoinFillN(slackVal, numRows + 2 * numCols, 0.0);
  
  for(i = 0; i < numRows; i++){
     
    if(1){//if(rowStatus[i] != 3) && (rowStatus[i] != 2)){
      if((rowsense[i] =='R') && (rowrange[i] < 0))
	slackVal[i] = - rhs[i] + rowActivity[i];
      else
	slackVal[i] = rhs[i] - rowActivity[i];
      
      slack = slackVal[i];
      
      switch(rowsense[i]){
      case 'L':
	if(slack > etol){
	  if(0) 
	    std::cout << "constraint " << i << " is not tight." << std::endl;
	  binding[i] = 0;
	}
	else if((slack > - etol) && (bindingCons < MAXCONS)){
	  if(0){
	    std::cout << "constraint " << i << " is tight." << std::endl;
	    std::cout << "status is " << rowStatus[i] << std::endl;
	  }
	  binding[i] = 1;
	  upper += rhs[i];
	  //upper += rhs[i] * solver->getRowPrice()[i];
	  bindingCons++;
	}
	break;
      case 'G':
	if(slack < - etol){
	  if(0) 
	    std::cout << "constraint " << i << " is not tight." << std::endl;
	  binding[i] = 0;
	}
	else if((slack < etol) && (bindingCons < MAXCONS)){
	  if(0){
	    std::cout << "constraint " << i << " is tight." << std::endl;
	    std::cout << "status is " << rowStatus[i] << std::endl;
	  }
	  binding[i] = 1;
	  upper -= rhs[i];
	  //upper -= rhs[i] * solver->getRowPrice()[i];
	  bindingCons++;
	}
	break;
      case 'E':
	if((slack < - etol) || (slack > etol)){
	  if(0) 
	    std::cout << "constraint " << i << " is not tight." << std::endl;
	  binding[i] = 0;
	}
	else if(bindingCons < MAXCONS){
	  if(0){
	    std::cout << "constraint " << i << " is tight." << std::endl;
	    std::cout << "status is " << rowStatus[i] << std::endl;
	  }
	  binding[i] = 1;
	  upper += rhs[i];
	  //upper += rhs[i] * solver->getRowPrice()[i];
	  bindingCons++;
	}
	break;
      case 'R':
	if(rowrange[i] < 0){
	  if((slack > rowrange[i] + etol) || (slack < - etol)){
	    //std::cout << "constraint " 
	    //	<< i << " is not tight." << std::endl;
	    binding[i] = 0;
	  }
	  else if(((slack > rowrange[i] - etol) || (slack < etol)) 
		  && (bindingCons < MAXCONS)){
	    if(0){
	      std::cout << "constraint " << i << " is tight." << std::endl;
	      std::cout << "status is " << rowStatus[i] << std::endl;
	    }	   
	    binding[i] = 1;
	    upper += rhs[i];
	    //upper += rhs[i] * solver->getRowPrice()[i];
	    bindingCons++;
	  }
	}
	else{
	  if((slack < rowrange[i] - etol) || (slack > etol)){
	    if(0)
	      std::cout << "constraint " 
			<< i << " is not tight." << std::endl;
	    binding[i] = 0;
	  }
	  else if(((slack < rowrange[i] + etol) || (slack > - etol)) 
		  && (bindingCons < MAXCONS)){
	    if(0){
	      std::cout << "constraint " << i << " is tight." << std::endl;
	      std::cout << "status is " << rowStatus[i] << std::endl;
	    }
	    binding[i] = 1;
	    upper += rhs[i];
	    //upper += rhs[i] * solver->getRowPrice()[i];
	    bindingCons++;
	  }
	}
	break;
      }
    }
  }

  double upSlack(0.0);
  double downSlack(0.0);

  for(i = 0; i < numCols; i++){
    if(0){
      std::cout << "status of column " << i << ": " 
		<< colStatus[i] << std::endl;
    }
    index = numRows + i;
   
    /*   
    if((fabs(colub[i] - collb[i]) < etol) && (colStatus[i] == 3)){
      if(1)
	std::cout << "col " << i << " is fixed with status." << colStatus[i] 
		  << std::endl;
      //binding[index] = 1;
      //upper += colub[i];
      //bindingCons++;
    }
    */
    //else{// 
    if(1){//if(colStatus[i] != 1){ 
      upSlack = colub[i] - sol[i];
      downSlack = sol[i] - collb[i];
      
      if((upSlack > - etol) && (upSlack < etol)){
	if(0){
	  std::cout << "upper bound " << i << " is tight." << std::endl;
	  std::cout << "status is " << colStatus[i] << std::endl;
	}
	binding[index] = 1;
	upper += colub[i];
	bindingCons++;
      }
      else if((downSlack > - etol) && (downSlack < etol)){
	if(0){
	  std::cout << "lower bound " << i << " is tight." << std::endl;
	  std::cout << "status is " << colStatus[i] << std::endl;
	}
	binding[index + numCols] = 1;
	upper -= collb[i];
	bindingCons++;
      }
    }
  }

  
  if(!bindingCons){
     std::cout << "No constraints are binding." << std::endl;
     //return 0;
     abort();
  }
  
  /** Set the upper bound for the current cut **/ 
  setCutUpperBound(upper);
  
  delete [] slackVal;  

  return binding;
  
}

//#############################################################################
int *
MibSCutGenerator::getBindingConsBasis()
{

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const char * rowsense = solver->getRowSense();
  const double * collb = solver->getColLower();
  const double * colub = solver->getColUpper();
  const double * rhs = solver->getRightHandSide();
  int i(0), index(0);

  /** Find binding constraints at current solution using basis information **/

  int bindingCons(0);
  double upper(0.0);

  int * binding = new int[numRows + 2 * numCols];
  CoinZeroN(binding, numRows + 2 * numCols);

  CoinWarmStartBasis::Status rowStatus;
  CoinWarmStartBasis * ws = 
     dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
  
  for(i = 0; i < numRows; i++){
    
    rowStatus = ws->getArtifStatus(i);
    
     switch(rowsense[i]){
      case 'L':
	if(rowStatus == CoinWarmStartBasis::basic){
	  binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper += rhs[i];
	  bindingCons++;
	}
	break;
      case 'G':
	if(rowStatus == CoinWarmStartBasis::basic){
	  binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper -= rhs[i];
	  bindingCons++;
	}
	break;
      case 'E':
	if(rowStatus == CoinWarmStartBasis::basic){
	   binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper += rhs[i];
	  bindingCons++;
	}
	break;
      case 'R':
	if(rowStatus == CoinWarmStartBasis::basic){
	   binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper += rhs[i];
	  bindingCons++;
	}
	break;
     }
  }
  
  /** Find binding bounds at current solution **/

  const double * sol = solver->getColSolution();
  double etol(localModel_->etol_);

  double upSlack(0.0);
  double downSlack(0.0);
  
  for(i = 0; i < numCols; i++){
    index = numRows + i;
    upSlack = colub[i] - sol[i];
    downSlack = sol[i] - collb[i];
    
    if((upSlack > - etol) && (upSlack < etol)){
      //std::cout << "upper bound " << i << " is tight." << std::endl;
      binding[index] = 1;
      upper += colub[i];
      bindingCons++;
     }
    else if((downSlack > - etol) && (downSlack < etol)){
      binding[index + numCols] = 1;
      upper -= collb[i];
      bindingCons++;
    }
  }

  
  if(!bindingCons){
    std::cout << "Invalid basis." << std::endl;
    //return 0;
    abort();
  }
  
  /** Set the upper bound for the current cut **/ 
  setCutUpperBound(upper);
  
  return binding;

}

//#############################################################################
void
MibSCutGenerator::getLowerMatrices(bool getLowerConstCoefMatrix,
				   bool getA2Matrix, bool getG2Matrix)
{
    bool isLowerConstCoefMatrixSet(false), isMatrixA2Set(false), isMatrixG2Set(false);
    int i, j;
    int index(0), numElements(0), pos(0);
    int numCols(localModel_->getNumCols());
    int uCols(localModel_->getUpperDim());
    int lCols(localModel_->getLowerDim());
    int lRows(localModel_->getLowerRowNum());
    int *uColInd(localModel_->getUpperColInd());
    int *lColInd(localModel_->getLowerColInd());
    int *lRowInd(localModel_->getLowerRowInd());
    char *rowSense = localModel_->getOrigRowSense();
    CoinPackedMatrix origMatrix = *localModel_->getOrigConstCoefMatrix();

    CoinShallowPackedVector origRow;
    CoinPackedVector row;
    CoinPackedVector rowA2;
    CoinPackedVector rowG2;
    CoinPackedMatrix *matrixA2G2 = 0;
    CoinPackedMatrix *matrixA2 = 0;
    CoinPackedMatrix *matrixG2 = 0;

    if(localModel_->getLowerConstCoefMatrix()){
	isLowerConstCoefMatrixSet = true;
    }

    if(localModel_->getA2Matrix()){
	isMatrixA2Set = true;
    }

    if(localModel_->getG2Matrix()){
	isMatrixG2Set =true;
    }
    
    origMatrix.reverseOrdering();
    if((getLowerConstCoefMatrix) && (!isLowerConstCoefMatrixSet)){
	matrixA2G2 = new CoinPackedMatrix(false, 0, 0);
	matrixA2G2->setDimensions(0, numCols);
    }
    if((getA2Matrix) && (!isMatrixA2Set)){
	matrixA2 = new CoinPackedMatrix(false, 0, 0);
        matrixA2->setDimensions(0, uCols);
    }
    if((getG2Matrix) && (!isMatrixG2Set)){
	matrixG2 = new CoinPackedMatrix(false, 0, 0);
	matrixG2->setDimensions(0, lCols);
    }
    
    //extracting matrices A2G2, A2 and G2
    for(i = 0; i < lRows; i++){
	index = lRowInd[i];
	origRow = origMatrix.getVector(index);
	if(rowSense[index] == 'G'){
	    row = -1 * origRow;
	}
	else{
	    row = origRow;
	}
	if(matrixA2G2){
	    matrixA2G2->appendRow(row);
	}
	numElements = row.getNumElements();
	const int *indices = row.getIndices();
	const double *elements = row.getElements();
	for(j = 0; j < numElements; j++){
	    pos = localModel_->binarySearch(0, uCols -1, indices[j], uColInd);
	    if(pos >= 0){
		rowA2.insert(pos, elements[j]);
	    }
	    else{
		pos = localModel_->binarySearch(0, lCols -1, indices[j], lColInd);
		rowG2.insert(pos, elements[j]);
	    }
	}
	if(matrixA2){
	    matrixA2->appendRow(rowA2);
	}
	if(matrixG2){
	    matrixG2->appendRow(rowG2);
	}
	rowA2.clear();
	rowG2.clear();
    }

    if((getLowerConstCoefMatrix) && (matrixA2G2)){
	localModel_->setLowerConstCoefMatrix(matrixA2G2);
    }
    if((getA2Matrix) && (matrixA2)){
	localModel_->setA2Matrix(matrixA2);
    }
    if((getG2Matrix) && (matrixG2)){
	localModel_->setG2Matrix(matrixG2);
    }

}
