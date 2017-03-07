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

#include "CoinPackedVector.hpp"
#include "OsiCbcSolverInterface.hpp"

#include "MibSBilevel.hpp"
#include "MibSModel.hpp"
#include "MibSTreeNode.hpp"
#include "MibSSolution.hpp"
#include "MibSHeuristic.hpp"
#include "MibSConfig.hpp"

#ifdef COIN_HAS_SYMPHONY
#include "symphony.h"
#include "SymConfig.h"
#include "OsiSymSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CPLEX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif

//#############################################################################
void 
MibSBilevel::createBilevel(CoinPackedVector* sol, 
			   MibSModel *mibs)
{
  
  /** Splits sol into two parts, upper- and lower-level components **/

  if(!mibs) return;
  
  model_ = mibs;

  heuristic_ = new MibSHeuristic(mibs);
  
  int i(0),j(0);
  int N(model_->numVars_);
  int lN(model_->lowerDim_); // lower-level dimension
  int uN(model_->upperDim_); // upper-level dimension
  double etol(model_->BlisPar()->entry(BlisParams::integerTol));

  

  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
      (model_->MibSPar_->entry(MibSParams::branchStrategy));

  int solveLowerWhenXYVarsInt(model_->MibSPar_->entry
			   (MibSParams::solveLowerWhenXYVarsInt));
  int solveLowerWhenXVarsInt(model_->MibSPar_->entry
			  (MibSParams::solveLowerWhenXVarsInt));
  int solveLowerWhenLinkVarsInt(model_->MibSPar_->entry
			  (MibSParams::solveLowerWhenLinkVarsInt));
  int solveLowerWhenLinkVarsFixed(model_->MibSPar_->entry
			    (MibSParams::solveLowerWhenLinkVarsFixed));
  int cutStrategy(model_->MibSPar_->entry
		  (MibSParams::cutStrategy));

  int saveSeenLinkingSols(model_->MibSPar_->entry
			   (MibSParams::saveSeenLinkingSols));
  
  assert(N == model_->solver()->getNumCols());

  model_->setNumRows(model_->solver()->getNumRows());
  model_->setUpperRowNum(model_->solver()->getNumRows() - model_->getLowerRowNum());
  model_->setUpperRowData();
  
  int *indices = sol->getIndices();
  double *values = sol->getElements();
  int numElements(sol->getNumElements()); // number of nonzero elements
  int * fixedInd = model_->fixedInd_; 
 
  if(!upperSolutionOrd_)
      upperSolutionOrd_ = new double[uN];
  if(!lowerSolutionOrd_)
      lowerSolutionOrd_ = new double[lN];
  if(!optUpperSolutionOrd_)
      optUpperSolutionOrd_ = new double[uN];
  if(!optLowerSolutionOrd_)
      optLowerSolutionOrd_ = new double[lN];
  
  CoinZeroN(upperSolutionOrd_, uN);
  CoinZeroN(lowerSolutionOrd_, lN);
  CoinZeroN(optUpperSolutionOrd_, uN);
  CoinZeroN(optLowerSolutionOrd_, lN);
  
  isIntegral_ = true;
  isUpperIntegral_ = true;
  isLinkVarsIntegral_ = true;
  LPSolStatus_ = MibSLPSolStatusUnknown;
  isLinkVarsFixed_ = true;
  shouldPrune_ = false;
  isLowerSolved_ = false;
  isUBSolved_ = false;
  isContainedInSetE_ = false;
  useBilevelBranching_ = true;
  isProvenOptimal_ = false;
  tagInSeenLinkingPool_ = MibSSetETagIsNotSet;
  haveHeurSolCand_ = false;

  model_->countIteration_ ++;
  /*std::cout << "countIteration = " << model_->countIteration_ << std::endl;
  if(model_->countIteration_ == 821){
      std::cout << "Stop here!" << std::endl;
      }*/
  
  int * lowerColInd = mibs->getLowerColInd();
  int * upperColInd = mibs->getUpperColInd();
  int index(0), uCount(0), lCount(0);

  const double *lower = model_->solver()->getColLower();
  const double *upper = model_->solver()->getColUpper();
  double value;
  for(i = 0; i < numElements; i ++){
    index = indices[i];
    value = CoinMax(values[i], lower[index]);
    value = CoinMin(value, upper[index]);
    if(binarySearch(0, uN - 1, index, upperColInd) >= 0){
       if(fabs(floor(value + 0.5) - value) > etol){
#if 1
	   if(fixedInd[index] == 1){
	       isLinkVarsIntegral_ = false;
	   }
	   if(mibs->solver()->isInteger(index)){
	       isUpperIntegral_ = false;
	       isIntegral_ = false;
	       LPSolStatus_ = MibSLPSolStatusInfeasible;
	   }
#endif
       }  
    }
    else{
       if(fabs(floor(value + 0.5) - value) > etol){
#if 1
	  //This check is failing when Blis has already declared the solution integral
	  //It's not really needed
	  if(mibs->solver()->isInteger(index)){
	     isIntegral_ = false;
	     LPSolStatus_ = MibSLPSolStatusInfeasible;
	  }
#endif
       }    
    }
  }

  for(i = 0; i < N; i ++){
      if(binarySearch(0, uN - 1, i, upperColInd) >= 0){
	  if((fixedInd[i] == 1) && (fabs(upper[i] - lower[i]) > etol)){
	      isLinkVarsFixed_ = false;
	      break;
	  }
      }
  }
  
  /* put the solution in order by integers first */

  int pos(0);

  for(i = 0; i < numElements; i++){
    index = indices[i];
    pos = binarySearch(0, lN - 1, index, lowerColInd);
    if(pos < 0){
      pos = binarySearch(0, uN - 1, index, upperColInd);
      if(mibs->solver()->isInteger(index)){
	  upperSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
      }
      else{
      upperSolutionOrd_[pos] = values[i];
      }
      optUpperSolutionOrd_[pos] = upperSolutionOrd_[pos];
    }
    else{
	if(mibs->solver()->isInteger(index)){
	    lowerSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	}
	else{
	    lowerSolutionOrd_[pos] = values[i];
	}
	optLowerSolutionOrd_[pos] = lowerSolutionOrd_[pos];	
    }
  }

  int solType(0);
  
  if(isLinkVarsIntegral_){
      std::vector<double> linkSol;
      for(i = 0; i < uN; i++){
	  index = upperColInd[i];
	  if(fixedInd[index] == 1){
	      linkSol.push_back(upperSolutionOrd_[i]);
	  }
      }
      if(model_->seenLinkingSolutions.find(linkSol) !=
	 model_->seenLinkingSolutions.end()){
	  isContainedInSetE_ = true;
	  solType = model_->seenLinkingSolutions.find(linkSol)->second.tag;
      }
  }
	 
  if(isContainedInSetE_){
      tagInSeenLinkingPool_ = static_cast<MibSSetETag>(solType);
  }
  if(tagInSeenLinkingPool_ == MibSSetETagVFIsInfeasible){
      LPSolStatus_ = MibSLPSolStatusInfeasible;
  }

  //steps 5-6
  if((isLinkVarsFixed_) && ((tagInSeenLinkingPool_ == MibSSetETagVFIsInfeasible) ||
			 (tagInSeenLinkingPool_ == MibSSetETagUBIsSolved))){
      useBilevelBranching_ = false;
      isProvenOptimal_ = false;
      shouldPrune_ = true;
  }

  //step 7
  if(!shouldPrune_){
      if(((tagInSeenLinkingPool_ == MibSSetETagVFIsFeasible)
	  || (tagInSeenLinkingPool_ == MibSSetETagUBIsSolved)) ||
	 ((!isContainedInSetE_) &&
	  (((branchPar == MibSBranchingStrategyLinking) &&
	    (isIntegral_) && (isLinkVarsFixed_)) ||
	   ((branchPar == MibSBranchingStrategyFractional)
	    && (isIntegral_)) ||
	   ((solveLowerWhenXYVarsInt == PARAM_ON) && (isIntegral_)) ||
	   ((solveLowerWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
	   ((solveLowerWhenLinkVarsInt == PARAM_ON) && (isLinkVarsIntegral_)) ||
	   ((solveLowerWhenLinkVarsFixed == PARAM_ON) && (isLinkVarsFixed_ ))))){
	  checkBilevelFeasiblity(mibs->isRoot_);
      }
  }
  if(cutStrategy == 1){
      useBilevelBranching_ = false;
  }
}


//#############################################################################
void
MibSBilevel::checkBilevelFeasiblity(bool isRoot)
{
    bool warmStartLL(model_->MibSPar_->entry
		     (MibSParams::warmStartLL));
    int maxThreadsLL(model_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(model_->MibSPar_->entry
		    (MibSParams::whichCutsLL));
    int probType(model_->MibSPar_->entry
		 (MibSParams::bilevelProblemType));
    std::string feasCheckSolver(model_->MibSPar_->entry
				(MibSParams::feasCheckSolver));
    MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
	(model_->MibSPar_->entry(MibSParams::branchStrategy));
    int computeBestUBWhenXVarsInt(model_->MibSPar_->entry
			      (MibSParams::computeBestUBWhenXVarsInt));
    int computeBestUBWhenLinkVarsInt(model_->MibSPar_->entry
			      (MibSParams::computeBestUBWhenLinkVarsInt));
    int computeBestUBWhenLinkVarsFixed(model_->MibSPar_->entry
			      (MibSParams::computeBestUBWhenLinkVarsFixed));
    int saveSeenLinkingSols(model_->MibSPar_->entry
			    (MibSParams::saveSeenLinkingSols));
    int lN(model_->lowerDim_); // lower-level dimension
    int uN(model_->upperDim_); // lower-level dimension
    int i(0), index(0), length(0), pos(0);
    int sizeFixedInd(model_->sizeFixedInd_);
    double etol(model_->etol_), objVal(0.0), lowerObj(0.0);
    int * fixedInd = model_->fixedInd_;
    int * lowerColInd = model_->getLowerColInd();
    int * upperColInd = model_->getUpperColInd();
    double *lowerSol = new double[lN];
    CoinFillN(lowerSol, lN, 0.0);

    std::vector<double> shouldStoreValues;

    const double * sol = model_->solver()->getColSolution();
    
    std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd_[i]);
	}
    }

    isProvenOptimal_ = true; 

    if(!isContainedInSetE_){
	
	//isProvenOptimal_ = true;
    
	if (warmStartLL && (feasCheckSolver == "SYMPHONY") && solver_){
	    solver_ = setUpModel(model_->getSolver(), false);
	}else{
	    if (solver_){
		delete solver_;
	    }
	    solver_ = setUpModel(model_->getSolver(), true);
	}
	
	OsiSolverInterface *lSolver = solver_;
	
	if(1)
	    lSolver->writeLp("lowerlevel");
	
	if (feasCheckSolver == "Cbc"){
	    dynamic_cast<OsiCbcSolverInterface *> 
		(lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
	}else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	    //dynamic_cast<OsiSymSolverInterface *> 
	    // (lSolver)->setSymParam("prep_level", -1);
	    
	    sym_environment *env = dynamic_cast<OsiSymSolverInterface *> 
		(lSolver)->getSymphonyEnvironment();
	    
	    if (warmStartLL){
		sym_set_int_param(env, "keep_warm_start", TRUE);
		if (probType == 1){ //Interdiction
		    sym_set_int_param(env, "should_use_rel_br", FALSE);
		    sym_set_int_param(env, "use_hot_starts", FALSE);
		    sym_set_int_param(env, "should_warmstart_node", TRUE);
		    sym_set_int_param(env, "sensitivity_analysis", TRUE);
		    sym_set_int_param(env, "sensitivity_bounds", TRUE);
		    sym_set_int_param(env, "set_obj_upper_lim", FALSE);
		}
	    }
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
	    lSolver->setHintParam(OsiDoReducePrint);
	    lSolver->messageHandler()->setLogLevel(0);
	    CPXENVptr cpxEnv = 
		dynamic_cast<OsiCpxSolverInterface*>(lSolver)->getEnvironmentPtr();
	    assert(cpxEnv);
	    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
	}
	
	//step 8
	if (warmStartLL && feasCheckSolver == "SYMPHONY"){
	    lSolver->resolve();
	    setWarmStart(lSolver->getWarmStart());
	}else{
	    lSolver->branchAndBound();
	}
  
	model_->counterVF_ ++;
	isLowerSolved_ = true;
  
	if(!lSolver->isProvenOptimal()){
	    LPSolStatus_ = MibSLPSolStatusInfeasible;
	    isProvenOptimal_ = false;
	    if(saveSeenLinkingSols){
	    //step 10
	    //Adding x_L to set E
		addSolutionToSetE(MibSSetETagVFIsInfeasible, shouldStoreValues, 0.0);
	    }
	}
	else{
	    //const double * sol = model_->solver()->getColSolution();
	    objVal = lSolver->getObjValue() * model_->getLowerObjSense();

	    objVal_ = objVal;

	    const double * values = lSolver->getColSolution();
	    
	    if(saveSeenLinkingSols){
		for(i = 0; i < lN; i++){
		    shouldStoreValues.push_back(values[i]);
		}
		
		//step 12
		//Adding x_L to set E  
		addSolutionToSetE(MibSSetETagVFIsFeasible, shouldStoreValues, objVal_);
		shouldStoreValues.clear();
	    }
	    else{
		for(i = 0; i < lN; i++){
		    lowerSol[i] = values[i];
		}
	    }
	    
	    MibSTreeNode * node = static_cast<MibSTreeNode *>(model_->activeNode_);
	    MibSTreeNode * parent =
		static_cast<MibSTreeNode *>(model_->activeNode_->getParent());
	    
	    if((!node->isBoundSet())
	       && (node->getIndex() != 0)){
		double parentBound = parent->getLowerUB();
		node->setLowerUB(parentBound);
		node->setIsBoundSet(true);
	    }
	    
	    if(objVal > node->getLowerUB()){
		node->setLowerUB(objVal);
		node->setIsBoundSet(true);
	    }
	}
    }

    //step 13
    if(((!saveSeenLinkingSols) && (isProvenOptimal_)) ||
       ((tagInSeenLinkingPool_ == MibSSetETagVFIsFeasible) ||
	(tagInSeenLinkingPool_ == MibSSetETagUBIsSolved))){

	OsiSolverInterface *UBSolver = 0;
	
	//double *lowerSol = new double[lN];
	//CoinFillN(lowerSol, lN, 0.0);

	if(saveSeenLinkingSols){
	    //get optimal value  of (VF) from solution pool
	    model_->it = model_->seenLinkingSolutions.find(linkSol);
	    objVal = model_->it->second.lowerObjVal1;
	    objVal_ = objVal;
	    for(i = 0; i < lN; i++){
		lowerSol[i] = model_->it->second.lowerSol1[i];
	    }  
	}
	lowerObj = getLowerObj(sol, model_->getLowerObjSense());
	
	if(isIntegral_){
	    assert((objVal - lowerObj) <= etol);
	}
	
	LPSolStatus_ = MibSLPSolStatusInfeasible;
	
	//step 15
	/** Current solution is bilevel feasible **/
	if((fabs(objVal - lowerObj) < etol) && (isIntegral_)){
	    LPSolStatus_ = MibSLPSolStatusFeasible;
	    useBilevelBranching_ = false;
	    shouldPrune_ = true;
	}
	if(!shouldPrune_){	
	    //step 18
	    if((tagInSeenLinkingPool_ != MibSSetETagUBIsSolved) &&
	       (((branchPar == MibSBranchingStrategyLinking) &&
		 (isIntegral_) && (isLinkVarsFixed_)) ||
		((computeBestUBWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
		((computeBestUBWhenLinkVarsInt == PARAM_ON)) ||
		((computeBestUBWhenLinkVarsFixed == PARAM_ON) && (isLinkVarsFixed_)))){  
		if(UBSolver){
		    delete UBSolver;
		}
		UBSolver = setUpUBModel(model_->getSolver(), objVal, true);
	    
#ifndef COIN_HAS_SYMPHONY
		dynamic_cast<OsiCbcSolverInterface *>
		    (UBSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
		dynamic_cast<OsiSymSolverInterface *>
		    (UBSolver)->setSymParam("prep_level", -1);
	    
		dynamic_cast<OsiSymSolverInterface *>
		    (UBSolver)->setSymParam("verbosity", -2);
	    
		dynamic_cast<OsiSymSolverInterface *>
		    (UBSolver)->setSymParam("max_active_nodes", 1);
#endif
		//step 19
		UBSolver->branchAndBound();
		model_->counterUB_ ++;
		isUBSolved_ = true;
		if(UBSolver->isProvenOptimal()){
		    haveHeurSolCand_ = true;
		    isProvenOptimal_ = true;
		    const double * valuesUB = UBSolver->getColSolution();
		    for(i = 0; i < uN + lN; i++){
			shouldStoreValues.push_back(valuesUB[i]);
			pos = binarySearch(0, uN - 1, i, upperColInd);
			if (pos >= 0){
			    if(UBSolver->isInteger(i)){
				optUpperSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
			    }
			    else{
				optUpperSolutionOrd_[pos] = (double) valuesUB[i];
			    }
			}
			else{
			    pos = binarySearch(0, lN - 1, i, lowerColInd);
			    if(UBSolver->isInteger(i)){
				optLowerSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
			    }
			    else{
				optLowerSolutionOrd_[pos] = (double) valuesUB[i];
			    }
			}
		    }
		    objVal = UBSolver->getObjValue() * model_->solver()->getObjSense();
		}
		else{
		    objVal = 10000000;
		    isProvenOptimal_ = false;
		}
		//step 22
		//Adding x_L to set E
		if(saveSeenLinkingSols){
		    addSolutionToSetE(MibSSetETagUBIsSolved, shouldStoreValues, objVal);
		}
		shouldStoreValues.clear();
	    
		//step 23
		if(isLinkVarsFixed_){
		    useBilevelBranching_ = false;
		    //isProvenOptimal_ = false;
		    shouldPrune_ = true;
		}	
	    }
	    else if((tagInSeenLinkingPool_ != MibSSetETagUBIsSolved) ||
		    ((!saveSeenLinkingSols) && (isUBSolved_))){
		haveHeurSolCand_ = true;
		for(i = 0; i < lN; i++){
		    index = lowerColInd[i];
		    if(model_->solver()->isInteger(index)){
			optLowerSolutionOrd_[i] = (double) floor(lowerSol[i] + 0.5);
		    }
		    else{
			optLowerSolutionOrd_[i] = (double) lowerSol[i];
		    }
		}
	    }
	}
    }
    delete [] lowerSol;
}
	    
//#############################################################################
void 
MibSBilevel::gutsOfDestructor()
{

  if(optUpperSolutionOrd_) delete [] optUpperSolutionOrd_;
  if(optLowerSolutionOrd_) delete [] optLowerSolutionOrd_;
  if(upperSolutionOrd_) delete [] upperSolutionOrd_;
  if(lowerSolutionOrd_) delete [] lowerSolutionOrd_;
  //delete heuristic_;
}


//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpUBModel(OsiSolverInterface * oSolver, double objValLL,
			      bool newOsi, const double *lpSol)
{
    int probType =
	model_->MibSPar_->entry(MibSParams::bilevelProblemType);

    bool warmStartLL =
	model_->MibSPar_->entry(MibSParams::warmStartLL);

    bool doDualFixing =
	model_->MibSPar_->entry(MibSParams::doDualFixing);

    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);

    OsiSolverInterface * nSolver;

    if (!lpSol){
	lpSol = oSolver->getColSolution();
    }

    int * fixedInd = model_->fixedInd_;

    int i(0), j(0), index1(0);
    double tmp(0.0);

    int uCols(model_->getUpperDim());
    int uRows(model_->getOrigUpperRowNum());
    int lRows(model_->getLowerRowNum());
    int lCols(model_->getLowerDim());
    double objSense(model_->getLowerObjSense());
    double uObjSense(1);
    int * uColIndices = model_->getUpperColInd();
    int * lColIndices = model_->getLowerColInd();
    int * uRowIndices = model_->getOrigUpperRowInd();
    int * lRowIndices = model_->getLowerRowInd();
    const double * origColLb = model_->getOrigColLb();
    const double * origColUb = model_->getOrigColUb();
    double * lObjCoeffs = model_->getLowerObjCoeffs();
    const double * uObjCoeffs = oSolver->getObjCoefficients();

    CoinPackedMatrix * matrix = NULL;
    matrix = model_->origConstCoefMatrix_;

    int rowNum(model_->getNumOrigCons() + 1);
    int colNum(model_->getNumOrigVars());
    double * rowUb = new double[rowNum];
    double * rowLb = new double[rowNum];
    double * colUb = new double[colNum];
    double * colLb = new double[colNum];

    CoinZeroN(colUb, colNum);
    CoinZeroN(colLb, colNum);

    /** Set the row bounds **/
    for(i = 0; i < rowNum-1; i++){
	rowLb[i] = model_->getOrigRowLb()[i];
	rowUb[i] = model_->getOrigRowUb()[i];
    }

    for(i = 0; i < uCols; i++){
	index1 = uColIndices[i];
	colLb[index1] = origColLb[index1];
	colUb[index1] = origColUb[index1];
	//if(isLinkVarsFixed_ == false){
	    if(fixedInd[index1] == 1){
		colLb[index1] = floor(lpSol[index1] + 0.5);
		colUb[index1] = colLb[index1];
	    }
	    //}
    }

    //sahar:To Do: check it again
    for(i = 0; i < lCols; i++){
	index1 = lColIndices[i];
	colLb[index1] = origColLb[index1];
	colUb[index1] = origColUb[index1];
    }
    
    if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"setUpModel", "MibsBilevel");
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"setUpModel", "MibsBilevel");
#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"setUpModel", "MibsBilevel");
    }
    

    int * integerVars = new int[colNum];
    double * objCoeffs = new double[colNum];
    double * newRow = new double[lCols];
    CoinFillN(integerVars, colNum, 0);
    CoinFillN(objCoeffs, colNum, 0.0);
    CoinFillN(newRow, lCols, 0.0);

    int intCnt(0);

    /** Fill in array of integer variables **/
    for(i = 0; i < colNum; i++){
	if(oSolver->isInteger(i)){
	    integerVars[intCnt] = i;
	    intCnt ++;
	}
    }

    CoinDisjointCopyN(uObjCoeffs, colNum, objCoeffs);

    CoinDisjointCopyN(lObjCoeffs, lCols, newRow);

    for(i = 0; i < lCols; i++){
	newRow[i] = newRow[i] * objSense;
    }

    rowUb[rowNum-1] = objValLL * objSense;
    rowLb[rowNum-1] = -1 * (oSolver->getInfinity());

    CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
    newMat->setDimensions(0, colNum);

    for(i = 0; i < rowNum - 1; i++){
	CoinPackedVector row;
	for(j = 0; j < colNum; j++){
	    tmp = matrix->getCoefficient(i, j);
	    row.insert(j, tmp);
	}
	newMat->appendRow(row);
    }

    CoinPackedVector row;
    for(i = 0; i < lCols; i++){
	index1 = lColIndices[i];
	tmp = newRow[i];
	row.insert(index1, tmp);
    }
    newMat->appendRow(row);

    nSolver->loadProblem(*newMat, colLb, colUb,
			 objCoeffs, rowLb, rowUb);

    for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
    }

    nSolver->setObjSense(uObjSense); //1 min; -1 max

    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    delete [] integerVars;


    return nSolver;

}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpModel(OsiSolverInterface * oSolver, bool newOsi,
			const double *lpSol)
{

  /** Create lower-level model with fixed upper-level vars **/

  int probType =
    model_->MibSPar_->entry(MibSParams::bilevelProblemType);

  bool warmStartLL =
    model_->MibSPar_->entry(MibSParams::warmStartLL);

  bool doDualFixing =
    model_->MibSPar_->entry(MibSParams::doDualFixing);

  std::string feasCheckSolver =
    model_->MibSPar_->entry(MibSParams::feasCheckSolver);

  OsiSolverInterface * nSolver;

  double etol(model_->etol_);
  int i(0), j(0), index1(0), index2(0);
  int uCols(model_->getUpperDim());
  int lRows(model_->getLowerRowNum());
  int lCols(model_->getLowerDim());
  int * uColIndices = model_->getUpperColInd();
  int * lColIndices = model_->getLowerColInd();
  int * lRowIndices = model_->getLowerRowInd();

  double objSense(model_->getLowerObjSense());  
  double * lObjCoeffs = model_->getLowerObjCoeffs();
     
  const CoinPackedMatrix * matrix = oSolver->getMatrixByRow();
  double coeff(0.0);
  if (!lpSol){
     lpSol = oSolver->getColSolution();   
  }
  const double * origRowLb = model_->getOrigRowLb();
  const double * origRowUb = model_->getOrigRowUb();
  const double * origColLb = model_->getOrigColLb();
  const double * origColUb = model_->getOrigColUb();
  double * rowUb = new double[lRows];
  double * rowLb = new double[lRows];
  double * colUb = new double[lCols];
  double * colLb = new double[lCols];
  //CoinZeroN(rowUb, lRows);
  //CoinZeroN(rowLb, lRows);
  //CoinZeroN(colUb, lCols);
  //CoinZeroN(colLb, lCols);

  /** Set the row bounds **/
  
  for(i = 0; i < lRows; i++){
      index1 = lRowIndices[i];
      rowLb[i] = oSolver->getRowLower()[index1];
      rowUb[i] = oSolver->getRowUpper()[index1];
  }
     
  for(i = 0; i < lCols; i++){
     colLb[i] = origColLb[lColIndices[i]];
     colUb[i] = origColUb[lColIndices[i]];
  }
  
  if (newOsi){
     if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
     }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"setUpModel", "MibsBilevel");
#endif
     }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"setUpModel", "MibsBilevel");
#endif
     }else{
	throw CoinError("Unknown solver chosen",
			"setUpModel", "MibsBilevel");
     }
	
     int * integerVars = new int[lCols];
     double * objCoeffs = new double[lCols];
     
     CoinFillN(integerVars, lCols, 0);
     //CoinZeroN(objCoeffs, lCols);
     
     int intCnt(0);
  
     /** Fill in array of lower-level integer variables **/
     
     for(i = 0; i < lCols; i++){
	index1 = lColIndices[i];
	if(oSolver->isInteger(index1)){
	   integerVars[intCnt] = i;
	   intCnt++;
	}
     }
     
     CoinDisjointCopyN(lObjCoeffs, lCols, objCoeffs);
     
     CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
     newMat->setDimensions(0, lCols);
     double tmp(0.0);
     
     /*
       for(i = 0; i < lRows; i++){
       CoinPackedVector row;
       index1 = lRowIndices[i];
       start = matStarts[index1];
       end = start + matrix->getVectorSize(index1);
       for(j = start; j < end; j++){
       index2 = matIndices[j];
       //tmp = findIndex(index, lCols, lColIndices);
       tmp = binarySearch(0, lCols - 1, index2, lColIndices);
       if(tmp > -1)
       row.insert(tmp, matElements[j]);
       }
       newMat->appendRow(row);
       }
     */
     for(i = 0; i < lRows; i++){
	CoinPackedVector row;
	index1 = lRowIndices[i];
	for(j = 0; j < lCols; j++){
	   index2 = lColIndices[j];
	   tmp = matrix->getCoefficient(index1, index2);
	   row.insert(j, tmp);
	}
	newMat->appendRow(row);
     }

     /*
       nSolver->assignProblem(newMat, colLb, colUb,
       objCoeffs, rowLb, rowUb);
     */
     
     nSolver->loadProblem(*newMat, colLb, colUb,
			  objCoeffs, rowLb, rowUb);
     
     for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
     }
     //nSolver->setInteger(integerVars, intCnt);
     
     nSolver->setObjSense(objSense); //1 min; -1 max
     
     nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

#if 0
     if(0){
	dynamic_cast<OsiCbcSolverInterface *> 
	   (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
     }
     else{
#ifdef COIN_HAS_SYMPHONY
	dynamic_cast<OsiSymSolverInterface *> 
	   (nSolver)->setSymParam("prep_level", -1);
	
	dynamic_cast<OsiSymSolverInterface *> 
	   (nSolver)->setSymParam("verbosity", -2);
	
	dynamic_cast<OsiSymSolverInterface *> 
	   (nSolver)->setSymParam("max_active_nodes", 1);
#endif
     }
#endif
     delete [] integerVars;

  }else{
     nSolver = solver_;
  }
  
#define SYM_VERSION_IS_WS strcmp(SYMPHONY_VERSION, "WS")  

#if SYMPHONY_VERSION_IS_WS
  if (feasCheckSolver == "SYMPHONY" && probType == 1 && warmStartLL &&
      !newOsi && doDualFixing){ //Interdiction

     /** Get upper bound from best known (feasible) lower level solution and 
	 try to fix additional variables by sensitivity analysis **/

     std::vector<std::pair<AlpsKnowledge*, double> > solutionPool;
     model_->getKnowledgeBroker()->
	getAllKnowledges(AlpsKnowledgeTypeSolution, solutionPool);

     const double * sol; 
     double objval, Ub(objSense*nSolver->getInfinity());
     BlisSolution* blisSol;
     std::vector<std::pair<AlpsKnowledge*, double> >::const_iterator si;
     for (si = solutionPool.begin(); si != solutionPool.end(); ++si){
	blisSol = dynamic_cast<BlisSolution*>(si->first);
	sol = blisSol->getValues();
	for (i = 0; i < uCols; i++){
	   if (lpSol[uColIndices[i]] > 1 - etol &&
	       sol[lColIndices[i]] > 1-etol){
	      break;
	   }
	}
	if (i == uCols && -objSense*blisSol->getQuality() < Ub){
	   Ub = -objSense*blisSol->getQuality();
	}
     }

     /** Figure out which variables get fixed by the upper level solution **/

     int *newUbInd = new int[uCols];
     int *newLbInd = new int[uCols];
     double *newUbVal = new double[uCols];
     double *newLbVal = new double[uCols];
     double newLb;
     for (i = 0; i < uCols; i++){
	newUbInd[i] = uColIndices[i];
	newLbInd[i] = uColIndices[i];
	newLbVal[i] = 0;
	if (lpSol[uColIndices[i]] > 1 - etol){
	   newUbVal[i] = 0;
	}else{
	   newUbVal[i] = 1;
	}
     }
	
     /** If we found an upper bound, then do the dual fixing **/

     if (Ub < objSense*nSolver->getInfinity()){
	sym_environment *env = 
	   dynamic_cast<OsiSymSolverInterface *>(nSolver)->getSymphonyEnvironment();

	for (i = 0; i < uCols; i++){
	   if (newUbVal[i] == 1){
	      // Try fixing it to zero
	      newUbVal[i] = 0; 
	      sym_get_lb_for_new_rhs(env, 0, NULL, NULL, uCols, newLbInd,
				     newLbVal, uCols, newUbInd, newUbVal,
				     &newLb);
	      if (objSense*newLb > Ub + etol){
		 //Victory! This variable can be fixed to 1 permanently
		 newLbVal[i] = 1;
	      }
	      //Set upper bound back to 1
	      newUbVal[i] = 1;
	      
	      if (newLbVal[i] == 0){
		 // Try fixing it to one
		 newLbVal[i] = 1;
		 sym_get_lb_for_new_rhs(env, 0, NULL, NULL, uCols, newLbInd,
					newLbVal, uCols, newUbInd, newUbVal,
					&newLb);
		 if (objSense*newLb > Ub + etol){
		    //Victory! This variable can be fixed to 0 permanently
		    newUbVal[i] = 0;
		 }
		 newLbVal[i] = 0;
	      }
	   }
	   
	}
     }
     
     /** Now set the row bounds to account for fixings **/
     
     /** This is probably very frgaile. Assuming interdiction 
	 rows come last. It would be better to set variable 
	 bounds directly, but this doesn't seem to work right now. **/
     int iRowStart = lRows-uCols; 

#if 1
     for(i = iRowStart; i < lRows; i++){
	nSolver->setRowLower(i, newLbVal[i-iRowStart]);
	nSolver->setRowUpper(i, newUbVal[i-iRowStart]);
     }
#else
     for(i = 0; i < uCols; i++){
	nSolver->setColLower(i, newLbVal[i]);
	nSolver->setColUpper(i, newUbVal[i]);
     }
#endif     
     delete[] newUbInd;
     delete[] newLbInd;
     delete[] newUbVal;
     delete[] newLbVal;
  }else{
#endif
     //FIXME: NEED TO GET ROW SENSE HERE
     
     /** Get contribution of upper-level columns **/
     
     double * upComp = new double[lRows];
     CoinFillN(upComp, lRows, 0.0);
     
     for(i = 0; i < lRows; i++){
	index1 = lRowIndices[i];
	for(j = 0; j < uCols; j++){
	   index2 = uColIndices[j];
	   coeff = matrix->getCoefficient(index1, index2);
	   if (coeff != 0){
	      upComp[i] += coeff * lpSol[index2];
	   }
	}
     }
     
     /** Correct the row bounds to account for fixed upper-level vars **/
     
     for(i = 0; i < lRows; i++){
	nSolver->setRowLower(i, rowLb[i] - upComp[i]);
	nSolver->setRowUpper(i, rowUb[i] - upComp[i]);
     }
     
     delete [] upComp;

#if SYMPHONY_VERSION_IS_WS
  }
#endif

  //I don't think this is needed
  //if(!getWarmStart())
  //  setWarmStart(nSolver->getWarmStart());

  return nSolver;

}

//#############################################################################
int
MibSBilevel::findIndex(int index, int size, int * indices)
{

  int i(0);
  int pos(-1);

  for(i = 0; i < size; i++){
    if(indices[i] == index)
       pos = i;
  }

  return pos;

}

//#############################################################################
int 
MibSBilevel::binarySearch(int start, int stop, int index, int * indexArray)
{
   int i(0);
   int pos(-1);

  //=========================================================================
  // If the list is short enough, finish with sequential search
  // Otherwise, call binary search recursively
  //=========================================================================

   //FIXME: CHANGED THIS 7/15
   //NOW USES ONLY SEQUENTIAL SEARCH
   //BINARY REQUIRES ORDERED ARRAY
   //OUR ARRAY IS ORDERED BY INTEGERS FIRST, NOT INDICES


   //if((stop - start) < 4){
   if(1){
      for(i = start; i < stop + 1; ++i){
	 if(indexArray[i] == index){
	    pos = i;
	    break;
	 }
      }
      
      return pos;
   }
   else{ 
      
      int midpoint((start + stop)/2);
      int val(indexArray[midpoint]);
      
      if(val == index){
	 pos = midpoint;
	 return pos;
      }
      else if(val > index){
	 pos = binarySearch(start, midpoint - 1, index, indexArray);
      }
      else{
	 pos = binarySearch(midpoint + 1, stop, index, indexArray);
      }
   }
   return pos;
}



//#############################################################################
double
MibSBilevel::getLowerObj(const double * sol, double objSense)
{

   int lCols(model_->getLowerDim());
   int * lColIndices = model_->getLowerColInd();
   double * lObjCoeffs = model_->getLowerObjCoeffs();

   int i(0), index(0);
   double objVal(0.0);


   for(i = 0; i < lCols; i++){
      index = lColIndices[i];
      if(0){
	std::cout << "sol[" << index << "]: " << sol[index] << std::endl;
	std::cout << "lObjCoeff[" << i << "]: " << lObjCoeffs[i] << std::endl;
      }      
      objVal += lObjCoeffs[i] * sol[index];
   }

   return objVal * objSense;

}
//#############################################################################
void
    MibSBilevel::addSolutionToSetE(MibSSetETag solTag, std::vector<double>
				   &shouldStoreValues, double objValue)
{
    int i(0),index(0);
    int uN(model_->upperDim_);
    int lN(model_->lowerDim_);
    int * upperColInd = model_->getUpperColInd();
    int * fixedInd = model_->fixedInd_;
    int solType = static_cast<int>(solTag);
    
    std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd_[i]);
	}
    }

    tagInSeenLinkingPool_ = solTag;
    model_->linkingSolution.lowerSol1.push_back(0);
    model_->linkingSolution.UBSol1.push_back(0);
    model_->linkingSolution.lowerSol1.clear();
    model_->linkingSolution.UBSol1.clear();
    
    switch(solTag){
    case MibSSetETagVFIsInfeasible:
	{
	    model_->linkingSolution.tag = solType;
	    model_->linkingSolution.lowerObjVal1 = 0.0;
	    model_->linkingSolution.UBObjVal1 = 0.0;
	    model_->linkingSolution.lowerSol1.push_back(0);
	    model_->linkingSolution.UBSol1.push_back(0);
	    model_->seenLinkingSolutions[linkSol] = model_->linkingSolution;
	    break;
	}
    case MibSSetETagVFIsFeasible:
	{
	    model_->linkingSolution.tag = solType;
	    model_->linkingSolution.lowerObjVal1 = objValue;
	    model_->linkingSolution.UBObjVal1 = 0.0;
	    for(i = 0; i < lN; i++){
		model_->linkingSolution.lowerSol1.push_back(shouldStoreValues[i]);
	    }
	    model_->linkingSolution.UBSol1.push_back(0);
	    model_->seenLinkingSolutions[linkSol] = model_->linkingSolution;
	    break;
	}
    case MibSSetETagUBIsSolved:
	{
	    //model_->it = model_->seenLinkingSolutions.find(linkSol);
	    //model_->it->second.tag = MibSSetETagUBIsSolved;
	    //model_->it->second.UBObjVal1 = objValue;
	    model_->seenLinkingSolutions[linkSol].tag = MibSSetETagUBIsSolved;
	    model_->seenLinkingSolutions[linkSol].UBObjVal1 = objValue;
	    if(isProvenOptimal_){
		//model_->it->second.UBSol1.clear();
		model_->seenLinkingSolutions[linkSol].UBSol1.clear();
		for(i = 0; i < uN + lN; i++){
		    //model_->it->second.UBSol1.push_back(shouldStoreValues[i]);
		    model_->seenLinkingSolutions[linkSol].UBSol1.push_back(shouldStoreValues[i]);
		}
	    }
	    break;
	}
    }
}
	    
	
	
	    






    
 
