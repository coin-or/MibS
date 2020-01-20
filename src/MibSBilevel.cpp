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

#include "CoinPackedVector.hpp"
#include "OsiCbcSolverInterface.hpp"

#include "MibSBilevel.hpp"
//#include "MibSModel.hpp"
#include "MibSTreeNode.hpp"
#include "MibSSolution.hpp"
//#include "MibSHeuristic.hpp"
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
//to run in parallel, uncomment this (should be modified later)
//#define _OPENMPMIBS 1
#ifdef _OPENMPMIBS
#include "omp.h"
#endif

//#############################################################################
MibSSolType
MibSBilevel::createBilevel(CoinPackedVector* sol, 
			   MibSModel *mibs)
{

  /** Splits sol into two parts, upper- and lower-level components **/

  if(!mibs) return MibSNoSol;
  
  model_ = mibs;

  int numScenarios(model_->getNumScenarios());

  if(numScenarios == 1){
      heuristic_ = new MibSHeuristic(mibs);
  }

  int i(0),j(0);
  int N(model_->numVars_);
  int truncLN(model_->getTruncLowerDim());
  int lN(model_->lowerDim_); // lower-level dimension
  int uN(model_->upperDim_); // upper-level dimension
  double etol(model_->BlisPar()->entry(BlisParams::integerTol));
  
  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
      (model_->MibSPar_->entry(MibSParams::branchStrategy));

  int solveSecondLevelWhenXYVarsInt(model_->MibSPar_->entry
			   (MibSParams::solveSecondLevelWhenXYVarsInt));
  int solveSecondLevelWhenXVarsInt(model_->MibSPar_->entry
			  (MibSParams::solveSecondLevelWhenXVarsInt));
  int solveSecondLevelWhenLVarsInt(model_->MibSPar_->entry
			  (MibSParams::solveSecondLevelWhenLVarsInt));
  int solveSecondLevelWhenLVarsFixed(model_->MibSPar_->entry
			    (MibSParams::solveSecondLevelWhenLVarsFixed));
  int cutStrategy(model_->MibSPar_->entry
		  (MibSParams::cutStrategy));

  int useLinkingSolutionPool(model_->MibSPar_->entry
			   (MibSParams::useLinkingSolutionPool));

  int maxActiveNodes(model_->MibSPar_->entry
		     (MibSParams::maxActiveNodes)); 

  MibSSolType storeSol(MibSNoSol);

  //saharSto: fix it later
  bool shouldStoreSolution(false);

  if(numScenarios == 1){
      shouldStoreSolution = true;
  }
  
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

  if(!optUpperSolutionOrd_)
      optUpperSolutionOrd_ = new double[uN];

  if(!optLowerSolutionOrd_)
      optLowerSolutionOrd_ = new double[lN];

  CoinZeroN(upperSolutionOrd_, uN);
  CoinZeroN(optUpperSolutionOrd_, uN);
  CoinZeroN(optLowerSolutionOrd_, lN);
  
  if(shouldStoreSolution == true){
      if(!lowerSolutionOrd_)
	  lowerSolutionOrd_ = new double[lN];
  
      CoinZeroN(lowerSolutionOrd_, lN);
  }
  
  isIntegral_ = true;
  isUpperIntegral_ = true;
  isLinkVarsIntegral_ = true;
  LPSolStatus_ = MibSLPSolStatusUnknown;
  isLinkVarsFixed_ = true;
  shouldPrune_ = false;
  isLowerSolved_ = false;
  isUBSolved_ = false;
  isContainedInLinkingPool_ = false;
  useBilevelBranching_ = true;
  isProvenOptimal_ = true;
  tagInSeenLinkingPool_ = MibSLinkingPoolTagIsNotSet;

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

  if((shouldStoreSolution == true) && (numScenarios == 1)){
      for (i = 0; i < numElements; i++){
	  index = indices[i];
          pos = binarySearch(0, lN - 1, index, lowerColInd);
	  if (pos < 0){
	      pos = binarySearch(0, uN - 1, index, upperColInd);
	      if ((mibs->solver()->isInteger(index)) &&
		  (((values[i] - floor(values[i])) < etol) ||
		   ((ceil(values[i]) - values[i]) < etol))){
		  upperSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	      }else{
		  upperSolutionOrd_[pos] = values[i];
	      }
	      optUpperSolutionOrd_[pos] = upperSolutionOrd_[pos];
	  }else{
	      if ((mibs->solver()->isInteger(index)) &&
		  (((values[i] - floor(values[i])) < etol) ||
		   ((ceil(values[i]) - values[i]) < etol))){
		  lowerSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	      }else{
		  lowerSolutionOrd_[pos] = values[i];
	      }
	      optLowerSolutionOrd_[pos] = lowerSolutionOrd_[pos];	
	  }
      }
  }
  else{
      for (i = 0; i < numElements; i++){
	  index = indices[i];
	  if(numScenarios == 1){
	      pos = binarySearch(0, lN - 1, index, lowerColInd);
	  }
	  else{
	      if(index < uN){
		  pos = -1;
	      }
	      else{
		  pos = index - uN;
	      }
	  }
	  if(pos < 0){
	      pos = binarySearch(0, uN - 1, index, upperColInd);
	      if ((mibs->solver()->isInteger(index)) &&
		  (((values[i] - floor(values[i])) < etol) ||
		   ((ceil(values[i]) - values[i]) < etol))){
		  upperSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	      }
	      else{
		  upperSolutionOrd_[pos] = values[i];
	      }
	      optUpperSolutionOrd_[pos] = upperSolutionOrd_[pos];
	  }
	  else{
	      if ((mibs->solver()->isInteger(index)) &&
		  (((values[i] - floor(values[i])) < etol) ||
		   ((ceil(values[i]) - values[i]) < etol))){
		  optLowerSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	      }
	      else{
		  optLowerSolutionOrd_[pos] = values[i];
	      }
	      if(shouldStoreSolution == true){
		  lowerSolutionOrd_[pos] = lowerSolutionOrd_[pos];
	      }
	  }
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
	  isContainedInLinkingPool_ = true;
	  solType = model_->seenLinkingSolutions.find(linkSol)->second.tag;
      }
  }
	 
  if(isContainedInLinkingPool_){
      tagInSeenLinkingPool_ = static_cast<MibSLinkingPoolTag>(solType);
  }
  if(tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsInfeasible){
      LPSolStatus_ = MibSLPSolStatusInfeasible;
  }

  //steps 5-6
  if((isLinkVarsFixed_) && ((tagInSeenLinkingPool_ ==
			     MibSLinkingPoolTagLowerIsInfeasible) ||
			 (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){
      useBilevelBranching_ = false;
      shouldPrune_ = true;
  }

#ifdef _OPENMPMIBS
  if((maxActiveNodes == 1) || (numScenarios == 1)){
    throw CoinError("Sequential execution with OpenMp!",
		    "createBilevel", "MibSBilevel");
  }
#else
  if((maxActiveNodes > 1) && (numScenarios > 1)){
    throw CoinError("No OpenMp for parallel execution",
		    "createBilevel", "MibSBilevel");
  }
#endif

  //step 7
  if(!shouldPrune_){
      if(((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible)
	  || (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved)) ||
	 ((!isContainedInLinkingPool_) &&
	  (((branchPar == MibSBranchingStrategyLinking) &&
	    (isIntegral_) && (isLinkVarsFixed_)) ||
	   ((branchPar == MibSBranchingStrategyFractional)
	    && (isIntegral_)) ||
	   ((solveSecondLevelWhenXYVarsInt == PARAM_ON) && (isIntegral_)) ||
	   ((solveSecondLevelWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
	   ((solveSecondLevelWhenLVarsInt == PARAM_ON) && (isLinkVarsIntegral_)) ||
	   ((solveSecondLevelWhenLVarsFixed == PARAM_ON) && (isLinkVarsFixed_ ))))){
#ifdef _OPENMPMIBS
	storeSol = checkBilevelFeasiblityParallel(mibs->isRoot_);
#else
	storeSol = checkBilevelFeasiblity(mibs->isRoot_);
#endif
      }
  }
  if(cutStrategy == 1){
      useBilevelBranching_ = false;
  }

  if(numScenarios == 1){
      heuristic_->findHeuristicSolutions();

      delete heuristic_;
  }

  return storeSol;
}


//#############################################################################
MibSSolType
MibSBilevel::checkBilevelFeasiblity(bool isRoot)
{
    OsiSolverInterface *oSolver = model_->getSolver();
    
    int numScenarios(model_->getNumScenarios());
    
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
    int computeBestUBWhenLVarsInt(model_->MibSPar_->entry
				     (MibSParams::computeBestUBWhenLVarsInt));
    int computeBestUBWhenLVarsFixed(model_->MibSPar_->entry
				    (MibSParams::computeBestUBWhenLVarsFixed));
    int useLinkingSolutionPool(model_->MibSPar_->entry
			    (MibSParams::useLinkingSolutionPool));
    double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));
    double remainingTime(0.0), startTimeVF(0.0), startTimeUB(0.0);
    MibSSolType storeSol(MibSNoSol);
    int truncLN(model_->truncLowerDim_);;
    int lN(model_->lowerDim_); // lower-level dimension
    int uN(model_->upperDim_); // lower-level dimension
    int i(0), j(0), k(0),  index(0), length(0), pos(0), begPos(0);
    int lpStat;
    int sizeFixedInd(model_->sizeFixedInd_);
    double etol(model_->etol_), objVal(0.0), lowerObj(0.0);
    int * fixedInd = model_->fixedInd_;
    int * lowerColInd = model_->getLowerColInd();
    int * upperColInd = model_->getUpperColInd();
    CoinPackedMatrix *truncMatrixG2 = NULL;
    double *multA2XOpt = NULL;
    OsiSolverInterface *UBSolver = 0;

    //saharSto: fix it later
    if(numScenarios == 1){
	truncLN = lN;
    }
	

    //saharSto: fix it later
    bool shouldStoreSolution(false);

    if(numScenarios == 1){
	shouldStoreSolution = true;
    }

    //saharSto: Check these changes
    double *lowerSol = NULL;
    std::vector<double> shouldStoreValuesLowerSol(1);
    std::vector<double> shouldStoreValuesUBSol(1);
    std::vector<double> shouldStoreObjValues;
    

    if(shouldStoreSolution == true){
        shouldStoreValuesLowerSol.resize(lN);
        shouldStoreValuesUBSol.resize(lN + uN);
    }

    const double * sol = model_->solver()->getColSolution();
    
    std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd_[i]);
	}
    }

    isProvenOptimal_ = true; 

    if(!isContainedInLinkingPool_){
	for(i = 0; i < numScenarios; i++){
            remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();
	    if(remainingTime <= etol){
		shouldPrune_ = true;
		storeSol = MibSNoSol;
		goto TERM_CHECKBILEVELFEAS;
	    }
	
	//isProvenOptimal_ = true;
    
	/*if (warmStartLL && (feasCheckSolver == "SYMPHONY") && solver_){
	    solver_ = setUpModel(model_->getSolver(), false);
	}else{
	    //if (solver_){
	//	delete solver_;
	   // }
	    solver_ = setUpModel(model_->getSolver(), true);
	}*/

	    if(warmStartLL && (feasCheckSolver == "SYMPHONY") && lSolver_){
		lSolver_ = setUpModel(oSolver, false, i);
	    }
	    else{
		if(lSolver_){
		    delete lSolver_;
		}
	    	lSolver_ = setUpModel(oSolver, true, i);
	     }
	    
	    OsiSolverInterface *lSolver = lSolver_;

	    if(0)
		lSolver->writeLp("lSolver");

	    remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();
	    if(remainingTime <= etol){
		shouldPrune_ = true;
	        storeSol = MibSNoSol;
	        goto TERM_CHECKBILEVELFEAS;
	    }
	
	    remainingTime = CoinMax(remainingTime, 0.00);

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
		lSolver->setHintParam(OsiDoReducePrint);
	        lSolver->messageHandler()->setLogLevel(0);
		CPXENVptr cpxEnv =
		    dynamic_cast<OsiCpxSolverInterface*>(lSolver)->getEnvironmentPtr();
	        assert(cpxEnv);
	        CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	        CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
		CPXsetintparam(cpxEnv, CPX_PARAM_CLOCKTYPE, 1);
		CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, remainingTime);
#endif
	    }

	    //step 8
	    if (warmStartLL && feasCheckSolver == "SYMPHONY"){
		lSolver->resolve();
	        setWarmStart(lSolver->getWarmStart());
	    }else{
		startTimeVF = model_->broker_->subTreeTimer().getTime();
	        lSolver->branchAndBound();
		if(0)
		    lSolver->writeLp("lSolverAfter");
	        model_->timerVF_ += model_->broker_->subTreeTimer().getTime() - startTimeVF;
	    }
  
	    model_->counterVF_ ++;

	    if(i == numScenarios - 1){
		isLowerSolved_ = true;
	    }
	
	    if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
		if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface *>
					     (lSolver)->getSymphonyEnvironment())){
		shouldPrune_ = true;
	        storeSol = MibSNoSol;
	        goto TERM_CHECKBILEVELFEAS;
		}
#endif
	    }
	    else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
		lpStat = CPXgetstat(dynamic_cast<OsiCpxSolverInterface*>
				    (lSolver)->getEnvironmentPtr(),
				    dynamic_cast<OsiCpxSolverInterface*>
				    (lSolver)->getLpPtr());
		if((lpStat == CPXMIP_TIME_LIM_FEAS) ||
		   (lpStat == CPXMIP_TIME_LIM_INFEAS)){
		    shouldPrune_ = true;
		    storeSol = MibSNoSol;
		    goto TERM_CHECKBILEVELFEAS;
		}
#endif
	    }
	    
	    if(!lSolver->isProvenOptimal()){
		LPSolStatus_ = MibSLPSolStatusInfeasible;
		isProvenOptimal_ = false;
	        if(useLinkingSolutionPool){
		    //step 10
		    //Adding x_L to set E
		    shouldStoreObjValues.push_back(0);
		    addSolutionToSeenLinkingSolutionPool
			(MibSLinkingPoolTagLowerIsInfeasible, shouldStoreValuesLowerSol,
			 shouldStoreObjValues);
		    shouldStoreObjValues.clear();
		}
	        if(isLinkVarsFixed_){
		    useBilevelBranching_ = false;
	            shouldPrune_ = true;
		}
	        break;
	    }
	    else{
		//const double * sol = model_->solver()->getColSolution();
		objVal = lSolver->getObjValue() * model_->getLowerObjSense();

	        //objVal_ = objVal;

		if((i == 0) && (!objValVec_.empty())){
		    objValVec_.clear();
		}

		objValVec_.push_back(objVal);
	        shouldStoreObjValues.push_back(objVal);

	        const double * values = lSolver->getColSolution();

		if(useLinkingSolutionPool){
		    begPos = i * truncLN;
		    if(shouldStoreSolution == true){
		        std::copy(values, values + lN, shouldStoreValuesLowerSol.begin() + begPos);
		    }
		    else{
			if(lowerSol == NULL){
			    lowerSol = new double[lN];
			    CoinFillN(lowerSol, lN, 0.0);
			}
			CoinDisjointCopyN(values, truncLN, lowerSol + begPos);
		    }

		    //step 12
		    //Adding x_L to set E
		    if(i == numScenarios - 1){
			addSolutionToSeenLinkingSolutionPool
			    (MibSLinkingPoolTagLowerIsFeasible, shouldStoreValuesLowerSol,
			     shouldStoreObjValues);
			shouldStoreValuesLowerSol.clear();
		        shouldStoreObjValues.clear();
		    }
		}
		else{
		    if(lowerSol == NULL){
			lowerSol = new double[lN];
		        CoinFillN(lowerSol, lN, 0.0);
		    }
		    //memcpy(lowerSol, values, sizeof(double) * lN);
		    begPos = i * truncLN;
		    CoinDisjointCopyN(values, truncLN, lowerSol + begPos);
		}

		if(i == numScenarios - 1){
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
	    /*if (!warmStartLL){
              delete solver_;
	    }*/
	}
    }

    //step 13
    if(((!useLinkingSolutionPool) && (isProvenOptimal_)) ||
       ((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible) ||
	(tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){
	
	//double *lowerSol = new double[lN];
	//CoinFillN(lowerSol, lN, 0.0);

	if(useLinkingSolutionPool){
	    //get optimal value  of (VF) from solution pool
	    //model_->it = seenLinkingSolutions.find(linkSol);
	    //objVal = model_->it->second.lowerObjVal1;
	    shouldStoreObjValues = model_->seenLinkingSolutions[linkSol].lowerObjValue; 
	    //objVal = seenLinkingSolutions.find(linkSol).
	    //objVal_ = objVal;
	    //objVal_ = shouldStoreObjValues[0];
	    if(!objValVec_.empty()){
		objValVec_.clear();
	    }
	    objValVec_ = model_->seenLinkingSolutions[linkSol].lowerObjValue;
          

	    if(shouldStoreSolution == true){
		if(lowerSol == NULL){
		    lowerSol = new double[lN];
		    CoinFillN(lowerSol, lN, 0.0);
		}
	    std::copy(model_->seenLinkingSolutions[linkSol].lowerSolution.begin(),
		      model_->seenLinkingSolutions[linkSol].lowerSolution.end(), lowerSol);
	    }
	}

	for(i = 0; i < numScenarios; i++){
	    lowerObj = getLowerObj(sol, model_->getLowerObjSense(), i);

	    objVal = shouldStoreObjValues[i];

	    if(isIntegral_){
		assert((objVal - lowerObj) <= etol);
	    }
	
	    LPSolStatus_ = MibSLPSolStatusInfeasible;

	    //step 15
	    /** Current solution is bilevel feasible **/
	    if((fabs(objVal - lowerObj) < etol) && (isIntegral_)){
		if(i == numScenarios - 1){
		    LPSolStatus_ = MibSLPSolStatusFeasible;
		    useBilevelBranching_ = false;
	            shouldPrune_ = true;
	            storeSol = MibSRelaxationSol;
		}
	    }
	    else{
		if(lowerSol != NULL){
		    memcpy(optLowerSolutionOrd_, lowerSol, sizeof(double) * lN);
		}
		break;
	    }
	}
    
	
	if(!shouldPrune_){	
	    //step 18
	    if((tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved) &&
	       (((branchPar == MibSBranchingStrategyLinking) &&
		 (isIntegral_) && (isLinkVarsFixed_)) ||
		((computeBestUBWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
		((computeBestUBWhenLVarsInt == PARAM_ON)) ||
		((computeBestUBWhenLVarsFixed == PARAM_ON) && (isLinkVarsFixed_)))){
		/*if(UBSolver_){
		    UBSolver_ = setUpUBModel(model_->getSolver(), shouldStoreObjValues,
					     false);
		}
		else{
		    UBSolver_ =setUpUBModel(model_->getSolver(), shouldStoreObjValues,
					    true);
		}

		UBSolver = UBSolver_;*/
	      double *valuesUB = new double[uN + lN];
	      objVal = 0.0;
	      bool isUBProvenOptimal(true);
	      //In the stochastic case, when all first-level variables are linking,
	      //we can decompose the giant problem UB to smaller MIPs
	      int numDecomposedProbs(1);
	      bool useUBDecompose(model_->MibSPar_->entry
				  (MibSParams::useUBDecompose));
	      if((numScenarios > 1) && (useUBDecompose == true) && (sizeFixedInd == uN)){
		numDecomposedProbs = numScenarios;
		const double *relaxSol = oSolver->getColSolution();
		const double * uObjCoeffs = oSolver->getObjCoefficients();
		double uObjSense = oSolver->getObjSense();
		CoinDisjointCopyN(relaxSol, uN, valuesUB);
		for(i = 0; i < uN; i++){
		    objVal += relaxSol[i] * uObjCoeffs[i] * uObjSense;
		}
	      }

	      for(i = 0; i < numDecomposedProbs; i++){
		  remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();
		  if(remainingTime <= etol){
		      shouldPrune_ = true;
		      storeSol = MibSNoSol;
		      goto TERM_CHECKBILEVELFEAS;
		  }
		  
		if(numDecomposedProbs == 1){
		  UBSolver = setUpUBModel(oSolver, shouldStoreObjValues, true);
		}
		else{
		    if(i == 0){
			CoinPackedMatrix coefMatrix = *model_->origConstCoefMatrix_;
			int truncNumCols = uN + truncLN;
			int uRowNum = model_->getOrigUpperRowNum();
			CoinPackedVector col1;
			CoinPackedVector col2;
			int *colInd = NULL;
			double *colElem = NULL;
			int colNumElem(0);
			truncMatrixG2 = new CoinPackedMatrix(true, 0, 0);
			truncMatrixG2->setDimensions(model_->getTruncLowerRowNum(), 0);
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
			int isA2Random(model_->MibSPar_->entry(MibSParams::isA2Random));
			if(isA2Random == PARAM_OFF){
			    multA2XOpt = new double[model_->getTruncLowerRowNum()];
			}
			else{
			    multA2XOpt = new double[model_->getLowerRowNum()];
			}
			const double *lpSol = oSolver->getColSolution();
			double *optUpperSol = new double[uN];
			CoinDisjointCopyN(lpSol, uN, optUpperSol);
			model_->getStocA2Matrix()->times(optUpperSol, multA2XOpt);
			delete [] optUpperSol;
		    }
		    UBSolver = setUpDecomposedUBModel(oSolver, shouldStoreObjValues, i,
						      truncMatrixG2, multA2XOpt);
		}
		
		if(0)
		    UBSolver->writeLp("UBSolver");

		remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();

		if(remainingTime <= etol){
		    shouldPrune_ = true;
		    storeSol = MibSNoSol;
		    goto TERM_CHECKBILEVELFEAS;
		}
		
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

		//step 19
		startTimeUB = model_->broker_->subTreeTimer().getTime(); 
		UBSolver->branchAndBound();
		model_->timerUB_ += model_->broker_->subTreeTimer().getTime() - startTimeUB;
		model_->counterUB_ ++;
		if(i == numDecomposedProbs - 1){
		  isUBSolved_ = true;
		}

		if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
		    if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface *>
						 (UBSolver)->getSymphonyEnvironment())){
			shouldPrune_ = true;
			storeSol = MibSNoSol;
			goto TERM_CHECKBILEVELFEAS;
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
			shouldPrune_ = true;
			storeSol = MibSNoSol;
			goto TERM_CHECKBILEVELFEAS;
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
		    objVal += UBSolver->getObjValue() * model_->solver()->getObjSense();
		}
		else{
		    isUBProvenOptimal = false;
		    //when all first-level variables are linking,
		    //UBSolver cannot be infeasible
		    if(numDecomposedProbs > 1){
			throw CoinError("When all first-level variables are linking, problem UB cannot be infeasible",
					"checkBilevelFeasiblity", "MibSBilevel");
		    }
		}
		if(UBSolver){
		    delete UBSolver;
		}
	      }

	      if(isUBProvenOptimal){
		    isProvenOptimal_ = true;
		    if(shouldStoreSolution == true){
			std::copy(valuesUB, valuesUB + uN + lN, shouldStoreValuesUBSol.begin());
		    }
		    for (i = 0; i < uN + lN; i++){
			pos = binarySearch(0, uN - 1, i, upperColInd);
			if (pos >= 0){
			    if ((oSolver->isInteger(i)) &&
				(((valuesUB[i] - floor(valuesUB[i])) < etol) ||
				 ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
				optUpperSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
			    }else{
				optUpperSolutionOrd_[pos] = (double) valuesUB[i];
			    }
			}else{
			    if(numScenarios == 1){
				pos = binarySearch(0, lN - 1, i, lowerColInd);
			    }
			    else{
				pos = i - uN;
			    }
			    if ((oSolver->isInteger(i)) &&
				(((valuesUB[i] - floor(valuesUB[i])) < etol) ||
				 ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
				optLowerSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
			    }else{
				optLowerSolutionOrd_[pos] = (double) valuesUB[i];
			    }
			}
		    }
		    shouldStoreObjValues.push_back(0);
		    shouldStoreObjValues.clear();
		    shouldStoreObjValues.push_back(objVal);
		    storeSol = MibSHeurSol;
		}else{
		    shouldStoreObjValues.push_back(0);
		    isProvenOptimal_ = false;
		}
		//step 22
		//Adding x_L to set E
		if(useLinkingSolutionPool){
		    addSolutionToSeenLinkingSolutionPool
			(MibSLinkingPoolTagUBIsSolved, shouldStoreValuesUBSol,
			 shouldStoreObjValues);
		}
		shouldStoreValuesUBSol.clear();
		shouldStoreObjValues.clear();
	    
		//step 23
		if(isLinkVarsFixed_){
		    useBilevelBranching_ = false;
		    //isProvenOptimal_ = false;
		    shouldPrune_ = true;
		}
		delete [] valuesUB;
		if(truncMatrixG2){
		    delete truncMatrixG2;
		}
		if(multA2XOpt){
		    delete [] multA2XOpt;
		}
	    }
	    else if ((tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved) ||
		     ((!useLinkingSolutionPool) && (isUBSolved_))){
		if(lowerSol != NULL){
		    for (i = 0; i < lN; i++){
			if(numScenarios == 1){
			    index = lowerColInd[i];
			}
			else{
			    index = uN + i;
			}
			if ((model_->solver()->isInteger(index)) &&
			    (((lowerSol[i] - floor(lowerSol[i])) < etol) ||
			     ((ceil(lowerSol[i]) - lowerSol[i]) < etol))){
			    optLowerSolutionOrd_[i] = (double) floor(lowerSol[i] + 0.5);
			}else{
			    optLowerSolutionOrd_[i] = (double) lowerSol[i];
			}
		    }
		    if(isUpperIntegral_){
			storeSol = MibSHeurSol;
		    }
		}
	    }
	}
    }

 TERM_CHECKBILEVELFEAS:
    
    delete [] lowerSol;

    return storeSol;
}
	    
//#############################################################################
MibSSolType
MibSBilevel::checkBilevelFeasiblityParallel(bool isRoot)
{
  
  OsiSolverInterface *oSolver = model_->getSolver();

  int clockType(model_->AlpsPar()->entry(AlpsParams::clockType));

  int maxActiveNodes(model_->MibSPar_->entry
		     (MibSParams::maxActiveNodes)); 

  int numScenarios(model_->getNumScenarios());

  int maxThreadsLL(model_->MibSPar_->entry
		   (MibSParams::maxThreadsLL));
  int whichCutsLL(model_->MibSPar_->entry
		  (MibSParams::whichCutsLL));
  std::string feasCheckSolver(model_->MibSPar_->entry
			      (MibSParams::feasCheckSolver));
  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
    (model_->MibSPar_->entry(MibSParams::branchStrategy));
  int computeBestUBWhenXVarsInt(model_->MibSPar_->entry
				(MibSParams::computeBestUBWhenXVarsInt));
  int computeBestUBWhenLVarsInt(model_->MibSPar_->entry
				(MibSParams::computeBestUBWhenLVarsInt));
  int computeBestUBWhenLVarsFixed(model_->MibSPar_->entry
				  (MibSParams::computeBestUBWhenLVarsFixed));
  int useLinkingSolutionPool(model_->MibSPar_->entry
			     (MibSParams::useLinkingSolutionPool));
  double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));
  double remainingTime(0.0), startTimeVF(0.0), startTimeUB(0.0);
  MibSSolType storeSol(MibSNoSol);
  int truncLN(model_->truncLowerDim_);;
  int lN(model_->lowerDim_); // lower-level dimension
  int uN(model_->upperDim_); // lower-level dimension
  int i(0), j(0), k(0), index(0), length(0), pos(0), begPos(0);
  int lpStat;
  int sizeFixedInd(model_->sizeFixedInd_);
  double etol(model_->etol_), objVal(0.0), lowerObj(0.0);
  int * fixedInd = model_->fixedInd_;
  int * lowerColInd = model_->getLowerColInd();
  int * upperColInd = model_->getUpperColInd();
  CoinPackedMatrix *truncMatrixG2 = NULL;
  double *multA2XOpt = NULL;
  //OsiSolverInterface *UBSolver = 0;
  if(clockType != AlpsClockTypeWallClock){
    throw CoinError("AlpsClockType should be wallclock in the parallel mode",
		    "checkBilevelFeasiblityParallel", "MibSBilevel");
  }

  //saharSto: Check these changes
  double *lowerSol = NULL;
  std::vector<double> shouldStoreValuesLowerSol(1);
  std::vector<double> shouldStoreValuesUBSol(1);
  std::vector<double> shouldStoreObjValues;


  const double * sol = model_->solver()->getColSolution();

  std::vector<double> linkSol;
  for(i = 0; i < uN; i++){
    index = upperColInd[i];
    if(fixedInd[index] == 1){
      linkSol.push_back(upperSolutionOrd_[i]);
    }
  }

  isProvenOptimal_ = true;

  if(!isContainedInLinkingPool_){
    std::vector<OsiSolverInterface *>lSolverVec(numScenarios);
    for(i = 0; i < numScenarios; i++){
      remainingTime = timeLimit - model_->broker_->subTreeTimer().getWallClock();
      if(remainingTime <= etol){
	shouldPrune_ = true;
	storeSol = MibSNoSol;
	goto TERM_CHECKBILEVELFEASPARAL;
      }
      lSolverVec[i] = setUpModel(oSolver, true, i);
    }

    for(i = 0; i < numScenarios; i++){
      remainingTime = timeLimit - model_->broker_->subTreeTimer().getWallClock();
      if(remainingTime <= etol){
	shouldPrune_ = true;
	storeSol = MibSNoSol;
	goto TERM_CHECKBILEVELFEASPARAL;
      }

      remainingTime = CoinMax(remainingTime, 0.00);

      if (feasCheckSolver == "Cbc"){
	dynamic_cast<OsiCbcSolverInterface *>
	  (lSolverVec[i])->getModelPtr()->messageHandler()->setLogLevel(0);
      }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	//dynamic_cast<OsiSymSolverInterface *>
	// (lSolver)->setSymParam("prep_level", -1);

	sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (lSolverVec[i])->getSymphonyEnvironment();
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
	lSolverVec[i]->setHintParam(OsiDoReducePrint);
	lSolverVec[i]->messageHandler()->setLogLevel(0);
	CPXENVptr cpxEnv =
	  dynamic_cast<OsiCpxSolverInterface*>(lSolverVec[i])->getEnvironmentPtr();
	assert(cpxEnv);
	CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
	CPXsetintparam(cpxEnv, CPX_PARAM_CLOCKTYPE, 2);
	CPXsetdblparam(cpxEnv, CPX_PARAM_TILIM, remainingTime);
#endif
      }
    }

    startTimeVF = model_->broker_->subTreeTimer().getWallClock();
#ifdef _OPENMPMIBS
    omp_set_num_threads(maxActiveNodes);
#pragma omp parallel for
#endif
    for(i = 0; i < numScenarios; i++){
      lSolverVec[i]->branchAndBound();
    }
    model_->timerVF_ += model_->broker_->subTreeTimer().getWallClock() - startTimeVF;
    model_->counterVF_ += numScenarios;
    isLowerSolved_ = true;

    for(i = 0; i < numScenarios; i++){
      if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface *>
				     (lSolverVec[i])->getSymphonyEnvironment())){
	  shouldPrune_ = true;
	  storeSol = MibSNoSol;
	  goto TERM_CHECKBILEVELFEASPARAL;
	}
#endif
      }
      else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	lpStat = CPXgetstat(dynamic_cast<OsiCpxSolverInterface*>
			    (lSolverVec[i])->getEnvironmentPtr(),
			    dynamic_cast<OsiCpxSolverInterface*>
			    (lSolverVec[i])->getLpPtr());
	if((lpStat == CPXMIP_TIME_LIM_FEAS) ||
	   (lpStat == CPXMIP_TIME_LIM_INFEAS)){
	  shouldPrune_ = true;
	  storeSol = MibSNoSol;
	  goto TERM_CHECKBILEVELFEASPARAL;
	}
#endif
      }

      if(!lSolverVec[i]->isProvenOptimal()){
	LPSolStatus_ = MibSLPSolStatusInfeasible;
	isProvenOptimal_ = false;
	if(useLinkingSolutionPool){
	  //step 10
	  //Adding x_L to set E
	  shouldStoreObjValues.push_back(0);
	  addSolutionToSeenLinkingSolutionPool
	    (MibSLinkingPoolTagLowerIsInfeasible, shouldStoreValuesLowerSol,
	     shouldStoreObjValues);
	  shouldStoreObjValues.clear();
	}
	if(isLinkVarsFixed_){
	  useBilevelBranching_ = false;
	  shouldPrune_ = true;
	}
	break;
      }
      else{
	//const double * sol = model_->solver()->getColSolution();
	objVal = lSolverVec[i]->getObjValue() * model_->getLowerObjSense();

	//objVal_ = objVal;

	if((i == 0) && (!objValVec_.empty())){
	  objValVec_.clear();
	}

	objValVec_.push_back(objVal);
	shouldStoreObjValues.push_back(objVal);

	const double * values = lSolverVec[i]->getColSolution();

	if(useLinkingSolutionPool){
	  begPos = i * truncLN;
	  if(lowerSol == NULL){
	    lowerSol = new double[lN];
	    CoinFillN(lowerSol, lN, 0.0);
	  }
	  CoinDisjointCopyN(values, truncLN, lowerSol + begPos);


	  //step 12
	  //Adding x_L to set E
	  if(i == numScenarios - 1){
	    addSolutionToSeenLinkingSolutionPool
	      (MibSLinkingPoolTagLowerIsFeasible, shouldStoreValuesLowerSol,
	       shouldStoreObjValues);
	    shouldStoreValuesLowerSol.clear();
	    shouldStoreObjValues.clear();
	  }
	}
	else{
	  if(lowerSol == NULL){
	    lowerSol = new double[lN];
	    CoinFillN(lowerSol, lN, 0.0);
	  }
	  //memcpy(lowerSol, values, sizeof(double) * lN);
	  begPos = i * truncLN;
	  CoinDisjointCopyN(values, truncLN, lowerSol + begPos);
	}

	if(i == numScenarios - 1){
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
    }
    for(i = 0; i < numScenarios; i++){
      if(lSolverVec[i]){
	delete lSolverVec[i];
      }
    }
  }

  //step 13
  if(((!useLinkingSolutionPool) && (isProvenOptimal_)) ||
     ((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible) ||
      (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){

    //double *lowerSol = new double[lN];
    //CoinFillN(lowerSol, lN, 0.0);

    if(useLinkingSolutionPool){
      //get optimal value  of (VF) from solution pool
      //model_->it = seenLinkingSolutions.find(linkSol);
      //objVal = model_->it->second.lowerObjVal1;
      shouldStoreObjValues = model_->seenLinkingSolutions[linkSol].lowerObjValue;
      //objVal = seenLinkingSolutions.find(linkSol).
      //objVal_ = objVal;
      //objVal_ = shouldStoreObjValues[0];
      if(!objValVec_.empty()){
	objValVec_.clear();
      }
      objValVec_ = model_->seenLinkingSolutions[linkSol].lowerObjValue;
    }

    for(i = 0; i < numScenarios; i++){
      lowerObj = getLowerObj(sol, model_->getLowerObjSense(), i);

       objVal = shouldStoreObjValues[i];

       if(isIntegral_){
	 assert((objVal - lowerObj) <= etol);
       }

       LPSolStatus_ = MibSLPSolStatusInfeasible;

       //step 15
       /** Current solution is bilevel feasible **/
       if((fabs(objVal - lowerObj) < etol) && (isIntegral_)){
	 if(i == numScenarios - 1){
	   LPSolStatus_ = MibSLPSolStatusFeasible;
	   useBilevelBranching_ = false;
	   shouldPrune_ = true;
	   storeSol = MibSRelaxationSol;
	 }
       }
       else{
	 if(lowerSol != NULL){
	   memcpy(optLowerSolutionOrd_, lowerSol, sizeof(double) * lN);
	 }
	 break;
       }
    }

    if(!shouldPrune_){
      //step 18
      if((tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved) &&
	 (((branchPar == MibSBranchingStrategyLinking) &&
	   (isIntegral_) && (isLinkVarsFixed_)) ||
	  ((computeBestUBWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
	  ((computeBestUBWhenLVarsInt == PARAM_ON)) ||
	  ((computeBestUBWhenLVarsFixed == PARAM_ON) && (isLinkVarsFixed_)))){
	double *valuesUB = new double[uN + lN];
	objVal = 0.0;
	bool isUBProvenOptimal(true);
	//In the stochastic case, when all first-level variables are linking,
	//we can decompose the giant problem UB to smaller MIPs
	int numDecomposedProbs(1);
	bool useUBDecompose(model_->MibSPar_->entry
			    (MibSParams::useUBDecompose));
	if((useUBDecompose == true) && (sizeFixedInd < uN)){
	  throw CoinError("UB decomposition can be used for stochastic problems in which all upper level variables are linking ",
			  "checkBilevelFeasiblityParallel", "MibsBilevel");
	}
	
	if((numScenarios > 1) && (useUBDecompose == true) && (sizeFixedInd == uN)){
	  numDecomposedProbs = numScenarios;
	  const double *relaxSol = oSolver->getColSolution();
	  const double * uObjCoeffs = oSolver->getObjCoefficients();
	  double uObjSense = oSolver->getObjSense();
	  CoinDisjointCopyN(relaxSol, uN, valuesUB);
	  for(i = 0; i < uN; i++){
	    objVal += relaxSol[i] * uObjCoeffs[i] * uObjSense;
	  }
	}
	std::vector<OsiSolverInterface *>UBSolverVec(numDecomposedProbs);  
	for(i = 0; i < numDecomposedProbs; i++){
	  remainingTime = timeLimit - model_->broker_->subTreeTimer().getWallClock();
	  if(remainingTime <= etol){
	    shouldPrune_ = true;
	    storeSol = MibSNoSol;
	    goto TERM_CHECKBILEVELFEASPARAL;
	  }

	  if(numDecomposedProbs == 1){
	    UBSolverVec[i] = setUpUBModel(oSolver, shouldStoreObjValues, true);
	  }
	  else{
	    if(i == 0){
	      CoinPackedMatrix coefMatrix = *model_->origConstCoefMatrix_;
	      int truncNumCols = uN + truncLN;
	      int uRowNum = model_->getOrigUpperRowNum();
	      CoinPackedVector col1;
	      CoinPackedVector col2;
	      int *colInd = NULL;
	      double *colElem = NULL;
	      int colNumElem(0);
	      truncMatrixG2 = new CoinPackedMatrix(true, 0, 0);
	      truncMatrixG2->setDimensions(model_->getTruncLowerRowNum(), 0);
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
	      int isA2Random(model_->MibSPar_->entry(MibSParams::isA2Random));
	      if(isA2Random == PARAM_OFF){
		multA2XOpt = new double[model_->getTruncLowerRowNum()];
	      }
	      else{
		multA2XOpt = new double[model_->getLowerRowNum()];
	      }
	      const double *lpSol = oSolver->getColSolution();
	      double *optUpperSol = new double[uN];
	      CoinDisjointCopyN(lpSol, uN, optUpperSol);
	      model_->getStocA2Matrix()->times(optUpperSol, multA2XOpt);
	      delete [] optUpperSol;
	    }
	    UBSolverVec[i] = setUpDecomposedUBModel(oSolver, shouldStoreObjValues, i,
					      truncMatrixG2, multA2XOpt);
	  }
	}

	for(i = 0; i < numDecomposedProbs; i++){   
	  remainingTime = timeLimit - model_->broker_->subTreeTimer().getWallClock();

	  if(remainingTime <= etol){
	    shouldPrune_ = true;
	    storeSol = MibSNoSol;
	    goto TERM_CHECKBILEVELFEASPARAL;
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

	//step 19
	startTimeUB = model_->broker_->subTreeTimer().getWallClock();
#ifdef _OPENMPMIBS
	omp_set_num_threads(maxActiveNodes); 
#pragma omp parallel for
#endif
	for(i = 0; i < numDecomposedProbs; i++){
	  UBSolverVec[i]->branchAndBound();
	}
	model_->timerUB_ += model_->broker_->subTreeTimer().getWallClock() - startTimeUB;
	model_->counterUB_ += numDecomposedProbs;
        isUBSolved_ = true;

	for(i = 0; i < numDecomposedProbs; i++){ 
	  if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	    if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface *>
					 (UBSolverVec[i])->getSymphonyEnvironment())){
	      shouldPrune_ = true;
	      storeSol = MibSNoSol;
	      goto TERM_CHECKBILEVELFEASPARAL;
	    }
#endif
	  }
	  else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    lpStat = CPXgetstat(dynamic_cast<OsiCpxSolverInterface*>
				(UBSolverVec[i])->getEnvironmentPtr(),
				dynamic_cast<OsiCpxSolverInterface *>
				(UBSolverVec[i])->getLpPtr());
		if((lpStat == CPXMIP_TIME_LIM_FEAS) ||
		   (lpStat == CPXMIP_TIME_LIM_INFEAS)){
		  shouldPrune_ = true;
		  storeSol = MibSNoSol;
		  goto TERM_CHECKBILEVELFEASPARAL;
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
	    objVal += UBSolverVec[i]->getObjValue() * model_->solver()->getObjSense();
	  }
	  else{
	    isUBProvenOptimal = false;
	    //when all first-level variables are linking,
	    //UBSolver cannot be infeasible
	    if(numDecomposedProbs > 1){
	      throw CoinError("When all first-level variables are linking, problem UB cannot be infeasible",
			      "checkBilevelFeasiblity", "MibSBilevel");
	    }
	  }
	  if(UBSolverVec[i]){
	    delete UBSolverVec[i];
	  }
	}

	if(isUBProvenOptimal){
	  isProvenOptimal_ = true;
	  for (i = 0; i < uN + lN; i++){
	    pos = binarySearch(0, uN - 1, i, upperColInd);
	    if (pos >= 0){
	      if ((oSolver->isInteger(i)) &&
		  (((valuesUB[i] - floor(valuesUB[i])) < etol) ||
		   ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
		optUpperSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
	      }else{
		optUpperSolutionOrd_[pos] = (double) valuesUB[i];
	      }
	    }else{
		pos = i - uN;
	      if ((oSolver->isInteger(i)) &&
		  (((valuesUB[i] - floor(valuesUB[i])) < etol) ||
		   ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
		optLowerSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
	      }else{
		optLowerSolutionOrd_[pos] = (double) valuesUB[i];
	      }
	    }
	  }
	  shouldStoreObjValues.push_back(0);
	  shouldStoreObjValues.clear();
	  shouldStoreObjValues.push_back(objVal);
	  storeSol = MibSHeurSol;
	}else{
	  shouldStoreObjValues.push_back(0);
	  isProvenOptimal_ = false;
	}
	//step 22
	//Adding x_L to set E
	if(useLinkingSolutionPool){
	  addSolutionToSeenLinkingSolutionPool
	    (MibSLinkingPoolTagUBIsSolved, shouldStoreValuesUBSol,
	     shouldStoreObjValues);
	}
	shouldStoreValuesUBSol.clear();
	shouldStoreObjValues.clear();

	//step 23
	if(isLinkVarsFixed_){
	  useBilevelBranching_ = false;
	  //isProvenOptimal_ = false;
	  shouldPrune_ = true;
	}
	delete [] valuesUB;
	if(truncMatrixG2){
	  delete truncMatrixG2;
	}
	if(multA2XOpt){
	  delete [] multA2XOpt;
	}
      }
      else if ((tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved) ||
	       ((!useLinkingSolutionPool) && (isUBSolved_))){
	if(lowerSol != NULL){
	  for (i = 0; i < lN; i++){
	    index = uN + i;
	    if ((model_->solver()->isInteger(index)) &&
		(((lowerSol[i] - floor(lowerSol[i])) < etol) ||
		 ((ceil(lowerSol[i]) - lowerSol[i]) < etol))){
	      optLowerSolutionOrd_[i] = (double) floor(lowerSol[i] + 0.5);
	    }else{
	      optLowerSolutionOrd_[i] = (double) lowerSol[i];
	    }
	  }
	  if(isUpperIntegral_){
	    storeSol = MibSHeurSol;
	  }
	}
      }
    }
  }

 TERM_CHECKBILEVELFEASPARAL:
  delete [] lowerSol;

  return storeSol;
}

//#############################################################################
void 
MibSBilevel::gutsOfDestructor()
{

  if(optUpperSolutionOrd_) delete [] optUpperSolutionOrd_;
  if(optLowerSolutionOrd_) delete [] optLowerSolutionOrd_;
  if(upperSolutionOrd_) delete [] upperSolutionOrd_;
  if(lowerSolutionOrd_) delete [] lowerSolutionOrd_;
  if(lSolver_) delete lSolver_;
  //delete heuristic_;
}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpDecomposedUBModel(OsiSolverInterface * oSolver,
				    std::vector<double> &objValuesVec,
				    int scenarioIndex,
				    CoinPackedMatrix *truncMatrixG2,
				    double *multA2XOpt)
{

    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);

    int isA2Random(model_->MibSPar_->entry(MibSParams::isA2Random));

    OsiSolverInterface * nSolver;

    const double *lpSol = oSolver->getColSolution();

    int i(0), j(0);
    int rowIndex(0), colIndex(0), tmpIndex(0);
    int numElements(0);
    int uCols(model_->getUpperDim());
    int origURows(model_->getOrigUpperRowNum());
    int lRows(model_->getTruncLowerRowNum());
    int numTruncRows(origURows + lRows);
    int numCols(model_->getTruncLowerDim());
    int numRows(lRows + 1);
    double lObjSense(model_->getLowerObjSense());
    double uObjSense(1);
    const double * uObjCoeffs(oSolver->getObjCoefficients());
    double * lObjCoeffs(model_->getLowerObjCoeffs());
    const double * origColLb(model_->getOrigColLb());
    const double * origColUb(model_->getOrigColUb());
    const double * origRowLb(model_->getOrigRowLb());
    const double * origRowUb(model_->getOrigRowUb());
    double * upComp = new double[lRows];
    if(isA2Random != PARAM_OFF){
	rowIndex = scenarioIndex * lRows;
    }
    CoinDisjointCopyN(multA2XOpt + rowIndex, lRows, upComp); 

    /*CoinPackedMatrix matrix = *model_->origConstCoefMatrix_;
    matrix.reverseOrdering();
    CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
    newMat->setDimensions(0, numCols);

    CoinShallowPackedVector origRow;
    CoinPackedVector row;
    CoinPackedVector rowG2;
    CoinPackedVector rowA2;

    rowIndex = origURows + scenarioIndex * lRows;
    for(i = 0; i < lRows; i++){
	origRow = matrix.getVector(rowIndex);
	row = origRow;
	numElements = row.getNumElements();
	const int *indices = row.getIndices();
	const double *elements = row.getElements();
	for(j = 0; j < numElements; j++){
	    tmpIndex = indices[j];
	    if(tmpIndex < uCols){
		upComp[i] += elements[j] * lpSol[tmpIndex];
	    }
	    else{
		colIndex = indices[j] - uCols - scenarioIndex * numCols;
		rowG2.insert(colIndex, elements[j]);
	    }
	}
	newMat->appendRow(rowG2);
	rowG2.clear();
	row.clear();
	rowIndex ++;
	}*/

    CoinPackedMatrix newMat;
    CoinPackedVector row;
    newMat.reverseOrderedCopyOf(*truncMatrixG2); 

    for(i = 0; i < numCols; i++){
	row.insert(i, lObjCoeffs[i] * lObjSense);
    }
    newMat.appendRow(row);
    row.clear();

    double * objCoeffs = new double[numCols];
    double * rowUb = new double[numRows];
    double * rowLb = new double[numRows];
    double * colUb = new double[numCols];
    double * colLb = new double[numCols];

    CoinDisjointCopyN(uObjCoeffs + uCols, numCols, objCoeffs);

    CoinDisjointCopyN(origColLb + uCols, numCols, colLb);
    CoinDisjointCopyN(origColUb + uCols, numCols, colUb);

    tmpIndex = origURows + scenarioIndex * lRows;
    CoinDisjointCopyN(origRowLb + tmpIndex, lRows, rowLb);
    CoinDisjointCopyN(origRowUb + tmpIndex, lRows, rowUb);
    rowUb[numRows - 1] = objValuesVec[scenarioIndex];
    rowLb[numRows - 1] = -1 * oSolver->getInfinity();

    /*double etol(model_->etol_);
    double optMultTol(10000000.0);
    double optObj = objValuesVec[scenarioIndex];
    optObj = trunc(optObj * trunc(optMultTol));
    if(optObj >= etol){
	optObj ++;
    }
    optObj = optObj/trunc(optMultTol);
    rowUb[numRows - 1] = optObj;*/

    int * integerVars = new int[numCols];
    int intCnt(0);
    colIndex = uCols + scenarioIndex * numCols;
    /** Fill in array of integer variables **/
    for(i = 0; i < numCols; i++){
	if(oSolver->isInteger(colIndex)){
	    integerVars[intCnt] = i;
	    intCnt ++;
	}
	colIndex ++;
    }

    if(feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
    }else if(feasCheckSolver == "SYMPHONY"){
	#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
	#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"setUpDecomposedUBModel", "MibsBilevel");
	#endif
    }else if(feasCheckSolver == "CPLEX"){
	#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
	#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"setUpDecomposedUBModel", "MibsBilevel");
	#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"setUpDecomposedUBModel", "MibsBilevel");
    }

    nSolver->loadProblem(newMat, colLb, colUb,
			 objCoeffs, rowLb, rowUb);

    for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
    }

    nSolver->setObjSense(uObjSense); //1 min; -1 max

    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    for(i = 0; i < lRows; i++){
	nSolver->setRowLower(i, floor(rowLb[i] - upComp[i] + 0.5));
	nSolver->setRowUpper(i, floor(rowUb[i] - upComp[i] + 0.5));
    }

    delete [] integerVars;
    delete [] rowUb;
    delete [] rowLb;
    delete [] colUb;
    delete [] colLb;
    delete [] objCoeffs;
    delete [] upComp;

    return nSolver;
}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpUBModel(OsiSolverInterface * oSolver,
			  std::vector<double> &objValuesVec,
			  bool newOsi, const double *lpSol)
{
    
    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);

    OsiSolverInterface * nSolver;

    if (!lpSol){
	lpSol = oSolver->getColSolution();
    }

    int numScenarios(model_->getNumScenarios());
    
    int * fixedInd = model_->fixedInd_;

    int i(0), j(0), index1(0), beg(0);
    double value(0.0);
    double minusInf(-1 * oSolver->getInfinity());

    int uCols(model_->getUpperDim());
    int uRows(model_->getOrigUpperRowNum());
    int lRows(model_->getLowerRowNum());
    int lCols(model_->getLowerDim());
    int truncLCols(model_->getTruncLowerDim());
    int rowNum(model_->getNumOrigCons() + numScenarios);
    int colNum(model_->getNumOrigVars());
    int truncColNum(uCols + truncLCols);
    int * uColIndices(model_->getUpperColInd());

    if(newOsi){
	double objSense(model_->getLowerObjSense());
	double uObjSense(1);
        int * lColIndices(model_->getLowerColInd());
	int * uRowIndices(model_->getOrigUpperRowInd());
	int * lRowIndices(model_->getLowerRowInd());
	const double * origColLb(model_->getOrigColLb());
	const double * origColUb(model_->getOrigColUb());
	const double * origRowLb(model_->getOrigRowLb());
	const double * origRowUb(model_->getOrigRowUb());
	double * lObjCoeffs(model_->getLowerObjCoeffs());
	const double * uObjCoeffs(oSolver->getObjCoefficients());
	int tmpRowNum(rowNum - numScenarios);

	CoinPackedMatrix matrix = *model_->origConstCoefMatrix_;
	//matrix.reverseOrdering();

	double * rowUb = new double[rowNum];
        double * rowLb = new double[rowNum];
        double * colUb = new double[colNum];
        double * colLb = new double[colNum];

        CoinZeroN(colUb, colNum);
        CoinZeroN(colLb, colNum);

	/** Set the row bounds **/
	memcpy(rowLb, origRowLb, sizeof(double) * tmpRowNum);
	memcpy(rowUb, origRowUb, sizeof(double) * tmpRowNum);

	//Set the col bounds
	if(numScenarios == 1){
	    memcpy(colLb, origColLb, sizeof(double) * colNum);
	    memcpy(colUb, origColUb, sizeof(double) * colNum);
	}
	else{
	    memcpy(colLb, origColLb, sizeof(double) * truncColNum);
	    memcpy(colUb, origColUb, sizeof(double) * truncColNum);
	    for(i = 1; i < numScenarios; i++){
		beg = uCols + i * truncLCols;
	        CoinDisjointCopyN(origColLb + uCols, truncLCols, colLb + beg);
	        CoinDisjointCopyN(origColUb + uCols, truncLCols, colUb + beg);
	    }
	}

	for(i = 0; i < uCols; i++){
	    index1 = uColIndices[i];
	    if(fixedInd[index1] == 1){
		colLb[index1] = floor(lpSol[index1] + 0.5);
	        colUb[index1] = colLb[index1];
	    }
	}
    
        if (feasCheckSolver == "Cbc"){
	    nSolver = new OsiCbcSolverInterface();
	}else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	    nSolver = new OsiSymSolverInterface();
#else
	    throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			    "setUpUBModel", "MibsBilevel");
#endif
	}else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    nSolver = new OsiCpxSolverInterface();
#else
	    throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			    "setUpUBModel", "MibsBilevel");
#endif
	}else{
	    throw CoinError("Unknown solver chosen",
			    "setUpUBModel", "MibsBilevel");
	}
    
        int * integerVars = new int[colNum];
        double * objCoeffs = new double[colNum];
        //double * newRow = new double[lCols];
        CoinFillN(integerVars, colNum, 0);
        CoinFillN(objCoeffs, colNum, 0.0);
        //CoinFillN(newRow, lCols, 0.0);

        int intCnt(0);

	/** Fill in array of integer variables **/
	for(i = 0; i < colNum; i++){
	    if(oSolver->isInteger(i)){
		integerVars[intCnt] = i;
	        intCnt ++;
	    }
	}

        CoinDisjointCopyN(uObjCoeffs, colNum, objCoeffs);
        //CoinDisjointCopyN(lObjCoeffs, lCols, newRow);

        //for(i = 0; i < lCols; i++){
	//  newRow[i] = newRow[i] * objSense;
	//}

        //rowUb[rowNum-1] = objValLL;
        //rowLb[rowNum-1] = -1 * (oSolver->getInfinity());

        //CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
        //newMat->setDimensions(0, colNum);
	CoinPackedMatrix newMat;
	newMat.reverseOrderedCopyOf(matrix);    

        //for(i = 0; i < rowNum - numScenarios; i++){
	//  newMat->appendRow(matrix.getVector(i));
	//}

        CoinPackedVector row;
	//for(i = 0; i < lCols; i++){
	//  index1 = lColIndices[i];
	//  row.insert(index1, newRow[i]);
	//}
        //newMat->appendRow(row);

	if(numScenarios == 1){
	    double * newRow = new double[lCols];
	    CoinFillN(newRow, lCols, 0.0);
	    CoinDisjointCopyN(lObjCoeffs, lCols, newRow);
	    for(i = 0; i < lCols; i++){
		newRow[i] = newRow[i] * objSense;
	    }
	    for(i = 0; i < lCols; i++){
		index1 = lColIndices[i];
		row.insert(index1, newRow[i]);
	    }
	    newMat.appendRow(row);
	    rowUb[rowNum-1] = objValuesVec[0];
	    rowLb[rowNum-1] = minusInf;
	    delete [] newRow;
	}
	else{
	    CoinFillN(rowLb + tmpRowNum, rowNum - tmpRowNum, minusInf);
	    for(i = 0; i < numScenarios; i++){
		index1 = uCols + i * truncLCols;
		for(j = 0; j < truncLCols; j++){
		    row.insert(index1, lObjCoeffs[j] * objSense);
		    index1 ++;
		}
		newMat.appendRow(row);
		row.clear();
		rowUb[tmpRowNum + i] = objValuesVec[i];
	    }
	    
	}

        nSolver->loadProblem(newMat, colLb, colUb,
			     objCoeffs, rowLb, rowUb);

        for(i = 0; i < intCnt; i++){
	    nSolver->setInteger(integerVars[i]);
	}

        nSolver->setObjSense(uObjSense); //1 min; -1 max

        nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);
    
        delete [] integerVars;
        delete [] rowUb;
        delete [] rowLb;
        delete [] colUb;
        delete [] colLb;
        delete [] objCoeffs;
        //delete newMat;
    }
    /*else{
	nSolver = UBSolver_;
	for(i = 0; i < numScenarios; i++){
	    index1 = rowNum - (numScenarios - i);
	    nSolver->setRowUpper(index1, objValuesVec[i]);
	}
	for(i = 0; i < uCols; i++){
	    index1 = uColIndices[i];
	    if(fixedInd[index1] == 1){
		value = floor(lpSol[index1] + 0.5);
		nSolver->setColLower(index1, value);
		nSolver->setColUpper(index1, value);
	    }
	}
    }*/

    return nSolver;

}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpModel(OsiSolverInterface * oSolver, bool newOsi,
			int scenarioIndex, const double *lpSol)
{

  /** Create lower-level model with fixed upper-level vars **/

  int numScenarios(model_->getNumScenarios());

  int probType =
    model_->MibSPar_->entry(MibSParams::bilevelProblemType);

  bool warmStartLL =
    model_->MibSPar_->entry(MibSParams::warmStartLL);

  bool doDualFixing =
    model_->MibSPar_->entry(MibSParams::doDualFixing);

  std::string feasCheckSolver =
    model_->MibSPar_->entry(MibSParams::feasCheckSolver);

  OsiSolverInterface * nSolver;

  bool signsChanged(false);
  double etol(model_->etol_);
  int i(0), j(0), index1(0), index2(0);
  double mult(0.0);
  int uCols(model_->getUpperDim());
  int origURows(model_->getOrigUpperRowNum());
  int lRows(0), lCols(0);
  int * uColIndices = model_->getUpperColInd();
  int * lColIndices = model_->getLowerColInd();
  int * lRowIndices = model_->getLowerRowInd();

  if(numScenarios == 1){
      lRows = model_->getLowerRowNum();
      lCols = model_->getLowerDim();
  }
  else{
      lRows = model_->getTruncLowerRowNum();
      lCols = model_->getTruncLowerDim();
  }

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
  const double * lColLbInLProb = model_->getLColLbInLProb();
  const double * lColUbInLProb = model_->getLColUbInLProb();
  const char *origRowSense = model_->getOrigRowSense();
  double * rowUb = new double[lRows];
  double * rowLb = new double[lRows];
  //CoinZeroN(rowUb, lRows);
  //CoinZeroN(rowLb, lRows);
  //CoinZeroN(colUb, lCols);
  //CoinZeroN(colLb, lCols);

  
  if (newOsi){
      /** Set the row bounds **/
      for(i = 0; i < lRows; i++){
	  if(numScenarios == 1){
	      index1 = lRowIndices[i];
	  }
	  else{
	      index1 = origURows + scenarioIndex * lRows + i;
	  }
	  rowLb[i] = floor(origRowLb[index1] + 0.5);
	  rowUb[i] = floor(origRowUb[index1] + 0.5);
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

     double * colUb = new double[lCols];
     double * colLb = new double[lCols];
     int * integerVars = new int[lCols];
     double * objCoeffs = new double[lCols];
     
     CoinFillN(integerVars, lCols, 0);
     //CoinZeroN(objCoeffs, lCols);

     // set the col bounds
     if(lColLbInLProb == NULL){
	 for(i = 0; i < lCols; i++){
	     colLb[i] = origColLb[lColIndices[i]];
	     //colUb[i] = origColUb[lColIndices[i]];
	 }
     }
     else{
	 memcpy(colLb, lColLbInLProb, sizeof(double) * lCols);
     }

     if(lColUbInLProb == NULL){
	 for(i = 0; i < lCols; i++){
	     colUb[i] = origColUb[lColIndices[i]];
	 }
     }
     else{
	 memcpy(colUb, lColUbInLProb, sizeof(double) * lCols);
     }
     
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
     /*for(i = 0; i < lCols; i++){
	 objCoeffs[i] = objCoeffs[i] * objSense;
	 }*/
     
     CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
     newMat->setDimensions(0, lCols);
     CoinShallowPackedVector origRow;
     CoinPackedVector row;
     CoinPackedVector rowG2;
     int numElements(0), pos(0);

     for(i = 0; i < lRows; i++){
	 index1 = lRowIndices[i];
	 origRow = matrix->getVector(index1);
	 row = origRow;
	 numElements = row.getNumElements();
	 const int *indices = row.getIndices();
	 const double *elements = row.getElements();
	 for(j = 0; j < numElements; j++){
	     pos = binarySearch(0, lCols -1, indices[j], lColIndices);
	     if(pos >= 0){
		 rowG2.insert(pos, elements[j]);
	     }
	 }
	 newMat->appendRow(rowG2);
	 rowG2.clear();
	 row.clear();
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

     /*#if 0
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
#endif*/
     delete [] colLb;
     delete [] colUb;
     delete [] integerVars;
     delete [] objCoeffs;
     delete newMat;

  }else{
     nSolver = lSolver_;
     
     /*signsChanged = false;
     if(feasCheckSolver == "SYMPHONY"){
	 for(i = 0; i < lRows; i++){
	     index1 = lRowIndices[i];
	     if(origRowSense[index1] == 'G'){
		     signsChanged = true;
		     break;
	     }
	 }
     }

     if(!signsChanged){
	 for(i = 0; i < lRows; i++){
	     if(numScenarios == 1){
		 index1 = lRowIndices[i];
	     }
	     else{
		 index1 = origURows + scenarioIndex * lRows + i;
	     }
	     rowLb[i] = origRowLb[index1];
	     rowUb[i] = origRowUb[index1];
	 }
     }
     else{
	 for(i = 0; i < lRows; i++){
	     index1 = lRowIndices[i];
	     if(numScenarios == 1){
		 index2 = index1;
	     }
	     else{
		 index2 = origURows + scenarioIndex * lRows + i;
	     }
	     if(origRowSense[index1] == 'G'){
		 rowLb[i] = -origRowUb[index2];
		 rowUb[i] = -origRowLb[index2];
	     }
	     else{
		 rowLb[i] = origRowLb[index2];
		 rowUb[i] = origRowUb[index2];
	     }
	 }
     }*/
     
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
	  mult = 1;
	  if(numScenarios == 1){
	      index1 = lRowIndices[i];
	      if((signsChanged) && (origRowSense[index1] == 'G')){
		  mult = -1;
	      }
	  }
	  else{
	      index1 = origURows + scenarioIndex * lRows + i;
	      if((signsChanged) && (origRowSense[origURows + i] == 'G')){
		  mult = -1;
	      }
	  }
	  
	  for(j = 0; j < uCols; j++){
	      index2 = uColIndices[j];
	      coeff = matrix->getCoefficient(index1, index2);
	      if (coeff != 0){
		  upComp[i] += mult * coeff * lpSol[index2];
	      }
	  }
      }
     
     /** Correct the row bounds to account for fixed upper-level vars **/
     
     for(i = 0; i < lRows; i++){
	 nSolver->setRowLower(i, floor(rowLb[i] - upComp[i] + 0.5));
	 nSolver->setRowUpper(i, floor(rowUb[i] - upComp[i] + 0.5));
     }
     
     delete [] upComp;

#if SYMPHONY_VERSION_IS_WS
  }
#endif

  //I don't think this is needed
  //if(!getWarmStart())
  //  setWarmStart(nSolver->getWarmStart());

  delete [] rowLb;
  delete [] rowUb;
  
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
MibSBilevel::getLowerObj(const double * sol, double objSense, double scenarioIndex)
{

   int numScenarios(model_->getNumScenarios());
   int truncLCols(model_->getTruncLowerDim());
   int lCols(model_->getLowerDim());
   int * lColIndices = model_->getLowerColInd();
   double * lObjCoeffs = model_->getLowerObjCoeffs();

   int i(0), index(0);
   double objVal(0.0);


   if(numScenarios == 1){
       for(i = 0; i < lCols; i++){
       index = lColIndices[i];
       if(0){
	   std::cout << "sol[" << index << "]: " << sol[index] << std::endl;
	   std::cout << "lObjCoeff[" << i << "]: " << lObjCoeffs[i] << std::endl;
       }
       objVal += lObjCoeffs[i] * sol[index];
       }
   }
   else{
       for(i = 0; i < truncLCols; i++){
	   index = scenarioIndex * truncLCols + lColIndices[i];
	   objVal += lObjCoeffs[i] * sol[index];
       }
   }

   return objVal * objSense;

}
//#############################################################################
void
    MibSBilevel::addSolutionToSeenLinkingSolutionPool(MibSLinkingPoolTag solTag, std::vector<double>
						      &shouldStoreValues,
						      std::vector<double> &shouldStoreObjValues)
{
    
    int i(0),index(0);
    int uN(model_->upperDim_);
    int lN(model_->lowerDim_);
    int * upperColInd = model_->getUpperColInd();
    int * fixedInd = model_->fixedInd_;
    int solType = static_cast<int>(solTag);

    LINKING_SOLUTION linkingSolution;
    
    std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd_[i]);
	}
    }

    tagInSeenLinkingPool_ = solTag;
    linkingSolution.lowerObjValue.push_back(0);
    linkingSolution.lowerObjValue.clear();
    linkingSolution.lowerSolution.push_back(0);
    linkingSolution.UBSolution.push_back(0);
    linkingSolution.lowerSolution.clear();
    linkingSolution.UBSolution.clear();
    
    switch(solTag){
    case MibSLinkingPoolTagLowerIsInfeasible:
	{
	    linkingSolution.tag = solType;
	    linkingSolution.lowerObjValue.push_back(0);
	    linkingSolution.UBObjValue = 0.0;
	    linkingSolution.lowerSolution.push_back(0);
	    linkingSolution.UBSolution.push_back(0);
	    model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    break;
	}
    case MibSLinkingPoolTagLowerIsFeasible:
	{
	    linkingSolution.tag = solType;
	    linkingSolution.lowerObjValue = shouldStoreObjValues;
	    linkingSolution.UBObjValue = 0.0;
	    linkingSolution.lowerSolution = shouldStoreValues;
	    linkingSolution.UBSolution.push_back(0);
	    model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    break;
	}
    case MibSLinkingPoolTagUBIsSolved:
	{
	    model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagUBIsSolved;
	    model_->seenLinkingSolutions[linkSol].UBObjValue = shouldStoreObjValues[0];
	    if(isProvenOptimal_){
		model_->seenLinkingSolutions[linkSol].UBSolution.clear();
		model_->seenLinkingSolutions[linkSol].UBSolution = shouldStoreValues;
	    }
	    break;
	}
    }
}
	    
	
	
	    






    
 
