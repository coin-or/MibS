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
MibSSolType
MibSBilevel::createBilevel(CoinPackedVector* sol, 
			   MibSModel *mibs)
{
  
  /** Splits sol into two parts, upper- and lower-level components **/

  if(!mibs) return MibSNoSol;
  
  model_ = mibs;

  heuristic_ = new MibSHeuristic(mibs);

  int i(0),j(0);
  int N(model_->numVars_);
  int lN(model_->lowerDim_); // lower-level dimension
  int uN(model_->upperDim_); // upper-level dimension
  double etol(model_->BlisPar()->entry(BlisParams::integerTol));
  
  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
      (model_->MibSPar_->entry(MibSParams::branchStrategy));

  bool solveSecondLevelEveryIteration(model_->MibSPar_->entry
                                     (MibSParams::solveSecondLevelEveryIteration) == PARAM_ON);
  bool solveSecondLevelWhenXYVarsInt(model_->MibSPar_->entry
                                     (MibSParams::solveSecondLevelWhenXYVarsInt) == PARAM_ON);
  bool solveSecondLevelWhenXVarsInt(model_->MibSPar_->entry
                                    (MibSParams::solveSecondLevelWhenXVarsInt) == PARAM_ON);
  bool solveSecondLevelWhenLVarsInt(model_->MibSPar_->entry
                                    (MibSParams::solveSecondLevelWhenLVarsInt) == PARAM_ON);
  bool solveSecondLevelWhenLVarsFixed(model_->MibSPar_->entry
                                      (MibSParams::solveSecondLevelWhenLVarsFixed) == PARAM_ON);
  int cutStrategy(model_->MibSPar_->entry
		  (MibSParams::cutStrategy));

  int useLinkingSolutionPool(model_->MibSPar_->entry
			   (MibSParams::useLinkingSolutionPool));

  MibSSolType storeSol(MibSNoSol);
  
  assert(N == model_->solver()->getNumCols());

  model_->setNumRows(model_->solver()->getNumRows());
  model_->setUpperRowNum(model_->solver()->getNumRows() - model_->getLowerRowNum());
  model_->setUpperRowData();
  
  int *indices = sol->getIndices();
  double *values = sol->getElements();
  int numElements(sol->getNumElements()); // number of nonzero elements
  int * fixedInd = model_->fixedInd_; 
 
  if(!upperSolutionOrd_){
      upperSolutionOrd_ = new double[uN];
  }
  if(!lowerSolutionOrd_){
      lowerSolutionOrd_ = new double[lN];
  }
  if(!optUpperSolutionOrd_){
      optUpperSolutionOrd_ = new double[uN];
  }
  if(!optLowerSolutionOrd_){
      optLowerSolutionOrd_ = new double[lN];
  }
  
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
  if((isLinkVarsFixed_) && 
    ((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsInfeasible) ||
      (tagInSeenLinkingPool_ == MibSLinkingPoolTagPesIsInfeasible) || // YX: pessimistic case 
      (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){
      useBilevelBranching_ = false;
      shouldPrune_ = true;
  }

  //step 7
  if(!shouldPrune_){
     if((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible ||
         tagInSeenLinkingPool_ == MibSLinkingPoolTagPesIsFeasible || // YX: pessimistic case
          tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved) ||
         (!isContainedInLinkingPool_ &&
	  ((branchPar == MibSBranchingStrategyLinking && isIntegral_ && isLinkVarsFixed_) ||
	   (branchPar == MibSBranchingStrategyFractional && isIntegral_) ||
	   (solveSecondLevelEveryIteration) ||
	   (solveSecondLevelWhenXYVarsInt && isIntegral_) ||
	   (solveSecondLevelWhenXVarsInt && isUpperIntegral_) ||
	   (solveSecondLevelWhenLVarsInt && isLinkVarsIntegral_) ||
	   (solveSecondLevelWhenLVarsFixed && isLinkVarsFixed_ )))){
	  storeSol = checkBilevelFeasibility(mibs->isRoot_);
      }
  }
  
  if(cutStrategy == 1){
      useBilevelBranching_ = false;
  }
  
  heuristic_->findHeuristicSolutions();

  delete heuristic_;

  return storeSol;
}


//#############################################################################
MibSSolType
MibSBilevel::checkBilevelFeasibility(bool isRoot)
{
    bool warmStartLL(model_->MibSPar_->entry
		     (MibSParams::warmStartLL));
    bool findPesSol(model_->MibSPar_->entry
		     (MibSParams::findPesSol)); // YX: pessimistic case
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

	// YX: target optimality gap for bounded rationality			
	double targetGap(model_->MibSPar()->entry(MibSParams::slTargetGap));
    double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));
    double remainingTime(0.0), startTimeVF(0.0), startTimePES(0.0), startTimeUB(0.0); // YX: (PES-MILP)
    MibSSolType storeSol(MibSNoSol);
    int lN(model_->lowerDim_); // lower-level dimension
    int uN(model_->upperDim_); // upper-level dimension
    int uRows(model_->getOrigUpperRowNum());
    int i(0), index(0), length(0), pos(0);
    int sizeFixedInd(model_->sizeFixedInd_);
    double etol(model_->etol_), objVal(0.0), lowerObj(0.0);
    double pesRF(0.0), lpPesVal(0.0); // YX: pes rf values
    int * fixedInd = model_->fixedInd_;
    int * lowerColInd = model_->getLowerColInd();
    int * upperColInd = model_->getUpperColInd();
    const double * valuesPes; // YX: (PES-MILP) solution (x,y) ptr 
    double *lowerSol = new double[lN];
    double *pesSol = new double[lN + uN]; // YX: processed pes solution
    CoinFillN(lowerSol, lN, 0.0);
    CoinFillN(pesSol, lN + uN, 0.0);

    std::vector<double> shouldStoreValuesLowerSol(lN);
    std::vector<double> shouldStoreValuesPesSol(lN + uN); // YX: pessimistic case
    std::vector<double> shouldStoreValuesUBSol(lN + uN);

    const double * sol = model_->solver()->getColSolution();
    double gap = (targetGap < etol) ? 0.0 : targetGap; // YX: added SL gap

    std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd_[i]);
	}
    }

    isProvenOptimal_ = true; 

    if(!isContainedInLinkingPool_){
	
	//isProvenOptimal_ = true;
    
	/*if (warmStartLL && (feasCheckSolver == "SYMPHONY") && solver_){
	    solver_ = setUpModel(model_->getSolver(), false);
	}else{
	    //if (solver_){
	//	delete solver_;
	   // }
	    solver_ = setUpModel(model_->getSolver(), true);
	}*/

	if(lSolver_){
	    lSolver_ = setUpModel(model_->getSolver(), false);
	}
	else{
	    lSolver_ = setUpModel(model_->getSolver(), true);
	}
	    
	OsiSolverInterface *lSolver = lSolver_;

    lSolver->writeLp("lSolver");
    model_->getSolver()->writeLp("oSolver");

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
	  startTimeVF = model_->broker_->subTreeTimer().getTime();
	  lSolver->branchAndBound();
	  model_->timerVF_ += model_->broker_->subTreeTimer().getTime() - startTimeVF;
	}
  
	model_->counterVF_ ++;
	isLowerSolved_ = true;
	
	if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					       (dynamic_cast<OsiSymSolverInterface *>
						(lSolver)->getSymphonyEnvironment()))){
	    shouldPrune_ = true;
	    storeSol = MibSNoSol;
	    goto TERM_CHECKBILEVELFEAS;
	}
	else if(!lSolver->isProvenOptimal()){
	    LPSolStatus_ = MibSLPSolStatusInfeasible;
	    isProvenOptimal_ = false;
	    if(useLinkingSolutionPool && isLinkVarsIntegral_){
	    //step 10
	    //Adding x_L to set E
		addSolutionToSeenLinkingSolutionPool
		    (MibSLinkingPoolTagLowerIsInfeasible, shouldStoreValuesLowerSol, 0.0);
	    }
	    if(isLinkVarsFixed_){
	      useBilevelBranching_ = false;
	      shouldPrune_ = true;
	    }
	}
	else{
	    //const double * sol = model_->solver()->getColSolution();
	    objVal = lSolver->getObjValue() * model_->getLowerObjSense();

	    objVal_ = objVal;

	    const double * values = lSolver->getColSolution();
	    
        // YX: round lSolver solution; sol (Relaxation soln) is not rounded?
        for (i = 0; i < lN; i++){			 
            if ((lSolver->isInteger(i)) &&
                (((values[i] - floor(values[i])) < etol) ||
                ((ceil(values[i]) - values[i]) < etol))){
                lowerSol[i] = (double) floor(values[i] + 0.5);
            }else{
                lowerSol[i] = (double) values[i];
            }
        }

	    if(useLinkingSolutionPool && isLinkVarsIntegral_){
		std::copy(lowerSol, lowerSol + lN, shouldStoreValuesLowerSol.begin());
		
		//step 12
		//Adding x_L to set E  
		addSolutionToSeenLinkingSolutionPool
		    (MibSLinkingPoolTagLowerIsFeasible, shouldStoreValuesLowerSol, objVal_);
		shouldStoreValuesLowerSol.clear();
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
        
        // YX: pessimistic case; enter only when isUpperIntegral_?
        if(findPesSol){
            if(pSolver_){
                pSolver_ = setUpPesModel(objVal, false);
            }else{
                pSolver_ = setUpPesModel(objVal, true);
            }
            
            OsiSolverInterface *pSolver = pSolver_;

            pSolver->writeLp("pSolverLoaded"); // YX: debug only

            // YX: add timer here
            remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();

            if(remainingTime <= etol){
                shouldPrune_ = true;
                goto TERM_CHECKBILEVELFEAS;
            }
            
            remainingTime = CoinMax(remainingTime, 0.00);

            if (feasCheckSolver == "Cbc"){
                dynamic_cast<OsiCbcSolverInterface *>
                (pSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
            }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
                //dynamic_cast<OsiSymSolverInterface *>
                // (pSolver)->setSymParam("prep_level", -1);
                sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
                (pSolver)->getSymphonyEnvironment();
                //Always uncomment for debugging!!
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
                pSolver->setHintParam(OsiDoReducePrint);
                pSolver->messageHandler()->setLogLevel(0);
                CPXENVptr cpxEnv =
                dynamic_cast<OsiCpxSolverInterface*>(pSolver)->getEnvironmentPtr();
                assert(cpxEnv);
                CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
                CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
            }

            startTimePES = model_->broker_->subTreeTimer().getTime(); 
            pSolver->branchAndBound();
            model_->timerPES_ += model_->broker_->subTreeTimer().getTime() - startTimePES;
            // YX: add tracker later? solved times to match SL-MILP counter?
            
            if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
                                (dynamic_cast<OsiSymSolverInterface *>
                                (pSolver)->getSymphonyEnvironment()))){
                shouldPrune_ = true;
                goto TERM_CHECKBILEVELFEAS;
            }

            if(!pSolver->isProvenOptimal()){
                isProvenOptimal_ = false;
                if(useLinkingSolutionPool && isLinkVarsIntegral_){
                    addSolutionToSeenLinkingSolutionPool
                        (MibSLinkingPoolTagPesIsInfeasible, shouldStoreValuesPesSol, 0.0);
                }
                if(isLinkVarsFixed_){
                    useBilevelBranching_ = false;
                    shouldPrune_ = true;
                }
            }else{
                pesRF = pSolver->getObjValue(); // YX: obtain max d^1y
                valuesPes = pSolver->getColSolution();
                
                // YX: sort pes solutions here with position shift;
                //  Then saved linking pool pes solutions are sorted and rounded 
                //  (！！which is different from UB solutions.)
                for(i = 0; i < uN + lN; i++){
                    pos = binarySearch(0, uN - 1, i, upperColInd);
                    if(pos >= 0){
                        if((pSolver->isInteger(i)) &&
                            (((valuesPes[i] - floor(valuesPes[i])) < etol) ||
                            ((ceil(valuesPes[i]) - valuesPes[i]) < etol))){
                            pesSol[pos] = (double) floor(valuesPes[i] + 0.5);
                        }else{
                            pesSol[pos] = (double) valuesPes[i];
                        }
                    }else{
                        pos = binarySearch(0, lN - 1, i, lowerColInd);
                        if((pSolver->isInteger(i)) &&
                            (((valuesPes[i] - floor(valuesPes[i])) < etol) ||
                            ((ceil(valuesPes[i]) - valuesPes[i]) < etol))){
                            pesSol[pos+uN] = (double) floor(valuesPes[i] + 0.5);
                        }else{
                            pesSol[pos+uN] = (double) valuesPes[i];
                        }
                    }
                }

                if(useLinkingSolutionPool && isLinkVarsIntegral_){
                    std::copy(pesSol, pesSol + uN + lN, shouldStoreValuesPesSol.begin());
                    addSolutionToSeenLinkingSolutionPool
                        (MibSLinkingPoolTagPesIsFeasible, shouldStoreValuesPesSol, pesRF);
                    shouldStoreValuesPesSol.clear();
                }
            }
        }
	}
	/*if (!warmStartLL){
	    delete solver_;
	}*/
    } // end of (solving SL-MILP and PES-MILP if not solved before)

    //step 13
	if(((!useLinkingSolutionPool || !isLinkVarsIntegral_) && (isProvenOptimal_)) ||
        ((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible) ||
        (tagInSeenLinkingPool_ == MibSLinkingPoolTagPesIsFeasible) ||
        (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){

	if(useLinkingSolutionPool && isLinkVarsIntegral_){
	    //get optimal value of (VF) from solution pool
	    objVal = model_->seenLinkingSolutions[linkSol].lowerObjValue; 
	    objVal_ = objVal;
	    std::copy(model_->seenLinkingSolutions[linkSol].lowerSolution.begin(),
		      model_->seenLinkingSolutions[linkSol].lowerSolution.end(), lowerSol);
        
        if(findPesSol){ // YX: pessimistic case; retrieve max d^1y
            pesRF = model_->seenLinkingSolutions[linkSol].pesRFValue;
            std::copy(model_->seenLinkingSolutions[linkSol].pesSolution.begin(),
              model_->seenLinkingSolutions[linkSol].pesSolution.end(), pesSol);
        }  
	}
	lowerObj = getLowerObj(sol, model_->getLowerObjSense());
	
	if(isIntegral_){
           assert((objVal - lowerObj) <= etol);
	}
	
	LPSolStatus_ = MibSLPSolStatusInfeasible;
	
	//step 15
    // YX: [d2^y - phi(A^x)]/abs(phi(A^x)), combined check when gap is set to 0;
    if(!findPesSol){
        if(((lowerObj - objVal) <= fabs(objVal) * gap/100) && (isIntegral_)){
        // if(((lowerObj - objVal) < etol) && (isIntegral_)){
            LPSolStatus_ = MibSLPSolStatusFeasible;
            useBilevelBranching_ = false;
            shouldPrune_ = true;
            storeSol = MibSRelaxationSol;            
        }else{
            memcpy(optLowerSolutionOrd_, lowerSol, sizeof(double) * lN);
            if(isUpperIntegral_){
                storeSol = MibSHeurSol;
            }
        }
    }else{
        lpPesVal = getRiskFuncVal(lowerSolutionOrd_);           
        // YX: check whether (LR) solution is a pessimistic solution
        if(((lowerObj - objVal) <= fabs(objVal) * gap/100) && (isIntegral_) &&
            (fabs(lpPesVal - pesRF) < etol)){
            useBilevelBranching_ = false;
            shouldPrune_ = true;
            storeSol = MibSRelaxationSol;
        }else{
            // YX: copy (PES-MILP) solution to optimal solution container
            memcpy(optUpperSolutionOrd_, pesSol, sizeof(double) * uN);
            memcpy(optLowerSolutionOrd_, pesSol + uN, sizeof(double) * lN); 
            storeSol = MibSHeurSol;
        }
    }		
	
	if(!shouldPrune_){	
	    //step 18
	    if((tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved) &&
	       (((branchPar == MibSBranchingStrategyLinking) &&
		 (isIntegral_) && (isLinkVarsFixed_)) ||
		((computeBestUBWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
		((computeBestUBWhenLVarsInt == PARAM_ON)  && (isLinkVarsIntegral_)) ||
		((computeBestUBWhenLVarsFixed == PARAM_ON) && (isLinkVarsFixed_)))){
        
        // YX: assume (PES-MILP) is solved and a solution is found;
        if((findPesSol)&&(((sizeFixedInd - uN) >= 0)||(uRows == 0))){
            // YX: all variables are linking variables; solution fixed
            isUBSolved_ = true; 
            
            // YX: update Upper bound value; UBSolution default is 0;
            objVal = getUpperObj(optLowerSolutionOrd_, optUpperSolutionOrd_);
            
        }else{
            if(UBSolver_){
                UBSolver_ = setUpUBModel(model_->getSolver(), objVal, false);
            }else{
                UBSolver_ = setUpUBModel(model_->getSolver(), objVal, true);
            }

        OsiSolverInterface * UBSolver = UBSolver_;

		remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();

		if(remainingTime <= etol){
		    shouldPrune_ = true;
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
		    UBSolver->setHintParam(OsiDoReducePrint);
		    UBSolver->messageHandler()->setLogLevel(0);
		    CPXENVptr cpxEnv =
			dynamic_cast<OsiCpxSolverInterface*>(UBSolver)->getEnvironmentPtr();
		    assert(cpxEnv);
		    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
		    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
		}

		//step 19
		startTimeUB = model_->broker_->subTreeTimer().getTime(); 
		UBSolver->branchAndBound();
		model_->timerUB_ += model_->broker_->subTreeTimer().getTime() - startTimeUB;
		model_->counterUB_ ++;
		
		if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
						       (dynamic_cast<OsiSymSolverInterface *>
							(UBSolver)->getSymphonyEnvironment()))){
		    shouldPrune_ = true;
		    goto TERM_CHECKBILEVELFEAS;
		}

		isUBSolved_ = true;

        // YX: pess UB has to have an optimal solution in case of feasible (PES-MILP)
        if(findPesSol){
            assert(UBSolver->isProvenOptimal());
        }

        if (UBSolver->isProvenOptimal()){
		    isProvenOptimal_ = true;
		    const double * valuesUB = UBSolver->getColSolution();
		    std::copy(valuesUB, valuesUB + uN + lN, shouldStoreValuesUBSol.begin());
		    for (i = 0; i < uN + lN; i++){
			pos = binarySearch(0, uN - 1, i, upperColInd);
			if (pos >= 0){
			    if ((UBSolver->isInteger(i)) &&
				(((valuesUB[i] - floor(valuesUB[i])) < etol) ||
				 ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
				optUpperSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
			    }else{
				optUpperSolutionOrd_[pos] = (double) valuesUB[i];
			    }
			}else{
			    pos = binarySearch(0, lN - 1, i, lowerColInd);
			    if ((UBSolver->isInteger(i)) &&
				(((valuesUB[i] - floor(valuesUB[i])) < etol) ||
				 ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
				optLowerSolutionOrd_[pos] = (double) floor(valuesUB[i] + 0.5);
			    }else{
				optLowerSolutionOrd_[pos] = (double) valuesUB[i];
			    }
			}
		    }
            if(!findPesSol){
		        objVal = UBSolver->getObjValue() * model_->solver()->getObjSense();
            }else{
                objVal = getUpperObj(optLowerSolutionOrd_, optUpperSolutionOrd_);
            }
		    storeSol = MibSHeurSol;
		}else{
		    isProvenOptimal_ = false;
		    storeSol = MibSNoSol;
		}
        
        } // end of (solving an UB problem in either optimistic or pessimistic case)
		
        //step 22
		//Adding x_L to set E
		if(useLinkingSolutionPool && isLinkVarsIntegral_){
		    addSolutionToSeenLinkingSolutionPool
			(MibSLinkingPoolTagUBIsSolved, shouldStoreValuesUBSol, objVal);
		    shouldStoreValuesUBSol.clear();
		}

		//step 23
		if(isLinkVarsFixed_){
		    useBilevelBranching_ = false;
		    //isProvenOptimal_ = false;
		    shouldPrune_ = true;
		}	
	    }
	} // end of (if not prune)
    } // end of (bilevel feasibility check and UB finding)

 TERM_CHECKBILEVELFEAS:
    
    delete [] lowerSol;
    delete [] pesSol;

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
  if(pSolver_) delete pSolver_; // YX: pessimistic case
  if(UBSolver_) delete UBSolver_;
  //delete heuristic_;
}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpUBModel(OsiSolverInterface * oSolver, double objValLL,
			      bool newOsi, const double *lpSol)
{
    
    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);
    double targetGap(model_->MibSPar_->entry(MibSParams::slTargetGap));
    bool findPes(model_->MibSPar_->entry(MibSParams::findPesSol)); // YX: add pes option

    OsiSolverInterface * nSolver;

    if (!lpSol){
	lpSol = oSolver->getColSolution();
    }

    int * fixedInd = model_->fixedInd_;

    int i(0), j(0), index1(0);
    double value(0.0), objUb(0.0);

    int uCols(model_->getUpperDim());
    int uRows(model_->getOrigUpperRowNum());
    int lRows(model_->getLowerRowNum());
    int lCols(model_->getLowerDim());
    int rowNum = findPes ? model_->getNumOrigCons() + 2: model_->getNumOrigCons() + 1;
    int colNum(model_->getNumOrigVars());
    int * uColIndices(model_->getUpperColInd());
    double gap = (targetGap < model_->getTolerance()) ? 0.0 : targetGap; // YX: added SL gap

    // YX: d2^y <= phi(A^x) + \delta * abs(phi(A^x))
    objUb = objValLL + fabs(objValLL) * gap/100;

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
    int tmpRowNum = findPes ? rowNum - 2: rowNum - 1;

	CoinPackedMatrix matrix = *model_->origConstCoefMatrix_;
	matrix.reverseOrdering();

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
	memcpy(colLb, origColLb, sizeof(double) * colNum);
	memcpy(colUb, origColUb, sizeof(double) * colNum);

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

        CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
        newMat->setDimensions(0, colNum);

        for(i = 0; i < rowNum - 1; i++){
	    newMat->appendRow(matrix.getVector(i));
	}
        /** Add pessimistic risk function constraint at (rowNum-2) **/
        if(findPes){
            CoinPackedVector row2;
            for(i = 0; i < lCols; i++){
                index1 = lColIndices[i];
                row2.insert(index1, uObjCoeffs[index1] * uObjSense);
                value += uObjCoeffs[index1] * optLowerSolutionOrd_[i];
                // YX: set up obj in pessimistic case cx+0y
                objCoeffs[index1] = 0.0; // !! verify obj index
            }
            newMat->appendRow(row2);
            rowUb[rowNum-2] = oSolver->getInfinity();
            rowLb[rowNum-2] = value * uObjSense;
        }

        /** Add value function constraint at (rowNum-1) **/
        CoinDisjointCopyN(lObjCoeffs, lCols, newRow);

        for(i = 0; i < lCols; i++){
            newRow[i] = newRow[i] * objSense;
        }

        rowUb[rowNum-1] = objUb; // YX: gap added to phi(A^2x)
        rowLb[rowNum-1] = -1 * (oSolver->getInfinity());
        
        CoinPackedVector row;
	for(i = 0; i < lCols; i++){
	    index1 = lColIndices[i];
	    row.insert(index1, newRow[i]);
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
        delete [] rowUb;
        delete [] rowLb;
        delete [] colUb;
        delete [] colLb;
        delete [] objCoeffs;
        delete [] newRow;
        delete newMat;
    }
    else{
	nSolver = UBSolver_;
	nSolver->setRowUpper(rowNum-1, objUb); // YX: gap added to phi(A^2x)
    if(findPes){
        nSolver->setRowLower(rowNum-2, getRiskFuncVal(optLowerSolutionOrd_));
    }
	for(i = 0; i < uCols; i++){
	    index1 = uColIndices[i];
	    if(fixedInd[index1] == 1){
		value = floor(lpSol[index1] + 0.5);
		nSolver->setColLower(index1, value);
		nSolver->setColUpper(index1, value);
	    }
	}
    }

    return nSolver;

}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpPesModel(double objValLL, bool newOsi, const double *lpSol)
{
   /** Setup pessimistic risk function model (PES-MILP) with fixed linking part **/

   std::string feasCheckSolver(model_->MibSPar()->entry
      (MibSParams::feasCheckSolver));
   double targetGap(model_->MibSPar()->entry(MibSParams::slTargetGap));

   OsiSolverInterface * nSolver;

   int i(0), idx(0);
   double value(0.0), objUb(0.0);

   int uCols(model_->getUpperDim());
   int lCols(model_->getLowerDim());
   int rowNum(model_->getNumOrigCons() + 1);
   int colNum(model_->getNumOrigVars());
   int *uColIndices(model_->getUpperColInd());
   int *fixedInd(model_->getFixedInd());
   const double *origRowLb(model_->getOrigRowLb());
   const double *origRowUb(model_->getOrigRowUb());
   double gap = (targetGap < model_->getTolerance()) ? 0.0 : targetGap;

   // YX: d2^y <= phi(A^x) + \delta * abs(phi(A^x))
   objUb = objValLL + fabs(objValLL) * gap/100;

   if (!lpSol){
      lpSol = model_->getSolver()->getColSolution();
   }

   if(newOsi){
      double lObjSense(model_->getLowerObjSense());
      double *lObjCoeffs(model_->getLowerObjCoeffs());
      int *lColIndices(model_->getLowerColInd());
      const double *origColLb(model_->getOrigColLb());
      const double *origColUb(model_->getOrigColUb());
      const double *uObjCoeffs(model_->getSolver()->getObjCoefficients());
      
      int intCnt(0);
      double uObjSense(1);

      CoinPackedVector row;
      CoinPackedMatrix matrix = *model_->origConstCoefMatrix_;
      matrix.reverseOrdering();

      CoinPackedMatrix *newMat = new CoinPackedMatrix(false, 0, 0);
      newMat->setDimensions(0, colNum);

      double *colUb = new double[colNum];
      double *colLb = new double[colNum];
      double *rowUb = new double[rowNum];
      double *rowLb = new double[rowNum];   
      int *integerVars = new int[colNum];
      double *objCoeffs = new double[colNum];
      double *newRow = new double[lCols];

      if(feasCheckSolver == "Cbc"){
         nSolver = new OsiCbcSolverInterface();
      }else if(feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
         nSolver = new OsiSymSolverInterface();
#else
         throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
               "setUpPesModel", "MibsBilevel");
#endif
      }else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
         nSolver = new OsiCpxSolverInterface();
#else
         throw CoinError("CPLEX chosen as solver, but it has not been enabled",
               "setUpPesModel", "MibsBilevel");
#endif
      }else{
         throw CoinError("Unknown solver chosen",
               "setUpPesModel", "MibsBilevel");
      }
      
      CoinZeroN(rowUb, rowNum);
      CoinZeroN(rowLb, rowNum);
      CoinZeroN(colUb, colNum);
      CoinZeroN(colLb, colNum);
      CoinZeroN(integerVars, colNum);
      CoinZeroN(objCoeffs, colNum);
      CoinZeroN(newRow, lCols);

      /** Set row bounds **/
      memcpy(rowLb, origRowLb, sizeof(double) * (rowNum - 1));
      memcpy(rowUb, origRowUb, sizeof(double) * (rowNum - 1));

      /** Set col bounds **/
      memcpy(colLb, origColLb, sizeof(double) * colNum);
      memcpy(colUb, origColUb, sizeof(double) * colNum);

      /** Fix linking variables **/
      for(i = 0; i < uCols; i++){
         idx = uColIndices[i];
         if(fixedInd[idx] == 1){
            colLb[idx] = floor(lpSol[idx] + 0.5);
            colUb[idx] = colLb[idx];
         }
      }

      /** Fill in array of integer variables **/
      for(i = 0; i < colNum; i++){
         if(model_->getSolver()->isInteger(i)){
            integerVars[intCnt] = i;
            intCnt ++;
         }
      }

      /** Prepare for value function constraint **/
      CoinDisjointCopyN(lObjCoeffs, lCols, newRow);

      rowUb[rowNum-1] = objUb;
      rowLb[rowNum-1] = -1 * (model_->getSolver()->getInfinity());

      /** Fill in new objective and constraint coeffs **/
      for(i = 0; i < lCols; i++){
         idx = lColIndices[i];
         objCoeffs[idx] = uObjCoeffs[idx]; // pess risk func 0x+d^1y
         row.insert(idx, newRow[i] * lObjSense); // val func d^2y
      }

      /** Set up constraints **/
      for(i = 0; i < rowNum - 1; i++){
         newMat->appendRow(matrix.getVector(i));
      }

      newMat->appendRow(row);  

      nSolver->loadProblem(*newMat, colLb, colUb,
            objCoeffs, rowLb, rowUb);

      for(i = 0; i < intCnt; i++){
         nSolver->setInteger(integerVars[i]);
      }
      
      /** Pessimistic risk function **/
      nSolver->setObjSense(-1 * uObjSense); //1 min; -1 max 
      
      // YX: still need this?
      nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);
   
      delete [] colUb;
      delete [] colLb;
      delete [] rowUb;
      delete [] rowLb;
      delete [] objCoeffs;
      delete [] integerVars;
      delete [] newRow;
      delete newMat;

   }else{

      nSolver = pSolver_;
      
      /** Update value function constraint bound **/
      nSolver->setRowUpper(rowNum-1, objUb); // YX: gap added to phi(A^2x)

      /** Update fixed linking variable bounds **/
      for(i = 0; i < uCols; i++){
         idx = uColIndices[i];
         if(fixedInd[idx] == 1){
            value = floor(lpSol[idx] + 0.5);
            nSolver->setColLower(idx, value);
            nSolver->setColUpper(idx, value);
         }
      }
   }

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

  bool signsChanged(false);
  double etol(model_->etol_);
  int i(0), j(0), index1(0), index2(0);
  double mult(0.0);
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
	  index1 = lRowIndices[i];
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
     
     signsChanged = false;
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
	     index1 = lRowIndices[i];
	     rowLb[i] = origRowLb[index1];
	     rowUb[i] = origRowUb[index1];
	 }
     }
     else{
	 for(i = 0; i < lRows; i++){
	     index1 = lRowIndices[i];
	     if(origRowSense[index1] == 'G'){
		 rowLb[i] = -origRowUb[index1];
		 rowUb[i] = -origRowLb[index1];
	     }
	     else{
		 rowLb[i] = origRowLb[index1];
		 rowUb[i] = origRowUb[index1];
	     }
	 }
     }
     
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
	  mult = 1;
	  if((signsChanged) && (origRowSense[index1] == 'G')){
	      mult = -1;
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
double
MibSBilevel::getRiskFuncVal(double *lowerSol)
{
   int lN(model_->getLowerDim());
   int *lColIndices(model_->getLowerColInd());
   const double uObjSense(model_->solver()->getObjSense());
   const double *uObjCoeffs(model_->solver()->getObjCoefficients());

   int i(0), idx(0);
   double rfVal(0.0);

   for(i = 0; i < lN; i++){
      idx = lColIndices[i];
      if(0){
         std::cout << "lowerSol[" << i << "]: " << lowerSol[i] << std::endl;
         std::cout << "uObjCoeff[" << idx << "]: " << uObjCoeffs[idx] << std::endl;
      }      
      rfVal += uObjCoeffs[idx] * lowerSol[i];
   }

   return rfVal * uObjSense;

}

//#############################################################################
double
MibSBilevel::getUpperObj(double *lowerSol, double *upperSol)
{  
   /** calculate upper-level objective value for (x, hat{y}) or (x, y') **/
   // YX: assuming upperSol passed in are ordered;
   int uN(model_->getUpperDim());
   int *uColIndices(model_->getUpperColInd());
   const double uObjSense(model_->solver()->getObjSense());
   const double *uObjCoeffs(model_->solver()->getObjCoefficients());
   
   if(!upperSol){
      const double *upperSol(model_->solver()->getColSolution());
   }

   int i(0), idx(0);
   double rfuncVal(0.0), uObjVal(0.0);

   for(i = 0; i < uN; i++){
      if(!upperSol){
         idx = uColIndices[i];
      }else{
         idx = i;
      } 
      if(0){
         std::cout << "upperSol[" << i << "]: " << upperSol[i] << std::endl;
         std::cout << "uObjCoeff[" << idx << "]: " << uObjCoeffs[idx] << std::endl;
      }      
      uObjVal += uObjCoeffs[idx] * upperSol[i] * uObjSense;
   }

   rfuncVal = getRiskFuncVal(lowerSol);
   uObjVal += rfuncVal;

   return uObjVal;
}

//#############################################################################
void
    MibSBilevel::addSolutionToSeenLinkingSolutionPool(MibSLinkingPoolTag solTag, std::vector<double>
				   &shouldStoreValues, double objValue)
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
    linkingSolution.lowerSolution.push_back(0);
    linkingSolution.pesSolution.push_back(0); // YX: pessimistic case
    linkingSolution.UBSolution.push_back(0);
    linkingSolution.lowerSolution.clear();
    linkingSolution.pesSolution.clear(); // YX: pessimistic case
    linkingSolution.UBSolution.clear();
    
    switch(solTag){
    case MibSLinkingPoolTagLowerIsInfeasible:
	{
	    linkingSolution.tag = solType;
	    linkingSolution.lowerObjValue = 0.0;
        linkingSolution.pesRFValue = 0.0;
	    linkingSolution.UBObjValue = 0.0;
	    linkingSolution.lowerSolution.push_back(0);
        linkingSolution.pesSolution.push_back(0); 
	    linkingSolution.UBSolution.push_back(0);
	    model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    break;
	}
    case MibSLinkingPoolTagLowerIsFeasible:
	{
	    linkingSolution.tag = solType;
	    linkingSolution.lowerObjValue = objValue;
        linkingSolution.pesRFValue = 0.0;
	    linkingSolution.UBObjValue = 0.0;
	    linkingSolution.lowerSolution = shouldStoreValues;
        linkingSolution.pesSolution.push_back(0);
	    linkingSolution.UBSolution.push_back(0);
	    model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    break;
	}
    case MibSLinkingPoolTagPesIsInfeasible: // YX: pessimistic case; following one tag above;
	{
        model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagPesIsInfeasible;
        break;
	}
    case MibSLinkingPoolTagPesIsFeasible: // YX: pessimistic case
	{
        model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagPesIsFeasible;
        model_->seenLinkingSolutions[linkSol].pesRFValue = objValue;
        model_->seenLinkingSolutions[linkSol].pesSolution.clear();
		model_->seenLinkingSolutions[linkSol].pesSolution = shouldStoreValues;
        break;
	}
    case MibSLinkingPoolTagUBIsSolved:
	{
	    model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagUBIsSolved;
	    model_->seenLinkingSolutions[linkSol].UBObjValue = objValue;
	    if(isProvenOptimal_){
		model_->seenLinkingSolutions[linkSol].UBSolution.clear();
		model_->seenLinkingSolutions[linkSol].UBSolution = shouldStoreValues;
	    }
	    break;
	}
    }
}
	    
	
	
	    






    
 
