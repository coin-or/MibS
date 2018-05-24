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

#include <map>
#include <time.h>
#include "OsiCbcSolverInterface.hpp"
#include "OsiSolverInterface.hpp"

#include "MibSHeuristic.hpp"
#include "MibSModel.hpp"
#include "MibSParams.hpp"
#include "MibSSolution.hpp"
#include "MibSSolTypes.hpp"
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
MibSHeuristic::MibSHeuristic(MibSModel * model)
{

  initializeData(model);

}

//#############################################################################
void 
MibSHeuristic::initializeData(MibSModel * model)
{

  MibSModel_ = model;
  numProblems_ = 0;
  etol_ = 1e-5;
  bestObjVal_ = 0.0; 
  bestSol_ = NULL;

}

//#############################################################################
MibSHeuristic::~MibSHeuristic()
{

  if(bestSol_) delete [] bestSol_;
  //if(MibSModel_) delete MibSModel_;
}

//#############################################################################
void 
MibSHeuristic::findHeuristicSolutions()
{

  MibSModel * model = MibSModel_;

  int msgLevel(model->AlpsPar_->entry(AlpsParams::msgLevel));

  bool useLowerObjHeuristic
    = model->MibSPar_->entry(MibSParams::useLowerObjHeuristic);

  bool useObjCutHeuristic
    = model->MibSPar_->entry(MibSParams::useObjCutHeuristic);

  bool useWSHeuristic
    = model->MibSPar_->entry(MibSParams::useWSHeuristic);

  bool useGreedyHeuristic
    = model->MibSPar_->entry(MibSParams::useGreedyHeuristic);

  int heurFrequency(100);
  int numCallsHeur(model->countIteration_ - 1);
  int numCallObjCutHeur(model->bS_->linkIntegralCount_);
  if(model->bS_->isLinkVarsIntegral_ == true){
      model->bS_->linkIntegralCount_ ++;
  }

  if(numCallsHeur%heurFrequency == 0){
    if(msgLevel > 100){
      if(useLowerObjHeuristic == true){
	std::cout << "lowerObj heuristic is on." << std::endl;
      }
      if(useWSHeuristic == true){
	std::cout << "ws heuristic is on." << std::endl;
      }
      if(useGreedyHeuristic == true){
	std::cout << "greedy heuristic is on." << std::endl;
      }
      std::cout << "Heuristic frequency = " << heurFrequency << std::endl;
    }

    if(useLowerObjHeuristic)
      lowerObjHeuristic();

    if(useWSHeuristic)
      weightedSumsHeuristic();

    if(useGreedyHeuristic)
      greedyHeuristic();
  }

  if((numCallObjCutHeur%heurFrequency == 0) && (useObjCutHeuristic == true) &&
     (model->bS_->isLinkVarsIntegral_ == true)){
      objCutHeuristic();
      if(msgLevel > 100){
	  std::cout << "objCut heuristic is on." << std::endl;
	  std::cout << "Heuristic frequency = " << heurFrequency << std::endl;
      }
  }
	  
      
}

//#############################################################################
void 
MibSHeuristic::lowerObjHeuristic()
{

  /* 
     optimize wrt to lower-level objective 
     over current feasible lp feasible region 
  */

  MibSModel * model = MibSModel_;

  std::string feasCheckSolver(model->MibSPar_->entry
			      (MibSParams::feasCheckSolver));

  int maxThreadsLL(model->MibSPar_->entry
		   (MibSParams::maxThreadsLL));

  int whichCutsLL(model->MibSPar_->entry
		  (MibSParams::whichCutsLL));

  double timeLimit(model->AlpsPar()->entry(AlpsParams::timeLimit));

  int useLinkingSolutionPool(model->MibSPar_->entry
			     (MibSParams::useLinkingSolutionPool));

  OsiSolverInterface * oSolver = model->getSolver();

  OsiSolverInterface * hSolver;

  if (feasCheckSolver == "Cbc"){
    hSolver = new OsiCbcSolverInterface();
  }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
    hSolver = new OsiSymSolverInterface();
#else
    throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
		    "lowerObjHeuristic", "MibSHeuristic");
#endif
  }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
    hSolver = new OsiCpxSolverInterface();
#else
    throw CoinError("CPLEX chosen as solver, but it has not been enabled",
		    "lowerObjHeuristic", "MibSHeuristic");
#endif
  }else{
    throw CoinError("Unknown solver chosen",
		    "lowerObjHeuristic", "MibSHeuristic");
  }

  int i(0), j(0);
  int index(0), solType(0);
  double value(0.0), lObjVal(0.0);
  double remainingTime(0.0), startTimeVF(0.0);
  bool isContainedInLinkingPool(false);
  MibSLinkingPoolTag tagInSeenLinkingPool(MibSLinkingPoolTagIsNotSet);
  double etol(etol_);
  double uObjSense(oSolver->getObjSense());
  double objSense(model->getLowerObjSense());
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int * lColIndices(model->getLowerColInd());
  int * uColIndices(model->getUpperColInd());
  double * lObjCoeffs(model->getLowerObjCoeffs());
  const double * uObjCoeffs(oSolver->getObjCoefficients());
  //CoinPackedMatrix origMatrix = *model->origConstCoefMatrix_;
  int * fixedInd = model->fixedInd_;
  MibSBilevel *bS = model->bS_;
  double * optUpperSolutionOrd = NULL;
  double * lSol = new double[lCols];

  int tCols(oSolver->getNumCols());

  double * nObjCoeffs = new double[tCols];
  CoinZeroN(nObjCoeffs, tCols);

  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    nObjCoeffs[index] = lObjCoeffs[i];
  } 

  //sahar: we use the info of node t like the first version of implementation
  //If we use the info of node t, hSolver may be infeasble.
  hSolver->loadProblem(*oSolver->getMatrixByCol(),
		       oSolver->getColLower(), oSolver->getColUpper(), nObjCoeffs,
		       oSolver->getRowLower(), oSolver->getRowUpper());
  
  for(j = 0; j < tCols; j++){
    if(oSolver->isInteger(j)){
      hSolver->setInteger(j);
    }
  }

  hSolver->setObjSense(objSense);
  
  double cutoff(model->getKnowledgeBroker()->getIncumbentValue());

  if(model->getNumSolutions()){

    CoinPackedVector objCon;
    double rhs = cutoff * uObjSense;

    for(i = 0; i < tCols; i++){
      objCon.insert(i, uObjCoeffs[i] * uObjSense);
    }
    hSolver->addRow(objCon, -hSolver->getInfinity(), rhs);
  }

  remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
  remainingTime = CoinMax(remainingTime, 0.00);

  if (feasCheckSolver == "Cbc"){
    dynamic_cast<OsiCbcSolverInterface *>
      (hSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
    sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (hSolver)->getSymphonyEnvironment();
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
    hSolver->setHintParam(OsiDoReducePrint);
    hSolver->messageHandler()->setLogLevel(0);
    CPXENVptr cpxEnv =
      dynamic_cast<OsiCpxSolverInterface*>(hSolver)->getEnvironmentPtr();
    assert(cpxEnv);
    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
  }

  hSolver->branchAndBound();

  if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					 (dynamic_cast<OsiSymSolverInterface *>
					  (hSolver)->getSymphonyEnvironment()))){
    
    goto TERM_LOWEROBJHEURISTIC;
  }

  if(hSolver->isProvenOptimal()){

    const double * sol = hSolver->getColSolution();

    optUpperSolutionOrd = new double[uCols];

    std::vector<double> linkSol;

    for(i = 0; i < uCols; i++){
      index = uColIndices[i];
      value = sol[index];
      if ((oSolver->isInteger(index)) &&
	  (((value - floor(value)) < etol) ||
	   ((ceil(value) - value) < etol))){
	optUpperSolutionOrd[i] = (double) floor(value + 0.5);
      }
      else{
	optUpperSolutionOrd[i] = value;
      }
      if(fixedInd[index] == 1){
	linkSol.push_back(optUpperSolutionOrd[i]);
      }
    }
  
    if(useLinkingSolutionPool){
      if(model->seenLinkingSolutions.find(linkSol) !=
	 model->seenLinkingSolutions.end()){
	isContainedInLinkingPool = true;
	solType = model->seenLinkingSolutions.find(linkSol)->second.tag;
      }
    }
  
    if(isContainedInLinkingPool == true){
      tagInSeenLinkingPool = static_cast<MibSLinkingPoolTag>(solType);
    }

    if(tagInSeenLinkingPool == MibSLinkingPoolTagUBIsSolved){
      goto TERM_LOWEROBJHEURISTIC;
    }

    //tag cannot be MibSLinkingPoolTagLowerIsInfeasible because
    //hSolver includes all integrality constraints.
    
    double upperObjVal(0.0);

    MibSSolution *mibSol = NULL;

    if(!isContainedInLinkingPool){

      if(bS->lSolver_){
	bS->lSolver_ = bS->setUpModel(oSolver, false, sol);
      }
      else{
	bS->lSolver_ = bS->setUpModel(oSolver, true, sol);
      }

      OsiSolverInterface *lSolver = bS->lSolver_; 

      remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
      remainingTime = CoinMax(remainingTime, 0.00);

      if (feasCheckSolver == "Cbc"){
	dynamic_cast<OsiCbcSolverInterface *>
	  (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
      }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (lSolver)->getSymphonyEnvironment();
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

      startTimeVF = model->broker_->subTreeTimer().getTime();
      lSolver->branchAndBound();
      model->timerVF_ += model->broker_->subTreeTimer().getTime() - startTimeVF; 
      model->counterVF_ ++;

      if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					     (dynamic_cast<OsiSymSolverInterface *>
					      (lSolver)->getSymphonyEnvironment()))){
	goto TERM_LOWEROBJHEURISTIC;
      }

      if (lSolver->isProvenOptimal()){
	lObjVal = lSolver->getObjValue() * objSense;
	const double * lSolTmp = lSolver->getColSolution();
	memcpy(lSol, lSolTmp, sizeof(double) * lCols);
	if(useLinkingSolutionPool){
	  //Adding x_L to set E
	  std::vector<double> shouldStoreValuesLowerSol(lCols);
	  std::copy(lSol, lSol + lCols, shouldStoreValuesLowerSol.begin());
	  addSolutionToSeenLinkingSolutionPoolHeur
	    (linkSol, shouldStoreValuesLowerSol, lObjVal);
	}
      }
      else{
	std::cout << "lSolver cannnot be infeasible! MibSHeuristic, lowerObjHeuristic" << std::endl;
        abort();
      }
    }
    //sahar:else means that linkingSol is in pool and tag is
    //MibSLinkingPoolTagLowerIsFeasible
    else{
      lObjVal = model->seenLinkingSolutions[linkSol].lowerObjValue;
      std::copy(model->seenLinkingSolutions[linkSol].lowerSolution.begin(),
		model->seenLinkingSolutions[linkSol].lowerSolution.end(), lSol);
    }

    double lowerObj = getLowerObj(sol, objSense);

    if(fabs(lObjVal - lowerObj) < etol){
      /** Current solution is bilevel feasible **/
      for(i = 0; i < tCols; i++)
	upperObjVal +=
	  sol[i] * uObjCoeffs[i];

      mibSol = new MibSSolution(hSolver->getNumCols(),
				sol, upperObjVal, model);

      model->storeSolution(BlisSolutionTypeHeuristic, mibSol);

      mibSol = NULL;
    }
    else{// the solution of hSolver is not bilevel feasible, but we
      //try to generate a heuristic solution by using the upper  part
      //of hSolver solution and the solution of lower level problem 
      int numElements(hSolver->getNumCols());
      double * lpSolution = new double[numElements];

      for(i = 0; i < uCols; i++){
	index = uColIndices[i];
        lpSolution[index] = optUpperSolutionOrd[i];
        upperObjVal +=
	  optUpperSolutionOrd[i] * uObjCoeffs[index];
      }

      for(i = 0; i < lCols; i++){
	index = lColIndices[i];
        value = lSol[i];
        if ((oSolver->isInteger(index)) &&
	    (((value - floor(value)) < etol) ||
	     ((ceil(value) - value) < etol))){
	  lpSolution[index] = (double) floor(value + 0.5);
	}
	else{
	  lpSolution[index] = value;
	}
	upperObjVal +=
	  lpSolution[index] * uObjCoeffs[index];
      }

      if(model->checkUpperFeasibility(lpSolution)){
	mibSol = new MibSSolution(numElements, lpSolution,
				  upperObjVal * oSolver->getObjSense(), model);
	model->storeSolution(BlisSolutionTypeHeuristic, mibSol);
	mibSol = NULL;
      }
      delete [] lpSolution;
    }
  }

 TERM_LOWEROBJHEURISTIC:
  delete hSolver;
  delete [] lSol;
  delete [] nObjCoeffs;
  if(optUpperSolutionOrd){
    delete [] optUpperSolutionOrd;
  }
}
      
//#############################################################################
void 
MibSHeuristic::objCutHeuristic()
{
  MibSModel * model = MibSModel_;
  
  std::string feasCheckSolver(model->MibSPar_->entry
			      (MibSParams::feasCheckSolver));

  int maxThreadsLL(model->MibSPar_->entry
		   (MibSParams::maxThreadsLL));

  int whichCutsLL(model->MibSPar_->entry
		  (MibSParams::whichCutsLL));

  double timeLimit(model->AlpsPar()->entry(AlpsParams::timeLimit));

  int useLinkingSolutionPool(model->MibSPar_->entry
			     (MibSParams::useLinkingSolutionPool));


  OsiSolverInterface * oSolver = model->getSolver();

  int i(0), j(0);
  int index(0), solType(0);
  double remainingTime(0.0), startTimeVF(0.0);
  double rhsAddedRow(0.0), value(0.0);
  double lObjCoeff(0.0), lObjVal(0.0), objCoeff(0.0);
  bool shouldSolveLowerProb(false);
  MibSBilevel *bS = model->bS_;
  double etol(etol_);
  int uCols(model->getUpperDim());
  int lCols(model->getLowerDim());
  int tCols(oSolver->getNumCols());  
  double objSense(model->getLowerObjSense());
  int * uColIndices = model->getUpperColInd();
  int * lColIndices = model->getLowerColInd();
  const double * uObjCoeffs = oSolver->getObjCoefficients();
  double * lObjCoeffs = model->getLowerObjCoeffs();
  int * fixedInd = model->fixedInd_;
  MibSLinkingPoolTag tagInSeenLinkingPool(bS->tagInSeenLinkingPool_);

  double *lSol = NULL;
  double *optUpperSolutionOrd = NULL; 
  CoinPackedVector objCon;

  std::vector<double> linkSol;
  std::vector<double> shouldStoreValuesLowerSol(lCols);

  OsiSolverInterface * hSolver = 0;
  OsiSolverInterface *lSolver = 0;

  if(bS->isLinkVarsIntegral_ == false){
    goto TERM_OBJCUTHEURISTIC;
  }

  if(useLinkingSolutionPool){
    if(tagInSeenLinkingPool == MibSLinkingPoolTagLowerIsInfeasible){
      goto TERM_OBJCUTHEURISTIC;
    }
    else if(tagInSeenLinkingPool == MibSLinkingPoolTagIsNotSet){
      shouldSolveLowerProb = true;
    }
  }
  else{
    if((bS->isLowerSolved_ == true) && (bS->isUBSolved_ == false) &&
       (bS->isProvenOptimal_ == false)){
      goto TERM_OBJCUTHEURISTIC;
    }
    else if(bS->isLowerSolved_ == false){
      shouldSolveLowerProb = true;
    }
  }
  
  lSol = new double[lCols];

  if(shouldSolveLowerProb){
    if(bS->lSolver_){
      bS->lSolver_ = bS->setUpModel(oSolver, false);
    }
    else{
      bS->lSolver_ = bS->setUpModel(oSolver, true);
    }

    lSolver = bS->lSolver_;

    remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
    remainingTime = CoinMax(remainingTime, 0.00);

    if (feasCheckSolver == "Cbc"){
      dynamic_cast<OsiCbcSolverInterface *>
	(lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
      sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	(lSolver)->getSymphonyEnvironment();
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

    startTimeVF = model->broker_->subTreeTimer().getTime();
    lSolver->branchAndBound();
    model->timerVF_ += model->broker_->subTreeTimer().getTime() - startTimeVF;
    model->counterVF_ ++;

    if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					   (dynamic_cast<OsiSymSolverInterface *>
					    (lSolver)->getSymphonyEnvironment()))){
      goto TERM_OBJCUTHEURISTIC;
    }

    if(lSolver->isProvenOptimal()){
      rhsAddedRow = lSolver->getObjValue() * objSense;
      if(useLinkingSolutionPool){
	//Adding x_L to set E
	double *upperSolutionOrd = bS->upperSolutionOrd_;
        for(i = 0; i < uCols; i++){
	  index = uColIndices[i];
	  if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd[i]);
	  }
	}
        const double * lSolTmp = lSolver->getColSolution();
        memcpy(lSol, lSolTmp, sizeof(double) * lCols);
        std::copy(lSol, lSol + lCols, shouldStoreValuesLowerSol.begin());
        addSolutionToSeenLinkingSolutionPoolHeur
	  (linkSol, shouldStoreValuesLowerSol, rhsAddedRow);
      }
    }
    else{
      goto TERM_OBJCUTHEURISTIC;
    }
  }

  if (feasCheckSolver == "Cbc"){
    hSolver = new OsiCbcSolverInterface();
  }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
    hSolver = new OsiSymSolverInterface();
#else
    throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
		    "objCutHeuristic", "MibSHeuristic");
#endif
  }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
    hSolver = new OsiCpxSolverInterface();
#else
    throw CoinError("CPLEX chosen as solver, but it has not been enabled",
		    "objCutHeuristic", "MibSHeuristic");
#endif
  }else{
    throw CoinError("Unknown solver chosen",
		    "objCutHeuristic", "MibSHeuristic");
  }

  hSolver->loadProblem(*oSolver->getMatrixByCol(),
		       oSolver->getColLower(), oSolver->getColUpper(),
		       uObjCoeffs, oSolver->getRowLower(),
		       oSolver->getRowUpper());

  for(j = 0; j < tCols; j++){
    if(oSolver->isInteger(j))
      hSolver->setInteger(j);
  }

  //CoinPackedVector objCon;

  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    objCon.insert(index, lObjCoeffs[i] * objSense);
  }

  if(shouldSolveLowerProb == false){
    double * optLowerSolutionOrd = bS->optLowerSolutionOrd_;
    for(i = 0; i < lCols; i++){
      rhsAddedRow += optLowerSolutionOrd[i] * lObjCoeffs[i] * objSense;
    }
  }

  hSolver->addRow(objCon, -hSolver->getInfinity(), rhsAddedRow);

  remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
  remainingTime = CoinMax(remainingTime, 0.00);

  if (feasCheckSolver == "Cbc"){
        dynamic_cast<OsiCbcSolverInterface *>
	  (hSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
        sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (hSolver)->getSymphonyEnvironment();
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
    hSolver->setHintParam(OsiDoReducePrint);
    hSolver->messageHandler()->setLogLevel(0);
        CPXENVptr cpxEnv =
	  dynamic_cast<OsiCpxSolverInterface*>(hSolver)->getEnvironmentPtr();
	assert(cpxEnv);
	CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
  }

  hSolver->branchAndBound();

  if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					 (dynamic_cast<OsiSymSolverInterface *>
					  (hSolver)->getSymphonyEnvironment()))){

    goto TERM_OBJCUTHEURISTIC;
  }

  if(hSolver->isProvenOptimal()){
    const double * sol = hSolver->getColSolution();

    optUpperSolutionOrd = new double[uCols];

    for(i = 0; i < uCols; i++){
      index = uColIndices[i];
      value = sol[index];
      if ((oSolver->isInteger(index)) &&
	  (((value - floor(value)) < etol) ||
	   ((ceil(value) - value) < etol))){
	optUpperSolutionOrd[i] = (double) floor(value + 0.5);
      }
      else{
	optUpperSolutionOrd[i] = value;
      }
      if(fixedInd[index] == 1){
	linkSol.push_back(optUpperSolutionOrd[i]);
      }
    }

    bool isContainedInLinkingPool(false);
    tagInSeenLinkingPool = MibSLinkingPoolTagIsNotSet;

    if(useLinkingSolutionPool){
      if(model->seenLinkingSolutions.find(linkSol) !=
	 model->seenLinkingSolutions.end()){
	isContainedInLinkingPool = true;
	solType = model->seenLinkingSolutions.find(linkSol)->second.tag;
      }
    }

    if(isContainedInLinkingPool == true){
      tagInSeenLinkingPool = static_cast<MibSLinkingPoolTag>(solType);
    }

    if(tagInSeenLinkingPool == MibSLinkingPoolTagUBIsSolved){
      goto TERM_OBJCUTHEURISTIC;
    }

    //tag cannot be MibSLinkingPoolTagLowerIsInfeasible because
    //hSolver includes all integrality constraints.

    double upperObjVal(0.0);

    MibSSolution *mibSol = NULL;

    if(!isContainedInLinkingPool){

      if(bS->lSolver_){
	bS->lSolver_ = bS->setUpModel(oSolver, false, sol);
      }
      else{
	bS->lSolver_ = bS->setUpModel(oSolver, true, sol);
      }

      lSolver = bS->lSolver_;

      remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
      remainingTime = CoinMax(remainingTime, 0.00);

      if (feasCheckSolver == "Cbc"){
	dynamic_cast<OsiCbcSolverInterface *>
	  (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
      }else if (feasCheckSolver == "SYMPHONY"){
	#if COIN_HAS_SYMPHONY
	sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (lSolver)->getSymphonyEnvironment();
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

      startTimeVF = model->broker_->subTreeTimer().getTime();
      lSolver->branchAndBound();
      model->timerVF_ += model->broker_->subTreeTimer().getTime() - startTimeVF;
      model->counterVF_ ++;

      if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					     (dynamic_cast<OsiSymSolverInterface *>
					      (lSolver)->getSymphonyEnvironment()))){
	goto TERM_OBJCUTHEURISTIC;
      }

      if(lSolver->isProvenOptimal()){
	lObjVal = lSolver->getObjValue() * objSense;
	const double * lSolTmp1 = lSolver->getColSolution();
	memcpy(lSol, lSolTmp1, sizeof(double) * lCols);
	if(useLinkingSolutionPool){
	  //Adding x_L to set E
	  std::copy(lSol, lSol + lCols, shouldStoreValuesLowerSol.begin());
	    addSolutionToSeenLinkingSolutionPoolHeur
	      (linkSol, shouldStoreValuesLowerSol, lObjVal);
	}
      }
      else{
	std::cout << "lSolver cannnot be infeasible! MibSHeuristic, objCutHeuristic" << std::endl;
	abort();
      }
    }
    //sahar:else means that linkingSol is in pool and tag is
    //MibSLinkingPoolTagLowerIsFeasible
    else{
      lObjVal = model->seenLinkingSolutions[linkSol].lowerObjValue;
      std::copy(model->seenLinkingSolutions[linkSol].lowerSolution.begin(),
		model->seenLinkingSolutions[linkSol].lowerSolution.end(), lSol);
    }

    double lowerObj = getLowerObj(sol, objSense);

    if(fabs(lObjVal - lowerObj) < etol){
      //Current solution is bilevel feasible
      for(i = 0; i < tCols; i++)
	upperObjVal +=
	  sol[i] * uObjCoeffs[i];

      mibSol = new MibSSolution(hSolver->getNumCols(),
				sol, upperObjVal, model);

      model->storeSolution(BlisSolutionTypeHeuristic, mibSol);

      mibSol = NULL;
    }
    else{// the solution of hSolver is not bilevel feasible, but we
      //try to generate a heuristic solution by using the upper  part
      //of hSolver solution and the solution of lower level problem
      int numElements(hSolver->getNumCols());
      double * lpSolution = new double[numElements];

      for(i = 0; i < uCols; i++){
	index = uColIndices[i];
	lpSolution[index] = optUpperSolutionOrd[i];
	upperObjVal += optUpperSolutionOrd[i] * uObjCoeffs[index];
      }

      for(i = 0; i < lCols; i++){
	index = lColIndices[i];
	value = lSol[i];
	if ((oSolver->isInteger(index)) &&
	    (((value - floor(value)) < etol) ||
	     ((ceil(value) - value) < etol))){
	  lpSolution[index] = (double) floor(value + 0.5);
	}
	else{
	  lpSolution[index] = value;
	}
	upperObjVal +=
	  lpSolution[index] * uObjCoeffs[index];
      }

      if(model->checkUpperFeasibility(lpSolution)){
	mibSol = new MibSSolution(numElements, lpSolution,
				  upperObjVal * oSolver->getObjSense(), model);
	model->storeSolution(BlisSolutionTypeHeuristic, mibSol);
	mibSol = NULL;
      }
      delete [] lpSolution;
    }
  }

 TERM_OBJCUTHEURISTIC:
  if(hSolver){
    delete hSolver;
  }
  if(lSol){
    delete [] lSol;
  }
  if(optUpperSolutionOrd){
    delete [] optUpperSolutionOrd;
  }
}

//#############################################################################
void 
MibSHeuristic::greedyHeuristic()
{

  MibSModel * model = MibSModel_;
  //OsiSolverInterface * oSolver = model->getSolver();
  OsiSolverInterface * oSolver = model->solver();
  
  double uObjSense(oSolver->getObjSense());
  double lObjSense(model->getLowerObjSense());  
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int * uColIndices = model->getUpperColInd();
  int * lColIndices = model->getLowerColInd();
  double * lObjCoeffs = model->getLowerObjCoeffs();
  double * intCost = model->getInterdictCost();
  double intBudget = model->getInterdictBudget();

  int tCols(uCols + lCols); 

  //assert(tCols == oSolver->getNumCols());

  int i(0), ind_min_wt(0);
  double usedBudget(0.0); 
  double * fixedVars = new double[lCols];
  double * testsol = new double[tCols];
  CoinZeroN(fixedVars, lCols);
  CoinZeroN(testsol, tCols);

  std::multimap<double, int> lObjCoeffsOrd;

  for(i = 0; i < lCols; i++)
    lObjCoeffsOrd.insert(std::pair<double, int>(lObjCoeffs[i] * lObjSense, i));
  
  if(!bestSol_)
    bestSol_ = new double[tCols];

  //initialize the best solution information
  //bestObjVal_ = model->getSolver()->getInfinity() * uObjSense;
  //CoinZeroN(bestSol_, tCols);

  std::multimap<double, int>::iterator iter;
  //std::multimap<double, int>::iterator first;
  //std::multimap<double, int>::iterator last;
  //int dist = std::distance(first, last);
  SRANDOM((unsigned) time(NULL));

  int randchoice(0); 
  if(0)
    std::cout << "randchoice " << randchoice << std::endl;


  double cost(0.0);

  //starting from the largest, fix corr upper-level variables
  //then, with these fixed, solve the lower-level problem
  //this yields a feasible solution

  iter = lObjCoeffsOrd.begin(); 
  
  while((usedBudget < intBudget) && (iter != lObjCoeffsOrd.end())){
    
    ind_min_wt = iter->second;
    cost = intCost[ind_min_wt];
    testsol[uColIndices[ind_min_wt]] = 1.0;
    double min_wt = iter->first;
    
    if(0){
      std::cout << "upper: " << ind_min_wt << " " 
		<< uColIndices[ind_min_wt] << " "  
		<< oSolver->getColUpper()[uColIndices[ind_min_wt]] << " " 
		<< oSolver->getColLower()[uColIndices[ind_min_wt]] << std::endl;
      
      std::cout << "lower: " << ind_min_wt << " " 
		<< lColIndices[ind_min_wt] << " "  
		<< oSolver->getColUpper()[lColIndices[ind_min_wt]] << std::endl;
    }

    //if((oSolver->getColUpper()[uColIndices[ind_min_wt]] == 1.0) 
       //&& (oSolver->getColUpper()[lColIndices[ind_min_wt]] > 0)){
    if(oSolver->getColUpper()[uColIndices[ind_min_wt]] > etol_){ 
      
      //if(((usedBudget + cost) <= intBudget) 
      // && checkLowerFeasibility(oSolver, testsol)){
      if((usedBudget + cost) <= intBudget){
	
	//FIXME: SHOULD BE CHECKING FOR CURRENT BOUNDS HERE  
	//fix the corresponding upper-level variable to 1
	 randchoice = RANDOM() % 2;
	if(0)
	  std::cout << "randchoice " << randchoice << std::endl;
	if(randchoice){
	  fixedVars[ind_min_wt] = 1.0;
	  usedBudget += intCost[ind_min_wt];
	}
      }
    }
    else{
      
      testsol[uColIndices[ind_min_wt]] = 0;
      //break;
      
    }
    iter++;
  }
  
  
   // now we find a feasible solution by fixing upper-level vars
  //  and solving the lower-level problem
  
  
  double * incumbentSol = new double[tCols];
  double * colsol = new double[tCols];

  CoinZeroN(colsol, tCols);

  for(i = 0; i < uCols; i++){
    colsol[uColIndices[i]] = fixedVars[i];
    if(fixedVars[i] == 1.0)
      if(0)
	std::cout << "fixed " << i << std::endl;
  }

  bool isTimeLimReached(false);
  bool foundFeasible(true);

  bfSol sol = getBilevelSolution(colsol, lObjSense * oSolver->getInfinity(),
				 isTimeLimReached, foundFeasible);

  if(sol.getColumnSol()){
    double incumbentObjVal = sol.getObjVal();
    CoinCopyN(sol.getColumnSol(), tCols, incumbentSol);
    
    MibSSolution * mibSol = new MibSSolution(tCols,
					     incumbentSol,
					     incumbentObjVal,
					     model);
    
    model->storeSolution(BlisSolutionTypeHeuristic, mibSol);
  }

  //bestObjVal_ = incumbentObjVal;
  //CoinCopyN(incumbentSol, tCols, bestSol_);

  delete [] incumbentSol;
  delete [] testsol;
  //delete [] colsol;
  //delete [] fixedVars;
  //delete sol;
}

//#############################################################################
void 
MibSHeuristic::weightedSumsHeuristic()
{

  std::queue< SolPair > probQueue;
  std::list< double > BETAS;
  std::map<double, mcSol> mcSolutions;
  //mcSol sol0 = mcSol();
  //mcSol sol1 = mcSol();
  //mcSol sol = mcSol();
  double beta(0.0);
  double slope(0.0);
  double etol(etol_);
  bool foundSolution(true);

  BETAS.push_back(1.0);
  BETAS.push_back(0.0);

  /** First solve the subproblem for beta = 1 and beta = 0 **/
  
  mcSol sol1 = solveSubproblem(1.0, foundSolution); // beta = 1
  numProblems_++;
  if(foundSolution == false){
    return;
  }
  mcSolutions.insert(std::make_pair(1.0, sol1));

  int i(0);
  //int tCols(MibSModel_->getLowerDim() + MibSModel_->getUpperDim());
  int tCols(MibSModel_->getSolver()->getNumCols());
  
  if(0){
    std::cout << "supported solution added." << std::endl;
    std::cout << "support obj value: " 
	      << sol1.getObjPair().second << std::endl;
    std::cout << "upper obj value: "  
	      << sol1.getObjPair().first << std::endl;
    for(i = 0; i < tCols; i++){
      std::cout << "sol1[" << i << "]: " 
		<< sol1.getColumnSol()[i] << std::endl;
    }
  }
  
  mcSol sol0 = solveSubproblem(0.0, foundSolution);// beta = 0
  numProblems_++;
  if(foundSolution == false){
    return;
  } 
  mcSolutions.insert(std::make_pair(0.0, sol0));

  if(0){
    std::cout << "supported solution added." << std::endl;
    std::cout << "support obj value: " 
	      << sol0.getObjPair().second << std::endl;
    std::cout << "upper obj value: "  
	      << sol0.getObjPair().first << std::endl;
     for(i = 0; i < tCols; i++){
      std::cout << "sol0[" << i << "]" 
		<< sol0.getColumnSol()[i] << std::endl;
    }
  }

  probQueue.push(SolPair(sol1.getObjPair(), sol0.getObjPair()));

  double delta_y(0.0);
  double delta_x(0.0);

  while(probQueue.size() > 0){

    /** Get the first pair in the queue **/
    SolPair sP = probQueue.front();
    /** Remove the first pair **/
    probQueue.pop();
   
    delta_y = sP.getFirst().second - sP.getSecond().second;
    delta_x = sP.getSecond().first - sP.getFirst().first;

    if((fabs(delta_y) > etol) && (fabs(delta_x) > etol)){
      
      slope = delta_y / delta_x;
      beta = slope/(1 + slope);
      
      /** should probably add a tolerance here **/
      
      if(find(BETAS.begin(), BETAS.end(), beta) == BETAS.end()){
	  //std::cout << "Solving with beta = " << beta << std::endl;
	
	BETAS.push_back(beta);
	mcSol sol = solveSubproblem(beta, foundSolution);
	if(foundSolution == false){
	 return;
	} 
	numProblems_++;
      
	/*******************************************/ 
	/* If the outcome does not already exist,  */ 
	/* add it to the solutions list and create */ 
	/* two new solution pairs to be analyzed   */
	/*******************************************/ 
	
	double pair1(sol.getObjPair().first - sP.getFirst().first);
	double pair2(sol.getObjPair().second - sP.getFirst().second);
	double pair3(sol.getObjPair().first - sP.getSecond().first);
	double pair4(sol.getObjPair().second - sP.getSecond().second);

	if(((fabs(pair1) > etol) && (fabs(pair2) > etol)) 
	   && ((fabs(pair3) > etol) && (fabs(pair4) > etol))){
	    
	    mcSolutions.insert(std::make_pair(beta, sol));
	    
	    if(0){
	      std::cout << "supported solution added." << std::endl;
	      std::cout << "support obj value: " 
			<< sol.getObjPair().second << std::endl;
	      std::cout << "upper obj value: "  
			<< sol.getObjPair().first << std::endl;
	      for(i = 0; i < tCols; i++){
		std::cout << "sol[" << i << "]: " 
		<< sol.getColumnSol()[i] << std::endl;
	      }
	    }


	    probQueue.push(SolPair(sP.getFirst(), sol.getObjPair()));
	    probQueue.push(SolPair(sol.getObjPair(), sP.getSecond()));
	    
	  }
      }
      /*else{
	std::cout << "Repeated beta value.  Skipping problem pair." <<std::endl;
	}*/
    }
  }

  createBilevelSolutions(mcSolutions);

}

//#############################################################################
void 
MibSHeuristic::createBilevelSolutions(std::map<double, mcSol> mcSolutions)
{

  MibSModel * model = MibSModel_;
  OsiSolverInterface * oSolver = model->getSolver(); 
  double uObjSense(oSolver->getObjSense());
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());

  //int tCols(uCols + lCols);
  int tCols(oSolver->getNumCols());

  double incumbentObjVal(model->getSolver()->getInfinity() * uObjSense);
  double * incumbentSol = new double[tCols];

  int i(0);
  bool shouldStoreSol(false), isTimeLimReached(false);

  if(!bestSol_)
    bestSol_ = new double[tCols];

  //initialize the best solution information
  //bestObjVal_ = model->getSolver()->getInfinity() * uObjSense;
  //CoinZeroN(bestSol_, tCols);

  std::map<double, mcSol >::iterator iter = mcSolutions.begin();

  for(iter = mcSolutions.begin(); iter != mcSolutions.end(); iter++){
    
    mcSol tmpsol = iter->second;
    const double * colsol = tmpsol.getColumnSol();
    double origLower = tmpsol.getObjPair().second;
    bool foundFeasible(true);
    
    bfSol sol = getBilevelSolution(colsol, origLower, isTimeLimReached,
				   foundFeasible);

    if(isTimeLimReached == true){
      break;
    }

    if(foundFeasible == true){
      if(sol.getObjVal() < incumbentObjVal){
	incumbentObjVal = sol.getObjVal();
	if(sol.getColumnSol()){
	  shouldStoreSol = true;
	  CoinCopyN(sol.getColumnSol(), tCols, incumbentSol);
	}
	else{
	  shouldStoreSol = false;
	}
      }
    }
  }
  
  if(shouldStoreSol == true){
      //double cutOff = model->getKnowledgeBroker()->getIncumbentValue();
      //if(incumbentObjVal < cutOff){
	  MibSSolution * mibSol = new MibSSolution(tCols,
						   incumbentSol,
						   incumbentObjVal,
						   model);
	  model->storeSolution(BlisSolutionTypeHeuristic, mibSol);
	  //}
  }

  //need to add this solution to mibssolutions instead
  //bestObjVal_ = incumbentObjVal;
  //CoinCopyN(incumbentSol, tCols, bestSol_);

  //delete bfSol;
  delete [] incumbentSol;

}

//#############################################################################
bool
MibSHeuristic::checkUpperFeasibility(double * solution)
{

  bool upperFeasible(true);
  MibSModel * model = MibSModel_;
  int * uRowIndices = model->getUpperRowInd();
  int uRows(model->getUpperRowNum());
  const double * origRowLb = model->getOrigRowLb();
  const double * origRowUb = model->getOrigRowUb();
  const CoinPackedMatrix * matrix = model->getSolver()->getMatrixByRow();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();


  double lhs(0.0);
  int i(0), j(0), index1(0), index2(0), start(0), end(0);

  for(i = 0; i < uRows; i++){
    index1 = uRowIndices[i];
    start = matStarts[index1];
    end = start + matrix->getVectorSize(index1);
    for(j = start; j < end; j++){
      index2 = matIndices[j];
      lhs += matElements[j] * solution[index2];
    }
    if((origRowLb[index1] > lhs) || (lhs > origRowUb[index1]))
      upperFeasible = false;
    lhs = 0.0;
  }

  return upperFeasible;
}

//#############################################################################
bool
MibSHeuristic::checkLowerFeasibility1(double * solution)
{

  bool lowerFeasible(true);
  MibSModel * model = MibSModel_;
  int * lRowIndices = model->getLowerRowInd();
  int lRows(model->getLowerRowNum());
  const double * origRowLb = model->getOrigRowLb();
  const double * origRowUb = model->getOrigRowUb();
  const CoinPackedMatrix * matrix = model->getSolver()->getMatrixByRow();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();


  double lhs(0.0);
  int i(0), j(0), index1(0), index2(0), start(0), end(0);

  for(i = 0; i < lRows; i++){
    index1 = lRowIndices[i];
    start = matStarts[index1];
    end = start + matrix->getVectorSize(index1);
    for(j = start; j < end; j++){
      index2 = matIndices[j];
      lhs += matElements[j] * solution[index2];
    }
    if((origRowLb[index1] > lhs) || (lhs > origRowUb[index1]))
      lowerFeasible = false;
    lhs = 0.0;
  }

  return lowerFeasible;
}

//#############################################################################
bool
MibSHeuristic::checkLowerFeasibility(OsiSolverInterface * si, 
				     double * solution)
{

  MibSModel * model = MibSModel_;
  OsiSolverInterface * lSolver = model->bS_->setUpModel(si, true, solution);

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

  if(lSolver->isProvenOptimal())
    return true;
  else
    return false;
}


//#############################################################################
bfSol
MibSHeuristic::getBilevelSolution(const double * sol, double origLower,
				  bool &isTimeLimReached, bool &foundFeasible)
{

  //Find a bilevel feasible solution by solving the LL problem
  //for a fixed UL solution, given by the UL portion of sol

  MibSModel * model = MibSModel_;
  MibSBilevel *bS = model->bS_;
  OsiSolverInterface * oSolver = model->getSolver();

  std::string feasCheckSolver(model->MibSPar_->entry
			      (MibSParams::feasCheckSolver));

  int maxThreadsLL(model->MibSPar_->entry
		   (MibSParams::maxThreadsLL));

  int whichCutsLL(model->MibSPar_->entry
		  (MibSParams::whichCutsLL));

  double timeLimit(model->AlpsPar()->entry(AlpsParams::timeLimit));

  int useLinkingSolutionPool(model->MibSPar_->entry
			     (MibSParams::useLinkingSolutionPool));

  int i(0);
  int index(0), solType(0);
  double value(0.0), roundedValue(0.0);
  double remainingTime(0.0), startTimeVF(0.0);
  double objVal(0.0), lowerObj(0.0);
  bool isContainedInLinkingPool(false);
  MibSLinkingPoolTag tagInSeenLinkingPool(MibSLinkingPoolTagIsNotSet);
  double etol(etol_);
  int uCols(model->getUpperDim());
  int lCols(model->getLowerDim());
  int tCols(oSolver->getNumCols());
  int *uColIndices = model->getUpperColInd();
  int *lColIndices = model->getLowerColInd();
  int *fixedInd = model->fixedInd_;
  const double *uObjCoeff = oSolver->getObjCoefficients();

  double * colSol = new double[tCols];
  double *lSol = new double[lCols];
  std::vector<double> linkSol; 
  //bfSol * bfsol;

  if(useLinkingSolutionPool){
    for(i = 0; i < uCols; i++){
      index = uColIndices[i];
      if(fixedInd[index] == 1){
	value = sol[index];
	if((oSolver->isInteger(index)) &&
	   (((value - floor(value)) < etol) ||
	    ((ceil(value) - value) < etol))){
	  roundedValue = (double) floor(value + 0.5);
	}
	else{
	  roundedValue = value;
	}
	linkSol.push_back(roundedValue);
      }
    }
    
    if(model->seenLinkingSolutions.find(linkSol) !=
       model->seenLinkingSolutions.end()){
      isContainedInLinkingPool = true;
      solType = model->seenLinkingSolutions.find(linkSol)->second.tag;
    }

    if(isContainedInLinkingPool == true){
      tagInSeenLinkingPool = static_cast<MibSLinkingPoolTag>(solType);
    }

    if(tagInSeenLinkingPool == MibSLinkingPoolTagUBIsSolved){
      objVal = model->seenLinkingSolutions[linkSol].UBObjValue;
      delete [] colSol;
      delete [] lSol;
      return bfSol(objVal, NULL, 0);
    }
  //tag cannot be MibSLinkingPoolTagLowerIsInfeasible because all
  //integrality constraints are included.
  }

  if(!isContainedInLinkingPool){
    if(bS->lSolver_){
      bS->lSolver_ = bS->setUpModel(oSolver, false, sol);
    }
    else{
      bS->lSolver_ = bS->setUpModel(oSolver, true, sol);
    }

    OsiSolverInterface *lSolver = bS->lSolver_;

    remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
    remainingTime = CoinMax(remainingTime, 0.00);

    if (feasCheckSolver == "Cbc"){
      dynamic_cast<OsiCbcSolverInterface *>
	(lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }else if (feasCheckSolver == "SYMPHONY"){
      #if COIN_HAS_SYMPHONY
      sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	(lSolver)->getSymphonyEnvironment();
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

    startTimeVF = model->broker_->subTreeTimer().getTime();
    lSolver->branchAndBound();
    model->timerVF_ += model->broker_->subTreeTimer().getTime() - startTimeVF;
    model->counterVF_ ++;

    if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					   (dynamic_cast<OsiSymSolverInterface *>
					    (lSolver)->getSymphonyEnvironment()))){
      delete [] colSol;
      delete [] lSol;
      isTimeLimReached = true;
      return bfSol(objVal, NULL, 0);
    }

    if(lSolver->isProvenOptimal()){
      lowerObj = lSolver->getObjValue();
      const double * lSolTmp = lSolver->getColSolution();
      memcpy(lSol, lSolTmp, sizeof(double) * lCols);
      if(useLinkingSolutionPool){
	//Adding x_L to set E
	std::vector<double> shouldStoreValuesLowerSol(lCols);
	std::copy(lSol, lSol + lCols, shouldStoreValuesLowerSol.begin());
	addSolutionToSeenLinkingSolutionPoolHeur
	  (linkSol, shouldStoreValuesLowerSol, lowerObj);
      }
    }
    else{
      std::cout << "lSolver cannnot be infeasible! MibSHeuristic, WSHeuristic" << std::endl;
      abort();
    }
  }
  //sahar:else means that linkingSol is in pool and tag is
  //MibSLinkingPoolTagLowerIsFeasible
  else{
    lowerObj = model->seenLinkingSolutions[linkSol].lowerObjValue;
    std::copy(model->seenLinkingSolutions[linkSol].lowerSolution.begin(),
	      model->seenLinkingSolutions[linkSol].lowerSolution.end(), lSol);
  }

  if(fabs(origLower - lowerObj) < etol){
    memcpy(colSol, sol, sizeof(double) * tCols);
  }
  else{
    memcpy(colSol, sol, sizeof(double) * tCols);
    for(i = 0; i < lCols; i++){
      index = lColIndices[i];
      colSol[index] = lSol[i];
    }
    foundFeasible = model->checkUpperFeasibility(colSol);
  }

  if(foundFeasible == true){
      for(i = 0; i < tCols; i++){
	  objVal += colSol[i] * uObjCoeff[i];
      }
  }

  delete [] lSol;
  return bfSol(objVal, colSol, tCols);
}

//#############################################################################
bfSol*
MibSHeuristic::getBilevelSolution1(const double * sol)
{

  //
  // Find a bilevel feasible solution by solving the LL problem
  // for a fixed UL solution, given by the UL portion of sol
  


  MibSModel * model = MibSModel_;
  OsiSolverInterface * oSolver = model->getSolver();
  OsiSolverInterface * lSolver = new OsiCbcSolverInterface(oSolver);
  
  //double uObjSense(model->getSolver()->getObjSense());
  double lObjSense(model->getLowerObjSense());  
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int * lColIndices = model->getLowerColInd();
  int uRowNum = model->getUpperRowNum();
  int lRowNum = model->getLowerRowNum();
  int * uRowIndices = model->getUpperRowInd();
  int * lRowIndices = model->getLowerRowInd();
  int * uColIndices = model->getUpperColInd();
  double * lObjCoeffs = model->getLowerObjCoeffs();

  int tCols(uCols + lCols); 
  int i(0), index(0);

  // delete the UL rows 
  lSolver->deleteRows(uRowNum, uRowIndices);

  // Fix the UL variables to their current value in sol 

  for(i = 0; i < uCols; i++){
    index = uColIndices[i];
    lSolver->setColLower(index, sol[index]);
    lSolver->setColUpper(index, sol[index]);
  }

  // Set the objective to the LL objective coefficients 

  double * nObjCoeffs = new double[tCols];
  CoinZeroN(nObjCoeffs, tCols);
      
  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    nObjCoeffs[index] = lObjCoeffs[i] * lObjSense;
  }
  
  lSolver->setObjective(nObjCoeffs);

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

  if(lSolver->isProvenOptimal()){

    double objVal(0.0);

    for(i = 0; i < tCols; i++)
      objVal += lSolver->getColSolution()[i] * oSolver->getObjCoefficients()[i];
    
    double * colsol = new double[tCols];
    CoinCopyN(lSolver->getColSolution(), tCols, colsol);
 
    bfSol *bfsol = new bfSol(objVal, colsol, tCols);
    delete lSolver;
    return bfsol;
  }
  else{
    delete lSolver;
    return NULL;
  }
}

//#############################################################################
mcSol 
MibSHeuristic::solveSubproblem(double beta, bool &foundSolution)
{
  //optimize wrt to weighted upper-level objective
  //over current feasible lp feasible region 

  MibSModel * model = MibSModel_;
  OsiSolverInterface * oSolver = model->getSolver();

  std::string feasCheckSolver(model->MibSPar_->entry
			      (MibSParams::feasCheckSolver));

  int maxThreadsLL(model->MibSPar_->entry
		   (MibSParams::maxThreadsLL));

  int whichCutsLL(model->MibSPar_->entry
		  (MibSParams::whichCutsLL));

  double timeLimit(model->AlpsPar()->entry(AlpsParams::timeLimit));

  OsiSolverInterface * sSolver = 0;

  foundSolution = true;

  if (feasCheckSolver == "Cbc"){
    sSolver = new OsiCbcSolverInterface();
  }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
    sSolver = new OsiSymSolverInterface();
#else
    throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
		    "solveSubproblem", "MibSHeuristic");
#endif
  }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
    sSolver = new OsiCpxSolverInterface();
#else
    throw CoinError("CPLEX chosen as solver, but it has not been enabled",
		    "solveSubproblem", "MibSHeuristic");
#endif
  }else{
    throw CoinError("Unknown solver chosen",
		    "solveSubproblem", "MibSHeuristic");
  }

  int i(0);
  int index(0);
  double remainingTime(0.0);
  double upperObjVal(0.0), lowerObjVal(0.0); 
  double etol(etol_);
  int tCols(oSolver->getNumCols());
  double uObjSense(oSolver->getObjSense());
  double lObjSense(model->getLowerObjSense());
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int *lColIndices = model->getLowerColInd();
  int *uColIndices = model->getUpperColInd();
  const double *uObjCoeffs = oSolver->getObjCoefficients();
  double *lObjCoeffs = model->getLowerObjCoeffs();

  double *nObjCoeffs = new double[tCols];

  double *colSol = new double[tCols]();

  sSolver->loadProblem(*oSolver->getMatrixByCol(),
		       oSolver->getColLower(), oSolver->getColUpper(),
		       uObjCoeffs, oSolver->getRowLower(),
		       oSolver->getRowUpper());

  for(i = 0; i < tCols; i++){
    if(oSolver->isInteger(i))
      sSolver->setInteger(i);
  }

  //Multiply all columns of the UL objective by beta
  for(i = 0; i < tCols; i++){
    if(fabs(uObjCoeffs[i]) > etol){
      nObjCoeffs[i] = beta * uObjCoeffs[i] * uObjSense;
    }
    else{
      nObjCoeffs[i] = 0.0;
    }
  }

  //Add the LL columns of the LL objective multiplied by (1 - beta)
  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    if(fabs(lObjCoeffs[i]) > etol){
      nObjCoeffs[index] += (1 - beta) * lObjCoeffs[i] * lObjSense;
    }
  }

  sSolver->setObjective(nObjCoeffs);

  remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
  remainingTime = CoinMax(remainingTime, 0.00);

  if (feasCheckSolver == "Cbc"){
        dynamic_cast<OsiCbcSolverInterface *>
	  (sSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
        sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (sSolver)->getSymphonyEnvironment();
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
    sSolver->setHintParam(OsiDoReducePrint);
    sSolver->messageHandler()->setLogLevel(0);
    CPXENVptr cpxEnv =
      dynamic_cast<OsiCpxSolverInterface*>(sSolver)->getEnvironmentPtr();
    assert(cpxEnv);
    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
  }

  sSolver->branchAndBound();

  if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					 (dynamic_cast<OsiSymSolverInterface *>
					  (sSolver)->getSymphonyEnvironment()))){
    foundSolution = false;
    //delete [] colSol;
    goto TERM_SOLVESUBPROBLEM;
  }

  if(sSolver->isProvenOptimal()){
    CoinCopyN(sSolver->getColSolution(), tCols, colSol);
    for(i = 0; i < tCols; i++){
	if(sSolver->isInteger(i)){
	    colSol[i] = (double) floor(colSol[i] + 0.5);
	}
	upperObjVal += colSol[i] * uObjCoeffs[i];
    }
    lowerObjVal = getLowerObj(colSol, lObjSense);

    CoinPackedVector objCon;

    if(fabs(beta - 1.0) < etol){
      //modify sSolver to fix upper-level objective to current value and
      //reoptimize wrt to lower-level objective

      CoinZeroN(nObjCoeffs, tCols);
      for(i = 0; i < lCols; i++){
	index = lColIndices[i];
        nObjCoeffs[index] = lObjCoeffs[i] * lObjSense;
      }

      sSolver->setObjective(nObjCoeffs);

      for(i = 0; i < tCols; i++){
	objCon.insert(i, uObjCoeffs[i] * uObjSense);
      }

      sSolver->addRow(objCon, upperObjVal, upperObjVal);
      objCon.clear();

      remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
      remainingTime = CoinMax(remainingTime, 0.00);
      if(feasCheckSolver == "SYMPHONY"){
	sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	  (sSolver)->getSymphonyEnvironment();
	sym_set_dbl_param(env, "time_limit", remainingTime);
      }

      sSolver->branchAndBound();

      if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					     (dynamic_cast<OsiSymSolverInterface *>
					      (sSolver)->getSymphonyEnvironment()))){
	foundSolution = false;
	//delete [] colSol;
	goto TERM_SOLVESUBPROBLEM;
      }

      if(sSolver->isProvenOptimal()){
	lowerObjVal = sSolver->getObjValue();
	CoinCopyN(sSolver->getColSolution(), tCols, colSol);
      }
    }
    else if(fabs(beta) < etol){
      //modify sSolver to fix lower-level objective to current value and
      //reoptimize wrt to upper-level objective

      CoinZeroN(nObjCoeffs, tCols);
      for(i = 0; i < tCols; i++){
	nObjCoeffs[i] = uObjCoeffs[i] * uObjSense;
      }
      sSolver->setObjective(nObjCoeffs);

      for(i = 0; i < lCols; i++){
	index = lColIndices[i];
	objCon.insert(index, lObjCoeffs[i] * lObjSense);
      }
      sSolver->addRow(objCon, lowerObjVal, lowerObjVal);

      remainingTime = timeLimit - model->broker_->subTreeTimer().getTime();
      remainingTime = CoinMax(remainingTime, 0.00);
      if(feasCheckSolver == "SYMPHONY"){
	        sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		  (sSolver)->getSymphonyEnvironment();
		sym_set_dbl_param(env, "time_limit", remainingTime);
      }

      sSolver->branchAndBound();

      if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
					     (dynamic_cast<OsiSymSolverInterface *>
					      (sSolver)->getSymphonyEnvironment()))){
	foundSolution = false;
	//delete [] colSol;
	goto TERM_SOLVESUBPROBLEM;
      }

      if(sSolver->isProvenOptimal()){
	upperObjVal = sSolver->getObjValue();
	CoinCopyN(sSolver->getColSolution(), tCols, colSol);
      }
    }
  }
  else{
    //delete [] colSol;
    foundSolution = false;
  }

 TERM_SOLVESUBPROBLEM:
  //if(foundSolution == true){
  //mcSol foundSol = mcSol(std::make_pair(upperObjVal, lowerObjVal), colSol, tCols);
    //foundSol.setObjPair(std::make_pair(upperObjVal, lowerObjVal));
    //foundSol.setColSol(colSol);
    //foundSol.setLen(tCols);
    //}
  
  delete sSolver;
  delete [] nObjCoeffs;
  //delete [] colSol;

  return mcSol(std::make_pair(upperObjVal, lowerObjVal), colSol, tCols);
}

//#############################################################################
double
MibSHeuristic::getLowerObj(const double * sol, double objSense)
{

  int lCols(MibSModel_->getLowerDim());
   int * lColIndices = MibSModel_->getLowerColInd();
   double * lObjCoeffs = MibSModel_->getLowerObjCoeffs();

   int i(0), index(0);
   double objVal(0.0);



   for(i = 0; i < lCols; i++){
      index = lColIndices[i];
      if(0){
	std::cout << "obj x sol: " << lObjCoeffs[i] 
		  << " x " << sol[index] << std::endl;
      }
      objVal += lObjCoeffs[i] * sol[index];
   }

   return (objVal * objSense);
   
}

//#############################################################################
void
MibSHeuristic::addSolutionToSeenLinkingSolutionPoolHeur(std::vector<double> &linkSol,
						  std::vector<double> &shouldStoreValues,
						  double objValue){
  MibSModel * model = MibSModel_;
  
  LINKING_SOLUTION linkingSolution;

  linkingSolution.tag = MibSLinkingPoolTagLowerIsFeasible;
  linkingSolution.lowerObjValue = objValue;
  linkingSolution.UBObjValue = 0.0;
  linkingSolution.lowerSolution = shouldStoreValues;
  linkingSolution.UBSolution.push_back(0);
  model->seenLinkingSolutions[linkSol] = linkingSolution;

}

  
