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

  bool useLowerObjHeuristic
    = model->MibSPar_->entry(MibSParams::useLowerObjHeuristic);

  bool useObjCutHeuristic
    = model->MibSPar_->entry(MibSParams::useObjCutHeuristic);

  bool useWSHeuristic
    = model->MibSPar_->entry(MibSParams::useWSHeuristic);

  bool useGreedyHeuristic
    = model->MibSPar_->entry(MibSParams::useGreedyHeuristic);

  if(useLowerObjHeuristic)
    lowerObjHeuristic();
  
  if(useObjCutHeuristic)
    objCutHeuristic();

  if(useWSHeuristic)
    weightedSumsHeuristic();

  if(useGreedyHeuristic)
    greedyHeuristic();

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

  /* Solve the LP relaxation with the new constraint d^2 y <= d^y* */

  MibSModel * model = MibSModel_;

  //OsiSolverInterface * oSolver = model->origLpSolver_;
  OsiSolverInterface * oSolver = model->getSolver();
#ifndef COIN_HAS_SYMPHONY
  OsiSolverInterface * hSolver = new OsiCbcSolverInterface();
#else
  OsiSolverInterface * hSolver = new OsiSymSolverInterface();
#endif
  
  double objSense(model->getLowerObjSense());  
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int tCols(lCols + uCols);
  int * lColIndices = model->getLowerColInd();
  int * uColIndices = model->getUpperColInd();
  double * lObjCoeffs = model->getLowerObjCoeffs();

  hSolver->loadProblem(*oSolver->getMatrixByCol(),
		       oSolver->getColLower(), oSolver->getColUpper(),
		       oSolver->getObjCoefficients(),
		       oSolver->getRowLower(), oSolver->getRowUpper());

  int j(0);
  for(j = 0; j < tCols; j++){
    if(oSolver->isInteger(j))
      hSolver->setInteger(j);
  }

  double * optLowerSolutionOrd = model->bS_->optLowerSolutionOrd_;

  CoinPackedVector objCon;
  int i(0), index(0);
  double rhs(0.0);

  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    objCon.insert(index, lObjCoeffs[i] * objSense);
    //should this be ordered? and should lObjCoeffs by at index?
    //rhs += optLowerSolutionOrd_[i] * lObjCoeffs[i] * objSense;
    rhs += optLowerSolutionOrd[i] * lObjCoeffs[i] * objSense;
  }

  //Hmm, I think this was wrong before...?
  //  hSolver->addRow(objCon, - hSolver->getInfinity(), rhs);
  hSolver->addRow(objCon, rhs, hSolver->getInfinity());

  /* optimize w.r.t. to the UL objective with the new row */
#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiCbcSolverInterface *> 
     (hSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
    dynamic_cast<OsiSymSolverInterface *> 
      (hSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (hSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (hSolver)->setSymParam("max_active_nodes", 1);
#endif

  hSolver->branchAndBound();

  if(0)
    hSolver->writeLp("objcutheuristic");

  if(hSolver->isProvenOptimal()){

    MibSSolution *mibSol = NULL;

    OsiSolverInterface * lSolver = model->bS_->setUpModel(hSolver, true);

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

    const double * sol = hSolver->getColSolution();
    double objVal(lSolver->getObjValue() * objSense);
    double etol(etol_);
    double lowerObj = getLowerObj(sol, objSense);  
    
    double * optUpperSolutionOrd = new double[uCols];
    double * optLowerSolutionOrd = new double[lCols];
    
    CoinZeroN(optUpperSolutionOrd, uCols);
    CoinZeroN(optLowerSolutionOrd, lCols);

    if(fabs(objVal - lowerObj) < etol){
    
      /** Current solution is bilevel feasible **/
     
      mibSol = new MibSSolution(hSolver->getNumCols(),
				hSolver->getColSolution(),
				hSolver->getObjValue(),
				model);

     model->storeSolution(BlisSolutionTypeHeuristic, mibSol);
     mibSol = NULL;
 
    }
    else{

      /* solution is not bilevel feasible, create one that is */

     const double * uSol = hSolver->getColSolution();
     const double * lSol = lSolver->getColSolution();
     //int numElements(lSolver->getNumCols());
     int numElements(hSolver->getNumCols());
     int i(0), pos(0), index(0);
     double * lpSolution = new double[numElements];
     double upperObj(0.0);

     //FIXME: problem is still here.  indices may be wrong.  
     //also is all this necessary, or can we just paste together uSol and lSol?
     
     for(i = 0; i < numElements; i++){
       //index = indices[i];
       pos = model->bS_->binarySearch(0, lCols - 1, i, lColIndices);
       if(pos < 0){
	 pos = model->bS_->binarySearch(0, uCols - 1, i, uColIndices);
	 //optUpperSolutionOrd[pos] = values[i];
	 //optUpperSolutionOrd[pos] = uSol[pos];
	 if (pos >= 0){
	    optUpperSolutionOrd[pos] = uSol[i];
	 }
       }
       else{
	 //optLowerSolutionOrd[pos] = lSol[i];
	 optLowerSolutionOrd[pos] = lSol[pos];
       }
     }

     for(i = 0; i < uCols; i++){
       index = uColIndices[i];
       lpSolution[index] = optUpperSolutionOrd[i];
       upperObj += 
	 optUpperSolutionOrd[i] * hSolver->getObjCoefficients()[index];
     }

     for(i = 0; i < lCols; i++){
       index = lColIndices[i];
       lpSolution[index] = optLowerSolutionOrd[i];
       upperObj += 
	 optLowerSolutionOrd[i] * hSolver->getObjCoefficients()[index];
     }

     if(model->checkUpperFeasibility(lpSolution)){
       mibSol = new MibSSolution(hSolver->getNumCols(),
				 lpSolution,
				 upperObj * hSolver->getObjSense(),
				 model);

       model->storeSolution(BlisSolutionTypeHeuristic, mibSol);
       mibSol = NULL;
     }
     delete [] lpSolution;
    }
    delete lSolver;
  }

  delete hSolver;

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
  
  /*
    now we find a feasible solution by fixing upper-level vars
    and solving the lower-level problem
  */
  
  double * incumbentSol = new double[tCols];
  double * colsol = new double[tCols];

  CoinZeroN(colsol, tCols);

  for(i = 0; i < uCols; i++){
    colsol[uColIndices[i]] = fixedVars[i];
    if(fixedVars[i] == 1.0)
      if(0)
	std::cout << "fixed " << i << std::endl;
  }

  bfSol * sol = getBilevelSolution(colsol, lObjSense * oSolver->getInfinity());

  if(sol){
    double incumbentObjVal = sol->getObjVal();
    CoinCopyN(sol->getColumnSol(), tCols, incumbentSol);
    
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
  mcSol sol0 = mcSol();
  mcSol sol1 = mcSol();
  mcSol sol = mcSol();
  double beta(0.0);
  double slope(0.0);
  double etol(etol_);   

  BETAS.push_back(1.0);
  BETAS.push_back(0.0);

  /** First solve the subproblem for beta = 1 and beta = 0 **/
  
  sol1 = solveSubproblem(1.0); // beta = 1
  numProblems_++;
  mcSolutions.insert(std::make_pair(1.0, sol1));

  int i(0);
  int tCols(MibSModel_->getLowerDim() + MibSModel_->getUpperDim());
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
  
  sol0 = solveSubproblem(0.0);// beta = 0
  numProblems_++;
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
	sol = solveSubproblem(beta);
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
  double uObjSense(model->getSolver()->getObjSense());
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());

  int tCols(uCols + lCols); 

  double incumbentObjVal(model->getSolver()->getInfinity() * uObjSense);
  double * incumbentSol = new double[tCols];

  int i(0);

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
    if(0){
      std::cout << "Candidate solution value: " << origLower << std::endl;
      for(i = 0; i < tCols; i++){
	std::cout << "colsol[" << i << "] :" 
		  << colsol[i] << std::endl;
      }
    }
    
    bfSol * sol = getBilevelSolution(colsol, origLower);

    if(sol){
      

      if(0){
	std::cout << "Returned solution: " << std::endl;
	for(i = 0; i < tCols; i++){
	  std::cout << "sol->getColumnSol[" << i << "]" 
		    << sol->getColumnSol()[i] << std::endl;
	}
	std::cout << "sol->getObjVal: " << sol->getObjVal() << std::endl;
      }
      
      if(sol->getObjVal() < incumbentObjVal){

	incumbentObjVal = sol->getObjVal();
	CoinCopyN(sol->getColumnSol(), tCols, incumbentSol);

	if(0){
	  std::cout << "New incumbent found." << std::endl;
	  for(i = 0; i < tCols; i++){
	    std::cout << "incumbentSol[" << i 
		      << "]: " << incumbentSol[i] << std::endl;
	  }
	  std::cout << "incumbentObjVal: " << incumbentObjVal << std::endl;
	}
	
      }
    }

  }


 
  if(0){
    std::cout << "This solution comes from MibSHeuristic.cpp:742" << std::endl;
  }
  MibSSolution * mibSol = new MibSSolution(tCols,
					   incumbentSol,
					   incumbentObjVal,
					   model);
  
  model->storeSolution(BlisSolutionTypeHeuristic, mibSol);

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
bfSol*
MibSHeuristic::getBilevelSolution(const double * sol, double origLower)
{

  /* 
     Find a bilevel feasible solution by solving the LL problem
     for a fixed UL solution, given by the UL portion of sol
  */

  MibSModel * model = MibSModel_;
  OsiSolverInterface * oSolver = model->getSolver();
  OsiSolverInterface * lSolver = model->bS_->setUpModel(oSolver, true, sol);
  
  //double uObjSense(model->getSolver()->getObjSense());
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int * lColIndices = model->getLowerColInd();
  int * uColIndices = model->getUpperColInd();
  double etol(etol_);

  int tCols(uCols + lCols); 
  int i(0);

  if(0){
     lSolver->writeLp("bilevelsolver");
     std::cout 
       << "Original Lower-level Solution Value: " << origLower << std::endl;
     for(i = 0; i < lCols; i++){
       std::cout << "lowsol[" << i << "]: " << sol[lColIndices[i]] << std::endl;
     }
  }
  
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
    double lowerObj(lSolver->getObjValue());
    double * colsol = new double[tCols];

    for(i = 0; i < uCols; i++){
      colsol[uColIndices[i]] = sol[uColIndices[i]];
    }

    if(0){
      std::cout << "candidate lower solution value: " << origLower << std::endl;
      std::cout << "actual lower solution value: " << lowerObj << std::endl;
    }

    if(fabs(origLower - lowerObj) < etol){
      //original solution was bilevel feasible
      if(0)
	std::cout << "Original solution was bilevel feasible:" << std::endl;
      for(i = 0; i < lCols; i++){
	if(0){
	  std::cout << "lowerportion[" 
		    << i << "]: " << sol[lColIndices[i]] << std::endl; 
	}	
	colsol[lColIndices[i]] = sol[lColIndices[i]];
      }
    }
    else{
      if(0){
	std::cout << "Not bilevel feasible." << std::endl;
      }
      for(i = 0; i < lCols; i++){
	if(0){
	  std::cout << "newportion[" 
		    << i << "]: " << lSolver->getColSolution()[i] << std::endl; 
	}
	colsol[lColIndices[i]] = lSolver->getColSolution()[i];
      }
    }


    for(i = 0; i < tCols; i++)
      objVal += colsol[i] * oSolver->getObjCoefficients()[i];

    bfSol * bfsol = 
      new bfSol(objVal, colsol);
    delete lSolver;
    return bfsol;
  }
  else{
    delete lSolver;
    return NULL;
  }

}

//#############################################################################
bfSol*
MibSHeuristic::getBilevelSolution1(const double * sol)
{

  /* 
     Find a bilevel feasible solution by solving the LL problem
     for a fixed UL solution, given by the UL portion of sol
  */


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

  /* delete the UL rows */
  lSolver->deleteRows(uRowNum, uRowIndices);

  /* Fix the UL variables to their current value in sol */

  for(i = 0; i < uCols; i++){
    index = uColIndices[i];
    lSolver->setColLower(index, sol[index]);
    lSolver->setColUpper(index, sol[index]);
  }

  /* Set the objective to the LL objective coefficients */

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
 
    bfSol * bfsol = 
      new bfSol(objVal, colsol);
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
MibSHeuristic::solveSubproblem(double beta)
{

  /* 
     optimize wrt to weighted upper-level objective 
     over current feasible lp feasible region 
  */

  MibSModel * model = MibSModel_;
  OsiSolverInterface * oSolver = model->getSolver();
#ifndef COIN_HAS_SYMPHONY
  OsiSolverInterface * sSolver = new OsiCbcSolverInterface();
#else
  OsiSolverInterface* sSolver = new OsiSymSolverInterface();
  //sSolver = oSolver->clone();
  //OsiSolverInterface * sSolver = tmpSolver;
  //OsiSolverInterface * tmpSolver = new OsiSolverInterface(oSolver);
#endif
  
  double uObjSense(oSolver->getObjSense());
  double lObjSense(model->getLowerObjSense());  
  int lCols(model->getLowerDim());
  int uCols(model->getUpperDim());
  int * lColIndices = model->getLowerColInd();
  int * uColIndices = model->getUpperColInd();
  double * lObjCoeffs = model->getLowerObjCoeffs();
  const double * uObjCoeffs = oSolver->getObjCoefficients();

  double etol(etol_);
  //int tCols(uCols + lCols);
  int tCols(oSolver->getNumCols());

  //assert(tCols == oSolver->getNumCols());


  sSolver->loadProblem(*oSolver->getMatrixByCol(),
		       oSolver->getColLower(), oSolver->getColUpper(),
		       oSolver->getObjCoefficients(),
		       oSolver->getRowLower(), oSolver->getRowUpper());

  int j(0);
  for(j = 0; j < tCols; j++){
    if(oSolver->isInteger(j))
      sSolver->setInteger(j);
  }


  double * nObjCoeffs = new double[tCols];
  int i(0), index(0);
  
  CoinZeroN(nObjCoeffs, tCols);
  
  /* Multiply the UL columns of the UL objective by beta */
  for(i = 0; i < uCols; i++){
    index = uColIndices[i];
    if(fabs(uObjCoeffs[index]) > etol)
      nObjCoeffs[index] = beta * uObjCoeffs[index] * uObjSense;
    else 
      nObjCoeffs[index] = 0.0;
  }
    
  /* Multiply the LL columns of the UL objective by beta */
  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    if(fabs(uObjCoeffs[index]) > etol)
      nObjCoeffs[index] = beta* uObjCoeffs[index] * uObjSense;
    else
      nObjCoeffs[index] = 0.0;
  }
  
  /* Add the LL columns of the LL objective multiplied by (1 - beta) */
  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    if(fabs(lObjCoeffs[i]) > etol)
      nObjCoeffs[index] += (1 - beta) * lObjCoeffs[i] * lObjSense;
  }
  
  sSolver->setObjective(nObjCoeffs);

  //int i(0);
  if(0){
    for(i = 0; i < sSolver->getNumCols(); i++){
      std::cout << "betaobj " << sSolver->getObjCoefficients()[i] << std::endl;
    }
  }

  if(0){
     sSolver->writeLp("afterbeta");
     //sSolver->writeMps("afterbeta");
  }
  
  if(0){  
    for(i = 0; i < sSolver->getNumCols(); i++){
      std::cout << "obj " << sSolver->getObjCoefficients()[i] << std::endl;
      std::cout << "upper " << sSolver->getColUpper()[i] << std::endl;
      std::cout << "lower " << sSolver->getColLower()[i] << std::endl;
    }
  }

#ifndef COIN_HAS_SYMPHONY
  dynamic_cast<OsiCbcSolverInterface *> 
     (sSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
  dynamic_cast<OsiSymSolverInterface *> 
     (sSolver)->setSymParam("prep_level", -1);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (sSolver)->setSymParam("verbosity", -2);
  
  dynamic_cast<OsiSymSolverInterface *> 
     (sSolver)->setSymParam("max_active_nodes", 1);
#endif

  //dynamic_cast<OsiSymSolverInterface *> (sSolver)->branchAndBound();

  sSolver->branchAndBound();

  if(sSolver->isProvenOptimal()){

    if(0){
      std::cout << "writing lp file." << std::endl;
      sSolver->writeLp("afterbeta");
      //sSolver->writeMps("afterbeta");
    }
    
    double upperObjVal(0.0);
    double lowerObjVal(0.0);
    

    for(i = 0; i < tCols; i++){
      upperObjVal += 
	sSolver->getColSolution()[i] * oSolver->getObjCoefficients()[i];
      if(0){
	std::cout << "sSolver->getColSolution()[" << i << "] :"
		  << sSolver->getColSolution()[i] << std::endl;
      }
    }
    lowerObjVal = getLowerObj(sSolver->getColSolution(), lObjSense);
    
    if(beta == 1.0){
      
      /*
	fix upper-level objective to current value and 
	reoptimize wrt to lower-level objective
      */
      
#ifndef COIN_HAS_SYMPHONY
       OsiSolverInterface * nSolver = new OsiCbcSolverInterface();
#else
       OsiSolverInterface * nSolver = new OsiSymSolverInterface();
#endif

       nSolver->loadProblem(*oSolver->getMatrixByCol(),
			    oSolver->getColLower(), oSolver->getColUpper(),
			    oSolver->getObjCoefficients(),
			    oSolver->getRowLower(), oSolver->getRowUpper());
       for(j = 0; j < tCols; j++){
	  if(oSolver->isInteger(j))
	     nSolver->setInteger(j);
       }
       

      CoinZeroN(nObjCoeffs, tCols);
      
      for(i = 0; i < lCols; i++){
	index = lColIndices[i];
	nObjCoeffs[index] = lObjCoeffs[i] * lObjSense;
      }
      
      nSolver->setObjective(nObjCoeffs);
      
      CoinPackedVector objCon;
      
      for(i = 0; i < tCols; i++){
	objCon.insert(i, uObjCoeffs[i] * uObjSense);
      }
      
      nSolver->addRow(objCon, upperObjVal, upperObjVal);
      nSolver->writeLp("beta1");

#ifndef COIN_HAS_SYMPHONY
      dynamic_cast<OsiCbcSolverInterface *> 
	 (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
      dynamic_cast<OsiSymSolverInterface *> 
	 (nSolver)->setSymParam("prep_level", -1);
      
      dynamic_cast<OsiSymSolverInterface *> 
	 (nSolver)->setSymParam("verbosity", -2);
      
      dynamic_cast<OsiSymSolverInterface *> 
	 (nSolver)->setSymParam("max_active_nodes", 1);
#endif

      nSolver->branchAndBound();
     

      double * colsol = new double[tCols];

      if(nSolver->isProvenOptimal()){
	lowerObjVal = nSolver->getObjValue();
	CoinCopyN(nSolver->getColSolution(), tCols, colsol);
      }
      else{
	//just take the current solution
	lowerObjVal = sSolver->getObjValue();
	CoinCopyN(sSolver->getColSolution(), tCols, colsol);
      }

      delete[] nObjCoeffs;
      nObjCoeffs = 0;
      delete sSolver;
      delete nSolver;
      return mcSol(std::make_pair(upperObjVal, lowerObjVal), colsol);
    }
    else if(beta == 0.0){
      
      /*
	fix lower-level objective to current value and 
	reoptimize wrt to upper-level objective
      */
      
#ifndef COIN_HAS_SYMPHONY
      OsiSolverInterface * nSolver = new OsiCbcSolverInterface();
#else
      OsiSolverInterface * nSolver = new OsiSymSolverInterface();
#endif

      nSolver->loadProblem(*oSolver->getMatrixByCol(),
			   oSolver->getColLower(), oSolver->getColUpper(),
			   oSolver->getObjCoefficients(),
			   oSolver->getRowLower(), oSolver->getRowUpper());
      for(j = 0; j < tCols; j++){
	if(oSolver->isInteger(j))
	  nSolver->setInteger(j);
      }
      
      CoinZeroN(nObjCoeffs, tCols);
	
      for(i = 0; i < tCols; i++)
	nObjCoeffs[i] = uObjCoeffs[i] * uObjSense;
      
      nSolver->setObjective(nObjCoeffs);
	
      CoinPackedVector objCon;
	
      for(i = 0; i < lCols; i++){
	index = lColIndices[i];
	objCon.insert(index, lObjCoeffs[i] * lObjSense);  
      }
      
      nSolver->addRow(objCon, lowerObjVal, lowerObjVal);
      
#ifndef COIN_HAS_SYMPHONY
      dynamic_cast<OsiCbcSolverInterface *> 
	 (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
#else
      dynamic_cast<OsiSymSolverInterface *> 
	 (nSolver)->setSymParam("prep_level", -1);
      
      dynamic_cast<OsiSymSolverInterface *> 
	 (nSolver)->setSymParam("verbosity", -2);
      
      dynamic_cast<OsiSymSolverInterface *> 
	 (nSolver)->setSymParam("max_active_nodes", 1);
#endif

      if(0)      
	nSolver->writeLp("nSolver");
      

      nSolver->branchAndBound();
	
      double * colsol = new double[tCols];
	
      if(nSolver->isProvenOptimal()){
	upperObjVal = nSolver->getObjValue();
	CoinCopyN(nSolver->getColSolution(), tCols, colsol);
      }
      else{
	upperObjVal = nSolver->getObjValue();
	CoinCopyN(nSolver->getColSolution(), tCols, colsol);
      }

      delete[] nObjCoeffs;
      nObjCoeffs = 0;
      delete sSolver;
      delete nSolver;
      return mcSol(std::make_pair(upperObjVal, lowerObjVal), colsol);
	
    }
    else{
      
      //no optimality cut needed here.  all solutions are supported.
      
      double * colsol = new double[tCols];
      CoinCopyN(sSolver->getColSolution(), tCols, colsol);	
      
      delete[] nObjCoeffs;
      nObjCoeffs = 0;
      delete sSolver;
      return mcSol(std::make_pair(upperObjVal, lowerObjVal), colsol);
      
    }
    
  }
  else{
    //FIXME:SHOULD JUST TAKE THIS OUT.  DELETE sSolver and remove it from above
    
    nObjCoeffs = 0;
    delete[] nObjCoeffs;
    delete sSolver;
    std::cout << "Subproblem is not proven optimal." << std::endl;
    //return NULL;
    //abort();
  }
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

  
