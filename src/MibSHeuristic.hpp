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

#ifndef MibSHeuristic_h_
#define MibSHeuristic_h_

#include "MibSModel.hpp"
#include "MibSSolTypes.hpp"

#include <map>

class MibSModel;

//#############################################################################

class MibSHeuristic {

 private:

  MibSModel * MibSModel_;  
  std::map<double, mcSol> mcSolutions_;
  int numProblems_;
  double etol_;
  double bestObjVal_;
  double * bestSol_;
  
 public: 
  MibSHeuristic();
  MibSHeuristic(MibSModel * model);
  ~MibSHeuristic();

  /** Initialize the data for this class **/
  void initializeData(MibSModel * model);

  /** Simple function that calls the correct solver **/
  void findHeuristicSolutions();
  
  /** Run the heuristic solver using weighted sums algorithm **/
  void weightedSumsHeuristic();

  /** Run the heuristic solver using greedy algorithm for interdiction**/
  void greedyHeuristic();

  void objCutHeuristic();
  void lowerObjHeuristic();

  /** Solve the mip subproblem for a given value of beta **/
  mcSol solveSubproblem(double beta, bool &foundSolution);
  mcSol solveSubproblem1(double beta);

  void createBilevelSolutions(std::map<double, mcSol>);

  bool checkUpperFeasibility(double * solution);

  bool checkLowerFeasibility1(double * solution);

  bool checkLowerFeasibility(OsiSolverInterface * si, double * solution);

  double getLowerObj(const double * sol, double objSense);

  void addSolutionToSeenLinkingSolutionPoolHeur(std::vector<double> &linkSol,
						std::vector<double> &shouldStoreValues,
						double objValue);

  bfSol getBilevelSolution(const double * sol, double origLower,
			   bool &isTimeLimReached, bool &foundFeasible);
  bfSol * getBilevelSolution1(const double * sol);


};


#endif
