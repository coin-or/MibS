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

#ifndef MibSBilevel_h_
#define MibSBilevel_h_

#include "OsiSolverInterface.hpp"
#include "CoinWarmStart.hpp"

#include "MibSModel.hpp"
#include "MibSHeuristic.hpp"

class MibSModel;
class MibSCutGenerator;
class MibSHeuristic;

//#############################################################################

class MibSBilevel {

    friend class MibSModel;
    friend class MibSCutGenerator;
    friend class MibSBranchStrategyMaxInf;
    friend class MibSHeuristic;

private:

    bool isIntegral_;
    bool isBilevelFeasible_;
    bool isUpperIntegral_;
    bool useBilevelBranching_;
    bool upperFixed_;
    bool isProvenOptimal_;

    /** Optimal value of LL objective **/
    double objVal_;
    
    double *upperSolution_;
    double *lowerSolution_;
    double *upperSolutionOrd_;
    double *lowerSolutionOrd_;
    double *optLowerSolution_;
    double *optLowerSolutionOrd_;
   
    MibSModel *model_;
    MibSHeuristic *heuristic_;
    OsiSolverInterface * solver_;
    OsiSolverInterface * solver2_;  
    CoinWarmStart * ws_;

   
public:
   
    MibSBilevel() : isIntegral_(false), isBilevelFeasible_(false),
		    isUpperIntegral_(false), useBilevelBranching_(false),
		    objVal_(0.0){
	upperSolution_ = 0;
	lowerSolution_ = 0;
	upperSolutionOrd_ = 0;
	lowerSolutionOrd_ = 0;
	optLowerSolution_ = 0;
	optLowerSolutionOrd_ = 0;
	model_ = 0;
	heuristic_= 0;
	solver_ = 0;
	solver2_ = 0;
	ws_ = 0;
    }
   
    ~MibSBilevel() {
	gutsOfDestructor();
    }
   
    void createBilevel(CoinPackedVector *sol,
		       MibSModel *mibs=0);
    void checkBilevelFeasiblity(bool isRoot, bool upperFixed_);
    void gutsOfDestructor();

private:
   
    int findIndex(int index, int size, int * indices);
    OsiSolverInterface * setUpRefineModel(OsiSolverInterface * solver, double objValLL,
					  bool newOsi, const double *sol = NULL);
    OsiSolverInterface * setUpModel(OsiSolverInterface * solver, double objValLL,
				    bool bestSolLL, bool newOsi, const double *sol = NULL);
    double getLowerObj(const double * sol, double objSense);
    int binarySearch(int index,int start, int stop, int * indexArray);
    CoinWarmStart * getWarmStart() {return ws_;}
    void setWarmStart(CoinWarmStart * ws) {ws_ = ws;}
    //void findHeuristicSolutions();
    //void objCutHeuristic();
    //void lowerObjHeuristic();
    //void weightedSumsHeuristic();
    //void greedyHeuristic();

};

#endif
