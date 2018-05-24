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
#include "MibSConstants.hpp"
#include "MibSHelp.hpp"

class MibSModel;
class MibSCutGenerator;
class MibSHeuristic;
class MibSTreeNode;

//#############################################################################

class MibSBilevel {

    friend class MibSModel;
    friend class MibSCutGenerator;
    friend class MibSBranchStrategyMaxInf;
    friend class MibSHeuristic;
    friend class MibSTreeNode;
    friend class MibSBranchStrategyPseudo;
    friend class MibSBranchStrategyStrong;

private:

    bool isIntegral_;
    bool isUpperIntegral_;
    bool isLinkVarsIntegral_;
    bool useBilevelBranching_;
    bool isLinkVarsFixed_;
    bool isProvenOptimal_;
    /** is lower-level problem solved or not **/
    bool isLowerSolved_;
    /** is problem UB solved or not **/
    bool isUBSolved_;
    /** should prune the node or not **/
    bool shouldPrune_;
    bool isContainedInLinkingPool_;
    MibSLinkingPoolTag tagInSeenLinkingPool_;

    MibSLPSolStatus LPSolStatus_;

    /** Optimal value of LL objective **/
    double objVal_;
    int linkIntegralCount_;

    double *upperSolutionOrd_;
    double *lowerSolutionOrd_;
    //double *optLowerSolution_;
    double *optUpperSolutionOrd_;// result of solving (UB)
    double *optLowerSolutionOrd_;
   
    MibSModel *model_;
    MibSHeuristic *heuristic_;
    OsiSolverInterface * lSolver_;
    OsiSolverInterface * UBSolver_;
    CoinWarmStart * ws_;
   
public:
   
    MibSBilevel() : isIntegral_(true), isUpperIntegral_(true),
		    isLinkVarsIntegral_(true), useBilevelBranching_(true),
		    isLinkVarsFixed_(true), isProvenOptimal_(false),
		    isLowerSolved_(false), isUBSolved_(false),
		    shouldPrune_(false), isContainedInLinkingPool_(false),
		    tagInSeenLinkingPool_(MibSLinkingPoolTagIsNotSet),
		    LPSolStatus_(MibSLPSolStatusUnknown), objVal_(0.0),
		    linkIntegralCount_(0){
	upperSolutionOrd_ = 0;
	lowerSolutionOrd_ = 0;
	//optLowerSolution_ = 0;
	optUpperSolutionOrd_ = 0;
	optLowerSolutionOrd_ = 0;
	model_ = 0;
	heuristic_= 0;
	lSolver_ = 0;
	UBSolver_ = 0;
	ws_ = 0;
    }
   
    ~MibSBilevel() {
	gutsOfDestructor();
    }
   
    MibSSolType createBilevel(CoinPackedVector *sol,
		       MibSModel *mibs=0);
    MibSSolType checkBilevelFeasiblity(bool isRoot);
    void gutsOfDestructor();

private:
   
    int findIndex(int index, int size, int * indices);
    OsiSolverInterface * setUpUBModel(OsiSolverInterface * solver, double objValLL,
					  bool newOsi, const double *sol = NULL);
    OsiSolverInterface * setUpModel(OsiSolverInterface * solver,
				    bool newOsi, const double *sol = NULL);
    double getLowerObj(const double * sol, double objSense);
    int binarySearch(int index,int start, int stop, int * indexArray);
    CoinWarmStart * getWarmStart() {return ws_;}
    void setWarmStart(CoinWarmStart * ws) {ws_ = ws;}
    void addSolutionToSeenLinkingSolutionPool(MibSLinkingPoolTag solTag, std::vector<double>
		      &shouldStoreValues, double objValue);
    //void findHeuristicSolutions();
    //void objCutHeuristic();
    //void lowerObjHeuristic();
    //void weightedSumsHeuristic();
    //void greedyHeuristic();

};

#endif
