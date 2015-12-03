/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*                                                                           */
/* Copyright (C) 2007-2015 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#include "MibSBranchStrategyPseudo.h"
#include "MibSModel.h"

#include "BlisObjectInt.h"
#include "BlisModel.h"
#include "BlisHelp.h"

//#############################################################################

struct BlisPseuoGreater
{
    bool operator()(double x, double y) const {
        return (x > y);
    }
};

//#############################################################################
MibSBranchStrategyPseudo::MibSBranchStrategyPseudo()
   : BlisBranchStrategyPseudo()
{


}

//#############################################################################
MibSBranchStrategyPseudo::MibSBranchStrategyPseudo(BlisModel *model, int rel)
   : BlisBranchStrategyPseudo(model, rel)
{

}

//#############################################################################
MibSBranchStrategyPseudo::~MibSBranchStrategyPseudo()
{

}

//#############################################################################

/** Create a set of candidate branching objects. */
int 
MibSBranchStrategyPseudo::createCandBranchObjects(int numPassesLeft, double ub)
{

   //FIXEME::ADD IN BILEVEL BRANCHING

    int bStatus = 0;
    int i, pass, colInd;

    int preferDir, saveLimit;
    int numFirsts  = 0;
    int numInfs = 0;
    int minCount = 0;
    int numLowerTightens = 0;
    int numUpperTightens = 0;
    double lpX, score, infeasibility, downDeg, upDeg, sumDeg = 0.0; 
    
    bool roundAgain, downKeep, downGood, upKeep, upGood;


    int *lbInd = NULL;
    int *ubInd = NULL;
    double *newLB = NULL;
    double *newUB = NULL;

    double *saveUpper = NULL;
    double *saveLower = NULL;
    double *saveSolution = NULL;

    BlisModel *model = dynamic_cast<BlisModel *>(model_);
    OsiSolverInterface *solver = model->solver();
    
    int numCols = model->getNumCols();
    int numObjects = model->numObjects();
    int aveIterations = model->getAveIterations();


    //std::cout <<  "aveIterations = " <<  aveIterations << std::endl;

     //------------------------------------------------------
    // Check if max time is reached or no pass is left.
    //------------------------------------------------------
    
    double timeLimit = model->AlpsPar()->entry(AlpsParams::timeLimit);
    AlpsKnowledgeBroker *broker = model->getKnowledgeBroker();
    bool maxTimeReached = (broker->timer().getTime() > timeLimit);
    bool selectNow = false;
    
    if (maxTimeReached || !numPassesLeft) {
        selectNow = true;
#ifdef BLIS_DEBUG
        printf("PSEUDO: CREATE: maxTimeReached %d, numPassesLeft %d\n", 
               maxTimeReached, numPassesLeft);
#endif
    }
    
    // Store first time objects.
    std::vector<BlisObjectInt *> firstObjects;

    // Store infeasible objects.
    std::vector<BlisObjectInt *> infObjects;

    // TODO: check if sorting is expensive.
    std::multimap<double, BcpsBranchObject*, BlisPseuoGreater> candObjects;

    double objValue = solver->getObjSense() * solver->getObjValue();

    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    saveSolution = new double[numCols];
    memcpy(saveSolution, solver->getColSolution(), numCols*sizeof(double));

    //--------------------------------------------------
    // Find the infeasible objects.
    // NOTE: we might go round this loop twice if we are feed in
    //       a "feasible" solution.
    //--------------------------------------------------
    
    for (pass = 0; pass < 2; ++pass) {
	
        numInfs = 0;

        BcpsObject * object = NULL;
        BlisObjectInt * intObject = NULL;
            
        infObjects.clear();
        firstObjects.clear();
        
        for (i = 0; i < numObjects; ++i) {
                
            object = model->objects(i);
            infeasibility = object->infeasibility(model, preferDir);
            
            if (infeasibility) {
                
                ++numInfs;
                intObject = dynamic_cast<BlisObjectInt *>(object);
                
                if (intObject) {
                    infObjects.push_back(intObject);
                    
                    if (!selectNow) {
                        minCount = 
                            ALPS_MIN(intObject->pseudocost().getDownCount(),
                                     intObject->pseudocost().getUpCount());
                        
                        if (minCount < 1) {
                            firstObjects.push_back(intObject);
                        }
                    }

#ifdef BLIS_DEBUG
                    if (intObject->columnIndex() == 40) {
                        std::cout << "x[40] = " << saveSolution[40] 
                                  << std::endl;
                    }
#endif

                    intObject = NULL;
                }
                else {
                    // TODO: currently all are integer objects.
#ifdef BLIS_DEBU
                    assert(0);
#endif
                }
                
            }
        }
            
        if (numInfs) {
#if 0
            std::cout << "PSEUDO: numInfs = " << numInfs
                      << std::endl;
#endif
            break;
        }
        else if (pass == 0) {
            // The first pass and is IP feasible.
            
#if 1
            std::cout << "ERROR: PSEUDO: given a integer feasible sol, no fraction variable" << std::endl;
            assert(0);
#endif      
            
            roundAgain = false;
            CoinWarmStartBasis * ws = 
                dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
            if (!ws) break;
            
            // Force solution values within bounds
            for (i = 0; i < numCols; ++i) {
                lpX = saveSolution[i];
                if (lpX < lower[i]) {
                    saveSolution[i] = lower[i];
                    roundAgain = true;
                    ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
                } 
                else if (lpX > upper[i]) {
                    saveSolution[i] = upper[i];
                    roundAgain = true;
                    ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
                } 
            }
            
            if (roundAgain) {
                // Need resolve and do the second round selection.
                solver->setWarmStart(ws);
                delete ws;
                
                // Resolve.
                solver->resolve();
		
                if (!solver->isProvenOptimal()) {
                    // Become infeasible, can do nothing. 
                    bStatus = -2;
                    goto TERM_CREATE;
                }
                else {
                    // Save new lp solution.
                    memcpy(saveSolution, solver->getColSolution(),
                           numCols * sizeof(double));
                    objValue = solver->getObjSense() * solver->getObjValue();
                }
            } 
            else {
                delete ws;
                break;
            }
        }
    } // EOF 2 pass

    //--------------------------------------------------
    // If we have a set of first time object, 
    // branch up and down to initialize pseudo-cost.
    //--------------------------------------------------
    
    numFirsts = static_cast<int> (firstObjects.size());
    //std::cout << "PSEUDO: numFirsts = " << numFirsts << std::endl;
    if (numFirsts > 0) {
        //std::cout << "PSEUDO: numFirsts = " << numFirsts << std::endl;
      
        //--------------------------------------------------
        // Backup solver status and mark hot start.
        //--------------------------------------------------
        saveLower = new double[numCols];
        saveUpper = new double[numCols];
        memcpy(saveLower, lower, numCols * sizeof(double));
        memcpy(saveUpper, upper, numCols * sizeof(double));

        CoinWarmStart * ws = solver->getWarmStart();
        solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
	aveIterations = ALPS_MIN(50, aveIterations);
        solver->setIntParam(OsiMaxNumIterationHotStart, aveIterations);
        
        solver->markHotStart();
        
        lbInd = new int [numFirsts];
        ubInd = new int [numFirsts];
            
        newLB = new double [numFirsts];
        newUB = new double [numFirsts];
            
        for (i = 0; i < numFirsts && bStatus != -2; ++i) {

            colInd = firstObjects[i]->columnIndex();
            
            lpX = saveSolution[colInd];
            
            BlisStrongBranch(model, objValue, colInd, lpX,
                             saveLower, saveUpper,
                             downKeep, downGood, downDeg,
                             upKeep, upGood, upDeg);
            
            if(!downKeep && !upKeep) {
                // Both branch can be fathomed
                bStatus = -2;
            }
            else if (!downKeep) {
                // Down branch can be fathomed.
                lbInd[numLowerTightens] = colInd;
                newLB[numLowerTightens++] = ceil(lpX);
            }
            else if (!upKeep) {
                // Up branch can be fathomed.
                ubInd[numUpperTightens] = colInd;
                newUB[numUpperTightens++] = floor(lpX);
            }
        }

        //--------------------------------------------------
        // Set new bounds in lp solver for resolving
        //--------------------------------------------------
        
        if (bStatus != -2) {
            if (numUpperTightens > 0) {
                bStatus = -1;
                for (i = 0; i < numUpperTightens; ++i) {
                    solver->setColUpper(ubInd[i], newUB[i]);
                }
            }
            if (numLowerTightens > 0) {
                bStatus = -1;
                for (i = 0; i < numLowerTightens; ++i) {
                    solver->setColLower(lbInd[i], newLB[i]);
                }
            }
        }
	
        //--------------------------------------------------
        // Unmark hotstart and recover LP solver.
        //--------------------------------------------------
        
        solver->unmarkHotStart();
        solver->setColSolution(saveSolution);
        solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
        solver->setWarmStart(ws);
        delete ws;
    }
    
    if (bStatus < 0) {
	goto TERM_CREATE;
    }
    else {
        // Create a set of candidate branching objects. 
        numBranchObjects_ = numInfs;
        branchObjects_ = new BcpsBranchObject* [numInfs];        
        
        // NOTE: it set model->savedLpSolution.
        
        sumDeg = 0.0;
	
        for (i = 0; i < numInfs; ++i) {

            if (infObjects[i]->pseudocost().getUpCost() < 
                infObjects[i]->pseudocost().getDownCost()) {
                preferDir = 1;
            }
            else {
                preferDir = -1;
            }
            
            branchObjects_[i] = infObjects[i]->createBranchObject(model,
                                                                  preferDir);
            score = infObjects[i]->pseudocost().getScore();
            branchObjects_[i]->setUpScore(score);
            sumDeg += score;
            

#ifdef BLIS_DEBUG_MORE
            std::cout << "col[" << infObjects[i]->columnIndex() << "]: score="
                      << score << ", dir=" << branchObjects_[i]->getDirection()
                      << ", up=" << infObjects[i]->pseudocost().getUpCost()
                      << ", down=" << infObjects[i]->pseudocost().getDownCost()
                      << std::endl;
#endif
        }
        
        model->setSolEstimate(objValue + sumDeg);
    }
    

 TERM_CREATE:
    
    //------------------------------------------------------
    // Cleanup.
    //------------------------------------------------------

    delete [] lbInd;
    delete [] ubInd;
    delete [] newLB;
    delete [] newUB;
    delete [] saveSolution;
    delete [] saveLower;
    delete [] saveUpper;

    return bStatus;

}
