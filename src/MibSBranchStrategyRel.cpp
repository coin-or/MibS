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

#include "MibSBranchStrategyRel.hpp"
#include "MibSModel.hpp"

#include "BlisObjectInt.h"
//#include "BlisModel.h"
#include "BlisHelp.h"

//#############################################################################

struct BlisPseuoGreater
{
    bool operator()(double x, double y) const {
        return (x > y);
    }
};

//#############################################################################
MibSBranchStrategyRel::MibSBranchStrategyRel()
   : BlisBranchStrategyRel()
{


}

//#############################################################################
MibSBranchStrategyRel::MibSBranchStrategyRel(BlisModel *model, int rel)
   : BlisBranchStrategyRel(model, rel)
{

}

//#############################################################################
MibSBranchStrategyRel::~MibSBranchStrategyRel()
{

}

//#############################################################################

/** Create a set of candidate branching objects. */
int 
MibSBranchStrategyRel::createCandBranchObjects(int numPassesLeft, double ub)
{

   //FIXME::ADD IN BILEVEL BRANCHING

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

    double * saveUpper = NULL;
    double * saveLower = NULL;
    double * saveSolution = NULL;

    
    BlisModel *model = dynamic_cast<BlisModel *>(model_);
    OsiSolverInterface * solver = model->solver();
    
    int numCols = model->getNumCols();
    int numObjects = model->numObjects();
    
    //int lookAhead = dynamic_cast<BlisParams*>
    //  (model->blisPar())->entry(BlisParams::lookAhead);    
    
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
        printf("REL: CREATE: maxTimeReached %d, numPassesLeft %d\n", 
               maxTimeReached, numPassesLeft);
#endif
    }

    
    // Store first time objects.
    std::vector<BlisObjectInt *> firstObjects;

    // Store infeasible objects.
    std::vector<BlisObjectInt *> infObjects;

    // TODO: check if sorting is expensive.
    std::multimap<double, BlisObjectInt*, BlisPseuoGreater> sortedObjects;

    double objValue = solver->getObjSense() * solver->getObjValue();

    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();

    int lookAhead = dynamic_cast<BlisParams*>
	(model->BlisPar())->entry(BlisParams::lookAhead);

    BlisObjectInt * intObject = NULL;

    //------------------------------------------------------
    // Backup solver status and mark hot start.
    //-----------------------------------------------------

    saveSolution = new double[numCols];
    memcpy(saveSolution, solver->getColSolution(), numCols*sizeof(double));
    saveLower = new double[numCols];
    saveUpper = new double[numCols];
    memcpy(saveLower, lower, numCols * sizeof(double));
    memcpy(saveUpper, upper, numCols * sizeof(double));

    //------------------------------------------------------
    // Find the infeasible objects.
    // NOTE: we might go round this loop twice if we are feed in
    //       a "feasible" solution.
    //------------------------------------------------------
    
    for (pass = 0; pass < 2; ++pass) {
	
        numInfs = 0;
        
        BcpsObject * object = NULL;

            
        infObjects.clear();
        firstObjects.clear();
        
        for (i = 0; i < numObjects; ++i) {
                
            object = model->objects(i);
            infeasibility = object->infeasibility(model, preferDir);
            
            if (infeasibility) {
                
                ++numInfs;
                intObject = dynamic_cast<BlisObjectInt *>(object);
                
                if (intObject) {
                        
                    //score = object->pseudocost().getScore();
                    //tempBO = object->createBranchObject(model, preferDir);
                    //candObjects.insert(std::make_pair(score, tempBO));
                    //tempBO = NULL;

                    infObjects.push_back(intObject);
                    
                    if (!selectNow) {
                        minCount = 
                            ALPS_MIN(intObject->pseudocost().getDownCount(),
                                     intObject->pseudocost().getUpCount());
                        
                        if (minCount < 1) {
                            firstObjects.push_back(intObject);
                        }
                    }

#ifdef BLIS_DEBUG_MORE
                    if (intObject->columnIndex() == 15) {
                        std::cout << "x[15] = " << saveSolution[15] 
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
#ifdef BLIS_DEBUG_MORE
            std::cout << "REL: numInfs = " << numInfs
                      << std::endl;
#endif
            break;
        }
        else if (pass == 0) {
            // The first pass and is IP feasible.
            
#ifdef BLIS_DEBUG
            std::cout << "REL: given a feasible sol" << std::endl;
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
    if (numFirsts > 0) {
	
        CoinWarmStart * ws = solver->getWarmStart();
        solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
	int maxIter = ALPS_MAX(model->getAveIterations(), 50);
        solver->setIntParam(OsiMaxNumIterationHotStart, maxIter);

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
                //break;
            }
            else if (!upKeep) {
                // Up branch can be fathomed.
                ubInd[numUpperTightens] = colInd;
                newUB[numUpperTightens++] = floor(lpX);
                // break;
            }
            
            // Update pseudocost.
            if(downGood) {
                firstObjects[i]->pseudocost().update(-1, downDeg, lpX);
                model->setSharedObjectMark(firstObjects[i]->getObjectIndex());
            }
            if(downGood) {
                firstObjects[i]->pseudocost().update(1, upDeg, lpX);
                model->setSharedObjectMark(firstObjects[i]->getObjectIndex());
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
    
    //std::cout << "REL: bStatus = " << bStatus << std::endl;

    if (bStatus < 0) {
        // Infeasible or monotone.
	goto TERM_CREATE;
    }
    else {
        // All object's pseudocost have been initialized.
        // Sort them, and do strong branch for the unreliable one
        // NOTE: it set model->savedLpSolution.
        
        sumDeg = 0.0;
        
        for (i = 0; i < numInfs; ++i) {
            score = infObjects[i]->pseudocost().getScore();
            sumDeg += score;
            
	    std::pair<const double, BlisObjectInt*> sa(score, infObjects[i]);
            sortedObjects.insert(sa);
            
#ifdef BLIS_DEBUG_MORE
            std::cout << "col[" << infObjects[i]->columnIndex() << "]="
                      << score << ", "<< std::endl;
#endif
        }

        int numNotChange = 0;

        std::multimap< double, BlisObjectInt*, BlisPseuoGreater>::iterator pos;
        
        CoinWarmStart * ws = solver->getWarmStart();
        solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
	int maxIter = ALPS_MAX(model->getAveIterations(), 50);
        solver->setIntParam(OsiMaxNumIterationHotStart, maxIter);
        solver->markHotStart();
        
        BlisObjectInt *bestObject = NULL;
        double bestScore = -10.0;
        
        for (pos = sortedObjects.begin(); pos != sortedObjects.end(); ++pos) {
	    
	    intObject  = pos->second;
            
            colInd = intObject->columnIndex();
            
#ifdef BLIS_DEBUG_MORE
            std::cout << "col[" << colInd << "]: " 
                      << "score=" << pos->first 
                      << ", upCount=" << intObject->pseudocost().getUpCount()
                      <<", downCount="<< intObject->pseudocost().getDownCount()
                      << std::endl;
#endif
            
            // Check if reliable.
	    int objRelibility=ALPS_MIN(intObject->pseudocost().getUpCount(),
                                       intObject->pseudocost().getDownCount());
	    
            if (objRelibility < relibility_) {
		// Unrelible object. Do strong branching. 
		
                
                lpX = saveSolution[colInd];
            
                BlisStrongBranch(model, objValue, colInd, lpX,
                                 saveLower, saveUpper,
                                 downKeep, downGood, downDeg,
                                 upKeep, upGood, upDeg);
                // Update pseudocost.
                if(downGood) {
                    intObject->pseudocost().update(-1, downDeg, lpX);
                }
                if(downGood) {
                    intObject->pseudocost().update(1, upDeg, lpX);
                }
	    }
	    
            // Compare with the best. 
            if (intObject->pseudocost().getScore() > bestScore) {
                bestScore = intObject->pseudocost().getScore();
                bestObject = intObject;
                // Reset 
                numNotChange = 0;
            }
	    else {
                // If best doesn't change for "lookAhead" comparisons, then 
                // the best is reliable.
		if (++numNotChange > lookAhead) {
                    if (bestObject->pseudocost().getUpCost() >
                        bestObject->pseudocost().getDownCost()) {
                        preferDir = 1;
                    }
                    else {
                        preferDir = -1;
                    }
                    break;
                }
	    }
        }
        
        solver->unmarkHotStart();
        solver->setColSolution(saveSolution);
        solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
        solver->setWarmStart(ws);
        delete ws;
        
        model->setSolEstimate(objValue + sumDeg);

        assert(bestObject != NULL);
        bestBranchObject_ = bestObject->createBranchObject(model, preferDir);
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
