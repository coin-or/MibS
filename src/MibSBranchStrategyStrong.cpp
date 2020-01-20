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

#include "MibSBranchStrategyStrong.hpp"
#include "MibSModel.hpp"

#include "BlisObjectInt.h"
//#include "BlisModel.h"
#include "BlisHelp.h"


//#############################################################################
MibSBranchStrategyStrong::MibSBranchStrategyStrong()
   : BlisBranchStrategyStrong()
{


}

//#############################################################################
MibSBranchStrategyStrong::MibSBranchStrategyStrong(BlisModel *model)
   : BlisBranchStrategyStrong(model)
{

}

//#############################################################################
MibSBranchStrategyStrong::~MibSBranchStrategyStrong()
{

}

//#############################################################################

/** Create a set of candidate branching objects. */
int 
MibSBranchStrategyStrong::createCandBranchObjects(int numPassesLeft, double ub)
{
    int bStatus = 0;
    int i, j, pass;
    int colInd, objInd;
    int numInfs = 0;
    int lastObj = -1;

    bool downKeep, downGood, upKeep, upGood, roundAgain(false);
    double lpX, downDeg, upDeg, sumDeg = 0.0;

    int numLowerTightens = 0;
    int numUpperTightens = 0;

    int *lbInd = NULL;
    int *ubInd = NULL;
    double *newLB = NULL;
    double *newUB = NULL;
    
    BlisObjectInt *intObject = NULL;
    
    BlisModel *model = dynamic_cast<BlisModel *>(model_);
    MibSModel *mibsmodel = dynamic_cast<MibSModel *>(model);
    MibSBilevel *bS = mibsmodel->bS_;
    OsiSolverInterface * solver = model->solver();
    
    int numCols = model->getNumCols();
    int numObjects = model->numObjects();
    bool beforeSolution = (model->getNumSolutions() == 0);
    
    int givenStrongLen = dynamic_cast<BlisParams*>
        (model->BlisPar())->entry(BlisParams::strongCandSize);
    int strongLen = givenStrongLen;
    int maxStrongLen = CoinMax(CoinMin(givenStrongLen, numObjects), 1);
    
    BlisStrong * candStrongs = new BlisStrong [maxStrongLen];

    double etol = mibsmodel->etol_;
    int uN = mibsmodel->upperDim_;
    int * upperColInd = mibsmodel->getUpperColInd();
    int * fixedInd = mibsmodel->fixedInd_;
    char * colType = mibsmodel->colType_;

    // If upper-level variable is fixed -> fixedVar = 1
    int *fixedVar = new int[numCols]();

    int *candidate = new int[numCols]();
    

    //------------------------------------------------------
    // Allocate memory and store current solution info.
    //------------------------------------------------------

    // Objective value. 
    double objValue = solver->getObjSense() * solver->getObjValue();
    
    // Column bounds.
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    double * saveUpper = new double[numCols];
    double * saveLower = new double[numCols];
    memcpy(saveLower, lower, numCols * sizeof(double));
    memcpy(saveUpper, upper, numCols * sizeof(double));
    
    // Primal solution.
    double * saveSolution = new double[numCols];
    memcpy(saveSolution, solver->getColSolution(),numCols * sizeof(double));

    lbInd = new int [numCols];
    ubInd = new int [numCols];    
    newLB = new double [numCols];
    newUB = new double [numCols];

    //------------------------------------------------------
    // Select a set of objects based on feasibility. 
    // NOTE: we might go round this loop twice if we are feed in
    //       a "feasible" solution.
    //------------------------------------------------------
    
    for (pass = 0; pass < 2; ++pass) {
	
	// Compute how many infeasible objects. 

	sumDeg = 0.0; 
	numInfs = 0;

	BcpsObject * object = NULL;

        int preferDir;        
	int leastFrac = 0;
        double infeasibility;
        double minInf = ALPS_ZERO;
        
	for (i = 0; i < maxStrongLen; ++i) {
	    candStrongs[i].bObject = NULL;
	}
	
	strongLen = 0;

	MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
	    (mibsmodel->MibSPar_->entry(MibSParams::branchStrategy));

	int index(0), found(0);
	double value(0.0);

	if(branchPar == MibSBranchingStrategyLinking){
	    for (i = 0; i < uN; ++i){
		index = upperColInd[i];
		if (fabs(lower[index]-upper[index])<=etol){
		    fixedVar[index]=1;
		}
	    }
	}

	for (i = 0; i < numCols; ++i) {
	    if(colType[i] == 'C'){
		candidate[i] = 0;
	    }
	    else{
		candidate[i] = 2;
		value = saveSolution[i];
		infeasibility = fabs(floor(value + 0.5) - value);
		if(branchPar == MibSBranchingStrategyLinking){
		    if((fixedInd[i] == 1) && (fabs(infeasibility) > etol)){
			found = 1;
			break;
		    }
		    if((bS->isLinkVarsFixed_ == true) && (bS->isIntegral_ == false)){
			if(fabs(infeasibility) > etol){
			    candidate[i] = 1;
			}
		    }
		    else if((fixedInd[i] == 1) && (((found == 0) &&
						    (fixedVar[i] != 1)) ||
						   (fabs(infeasibility) > etol))){
			candidate[i] = 1;
		    }
		}
		else if(infeasibility > etol){
		    candidate[i] = 1;
		}
	    }
	}

	//*******
	
	//for (i = 0; i < numObjects; ++i) {
            
            // TODO: currently all integer object.
	    //intObject = dynamic_cast<BlisObjectInt *>(model->objects(i));
	    //infeasibility = intObject->infeasibility(model, preferDir);
            
	    //if (infeasibility) {

	index = -1;
	for (i = 0; i < numCols; ++i) {
	    if(candidate[i] != 0){
		index ++;
		object = model->objects(index);
	    }

	    if (candidate[i] == 1) {
		++numInfs;
                intObject = dynamic_cast<BlisObjectInt *>(object);
		infeasibility = intObject->infeasibility(model, preferDir);
		
		// Increase estimated degradation to solution
		sumDeg += intObject->pseudocost().getScore();
		
		// Check for suitability based on infeasibility.
                if ((infeasibility > minInf) ||
		    ((branchPar == MibSBranchingStrategyLinking) &&
		     (found == 0))){
		    
                    if (candStrongs[leastFrac].bObject) {
                        // The slot already has one, free it.
                        delete candStrongs[leastFrac].bObject;
                    }
                    
                    // Create new branching object.
		    candStrongs[leastFrac].bObject = 
			intObject->createBranchObject(model, preferDir);

		    candStrongs[leastFrac].bObject->setUpScore(infeasibility);
                    candStrongs[leastFrac].bObject->setObjectIndex(i);
                    candStrongs[leastFrac].objectIndex = i;
                    
		    strongLen = CoinMax(strongLen, leastFrac + 1);
                    
                    // Find an empty or the worst slot.
		    leastFrac = -1;
		    minInf = ALPS_INFINITY;
                    
		    for(j = 0; j < maxStrongLen; ++j) {
                        if (!candStrongs[j].bObject) {
                            // j is an empty slots.
                            minInf = ALPS_ZERO;
                            leastFrac = j;
                            break;
                        }
			//sahar:To Do:If there is no fractional valued var, we consider
			//the first "maxStrongLen" objects 
			else if(((branchPar == MibSBranchingStrategyFractional)
				 || (found == 1)) &&
				(candStrongs[j].bObject->getUpScore() < minInf)){
			    minInf = candStrongs[j].bObject->getUpScore();
			    leastFrac = j;
			}
		    }
		}
	    }
	}

	if (numInfs) {
#ifdef BLIS_DEBUG_MORE
            std::cout << "STRONG: numInfs = " << numInfs
                      << std::endl;
#endif
            model->setSolEstimate(objValue + sumDeg);
            
	    break;
	}
	else if (pass == 0) {
	    // The first pass and ip feasible
            
#ifdef BLIS_DEBUG_MORE
	    std::cout << "STRONG: given a feasible sol" << std::endl;
#endif
            
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
		// Need resolve, then do the second round selection.
		solver->setWarmStart(ws);
		delete ws;

		// Resolve.
		solver->resolve();
		
		// Save new lp solution.
		memcpy(saveSolution, 
		       solver->getColSolution(),
		       numCols * sizeof(double));

		if (!solver->isProvenOptimal()) {
		    // Become infeasible, can do nothing. 
		    bStatus = -2;
		    break;
		}
	    } 
	    else {
		delete ws;
		break;
	    }
	}
    }

    //------------------------------------------------------
    // Now we have a set of candidate branching object, evaluate them.
    //------------------------------------------------------
    
    int bestBO = 0;
    double bestScore = 0.0;
    
    for (i = 0; i < strongLen; ++i) {
	candStrongs[i].numIntInfUp = numInfs;
	candStrongs[i].numIntInfDown = numInfs;
	if ( !model->objects(candStrongs[i].bObject->getObjectIndex())->
             boundBranch(model) ) {
            // Weild branching: not by modifying variable bounds.
            // Exit selection routine.
	    strongLen = 0;
        }
	// Find the most fractional object in case of doing simple branch
	if (candStrongs[i].bObject->getUpScore() > bestScore) {
	    bestScore = candStrongs[i].bObject->getUpScore();
	    bestBO = i;
	}
    }
    
    // If we have hit max time don't do strong branching
    double timeLimit = model->AlpsPar()->entry(AlpsParams::timeLimit);
    AlpsKnowledgeBroker *broker = model->getKnowledgeBroker();
    bool maxTimeReached = (broker->timer().getTime() > timeLimit);
    
#ifdef BLIS_DEBUG_MORE
    printf("1. strongLen = %d, maxTimeReached %d, numPassesLeft %d\n", 
	   strongLen, maxTimeReached, numPassesLeft);
#endif
    
    if (strongLen <= 0 || maxTimeReached) {
        
        //--------------------------------------------------
        // Simple max infeasibility branching.
        //--------------------------------------------------
        
#ifdef BLIS_DEBUG
        std::cout << "NOT STRONG: maxTimeReached=" << maxTimeReached
                  << "; numPassesLeft=" << numPassesLeft 
                  << std::endl;
#endif

        // Create candidate set.
        numBranchObjects_ = 1;
        branchObjects_ = new BcpsBranchObject * [1];
        branchObjects_[0] = candStrongs[bestBO].bObject;
        candStrongs[bestBO].bObject = NULL;
        strongLen = 0; // Stop here.
    }
    else {
        
        //--------------------------------------------------------
        // Strong branching
        //--------------------------------------------------------
        
        // set true to say look at all even if some fixed (experiment)
	bool solveAll = false;
	int saveLimit;
        
	CoinWarmStart * ws = solver->getWarmStart();
	solver->getIntParam(OsiMaxNumIterationHotStart, saveLimit);
	if (beforeSolution) {
	    solver->setIntParam(OsiMaxNumIterationHotStart, 100000);
        }
        else {
            solver->setIntParam(OsiMaxNumIterationHotStart, 1000);
        }
        
        // Mark hot start
        solver->markHotStart();
        
#if 0
	printf("BEFORE LOOP: strongLen = %d\n",strongLen);
#endif
        
	for (i = 0; i < strongLen; ++i) {

#if 1
            objInd = candStrongs[i].bObject->getObjectIndex();
            colInd = model->getIntColIndices()[objInd];
            lpX = saveSolution[colInd];
            
            BlisStrongBranch(model, objValue, colInd, lpX,
                             saveLower, saveUpper,
                             downKeep, downGood, downDeg,
                             upKeep, upGood, upDeg);
            
            candStrongs[i].bObject->setDownScore(downDeg);
            candStrongs[i].bObject->setUpScore(upDeg);

            if(!downKeep && !upKeep) {
                // This node can be fathomed
                bStatus = -2;
                break;
            }
            else if (!downKeep) {
                // Down branch can be fathomed.
                lbInd[numLowerTightens] = colInd;
                newLB[numLowerTightens++] = ceil(lpX);
                lastObj = i;
                if (numPassesLeft == 0) break;
            }
            else if (!upKeep) {
                // Up branch can be fathomed.
                ubInd[numUpperTightens] = colInd;
                newUB[numUpperTightens++] = floor(lpX);
                lastObj = i;
                if (numPassesLeft == 0) break;
            }
#else
	    double objChange;
	    double newObjValue = ALPS_DBL_MAX;
            
	    // status is 0 finished, 1 infeasible and other
	    int lpStatus;            

            //----------------------------------------------
	    // Branching down.
            //----------------------------------------------
            
            candStrongs[i].bObject->setDirection(-1);
            candStrongs[i].bObject->branch();
            solver->solveFromHotStart();

            if (solver->isProvenOptimal()) {
                lpStatus = 0; // optimal
            }
            else if (solver->isIterationLimitReached()
                     &&!solver->isDualObjectiveLimitReached()) {
                lpStatus = 2; // unknown 
            }
            else {
                lpStatus = 1; // infeasible
            }
            
            newObjValue = solver->getObjSense() * solver->getObjValue();
            objChange = newObjValue-objValue ;
            
#ifdef BLIS_DEBUG_MORE
	    std::cout << "Down: lpStatus = " << lpStatus << std::endl;
#endif

	    if (lpStatus == 0) {

                // Update pseudocost
                ind = candStrongs[i].objectIndex;
                intObject = dynamic_cast<BlisObjectInt *>(model->objects(ind));
                col = intObject->columnIndex();
                lpX = saveSolution[col];
                intObject->pseudocost().update(-1, objChange, lpX);
                model->setSharedObjectMark(intObject->getObjectIndex());

		candStrongs[i].finishedDown = true ;
		if (newObjValue > model->getCutoff()) {
		    objChange = ALPS_DBL_MAX; // say infeasible
		} 
                else {    
		    // See if integer solution
                    BlisSolution* ksol = 
                        model->feasibleSolution(candStrongs[i].numIntInfDown,
                                                candStrongs[i].numObjInfDown);
                    
		    if (ksol) {
#ifdef BLIS_DEBUG_MORE
			printf("STRONG:down:found a feasible solution\n");
#endif
			
			model->storeSolution(BlisSolutionTypeStrong, ksol);
			objChange = ALPS_DBL_MAX ;
		    }
		}
	    }
	    else if (lpStatus == 1) {
	      objChange = ALPS_DBL_MAX ;
	    } 
	    else {
		// Can't say much as we did not finish
		candStrongs[i].finishedDown = false ;
	    }
            
	    candStrongs[i].bObject->setDownScore(objChange);
	    
	    // restore bounds
            int numDiff = 0;
            for (j = 0; j < numCols; ++j) {
                if (saveLower[j] != lower[j]) {
                    solver->setColLower(j, saveLower[j]);
                    ++numDiff;
                }
                if (saveUpper[j] != upper[j]) {
                    solver->setColUpper(j, saveUpper[j]);
                    ++numDiff;
                }
            }

#ifdef BLIS_DEBUG_MORE
            std::cout << "numDiff = " << numDiff << std::endl;
#endif	    
            
            //----------------------------------------------
	    // Branching up.
            //----------------------------------------------

            candStrongs[i].bObject->branch();
            solver->solveFromHotStart();

            if (solver->isProvenOptimal()) {
                lpStatus = 0; // optimal
            }
            else if (solver->isIterationLimitReached()
                     &&!solver->isDualObjectiveLimitReached()) {
                lpStatus = 2; // unknown 
            }
            else {
                lpStatus = 1; // infeasible
            }
            
            newObjValue = solver->getObjSense() * solver->getObjValue();
            objChange = newObjValue - objValue ;
            
#ifdef BLIS_DEBUG_MORE
	    std::cout << "STRONG: Up: lpStatus = " << lpStatus << std::endl;
#endif      

	    if (lpStatus == 0) {
                // Update pseudocost
                ind = candStrongs[i].objectIndex;
                intObject = dynamic_cast<BlisObjectInt *>(model->objects(ind));
                col = intObject->columnIndex();
                lpX = saveSolution[col];
                intObject->pseudocost().update(1, objChange, lpX);
                model->setSharedObjectMark(intObject->getObjectIndex());
                
		candStrongs[i].finishedUp = true ;
		if (newObjValue > model->getCutoff()) {
		    objChange = ALPS_DBL_MAX; // Cutoff
		} 
		else {
		    // See if integer solution
                    BlisSolution* ksol = 
                        model->feasibleSolution(candStrongs[i].numIntInfDown,
                                                candStrongs[i].numObjInfDown);
                    
		    if (ksol) {
#ifdef BLIS_DEBUG_MORE
			printf("STRONG:Up:found a feasible solution\n");
#endif
                        
			model->storeSolution(BlisSolutionTypeStrong, ksol);
			objChange = ALPS_DBL_MAX ;
		    }
		}
	    } 
	    else if (lpStatus == 1) {
		objChange = ALPS_DBL_MAX;
	    } 
            else {
		// Can't say much as we did not finish
		candStrongs[i].finishedUp = false ;
	    }
	    candStrongs[i].bObject->setUpScore(objChange);
	    
	    // restore bounds
            for (j = 0; j < numCols; ++j) {
                if (saveLower[j] != lower[j]) {
                    solver->setColLower(j,saveLower[j]);
                }
                if (saveUpper[j] != upper[j]) {
                    solver->setColUpper(j,saveUpper[j]);
                }
	    }
            
            //----------------------------------------------
            // End of evaluation for this branching object.
            // Possibilities are:
	    // 1) Both sides below cutoff; this variable is a 
            //    candidate for branching.
	    // 2) Both sides infeasible or above the obj cutoff: 
	    //    no further action here. Break from the evaluation loop and 
	    //    assume the node will be purged by the caller.
	    // 3) One side below cutoff: Install the branch (i.e., fix the 
	    //    variable). Break from the evaluation loop and assume the 
	    //    node will be reoptimised by the caller.
            //----------------------------------------------
            
	    if (candStrongs[i].bObject->getUpScore() < ALPS_INFINITY) {
		if(candStrongs[i].bObject->getDownScore() < ALPS_INFINITY) {
		    // feasible - no action
		} 
                else {
		    // up feasible, down infeasible
		    bStatus = -1;
		    if (!solveAll) {
			candStrongs[i].bObject->setDirection(1);
			candStrongs[i].bObject->branch();
			break;
		    }
		}
	    } 
            else {
		if(candStrongs[i].bObject->getDownScore() < ALPS_INFINITY) {
		    // down feasible, up infeasible
		    bStatus = -1;
		    if (!solveAll) {
			candStrongs[i].bObject->setDirection(-1);
			candStrongs[i].bObject->branch();
			break;
		    }
		} 
                else {
		    // neither side feasible
		    bStatus = -2;
		    break;
		}
	    }
#endif

	}// EOF the loop for checking each candiate
        
        //--------------------------------------------------
        // Set new bounds in lp solver for resolving
        //--------------------------------------------------
        
        if (bStatus != -2 && numPassesLeft > 0) {
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
        // Unmark hotstart and reset lp solver.
        //--------------------------------------------------
        
        solver->unmarkHotStart();
        solver->setColSolution(saveSolution);
        solver->setIntParam(OsiMaxNumIterationHotStart, saveLimit);
        solver->setWarmStart(ws);
        delete ws;
    }

    //------------------------------------------------------
    // Record the list of candidates.
    //------------------------------------------------------

    if (bStatus >= 0 && strongLen > 0) {
        if (numPassesLeft == 0 && lastObj > -1) {
            numBranchObjects_ = 1;
            branchObjects_ = new BcpsBranchObject* [1]; 
            for (i = 0; i < strongLen; ++i) {
                if (i == lastObj) {
                    branchObjects_[0] = candStrongs[i].bObject;
                    candStrongs[i].bObject = NULL;
                    //std::cout << "strong: one direction is good\n";
                }
            }
        }
        else {
            // Store the set of candidate branching objects. 
            numBranchObjects_ = strongLen;
            branchObjects_ = new BcpsBranchObject* [strongLen];        
            for (i = 0; i < strongLen; ++i) {
                branchObjects_[i] = candStrongs[i].bObject;
                candStrongs[i].bObject = NULL;
            }
        }
    }  
    
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
    delete [] fixedVar;
    delete [] candidate;
    for (i = 0; i < maxStrongLen; ++i) {
        if (candStrongs[i].bObject) {
            delete candStrongs[i].bObject;
        }
    }
    delete [] candStrongs;

    return bStatus;
}
