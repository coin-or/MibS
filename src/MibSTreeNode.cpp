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

#include "MibSModel.hpp"
#include "MibSTreeNode.hpp"
#include "BlisObjectInt.h"
#include "BlisConstraint.h"

//#############################################################################
MibSTreeNode::MibSTreeNode()
   : BlisTreeNode()
{

   lowerUpperBound_ = - ALPS_DBL_MAX;
   boundSet_ = false;

}

//#############################################################################
MibSTreeNode::MibSTreeNode(AlpsNodeDesc *&desc)
   :  BlisTreeNode(desc)
{

   lowerUpperBound_ = - ALPS_DBL_MAX;
   boundSet_ = false;

}

//#############################################################################
MibSTreeNode::MibSTreeNode(BlisModel *m)
   :  BlisTreeNode(m)
{

   lowerUpperBound_ = - ALPS_DBL_MAX;
   boundSet_ = false;

}

//#############################################################################
MibSTreeNode::~MibSTreeNode()
{

}

//#############################################################################
AlpsTreeNode*
MibSTreeNode::createNewTreeNode(AlpsNodeDesc *&desc) const
{
    double estimate = solEstimate_;

    // Set solution estimate for this nodes.
    // double solEstimate = quality_ + sum_i{min{up_i, down_i}}
    int branchDir = dynamic_cast<BlisNodeDesc *>(desc)->getBranchedDir();
    int branchInd = dynamic_cast<BlisNodeDesc *>(desc)->getBranchedInd();
    double lpX = dynamic_cast<BlisNodeDesc *>(desc)->getBranchedVal();
    double f = lpX - floor(lpX);
    //assert(f > 0.0);
    
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());
    int objInd = model->getIntObjIndices()[branchInd];
    BlisObjectInt *obj = dynamic_cast<BlisObjectInt *>(model->objects(objInd));
    
    if (branchDir == -1) {
        estimate -= (1.0-f) * obj->pseudocost().getUpCost();
    }
    else {
	estimate -= f * obj->pseudocost().getDownCost();
    }
    
#ifdef BLIS_DEBUG_MORE
    printf("BLIS:createNewTreeNode: quality=%g, solEstimate=%g\n",
           quality_, solEstimate_);
#endif

    // Create a new tree node
    BlisTreeNode *node = new MibSTreeNode(desc);    
    desc = NULL;
    
    return node;
}

//#############################################################################

// NOTE: if rampup,
// - parent must be explicit if not NULL,
// - this node is explicit.

int
MibSTreeNode::process(bool isRoot, bool rampUp)
{
    BlisReturnStatus returnStatus = BlisReturnStatusUnknown;
    BlisLpStatus lpStatus = BlisLpStatusUnknown;
    int j, k = -1;
    int numCols, numRows, numCoreCols, numCoreRows;
    int numStartRows, origNumStartRows;
    int maxNumCons, tempNumCons = 0;

    int numIntInfs = 0;
    int numObjInfs = 0;

    int voilatedNumCons = 0;
    int origNumOldCons = 0;
    int currNumOldCons = 0;
    int newNumCons = 0;
    int maxNewNumCons = 0;

    int numAppliedCons = 0;
    int cutStrategy;
    
    // Only autmatic stategy has depth limit.
    int maxConstraintDepth = 20;

    int numPassesLeft = 0;      
    int bStatus = -1;
 
    double cutoff = ALPS_INC_MAX;
    double parentObjValue = getQuality();
    double preObjValue = -ALPS_OBJ_MAX;
    double improvement = 100.0;

    double  heurObjValue;
    double *heurSolution = NULL;
    double *currLpSolution = NULL;

    bool keepOn = true;
    bool needBranch = false;
    bool lpFeasible = false;
    bool foundSolution = false;
    bool genConsHere = false;
    bool shareCon = false;
    bool shareVar = false;

    CoinWarmStartBasis::Status rowStatus;
    BlisConstraint *aCon = NULL;
    BcpsObject **newConstraints = NULL;

    BlisSolution *ipSol = NULL;
    
    int numDelRows = 0;
    int *delRow = NULL;
    int *oldConsPos = NULL;
    
    std::vector<int> delIndices;
    
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());
    MibSModel *mibsModel = dynamic_cast<MibSModel *>(model);
    MibSBilevel *bS = mibsModel->bS_;

    AlpsPhase phase = knowledgeBroker_->getPhase();

    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);
    int hubMsgLevel = model->AlpsPar()->entry(AlpsParams::hubMsgLevel);
    int workerMsgLevel = model->AlpsPar()->entry(AlpsParams::workerMsgLevel);

    BlisParams * BlisPar = model->BlisPar();

    int maxPass = BlisPar->entry(BlisParams::cutPass);
    double tailOffTol = BlisPar->entry(BlisParams::tailOff);

    MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
	(mibsModel->MibSPar_->entry(MibSParams::branchStrategy));

    //tailOffTol = 1e-7;

    /*if(bS->useBilevelBranching_ == false){
	tailOffTol = -1000;
	}*/

    if (maxPass < ALPS_INT_MAX) {
	++maxPass;
    }
    
    shareCon = BlisPar->entry(BlisParams::shareConstraints);
    shareVar = BlisPar->entry(BlisParams::shareVariables);    

    cutoff = model->getCutoff();

    numPassesLeft = model->getNumBranchResolve();
    
    // std::cout << "numPassesLeft = " << numPassesLeft << std::endl;
    
    //------------------------------------------------------
    // Check if this can be fathomed by objective cutoff.
    //------------------------------------------------------
    
#if 0
    std::cout << "parentObjValue = " << parentObjValue 
              << "cutoff = " << cutoff << std::endl;
#endif

    if (parentObjValue > cutoff) {
	setStatus(AlpsNodeStatusFathomed);
        //std::cout << "fathom!" <<std::endl;
	goto TERM_PROCESS;
    }
    else {
        model->setActiveNode(this);
        model->addNumNodes();
    }
    
    //------------------------------------------------------
    // Get model information and parameters.
    //------------------------------------------------------

    numCols = model->solver()->getNumCols();
    numRows = model->solver()->getNumRows();
    numCoreCols = model->getNumCoreVariables();
    numCoreRows = model->getNumCoreConstraints();
    numAppliedCons =  numRows - numCoreRows;
    
    maxNumCons = model->getMaxNumCons();

    heurSolution = new double [numCols];
    currLpSolution = new double [numCols];

    //------------------------------------------------------
    // Decides if can generate constraints.
    //------------------------------------------------------

    // Mark if this node is root or not.
    model->isRoot_ = isRoot;

    genConsHere = false;

    cutStrategy = model->getCutStrategy();

    assert(cutStrategy != BlisCutStrategyNotSet);
    
    if (cutStrategy == BlisCutStrategyNone) {
	genConsHere = false;
    }
    else if (cutStrategy == BlisCutStrategyRoot) {
	// The original root only
	if (isRoot && (index_ == 0)) genConsHere = true;
    }
    else if (cutStrategy == BlisCutStrategyAuto) {
	if (depth_ < maxConstraintDepth) {
           genConsHere = true;
	}
    }
    else if (cutStrategy == BlisCutStrategyPeriodic) {
	genConsHere = true;
    }
    else {
	genConsHere = true;
    }

    if (genConsHere && (phase == AlpsPhaseRampup)) {
	if (!(BlisPar->entry(BlisParams::cutRampUp))) {
	    genConsHere = false;
	}
    }

    //======================================================
    // Restore, load and solve the subproblem.
    // (1) LP infeasible
    //     a. set status to be fathom.
    // (2) LP feasible
    //     a. MILP feasible. Check whether need update incumbent.
    //     b. LP feasible but not MIP feasible. Check whether can be 
    //        fathomed, if not, choose a branch variable.
    //======================================================

    //------------------------------------------------------
    // Extract info from this node and load subproblem into lp solver.
    //------------------------------------------------------
    
    installSubProblem(model);

    //------------------------------------------------------
    // Bounding, heuristic searching, constraint generating.
    //------------------------------------------------------

    numStartRows = model->solver()->getNumRows();
    origNumStartRows = numStartRows;

    origNumOldCons = numStartRows - numCoreRows;
    currNumOldCons = origNumOldCons;
    
#if 0
    std::cout << "PROCESS: genConsHere =" << genConsHere
	      << ", cut strategy =" << model->getCutStrategy()
	      << ", numCoreRows =" << numCoreRows
	      << ", numStartRows =" << numStartRows
	      << ", currNumOldCons =" << currNumOldCons 
              << ", index_ = " << index_ << std::endl;
#endif
    
    if (currNumOldCons > 0) {
	oldConsPos = new int [currNumOldCons];
	for (k = 0; k < currNumOldCons; ++k) {
	    oldConsPos[k] = k;
	}
    }

    if (genConsHere) {
	if (maxNumCons > ALPS_INT_MAX - 100) {
	    maxNewNumCons = 10000;
	}
	else {
	    maxNewNumCons = maxNumCons;
	}
        newConstraints = new BcpsObject* [maxNewNumCons];    
    }

    model->boundingPass_ = 0;
    while (keepOn && (model->boundingPass_ < maxPass)) {
        ++(model->boundingPass_);
        keepOn = false;
        
        //--------------------------------------------------
        // Bounding to get the quality of this node.
        //--------------------------------------------------
        
        if ( (knowledgeBroker_->getProcType() == AlpsProcessTypeMaster) && 
             isRoot && (model->boundingPass_ == 1) ) {
            if (msgLevel > 0) {
                model->blisMessageHandler()->message(BLIS_ROOT_PROCESS, 
                                                     model->blisMessages())
                    << model->getNumRows()
                    << model->getNumCols()
                    << CoinMessageEol;
                model->solver()->messageHandler()->setLogLevel(1);
            }
            else {
                model->solver()->messageHandler()->setLogLevel(-1);
            }
            getKnowledgeBroker()->tempTimer().start();
        }
	
        lpStatus = static_cast<BlisLpStatus> (bound(model));

	if (model->boundingPass_ == 1) {
	    int iter = model->solver()->getIterationCount();
	    model->addNumIterations(iter);
            if (isRoot) {
                getKnowledgeBroker()->tempTimer().stop();
                if ((knowledgeBroker_->getProcType() == AlpsProcessTypeMaster)
		    && (msgLevel > 0)) {
                    model->solver()->messageHandler()->setLogLevel(0);
                    model->blisMessageHandler()->message(BLIS_ROOT_TIME, 
							 model->blisMessages())
                        << getKnowledgeBroker()->tempTimer().getCpuTime() 
                        << CoinMessageEol;
                }
            }
	}
        
        switch(lpStatus) {
        case BlisLpStatusOptimal:
	  //can we put preprocessor here?

	  //MibSModel * mibs = dynamic_cast<MibSModel *>(model);
	  //bool usePreprocessor =
	  // mibs->MibSPar_->entry(MibSParams::usePreprocessor);

	  //if((model->boundingPass_ == 1) && usePreprocessor){
	  // mibs->runPreprocessor();
	  //}

	  // Check if IP feasible 
	  ipSol = model->feasibleSolution(numIntInfs, numObjInfs);

	  //if((bS->useBilevelBranching_ == false) &&
	  // (bS->LPSolStatus_ != MibSLPSolStatusFeasible)){
	  if((((branchPar == MibSBranchingStrategyLinking) && (bS->isLinkVarsFixed_)) ||
	      (branchPar == MibSBranchingStrategyFractional)) && (bS->isIntegral_)){
	      tailOffTol = -1000;
	  }
	  else{
	      tailOffTol = BlisPar->entry(BlisParams::tailOff);
	  }

	  if((!ipSol) && (bS->shouldPrune_)){
	      setStatus(AlpsNodeStatusFathomed);
	      goto TERM_PROCESS;
	  }
	      
	  if (ipSol) {
		// IP feasible
		model->storeSolution(BlisSolutionTypeHeuristic, ipSol);
		//model->storeSolution(BlisSolutionTypeRounding, ipSol);
                // Update cutoff
                cutoff = model->getCutoff();
                setStatus(AlpsNodeStatusFathomed);
		goto TERM_PROCESS;
            }
            else {
                if (quality_ > cutoff) {
                    setStatus(AlpsNodeStatusFathomed);
                    goto TERM_PROCESS;
                }
                needBranch = true;
                //reducedCostFix(model);
                
                //------------------------------------------
                // Check if tailoff
                //------------------------------------------

                if (model->boundingPass_ > 1) {
                    improvement = quality_ - preObjValue;
                    if (improvement > tailOffTol) {
                        // NOTE: still need remove slacks, although
                        //       tailoff.
                        keepOn = true;
                    }
                    
#if 0
                    std::cout << "PROCESS: boundingPass_["
                              << model->boundingPass_ << "], improvement=" 
                              << improvement << ", tailOffTol=" << tailOffTol
                              << ", preObj=" << preObjValue 
                              << ", newObj=" << quality_
                              << std::endl;
#endif
                }
                else {
                    keepOn = true;
                }
                // Update previous objective value.
                preObjValue = quality_;

                //------------------------------------------
                // Remove non-core slack constraints. 
                //------------------------------------------

                numRows = model->getNumRows();
                
                if ( genConsHere &&
                     //(improvement > tailOffTol) && 
                     //(numRows > numStartRows) ) {
                     (numRows > numCoreRows) ) {   
                 
#if 1                  
                    if ( (numStartRows + newNumCons != numRows) ||
                         (numCoreRows+currNumOldCons +newNumCons != numRows) ){
                        
                        std::cout << "ERROR: numRows=" << numRows
                                  << "; numCoreRows=" << numCoreRows
                                  << "; numStartRows=" << numStartRows
                                  << "; newNumCons=" << newNumCons
                                  << "; currNumOldCons=" << currNumOldCons
                                  << std::endl;
                        
                        assert(numRows - numStartRows == newNumCons);
                    }
#endif
                    
                    int *oldDelMark = NULL;
                    if (currNumOldCons > 0) {
                        oldDelMark = new int [currNumOldCons];
                        CoinZeroN(oldDelMark, currNumOldCons);
                    }
                    int *newDelMark = NULL;
                    if (newNumCons > 0) {
                        newDelMark = new int [newNumCons];
                        CoinZeroN(newDelMark, newNumCons);
                    }
                    
                    const CoinWarmStartBasis* ws= 
                        dynamic_cast<CoinWarmStartBasis*>
                        (model->solver()->getWarmStart());
                    
                    // Make sure delIndices is empty.
                    assert(delIndices.size()==0);

#if REMOVE_SLACK
                    for (k = numCoreRows; k < numRows; ++k) {
                        rowStatus = ws->getArtifStatus(k);

                        if (rowStatus == CoinWarmStartBasis::basic) {
			    int count;
                            if (k < numStartRows) {
				BlisConstraint *tmpCon = 
				    model->oldConstraints()[(k-numCoreRows)];
				count = tmpCon->getNumInactive() + 1;
				tmpCon->setNumInactive(count);
				if (tmpCon->getNumInactive() > BLIS_SLACK_MAX){
				    oldDelMark[(k-numCoreRows)] = 1;
				    delIndices.push_back(k);
				}
                            }
                            else {
				BcpsObject* tmpCon = 
				    newConstraints[(k-numStartRows)];
				count = tmpCon->getNumInactive() + 1;
				tmpCon->setNumInactive(count);
				if (tmpCon->getNumInactive() > BLIS_SLACK_MAX){
				    newDelMark[(k-numStartRows)] = 1;
				    delIndices.push_back(k);
				}
                            }
                        }
                    }
#endif
                    numDelRows = static_cast<int> (delIndices.size());
		    
		    if (msgLevel > 100) {
			std::cout << "PROCESS: new cuts=" << newNumCons 
				  << ", delete slack cuts=" << numDelRows 
				  << ", numRows=" << numRows 
				  << ", numStartRows=" <<numStartRows 
				  << std::endl;
		    }
    
                    
                    if (numDelRows > 0) {
                        delRow = new int [numDelRows];
                        for (k = 0; k < numDelRows; ++k) {
                            delRow[k] = delIndices[k];
#ifdef BLIS_DEBUG
			    std::cout << "REMOVE: slack row " << delRow[k] 
				      << std::endl;
#endif
                        }
                        
                        //----------------------------------
                        // Delete from lp solver.
                        //----------------------------------
                        
                        model->solver()->deleteRows(numDelRows, delRow);
                        
                        delete [] delRow;
                        delRow = NULL;
                        delIndices.clear();

                        //----------------------------------
                        // Delete from the old cut position array.
                        //----------------------------------

                        tempNumCons = 0;
                        for (k = 0; k < currNumOldCons; ++k) {
                            if (oldDelMark[k] != 1) {
                                // Survived
                                oldConsPos[tempNumCons++] = oldConsPos[k];    
                            }
                        }
                        currNumOldCons = tempNumCons;
                        numStartRows = numCoreRows + currNumOldCons;
                        
                        //----------------------------------
                        // Delete from new cut vector.
                        //----------------------------------

                        //std::cout << std::endl;
                        tempNumCons = 0;
                        for (k = 0; k < newNumCons; ++k) {
                            if (newDelMark[k] == 1) {
                                // Deleted
#ifdef BLIS_DEBUG_MORE
                                std::cout << "delete cut " << k 
                                          << ", size=" 
                                          << dynamic_cast<BlisConstraint*>(newConstraints[k])->getSize()
                                          << std::endl;
#endif
                                
                                delete newConstraints[k];
                                newConstraints[k] = NULL;
                            }
                            else {
                                // Survived
                                newConstraints[tempNumCons++] = newConstraints[k];
                            }
                        }
                        //assert(tempNumCons + numDelRows == newNumCons);
                        numAppliedCons -= newNumCons;
                        numAppliedCons += tempNumCons;
                        newNumCons = tempNumCons;                        
                        
                        //----------------------------------
                        // Resolve to update solution info in lp solver.
                        //----------------------------------
                        
                        int easy = 2;
                        model->solver()->setHintParam(OsiDoInBranchAndCut,
                                                      true, OsiHintDo, &easy);
                        model->solver()->resolve();
                        model->solver()->setHintParam(OsiDoInBranchAndCut,
                                                      true, OsiHintDo, NULL) ;
#ifdef BLIS_DEBUG
                        if (model->solver()->getIterationCount() != 0) {
                            // TODO: maybe some cuts become slack again
#ifdef BLIS_DEBUG
                            std::cout << "SLACK: resolve changed solution!"
                                      << ", iter=" 
				      << model->solver()->getIterationCount()
				      << std::endl;
#endif
                        }
			else {
#ifdef BLIS_DEBUG
			    std::cout<<"SLACK: resolve don't changed solution!"
                                     << std::endl;
#endif
			}
			
                        double tempOV = model->solver()->getObjValue();
                        double ovDiff = fabs(quality_ - tempOV);
                        
                        if (ovDiff /(1.0 + tempOV) > 1.0e-3) {
                            std::cout << "ERROR: SLACK: quality_("<<quality_ 
                                      << ") != tempOV(" << tempOV
                                      << ")" << std::endl;
                            assert(0);
                        }
                        else {
                            std::cout << "AFTER SLACK: quality_("<<quality_ 
                                      << ") == tempOV(" << tempOV
                                      << ")" << std::endl;
                        }
#endif
                    }
                    
                    delete ws;
                    delete [] newDelMark;
                    delete [] oldDelMark;
                }
            }
            
            break;
        case BlisLpStatusAbandoned:
            assert(0);
	    returnStatus = BlisReturnStatusErrLp;
            goto TERM_PROCESS;
        case BlisLpStatusDualInfeasible:
            // FIXME: maybe also primal infeasible
#ifdef BLIS_DEBUG
	    assert(0);
#endif
	    returnStatus = BlisReturnStatusUnbounded;
            goto TERM_PROCESS;
        case BlisLpStatusPrimalInfeasible:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
	    returnStatus = BlisReturnStatusInfeasible;
            goto TERM_PROCESS;
        case BlisLpStatusDualObjLim:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
	    returnStatus = BlisReturnStatusOverObjLim;
            goto TERM_PROCESS;
        case BlisLpStatusPrimalObjLim:
        case BlisLpStatusIterLim:
            /* Can't say much, need branch */
            needBranch = true;
#ifdef BLIS_DEBUG
            assert(0);
#endif
	    returnStatus = BlisReturnStatusBranch;
            goto TERM_BRANCH;
            break;
        default:
#ifdef BLIS_DEBUG
            std::cout << "PROCESS: unknown status "  <<  status << std::endl;
            assert(0);
#endif
            break;
        }

        //--------------------------------------------------
        // Call heuristics.
        //--------------------------------------------------
        
        if (keepOn && (model->heurStrategy_ != BlisHeurStrategyNone)) {
            heurObjValue = getKnowledgeBroker()->getIncumbentValue();

            for (k = 0; k < model->numHeuristics(); ++k) {
                int heurStrategy = model->heuristics(k)->strategy();
                
                if (heurStrategy != BlisHeurStrategyNone) {

                    getKnowledgeBroker()->tempTimer().start();
                    foundSolution = false;
                    foundSolution = 
                        model->heuristics(k)->searchSolution(heurObjValue,
                                                             heurSolution);
                    getKnowledgeBroker()->tempTimer().stop();

                    model->heuristics(k)->
                        addTime(getKnowledgeBroker()->tempTimer().getCpuTime());
                    model->heuristics(k)->addCalls(1);

                    if (foundSolution) {
                        ipSol = model->feasibleSolution(numIntInfs,numObjInfs);
                    }
                    
                    if (ipSol) {
                        model->heuristics(k)->addNumSolutions(1);
                        int noSols = model->heuristics(k)->noSolCalls();
                        model->heuristics(k)->addNoSolCalls(-noSols);

                        model->storeSolution(BlisSolutionTypeBounding, ipSol);
                        cutoff = model->getCutoff();
                        if (quality_ > cutoff) {
                            setStatus(AlpsNodeStatusFathomed);
                            goto TERM_PROCESS;
                        }
                    }
                    else {
                        model->heuristics(k)->addNoSolCalls(1);
                    }
                    
                } // EOF heurStrategy
            }    
        }
        
        //--------------------------------------------------
        // Generate constraints.
        //--------------------------------------------------

#if 0
        std::cout << "keepOn = " << keepOn
                  << ", geneConsHere = " << genConsHere
                  << ", numAppliedCons = " << numAppliedCons
                  << ", maxNumCons = " << maxNumCons
                  << std::endl;
#endif

        if ( keepOn && genConsHere && 
             (numAppliedCons < maxNumCons) && 
             (model->boundingPass_ < maxPass) ) {
            
            //OsiCuts newOsiCuts;
            BcpsConstraintPool newConPool;

            memcpy(currLpSolution, 
                   model->getLpSolution(),
                   numCols * sizeof(double));
            
	    // Get violated constraints that are from other processes.
            tempNumCons = newConPool.getNumConstraints();
	    getViolatedConstraints(model, currLpSolution, 
				   *(model->constraintPoolReceive()));
	    voilatedNumCons = newConPool.getNumConstraints() - tempNumCons;
	    
	    // Generate constraints (only if no violated).
	    if ((voilatedNumCons == 0) && (bS->LPSolStatus_ ==
					   MibSLPSolStatusInfeasible)) {
		lpStatus = static_cast<BlisLpStatus> 
		    (generateConstraints(model, newConPool));

		if(bS->shouldPrune_){
		    setStatus(AlpsNodeStatusFathomed);
		    quality_ = -ALPS_OBJ_MAX;
		    goto TERM_PROCESS;
		}
		    
            
		if (lpStatus != BlisLpStatusOptimal) {
		    setStatus(AlpsNodeStatusFathomed);
		    quality_ = -ALPS_OBJ_MAX; // Remove it as soon as possilbe
		    goto TERM_PROCESS;
		}
	    }
            
	    tempNumCons = newConPool.getNumConstraints();
            
            if (tempNumCons > 0) {
		// Select and install new constraints
		applyConstraints(model, currLpSolution, newConPool);
		
		// Some weak/parallel/dense constraints might be discarded.
		tempNumCons = newConPool.getNumConstraints();
		if (tempNumCons > 0) {
		    keepOn = true;
		}
		else {
		    keepOn = false;
		}

                // Move cuts from pool to array newConstraints.
                for (k = 0; k < tempNumCons; ++k) {
                    aCon = dynamic_cast<BlisConstraint *>
			(newConPool.getConstraint(k));
                    newConstraints[newNumCons++] = aCon;
                    if (newNumCons >= maxNewNumCons) {
                        // No space, need resize
#ifdef BLIS_DEBUG
                        std::cout << "NEWCUT: resize, maxNewNumCons = " 
                                  << maxNewNumCons << std::endl;
#endif
                        maxNewNumCons *= 2;
                        BcpsObject **tempNews = new BcpsObject* [maxNewNumCons];
                        memcpy(tempNews, 
                               newConstraints,
                               newNumCons * sizeof(BcpsObject *));
                        delete [] newConstraints;
                        newConstraints = tempNews;
                    }

		    // Make a copy to send pool if share
		    if (shareCon && (voilatedNumCons == 0)) {
			if (aCon->getValidRegion() == BcpsValidGlobal) {
			    model->constraintPoolSend()->
				addConstraint(new BlisConstraint(*aCon));
			}
#if 0
			std::cout << "+++ Num of send new constraint = " 
				  << model->constraintPoolSend()->getNumConstraints()
				  << std::endl;
#endif
		    }
                }
		
		newConPool.clear();
                numAppliedCons += tempNumCons;
            }
	    else { // Didn't generate any new constraints.
		keepOn = false;
	    }
        }
        else { // Don't allow to generate constraints.
            keepOn = false;
        }
    } // EOF bounding/cutting/heuristics loop
    
    //------------------------------------------------------
    // Select branching object
    //------------------------------------------------------
    
 TERM_BRANCH:

#ifdef BLIS_DEBUG_MORE
    printf("needBranch = %d\n", needBranch);
#endif
    
    if (needBranch) { 
	
        bStatus = -1;
        
        while (bStatus == -1) { 
            foundSolution = false;
            if(getKnowledgeBroker()->getProcRank() == -1) {
                std::cout << "*** I AM RANK ONE: before choose:bStatus = " 
                          << bStatus << std::endl;
            }
            bStatus = selectBranchObject(model, 
                                         foundSolution, 
                                         numPassesLeft);
            --numPassesLeft;
            
            if (bStatus == -1) { 
                lpFeasible = model->resolve();
                
                //resolved = true ;
#ifdef BLIS_DEBUG_MORE
                printf("Resolve since some col fixed, Obj value %g, numRows %d, cutoff %g\n",
                       model->solver()->getObjValue(),
                       model->solver()->getNumRows(),
                       cutoff);
#endif
                if (lpFeasible) {
		    // Update new quality.
                    setQuality(model->solver()->getObjValue());
                    if (getQuality() > cutoff) {
                        bStatus = -2;
                    }
                    // Check if feasible at the other branch due to random LP 
                    ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
                    if (ipSol) {         
                        // IP feasible
                        model->storeSolution(BlisSolutionTypeBounding, ipSol);
                        // Update cutoff
                        cutoff = model->getCutoff();
                        setStatus(AlpsNodeStatusFathomed);
                        goto TERM_PROCESS;
                    }
                }
                else {
                    // Should not happen. No, it will happen when other
                    // branch is ip feasible, and cause this branch to fathom
                    // when resolving. Test enigma.
                    // assert(0);
                    bStatus = -2;
                    setStatus(AlpsNodeStatusFathomed);
                }
            }
            
            if(getKnowledgeBroker()->getProcRank() == -1) {
                std::cout << "*** I AM RANK ONE: bStatus = " << bStatus
                          << std::endl;
            }
        }
        
        assert(bStatus != -1);
        
        //----------------------------------------------------
        // If found a branching object:
        // 1. Record basis
        // 2. Record soft var bound difference. 
        // 3. Record add/del constraints.
        // NOTE: Hard var bound differences have been recorded when branch().
        //       startXXXXX have branching bounds for this node
        //----------------------------------------------------
        
        if (bStatus >= 0) {
            
#ifdef BLIS_DEBUG_MORE
            BlisBranchObjectInt *branchObject =
                dynamic_cast<BlisBranchObjectInt *>(branchObject_);
            std::cout << "SetPregnant: branchedOn = " 
                      << model->getIntVars()[branchObject->variable()]
                      << std::endl;
#endif
            //--------------------------------------------------
            // Mark as pregnant.
            //--------------------------------------------------

            setStatus(AlpsNodeStatusPregnant);
	    
            //--------------------------------------------------
            // Save basis.
            //--------------------------------------------------

            CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
                (model->solver()->getWarmStart());
            BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc *>(desc_);
            desc->setBasis(ws);

            //----------------------------------------------
            // Save variable/constraint bound, non-core constraints
	    // and non-core variable.
	    // The sizes of hard variable bound vectors are numCols.
	    // The sizes of soft variable bound vectors are number
	    // of modified.
            //----------------------------------------------

	    int *tempVarLBPos = model->tempVarLBPos();
	    int *tempVarUBPos = model->tempVarUBPos();
	    //int *tempConLBPos = model->tempConLBPos();
	    //int *tempConUBPos = model->tempConUBPos();
	    
	    int numModSoftColLB = 0;
	    int numModSoftColUB = 0;
	    const double *currColLB = model->solver()->getColLower();
	    const double *currColUB = model->solver()->getColUpper();
	    //const double *currRowLB = model->solver()->getRowLower();
	    //const double *currRowUB = model->solver()->getRowUpper();
	    
	    double *startColLB = model->startVarLB();
	    double *startColUB = model->startVarUB();
	    //double *startRowLB = model->startConLB();
	    //double *startRowUB = model->startConUB();

	    //BlisConstraint **tempCons = NULL;


#ifdef BLIS_DEBUG_MORE
	    // Debug survived old constraints.
	    for (k = 0; k < currNumOldCons; ++k) {
		int oldPos = oldConsPos[k];
		BlisConstraint *aCon = model->oldConstraints()[oldPos];
		assert(aCon);
		std::cout << "SAVE: DBG: oldPos=" << oldPos
			  << ", k=" << k << ", len=" << aCon->getSize()
			  << ", node=" << index_ << std::endl;
	    }
#endif

	    //----------------------------------------------
	    // Decide if save explicit decription.
	    //----------------------------------------------
	    
            BlisParams * BlisPar = model->BlisPar();
            int difference = BlisPar->entry(BlisParams::difference);
    
            if (difference == -1) {
                if (depth_ % 30 == 0 || isRoot || (phase == AlpsPhaseRampup)) {
                    explicit_ = 1;
                    //std::cout << "SAVE: node "<< index_ <<" explicitly, "
                    //  << "depth=" << depth_ << std::endl;
                }
                else {
                    explicit_ = 0;
                    //std::cout << "SAVE: node "<< index_ <<" relatively, "
                    //  << "depth=" << depth_ << std::endl;
                }
            }
            else if (difference == 0) {
                explicit_ = 1;
            }
            else {
                explicit_ = 0;
            }
	    
	    //explicit_ = 1;

	    if (explicit_ || (phase == AlpsPhaseRampup) ) {
		// NOTE: full hard bound has been stored. 
		
		int index;
		int numModify = 0;
		int numSoftVarLowers = 0;
		int numSoftVarUppers = 0;
		double value;
		
		double *fVarHardLB = new double [numCols];
		double *fVarHardUB = new double [numCols];
		int *fVarHardLBInd = new int [numCols];
		int *fVarHardUBInd = new int [numCols];
		
		double *fVarSoftLB = new double [numCols];
		double *fVarSoftUB = new double [numCols];
		int *fVarSoftLBInd = new int [numCols];
		int *fVarSoftUBInd = new int [numCols];
		
		for (k = 0; k < numCols; ++k) {
		    fVarSoftLB[k] = ALPS_DBL_MAX;
		    fVarSoftUB[k] = -ALPS_DBL_MAX;
		    fVarHardLB[k] = ALPS_DBL_MAX;
		    fVarHardUB[k] = -ALPS_DBL_MAX;
		    fVarHardLBInd[k] = k;
		    fVarHardUBInd[k] = k;
		}

		//------------------------------------------
		// Build the path to an explicit node.
		//------------------------------------------
		
		AlpsTreeNode *parent = parent_;

		// First push this node since it has branching hard bounds.
		model->leafToRootPath.push_back(this);
		BlisNodeDesc* pathDesc = NULL;
		
		if (phase != AlpsPhaseRampup) {
		    while(parent) {
			model->leafToRootPath.push_back(parent);
			if (parent->getExplicit()) {
			    // Reach an explicit node, then stop.
			    break;
			}
			else {
			    parent = parent->getParent();
			}
		    }
		}

#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: EXP: path len = "<<model->leafToRootPath.size()
			  << std::endl;
#endif
		//------------------------------------------
		// Summarize bounds.
		//------------------------------------------
		
		for(j = static_cast<int> (model->leafToRootPath.size() - 1); j > -1; --j) {
		    
		    pathDesc = dynamic_cast<BlisNodeDesc*>
			((model->leafToRootPath.at(j))->getDesc());
		    
		    //--------------------------------------
		    // Full variable hard bounds.
		    //--------------------------------------
		    
		    numModify = pathDesc->getVars()->lbHard.numModify;
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->lbHard.posModify[k];
			value = pathDesc->getVars()->lbHard.entries[k];
			fVarHardLB[index] = value;
		    }
		    
		    numModify = pathDesc->getVars()->ubHard.numModify;
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->ubHard.posModify[k];
			value = pathDesc->getVars()->ubHard.entries[k];
			fVarHardUB[index] = value;
		    }
		    
		    //--------------------------------------
		    // Full variable soft bounds.
		    //--------------------------------------
		    
		    numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG_MORE
		    std::cout << "SAVE: EXP: j=" << j << ", numModify soft lb="
			      << numModify << std::endl;
#endif
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->lbSoft.posModify[k];
			value = pathDesc->getVars()->lbSoft.entries[k];
			fVarSoftLB[index] = value;
		    }
		    
		    numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef BLIS_DEBUG_MORE
		    std::cout << "SAVE: EXP: j=" << j << ", numModify soft ub="
			      << numModify << std::endl;
#endif
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->ubSoft.posModify[k];
			value = pathDesc->getVars()->ubSoft.entries[k];
			fVarSoftUB[index] = value;
		    }

		} // EOF of for(path)

		//------------------------------------------
		// Collect modified soft bounds at this node.
		// NOTE: Do this after collecting previous soft bounds.
		//------------------------------------------
		
		numModSoftColLB = 0;
		numModSoftColUB = 0;
		for (k = 0; k < numCoreCols; ++k) {
		    if (currColLB[k] != startColLB[k]) {
			fVarSoftLB[k] = currColLB[k];
			++numModSoftColLB;
#ifdef BLIS_DEBUG_MORE
			printf("Col %d, soft lb change, start %g, curr %g\n",
			       k, startColLB[k], currColLB[k]);
#endif
			
		    }
		    if (currColUB[k] != startColUB[k]) {
			fVarSoftUB[k] = currColUB[k];
			++numModSoftColUB;
		    }
		}
		
#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: EXP: THIS: numModSoftColLB = "<<numModSoftColLB 
			  << ", numModSoftColUB = " << numModSoftColUB << std::endl;
#endif
		
		//--------------------------------------
		// Debug if bounds are consistant.
		//--------------------------------------

#ifdef BLIS_DEBUG
		for (k = 0; k < numCols; ++k) {
		    
		    //std::cout << "EXP: COL[" << k <<"]: " 
		    //      <<"hardLB=" << fVarHardLB[k] 
		    //      <<", hardUB=" << fVarHardUB[k] << std::endl;
		    
		    // Hard lower bound should not greater than
		    // hard upper bound.
		    if (fVarHardLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
			printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
			       k, fVarHardLB[k], fVarHardUB[k]);
			assert(0);
		    }
		    
		    if (fVarSoftLB[k] < ALPS_BND_MAX) {
			// Soft lower changed, and should not greater 
			// than its hard upper bound.
			if (fVarSoftLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
			    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
				   k, fVarSoftLB[k], fVarHardUB[k]);
			    assert(0);	
			}
		    }
		    
		    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
			// Soft upper changed, and should not less 
			// than its hard lower bound.
			if (fVarSoftUB[k] < fVarHardLB[k] - ALPS_GEN_TOL) {
			    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
				   k, fVarHardLB[k], fVarSoftUB[k]);
			    assert(0);
			}
		    }
		}
#endif
		
		//------------------------------------------
		// Record hard variable bounds. FULL set.
		//------------------------------------------
		
		desc->assignVarHardBound(numCols,
					 fVarHardLBInd,
					 fVarHardLB,
					 numCols,
					 fVarHardUBInd,
					 fVarHardUB);
		
		//------------------------------------------
		// Recode soft variable bound. Modified.
		//------------------------------------------
		
		for (k = 0; k < numCols; ++k) {
		    if (fVarSoftLB[k] < ALPS_BND_MAX) {
			fVarSoftLBInd[numSoftVarLowers] = k;
			fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
		    }
		    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
			fVarSoftUBInd[numSoftVarUppers] = k;
			fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
		    }
		}

		
#ifdef BLIS_DEBUG_MORE
		// Print soft bounds.
		std::cout << "SAVE: EXP: numSoftVarLowers=" << numSoftVarLowers
			  << ", numSoftVarUppers=" << numSoftVarUppers
			  << std::endl;
		for (k = 0; k < numSoftVarLowers; ++k) {
		    std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
			      << fVarSoftLB[k] << std::endl;		    
		}
		std::cout << "------------------" << std::endl;
		for (k = 0; k < numSoftVarUppers; ++k) {
		    std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
			      << fVarSoftUB[k] << std::endl;
		}
		std::cout << "------------------" << std::endl << std::endl;
#endif

		//if ( (numSoftVarUppers > 0) || (numSoftVarLowers > 0) ) {
                
                // Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
                desc->assignVarSoftBound(numSoftVarLowers, 
                                         fVarSoftLBInd,
                                         fVarSoftLB,
                                         numSoftVarUppers, 
                                         fVarSoftUBInd,
                                         fVarSoftUB);

                //------------------------------------------
                // Full set of active non-core constraints.
                //------------------------------------------
                // old constraint: model->oldConstraints_, currNumOldCons.
                // new constraint: newConstraints, newNumCons.

                BcpsObject **toAddCons = new BcpsObject * [currNumOldCons + 
                                                           newNumCons];
                if (currNumOldCons > 0) {
                    // Hard copy of the survived old constraints.
                    for (k = 0; k < currNumOldCons; ++k) {
			int oldPos = oldConsPos[k];
                        BlisConstraint *aCon = model->oldConstraints()[oldPos];
                        assert(aCon);
#ifdef BLIS_DEBUG
			std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons 
				  << ", k=" << k << ", len=" << aCon->getSize()
                                  << ", node=" << index_ << std::endl;
#endif
			
                        BlisConstraint *newCon = new BlisConstraint(*aCon);
                        toAddCons[k] = newCon;
                    }
                }
                
		if (newNumCons > 0) {
		    // Include new constraints.
		    memcpy(toAddCons + currNumOldCons, 
			   newConstraints,
			   newNumCons * sizeof(BcpsObject *));
		}

		//------------------------------------------
                // Save in description. Add first delete exiting, then add.
		// Pointers in model->oldConstraints_ can be dangling.
		// It is not safe to use model->oldConstraints_ after adding.
		
		// If this node is the root of a subtree, and before processing
		// it has a list of cuts, then model->oldConstraints_ 
		// stores pointer to the cuts when installing.
		
		// Need update model->constraints_ here OR do not be smart!
		// 1/6/06: Choose to be dumn.
		//------------------------------------------

		//------------------------------------------
		// Generating constraints,
		// also means that slack ones might been removed.
		//------------------------------------------

		int numTotal = currNumOldCons + newNumCons;
                desc->setAddedConstraints(numTotal, toAddCons);
                
#ifdef BLIS_DEBUG
		std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons
			  << ", newNumCons=" << newNumCons
			  << std::endl;
#endif	
		
                //------------------------------------------
                // Full set of active non-core variables.
                //------------------------------------------
                // Todo.
                
		//------------------------------------------
		// Clear path vector.
		//------------------------------------------
		
		model->leafToRootPath.clear();
		assert(model->leafToRootPath.size() == 0);
	    }
	    else { // Relative.


	      //THIS CODE IS WHAT IS RUNNING NOW

		//------------------------------------------
		// Record soft bound changes for core vars.
		//------------------------------------------

		// Variable bound change 
		numModSoftColLB = 0;
		numModSoftColUB = 0;
		for (k = 0; k < numCoreCols; ++k) {
		    if (currColLB[k] != startColLB[k]) {
			tempVarLBPos[numModSoftColLB] = k;
			/* startColLB as a temporary storage vector */
			startColLB[numModSoftColLB] = currColLB[k];
			
#ifdef BLIS_DEBUG_MORE
			printf("Col %d, soft lb change, start %g, curr %g\n",
			       k, startColLB[k], currColLB[k]);
#endif
			
			++numModSoftColLB;
		    }
		    if (currColUB[k] != startColUB[k]) {
			tempVarUBPos[numModSoftColUB] = k;
			startColUB[numModSoftColUB] = currColUB[k];
			++numModSoftColUB;
		    }
		}

#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: REL: numModSoftColLB = " 
                          << numModSoftColLB 
			  << ", numModSoftColUB = " 
                          << numModSoftColUB 
			  << std::endl;
#endif
            
		if (numModSoftColLB > 0 || numModSoftColUB > 0) {
#ifdef BLIS_DEBUG
		    //assert(0);
#endif
		    desc->setVarSoftBound(numModSoftColLB, 
					  tempVarLBPos,
					  startColLB,
					  numModSoftColUB, 
					  tempVarUBPos,
					  startColUB);
		}
            
		//------------------------------------------
		// TODO: Constraint bounds change.
		//------------------------------------------

#if 0
		for (k = 0; k < numCoreRows; ++k) {
		    if (currRowLB[k] != startRowLB[k]) {
			tempConLBPos[numModSoftRowLB] = k;
			startRowLB[numModSoftRowLB] = currRowLB[k];
			++numModSoftRowLB;
		    }
		    if (currRowUB[k] != startRowUB[k]) {
			tempConUBPos[numModSoftRowUB] = k;
			startRowUB[numModSoftRowUB] = currRowUB[k];
			++numModSoftRowUB;
		    }
		}
		if (numModSoftRowLB > 0 || numModSoftRowUB > 0) {
		    desc->setConSoftBound(numModSoftRowLB, 
					  tempConLBPos,
					  startRowLB,
					  numModSoftRowUB, 
					  tempConUBPos,
					  startRowUB);
		}
#endif	
            
		if (genConsHere) {
		    // NOTE: explicit_ can NOT do this if, since genConsHere maybe
		    //       false here, but there are maybe cons from parents.

		    //--------------------------------------
		    // Record add constraints.
		    //--------------------------------------

		    // Apend will copy old, then add new. 
		    // If this node has a list of cuts before pointers in 
		    // model->oldConstraints() will be kept. Safe!
		    if (newNumCons > 0) {
			desc->appendAddedConstraints(newNumCons,
                                                     newConstraints);
                    }
		    
		    //--------------------------------------
		    // Record deleted constraint positions.
		    //--------------------------------------
		    
		    int *oldLeft = new int [origNumOldCons];
		    int leftCon;
		    CoinZeroN(oldLeft, origNumOldCons);
		    
		    for (k = 0; k < currNumOldCons; ++k) {
			leftCon = oldConsPos[k];
			assert(leftCon >= 0 && leftCon < origNumOldCons);
			oldLeft[leftCon] = 1;
		    }
		    // 
		    leftCon = 0;
		    for (k = 0; k < origNumOldCons; ++k) {
			if (oldLeft[k] == 0) {
			    // Deleted. Now oldLeft stores delete position.
			    oldLeft[leftCon++] = k;
			}
                        //FIXME: clean
                        //assert(k < 15196924);
		    }
		    desc->delConstraints(leftCon, oldLeft);
		    
#ifdef BLIS_DEBUG
		    std::cout << "PROCESS: ADD: new cuts=" << newNumCons 
			      << ", numRows=" << model->solver()->getNumRows() 
			      << ", numStartRows="<< numStartRows 
                              << ", origNumStartRows="<< origNumStartRows 
                              << ", num removed=" << leftCon << std::endl;                    
#endif
		    
		}// EOF of if(genConsHere)
	    } // EOF of relative
        }
        else if (bStatus == -2) {
            
#if 0
            std::cout << "bStatus = -2, fathom this node!" << std::endl;
#endif
            //branchObject->getDown()[0], branchObject->getDown()[1]);
	    
            setStatus(AlpsNodeStatusFathomed);
        }
        else {
            throw CoinError("No branch object found", "process", 
                            "BlisTreeNode");
        }
    }
    
    //------------------------------------------------------
    // End of process()
    //------------------------------------------------------
    
 TERM_PROCESS:

    if(0)
      model->solver()->writeLp("termprocess");

    bool printCutStat = false;
    if (genConsHere) {
        if ( (getKnowledgeBroker()->getProcType() == AlpsProcessTypeMaster) &&
             (msgLevel > 0) ) {
            printCutStat = true;
        }
        else if ( (getKnowledgeBroker()->getProcType()==AlpsProcessTypeHub) &&
                  (hubMsgLevel > 0) ) {
            printCutStat = true;
        }
        else if ((getKnowledgeBroker()->getProcType()==AlpsProcessTypeWorker)&&
		 (workerMsgLevel > 0)) {
            printCutStat = true;
        }
    }

    if (printCutStat == true) {
        printCutStat = false;
        if ( (msgLevel > 100) || (index_ == 0) ) {
            printCutStat = true;
        }
    }

    if (printCutStat) {
        int numT = model->numCutGenerators();
        for (k = 0; k < numT; ++k) {
            if ( model->cutGenerators(k)->calls() > -1) {
                model->blisMessageHandler()->message(BLIS_CUT_STAT_NODE,
                                                     model->blisMessages())
                    << index_
                    << model->cutGenerators(k)->name()
                    << model->cutGenerators(k)->calls()
                    << model->cutGenerators(k)->numConsGenerated()
                    << model->cutGenerators(k)->time()
                    << model->cutGenerators(k)->strategy()
                    << CoinMessageEol;
            }   
        }

        numT = model->numHeuristics();
        for (k = 0; k < numT; ++k) {
            if ( model->heuristics(k)->calls() > -1) {
                model->blisMessageHandler()->message(BLIS_HEUR_STAT_NODE,
                                                     model->blisMessages())
                    << index_
                    << model->heuristics(k)->name()
                    << model->heuristics(k)->calls()
                    << model->heuristics(k)->numSolutions()
                    << model->heuristics(k)->time()
                    << model->heuristics(k)->strategy()
                    << CoinMessageEol;
            }   
        }
    }

    delete [] heurSolution;
    delete [] currLpSolution;

    if (status_ == AlpsNodeStatusFathomed) {
	// Delete new cuts since no use anymore.
        for (k = 0; k < newNumCons; ++k) {
            delete newConstraints[k];
        }
    }
    delete [] newConstraints;
    delete [] oldConsPos;
    
    model->isRoot_ = false;

#ifdef BLIS_DEBUG_MORE
    // Debug survived old constraints.
    //int currNumOldCons = model->getNumOldConstraints();
    for (k = 0; k < currNumOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	std::cout << "SAVE: DBG: TERM: "
		  << "k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
    }
#endif
    
    return returnStatus;
}

//#############################################################################
