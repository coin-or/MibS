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

#include "MibSBranchStrategyMaxInf.hpp"
#include "MibSModel.hpp"
#include "MibSObjectInt.hpp"

//#include "BlisObjectInt.h"
//#include "BlisModel.h"

//#############################################################################
MibSBranchStrategyMaxInf::MibSBranchStrategyMaxInf()
   : BlisBranchStrategyMaxInf()
{


}

//#############################################################################
MibSBranchStrategyMaxInf::MibSBranchStrategyMaxInf(BlisModel *model)
   : BlisBranchStrategyMaxInf(model)
{

}

//#############################################################################
MibSBranchStrategyMaxInf::~MibSBranchStrategyMaxInf()
{

}

//#############################################################################

/** Create a set of candidate branching objects. */
int 
MibSBranchStrategyMaxInf::createCandBranchObjects(int numPassesLeft, double ub)
{

    int numInfs = 0;
    
    int i, col, preferDir, maxInfDir, maxScoreDir;
    
    double score, maxScore = 0.0;
    double infeasibility, maxInf = 0.0;
    
    BlisModel *model = dynamic_cast<BlisModel *>(model_);
    
    BlisObjectInt * intObject = 0;
    BlisObjectInt * maxInfIntObject = 0;
    BlisObjectInt * maxScoreIntObject = 0;

    //MibSObjectInt * intObject = 0;
    //MibSObjectInt * maxInfIntObject = 0;
    //MibSObjectInt * maxScoreIntObject = 0;
    
    int numObjects = model->numObjects();
    
    double *objCoef = model->getObjCoef();

    MibSModel *mibsmodel = dynamic_cast<MibSModel *>(model);
    bool bilevelBranch = mibsmodel->bS_->useBilevelBranching_;

    if(bilevelBranch){    
       //create branch based on bilevel infeasibility
       
       std::cout << "Using Bilevel Branching." << std::endl;



    }
    else{
       //use Blis MaxInf branching
       
       for (i = 0; i < numObjects; ++i) {
	  
	  // TODO: currently all integer object.
	  intObject = dynamic_cast<BlisObjectInt *>(model->objects(i));
	  //intObject = dynamic_cast<MibSObjectInt *>(model->objects(i));
	  infeasibility = intObject->infeasibility(model, preferDir);
	  
	  if (infeasibility) {
	     ++numInfs;
	  
	     if (infeasibility > maxInf) {
		maxInfIntObject = intObject;
		maxInfDir = preferDir;
		maxInf = infeasibility;
	     }
	     
	     col = intObject->columnIndex();
	     score = ALPS_FABS(objCoef[col] * infeasibility);
	     
	     if (score > maxScore) {
		maxScoreIntObject = intObject;
		maxScoreDir = preferDir;
		maxScore = score;
	     }
	  }
       }
       
       assert(numInfs > 0);
       
       if (maxScoreIntObject) {
	  maxInfIntObject = maxScoreIntObject;
	  maxInfDir = maxScoreDir;
       }
    }


    numBranchObjects_ = 1;

    //FIXME: THINK I NEED TO DERIVE MY OWN BRANCHING OBJECT CLASS
    branchObjects_ = new BcpsBranchObject* [1];
    branchObjects_[0] = maxInfIntObject->createBranchObject(model,
                                                            maxInfDir);
    
    return 0;
}
