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

#ifndef MibSBranchObjectInt_h_
#define MibSBranchObjectInt_h_

#include "BlisBranchObjectInt.h"
#include "BcpsModel.h"

//#############################################################################

class MibSBranchObjectInt : public BlisBranchObjectInt {

 private:
   
  
 public:
  
   //   MibSBranchObjectInt();

  ~MibSBranchObjectInt();



  /** Default constructor. */
  MibSBranchObjectInt()    
    : 
    BlisBranchObjectInt()
    {
      type_ = BlisBranchingObjectTypeInt;
      down_[0] = 0.0;
      down_[1] = 0.0;
      up_[0] = 0.0;
      up_[1] = 0.0;
      
    }

  /** Construct a branching object, which branching on variable varInd. 
      \param varInd     the index of integer variable in object set
      \param direction  the direction of first branching: 1(up), -1(down)
      \param value      the fractional solution value of variable varInd 
  */
  MibSBranchObjectInt(BlisModel * model, 
		      int varInd,
		      int direction,
		      double value)
     : 
     BlisBranchObjectInt(model, varInd, direction, value)
     {
	type_ = BlisBranchingObjectTypeInt;
	int iColumn = model->getIntColIndices()[objectIndex_];
	down_[0] = model->solver()->getColLower()[iColumn];
	down_[1] = floor(value_);
	up_[0] = ceil(value_);
	up_[1] = model->getColUpper()[iColumn];
     }
  
  /** Construct a branching object, which branching on variable varInd. 
      \param varInd     the index of integer variable in object set
      \param intScore   the integer score/goodness
      \param dblScore   the double score/goodness
      \param direction  the direction of first branching: 1(up), -1(down)
      \param value      the fractional solution value of variable varInd 
  */
  MibSBranchObjectInt(BlisModel * model,
		      int varInd,
		      int intScore,
		      double dblScore,
		      int direction,
		      double value)
     : 
     BlisBranchObjectInt(model, varInd, intScore, dblScore, direction, value)
     {
	type_ = BlisBranchingObjectTypeInt;
	int iColumn = model->getIntColIndices()[objectIndex_];
	down_[0] = model->solver()->getColLower()[iColumn];
	down_[1] = floor(value_);
	up_[0] = ceil(value_);
	up_[1] = model->getColUpper()[iColumn];
     }
    
  /** Create a degenerate branching object.
      Specifies a `one-direction branch'. Calling branch() for this 
      object will always result in lowerValue <= x <= upperValue. 
      Used to fix a variable when lowerValue = upperValue.
  */
  MibSBranchObjectInt(BlisModel * model,
		      int varInd, 
		      int direction,
		      double lowerValue, 
		      double upperValue)
     :
     BlisBranchObjectInt(model, varInd, direction, lowerValue)
     {
	type_ = BlisBranchingObjectTypeInt;
	numBranchesLeft_ = 1;
	down_[0] = lowerValue;
	down_[1] = upperValue;
	up_[0] = lowerValue;
	up_[1] = upperValue;
     }
  

  
  
};

#endif
