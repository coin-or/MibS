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

#ifndef MibSBranchStrategyStrong_h_
#define MibSBranchStrategyStrong_h_

#include "BlisBranchStrategyStrong.h"
#include "BlisModel.h"

//#############################################################################

class MibSBranchStrategyStrong : public BlisBranchStrategyStrong {

 private:
   
  
 public:
  
  MibSBranchStrategyStrong();

  MibSBranchStrategyStrong(BlisModel *model);

  ~MibSBranchStrategyStrong();

  int createCandBranchObjects(int numPassesLeft, double ub);

};

#endif
