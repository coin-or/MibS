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

#ifndef MibSBranchStrategyMaxInf_h_
#define MibSBranchStrategyMaxInf_h_

#include "BlisBranchStrategyMaxInf.h"
#include "BlisModel.h"

//#############################################################################

class MibSBranchStrategyMaxInf : public BlisBranchStrategyMaxInf {

 private:
   
  
 public:
  
  MibSBranchStrategyMaxInf();

  MibSBranchStrategyMaxInf(BlisModel *model);

  ~MibSBranchStrategyMaxInf();

  int createCandBranchObjects(int numPassesLeft, double ub);

};

#endif
