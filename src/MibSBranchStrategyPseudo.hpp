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

#ifndef MibSBranchStrategyPseudo_h_
#define MibSBranchStrategyPseudo_h_

#include "BlisBranchStrategyPseudo.h"
#include "BlisModel.h"

//#############################################################################

class MibSBranchStrategyPseudo : public BlisBranchStrategyPseudo {

 private:
   
  
 public:
  
  MibSBranchStrategyPseudo();

  MibSBranchStrategyPseudo(BlisModel *model, int rel);

  ~MibSBranchStrategyPseudo();

  int createCandBranchObjects(int numPassesLeft, double ub);

};

#endif
