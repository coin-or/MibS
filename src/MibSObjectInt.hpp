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

#ifndef MibSObjectInt_h_
#define MibSObjectInt_h_

#include "BlisObjectInt.h"
#include "BcpsModel.h"

//#############################################################################

class MibSObjectInt : public BlisObjectInt {

 private:
   
  
 public:
  
  MibSObjectInt()
     {
     }
  
  ~MibSObjectInt();

  virtual BcpsBranchObject * createBranchObject(BcpsModel *m,
						int direction) const;

};

#endif
