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

#ifndef MibSSolution_h_
#define MibSSolution_h_

#include "BlisSolution.h"

#include "MibSModel.hpp"

//#############################################################################

class MibSSolution : public BlisSolution {

 private:
   
   MibSModel * localModel_;

 public:
   
   MibSSolution();
   MibSSolution(int s,
		const double *values,
		double objVal,
		MibSModel *mibs=0);
   ~MibSSolution();

   void print(std::ostream& os) const;
	      
};

#endif
