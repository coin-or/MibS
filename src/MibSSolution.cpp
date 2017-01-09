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

#include <cmath>

#include "MibSSolution.hpp"

//#############################################################################
MibSSolution::MibSSolution()
   : BlisSolution()
{

   localModel_ = NULL;

}

//#############################################################################
MibSSolution::MibSSolution(int s, 
			   const double *values, 
			   double objVal, 
			   MibSModel *mibs)
   : BlisSolution(s, values, objVal)
{

  /** Print out the solution **/
   
   if(!mibs) return;
   
   localModel_ = mibs;

   int msgLevel = mibs->AlpsPar()->entry(AlpsParams::msgLevel);
   int uN(mibs->getUpperDim());
   int lN(mibs->getLowerDim());
   int index(0);
   double etol(mibs->getTolerance());
   int * upperColInd = mibs->getUpperColInd();
   int * lowerColInd = mibs->getLowerColInd();

   if (msgLevel > 5){
      std::cout << std::endl << "Bilevel feasible solution found:" 
		<< std::endl;
      std::cout << "Upper-level objective value: " << objVal 
		<< std::endl << std::endl;

      int i(0);
   
      for(i = 0; i < uN; i++){
	 index = upperColInd[i];
	 if((values[index] > etol) || (values[index] < - etol)) 
	    std::cout << "UL[" << i << "]: " << values[index] << std::endl;
      }
      
      for(i = 0; i < lN; i++){
	 index = lowerColInd[i];
	 if((values[index] > etol) || (values[index] < - etol)) 
	    std::cout << "LL[" << i << "]: " << values[index] << std::endl;
      }
   }

}

//#############################################################################
MibSSolution::~MibSSolution()
{

}

//#############################################################################
void
MibSSolution::print(std::ostream& os) const 
{
   
   double nearInt = 0.0;
   //int size(localModel_->getNumOrigVars());
   int j(0); 


   int uN(localModel_->getUpperDim());
   int lN(localModel_->getLowerDim());
   int index(0);
   int * upperColInd = localModel_->getUpperColInd();
   int * lowerColInd = localModel_->getLowerColInd();

   /*   
   for(j = 0; j < size; ++j) {
      if (values_[j] > 1.0e-15 || values_[j] < -1.0e-15) {
	 nearInt = floor(values_[j] + 0.5);
	 if (ALPS_FABS(nearInt - values_[j]) < 1.0e-6) {
	    os << "x[" << j << "] = " << nearInt << std::endl;
	 }
	 else {
	    os << "x[" << j << "] = " << values_[j] << std::endl;
	 }
      }
   }
   */

   for(j = 0; j < uN; ++j) {
      index = upperColInd[j];      
      if (values_[index] > 1.0e-15 || values_[index] < -1.0e-15) {
	 nearInt = floor(values_[index] + 0.5);
	 if (ALPS_FABS(nearInt - values_[index]) < 1.0e-6) {
	    os << "x[" << j << "] = " << nearInt << std::endl;
	 }
	 else {
	    os << "x[" << j << "] = " << values_[index] << std::endl;
	 }
      }
   }

   for(j = 0; j < lN; ++j) {
      index = lowerColInd[j];      
      if (values_[index] > 1.0e-15 || values_[index] < -1.0e-15) {
	 nearInt = floor(values_[index] + 0.5);
	 if (ALPS_FABS(nearInt - values_[index]) < 1.0e-6) {
	    os << "y[" << j << "] = " << nearInt << std::endl;
	 }
	 else {
	    os << "y[" << j << "] = " << values_[index] << std::endl;
	 }
      }
   }


}

