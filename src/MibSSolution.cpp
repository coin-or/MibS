/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*          Sahar Tahernajad, Lehigh University                              */
/*                                                                           */
/* Copyright (C) 2007-2023 Lehigh University, Scott DeNegre, and Ted Ralphs. */
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

   std::string inputFormat(localModel_->MibSPar_->entry
			    (MibSParams::inputFormat));
   int isInterdict = localModel_->MibSPar_->entry(MibSParams::bilevelProblemType);
   double nearInt = 0.0;
   //int size(localModel_->getNumOrigVars());
   int j(0); 

   int uN(localModel_->getUpperDim());
   int lN(localModel_->getLowerDim());
   int index(0);
   int * upperColInd = localModel_->getUpperColInd();
   int * lowerColInd = localModel_->getLowerColInd();
   std::string * columnName = localModel_->getColumnName();

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

   if (isInterdict){
      os << "First stage (interdiction) variable values:" << std::endl;
   }else{
      os << "First stage (upper level) variable values:" << std::endl;
   }
   for(j = 0; j < uN; ++j) {
      index = upperColInd[j];      
      if (values_[index] > 1.0e-15 || values_[index] < -1.0e-15) {
	 nearInt = floor(values_[index] + 0.5);
	 if (ALPS_FABS(nearInt - values_[index]) < 1.0e-6) {
	     if (!columnName) {
		 os << "x[" << j << "] = " << nearInt << std::endl;
	     }
	     else {
                if (isInterdict){
                   os << columnName[index-uN] << " = " << nearInt << std::endl;
                }else{
                   os << columnName[index] << " = " << nearInt << std::endl;
                }
	     }
	 }
	 else {
	     if (!columnName) {
		 os << "x[" << j << "] = " << values_[index] << std::endl;
	     }
	     else {
		 os << columnName[index] << " = " << values_[index] << std::endl;
	     }
	 }
      }
   }

   os << "Second stage (lower level) variable values:" << std::endl;
   for(j = 0; j < lN; ++j) {
      index = lowerColInd[j];      
      if (values_[index] > 1.0e-15 || values_[index] < -1.0e-15) {
	 nearInt = floor(values_[index] + 0.5);
	 if (ALPS_FABS(nearInt - values_[index]) < 1.0e-6) {
	     if (!columnName) {
		 os << "y[" << j << "] = " << nearInt << std::endl;
	     }
	     else {
		 os << columnName[index] << " = " << nearInt << std::endl;
	     }
	 }
	 else {
             if(!columnName){
		 os << "y[" << j << "] = " << values_[index] << std::endl;
	     }
	     else{
		 os << columnName[index] << " = " << values_[index] << std::endl;
	     }
	 }
      }
   }

  int useBendersInterdictionCut = 
     localModel_->MibSPar_->entry(MibSParams::useBendersInterdictionCut);
  int useImprovingSolutionIC =
     localModel_->MibSPar_->entry(MibSParams::useImprovingSolutionIC);
  int useImprovingDirectionIC =
     localModel_->MibSPar_->entry(MibSParams::useImprovingDirectionIC);
  int useHypercubeIC =
     localModel_->MibSPar_->entry(MibSParams::useHypercubeIC);
  int useGeneralizedNoGoodCut = 
     localModel_->MibSPar_->entry(MibSParams::useGeneralizedNoGoodCut);
  int useBendersBinaryCut
     = localModel_->MibSPar_->entry(MibSParams::useBendersBinaryCut);
  int useIntegerNoGoodCut
     = localModel_->MibSPar_->entry(MibSParams::useIntegerNoGoodCut);

   os << "Number of problems (VF) solved = " << localModel_->counterVF_
      << std::endl
      << "Number of problems (UB) solved = " << localModel_->counterUB_
      << std::endl
      << "Time for solving problem (VF) = " << localModel_->timerVF_
      << std::endl
      << "Time for solving problem (UB) = " << localModel_->timerUB_
      << std::endl;
   if (useBendersInterdictionCut){
      os << "Number of Benders Interdiction Cuts Generated: "
         << localModel_->counterBendersInterdict_ << std::endl;
   }
   if (useHypercubeIC){
      os << "Number of Hypercube Intersection Cuts Generated: "
         << localModel_->counterHypercubeIC_ << std::endl;
   }
   if (useGeneralizedNoGoodCut){
      os << "Number of Generalized No Good Cuts Generated: "
         << localModel_->counterGeneralizedNoGood_ << std::endl;
   }
   if (useBendersBinaryCut){
      os << "Number of Benders Binary Cuts Generated: "
         << localModel_->counterBendersBinary_ << std::endl;
   }
   if (useIntegerNoGoodCut){
      os << "Number of Integer No Good Cuts Generated: "
         << localModel_->counterIntegerNoGood_ << std::endl;
   }
   if (useImprovingDirectionIC){
      os << "Number of IDICs Generated:" << std::endl
         << "   Full Int IDIC:                " << localModel_->counterXYIntIDIC_
         << std::endl
         << "   Linking Int IDIC:             "  << localModel_->counterLIntIDIC_
         << std::endl
         << "   Second-level Int IDIC:        "  << localModel_->counterYIntIDIC_
         << std::endl
         << "   Fractional IDIC:              "  << localModel_->counterFracIDIC_
         << std::endl
         << "Number of IDIC Generation Failures:"  << std::endl
         << "   Full Int IDIC (Fail):         "  << localModel_->counterXYIntIDICFail_
         << std::endl
         << "   Linking Int IDIC (Fail):      "  << localModel_->counterLIntIDICFail_
         << std::endl
         << "   Second-level Int IDIC (Fail): "  << localModel_->counterYIntIDICFail_
         << std::endl
         << "   Fractional IDIC (Fail):       "  << localModel_->counterFracIDICFail_
         << std::endl;
   }
   if (useImprovingSolutionIC){
      os << "Number of ISICs Generated:" << std::endl
         << "   Full Int ISIC:                " << localModel_->counterXYIntISIC_
         << std::endl
         << "   Linking Int ISIC:             "  << localModel_->counterLIntISIC_
         << std::endl
         << "   Second-level Int ISIC:        "  << localModel_->counterYIntISIC_
         << std::endl
         << "   Fractional:                   "  << localModel_->counterFracISIC_
         << std::endl
         << "Number of ISIC Generation Failures:"  << std::endl
         << "   Full Int ISIC (Fail):         "  << localModel_->counterXYIntISICFail_
         << std::endl
         << "   Linking Int ISIC (Fail):      "  << localModel_->counterLIntISICFail_
         << std::endl
         << "   Second-level Int ISIC (Fail): "  << localModel_->counterYIntISICFail_
         << std::endl
         << "   Fractional ISIC (Fail):       "  << localModel_->counterFracISICFail_
         << std::endl;
   }
}

