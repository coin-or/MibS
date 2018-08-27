/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*                                                                           */
/* Copyright (C) 2007-2015 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#include "BlisModel.h"
#include "BlisConstraint.h"
#include "BlisVariable.h"
#include "BlisNodeDesc.h"
#include "BlisConfig.h"

#include "BcpsConfig.h"

#include "BlisBranchStrategyMaxInf.h"
#include "BlisBranchStrategyPseudo.h"
#include "BlisBranchStrategyRel.h"
#include "BlisBranchStrategyStrong.h"

#include "BlisHeurRound.h"

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"

#include "OsiCbcSolverInterface.hpp"

#include "MibSModel.hpp"
#include "MibSSolution.hpp"
#include "MibSCutGenerator.hpp"
#include "MibSBilevel.hpp"
#include "MibSTreeNode.hpp"
#include "MibSConstants.hpp"	

#include "MibSBranchStrategyMaxInf.hpp"
#include "MibSBranchStrategyPseudo.hpp"
#include "MibSBranchStrategyStrong.hpp"

//FIXME::RELIABILITY BRANCHING DOESNT WORK
//NECESSARY DATA MEMBERS ARE DESIGNATED AS PRIVATE
//IN PARENT CLASS.  DIDNT WANT TO ALTER BLIS CODE
//#include "MibSBranchStrategyRel.h"
#include "MibSBranchStrategyStrong.hpp"

//#############################################################################
MibSModel::MibSModel()
{
  initialize();
}

//#############################################################################
MibSModel::~MibSModel()
{
    
  if(upperColInd_) delete [] upperColInd_;
  if(lowerColInd_) delete [] lowerColInd_;
  if(upperRowInd_) delete [] upperRowInd_;
  if(origUpperRowInd_) delete [] origUpperRowInd_;
  if(lowerRowInd_) delete [] lowerRowInd_;
  if(structRowInd_) delete [] structRowInd_;
  if(fixedInd_) delete [] fixedInd_;
  if(interdictCost_) delete [] interdictCost_;
  if(origColLb_) delete [] origColLb_;
  if(origColUb_) delete [] origColUb_;
  if(lColLbInLProb_) delete [] lColLbInLProb_;
  if(lColUbInLProb_) delete [] lColUbInLProb_;
  if(origRowLb_) delete [] origRowLb_;
  if(origRowUb_) delete [] origRowUb_;
  if(origRowSense_) delete [] origRowSense_;
  if(lowerObjCoeffs_) delete [] lowerObjCoeffs_;
  if(columnName_) delete [] columnName_;
  if(rowName_) delete [] rowName_;
  if(MibSPar_) delete MibSPar_;
  if(lowerConstCoefMatrix_) delete lowerConstCoefMatrix_;
  if(A2Matrix_) delete A2Matrix_;
  if(G2Matrix_) delete G2Matrix_;
  if(bS_) delete bS_;
    
}

//#############################################################################
void 
MibSModel::initialize()
{

  llDataFile_ = "";
  ulDataFile_ = "";
  ulAmplModelFile_ = "";
  ulAmplDataFile_ = "";
  etol_ = 1e-5;
  numVars_ = 0;
  numOrigVars_ = 0;
  numCons_ = 0;
  numOrigCons_ = 0;
  objSense_ = 0.0;
  lowerDim_ = 0;
  lowerObjSense_ = 0.0;
  upperDim_ = 0;
  leftSlope_ = 0;
  rightSlope_ = 0;
  lowerRowNum_ = 0;
  upperRowNum_ = 0;
  origUpperRowNum_ =0;
  structRowNum_ = 0;
  sizeFixedInd_ = 0;
  counterVF_ = 0;
  counterUB_ = 0;
  timerVF_ = 0.0;
  timerUB_ = 0.0;
  countIteration_ = 0;
  isInterdict_ = false;
  isPureInteger_ = true;
  isUpperCoeffInt_ = true;
  isLowerCoeffInt_ = true;
  allUpperBin_ = true;
  allLowerBin_ = true;
  allLinkingBin_ = true;
  positiveA1_ = true;
  positiveA2_ = true;
  positiveG1_ = true;
  positiveG2_ = true;
  upperColInd_ = NULL;
  lowerColInd_ = NULL;
  upperRowInd_ = NULL;
  origUpperRowInd_ = NULL;
  lowerRowInd_ = NULL;
  structRowInd_ = NULL;
  fixedInd_ = NULL;
  origColLb_ = NULL;
  origColUb_ = NULL;
  lColLbInLProb_ = NULL;
  lColUbInLProb_ = NULL;
  origRowLb_ = NULL;
  origRowUb_ = NULL;
  origRowSense_ = NULL;
  lowerObjCoeffs_ = NULL;
  columnName_ = NULL;
  rowName_ = NULL;
  interdictCost_ = NULL;
  origConstCoefMatrix_ = NULL;
  lowerConstCoefMatrix_ = NULL;
  A2Matrix_ = NULL;
  G2Matrix_ = NULL;
  bS_ = new MibSBilevel();
  //simpleCutOnly_ = true; //FIXME: should make this a parameter
  //bindingMethod_ = "BLAND"; //FIXME: should make this a parameter
  //bindingMethod_ = "BASIS"; //FIXME: should make this a parameter
  MibSPar_ = new MibSParams;
  //maxAuxCols_ = 0; //FIXME: should make this a parameter

  MibSCutGenerator *cg = new MibSCutGenerator(this);

  cg->setStrategy(BlisCutStrategyPeriodic);
  cg->setCutGenerationFreq(1);  // Generate cuts at every node

  addCutGenerator(cg);

  setBlisParameters();
}

//#############################################################################

/** Read in parameters. */
void 
MibSModel::readParameters(const int argnum, const char * const * arglist)
{
     //std::cout << "Reading parameters ..." << std::endl;
    AlpsPar_->readFromArglist(argnum, arglist);
    BlisPar_->readFromArglist(argnum, arglist);
    MibSPar_->readFromArglist(argnum, arglist);
}

//#############################################################################
void 
MibSModel::readInstance(const char* dataFile)
{
   
#if 0
   std::ifstream data_stream(dataFile);
   
   if (!data_stream){
      std::cout << "Error opening input data file. Aborting.\n";
      abort();
   }
   
   std::string key;
   std::string file;
   
   while (data_stream >> key){
      if(key == "UPPER"){
	 data_stream >> file;
	 setUpperFile(file);
      }
      else if(key == "LOWER"){
	 data_stream >> file;
	 setLowerFile(file);
      }
      else if(key == "AMPL_MODEL"){
	 data_stream >> file;
	 setUpperAmplModelFile(file);
      }
      else if(key == "AMPL_DATA"){
	 data_stream >> file;
	 setUpperAmplDataFile(file);
      }
   }
   data_stream.close();
#endif

  setUpperFile(dataFile);
  //readAuxiliaryData(); // reads in lower-level vars, rows, obj coeffs
  readProblemData(); // reads in max c^1x + d^1y s.t. (x,y) in Omega^I
}

//#############################################################################
void 
MibSModel::setBlisParameters()
{

  int bliscuts(MibSPar_->entry(MibSParams::blisCutStrategy));
  int blisbranch(MibSPar_->entry(MibSParams::blisBranchStrategy));

  /* Set Blis Parameters to keep cutting until no cut is found */
  BlisPar()->setEntry(BlisParams::cutFactor, ALPS_DBL_MAX);
  BlisPar()->setEntry(BlisParams::cutPass, ALPS_INT_MAX);
  BlisPar()->setEntry(BlisParams::tailOff, -10000);
  BlisPar()->setEntry(BlisParams::denseConFactor, ALPS_DBL_MAX);
 
  /* Set cut generation frequency to 1 */
  BlisPar()->setEntry(BlisParams::cutGenerationFrequency, 1);

  /* Set Blis cut strategy using MibS parameters blisCutStrategy */
  BlisPar()->setEntry(BlisParams::cutStrategy, bliscuts);
  /* Set Blis branch strategy using MibS parameters blisBranchStrategy */
  BlisPar()->setEntry(BlisParams::cutStrategy, blisbranch);
  
}

//#############################################################################
void 
MibSModel::readAuxiliaryData(int numCols, int numRows)
{

  std::string fileName = getLowerFile();
  fileCoinReadable(fileName);
  std::ifstream data_stream(fileName.c_str());

  std::string inputFormat(MibSPar_->entry
			  (MibSParams::inputFormat));
  
  if (!data_stream){
      std::cout << "Error opening input data file. Aborting.\n";
      abort();
  }
  
  std::string key;
  int iValue(0);
  std::string cValue;
  double dValue(0.0);
  int i(0), j(0), k(0), m(0), p(0), pos(0);
  int lowerColNum(0), lowerRowNum(0);

  while (data_stream >> key){
     if(key == "N"){ 
	   data_stream >> iValue;
	   setLowerDim(iValue);
     }
     else if(key == "M"){
	   data_stream >> iValue;
	   setLowerRowNum(iValue);
     }
     else if(key == "LC"){
       if(!getLowerColInd())
	 lowerColInd_ = new int[getLowerDim()];

       if(inputFormat == "indexBased"){
	data_stream >> iValue;
	lowerColInd_[i] = iValue;
       }
       else{
	   data_stream >> cValue;
	   for(p = 0; p < numCols; ++p){
	       if(columnName_[p] == cValue){
		   pos = p;
		   break;
	       }
	   }
	   lowerColInd_[i] = pos;
       }
       
	i++;
     }
     else if(key == "LR"){
       if(!getLowerRowInd())
	 lowerRowInd_ = new int[getLowerRowNum()];
       
       if(inputFormat == "indexBased"){
	   data_stream >> iValue;
	   lowerRowInd_[j] = iValue;
       }
       else{
	   data_stream >> cValue;
	   for(p = 0; p < numRows; ++p){
	       if(rowName_[p] == cValue){
		   pos = p;
		   break;
	       }
	   }
	   lowerRowInd_[j] = pos;
       }
	
	j++;
     }
     else if(key == "LO"){
       if(!getLowerObjCoeffs())
	 lowerObjCoeffs_ = new double[getLowerDim()];
	
	data_stream >> dValue;
	lowerObjCoeffs_[k] = dValue;
	k++;
     }
     else if(key == "OS"){
	data_stream >> dValue;
	lowerObjSense_ = dValue; //1 min; -1 max
     }
     else if(key == "IC"){
       if(!getInterdictCost()){
	 //FIXME: ALLOW MORE THAN ONE ROW
	 interdictCost_ = new double[getLowerDim()];
       }
       
       data_stream >> dValue;
       interdictCost_[m] = dValue;
       m++;
     }
     else if(key == "IB"){
	 isInterdict_ = true;
       //FIXME: ALLOW MORE THAN ONE ROW
	data_stream >> dValue;
	interdictBudget_ = dValue;
     }
     else if(key == "@VARSBEGIN"){
	 pos = -1;
	 lowerColNum = getLowerDim();
	 if(!getLowerColInd())
	     lowerColInd_ = new int[lowerColNum];
	 if(!getLowerObjCoeffs())
	     lowerObjCoeffs_ = new double[lowerColNum];
	 for(i = 0; i < lowerColNum; i++){
	     data_stream >> cValue;
	     data_stream >> dValue;
	     for(p = 0; p < numCols; ++p){
		 if(columnName_[p] == cValue){
		     pos = p;
		     break;
		 }
	     }
	     if(pos < 0){
		 std::cout << cValue << " does not belong to the list of variables." << std::endl;
	     }
	     else{
		 lowerObjCoeffs_[pos] = dValue;
		 lowerColInd_[i] = pos;
	     }
	     pos = -1;
	 }
     }
     else if(key == "@CONSTSBEGIN"){
	 pos = -1;
	 lowerRowNum = getLowerRowNum();
	 if(!getLowerRowInd())
	     lowerRowInd_ = new int[lowerRowNum];
	 for(i = 0; i < lowerRowNum; i++){
	     data_stream >> cValue;
	     for(p = 0; p < numRows; ++p){
		 if(rowName_[p] == cValue){
		     pos = p;
		     break;
		 }
	     }
	     if(pos < 0){
		 std::cout << cValue << " does not belong to the list of constraints." << std::endl;
	     }
	     else{
		 lowerRowInd_[i] = pos;
	     }
	     pos = -1;
	 }
     }
  }

  data_stream.close();
  
  std::cout << "LL Data File: " << getLowerFile() << "\n";
  std::cout << "Number of LL Variables:   " 
	    << getLowerDim() << "\n\n";
}



//#############################################################################
void 
MibSModel::loadAuxiliaryData(int lowerColNum, int lowerRowNum,
			     const int *lowerColInd,
			     const int *lowerRowInd,
			     double lowerObjSense,
			     const double *lowerObjCoef,
			     int upperColNum, int upperRowNum,
			     const int *upperColInd,
			     const int *upperRowInd,
			     int structRowNum, 
			     const int *structRowInd,
			     double interdictBudget, 
			     const double *interdictCost,
			     const double *lColLbInLProb,
			     const double *lColUbInLProb)
{
   int *copyLowerColInd = new int[lowerColNum];
   int *copyLowerRowInd = new int[lowerRowNum];
   double *copyLowerObjCoef = new double[lowerColNum];
   int *copyUpperColInd = NULL;
   int *copyUpperRowInd = NULL;
   int *copyStructRowInd = NULL;
   double *copyInterdictCost = NULL;
   double *copyLColLbInLProb = NULL;
   double *copyLColUbInLProb = NULL;

   if (upperColInd != NULL){
      copyUpperColInd = new int[upperColNum];
   }
   if (upperRowInd != NULL){
      copyUpperRowInd = new int[upperRowNum];
   }
   if (structRowInd != NULL){
      copyStructRowInd = new int[structRowNum];
   }
   if (interdictCost != NULL){
      copyInterdictCost = new double[lowerColNum];
   }

   if (lColLbInLProb != NULL){
       copyLColLbInLProb = new double[lowerColNum];
   }

   if (lColUbInLProb != NULL){
       copyLColUbInLProb = new double[lowerColNum];
   }

   CoinDisjointCopyN(lowerColInd, lowerColNum, copyLowerColInd);
   CoinDisjointCopyN(lowerRowInd, lowerRowNum, copyLowerRowInd);
   CoinDisjointCopyN(lowerObjCoef, lowerColNum, copyLowerObjCoef);
   if (upperColInd != NULL){
      CoinDisjointCopyN(upperColInd, upperColNum, copyUpperColInd);
   }
   if (upperRowInd != NULL){
      CoinDisjointCopyN(upperRowInd, upperRowNum, copyUpperRowInd);
   }
   if (structRowInd != NULL){
      CoinDisjointCopyN(structRowInd, structRowNum, copyStructRowInd);
   }
   if (interdictCost != NULL){
      CoinDisjointCopyN(interdictCost, lowerColNum, copyInterdictCost);
   }
   if (lColLbInLProb != NULL){
       CoinDisjointCopyN(lColLbInLProb, lowerColNum, copyLColLbInLProb);
   }
   if (lColUbInLProb != NULL){
       CoinDisjointCopyN(lColUbInLProb, lowerColNum, copyLColUbInLProb);
   }

   setLowerDim(lowerColNum);
   setLowerRowNum(lowerRowNum);
   setLowerColInd(copyLowerColInd);
   setLowerRowInd(copyLowerRowInd);
   setLowerObjSense(lowerObjSense);
   setLowerObjCoeffs(copyLowerObjCoef);
   if (upperColInd != NULL){
      setUpperDim(upperColNum);
      setUpperColInd(copyUpperColInd);
   }
   if (upperRowInd != NULL){
      setUpperRowNum(upperRowNum);
      setUpperRowInd(copyUpperRowInd);
   }
   if (structRowInd != NULL){
      setStructRowNum(structRowNum);
      setStructRowInd(copyStructRowInd);
   }
   if (interdictCost != NULL){
       isInterdict_ = true;
      setInterdictBudget(interdictBudget);
      setInterdictCost(copyInterdictCost);
   }
   if (lColLbInLProb != NULL){
       setLColLbInLProb(copyLColLbInLProb);
   }
   if (lColUbInLProb != NULL){
       setLColUbInLProb(copyLColUbInLProb);
   }
   
}

//#############################################################################
void 
MibSModel::readProblemData()
{

   int j(0);
   
   int msgLevel(AlpsPar_->entry(AlpsParams::msgLevel));
   
   //------------------------------------------------------
   // Read in data from MPS or AMPL/GMPL file.
   // AMPL/GMPL files often have a separate data file
   // sepUpperData indicates if this is the case
   //------------------------------------------------------
   
   int rc(-1);
   int format(MibSPar_->entry(MibSParams::upperFileFormat));

   CoinMpsIO *mps = new CoinMpsIO;

   switch(format){
     
   case 0: // mps
     {
       rc = mps->readMps(getUpperFile().c_str(), "");
       break;
     }
   case 1: // ampl/gmpl
     {
       rc = mps->readGMPL(getUpperAmplModelFile().c_str(), 
			  getUpperAmplDataFile().c_str());
       break;
     }
   }
   
   //int rc = mps->readMps(getUpperFile().c_str(), "");
   
   if(rc) {
      delete mps;
      throw CoinError("Unable to read in instance",
		      "readInstance",
		      "MibSModel");
   }

   
   mps->messageHandler()->setLogLevel(msgLevel);
   
   CoinPackedMatrix matrix = *(mps->getMatrixByCol());

   double objSense(1.0);
   
   char * colType = NULL;

   int numCols = mps->getNumCols(); 
   int numRows = mps->getNumRows();

   std::string tmpString;
   
   if(!getColumnName())
       columnName_ = new std::string [numCols];
       
   for(j = 0; j < numCols; ++j){
       std::string tmpString(mps->columnName(j));
       columnName_[j] = tmpString;
   }

   if(!getRowName())
       rowName_ = new std::string [numRows];

   for(j = 0; j < numRows; ++j){
       std::string tmpString(mps->rowName(j));
       rowName_[j] = tmpString;
   }
   
   double *varLB = new double [numCols];
   double *varUB = new double [numCols];
   
   double *conLB = new double [numRows];
   double *conUB = new double [numRows];
   
   memcpy(varLB, mps->getColLower(), sizeof(double) * numCols);
   memcpy(varUB, mps->getColUpper(), sizeof(double) * numCols);
   
   memcpy(conLB, mps->getRowLower(), sizeof(double) * numRows);
   memcpy(conUB, mps->getRowUpper(), sizeof(double) * numRows);
   
   //------------------------------------------------------
   // Set colType_
   //------------------------------------------------------
   
   colType = new char [numCols];   
   
   for(j = 0; j < numCols; ++j) {
      if (mps->isContinuous(j)) {
	 colType[j] = 'C';
      }
      else {
	 if (varLB[j] == 0 && varUB[j] == 1.0) {
	    colType[j] = 'B';
	 }
	 else {
	    colType[j] = 'I';
	 }
      }
   }
   
   CoinPackedMatrix colMatrix = *(mps->getMatrixByCol());

   //FIXME: MPS is always minimization, but should we be able to override?
   //FIXME: In previous version of code, objSense was only set to -1
   //       for interdiction problems...
   //objSense = BlisPar_->entry(BlisParams::objSense);

   double *objCoef = new double [numCols];
   
   const double *mpsObj =  mps->getObjCoefficients();

   memcpy(objCoef, mpsObj, sizeof(double) * numCols);

   const char* rowSense = mps->getRowSense();

   readAuxiliaryData(numCols, numRows); // reads in lower-level vars, rows, obj coeffs 
   
   loadProblemData(matrix, varLB, varUB, objCoef, conLB, conUB, colType, 
		   objSense, mps->getInfinity(), rowSense);

   delete [] colType;
   delete [] varLB;
   delete [] varUB;
   delete [] conLB;
   delete [] conUB;
   delete [] objCoef;

   delete mps;
}

//#############################################################################
void
MibSModel::loadProblemData(const CoinPackedMatrix& matrix,
			   const double* colLB, const double* colUB,   
			   const double* obj, const double* rowLB,
			   const double* rowUB, const char *types,
			   double objSense, double infinity, const char *rowSense)
{
   //FIXME: THIS ISN'T TRUE IF WE LOAD AN INTERDICTION PROBLEM 
   //AS A "GENERAL" PROBLEM.  FOR NOW, IT'S OK SINCE WE ONLY
   //DO THIS FROM KNAP SOLVER, WHICH SHOULD SET THIS ITSELF.

   int problemType(MibSPar_->entry(MibSParams::bilevelProblemType));

   int i(0);
   
   if(isInterdict_ == true){
       if(problemType == PARAM_NOTSET){
	   MibSPar()->setEntry(MibSParams::bilevelProblemType, INTERDICT);
       }
       else if(problemType == GENERAL){
	   std::cout << "Wrong value for MibSProblemType. Automatically modifying its value." << std::endl;
	   MibSPar()->setEntry(MibSParams::bilevelProblemType, INTERDICT);
       }
   }
   else{
       if(problemType == PARAM_NOTSET){
	   MibSPar()->setEntry(MibSParams::bilevelProblemType, GENERAL);
       }
   else if(problemType == INTERDICT){
   	   std::cout << "Wrong value for MibSProblemType. Automatically modifying its value." << std::endl;
   	   MibSPar()->setEntry(MibSParams::bilevelProblemType, GENERAL);
     }
   }

   problemType = MibSPar_->entry(MibSParams::bilevelProblemType);
  
   int j(0);
   int beg(0);

   int numRows = matrix.getNumRows();
   int numCols = matrix.getNumCols();

   double *varLB(NULL);  
   double *varUB(NULL);  
   double *conLB(NULL);  
   double *conUB(NULL);  
   double *objCoef(NULL);
   char   *colType(NULL);

   CoinPackedMatrix *newMatrix = NULL;

   //switch (problemType){
      
       //case 0:

   if((problemType == 0) || (interdictCost_ == NULL)){

      if (!structRowInd_){
	 structRowInd_ = new int[numRows];
	 CoinIotaN(structRowInd_, numRows, 0);
	 structRowNum_ = numRows;
      }
      
      //Make copies of the data
      newMatrix = new CoinPackedMatrix();
      *newMatrix = matrix;
      
      varLB = new double [numCols];
      varUB = new double [numCols];
      conLB = new double [numRows];
      conUB = new double [numRows];
      objCoef = new double [numCols];
      colType = new char [numCols];
 
      CoinDisjointCopyN(colLB, numCols, varLB);
      CoinDisjointCopyN(colUB, numCols, varUB);
      CoinDisjointCopyN(rowLB, numRows, conLB);
      CoinDisjointCopyN(rowUB, numRows, conUB);
      CoinDisjointCopyN(obj, numCols, objCoef);
      memcpy(colType, types, numCols);

      origRowSense_ = new char[numRows];
      memcpy(origRowSense_, rowSense, numRows);

   }
   else{    
      //------------------------------------------------------
      // Add interdict vars and aux UL rows
      //------------------------------------------------------
      
      double * intCosts = getInterdictCost();
      int numInterdictNZ(0);
      //FIXME: ALLOW MORE THAN ONE ROW
      int auxULRows(1);
      int interdictRows(numCols);
      int auxRows(auxULRows + interdictRows);
      int numTotalCols(0), numTotalRows(0);
      //int numAuxCols(2 * numCols);
      int numAuxCols(1);
      //int numAuxCols(0);//this breaks orig interdiction cut
      int i(0);
      
      //FIXME:  NEED TO CHANGE THIS AROUND
      //maxAuxCols_ = numAuxCols;
      
      for(i = 0; i < numCols; i++){
	 if((intCosts[i] > etol_) || (intCosts[i] < - etol_)){
	    numInterdictNZ++;
	 }
      }
      
      numTotalCols = 2 * numCols + numAuxCols;
      numTotalRows = numRows + auxRows;
      
      int structRows(numTotalRows - interdictRows);
      structRowInd_ = new int[structRows];
      CoinIotaN(structRowInd_, structRows, 0);
      structRowNum_ = structRows;
      
      varLB = new double [numTotalCols];
      varUB = new double [numTotalCols];
      
      conLB = new double [numTotalRows];
      conUB = new double [numTotalRows];
      
      CoinDisjointCopyN(colLB, numCols, varLB + numCols);
      CoinDisjointCopyN(colUB, numCols, varUB + numCols);

      CoinFillN(varLB, numCols, 0.0); 
      CoinFillN(varUB, numCols, 1.0); 
      
      CoinFillN(varLB + 2 * numCols, numAuxCols, 0.0); 
      CoinFillN(varUB + 2 * numCols, numAuxCols, 1.0); 
      
      CoinDisjointCopyN(rowLB, numRows, conLB + auxULRows);
      CoinDisjointCopyN(rowUB, numRows, conUB + auxULRows);
      
      /* Add interdiction budget row */
      CoinFillN(conLB, auxULRows, - 1 * infinity);
      CoinFillN(conUB, auxULRows, getInterdictBudget());
      
      /* Add VUB rows */
      CoinFillN(conLB + (numTotalRows - interdictRows), 
		interdictRows, - 1 * infinity);
      CoinDisjointCopyN(colUB, interdictRows, conUB + (numTotalRows - interdictRows));
      
      objCoef = new double [numTotalCols];
      CoinZeroN(objCoef, numTotalCols);
      //This is a work-around because the MPS files in our test set have the lower-level
      //objective instead of the upper level one
      for (j = 0; j < numCols; j++){ 
	 objCoef[j + numCols] = -obj[j];
      }
      
      //------------------------------------------------------
      // Set colType_
      //------------------------------------------------------
      
      colType = new char [numTotalCols];   
      
      for(j = 0; j < numCols; ++j) {
	 colType[j] = 'B';
      }
      
      CoinDisjointCopyN(types, numCols, colType + numCols);
      
      /* Auxilliary indicator columns, used later*/
      
      for(j = 0; j < numAuxCols; ++j) {
	 colType[j + 2 * numCols] = 'B';
      }
       
      CoinPackedMatrix rowMatrix;
      rowMatrix = matrix;
      rowMatrix.reverseOrdering();
      const double * matElements = rowMatrix.getElements();
      const int * matIndices = rowMatrix.getIndices();
      const int * matStarts = rowMatrix.getVectorStarts();
      
      newMatrix = new CoinPackedMatrix(false, 0, 0);
      newMatrix->setDimensions(0, numTotalCols);
      int start(0), end(0), tmp(0), index(0);
      
      /* Add interdiction budget row */
      
      for(i = 0; i < auxULRows; i++){
	 CoinPackedVector row;
	 for(j = 0; j < numCols; j++){
	    row.insert(j, intCosts[j]);
	 }
	 newMatrix->appendRow(row);
      }
      
      /* lower-level rows */
      
      for (i = 0; i < numRows; i++){
	 CoinPackedVector row;
	 start = matStarts[i];
	 end = start + rowMatrix.getVectorSize(i);
	 for(j = start; j < end; j++){
	    index = matIndices[j] + numCols;
	    row.insert(index, matElements[j]);
	 }
	 newMatrix->appendRow(row);
      }
      
      /* Add VUB rows */
      
      for (i = 0; i < numCols; i++){
	 CoinPackedVector row;
	 row.insert(i, colUB[i]);
	 row.insert(i + numCols, 1.0);
	 newMatrix->appendRow(row);
      }
      
      newMatrix->reverseOrdering();
      
      setUpperDim(numCols);
      setUpperRowNum(1);

      int *upperColInd = new int[numCols];
      int *upperRowInd = new int[1];      
      CoinIotaN(upperColInd, numCols, 0);
      upperRowInd[0] = 0;

      origRowSense_ = new char[numTotalRows];
      CoinDisjointCopyN(rowSense, numRows, origRowSense_ + auxULRows);
      CoinFillN(origRowSense_, auxULRows, 'L');
      CoinFillN(origRowSense_ + (numTotalRows - interdictRows),
		interdictRows, 'L');

      setUpperColInd(upperColInd);
      setUpperRowInd(upperRowInd);
      
      // store the indices of the structural constraints
      //for(i = 0; i < interdictRows; i++)
      //	vubRowInd_[i] = lowerRowInd_[numTotalRows - interdictRows + i];
      
      //break;
      
   }

   setColMatrix(newMatrix);   
   
   numCoreConstraints_ = numCons_ = newMatrix->getMinorDim();
   numCoreVariables_   = numVars_ = newMatrix->getMajorDim();
   setNumCons(numCons_);
   setNumVars(numVars_);
   
   setConLb(conLB);
   setConUb(conUB);
   
   setVarLb(varLB);
   setVarUb(varUB);
   
   setObjCoef(objCoef);
   
   setColType(colType);

   BlisPar_->setEntry(BlisParams::objSense, objSense);

   //------------------------------------------------------
   // Create variables and constraints.
   //------------------------------------------------------

   const double *elements = newMatrix->getElements();
   const int *indices = newMatrix->getIndices();
   const int *lengths = newMatrix->getVectorLengths();
   const CoinBigIndex *starts = newMatrix->getVectorStarts();
      
   for (j = 0; j < numVars_; ++j) {
      
      beg = starts[j];
      
      BlisVariable * var = new BlisVariable(varLB[j],
					    varUB[j], 
					    varLB[j], 
					    varUB[j],
					    objCoef[j], 
					    lengths[j],
					    indices + beg,
					    elements + beg);
      
      var->setObjectIndex(j);
      var->setRepType(BCPS_CORE);
      var->setStatus(BCPS_NONREMOVALBE);
      var->setIntType(colType_[j]);
      variables_.push_back(var);
      var = NULL;
   }
   
   for (j = 0; j < numCons_; ++j) {
      BlisConstraint *con = new BlisConstraint(conLB[j], 
					       conUB[j], 
					       conLB[j], 
					       conUB[j]);
      con->setObjectIndex(j);
      con->setRepType(BCPS_CORE);
      con->setStatus(BCPS_NONREMOVALBE);
      constraints_.push_back(con);
      con = NULL;        
   }

   numOrigVars_ = numVars_;
   numOrigCons_ = numCons_;

   if(countIteration_ == 0){
       //origConstCoefMatrix_ = new CoinPackedMatrix();
   origConstCoefMatrix_ = newMatrix;
   }

   setUpperColData();
   setUpperRowData();
   setBounds(); // stores the original column and row bounds
   //checkProblemType(); // checks if MibS can solve problem entered
   setProblemType(); //determine the type of MIBLP
   //determine the list of first-stage variables participate in second-stage constraints
   setRequiredFixedList(newMatrix);
   instanceStructure(newMatrix, conLB, conUB, rowSense);
}

//#############################################################################
void 
MibSModel::setUpperColData()
{

   int lowerColNum(lowerDim_);
   if (!upperDim_){
      upperDim_ = numVars_ - lowerDim_;
      int * lowerColInd = getLowerColInd();
      
      if(!getUpperColInd())
	 upperColInd_ = new int[upperDim_];
      
      int i(0), cnt(0);

      for (i = 0; i < lowerDim_ + upperDim_; i++){
	 if((!findIndex(i, lowerColNum, lowerColInd)) 
	    && (colType_[i] != 'C')){
	    upperColInd_[cnt] = i;
	    cnt++;
	    if(0)
	       std::cout << "i: " << i << std::endl;
	 }
      }

      for(i = 0; i < lowerDim_ + upperDim_; i++){
	 if((!findIndex(i, lowerColNum, lowerColInd)) 
	    && (colType_[i] == 'C')){
	    upperColInd_[cnt] = i;
	    cnt++;
	 }
      }
   
      assert(cnt == upperDim_);
   }
   //numOrigVars_ = lowerDim_ + upperDim_;
}

//#############################################################################
void 
MibSModel::setUpperRowData()
{

   //FIXME: MAKE THIS MORE SIMPLE

   int lowerRowNum(lowerRowNum_);
   upperRowNum_ = numCons_ - lowerRowNum_;
   int * lowerRowInd = getLowerRowInd();

   if(upperRowInd_ != NULL){
       delete [] upperRowInd_;
       upperRowInd_ = NULL;
   }

   if(upperRowNum_ > 0){
       upperRowInd_ = new int[upperRowNum_];
   }

   int i(0), cnt(0);

   for(i = 0; i < lowerRowNum_ + upperRowNum_; i++){
       if(!findIndex(i, lowerRowNum, lowerRowInd)){
	   upperRowInd_[cnt] = i;
	   cnt++;
       }
   }
   assert(cnt == upperRowNum_);

   if(countIteration_ == 0){
       origUpperRowNum_ = upperRowNum_;

       if(origUpperRowNum_ > 0){
	   origUpperRowInd_ = new int[origUpperRowNum_];
       }
       
       for(i = 0; i < origUpperRowNum_; i++){
	   origUpperRowInd_[i] = upperRowInd_[i];
       }
   }
   
   //numOrigCons_ = lowerRowNum_ + upperRowNum_;
}

//#############################################################################
void 
MibSModel::findIntegers()
{

}

//#############################################################################
AlpsTreeNode * 
MibSModel::createRoot() {
    
  //-------------------------------------------------------------
  // NOTE: Root will be deleted by ALPS. Root is an explicit node.
  //-------------------------------------------------------------
  
  //determine value function slopes for valfunc cut, if necessary
  bool useValFuncCut 
    = MibSPar_->entry(MibSParams::useValFuncCut);
  
  if(useValFuncCut)
    setValFuncSlopes();
  
  MibSTreeNode * root = new MibSTreeNode;
  
  BlisNodeDesc* desc = new BlisNodeDesc(this);
  root->setDesc(desc);
  
  //-------------------------------------------------------------
  // NOTE: Although original data are stored in model when reading. 
  //   Root desc still store a full copy of col and row bounds when creating.
  //   It will store soft differences after finding a branching object. 
  //   The soft difference are due to reduced cost fixing and probing.
  //   Also the added cols and rows will be stored.
  //-------------------------------------------------------------
  int k;
    
  std::vector<BcpsVariable *> vars = getVariables();
  std::vector<BcpsConstraint *> cons = getConstraints();
  
  int numVars = static_cast<int> (vars.size());
  int numCons = static_cast<int> (cons.size());
  
#ifdef BLIS_DEBUG
  std::cout << "BLIS: createRoot(): numVars=" << numVars
	    << ", numCons=" << numCons 
	    << "; numCoreVariables_=" << numCoreVariables_
	    << ", numCoreConstraints_=" << numCoreConstraints_ << std::endl;
#endif
  
  int *varIndices1 = new int [numVars];
  int *varIndices2 = new int [numVars];
  int *varIndices3 = NULL;
  int *varIndices4 = NULL;
  double *vlhe = new double [numVars];
  double *vuhe = new double [numVars];
  double *vlse = NULL;
  double *vuse = NULL;
  
  int *conIndices1 = new int [numCons];
  int *conIndices2 = new int [numCons];
  int *conIndices3 = NULL;
  int *conIndices4 = NULL;
  double *clhe = new double [numCons];
  double *cuhe = new double [numCons];
  double *clse = NULL;
  double *cuse = NULL;
  
  //-------------------------------------------------------------
  // Get var bounds and indices.
  //-------------------------------------------------------------
  
  for (k = 0; k < numVars; ++k) {
    vlhe[k] = vars[k]->getLbHard();
    vuhe[k] = vars[k]->getUbHard();
    varIndices1[k] = k;
    varIndices2[k] = k;
    
#ifdef BLIS_DEBUG_MORE
    std::cout << "BLIS: createRoot(): var "<< k << ": hard: lb=" << vlhe[k]
	      << ", ub=" << vuhe[k] << std::endl;
#endif  
    }
  
  //-------------------------------------------------------------  
  // Get con bounds and indices.
  //-------------------------------------------------------------
  
  for (k = 0; k < numCons; ++k) {
    clhe[k] = cons[k]->getLbHard();
    cuhe[k] = cons[k]->getUbHard();
    conIndices1[k] = k;
    conIndices2[k] = k;
  }
  
  int *tempInd = NULL;
  BcpsObject **tempObj = NULL;
  
  desc->assignVars(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numVars, varIndices1, vlhe, /*Var hard lb*/
		   false, numVars, varIndices2, vuhe, /*Var hard ub*/
		   false, 0, varIndices3, vlse,       /*Var soft lb*/
		   false, 0, varIndices4, vuse);      /*Var soft ub*/
  desc->assignCons(0 /*numRem*/, tempInd,
		   0 /*numAdd*/, tempObj,
		   false, numCons, conIndices1, clhe, /*Con hard lb*/
		   false, numCons, conIndices2, cuhe, /*Con hard ub*/
		   false, 0,conIndices3,clse,         /*Con soft lb*/
		   false, 0,conIndices4,cuse);        /*Con soft ub*/
    
  //-------------------------------------------------------------  
  // Mark it as an explicit node.
  //-------------------------------------------------------------

  root->setExplicit(1);
  
  return root;
}

//############################################################################ 

/** Do necessary work to make model usable and presolve.
    Return success or not. 
    This function is called when constructing knowledge broker. */
bool 
MibSModel::setupSelf()
{
   
   int j;
   
   bcpsMessageHandler_->setLogLevel(broker_->getMsgLevel());
   blisMessageHandler_->setLogLevel(broker_->getMsgLevel());
   
   if (broker_->getMsgLevel() > 0) {
      
      //std::cout << "**** getProcType = " << broker_->getProcType() << std::endl;
      
      if (broker_->getMsgLevel() > 0) {
	 if (broker_->getProcRank() == broker_->getMasterRank()) {
	    if (strcmp(BCPS_VERSION, "trunk")){
	       std::cout << "==  Bcps Version: " << BCPS_VERSION << std::endl;
	    }else{
	       std::cout << "==  Bcps Version: Trunk (unstable)" << std::endl;
	    }
#ifdef BCPS_SVN_REV
	    std::cout << "==  Bcps Revision Number: " << BCPS_SVN_REV
	        << std::endl;
#endif
	    if (strcmp(BLIS_VERSION, "trunk")){
	       std::cout << "==  Blis Version: " << BLIS_VERSION << std::endl;
	    }else{
	       std::cout << "==  Blis Version: Trunk (unstable)" << std::endl;
	    }
#ifdef BLIS_SVN_REV
	    std::cout << "==  Blis Revision Number: " << BLIS_SVN_REV
	        << std::endl;
#endif
	    std::cout << std::endl;
	 }
      }
      
   }
   
   //------------------------------------------------------
   // Set numIntObjects_, intColIndices_, intObjectIndices_ 
   //------------------------------------------------------
   
   // TODO: now integer, later other objects
   intObjIndices_ = new int [numCols_];
   memset(intObjIndices_, 0, sizeof(int) * numCols_);
   
   numIntObjects_ = 0;  
   intColIndices_ = new int [numCols_];
   
   for(j = 0; j < numCols_; ++j) {
	if (colType_[j] == 'I' || colType_[j] == 'B') {
	   intColIndices_[numIntObjects_++] = j;
	}
   }
   if (numIntObjects_ == 0) {
      if (broker_->getMsgLevel() > 0) {
	 bcpsMessageHandler_->message(BLIS_W_LP, blisMessages())
	    << CoinMessageEol;
      }
   }   
   
   //------------------------------------------------------
   // Load data to LP solver.
   //------------------------------------------------------

   if (!lpSolver_) {
      // preprocessing causes this check.
      lpSolver_ = origLpSolver_;
   }
   
   //colMatrix_->dumpMatrix();

   lpSolver_->loadProblem(*colMatrix_,
			  varLB_, varUB_, 
			  objCoef_,
			  conLB_, conUB_);
   
   lpSolver_->setObjSense(BlisPar_->entry(BlisParams::objSense));
   lpSolver_->setInteger(intColIndices_, numIntObjects_);

   //------------------------------------------------------
   // Preprocess model
   //------------------------------------------------------

   int usePreprocessor =
     MibSPar_->entry(MibSParams::usePreprocessor);

   if(usePreprocessor == PARAM_ON){
     runPreprocessor();
   }
   //------------------------------------------------------
   // Create integer objects and analyze objective coef.
   //------------------------------------------------------
   
   createIntgerObjects(true);   

   // Do this after loading LP.
   analyzeObjective();
   
   //------------------------------------------------------
   // Allocate memory.
   //------------------------------------------------------
   
   startVarLB_ = new double [numCols_];
   startVarUB_ = new double [numCols_];
   
   startConLB_ = new double [numRows_];
   startConUB_ = new double [numRows_];
   
   tempVarLBPos_ = new int [numCols_];
   tempVarUBPos_ = new int [numCols_];
    
   tempConLBPos_ = new int [numRows_];
   tempConUBPos_ = new int [numRows_];
   
   //------------------------------------------------------
   // Get parameters.
   //------------------------------------------------------
   
   integerTol_ = BlisPar_->entry(BlisParams::integerTol);
   optimalRelGap_ = BlisPar_->entry(BlisParams::optimalRelGap);
   optimalAbsGap_ = BlisPar_->entry(BlisParams::optimalAbsGap);
   
   cutoff_ =  BlisPar_->entry(BlisParams::cutoff);
   
   //------------------------------------------------------
   // Modify parameters.
   //------------------------------------------------------
   
   // Disable Alps message
   // AlpsPar()->setEntry(AlpsParams::msgLevel, 1);
   
#ifdef BLIS_DEBUG_MORE
   std::string problemName;
   lpSolver_->getStrParam(OsiProbName, problemName);
   printf("BLIS: setupSelf: Problem name - %s\n", problemName.c_str());
   lpSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif

   //------------------------------------------------------
   // Set branch strategy.
   //------------------------------------------------------
   
   int brStrategy = BlisPar_->entry(BlisParams::branchStrategy);

   /*   
   if (brStrategy == BlisBranchingStrategyMaxInfeasibility) {
      // Max inf
      branchStrategy_ = new BlisBranchStrategyMaxInf(this);
   }
   else if (brStrategy == BlisBranchingStrategyPseudoCost) {
      // Pseudocost
      branchStrategy_ = new BlisBranchStrategyPseudo(this, 1);
   }
   else if (brStrategy == BlisBranchingStrategyReliability) {
      // Relibility
      branchStrategy_ = new BlisBranchStrategyRel(this, relibility);
   }
   else if (brStrategy == BlisBranchingStrategyStrong) {
      // Strong
      branchStrategy_ = new BlisBranchStrategyStrong(this);
    }
   else {
      throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
   }
   
   brStrategy = BlisPar_->entry(BlisParams::branchStrategyRampUp);
   
   if (brStrategy == BlisBranchingStrategyMaxInfeasibility) {
      // Max inf
      rampUpBranchStrategy_ = new BlisBranchStrategyMaxInf(this);
   }
   else if (brStrategy == BlisBranchingStrategyPseudoCost) {
      // Pseudocost
      rampUpBranchStrategy_ = new BlisBranchStrategyPseudo(this, 1);
   }
   else if (brStrategy == BlisBranchingStrategyReliability) {
      // Relibility
      rampUpBranchStrategy_ = new BlisBranchStrategyRel(this, relibility);
   }
   else if (brStrategy == BlisBranchingStrategyStrong) {
        // Strong
      rampUpBranchStrategy_ = new BlisBranchStrategyStrong(this);
   }
   else {
      throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
   }
   */

   if (brStrategy == BlisBranchingStrategyMaxInfeasibility) {
      // Max inf
      branchStrategy_ = new MibSBranchStrategyMaxInf(this);
   }
   else if (brStrategy == BlisBranchingStrategyPseudoCost) {
      // Pseudocost
      branchStrategy_ = new MibSBranchStrategyPseudo(this, 1);
   }
   else if (brStrategy == BlisBranchingStrategyReliability) {
      // Relibility
      std::cout << "Reliability branching doesn't work. Sorry." << std::endl;
      abort();
      //branchStrategy_ = new MibSBranchStrategyRel(this, relibility);
   }
   else if (brStrategy == BlisBranchingStrategyStrong) {
      // Strong
      branchStrategy_ = new MibSBranchStrategyStrong(this);
    }
   else {
      throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
   }
   
   brStrategy = BlisPar_->entry(BlisParams::branchStrategyRampUp);
   
   if (brStrategy == BlisBranchingStrategyMaxInfeasibility) {
      // Max inf
      rampUpBranchStrategy_ = new MibSBranchStrategyMaxInf(this);
   }
   else if (brStrategy == BlisBranchingStrategyPseudoCost) {
      // Pseudocost
      rampUpBranchStrategy_ = new MibSBranchStrategyPseudo(this, 1);
   }
   else if (brStrategy == BlisBranchingStrategyReliability) {
      // Relibility
      std::cout 
	 << "Reliability branching doesn't work. Sorry." << std::endl;
      abort();
      //rampUpBranchStrategy_ = new MibSBranchStrategyRel(this, relibility);
   }
   else if (brStrategy == BlisBranchingStrategyStrong) {
        // Strong
      rampUpBranchStrategy_ = new MibSBranchStrategyStrong(this);
   }
   else {
      throw CoinError("Unknown branch strategy.", "setupSelf","MibSModel");
   }


   //------------------------------------------------------
   // Add heuristics.
   //------------------------------------------------------
   
   heurStrategy_ = 
      static_cast<BlisHeurStrategy> (BlisPar_->entry(BlisParams::heurStrategy));
   heurCallFrequency_ = 
      static_cast<BlisHeurStrategy> (BlisPar_->entry(BlisParams::heurCallFrequency));
   BlisHeurStrategy heurRoundStrategy = 
      static_cast<BlisHeurStrategy> (BlisPar_->entry(BlisParams::heurRoundStrategy)); 
   int callFreq;
   
   if (heurRoundStrategy == BlisHeurStrategyNotSet) {
      heurRoundStrategy = heurStrategy_;
      callFreq = heurCallFrequency_;
   }
   if (heurRoundStrategy != BlisHeurStrategyNone && 
	heurRoundStrategy != BlisHeurStrategyNotSet) {
      // Add rounding heuristic
      BlisHeurRound *heurRound = new BlisHeurRound(this, 
						   "Rounding",
						   heurRoundStrategy,
						   callFreq);
      addHeuristic(heurRound);
   }
   
    // Adjust heurStrategy
   for (j = 0; j < numHeuristics_; ++j) {
      if (heuristics_[j]->strategy() != BlisHeurStrategyNone) {
	 // Doesn't matter what's the strategy, we just want to 
	 // call heuristics.
	 heurStrategy_ = BlisHeurStrategyAuto;
	 break;
      }
   }
   
   //------------------------------------------------------
   // Cut generators settings.
   //------------------------------------------------------
   
   // Compute dense cutoff.
    
   const CoinPackedMatrix * rowMatrix = lpSolver_->getMatrixByRow();
   const int * rowLen = rowMatrix->getVectorLengths();
   double maxLen = 0.0, minLen = ALPS_DBL_MAX, sumLen = 0.0;
   double aveLen, diffLen, stdLen;
   double denseConFactor = BlisPar_->entry(BlisParams::denseConFactor);
   
   for (j = 0; j < numRows_; ++j) {
      if (rowLen[j] > maxLen) maxLen = rowLen[j];
      if (rowLen[j] < minLen) minLen = rowLen[j];
      sumLen += rowLen[j];
   }
   assert(numRows_ > 0);
   aveLen = sumLen / numRows_;
   sumLen = 0.0;
   
   for (j = 0; j < numRows_; ++j) {
      diffLen = rowLen[j] - aveLen;
      diffLen *= diffLen;
      sumLen += diffLen;
   }
   stdLen = sqrt(sumLen/numRows_);
   if (denseConFactor > 10e5) {
      denseConCutoff_ = ALPS_INT_MAX;
   }
   else {
      denseConCutoff_ = static_cast<int>(aveLen + denseConFactor*stdLen);    
      denseConCutoff_ = ALPS_MIN(numCols_/2, denseConCutoff_);
      denseConCutoff_ = ALPS_MAX(100, denseConCutoff_);
   }
    
#ifdef BLIS_DEBUG
   std::cout << "aveLen=" << aveLen << ", minLen=" << minLen
	     << ", maxLen=" << maxLen << ", stdLen=" << stdLen
	     << ", denseConCutoff_=" << denseConCutoff_ << std::endl;
#endif
   
   // NOTE: maxNumCons is valid only for automatic strategy.
   double cutFactor = BlisPar_->entry(BlisParams::cutFactor);
   
   if (cutFactor > 1.0e5) {
      maxNumCons_ = ALPS_INT_MAX;
      // old constraint size will automatically double if no space.
      oldConstraintsSize_ = 10000;
   }
   else {
      maxNumCons_ = (int)((cutFactor - 1.0) * numRows_);
      oldConstraintsSize_ = maxNumCons_;
   }
   
   constraintPool_ = new BcpsConstraintPool();
   constraintPoolReceive_ = new BcpsConstraintPool();
   constraintPoolSend_ = new BcpsConstraintPool();
   oldConstraints_ = new BlisConstraint* [oldConstraintsSize_];
   
   cutStrategy_ = static_cast<BlisCutStrategy> 
      (BlisPar_->entry(BlisParams::cutStrategy)); 
   cutGenerationFrequency_ = static_cast<BlisCutStrategy> 
      (BlisPar_->entry(BlisParams::cutGenerationFrequency));
   
   if (cutGenerationFrequency_ < 1) {
      std::cout << "WARNING: Input cut generation frequency is " 
		<< cutGenerationFrequency_ 
		<< ", which is not allowed. Changed it to 1" << std::endl;
      cutGenerationFrequency_ = 1;
   }
#if 0
   std::cout << "Initially, cutStrategy_ = " << cutStrategy_ 
	     << "; freq = " << cutGenerationFrequency_
	     << std::endl;
#endif
    
   BlisCutStrategy cliqueStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutCliqueStrategy));
   BlisCutStrategy fCoverStrategy = static_cast<BlisCutStrategy>
      (BlisPar_->entry(BlisParams::cutFlowCoverStrategy));
   BlisCutStrategy gomoryStrategy = static_cast<BlisCutStrategy> 
      (BlisPar_->entry(BlisParams::cutGomoryStrategy)); 
   BlisCutStrategy knapStrategy = static_cast<BlisCutStrategy>
      (BlisPar_->entry(BlisParams::cutKnapsackStrategy)); 
   BlisCutStrategy mirStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutMirStrategy)); 
   BlisCutStrategy oddHoleStrategy = static_cast<BlisCutStrategy> 
      (BlisPar_->entry(BlisParams::cutOddHoleStrategy));
   BlisCutStrategy probeStrategy = static_cast<BlisCutStrategy> 
      (BlisPar_->entry(BlisParams::cutProbingStrategy));
   BlisCutStrategy twoMirStrategy = static_cast<BlisCutStrategy> 
      (BlisPar_->entry(BlisParams::cutTwoMirStrategy)); 
   
   int cliqueFreq = BlisPar_->entry(BlisParams::cutCliqueFreq);
   int fCoverFreq = BlisPar_->entry(BlisParams::cutFlowCoverFreq);
   int gomoryFreq = BlisPar_->entry(BlisParams::cutGomoryFreq); 
   int knapFreq = BlisPar_->entry(BlisParams::cutKnapsackFreq); 
   int mirFreq = BlisPar_->entry(BlisParams::cutMirFreq); 
   int oddHoleFreq = BlisPar_->entry(BlisParams::cutOddHoleFreq);
   int probeFreq = BlisPar_->entry(BlisParams::cutProbingFreq);
   int twoMirFreq = BlisPar_->entry(BlisParams::cutTwoMirFreq); 
   
   //------------------------------------------------------
   // Add cut generators.
   //------------------------------------------------------
   
   
   //----------------------------------
   // Add probe cut generator.
   //----------------------------------
   
   if (probeStrategy == BlisCutStrategyNotSet) {
      // Disable by default
      if (cutStrategy_ == BlisCutStrategyNotSet) {
	 probeStrategy = BlisCutStrategyNone;
      }
      else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	 probeStrategy = cutStrategy_;
	 probeFreq = cutGenerationFrequency_;
      }
	else {
	   probeStrategy = cutStrategy_;
	}
   }
   if (probeStrategy != BlisCutStrategyNone) {
      CglProbing *probing = new CglProbing;
      probing->setUsingObjective(true);
      probing->setMaxPass(1);
      probing->setMaxPassRoot(5);
      // Number of unsatisfied variables to look at
      probing->setMaxProbe(10);
      probing->setMaxProbeRoot(1000);
      // How far to follow the consequences
      probing->setMaxLook(50);
      probing->setMaxLookRoot(500);
      // Only look at rows with fewer than this number of elements
      probing->setMaxElements(200);
      probing->setRowCuts(3);
      addCutGenerator(probing, "Probing", probeStrategy, probeFreq);
   }
   
   //----------------------------------
   // Add clique cut generator.
   //----------------------------------
   
    if (cliqueStrategy == BlisCutStrategyNotSet) {
       // Only at root by default
       if (cutStrategy_ == BlisCutStrategyNotSet) {
	  cliqueStrategy = BlisCutStrategyRoot;
       }
       else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	  cliqueFreq = cutGenerationFrequency_;
	  cliqueStrategy = BlisCutStrategyPeriodic;
       }
       else { // Root or Auto
	  cliqueStrategy = cutStrategy_;
       }
    }
    if (cliqueStrategy != BlisCutStrategyNone) {
       CglClique *cliqueCut = new CglClique ;
       cliqueCut->setStarCliqueReport(false);
       cliqueCut->setRowCliqueReport(false);
       addCutGenerator(cliqueCut, "Clique", cliqueStrategy, cliqueFreq);
    }

    //----------------------------------
    // Add odd hole cut generator.
    //----------------------------------

    if (oddHoleStrategy == BlisCutStrategyNotSet) {
       if (cutStrategy_ == BlisCutStrategyNotSet) {
	  // Disable by default
	  oddHoleStrategy = BlisCutStrategyNone;
       }
       else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	  oddHoleStrategy = BlisCutStrategyPeriodic;
	  oddHoleFreq = cutGenerationFrequency_;
       }
	else {
	   oddHoleStrategy = cutStrategy_;
	}
    }
    if (oddHoleStrategy != BlisCutStrategyNone) {
       CglOddHole *oldHoleCut = new CglOddHole;
       oldHoleCut->setMinimumViolation(0.005);
       oldHoleCut->setMinimumViolationPer(0.00002);
       // try larger limit
       oldHoleCut->setMaximumEntries(200);
       addCutGenerator(oldHoleCut, "OddHole",oddHoleStrategy,oddHoleFreq);
    }

    //----------------------------------
    // Add flow cover cut generator.
    //----------------------------------

    if (fCoverStrategy == BlisCutStrategyNotSet) {
       if (cutStrategy_ == BlisCutStrategyNotSet) {
	  fCoverStrategy = BlisCutStrategyAuto;
	  fCoverFreq = cutGenerationFrequency_;
       }
       else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	  fCoverStrategy = cutStrategy_;
	  fCoverFreq = cutGenerationFrequency_;
       }
       else {
	  fCoverStrategy = cutStrategy_;
       }
    }
    if (fCoverStrategy != BlisCutStrategyNone) {
       CglFlowCover *flowGen = new CglFlowCover;
       addCutGenerator(flowGen, "Flow Cover", fCoverStrategy, fCoverFreq);
    }
    
    //----------------------------------
    // Add knapsack cut generator.
    //----------------------------------

    if (knapStrategy == BlisCutStrategyNotSet) {
       if (cutStrategy_ == BlisCutStrategyNotSet) {
	  // Only at root by default
	  knapStrategy = BlisCutStrategyRoot;
       }
       else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	  knapStrategy = cutStrategy_;
	  knapFreq = cutGenerationFrequency_;
       }
       else {
	  knapStrategy = cutStrategy_;
       }
    }
    if (knapStrategy != BlisCutStrategyNone) {
       CglKnapsackCover *knapCut = new CglKnapsackCover;
       addCutGenerator(knapCut, "Knapsack", knapStrategy, knapFreq);
    }
    
    //----------------------------------
    // Add MIR cut generator.
    //----------------------------------

    if (mirStrategy == BlisCutStrategyNotSet) {
       if (cutStrategy_ == BlisCutStrategyNotSet) {
	  // Disable by default
	  mirStrategy = BlisCutStrategyNone;
       }
       else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	  mirStrategy = cutStrategy_;
	  mirFreq = cutGenerationFrequency_;
       }
       else {
	  mirStrategy = cutStrategy_;
       }
    }
    if (mirStrategy != BlisCutStrategyNone) {
       CglMixedIntegerRounding2 *mixedGen = new CglMixedIntegerRounding2;
       addCutGenerator(mixedGen, "MIR", mirStrategy, mirFreq);
    }
    
    //----------------------------------
    // Add Gomory cut generator.
    //----------------------------------

    if (gomoryStrategy == BlisCutStrategyNotSet) {
       if (cutStrategy_ == BlisCutStrategyNotSet) {
	  // Only at root by default
	  gomoryStrategy = BlisCutStrategyRoot;	    
       }
       else if (cutStrategy_ == BlisCutStrategyPeriodic) {
	  gomoryStrategy = cutStrategy_;
	  gomoryFreq = cutGenerationFrequency_;
       }
       else {
	  gomoryStrategy = cutStrategy_;
       }
    }
    if (gomoryStrategy != BlisCutStrategyNone) {
       CglGomory *gomoryCut = new CglGomory;
       // try larger limit
       gomoryCut->setLimit(300);
       addCutGenerator(gomoryCut, "Gomory", gomoryStrategy, gomoryFreq);
    }
    
    //----------------------------------
    // Add Tow MIR cut generator.
    //----------------------------------

    // Disable forever, not useful.
    twoMirStrategy = BlisCutStrategyNone;
    if (twoMirStrategy != BlisCutStrategyNone) {
       CglTwomir *twoMirCut =  new CglTwomir;
       addCutGenerator(twoMirCut, "Two MIR", twoMirStrategy, twoMirFreq);
    }
    
    //--------------------------------------------
    // Adjust cutStrategy_ according to the strategies of
    // each cut generators.
    //--------------------------------------------

    if (numCutGenerators_ > 0) {
       BlisCutStrategy strategy0 = cutGenerators(0)->strategy();
       BlisCutStrategy strategy1;
       for (j = 1; j < numCutGenerators_; ++j) {
	  strategy1 = cutGenerators(j)->strategy();
	  if (strategy1 != strategy0) {
	     // A generator has different strategy.
	     break;
	  }
	}
       if (j == numCutGenerators_) {
	  // All cut generators has same strategy.
	  cutStrategy_ = strategy0;
       }
       else {
	  // Assume to generate cons at each node since generators 
	  // has various strategies.
	  cutStrategy_ = BlisCutStrategyPeriodic;
	  cutGenerationFrequency_ = 1;
       }
    }
    
#if 0
    std::cout << "AFTER: cutStrategy_ = " << cutStrategy_ << std::endl;
#endif

    //--------------------------------------------
    // Random vector
    //--------------------------------------------

    // TODO, if generating variable, then need more space.
    conRandoms_ = new double [numCols_];
    double deseed = 12345678.0;
    double DSEED2 = 2147483647.0;
        
    for (j = 0; j < numCols_; ++j) {
       deseed *= 16807.0;
       int jseed = (int) (deseed / DSEED2);
       deseed -= (double) jseed * DSEED2;
       double random = deseed / DSEED2;
       conRandoms_[j] = random;
       //std::cout << "conRandoms_[" << j << "]="
       //      <<conRandoms_[j]<< std::endl;
    }
    
    
    return true;
}

//#############################################################################
bool
MibSModel::findIndex(int index, int size, int * indices)
{

  int i(0);
  bool found(false);

  for(i = 0; i < size; i++){
    if(indices[i] == index)
      found = true;
  }

  return found;


}

//#############################################################################
BlisSolution * 
MibSModel::userFeasibleSolution(const double * solution, bool &userFeasible)
{

  /** User's criteria for a feasible solution.
   *  If solution is feasible then need to
   *     1) set userFeasible to true, and
   *     2) return a non-null MIBS solution.
   *  If solution is infeasible then need
   *     1) set userFeasible to false, and
   *     2) return a null.
   **/

  if(0)
    printCurSol();
  
  CoinPackedVector *sol = getSolution();

  MibSSolution *mibSol = NULL;
  
  if(0)
    solver()->writeLp("userfeasible1");

  int i(0), index(0);
  double upperObj(0.0);
  bool isHeurSolution(true);
  int * upperColInd = getUpperColInd();
  int * lowerColInd = getLowerColInd();
  double * lpSolution = new double[getNumCols()];
  CoinFillN(lpSolution, getNumCols(), 0.0);
  MibSSolType solType;
  

  //bool bilevelbranch = MibSPar_->entry(MibSParams::isBilevelBranchProb);
  bool bilevelbranch = false;

  //sahar:check
  if(bilevelbranch){
    for(i = 0; i < solver()->getNumRows(); i++){
      if(solver()->getRowSense()[i] == 'R'){
	double rhs(solver()->getRightHandSide()[i]);
	double range(solver()->getRowRange()[i]);
	//double lower = - solver()->getInfinity();
	//double upper = solver()->getRightHandSide()[i];
	solver()->setRowType(i, 'L', rhs, range);
      }
    }
  }

  if(0)
    solver()->writeLp("userfeasible");

  solType = createBilevel(sol);

  if(solType != MibSNoSol){
      for(i = 0; i < upperDim_; i++){
	  index = upperColInd[i];
	  lpSolution[index] = bS_->optUpperSolutionOrd_[i];
	  upperObj +=
	      bS_->optUpperSolutionOrd_[i] * solver()->getObjCoefficients()[index];
      }
      for(i = 0; i < lowerDim_; i++){
	  index = lowerColInd[i];
	  lpSolution[index] = bS_->optLowerSolutionOrd_[i];
	  upperObj +=
	      bS_->optLowerSolutionOrd_[i] * solver()->getObjCoefficients()[index];
      }
  }

  userFeasible = false;
  if(solType == MibSRelaxationSol){
      userFeasible = true;
  }

  if(userFeasible == true){
      mibSol = new MibSSolution(getNumCols(),
				lpSolution,
				upperObj,
				this);
  }
  else if(solType == MibSHeurSol){
      if(!bS_->isUBSolved_){
	  isHeurSolution = checkUpperFeasibility(lpSolution);
      }
      if(isHeurSolution == true){
	  mibSol = new MibSSolution(getNumCols(),
				    lpSolution,
				    upperObj,
				    this);
	  storeSolution(BlisSolutionTypeHeuristic, mibSol);
	  mibSol = NULL;
      }
  }

  delete sol;
  delete [] lpSolution;
  return mibSol;
}

//#############################################################################
bool
MibSModel::checkUpperFeasibility(double * solution)
{

  bool upperFeasible(true);
  int * uRowIndices = getUpperRowInd();
  int uRows(getUpperRowNum());
  const double * RowLb = getSolver()->getRowLower();
  const double * RowUb = getSolver()->getRowUpper();
  const CoinPackedMatrix * matrix = getSolver()->getMatrixByRow();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();

  double lhs(0.0);
  int i(0), j(0), index1(0), index2(0), start(0), end(0);

  if(!bS_->isUpperIntegral_){
      upperFeasible = false;
  }
  else{
      for(i = 0; i < uRows; i++){
	  index1 = uRowIndices[i];
	  start = matStarts[index1];
	  end = start + matrix->getVectorSize(index1);
	  for(j = start; j < end; j++){
	      index2 = matIndices[j];
	      lhs += matElements[j] * solution[index2];
	  }
	  if(((RowLb[index1] - lhs) > etol_)
	     || ((lhs - RowUb[index1]) > etol_)){
	      upperFeasible = false;
	  }
	  lhs = 0.0;
      }
  }
  
  return upperFeasible;
}

//#############################################################################
double
MibSModel::getObjectiveBound()
{

  /* just one method for now.  may want more later */

  double bound(0.0);
  int strategy(MibSPar_->entry(MibSParams::objBoundStrategy));
  
  switch(strategy){
  case 0: //get bound by optimizing wrt to LL objective
    bound = lowerObjectiveBound();
    break;
  case 1: //for interdiction problems, use reduced costs for bound
    bound = interdictionBound();
    break;
  }

  return bound;

}

//#############################################################################
double
MibSModel::interdictionBound()
{

  /* solve the linear relaxation of the lower-level of an
     interdiction problem and add variables to the 
     upper-level solution using greedy heuristic based on
     lower-level reduced costs. then find a feasible soln.
  */

  OsiSolverInterface * bSolver = new OsiClpSolverInterface();
  double bound(solver()->getInfinity());

  double objSense(getLowerObjSense());  
  int lCols(getLowerDim());
  int uCols(getUpperDim());
  //int tCols(uCols + lCols);
  int lRows(getLowerRowNum());
  int * lColIndices = getLowerColInd();
  int * uColIndices = getUpperColInd();
  int * lRowIndices = getLowerRowInd();
  double * lObjCoeffs = getLowerObjCoeffs();

  int i(0), j(0), index1(0), index2(0);

  double * objCoeffs = new double[lCols];
  double * varUb = new double[lCols];
  double * varLb = new double[lCols];
  double * rowUb = new double[structRowNum_];
  double * rowLb = new double[structRowNum_];
  const CoinPackedMatrix * origMat = solver()->getMatrixByRow();
  const double * matElements = origMat->getElements();
  const int * matIndices = origMat->getIndices();
  const int * matStarts = origMat->getVectorStarts();

  CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);

  //solver()->writeLp("orig");

  CoinZeroN(objCoeffs, lCols);

  for(i = 0; i < lCols; i++)
    objCoeffs[i] = lObjCoeffs[i] * objSense;

  for(i = 0; i < lCols; i++){
    varLb[i] = solver()->getColLower()[lColIndices[i]];
    varUb[i] = solver()->getColUpper()[lColIndices[i]];
  }

  int pos(0), start(0), end(0);

  for(i = 0; i < structRowNum_; i++){
    pos = binarySearch(0, lRows - 1, i, lRowIndices);
    if(pos > -1){
      rowLb[pos] = solver()->getRowLower()[structRowInd_[i]];
      rowUb[pos] = solver()->getRowUpper()[structRowInd_[i]];
    }
  }

  /*
  for(i = 0; i < structRowNum_; i++){
    index = structRowInd_[i];
    if(binarySearch(0, lRows - 1, index, lRowIndices) >= 0)
      newMat->appendMajorVector(origMat->getVector(i));
  }
  */

  for(i = 0; i < structRowNum_; i++){
    index1 = structRowInd_[i];
    if(binarySearch(0, lRows - 1, index1, lRowIndices) > -1){
      CoinPackedVector row;
      start = matStarts[index1];
      end = start + origMat->getVectorSize(index1);
      for(j = start; j < end; j++){
	index2 = matIndices[j];
	pos = binarySearch(0, lCols - 1, index2, lColIndices);
	if(pos > -1)
	  row.insert(pos, matElements[j]);
      }
      newMat->appendRow(row);
     }
  }
    
  bSolver->assignProblem(newMat,
			 varLb, varUb, 
			 objCoeffs,
			 rowLb, rowUb);
  
  bSolver->setObjSense(objSense);
  
  for(i = 0; i < lCols; i++){
    if(solver()->isInteger(lColIndices[i]))
      bSolver->setInteger(i);
  }

  bSolver->initialSolve();

  //bSolver->writeLp("knapsack");

  if(bSolver->isProvenOptimal()){

    /*
      add variables to the interdiction set with the largest reduced costs
      while the budget constraint is still satisfied
    */

    double usedBudget(0.0);
    double * intCost = getInterdictCost();
    double intBudget(getInterdictBudget());
    const double * reducedCost = bSolver->getReducedCost();

    double min_red(bSolver->getInfinity());
    int ind_min_red(-1);
    double cost(0.0);

    int uInd(-1);
    int * fixedVars = new int[lCols];
    CoinZeroN(fixedVars, lCols);

    while(usedBudget < intBudget){
      
      //find the variable with largest reduced cost
      for(i = 0; i < lCols; i++){

	//check this
	uInd = uColIndices[i];
	
	cost = reducedCost[i];

	if(cost < min_red && !fixedVars[uInd]){
	  min_red = cost;
	  ind_min_red = i;
	}

      }

      if((usedBudget + intCost[ind_min_red]) <= intBudget){
	
	//fix the corresponding upper-level variable to 1
	fixedVars[uColIndices[ind_min_red]] = 1;
	usedBudget += intCost[ind_min_red];

	std::cout << "UL var " << uColIndices[ind_min_red] << "fixed." <<std::endl;
	min_red = bSolver->getInfinity();
	ind_min_red = -1;

      }
      else
	break;

    }

    /*
      now we find a feasible solution by fixing upper-level vars
      and solving the lower-level problem
    */

    OsiSolverInterface * lSolver = setUpModel(fixedVars);

    //lSolver->writeLp("findfeas");

    lSolver->branchAndBound();

    double * optUpperSolutionOrd = new double[lCols];
    double * optLowerSolutionOrd = new double[lCols];
    
    CoinZeroN(optUpperSolutionOrd, lCols);
    CoinZeroN(optLowerSolutionOrd, lCols);

    const double * lSol = lSolver->getColSolution();
    double upperObj(0.0);
      
    for(i = 0; i < uCols; i++){
      index1 = uColIndices[i];
      //lpSolution[index] = optUpperSolutionOrd[i];
      upperObj += fixedVars[index1] * 
	solver()->getObjCoefficients()[index1] * solver()->getObjSense();
    }

    for(i = 0; i < lCols; i++){
      index1 = lColIndices[i];
      //lpSolution[index] = optLowerSolutionOrd[i];
      upperObj += lSol[i] * 
	solver()->getObjCoefficients()[index1] * solver()->getObjSense();
    }
   
    bound = upperObj; 
  }

  std::cout << "bound " << bound << std::endl;
  //NEED TO CLEAR A BUNCH OF MEMORY HERE

  return bound;
}

//#############################################################################
double
MibSModel::lowerObjectiveBound()
{

  /* 
     optimize wrt to lower-level objective 
     over current feasible lp feasible region 
  */

  OsiSolverInterface * oSolver = getSolver();
  //OsiSolverInterface * oSolver = solver();
  OsiSolverInterface * hSolver = new OsiCbcSolverInterface(oSolver);

  double bound(oSolver->getInfinity());
  double objSense(getLowerObjSense());  
  int lCols(getLowerDim());
  int uCols(getUpperDim());
  int * lColIndices = getLowerColInd();
  double * lObjCoeffs = getLowerObjCoeffs();

  int tCols(lCols + uCols);

  assert(tCols == oSolver->getNumCols());

  double * nObjCoeffs = new double[tCols];
  int i(0), index(0);
  
  CoinZeroN(nObjCoeffs, tCols);

  for(i = 0; i < lCols; i++){
    index = lColIndices[i];
    nObjCoeffs[index] = lObjCoeffs[i] * objSense;
  }

  hSolver->setObjective(nObjCoeffs);

  //hSolver->writeLp("bound");

  dynamic_cast<OsiCbcSolverInterface *> 
  (hSolver)->getModelPtr()->messageHandler()->setLogLevel(0);

  hSolver->branchAndBound();

  if(hSolver->isProvenOptimal()){

    double upperObjVal(0.0);
    //MibSSolution *mibSol = NULL;
    
    for(i = 0; i < tCols; i++)
      upperObjVal += 
	hSolver->getColSolution()[i] * oSolver->getObjCoefficients()[i];

    //mibSol = new MibSSolution(hSolver->getNumCols(),
    //		      hSolver->getColSolution(),
    //		      upperObjVal,
    //		      this);
    
    //storeSolution(BlisSolutionTypeHeuristic, mibSol);
    //mibSol = NULL;
    bound = upperObjVal;
  }

  return bound;

}

//#############################################################################
void
MibSModel::runPreprocessor()
{

  int uCols(getUpperDim());   // for now, just fixing upper level columns 
  int * upperColInd = getUpperColInd();

  OsiSolverInterface * oSolver = 
    dynamic_cast<OsiClpSolverInterface*>(solver());
  
  oSolver->initialSolve();

  if(oSolver->isProvenOptimal()){

    std::vector<BcpsVariable *> vars = getVariables();

    double lpObjVal = oSolver->getObjValue();
    double objBound = getObjectiveBound();
    setCutoff(objBound * oSolver->getObjSense());

    double objDiff = objBound - lpObjVal;
 
    const double * reducedCost = oSolver->getReducedCost();

    CoinWarmStartBasis * ws = 
       dynamic_cast<CoinWarmStartBasis*>(oSolver->getWarmStart());
    
    const double * sol = oSolver->getColSolution();
    
    int i(0), index(0);
    
    for(i = 0; i < uCols; i++){

      index = upperColInd[i];
      
      if(fabs(reducedCost[index]) >= objDiff){
	
	bool getStatus = true;

	if(oSolver->isInteger(index)){
	  
	  double infeas = fabs(floor(sol[index] + 0.5) - sol[index]);
	  if(infeas > etol_)
	    getStatus = false;
	}

	if(getStatus){

	  if(ws->getStructStatus(index) 
	     == CoinWarmStartBasis::atLowerBound){
	    vars[index]->setUbHard(vars[index]->getLbHard());
	    vars[index]->setUbSoft(vars[index]->getLbSoft());
	  }
	  else if(ws->getStructStatus(index) 
		  == CoinWarmStartBasis::atUpperBound){
	    vars[index]->setLbHard(vars[index]->getUbHard());
	    vars[index]->setLbSoft(vars[index]->getUbSoft());
	  }

	}
      }
    }
  }

}

//#############################################################################
void
MibSModel::runPreprocessor1()
{


  int uCols(getUpperDim());
  int lCols(getLowerDim());
  int * upperColInd = getUpperColInd();
  int * lowerColInd = getLowerColInd();

  OsiSolverInterface * oSolver = 
    dynamic_cast<OsiClpSolverInterface*>(solver());
  
  oSolver->initialSolve();

  if(oSolver->isProvenOptimal()){

    std::vector<BcpsVariable *> vars = getVariables();

    double lpObjVal = oSolver->getObjValue();
    double objBound = getObjectiveBound();
    setCutoff(objBound * oSolver->getObjSense());

    double objDiff = objBound - lpObjVal;
 
    const double * reducedCost = oSolver->getReducedCost();

    CoinWarmStartBasis * ws = 
      dynamic_cast<CoinWarmStartBasis*>(oSolver->getWarmStart());
    
    const double * sol = oSolver->getColSolution();
    
    int i(0), index(0);
    
    for(i = 0; i < uCols; i++){

      index = upperColInd[i];
      
      if(fabs(reducedCost[index]) >= objDiff){
	
	bool getStatus = true;

	if(oSolver->isInteger(index)){
	  
	  double infeas = fabs(floor(sol[index] + 0.5) - sol[index]);
	  if(infeas > etol_)
	    getStatus = false;
	}

	if(getStatus){

	  if(ws->getStructStatus(index) 
	     == CoinWarmStartBasis::atLowerBound){
	    //vars[index]->setUbHard(vars[index]->getLbHard());
	    //vars[index]->setUbSoft(vars[index]->getLbSoft());
	    lpSolver_->setColUpper(index, solver()->getColLower()[index]);
	    varUB_[index] = varLB_[index];
	    std::cout << "variable " << index << "fixed to LB" << std::endl;
	  }
	  else if(ws->getStructStatus(index) 
		  == CoinWarmStartBasis::atUpperBound){
	    //vars[index]->setLbHard(vars[index]->getUbHard());
	    //vars[index]->setLbSoft(vars[index]->getUbSoft());
	    lpSolver_->setColLower(index, solver()->getColUpper()[index]);
	    varLB_[index] =  varUB_[index];
	    std::cout << "variable " << index << "fixed to UB" << std::endl;
	  }

	}
      }
    }
    
    for(i = 0; i < lCols; i++){

      index = lowerColInd[i];
      
      if(fabs(reducedCost[index]) >= objDiff){
	
	bool getStatus = true;

	if(oSolver->isInteger(index)){
	  
	  double infeas = fabs(floor(sol[index] + 0.5) - sol[index]);
	  if(infeas > etol_)
	    getStatus = false;
	}

	//check this - vars arent getting fixed correctly
	if(getStatus){

	  if(ws->getStructStatus(index) 
	     == CoinWarmStartBasis::atLowerBound){
	    //vars[index]->setUbHard(vars[index]->getLbHard());
	    //vars[index]->setUbSoft(vars[index]->getLbSoft());
	    lpSolver_->setColUpper(index, solver()->getColLower()[index]);
	    varUB_[index] = varLB_[index];
	    std::cout << "variable " << index << "fixed to LB" << std::endl;
	  }
	  else if(ws->getStructStatus(index) 
		  == CoinWarmStartBasis::atUpperBound){
	    //vars[index]->setLbHard(vars[index]->getUbHard());
	    //vars[index]->setLbSoft(vars[index]->getUbSoft());
	    lpSolver_->setColLower(index, solver()->getColUpper()[index]);
	    varLB_[index] = varUB_[index];
	    std::cout << "variable " << index << "fixed to LB" << std::endl;
	  }

	}
      }
    }
  }

  //lpSolver_->writeLp("preprocessed");

}

//#############################################################################
CoinPackedVector * 
MibSModel::getSolution()
{

  int varnum = solver()->getNumCols();
  const double *sol = solver()->getColSolution();
  double etol = etol_;
  int *indices = new int[varnum];
  double *values = new double[varnum]; /* n */
  int i(0), cnt(0);
 	
  CoinZeroN(indices, varnum);
  CoinZeroN(values, varnum);

  int index(0);
  int udim(getUpperDim());
  int ldim(getLowerDim());

  int * uIndices = getUpperColInd();
  int * lIndices = getLowerColInd();

  for(i = 0; i < udim; i++){
     index = uIndices[i];
    if (sol[index] > etol || sol[index] < -etol){
      indices[cnt] = index;
      values[cnt++] = sol[index];
    }
  }

  for(i = 0; i < ldim; i++){
     index = lIndices[i];
    if (sol[index] > etol || sol[index] < -etol){
      indices[cnt] = index;
      values[cnt++] = sol[index];
    }
  }

  /*
  for (i = 0; i < varnum; i++){
    if (sol[i] > etol || sol[i] < -etol){
      indices[cnt] = i;
      values[cnt++] = sol[i];
    }
  }
  */

  return(new CoinPackedVector(varnum, cnt, indices, values, false));
}

//#############################################################################
MibSSolType
MibSModel::createBilevel(CoinPackedVector *vec)
{

   MibSSolType solType = bS_->createBilevel(vec, this);

   return solType;
}

//#############################################################################
void 
MibSModel::setBounds()
{

  double * varlower = varLB();
  double * varupper = varUB();
  double * rowlower = conLB();
  double * rowupper = conUB();

  int i(0);

  origColLb_ = new double[numOrigVars_];
  origColUb_ = new double[numOrigVars_];
  origRowLb_ = new double[numOrigCons_];
  origRowUb_ = new double[numOrigCons_];

  for(i = 0; i < numOrigVars_; i++){
    origColLb_[i] = varlower[i];
    origColUb_[i] = varupper[i];
  }

  for(i = 0; i < numOrigCons_; i++){
    origRowLb_[i] = rowlower[i];
    origRowUb_[i] = rowupper[i];
  }

}

//#############################################################################
void 
MibSModel::checkProblemType()
{

  //OsiSolverInterface * lpSolver = getSolver();
  char * colType = colType_;
  int type(MibSPar_->entry(MibSParams::bilevelProblemType));
  int * lColIndices = getLowerColInd();
  int * uColIndices = getUpperColInd();

  int uCols(getUpperDim());
  int lCols(getLowerDim());

  int i(0), ind(0);

  switch(type){
    
  case 0: // general
    {
   
      /********************************************************/
      /* general problem type:                                */
      /* for now, this means the problem must be pure integer */ 
      /********************************************************/
      
      for(i = 0; i < uCols; i++){
	ind = uColIndices[i];
	char type = colType[ind];
	if(type != 'I'){
	  throw CoinError("Instance must be a pure integer problem.",
			  "checkproblem",
			  "MibSModel");
	}
      }
      for(i = 0; i < lCols; i++){
	ind = lColIndices[i];
	char type = colType[ind];
	if(type != 'I'){
	  throw CoinError("Instance must be a pure integer problem.",
			  "checkproblem",
			  "MibSModel");
	}
      }
      std::cout << "Pure integer problem accepted by MibS." << std::endl;
      break;
    }
  case 1: // interdiction
    {

      /********************************************************************/
      /* interdiction problem type:                                       */
      /* for now, this means the lower-level problem must be pure integer */
      /********************************************************************/

      for(i = 0; i < uCols; i++){
	ind = uColIndices[i];
	char type = colType[ind];
	if(type != 'B'){
	  throw CoinError("All upper-level variables must be binary.",
			  "checkproblem",
			  "MibSModel");
	}
      }
      for(i = 0; i < lCols; i++){
	ind = lColIndices[i];
	char type = colType[ind];
	if((type != 'I') && (type != 'B')){
	  throw CoinError("All lower-level variables must be integer.",
			  "checkproblem",
			  "MibSModel");
	}
      }
      std::cout << "Interdiction problem accepted by MibS." << std::endl;
      break;
    }
  }


}

//#############################################################################
void 
MibSModel::setProblemType()
{
  /** Set the type of MIBLP:
   ** UILI - IBLP
   ** UMLM - MIBLP
   ** UBLB - Binary IBLP 
   ** UBLI - upper binary, lower integer
   ** INTD - Interdiction
   **/


}

//#############################################################################
void 
MibSModel::setValFuncSlopes()
{

  //FIXME:  THIS WILL ONLY WORK FOR ONE ROW 

  /** Set the slopes of the lower-level value function **/

  const CoinPackedMatrix * matrix = getSolver()->getMatrixByRow();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();

  int lCols(getLowerDim());  
  int lRows(getLowerRowNum());
  int * lColIndices = getLowerColInd();
  int * lRowIndices = getLowerRowInd();
  double * objCoeffs = getLowerObjCoeffs();
  int i(0), j(0), start(0), end(0), index1(0), index2(0);
  

  std::vector<double> posRowCoeffs; // positive row coeffs  
  std::vector<double> negRowCoeffs; // negative row coeffs  
  std::vector<double> posObjCoeffs; // obj coeffs associated with posRowCoeffs
  std::vector<double> negObjCoeffs; // obj coeffs associated with negRowCoeffs


  // negative and positive slopes of the val func
  double leftSlope(-getSolver()->getInfinity()); 
  double rightSlope(getSolver()->getInfinity());

  int tmp(0);

  for(i = 0; i < lRows; i++){
    //CoinPackedVector row;
     index1 = lRowIndices[i];
     start = matStarts[index1];
     end = start + matrix->getVectorSize(index1);
     for(j = start; j < end; j++){
	index2 = matIndices[j];
	tmp = binarySearch(0, lCols - 1, index2, lColIndices);
      if(tmp > -1)

	//NEED TO CHECK IF VAR IS CONTINUOUS
	if(getSolver()->isContinuous(index2)){
	  if(matElements[j] > 0){ 
	    // add to C+ vector
	    posRowCoeffs.push_back(matElements[j]);
	    posObjCoeffs.push_back(objCoeffs[tmp]);
	  }
	  else{
	    // add to C- vector
	    negRowCoeffs.push_back(matElements[j]);
	    negObjCoeffs.push_back(objCoeffs[tmp]);
	  }
	}
     }
  }

  assert(posRowCoeffs.size() == posObjCoeffs.size());
  assert(negRowCoeffs.size() == negObjCoeffs.size());

  int posSize(posRowCoeffs.size());
  int negSize(negRowCoeffs.size());
  double tmpRatio(0.0);

  for(i = 0; i < posSize; i++){
    tmpRatio = (double) (posObjCoeffs.at(i) / posRowCoeffs.at(i));
    if(tmpRatio < rightSlope)
      rightSlope = tmpRatio;
  }

  for(i = 0; i < negSize; i++){
    tmpRatio = (double) (negObjCoeffs.at(i) / negRowCoeffs.at(i));
    if(tmpRatio > leftSlope)
      leftSlope = tmpRatio;
  }

  rightSlope_ = rightSlope;
  leftSlope_ = leftSlope;

  if(0){
    std::cout << "leftSlope " << leftSlope_ << std::endl;
    std::cout << "rightSlope " << rightSlope_ << std::endl;
  }

}

//#############################################################################
void 
MibSModel::printCurSol()
{

   int uN(upperDim_);
   int lN(lowerDim_);
   int index(0);
   double etol(etol_);
   //double objVal(solver()->getObjValue()); //should add objsense here
   int * upperColInd = getUpperColInd();
   int * lowerColInd = getLowerColInd();
   const double * sol = solver()->getColSolution();

   std::cout << "Nonzero values in current solution" << std::endl;
   std::cout << "**********************************" << std::endl;

   int i(0);

   for(i = 0; i < uN; i++){
      index = upperColInd[i];
      if((sol[index] > etol) || (sol[index] < - etol)) 
	 std::cout << "x[" << i << "]" << sol[index] << std::endl;
   }

   for(i = 0; i < lN; i++){
      index = lowerColInd[i];
      if((sol[index] > etol) || (sol[index] < - etol)) 
	 std::cout << "y[" << i << "]" << sol[index] << std::endl;
   }

   //std::cout << "Upper-level objective value: " << objVal << std::endl;

}


//#############################################################################
int 
MibSModel::binarySearch(int start, int stop, int index, int * indexArray)
{
   int i(0);
   int pos(-1);

  //=========================================================================
  // If the list is short enough, finish with sequential search
  // Otherwise, call binary search recursively
  //=========================================================================

   //FIXME: CHANGED THIS 7/15
   //NOW USES ONLY SEQUENTIAL SEARCH
   //BINARY REQUIRES ORDERED ARRAY
   //OUR ARRAY IS ORDERED BY INTEGERS FIRST, NOT INDICES


   //if((stop - start) < 4){
   if(1){
      for(i = start; i < stop + 1; ++i){
	 if(indexArray[i] == index){
	    pos = i;
	    break;
	 }
      }
      
      return pos;
   }
   else{ 
      
      int midpoint((start + stop)/2);
      int val(indexArray[midpoint]);
      
      if(val == index){
	 pos = midpoint;
	 return pos;
      }
      else if(val > index){
	 pos = binarySearch(start, midpoint - 1, index, indexArray);
      }
      else{
	 pos = binarySearch(midpoint + 1, stop, index, indexArray);
      }
   }
   return pos;
}

//#############################################################################
OsiSolverInterface *
MibSModel::setUpModel(int * fixed)
{

  /** Create lower-level model with fixed upper-level vars **/

  OsiSolverInterface * oSolver = solver();
  OsiSolverInterface * nSolver = new OsiCbcSolverInterface();

  const double * origColLb = getOrigColLb();
  const double * origColUb = getOrigColUb();
  const double * origRowLb = getOrigRowLb();
  const double * origRowUb = getOrigRowUb();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByRow();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();
  //const double * sol = oSolver->getColSolution();   

  double objSense(getLowerObjSense());  
  int lCols(getLowerDim());
  int lRows(getLowerRowNum());
  int * lColIndices = getLowerColInd();
  int * lRowIndices = getLowerRowInd();
  double * lObjCoeffs = getLowerObjCoeffs();
  int i(0), j(0), start(0), end(0), index1(0), index2(0);
  
  double * colUb = new double[lCols];
  double * colLb = new double[lCols];
  double * rowUb = new double[lRows];
  double * rowLb = new double[lRows];
  double * objCoeffs = new double[lCols];
  
  int intCnt(0);
  int * integerVars = new int[lCols];
  CoinFillN(integerVars, lCols, 0);
  
  /** Fill in array of lower-level integer variables **/
  
  for(i = 0; i < lCols; i++){
    index1 = lColIndices[i];
    if(oSolver->isInteger(i)){
      integerVars[intCnt] = i;
      intCnt++;
    }
  }

  //FIXME: NEED TO GET ROW SENSE HERE
  
  double * upComp = new double[lRows];
  CoinFillN(upComp, lRows, 0.0);
  
  for(i = 0; i < lCols; i++){
    colLb[i] = origColLb[lColIndices[i]];
    colUb[i] = origColUb[lColIndices[i]];
  }
  
  CoinDisjointCopyN(lObjCoeffs, lCols, objCoeffs);
 
  /** Set the row bounds **/
  
  for(i = 0; i < lRows; i++){
    rowLb[i] = origRowLb[lRowIndices[i]];
    rowUb[i] = origRowUb[lRowIndices[i]];
  }
  
  /** Get contribution of upper-level columns **/
  
  for(i = 0; i < lRows; i++){
    index1 = lRowIndices[i];
    //index1 = i;
     start = matStarts[index1];
     end = start + matrix->getVectorSize(index1);
     for(j = start; j < end; j++){
	index2 = matIndices[j];
	//if(findIndex(index, lCols, lColIndices) < 0)
	//upComp[i] += matElements[j] * fixed[index2];
	if(binarySearch(0, lCols - 1, index2, lColIndices) < 0)
	  upComp[i] += matElements[j] * fixed[index2];
     }
  }
  
  /** Correct the row bounds to account for fixed upper-level vars **/
  
  for(i = 0; i < lRows; i++){
    rowLb[i] -= upComp[i];
    rowUb[i] -= upComp[i];
  }
  
  
  CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
  newMat->setDimensions(0, lCols);
  int tmp(0);

  for(i = 0; i < lRows; i++){
     CoinPackedVector row;
     index1 = lRowIndices[i];
     start = matStarts[index1];
     end = start + matrix->getVectorSize(index1);
     for(j = start; j < end; j++){
	index2 = matIndices[j];
	//tmp = findIndex(index, lCols, lColIndices);
	tmp = binarySearch(0, lCols - 1, index2, lColIndices);
	if(tmp > -1)
	  row.insert(tmp, matElements[j]);
     }
     newMat->appendRow(row);
  }

  nSolver->assignProblem(newMat, colLb, colUb,
			 objCoeffs, rowLb, rowUb);
  
  
  nSolver->setInteger(integerVars, intCnt);

  nSolver->setObjSense(objSense); //1 min; -1 max

  nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);
  dynamic_cast<OsiCbcSolverInterface *> 
    (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  
  delete [] upComp;
  delete [] integerVars;
  
  return nSolver;

}

//#############################################################################

AlpsReturnStatus 
MibSModel::encodeMibS(AlpsEncoded *encoded) const
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    MibSPar_->pack(*encoded);
    
    return status;
}

//#############################################################################

AlpsReturnStatus 
MibSModel::decodeMibS(AlpsEncoded &encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    MibSPar_->unpack(encoded);
    
    return status;
}

//#############################################################################

AlpsEncoded* 
MibSModel::encode() const 
{ 
    AlpsReturnStatus status = AlpsReturnStatusOk;

    // NOTE: "AlpsKnowledgeTypeModel" is the type name.
    AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeModel);

    status = encodeAlps(encoded);
    status = encodeBcps(encoded);
    status = encodeBlis(encoded);
    status = encodeMibS(encoded);

    return encoded;
}

//#############################################################################

void
MibSModel::decodeToSelf(AlpsEncoded& encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    status = decodeAlps(encoded);
    status = decodeBcps(encoded);
    status = decodeBlis(encoded);
    status = decodeMibS(encoded);
}

//#############################################################################
void
MibSModel::setRequiredFixedList(const CoinPackedMatrix *newMatrix)
{
    int uCols(upperDim_);
    int lRows(lowerRowNum_);
    int * upperColInd = getUpperColInd();
    int * lowerRowInd = getLowerRowInd();

    const double * matElements = newMatrix->getElements();
    const int * matIndices = newMatrix->getIndices();
    const int * matStarts = newMatrix->getVectorStarts();

    int index1, rowIndex, posRow, start, end;
    int i, j;
    int num(0);

    if(!fixedInd_){
	fixedInd_ = new int[numVars_]();
    }

    for(i = 0; i < numVars_; i++){
	fixedInd_[i] = 0;
	if(binarySearch(0, uCols - 1, i, upperColInd) >= 0){
	    start = matStarts[i];
	    end = start + newMatrix->getVectorSize(i);
	    for(j = start; j < end; j++){
		rowIndex = matIndices[j];
		posRow = binarySearch(0, lRows - 1, rowIndex, lowerRowInd);
		if(posRow >= 0){
		    fixedInd_[i] = 1;
		    sizeFixedInd_ ++;
		    break;
		}
	    }
	}
    }

}

//#############################################################################
void                                                                                                                                                       
MibSModel::instanceStructure(const CoinPackedMatrix *newMatrix, const double* rowLB,
			     const double* rowUB, const char *rowSense)
{

    bool printProblemInfo(MibSPar_->entry(MibSParams::printProblemInfo));
    
    /** Determines the properties of instance **/
    if(printProblemInfo == true){
	std::cout<<"======================================="<<std::endl;
	std::cout<<"             Problem Structure          "<<std::endl;
	std::cout<<"======================================="<<std::endl;
        std::cout<<"Number of UL Variables: "<<upperDim_<<std::endl;
        std::cout<<"Number of LL Variables: "<<lowerDim_<<std::endl;
        std::cout<<"Number of UL Rows: "<<upperRowNum_<<std::endl;                                                                                         
        std::cout<<"Number of LL Rows: "<<lowerRowNum_<<std::endl;
    }
    
    int i(0),j(0),index(0), mult(0);
    int numCols(numVars_);
    int numRows(numCons_);
    int uCols(upperDim_);
    int lCols(lowerDim_);
    int uRows(upperRowNum_);
    int lRows(lowerRowNum_);
    int numUpperInt(0);
    int numLowerInt(0);
    int * uColIndices = getUpperColInd();         
    int * lColIndices = getLowerColInd();
    int * lRowIndices = getLowerRowInd();
    char * newRowSense = new char[numRows];

    if (isInterdict_ == false){
	for(i = 0; i < uCols; i++){
	    index = uColIndices[i];
	    if(colType_[index] != 'C'){
		numUpperInt ++;
	    }
	}
	if(printProblemInfo == true){
	    std::cout <<"Number of integer UL Variables: " << numUpperInt << std::endl;
	}
	for(i = 0; i < lCols; i++){
	    index = lColIndices[i];
	    if(colType_[index] != 'C'){
		numLowerInt++;
	    }
	}
	if(printProblemInfo == true){
	    std::cout <<"Number of integer LL Variables: " << numLowerInt << std::endl;
	}
    }


    if (isInterdict_ == false){
	CoinDisjointCopyN(rowSense, numRows, newRowSense);
    }
    else{
	for(i = 0; i < numRows; i++){
	    newRowSense[i] = 'L';
	}
	for(i = 0; i < lRows-lCols; i++){
	    newRowSense[i+uRows] = rowSense[i];
	}
    }
    
    //Checks general or interdiction 
    if(isInterdict_ == true){
	if(printProblemInfo == true){
	    std::cout << "This instance is an interdiction bilevel optimization problem." << std::endl;
	}
    }

    //Checks type of variables
    if(isInterdict_ == false){
	for(i = 0; i < numCols; i++){
	    if(colType_[i] == 'C'){
		isPureInteger_ = false;
		break;
	    }
	}

	for(i = 0; i < numCols; i++){
	    if (colType_[i] != 'B'){
		if(binarySearch(0, lCols - 1, i, lColIndices) < 0){
		    if(fixedInd_[i] == 1){
			allLinkingBin_ = false;
		    }
		    allUpperBin_ = false;
		}
		else{
		    allLowerBin_ = false;
		}
		if((!allUpperBin_) && (!allLowerBin_) && (!allLinkingBin_)){
		    break;
		}
	    }
	}
    }

    if(printProblemInfo == true){
	if(isPureInteger_ == true){
	    std::cout << "This instance is a pure integer bilevel optimization problem" << std::endl;
	}
	else{
	    std::cout << "This instance is a mixed_integer bilevel optimization problem" << std::endl;
	}
    
        if(allUpperBin_ == true){
	    std::cout << "All of UL varibles are binary." << std::endl;
	}
    
        if(allLowerBin_ == true){
	    std::cout << "All of LL varibles are binary." << std::endl; 
	}
    }
    
    int nonZero (newMatrix->getNumElements());
    int counterStart, counterEnd;
    int rowIndex, posRow, posCol;
    double rhs(0.0);
    
    const double * matElements = newMatrix->getElements();
    const int * matIndices = newMatrix->getIndices();           
    const int * matStarts = newMatrix->getVectorStarts();
    //Signs of matrices A1, A2, G1 and G2 are determined
    //based on converting the row senses to 'L'.  
    for(i = 0; i < numCols; i++){
        counterStart = matStarts[i];
        counterEnd = matStarts[i+1];
        for(j = counterStart; j < counterEnd; j++){                                                                  
	    rowIndex = matIndices[j];
	    if(newRowSense[rowIndex] == 'L'){
		mult = 1;
	    }
	    else{
		mult = -1;
	    }
	    posRow = binarySearch(0, lRows - 1, rowIndex, lRowIndices);       
	    posCol = binarySearch(0, lCols - 1, i, lColIndices);
	    if((fabs(matElements[j] - floor(matElements[j])) > etol_) &&
	       (fabs(matElements[j] - ceil(matElements[j])) > etol_)){
		if(posRow < 0){
		    isUpperCoeffInt_ = false;
		}
		else{
		    isLowerCoeffInt_ = false;
		}
	    }
	    if(mult * matElements[j] < -etol_){
		if(posRow < 0){
		    if(posCol < 0){
			positiveA1_ = false;
		    }  
		    else{
                        positiveG1_ = false;
                    }  
                }   
                else{  
                    if(posCol < 0){
                        positiveA2_ = false;
                    }  
                    else{
                        positiveG2_ = false;                                                                                                               
                    }
		}
	    }
	}
    }

    if((isUpperCoeffInt_ == true) || (isLowerCoeffInt_ == true)){
	for(i = 0; i < numRows; i++){
	    switch(newRowSense[i]){
	    case 'L':
		rhs = rowUB[i];
		break;
	    case 'G':
		rhs = rowLB[i];
		break;
	    case 'E':
		std::cout
		    << "MibS cannot currently handle equality constraints." << std::endl;
		abort();
		break;
	    case 'R':
		std::cout
		    << "MibS cannot currently handle range constraints." << std::endl;
		abort();
		break;
	    }
	    if((fabs(rhs - floor(rhs)) > etol_) &&
	       (fabs(rhs - ceil(rhs)) > etol_)){
		posRow = binarySearch(0, lRows - 1, rowIndex, lRowIndices);
		if(posRow < 0){
		    isUpperCoeffInt_ = false;
		}
		else{
		    isLowerCoeffInt_ = false;
		}
	    }
	}
    }

    for(i = 0; i < numCols; i++){
	if(fixedInd_[i] == 1){
	    if(colType_[i] == 'C'){
		std::cout << "All linking variables should be discrete" << std::endl;
		i = -1;
		assert(i > 0);
	    }
	}
    }

    if(printProblemInfo == true){
	if(positiveA1_ == true){
	    std::cout << "Matrix A1 is positive." << std::endl;
	}
    
        if(positiveG1_ == true){
	    std::cout << "Matrix G1 is positive." << std::endl;
	}
    
        if(positiveA2_ == true){
	    std::cout << "Matrix A2 is positive." << std::endl;
	}

	if(positiveG2_ == true){
	    std::cout << "Matrix G2 is positive." << std::endl;
	}
    }

    //set the default value of time limit if it is not set by user
    //Alps sets the default value of this parameter to a very large value
    //but it is very large and creats memory issues for setting the time limit
    //when an extra mip should be solved. So, we set it to a smaller large value.
    if(AlpsPar()->entry(AlpsParams::timeLimit) > 1.0e200){
	AlpsPar()->setEntry(AlpsParams::timeLimit, 1.0e200);
    }
    
    int paramValue(0);
    MibSIntersectionCutType cutType(MibSIntersectionCutTypeNotSet);
    //bool isHypercubeOn(false);

    bool turnOffOtherCuts(MibSPar_->entry(MibSParams::turnOffOtherCuts));
    bool defaultCutIsOn(false);

    /*//Check if hypercube IC is on or off
    paramValue = MibSPar_->entry(MibSParams::useIntersectionCut);
    cutType =  MibSPar_->entry(MibSParams::intersectionCutType);

    if((paramValue == PARAM_ON) && (cutType == 2)){
	isHypercubeOn = true;
    }*/
    
    //Param: "MibS_usePreprocessor" 
    paramValue = MibSPar_->entry(MibSParams::usePreprocessor);
    
    if(paramValue == PARAM_NOTSET)
	MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
    else if(paramValue == PARAM_ON){
	std::cout << "The preprocessor is not currently functional. Automatically disabling your parameter choice." << std::endl;
	MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
    }

    //Param: "MibS_useGreedyHeuristic"
    paramValue = MibSPar_->entry(MibSParams::useGreedyHeuristic);

    if(paramValue == PARAM_NOTSET)
	MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);
    else if(paramValue == PARAM_ON){
	if(isInterdict_ == false){
	    MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);
	    std::cout << "This heuristic only works for interdiction problems. Automatically disabling this heuristic." << std::endl;
	}
    }

    //Param: "MibS_useIncObjCut"
    if((turnOffOtherCuts == true) &&
       (MibSPar_->entry(MibSParams::useIncObjCut) == PARAM_NOTSET)){
	MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useIncObjCut);

    if(paramValue == PARAM_NOTSET)
	MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    else if(paramValue == PARAM_ON){
	if((allLinkingBin_ == false) || (positiveA2_ == false)){
	    std::cout << "The increasing objective cut is not valid for this problem. Automatically disabling this cut." << std::endl;
	    MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
	}
    }

    //Param: "MibS_useInterSectionCut" (MibSIntersectionCutTypeIC)
    paramValue = MibSPar_->entry(MibSParams::useIntersectionCut);
    cutType = static_cast<MibSIntersectionCutType>
	(MibSPar_->entry(MibSParams::intersectionCutType));

    if((paramValue == PARAM_ON) && (cutType == MibSIntersectionCutTypeIC)){
	if((isPureInteger_ == false) || (isLowerCoeffInt_ == false)){
	    std::cout << "The intersection cut IC is not valid problem.  Automatically disabling this cut."
		      << std::endl;
	    MibSPar()->setEntry(MibSParams::useIntersectionCut, PARAM_NOTSET);
	}
    }    

    //Param: "MibS_useNoGoodCut"
    paramValue = MibSPar_->entry(MibSParams::useNoGoodCut);

    if(paramValue == PARAM_NOTSET){
	MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    }
    else if((paramValue == PARAM_ON) && (allUpperBin_ == false)){
	std::cout << "The no-good cut is not valid for this problem. Automatically disabling this cut." << std::endl;
	MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    }

    //Param: "MibS_useBendersCut"
    if((turnOffOtherCuts == true) &&
       (MibSPar_->entry(MibSParams::useBendersCut) == PARAM_NOTSET)){
	MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useBendersCut);

    if(paramValue == PARAM_NOTSET){
	if((isInterdict_ == false) || (positiveG2_ == false)){
	    MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	}else{
	    MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_ON);
	    MibSPar()->setEntry(MibSParams::bendersCutType, MibSBendersCutTypeJustOneCut);
	}
    }
    else if(paramValue == PARAM_ON){
	if((isInterdict_ == false) || (positiveG2_ == false)){
	    std::cout << "The benders cut is not valid for this problem. Automatically disabling this cut." << std::endl;
	    MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	}
    }
    if(MibSPar_->entry(MibSParams::useBendersCut) == PARAM_ON){
	    defaultCutIsOn = true;
	}

    //Param: "MibS_useGeneralNoGoodCut"
    if((turnOffOtherCuts == true) &&
       (MibSPar_->entry(MibSParams::useNoGoodCut) == PARAM_NOTSET)){
	MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useGeneralNoGoodCut);

    if(paramValue == PARAM_NOTSET){
	if(allLinkingBin_ == false){
	    MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
	}
	else if(defaultCutIsOn == false){
	    MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_ON);
	}
    }
    else if((paramValue == PARAM_ON) && (allLinkingBin_ == false)){
	std::cout << "Generalized no-good cut is not valid for this problem. Automatically disabling this cut." << std::endl;
	MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    }
    if(MibSPar_->entry(MibSParams::useGeneralNoGoodCut) == PARAM_ON){
	defaultCutIsOn = true;
    }

    //Param: "MibS_usePureIntegerCut"
    if((turnOffOtherCuts == true) &&
       (MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_NOTSET)){
	MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::usePureIntegerCut);

    if(paramValue == PARAM_NOTSET){
	if((isPureInteger_ == false) || (isUpperCoeffInt_ == false)
	   || (isLowerCoeffInt_ == false)){
	    MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
	}
	else if(defaultCutIsOn == false){
	    MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_ON);
	}
    }
    else if((paramValue == PARAM_ON) && ((isPureInteger_ == false) || (isUpperCoeffInt_ == false)
					 || (isLowerCoeffInt_ == false))){
	std::cout << "The pure integer cut is not valid for this problem. Automatically disabling this cut." << std::endl;
	MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    }
    if(MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_ON){
	defaultCutIsOn = true;
    }

    //Param: "MibS_useInterSectionCut" (MibSIntersectionCutTypeHypercubeIC) 
    if((turnOffOtherCuts == true) &&
       (MibSPar_->entry(MibSParams::useIntersectionCut) == PARAM_NOTSET)){
	MibSPar()->setEntry(MibSParams::useIntersectionCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useIntersectionCut);

    if(paramValue == PARAM_NOTSET){
	if(defaultCutIsOn == false){
	    MibSPar()->setEntry(MibSParams::useIntersectionCut, PARAM_ON);
	    MibSPar()->setEntry(MibSParams::intersectionCutType,
				MibSIntersectionCutTypeHypercubeIC);
	}
    }
    
    //Param: "MibS_branchProcedure"
    MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
	  (MibSPar_->entry(MibSParams::branchStrategy));
    if(branchPar == MibSBranchingStrategyNotSet){
	if((isInterdict_ == true) || (numUpperInt <= numLowerInt)){
	    MibSPar()->setEntry(MibSParams::branchStrategy,
				MibSBranchingStrategyLinking);
	    if(printProblemInfo == true){
		std::cout <<
		    "Since no branching procedure is selected, it is set to 'MibSBranchingStrategyLinking' automatically."
			  << std::endl;
	    }
	}
	else{
	    MibSPar()->setEntry(MibSParams::branchStrategy,
				MibSBranchingStrategyFractional);
	    if(printProblemInfo == true){
		std::cout <<
		    "Since no branching procedure is selected, it is set to 'MibSBranchingStrategyFractional' automatically."
			  << std::endl;
	    }
	}
    }


    if(printProblemInfo == true){
	if(MibSPar_->entry(MibSParams::useIncObjCut) == PARAM_ON){
	    std::cout << "'Increasing objective cut' genertor is on." << std::endl;
	}

	if(MibSPar_->entry(MibSParams::useBendersCut) == PARAM_ON){
	    if(MibSPar_->entry(MibSParams::bendersCutType) ==
	       MibSBendersCutTypeJustOneCut){
		std::cout << "'Benders cut' genertor (just one cut) is on." << std::endl;
	    }
	    else{
		std::cout << "'Benders cut' genertor (multiple cuts) is on." << std::endl;
	    }
	}

        if(MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_ON){
	    std::cout << "'Pure integer cut' genertor is on."<< std::endl;
	}

	if(MibSPar_->entry(MibSParams::useNoGoodCut) == PARAM_ON){
	    std::cout << "'No-good cut' genertor is on."<< std::endl;
	}

	if(MibSPar_->entry(MibSParams::useGeneralNoGoodCut) == PARAM_ON){
	    std::cout << "'Generalized no-good cut' genertor is on."<< std::endl;
	}

        if(MibSPar_->entry(MibSParams::useIntersectionCut) == PARAM_ON){
	    if(MibSPar_->entry(MibSParams::intersectionCutType) ==
	       MibSIntersectionCutTypeIC){
		std::cout << "'Intersection cut IC' generator is on." << std::endl;
	    }
	    if(MibSPar_->entry(MibSParams::intersectionCutType) ==
	       MibSIntersectionCutTypeWatermelon){
		std::cout << "'watermelon IC' generator is on." << std::endl;
	    }
	    if(MibSPar_->entry(MibSParams::intersectionCutType) ==
	       MibSIntersectionCutTypeHypercubeIC){
		std::cout << "'hypercube IC' generator is on." << std::endl;
	    }
	}

        if(MibSPar_->entry(MibSParams::branchStrategy) ==
	   MibSBranchingStrategyLinking){
	    std::cout << "Branching procedure is set to 'MibSBranchingStrategyLinking'."
		      << std::endl;
	}
	else if(MibSPar_->entry(MibSParams::branchStrategy) ==
		MibSBranchingStrategyFractional){
	    std::cout << "Branching procedure is set to 'MibSBranchingStrategyFractional'."
		      << std::endl;
	}
    }
	
    //Setting parameters of solving (VF) and (UB)
    int solveSecondLevelWhenXYVarsInt(MibSPar_->entry(MibSParams::
					    solveSecondLevelWhenXYVarsInt));
    int solveSecondLevelWhenXVarsInt(MibSPar_->entry(MibSParams::
					   solveSecondLevelWhenXVarsInt));
    int solveSecondLevelWhenLVarsInt(MibSPar_->entry(MibSParams::
					      solveSecondLevelWhenLVarsInt));
    int solveSecondLevelWhenLVarsFixed(MibSPar_->entry(MibSParams::
						solveSecondLevelWhenLVarsFixed));
    int computeBestUBWhenXVarsInt(MibSPar_->entry(MibSParams::
					      computeBestUBWhenXVarsInt));
    int computeBestUBWhenLVarsInt(MibSPar_->entry(MibSParams::
						  computeBestUBWhenLVarsInt));
    int computeBestUBWhenLVarsFixed(MibSPar_->entry(MibSParams::
						    computeBestUBWhenLVarsFixed));

    if((solveSecondLevelWhenXYVarsInt == PARAM_NOTSET) &&
       (solveSecondLevelWhenXVarsInt == PARAM_NOTSET) &&
       (solveSecondLevelWhenLVarsInt == PARAM_NOTSET) &&
       (solveSecondLevelWhenLVarsFixed == PARAM_NOTSET) &&
       (computeBestUBWhenXVarsInt == PARAM_NOTSET) &&
       (computeBestUBWhenLVarsInt == PARAM_NOTSET) &&
       (computeBestUBWhenLVarsFixed == PARAM_NOTSET)){
	MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt, PARAM_ON);
	MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed, PARAM_ON);
	MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed, PARAM_ON);
    }

    if(printProblemInfo == true){
	if(MibSPar_->entry(MibSParams::solveSecondLevelWhenXYVarsInt) == PARAM_ON){
	    std::cout << "'solveSecondLevelWhenXYVarsInt' is true." << std::endl;
	}
    
        if(MibSPar_->entry(MibSParams::solveSecondLevelWhenXVarsInt) == PARAM_ON){
	    std::cout << "'solveSecondLevelWhenXVarsInt' is true." <<std::endl;
	}

        if(MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsInt) == PARAM_ON){
	    std::cout << "'solveSecondLevelWhenLVarsInt' is true." <<std::endl;
	}

        if(MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsFixed) == PARAM_ON){
	    std::cout << "'solveSecondLevelWhenLVarsFixed' is true." <<std::endl;
	}

        if(MibSPar_->entry(MibSParams::computeBestUBWhenXVarsInt) == PARAM_ON){
	    std::cout << "'computeBestUBWhenXVarsInt' is true." <<std::endl;
	}

        if(MibSPar_->entry(MibSParams::computeBestUBWhenLVarsInt) == PARAM_ON){
	    std::cout << "'computeBestUBWhenLVarsInt' is true." <<std::endl;
	}

        if(MibSPar_->entry(MibSParams::computeBestUBWhenLVarsFixed) == PARAM_ON){
	    std::cout << "'computeBestUBWhenLVarsFixed' is true." <<std::endl;
	}
    }

    //Setting "saveSeenLinkingSols" parameter
    
    if(MibSPar_->entry(MibSParams::useLinkingSolutionPool) == PARAM_NOTSET){
	MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
    }

    if(printProblemInfo == true){
	if(MibSPar_->entry(MibSParams::useLinkingSolutionPool) == PARAM_ON){
	    std::cout << "'seenLinkingSolutions' pool is used." << std::endl;
	}
	else{
	    std::cout << "'seenLinkingSolutions' pool is not used." << std::endl;
	
	}
    }
	
    delete [] newRowSense;
}
