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

//#include "BlisModel.h"
#include "MibSConfig.hpp"
#include "MibSModel.hpp"
#ifdef COIN_HAS_SYMPHONY
#include "symphony.h"
#include "SymConfig.h"
#include "OsiSymSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CPLEX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif

#include "MibSSolution.hpp"
#include "MibSCutGenerator.hpp"
//#include "MibSBilevel.hpp"
#include "MibSTreeNode.hpp"
#include "MibSConstants.hpp"	

#include "MibSBranchStrategyMaxInf.hpp"
#include "MibSBranchStrategyPseudo.hpp"
#include "MibSBranchStrategyStrong.hpp"

//FIXME::RELIABILITY BRANCHING DOESNT WORK
//NECESSARY DATA MEMBERS ARE DESIGNATED AS PRIVATE
//IN PARENT CLASS.  DIDNT WANT TO ALTER BLIS CODE
//#include "MibSBranchStrategyRel.h"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

#ifdef COIN_HAS_C11
#include <random>
#endif

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
  if(stocA2Matrix_) delete stocA2Matrix_;
    
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
  truncLowerDim_ = 0;
  lowerObjSense_ = 0.0;
  upperDim_ = 0;
  leftSlope_ = 0;
  rightSlope_ = 0;
  lowerRowNum_ = 0;
  truncLowerRowNum_ = 0;
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
  boundProbRoot_ = NULL;
  boundProbLeafNum_ = 0;
  numScenarios_ = 0;
  stocA2Matrix_ = NULL;
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
MibSModel::readAuxiliaryData(const CoinPackedMatrix& rowMatrix,
			     const double* conLB, const double* conUB,
			     int numCols, int numRows, double infinity,
			     const char *rowSense)
{

 std::string stochasticityType(MibSPar_->entry
			       (MibSParams::stochasticityType));

 int isA2Random(MibSPar_->entry(MibSParams::isA2Random));

 int isSMPSFormat(MibSPar_->entry(MibSParams::isSMPSFormat));

 if((stochasticityType == "stochasticWithoutSAA") && (isSMPSFormat == PARAM_OFF)){
     throw CoinError("When the SAA method is not used(for the stochastic case), the input format should be SMPS",
		     "readAuxiliaryData",
		     "MibSModel");
 }

 if(stochasticityType == "deterministic"){
     numScenarios_ = 1;
     std::string fileName = getLowerFile();
     if (fileName == ""){
	 std::string mpsFile(getUpperFile());
	 int length = mpsFile.length();
	 char *tmpArr = new char[length + 1];
	 fileName = mpsFile.erase(length - 3, 3);
	 fileName.append("aux");
	 if (fileCoinReadable(fileName)){
	     fileName.copy(tmpArr, length);
	     tmpArr[length]='\0';
	     MibSPar()->setEntry(MibSParams::auxiliaryInfoFile, tmpArr);
	     std::cout << "Warning: The auxiliary file was not specified. ";
	     std::cout << std::endl;
	     std::cout << "MibS used " <<  fileName << " automatically.";
	     std::cout << std::endl;
	 }else{
	     fileName = mpsFile.erase(length - 3, 3);
	     fileName.append("txt");
	     if (fileCoinReadable(fileName)){
		 fileName.copy(tmpArr, length);
		 tmpArr[length]='\0';
		 MibSPar()->setEntry(MibSParams::auxiliaryInfoFile, tmpArr);
		 std::cout << "Warning: The auxiliary file is not specified. ";
		 std::cout << "MibS selected " <<  fileName << " automatically.";
		 std::cout << std::endl;
	     }else{
		 std::cout << "Error: The auxiliary file is not specified. ";
		 std::cout << "Aborting." << std::endl;
		 abort();
	     }
	 }
	 delete [] tmpArr;
     }
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
	     setTruncLowerDim(iValue);
	 }
         else if(key == "M"){
	     data_stream >> iValue;
	     setLowerRowNum(iValue);
	     setTruncLowerRowNum(iValue);
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
 else if(stochasticityType == "stochasticWithoutSAA"){//the format is smps necessarily
     std::string timFileName = getTimFile();
     fileCoinReadable(timFileName);

     std::string stoFileName = getStoFile();
     fileCoinReadable(stoFileName);

     std::ifstream timData_stream(timFileName.c_str());
     std::ifstream stoData_stream(stoFileName.c_str());

     if(!timData_stream){
	 std::cout << "Error opening input time data file. Aborting.\n";
	 abort();
     }

     if(!stoData_stream){
	 std::cout << "Error opening input stochastic data file. Aborting.\n";
	 abort();
     }

     int i(0), j(0), p(0);
     double etol(etol_);
     std::string key, fLColName, fLRowName;
     int iValue(0), posRow(-1), posCol(-1), index(0), rowNumElem(0), scCount(0), start(0);
     double dValue(0.0), element(0.0);
     int uColNum(0), lColNum(0), uRowNum(0), lRowNum(0), totalRowNum(0);
     int truncRowNum(0), truncLColNum(0), truncLRowNum(0);
     bool isFirstScenario(true);
     CoinPackedVector appendedRow;
     CoinPackedVector row;
     int *rowInd = NULL;
     double *rowElem = NULL;
     CoinPackedVector col;
     std::vector<double> stocLowerRhs;

     //read time file
     //ignoring first three lines
     for(i = 0; i < 7; i++){
	 timData_stream >> key;
     }

     //name of first ll variable
     timData_stream >> key;
     fLColName = key;

     //name of first ll row
     timData_stream >> key;
     fLRowName = key;

     for(p = 0; p < numCols; p++){
	 if(columnName_[p] == fLColName){
	     posCol = p;
	     break;
	 }
     }
     if(posCol < 0){
	 std::cout << fLColName << " does not belong to the list of variables." << std::endl;
	 abort();
     }
     else{
	 uColNum = p;
	 truncLColNum = numCols - uColNum;
	 setTruncLowerDim(truncLColNum);
     }

     posRow = -1;
     for(p = 0; p < numRows; p++){
	 if(rowName_[p] == fLRowName){
	     posRow = p;
	     break;
	 }
     }
     if(posRow < 0){
	 std::cout << fLRowName << " does not belong to the list of constraints." << std::endl;
	 abort();
     }
     else{
	 uRowNum = p;
	 truncLRowNum = numRows - uRowNum;
	 setTruncLowerRowNum(truncLRowNum);
     }
     timData_stream >> key;

     //set lower objective function coeffs and sense
     lowerObjCoeffs_ = new double[truncLColNum];
     for(i = 0; i < truncLColNum; i++){
	 timData_stream >> key;
	 timData_stream >> dValue;
	 lowerObjCoeffs_[i] = dValue;
     }
     timData_stream >> key;
     timData_stream >> iValue;
     lowerObjSense_ = iValue;//1 min; -1 max

     //saharSto: do we need fill lowerColInd_ and lowerRowInd_?
     //set lower col and row indices
     lowerColInd_ = new int[truncLColNum];
     lowerRowInd_ = new int[truncLRowNum];
     CoinIotaN(lowerColInd_, truncLColNum, uColNum);
     CoinIotaN(lowerRowInd_, truncLRowNum, uRowNum);

     timData_stream.close();

     //read stochastic data
     //find number of scenarios
     while(stoData_stream >> key){
	 if(key == "SC"){
	     numScenarios_ ++;
	 }
     }
     stoData_stream.close();

     lColNum = truncLColNum * numScenarios_;
     lRowNum = truncLRowNum * numScenarios_;
     totalRowNum = uRowNum + lRowNum;
     truncRowNum = uRowNum + truncLRowNum;
     setLowerDim(lColNum);
     setLowerRowNum(lRowNum);
     //saharSto3: check it
     //stoData_stream(stoFileName.c_str());
     stoData_stream.open(stoFileName.c_str());

     origRowLb_ = new double[totalRowNum];
     origRowUb_ = new double[totalRowNum];
     CoinFillN(origRowLb_, totalRowNum, -1 * infinity);
     CoinFillN(origRowUb_, totalRowNum, infinity);
     
     CoinPackedMatrix *stocMatrixA2 = 0;
     stocMatrixA2 = new CoinPackedMatrix(false, 0, 0);
     stocMatrixA2->setDimensions(0, uColNum);

     //saharSto3: check it
     //double **fullMatrixA2 = NULL;
     //double **fullMatrixA2Copy = NULL;
     std::vector<std::vector<double> > fullMatrixA2;
     std::vector<std::vector<double> > fullMatrixA2Copy;

     if(isA2Random != PARAM_OFF){
	 /*fullMatrixA2 = new double*[truncLRowNum];
	 fullMatrixA2Copy = new double*[truncLRowNum];
	 for(i = 0; i < truncLRowNum; i++){
	     //saharSto3: check if they are filled with zeros
	     fullMatrixA2[i] = new double[uColNum]();
	     fullMatrixA2Copy[i] = new double[uColNum]();
	     }*/
	 fullMatrixA2.resize(truncLRowNum, std::vector<double>(uColNum, 0.0));
	 fullMatrixA2Copy.resize(truncLRowNum, std::vector<double>(uColNum, 0.0));
     }

     for(i = uRowNum; i < numRows; i++){
	 row = rowMatrix.getVector(i);
	 rowInd = row.getIndices();
	 rowElem = row.getElements();
	 rowNumElem = row.getNumElements();
	 for(j = 0; j < rowNumElem; j++){
	     index = rowInd[j];
	     element = rowElem[j];
	     if(index < uColNum){
		 appendedRow.insert(index, element);
		 if(isA2Random != PARAM_OFF){
		     fullMatrixA2[i - uRowNum][index] = element;
		     fullMatrixA2Copy[i - uRowNum][index] = element;
		 }
	     }
	 }
	 //saharSto: when A2 is not random, we just store matrix A2
	 //as stocMatrixA2 because we know that it is the same
	 //for all scenarios
	 if(isA2Random == PARAM_OFF){
	     stocMatrixA2->appendRow(appendedRow);
	 }
	 appendedRow.clear();
     }

     //ignore first two lines
     for(i = 0; i < 4; i++){
	 stoData_stream >> key;
     }

     while(stoData_stream >> key){
	 if(key == "SC"){
	     stoData_stream >> key;
	     stoData_stream >> key;
	     stoData_stream >> dValue;
	     scenarioProb_.push_back(dValue);
	     stoData_stream >> key;
	     start = uRowNum + scCount * truncLRowNum;
	     CoinDisjointCopyN(conLB + uRowNum, truncLRowNum, origRowLb_ + start);
	     CoinDisjointCopyN(conUB + uRowNum, truncLRowNum, origRowUb_ + start);
	     if((isA2Random != PARAM_OFF) && (isFirstScenario == false)){
		 //start = scCount * truncLRowNum;
		 for(i = 0; i < truncLRowNum; i++){
		     for(j = 0; j < uColNum; j++){
			 element = fullMatrixA2Copy[i][j];
			 if(fabs(element) > etol){
			     appendedRow.insert(j,element); 
			 }
			 fullMatrixA2Copy[i][j] = fullMatrixA2[i][j];
		     }
		     stocMatrixA2->appendRow(appendedRow);
		     appendedRow.clear();
		 }
	     }
	     scCount ++;
	     if(isFirstScenario == true){
		 isFirstScenario = false;
	     }
	 }
	 else if(key == "ENDATA"){
	     if(isA2Random != PARAM_OFF){
		 for(i = 0; i < truncLRowNum; i++){
		     for(j = 0; j < uColNum; j++){
			 element = fullMatrixA2Copy[i][j];
			 if(fabs(element) > etol){
			     appendedRow.insert(j,element);
			 }
		     }
		     stocMatrixA2->appendRow(appendedRow);
		     appendedRow.clear();
		 }
	     }
	     break;
	 }
	 else if(key == "RHS"){
	     stoData_stream >> key;
	     stoData_stream >> dValue;
	     posRow = -1;
	     for(p = 0; p < numRows; p++){
		 if(rowName_[p] == key){
		     posRow = p;
		     break;
		 }
	     }
	     if(posRow < 0){
		 std::cout << key << " does not belong to the list of constraints." << std::endl;
		 throw CoinError("Wrong constriant name",
				 "readAuxiliaryData",
				 "MibSModel");
	     }
	     else if(posRow < uRowNum){
		 throw CoinError("Upper-level RHS cannot be random",
				 "readAuxiliaryData",
				 "MibSModel");
	     }
	     else{
		 index = (scCount - 1) * truncLRowNum + posRow;
		 if(rowSense[posRow] == 'L'){
		     origRowUb_[index] = dValue;
		 }
		 else{
		     origRowLb_[index] = dValue;
		 }
	     }  
	 }
	 else{
	     posCol = -1;
	     for(p = 0; p < numCols; p++){
		 if(columnName_[p] == key){
		     posCol = p;
		     break;
		 }
	     }
	     if(posCol < 0){
		 std::cout << key << " does not belong to the list of variables." << std::endl;
		 throw CoinError("Wrong variable name",
				 "readAuxiliaryData",
				 "MibSModel");
	     }
	     else if(posCol >= uColNum){
		 throw CoinError("Lower-level coefficients cannot be random",
				 "readAuxiliaryData",
				 "MibSModel");
	     }
	     else{
		 stoData_stream >> key;
		 stoData_stream >> dValue;
		 posRow = -1;
		 for(p = 0; p < numRows; p++){
		     if(rowName_[p] == key){
			 posRow = p;
			 break;
		     }
		 }
		 if(posRow < 0){
		     std::cout << key << " does not belong to the list of constraints." << std::endl;
		     throw CoinError("Wrong constriant name",
				     "readAuxiliaryData",
				     "MibSModel");
		 }
		 else if(posRow < uRowNum){
		     throw CoinError("Coefficients of lower-level constraints cannot be random",
				     "readAuxiliaryData",
				     "MibSModel");
		 }
		 else{
		     fullMatrixA2Copy[posRow - uRowNum][posCol] = dValue;
		 }
	     }
	 }
     }

     setStocA2Matrix(stocMatrixA2);

     stoData_stream.close();

     std::cout << "Time Data File: " << getTimFile() << "\n";
     std::cout << "Stochastic Data File: " << getStoFile() << "\n";
     std::cout << "Number of LL Variables:   "
	       << truncLColNum << "\n\n";
 }
 else{
     int iValue(0), pos(-1);
     double dValue(0.0);
     std::string key, fLColName, fLRowName;;
     int i(0), k(0), p(0);
     int uColNum(0), uRowNum(0);
     int truncLColNum(0), truncLRowNum(0);

     if(isSMPSFormat != PARAM_ON){
	 std::string fileName = getLowerFile();
         fileCoinReadable(fileName);
         std::ifstream data_stream(fileName.c_str());

         if(!data_stream){
	     std::cout << "Error opening input data file. Aborting.\n";
	     abort();
	 }

         while(data_stream >> key){
	     if(key == "N"){
		 data_stream >> iValue;
	         truncLColNum = iValue;
	         setTruncLowerDim(iValue);
	     }
	     else if(key == "M"){
		 data_stream >> iValue;
	         truncLRowNum = iValue;
	         setTruncLowerRowNum(iValue);
	     }
	     else if(key == "LC"){
		 data_stream >> iValue;
	     }
	     else if(key == "LR"){
		 data_stream >> iValue;
	     }
	     else if(key == "LO"){
		 if(!getLowerObjCoeffs())
		     lowerObjCoeffs_ = new double[getTruncLowerDim()];
		 data_stream >> dValue;
	         lowerObjCoeffs_[k] = dValue;
	         k++;
	     }
	     else if(key == "OS"){
		 data_stream >> dValue;
	         lowerObjSense_ = dValue; //1 min; -1 max
	     }
	 }
	 uColNum = numCols - truncLColNum;
	 uRowNum = numRows - truncLRowNum;
	 
	 data_stream.close();
     }
     else{
	 std::string timFileName = getTimFile();
	 fileCoinReadable(timFileName);

	 std::ifstream timData_stream(timFileName.c_str());

	 if(!timData_stream){
	     std::cout << "Error opening input time data file. Aborting.\n";
	     abort();
	 }

	 //read time file
	 //ignoring first three lines
	 for(i = 0; i < 7; i++){
	     timData_stream >> key;
	 }

	 //name of first ll variable
	 timData_stream >> key;
	 fLColName = key;

	 //name of first ll row
	 timData_stream >> key;
	 fLRowName = key;

	 for(p = 0; p < numCols; p++){
	     if(columnName_[p] == fLColName){
		 pos = p;
		 break;
	     }
	 }
	 if(pos < 0){
	     std::cout << fLColName << " does not belong to the list of variables." << std::endl;
	     abort();
	 }
	 else{
	     uColNum = p;
	     truncLColNum = numCols - uColNum;
	     setTruncLowerDim(truncLColNum);
	 }

	 pos = -1;
	 for(p = 0; p < numRows; p++){
	     if(rowName_[p] == fLRowName){
		 pos = p;
		 break;
	     }
	 }
	 if(pos < 0){
	     std::cout << fLRowName << " does not belong to the list of constraints." << std::endl;
	     abort();
	 }
	 else{
	     uRowNum = p;
	     truncLRowNum = numRows - uRowNum;
	     setTruncLowerRowNum(truncLRowNum);
	 }
	 timData_stream >> key;
	 
	 //set lower objective function coeffs and sense
	 if(!getLowerObjCoeffs())
	     lowerObjCoeffs_ = new double[truncLColNum];
     
	 for(i = 0; i < truncLColNum; i++){
	     timData_stream >> key;
	     timData_stream >> dValue;
	     lowerObjCoeffs_[i] = dValue;
	 }
	 timData_stream >> key;
	 timData_stream >> iValue;
	 lowerObjSense_ = iValue;//1 min; -1 max

	 timData_stream.close();
     }

     //saharSto:in the saa format, it is assumed that
     //The indices of upper-level variables are 0,...,n1-1 and
     //the indices of lower-level variables are n1,...,n1+n2-1
     //and the same for upper- and lower-level constraints
     //saharSto: we do not need LC and LR, but I considered them in case we
     //want to use the input files for deterministic problems
     if(!getLowerColInd())
	 lowerColInd_ = new int[truncLColNum];
     CoinIotaN(lowerColInd_, truncLColNum, uColNum);

     if(!getLowerRowInd())
	 lowerRowInd_ = new int[truncLRowNum];
     CoinIotaN(lowerRowInd_, truncLRowNum, uRowNum);

     if(isSMPSFormat == PARAM_ON){
	 std::cout << "Time Data File: " << getTimFile() << "\n";
	 std::cout << "Stochastic Data File: " << getStoFile() << "\n";
	 std::cout << "Number of LL Variables:   "
		   << truncLColNum << "\n\n";
	 
     }
     else{
	 std::cout << "LL Data File: " << getLowerFile() << "\n";
	 std::cout << "Number of LL Variables:   "
		   << getTruncLowerDim() << "\n\n";
     }
 }
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
   CoinPackedMatrix rowMatrix = *(mps->getMatrixByRow());

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

   //readAuxiliaryData(numCols, numRows); // reads in lower-level vars, rows, obj coeffs
   readAuxiliaryData(rowMatrix, conLB, conUB, numCols, numRows,
		     mps->getInfinity(), rowSense); // reads in lower-level vars, rows, obj coeffs

   std::string stochasticityType(MibSPar_->entry
				 (MibSParams::stochasticityType));

   int useProgresHedg(MibSPar_->entry(MibSParams::useProgresHedg));

   if(stochasticityType == "stochasticWithSAA"){
       setupSAA(matrix, rowMatrix, varLB, varUB, objCoef, conLB, conUB, colType,
		objSense, numCols, numRows, mps->getInfinity(), rowSense);
   }

   if(useProgresHedg == true){
       int uColNum = numCols - getTruncLowerDim();
       for(j = 0; j < uColNum; j++){
	   if(colType[j] != 'B'){
	       throw CoinError("The progressive hedging heuristic can be used only when all UL variables are binary",
			       "readProblemData", "MibSModel");
	   }
       }

       if(stochasticityType != "stochasticWithoutSAA"){
	   throw CoinError("The progressive hedging heuristic can be used only when the parameter stochasticityType is set to stochasticWithoutSAA",
			   "readProblemData", "MibSModel");
       }
       else{
	   setupProgresHedg(matrix, rowMatrix, varLB, varUB, objCoef, conLB, conUB, colType,
			    objSense, numCols, numRows, mps->getInfinity(), rowSense);
       }
   }
   
   
   loadProblemData(matrix, rowMatrix, varLB, varUB, objCoef, conLB, conUB, colType, 
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
//this function is used when user selects to use progressive hedging heuristic
//it is assumed that all subproblems are solved untile finding at least one feasible solution
void
MibSModel::setupProgresHedg(const CoinPackedMatrix& matrix,
			    const CoinPackedMatrix& rowMatrix,
			    const double* varLB, const double* varUB,
			    const double* objCoef, const double* conLB,
			    const double* conUB, const char * colType,
			    double objSense, int truncNumCols, int truncNumRows,
			    double infinity, const char *rowSense)
{
    const int clockType = AlpsPar()->entry(AlpsParams::clockType);
    broker_->timer().setClockType(clockType);
    broker_->subTreeTimer().setClockType(clockType);
    broker_->tempTimer().setClockType(clockType);
    broker_->timer().start();
    
    //saharPH: define these as parameter later
    int maxPHIteration(MibSPar_->entry(MibSParams::iterationLimitPH));
    //double rho(1.0);
    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));

    int i(0), j(0), k(0);
    int cntIter(0), index(0), rowNumElem(0);
    double element(0.0), optObj(0.0), value(0.0), removedObj(0.0);
    double maxVal(0.0), minVal(0.0);
    double etol(etol_);
    bool shouldTerminate(false), isTimeLimReached(false);
    bool isImplementable(true), isFixed(true), solveRestrictedProb(false);
    bool maxFound(false), minFound(false);
    int truncLColNum(getTruncLowerDim());
    int truncLRowNum(getTruncLowerRowNum());
    int uColNum(truncNumCols - truncLColNum);
    int uRowNum(truncNumRows - truncLRowNum);
    int numScenarios(getNumScenarios());
    std::vector<double> scenarioProb(getScenarioProb());
    CoinPackedVector row;
    CoinPackedVector appendedRow;
    int *rowInd = NULL;
    double *rowElem = NULL;
    double *optSol = NULL;

    double *rho = new double[uColNum];
    CoinZeroN(rho, uColNum);
    

    double **wArrs = new double *[numScenarios];
    for(i = 0; i < numScenarios; i++){
	wArrs[i] = new double[uColNum]();
    }

    double *tmpWArr = new double[uColNum];
    CoinZeroN(tmpWArr, uColNum);

    double *implemSol = new double[uColNum];
    CoinZeroN(implemSol, uColNum);

    double *implemSolBack = new double[uColNum];
    CoinZeroN(implemSolBack, uColNum);

    int *newIndexUL = new int[uColNum];
    CoinZeroN(newIndexUL, uColNum);
    int numFixed(0);

    double **allScenarioSols = new double *[numScenarios];
    for(i = 0; i < numScenarios; i++){
	allScenarioSols[i] = new double[uColNum]();
    }

    double *uLUpper = NULL;
    double *uLLower = NULL;
    double *scenarioSol = NULL;

    double *copyObjCoef = new double[truncNumCols];
    if(objSense == 1){
	memcpy(copyObjCoef, objCoef, sizeof(double) * truncNumCols);
    }
    else{
	for(i = 0; i < truncNumCols; i++){
	    copyObjCoef[i] = -1 * objCoef[i];
	}
    }

    CoinPackedMatrix *stocA2Matrix(getStocA2Matrix());
    CoinPackedMatrix *filledCoefMatrix = 0;
    if(isA2Random == PARAM_OFF){
	filledCoefMatrix = new CoinPackedMatrix(false, 0, 0);
	filledCoefMatrix->setDimensions(0, truncNumCols);
	for(i = 0; i < uRowNum; i++){
	    filledCoefMatrix->appendRow(rowMatrix.getVector(i));
	}
	for(i = 0; i < truncLRowNum; i++){
	    appendedRow = stocA2Matrix->getVector(i);
	    row = rowMatrix.getVector(i + uRowNum);
	    rowInd = row.getIndices();
	    rowElem = row.getElements();
	    rowNumElem = row.getNumElements();
	    for(j = 0; j < rowNumElem; j++){
		index = rowInd[j];
		if(index >= uColNum){
		    appendedRow.insert(index, rowElem[j]);
		}
	    }
	    filledCoefMatrix->appendRow(appendedRow);
	    appendedRow.clear();
	}
    }

    while(shouldTerminate == false){
	memcpy(implemSolBack, implemSol, sizeof(double) * uColNum);
	//CoinZeroN(foundSol, numScenarios);
	CoinZeroN(implemSol, uColNum);
	isImplementable = true;
	//numFoundSol = 0;
	if(cntIter > 0){
	    //compute rho
	    for(j = 0; j < uColNum; j++){
		maxVal = 0.0;
		minVal = 1.0;
		maxFound = false;
		minFound = false;
		for(k = 0; k < numScenarios; k++){
		    value = allScenarioSols[k][j];
		    if(maxFound == false){
			if(fabs(value - maxVal) > 0.5){
			    maxVal = 1.0;
			    maxFound = true;
			}
		    }
		    if(minFound == false){
			if(fabs(minVal - value) > 0.5){
			    minVal = 0.0;
			    minFound = true;
			}
		    }
		    if((minFound == true) && (maxFound == true)){
			break;
		    }
		}
		rho[j] = copyObjCoef[j]/(maxVal - minVal + 1);
	    }
	}
	for(i = 0; i < numScenarios; i ++){
	    isTimeLimReached = false;
	    if(cntIter > 0){
		//compute w
		for(j = 0; j < uColNum; j++){
		    tmpWArr[j] = wArrs[i][j] + rho[j] * (allScenarioSols[i][j] - implemSolBack[j]);
		    wArrs[i][j] = tmpWArr[j];
		}
	    }

	    std::cout << "Solving scenario " << i << " of iteration " << cntIter << std::endl;
	    scenarioSol = solvePHProb(rowMatrix,varLB, varUB, copyObjCoef, conLB,
				      conUB, colType, objSense, truncNumCols, truncNumRows,
				      infinity, rowSense, i, cntIter, isTimeLimReached,
				      tmpWArr, implemSolBack, rho, filledCoefMatrix);

	    if(isTimeLimReached == true){
		std::cout << "Time limit is reached" << std::endl;
		goto TERM_SETUPPH;
	    }

	    if(scenarioSol == NULL){
		throw CoinError("It is assumed that all subproblems are solved until finding at least one feasible solution",
				"setupProgresHedg", "MibSModel");
	    }

	    //if(scenarioSol != NULL){
	    //foundSol[i] = 1;
	    //numFoundSol ++;
	    for(j = 0; j < uColNum; j++){
		value = scenarioSol[j];
		allScenarioSols[i][j] = value;
		implemSol[j] += scenarioProb[i] * value;
	    }
	    delete [] scenarioSol;
		//}
	}

	/*for(j = 0; j < uColNum; j++){
	    implemSol[j] = implemSol[j]/numFoundSol;
	    }*/

	cntIter ++;
	//if(numFoundSol == numScenarios){
	for(j = 0; j < uColNum; j++){
	    value = implemSol[j];
	    for(i = 0; i < numScenarios; i++){
		if(fabs(allScenarioSols[i][j] - value) > etol){
		    isImplementable = false;
		    break;
		}
	    }
	    if(isImplementable == false){
		break;
	    }
	}
	    //}
	    //else{
	    //isImplementable = false;
	    //}

	if(isImplementable == true){
	    optObj = 0;
	    for(i = 0; i < uColNum; i++){
		value = floor(implemSol[i] + 0.5);
		//saharPH: correct it
		optObj += value * objCoef[i];
	    }
	    std::cout << "CPU time for solving subproblems: " << broker_->timer().getCpuTime() << std::endl;
	    std::cout << "Wallclock time for solving subproblems: " << broker_->timer().getWallClock() << std::endl;
	    printSolutionPH(implemSol, optObj, cntIter, uColNum);
	    shouldTerminate = true;
	}
	else if(cntIter == maxPHIteration){
	    solveRestrictedProb = true;
	    shouldTerminate = true;
	    if(uLUpper == NULL){
		uLUpper = new double[uColNum];
	    }
	    if(uLLower == NULL){
		uLLower = new double[uColNum];
	    }
	    //store the fixed part
	    removedObj = 0.0;
	    numFixed = 0;
	    for(j = 0; j < uColNum; j++){
		value = implemSol[j];
		isFixed = true;
		for(i = 0; i < numScenarios; i++){
		    //if(foundSol[i] == 1){
		      if(fabs(allScenarioSols[i][j] - value) > etol){
			  uLUpper[j] = 1;
			  uLLower[j] = 0;
			  isFixed = false;
			  newIndexUL[j] = j - numFixed;
			  break;
		      }
			//}
		}
		/*if(isFixed == false){
  if((implemSol[j] - 0.45) < etol){
    isFixed = true;
    value = 0.0;
  }
  else if((0.55 - implemSol[j]) < etol){
    isFixed = true;
    value = 1.0;
  }
  }*/
		if(isFixed == true){
		    newIndexUL[j] = -1;
		    numFixed ++;
		    value = floor(value + 0.5);
		    uLUpper[j] = value;
		    uLLower[j] = value;
		    removedObj += copyObjCoef[j] * value;
		    std:: cout << "x[" << j << "] is fixed to " << value << std::endl;
		}
	    }//
	}
    }

 TERM_SETUPPH:
    delete [] copyObjCoef;
    delete [] implemSol;
    delete [] implemSolBack;
    delete [] rho;
    delete [] tmpWArr;
    for(i = 0; i < numScenarios; i++){
	delete [] wArrs[i];
	delete [] allScenarioSols[i];
    }
    delete [] wArrs;
    delete [] allScenarioSols;
    if(filledCoefMatrix){
	delete filledCoefMatrix;
    }
    if(solveRestrictedProb == true){
	std::cout << "Removed obj = " << removedObj << std::endl;
	std::cout << "CPU time for solving subproblems: " << broker_->timer().getCpuTime() << std::endl;
	std::cout << "Wallclock time for solving subproblems: " << broker_->timer().getWallClock() << std::endl;
	solveRestrictedPH(matrix, rowMatrix, varLB, varUB, objCoef,
			  conLB, conUB, colType, objSense, truncNumCols,
			  truncNumRows, infinity, rowSense, uLUpper, uLLower, numFixed, newIndexUL);
	std::cout << "Total CPU time: " << broker_->timer().getCpuTime() << std::endl;
	std::cout << "Total wallclock time: " << broker_->timer().getWallClock() << std::endl;
    }

    exit(0);

}

//#############################################################################
void
MibSModel::solveRestrictedPH(const CoinPackedMatrix& matrix,
			     const CoinPackedMatrix& rowMatrix,
			     const double* varLB, const double* varUB,
			     const double* objCoefOrig, const double* conLB,
			     const double* conUB, const char * colTypeOrig,
			     double objSense, int truncNumColsOrig, int truncNumRows,
			     double infinity, const char *rowSense,
			     double *uLUpper, double *uLLower, int numFixed, int* newIndexUL)
{

    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));

    std::string feasCheckSolver
	= MibSPar_->entry(MibSParams::feasCheckSolver);

    double timeLimit(AlpsPar()->entry(AlpsParams::timeLimit));

    int clockType(AlpsPar()->entry(AlpsParams::clockType));

    bool useUBDecompose(MibSPar_->entry(MibSParams::useUBDecompose));

    int i(0), j(0), k(0);
    int index(0);
    double value(0.0), rowNumElem(0.0);
    double etol(etol_);
    int truncLColNum(getTruncLowerDim());
    int truncLRowNum(getTruncLowerRowNum());
    int uColNumOrig(truncNumColsOrig - truncLColNum);
    int uColNum(uColNumOrig - numFixed);
    int uRowNum(truncNumRows - truncLRowNum);
    int lColNum(truncLColNum * numScenarios_);
    int lRowNum(truncLRowNum * numScenarios_);
    int truncNumCols(uColNum + truncLColNum);
    int numCols(uColNum + lColNum);
    int numRows(uRowNum + lRowNum);
    CoinPackedVector row;
    CoinPackedVector appendedRow;
    int *rowInd = NULL;
    double *rowElem = NULL;
    double *optSol = NULL;

    double remainingTime = timeLimit - broker_->timer().getTime();

    if(remainingTime < etol){
	std::cout << "Time limit is reached" << std::endl;
	return;
    }

    double *colLower = new double[truncNumCols];
    double *colUpper = new double[truncNumCols];
    double *objCoef = new double[truncNumCols];
    char *colType = new char[truncNumCols];


    //CoinDisjointCopyN(uLLower, uColNum, colLower);
    //CoinDisjointCopyN(uLUpper, uColNum, colUpper);

    if(numFixed > 0){
    for(i = 0; i < uColNumOrig; i++){
	index = newIndexUL[i];
	if(index >= 0){
	    colLower[index] = uLLower[i];
	    colUpper[index] = uLUpper[i];
	    objCoef[index] = objCoefOrig[i];
	    colType[index] = colTypeOrig[i];
	}
    }
    CoinDisjointCopyN(varLB + uColNumOrig, truncLColNum, colLower + uColNum);
    CoinDisjointCopyN(varUB + uColNumOrig, truncLColNum, colUpper + uColNum);

    CoinDisjointCopyN(objCoefOrig + uColNumOrig, truncLColNum, objCoef + uColNum);

    CoinDisjointCopyN(colTypeOrig + uColNumOrig, truncLColNum, colType + uColNum);
    }
    else{
	memcpy(colLower, varLB, sizeof(double) * truncNumCols);
	memcpy(colUpper, varUB, sizeof(double) * truncNumCols);
	memcpy(objCoef, objCoefOrig, sizeof(double) * truncNumCols);
	memcpy(colType, colTypeOrig, sizeof(char) * truncNumCols);
    }

    //Create new MibS model to solve bilevel
    MibSModel *modelResPH = new MibSModel();

    modelResPH->numScenarios_ = getNumScenarios();
    modelResPH->setTruncLowerDim(truncLColNum);
    modelResPH->setTruncLowerRowNum(truncLRowNum);
    modelResPH->lowerObjCoeffs_ = new double[truncLColNum];
    memcpy(modelResPH->lowerObjCoeffs_, getLowerObjCoeffs(),
	   sizeof(double) * truncLColNum);

    modelResPH->setLowerObjSense(getLowerObjSense());

    modelResPH->lowerColInd_ = new int[truncLColNum];
    modelResPH->lowerRowInd_ = new int[truncLRowNum];
    CoinIotaN(modelResPH->lowerColInd_, truncLColNum, uColNum);
    CoinIotaN(modelResPH->lowerRowInd_, truncLRowNum, uRowNum);
    //modelResPH->lowerColInd_ = lowerColInd_;
    //modelResPH->lowerRowInd_ = lowerRowInd_;

    modelResPH->setLowerDim(lColNum);
    modelResPH->setLowerRowNum(lRowNum);

    modelResPH->origRowLb_ = new double[numRows];
    modelResPH->origRowUb_ = new double[numRows];
    memcpy(modelResPH->origRowLb_, origRowLb_, sizeof(double) * numRows);
    memcpy(modelResPH->origRowUb_, origRowUb_, sizeof(double) * numRows);
    //modelResPH->origRowLb_ = origRowLb_;
    //modelResPH->origRowUb_ = origRowUb_;

    CoinPackedMatrix *stocMatrixA2 = new CoinPackedMatrix(false, 0, 0);
    if(numFixed > 0){
        if(isA2Random != PARAM_OFF){
	    for(i = 0; i < lRowNum; i++){
		row = stocA2Matrix_->getVector(i);
	        rowInd = row.getIndices();
	        rowElem = row.getElements();
	        rowNumElem = row.getNumElements();
	        for(j = 0; j < rowNumElem; j++){
		    index = newIndexUL[rowInd[j]];
		    if(index >= 0){
			appendedRow.insert(index, rowElem[j]);
		    }
		    else{
			index = i + uRowNum;
		        value = rowElem[j] * uLLower[rowInd[j]];
		        if(rowSense[index] == 'G'){
			    modelResPH->origRowLb_[index] = origRowLb_[index] - value;
			}
		        else{
			    modelResPH->origRowUb_[index] = origRowUb_[index] - value;
			}
		    }
		}
	        stocMatrixA2->appendRow(appendedRow);
	        appendedRow.clear();
	    }
	}
        else{
	    for(i = 0; i < truncLRowNum; i++){
		row = stocA2Matrix_->getVector(i);
	        rowInd = row.getIndices();
	        rowElem = row.getElements();
	        rowNumElem = row.getNumElements();
	        for(j = 0; j < rowNumElem; j++){
		    index = newIndexUL[rowInd[j]];
		    if(index >= 0){
			appendedRow.insert(index, rowElem[j]);
		    }
		    else{
			value = rowElem[j] * uLLower[rowInd[j]];
		        if(rowSense[i + uRowNum] == 'G'){
			    for(k = 0; k < numScenarios_; k++){
				index = i + uRowNum + k * truncLRowNum;
			        modelResPH->origRowLb_[index] = origRowLb_[index] - value;
			    }
			}
		        else{
			    for(k = 0; k < numScenarios_; k++){
				index = i + uRowNum + k * truncLRowNum;
			        modelResPH->origRowUb_[index] = origRowUb_[index] - value;
			    }
			}
		    }
		}
	        stocMatrixA2->appendRow(appendedRow);
	        appendedRow.clear();
	    }
	}

	modelResPH->setStocA2Matrix(stocMatrixA2);
    }
    else{
	modelResPH->setStocA2Matrix(stocA2Matrix_);
    }

    CoinPackedMatrix *coefMatrix = new CoinPackedMatrix(false, 0, 0);
    CoinPackedMatrix *coefRowMatrix = new CoinPackedMatrix(false, 0, 0);
    if(numFixed > 0){
        coefMatrix->setDimensions(0, truncNumCols);
        coefRowMatrix->setDimensions(0, truncNumCols);

        for(i = 0; i < uRowNum; i++){
	    row = rowMatrix.getVector(i);
	    rowInd = row.getIndices();
	    rowElem = row.getElements();
	    rowNumElem = row.getNumElements();
	    for(j = 0; j < rowNumElem; j++){
		index = newIndexUL[rowInd[j]];
	        if(index >= 0){
		    appendedRow.insert(index, rowElem[j]);
		}
	        else{
		    value = rowElem[j] * uLLower[rowInd[j]];
		    if(rowSense[i] == 'G'){
			modelResPH->origRowLb_[i] = origRowLb_[i] - value;
		    }
		    else{
			modelResPH->origRowUb_[i] = origRowUb_[i] - value;
		    }
		}
	    }
	    coefMatrix->appendRow(appendedRow);
	    coefRowMatrix->appendRow(appendedRow);
	    appendedRow.clear();
	}


        for(i = 0; i < truncLRowNum; i++){
	    appendedRow = stocMatrixA2->getVector(i);
	    row = rowMatrix.getVector(i + uRowNum);
	    rowInd = row.getIndices();
	    rowElem = row.getElements();
	    rowNumElem = row.getNumElements();
	    for(j = 0; j < rowNumElem; j++){
		index = rowInd[j];
	        if(index >= uColNumOrig){
		    appendedRow.insert(index - numFixed, rowElem[j]);
		}
	    }
	    coefMatrix->appendRow(appendedRow);
	    coefRowMatrix->appendRow(appendedRow);
	    appendedRow.clear();
	}
        coefMatrix->reverseOrdering();
    }

    //stocMatrixA2->setDimensions(0, uColNum);
    //for(i = 0; i < lRowNum; i++){
    //stocA2Matrix->(stocA2Matrix_->getVector(i));
    //}

    modelResPH->scenarioProb_ = getScenarioProb();

    remainingTime = timeLimit - broker_->timer().getTime();

    if(remainingTime < etol){
	std::cout << "Time limit is reached" << std::endl;
	return;
    }
    remainingTime = CoinMax(remainingTime, 0.0);

    //Set up lp solver
    OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);
    modelResPH->setSolver(&lpSolver);

    modelResPH->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
    //modelResPH->AlpsPar()->setEntry(AlpsParams::nodeLimit, 150);
    modelResPH->BlisPar()->setEntry(BlisParams::heurStrategy, 0);
    modelResPH->MibSPar()->setEntry(MibSParams::feasCheckSolver, feasCheckSolver.c_str());
    //modelResPH->MibSPar()->setEntry(MibSParams::printProblemInfo, false);
    modelResPH->MibSPar()->setEntry(MibSParams::useBoundCut, false);
    modelResPH->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyLinking);
    modelResPH->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
    modelResPH->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useTypeWatermelon, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_ON);
    modelResPH->MibSPar()->setEntry(MibSParams::useTypeTenderIC, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useTypeHybridIC, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::useNewPureIntCut, false);
    modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt, PARAM_ON);
    modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXVarsInt, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsInt, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed, PARAM_ON);
    modelResPH->MibSPar()->setEntry(MibSParams::computeBestUBWhenXVarsInt, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsInt, PARAM_OFF);
    modelResPH->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed, PARAM_ON);
    modelResPH->MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);

    modelResPH->isInterdict_ = false;
    modelResPH->MibSPar()->setEntry(MibSParams::stochasticityType, "stochasticWithoutSAA");
    modelResPH->MibSPar()->setEntry(MibSParams::isA2Random, isA2Random);
    modelResPH->MibSPar()->setEntry(MibSParams::useUBDecompose, useUBDecompose);
    modelResPH->AlpsPar()->setEntry(AlpsParams::clockType, clockType);

    if(numFixed > 0){
	modelResPH->loadProblemData(*coefMatrix, *coefRowMatrix, colLower, colUpper, objCoef,
				    conLB, conUB, colType, objSense, infinity, rowSense);
    }
    else{
	modelResPH->loadProblemData(matrix, rowMatrix, colLower, colUpper, objCoef,
				    conLB, conUB, colType, objSense, infinity, rowSense);
    }

    int argc = 1;
    char** argv = new char* [1];
    argv[0] = (char *) "mibs";

#ifdef  COIN_HAS_MPI
    AlpsKnowledgeBrokerMPI brokerResPH(argc, argv, *modelResPH);
#else
    AlpsKnowledgeBrokerSerial brokerResPH(argc, argv, *modelResPH);
#endif

    brokerResPH.search(modelResPH);

    brokerResPH.printBestSolution();

    delete [] colLower;
    delete [] colUpper;
    delete [] objCoef;
    delete [] colType;
    if(coefMatrix){
	delete coefMatrix;
    }
    if(coefRowMatrix){
	delete coefRowMatrix;
    }
    
    delete modelResPH;
}


//#############################################################################
/*void
MibSModel::solveRestrictedPH(const CoinPackedMatrix& matrix,
     const CoinPackedMatrix& rowMatrix,
     const double* varLB, const double* varUB,
     const double* objCoef, const double* conLB,
     const double* conUB, const char * colType,
     double objSense, int truncNumCols, int truncNumRows,
     double infinity, const char *rowSense,
     double *uLUpper, double *uLLower, int numFixed, int *isFixedUL)
{
  int isA2Random(MibSPar_->entry(MibSParams::isA2Random));
  
  std::string feasCheckSolver
    = MibSPar_->entry(MibSParams::feasCheckSolver);
  double timeLimit(AlpsPar()->entry(AlpsParams::timeLimit));
  int clockType(AlpsPar()->entry(AlpsParams::clockType));
  bool useUBDecompose(MibSPar_->entry(MibSParams::useUBDecompose));
  double etol(etol_);
  int truncLColNum(getTruncLowerDim());
  int truncLRowNum(getTruncLowerRowNum());
  int uColNum(truncNumCols - truncLColNum);
  int uRowNum(truncNumRows - truncLRowNum);
  int lColNum(truncLColNum * numScenarios_);
  int lRowNum(truncLRowNum * numScenarios_);
  int numCols(uColNum + lColNum);
  int numRows(uRowNum + lRowNum);
  double remainingTime = timeLimit - broker_->subTreeTimer().getTime();
  if(remainingTime < etol){
  std::cout << "Time limit is reached" << std::endl; 
    return;
  }
  double *colLower = new double[truncNumCols];
  double *colUpper = new double[truncNumCols];
  CoinDisjointCopyN(uLLower, uColNum, colLower);
  CoinDisjointCopyN(uLUpper, uColNum, colUpper);
  CoinDisjointCopyN(varLB + uColNum, truncLColNum, colLower + uColNum);
  CoinDisjointCopyN(varUB + uColNum, truncLColNum, colUpper + uColNum);
  
  // Create new MibS model to solve bilevel 
  MibSModel *modelResPH = new MibSModel();
  modelResPH->numScenarios_ = getNumScenarios();
  modelResPH->setTruncLowerDim(truncLColNum);
  modelResPH->setTruncLowerRowNum(truncLRowNum);
  modelResPH->lowerObjCoeffs_ = new double[truncLColNum];
  memcpy(modelResPH->lowerObjCoeffs_, getLowerObjCoeffs(),
 sizeof(double) * truncLColNum);
  
  modelResPH->setLowerObjSense(getLowerObjSense());
  //modelResPH->lowerColInd_ = new int[truncLColNum];
  //modelResPH->lowerRowInd_ = new int[truncLRowNum];
  //CoinIotaN(modelResPH->lowerColInd_, truncLColNum, uColNum);
  //CoinIotaN(modelResPH->lowerRowInd_, truncLRowNum, uRowNum);
  modelResPH->lowerColInd_ = lowerColInd_;
  modelResPH->lowerRowInd_ = lowerRowInd_;
  
  modelResPH->setLowerDim(lColNum);
  modelResPH->setLowerRowNum(lRowNum);
  //modelResPH->origRowLb_ = new double[numRows];
  //modelResPH->origRowUb_ = new double[numRows];
  //memcpy(modelResPH->origRowLb_, origRowLb_, sizeof(double) * numRows);
  //memcpy(modelResPH->origRowUb_, origRowUb_, sizeof(double) * numRows);
  modelResPH->origRowLb_ = origRowLb_;
  modelResPH->origRowUb_ = origRowUb_; 
  //CoinPackedMatrix *stocMatrixA2 = new CoinPackedMatrix(false, 0, 0);
  //stocMatrixA2->setDimensions(0, uColNum);
  //for(i = 0; i < lRowNum; i++){
  //stocA2Matrix->(stocA2Matrix_->getVector(i));
  //}
  modelResPH->setStocA2Matrix(getStocA2Matrix());
  modelResPH->scenarioProb_ = getScenarioProb();
  remainingTime = timeLimit - broker_->subTreeTimer().getTime();
  if(remainingTime < etol){
  std::cout << "Time limit is reached" << std::endl;
    return;
  }
  remainingTime = CoinMax(remainingTime, 0.0);
  
  // Set up lp solver 
  OsiClpSolverInterface lpSolver;
  lpSolver.getModelPtr()->setDualBound(1.0e10);
  lpSolver.messageHandler()->setLogLevel(0);
  modelResPH->setSolver(&lpSolver);
  modelResPH->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
  modelResPH->BlisPar()->setEntry(BlisParams::heurStrategy, 0);
  modelResPH->MibSPar()->setEntry(MibSParams::feasCheckSolver, feasCheckSolver.c_str());
  modelResPH->MibSPar()->setEntry(MibSParams::printProblemInfo, false);
  modelResPH->MibSPar()->setEntry(MibSParams::useBoundCut, false);
  modelResPH->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyLinking);
  modelResPH->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
  modelResPH->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useTypeWatermelon, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_ON);
  modelResPH->MibSPar()->setEntry(MibSParams::useTypeTenderIC, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useTypeHybridIC, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::useNewPureIntCut, false);
  modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt, PARAM_ON);
  modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXVarsInt, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsInt, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed, PARAM_ON);
  modelResPH->MibSPar()->setEntry(MibSParams::computeBestUBWhenXVarsInt, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsInt, PARAM_OFF);
  modelResPH->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed, PARAM_ON);
  modelResPH->MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
  modelResPH->isInterdict_ = false;
  modelResPH->MibSPar()->setEntry(MibSParams::stochasticityType, "stochasticWithoutSAA");
  modelResPH->MibSPar()->setEntry(MibSParams::isA2Random, isA2Random);
  modelResPH->MibSPar()->setEntry(MibSParams::useUBDecompose, useUBDecompose);
  modelResPH->AlpsPar()->setEntry(AlpsParams::clockType, clockType);
  modelResPH->loadProblemData(matrix, rowMatrix, colLower, colUpper, objCoef,
      conLB, conUB, colType, objSense, infinity, rowSense);
  int argc = 1;
  char** argv = new char* [1];
  argv[0] = (char *) "mibs";
#ifdef  COIN_HAS_MPI
  AlpsKnowledgeBrokerMPI brokerResPH(argc, argv, *modelResPH);
#else
  AlpsKnowledgeBrokerSerial brokerResPH(argc, argv, *modelResPH);
#endif
  brokerResPH.search(modelResPH);
  brokerResPH.printBestSolution();
  delete modelResPH;
  }*/

//#############################################################################
void
MibSModel::printSolutionPH(double *optSol, double optObj, int numIter,
			   int uColNum)
{

    int i(0);
    double nearInt(0.0);
    std::cout <<  "============================================" << std::endl;
    std::cout << "Number of PH heuristic iteration: " << numIter << std::endl;
    std::cout << "Best solution: " << std::endl;
    std::cout << "Cost = " << optObj << std::endl;
    for(i = 0; i < uColNum; i++){
	if(optSol[i] > 1.0e-15 || optSol[i] < -1.0e-15){
	    nearInt = floor(optSol[i] + 0.5);
	    if(ALPS_FABS(nearInt - optSol[i]) < 1.0e-6){
		std::cout << "x[" << i << "] = " << nearInt << std::endl;
	    }
	    else{
		std::cout << "x[" << i << "] = " << optSol[i] << std::endl;
	    }
	}
    }

}

//#############################################################################
double *
MibSModel::solvePHProb(const CoinPackedMatrix& rowMatrix, const double *varLB,
		       const double *varUB, double *origObjCoef,
		       const double *conLB, const double *conUB, const char *colType,
		       double objSense, int numCols, int numRows, double infinity,
		       const char *rowSense, int scenarioIndex, int iterIndex,
		       bool &isTimeLimReached, double *wArr, double *implemSol, double *rho,
		       CoinPackedMatrix *filledCoefMatrix)
{

    int nodeLimit(MibSPar_->entry(MibSParams::nodeLimitPHSubprob));

    double gapLimit(MibSPar_->entry(MibSParams::optimalRelGapLimitPHSubprob));

    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));

    double timeLimit(AlpsPar()->entry(AlpsParams::timeLimit));

    std::string feasCheckSolver =
	MibSPar_->entry(MibSParams::feasCheckSolver);

    int clockType(AlpsPar()->entry(AlpsParams::clockType));

    int i(0), j(0);
    int index(0), rowNumElem(0);
    double remainingTime(0.0);
    double etol(etol_);
    int lColNum(getTruncLowerDim());
    int lRowNum(getTruncLowerRowNum());
    int uColNum(numCols - lColNum);
    int uRowNum(numRows - lRowNum);
    CoinPackedMatrix *stocA2Matrix(getStocA2Matrix());
    CoinPackedVector row;
    CoinPackedVector appendedRow;
    int *rowInd = NULL;
    double *rowElem = NULL;
    double *optSol = NULL;

    double *objCoef = new double[numCols];
    CoinZeroN(objCoef, numCols);

    double *rowLower = new double[numRows];
    double *rowUpper = new double[numRows];
    memcpy(rowLower, conLB, sizeof(double) * numRows);
    memcpy(rowUpper, conUB, sizeof(double) * numRows);

    int *uColInd = new int[uColNum];
    int *uRowInd = new int[uRowNum];
    int *structRowInd = new int[numRows];
    CoinIotaN(uColInd, uColNum, 0);
    CoinIotaN(uRowInd, uRowNum, 0);
    CoinIotaN(structRowInd, numRows, 0);

    CoinPackedMatrix *coefMatrix = 0;
    CoinPackedMatrix *coefRowMatrix = 0;
    CoinPackedMatrix copyCoefMatrix;
    CoinPackedMatrix copyCoefRowMatrix;

    if(isA2Random != PARAM_OFF){
	coefMatrix = new CoinPackedMatrix(false, 0, 0);
        coefMatrix->setDimensions(0, numCols);
        coefRowMatrix = new CoinPackedMatrix(false, 0, 0);
        coefRowMatrix->setDimensions(0, numCols);

        for(i = 0; i < uRowNum; i++){
	    coefMatrix->appendRow(rowMatrix.getVector(i));
	    coefRowMatrix->appendRow(rowMatrix.getVector(i));
	}

        for(i = 0; i < lRowNum; i++){
	    if(isA2Random != PARAM_OFF){
		appendedRow = stocA2Matrix->getVector(scenarioIndex * lRowNum + i);
	    }
	    row = rowMatrix.getVector(i + uRowNum);
	    rowInd = row.getIndices();
	    rowElem = row.getElements();
	    rowNumElem = row.getNumElements();
	    for(j = 0; j < rowNumElem; j++){
		index = rowInd[j];
		if(index >= uColNum){
		    appendedRow.insert(index, rowElem[j]);
		}
	    }
	    coefMatrix->appendRow(appendedRow);
	    coefRowMatrix->appendRow(appendedRow);
	    appendedRow.clear();
	}
        coefMatrix->reverseOrdering();
    }
    else{
	copyCoefMatrix.reverseOrderedCopyOf(*filledCoefMatrix);
	copyCoefRowMatrix.copyOf(*filledCoefMatrix);
    }

    //There is no penalty in the first iteration
    double rhoVal(0.0);
    memcpy(objCoef, origObjCoef, sizeof(double) * numCols);
    if(iterIndex > 0){
	for(i = 0; i < uColNum; i++){
	    //saharPH: check rho/2.0
	    rhoVal = rho[i];
	    objCoef[i] += wArr[i] - rhoVal * implemSol[i] + (rhoVal/2.0);
	}
    }

    //setting row bounds
    for(i = uRowNum; i < numRows; i++){
	index = scenarioIndex * lRowNum + i;
	if(rowSense[i] == 'L'){
	    rowLower[i] = -1 * infinity;
	    rowUpper[i] = origRowUb_[index];
	}
	else{
	    rowLower[i] = origRowLb_[index];
	    rowUpper[i] = infinity;
	}
    }

    /** Set up lp solver **/
    OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);

    //Create new MibS model to solve bilevel
    MibSModel *pHModel = new MibSModel();

    remainingTime = timeLimit - broker_->timer().getTime();
    if(remainingTime <= etol){
	isTimeLimReached = true;
	//saharPH: free the memory
	return NULL;
    }
    remainingTime = CoinMax(remainingTime, 0.00);
    pHModel->setSolver(&lpSolver);
    pHModel->AlpsPar()->setEntry(AlpsParams::msgLevel, -1);
    pHModel->BlisPar()->setEntry(BlisParams::optimalRelGap, gapLimit);
    pHModel->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
    //pHModel->AlpsPar()->setEntry(AlpsParams::timeLimit, 5);
    pHModel->AlpsPar()->setEntry(AlpsParams::nodeLimit, nodeLimit);
    pHModel->BlisPar()->setEntry(BlisParams::heurStrategy, 0);
    pHModel->MibSPar()->setEntry(MibSParams::feasCheckSolver, feasCheckSolver.c_str());
    pHModel->MibSPar()->setEntry(MibSParams::printProblemInfo, false);
    pHModel->MibSPar()->setEntry(MibSParams::useBoundCut, false);
    //saharPH: it is dependent on the problem
    //pHModel->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyLinking);
    pHModel->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
    pHModel->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useTypeWatermelon, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_ON);
    pHModel->MibSPar()->setEntry(MibSParams::useTypeTenderIC, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useTypeHybridIC, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::useNewPureIntCut, false);
    pHModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt, PARAM_ON);
    pHModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXVarsInt, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsInt, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed, PARAM_ON);
    pHModel->MibSPar()->setEntry(MibSParams::computeBestUBWhenXVarsInt, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsInt, PARAM_OFF);
    pHModel->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed, PARAM_ON);
    pHModel->MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
    pHModel->isInterdict_ = false;
    pHModel->MibSPar()->setEntry(MibSParams::stochasticityType, "deterministic");
    pHModel->AlpsPar()->setEntry(AlpsParams::clockType, clockType);

    pHModel->loadAuxiliaryData(lColNum, lRowNum, getLowerColInd(), getLowerRowInd(),
			       getLowerObjSense(), getLowerObjCoeffs(), uColNum,
			       uRowNum, uColInd, uRowInd, numRows, structRowInd, 0,
			       NULL, NULL, NULL);

    if(isA2Random != PARAM_OFF){
	pHModel->loadProblemData(*coefMatrix, *coefRowMatrix, varLB, varUB, objCoef,
			         rowLower, rowUpper, colType, 1, infinity,
			         rowSense);
    }
    else{
	pHModel->loadProblemData(copyCoefMatrix, copyCoefRowMatrix, varLB, varUB, objCoef,
				 rowLower, rowUpper, colType, 1, infinity,
				 rowSense);
    }

    int argc = 1;
    char** argv = new char* [1];
    argv[0] = (char *) "mibs";

#ifdef  COIN_HAS_MPI
    AlpsKnowledgeBrokerMPI brokerPH(argc, argv, *pHModel);
#else
    AlpsKnowledgeBrokerSerial brokerPH(argc, argv, *pHModel);
#endif

    brokerPH.search(pHModel);
    brokerPH.printBestSolution();

    MibSSolution *mibsSol = NULL;

    /*
  if(brokerPH.getSolStatus() == AlpsExitStatusOptimal){
    AlpsSolution *sol = dynamic_cast<AlpsSolution* >
      (brokerPH.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
    mibsSol = dynamic_cast<MibSSolution* >(sol);
  }
  
  if(mibsSol != NULL){
    optSol = new double[uColNum];
    //we know that all UL variables are binary
    for(i = 0; i < uColNum; i++){
      optSol[i] = floor(mibsSol->getValues()[i] + 0.5);
    }
  }
  else if(brokerPH.getSolStatus() == AlpsExitStatusTimeLimit){
    isTimeLimReached = true;
  }
  else if(brokerPH.getNumKnowledges(AlpsKnowledgeTypeSolution) <= 0){
    throw CoinError("Problem PH is infeasible",
    "solvePHProb", "MibSModel");
  }
    */
    if(timeLimit - broker_->timer().getTime() <= 0){
	isTimeLimReached = true;
    }

    if(brokerPH.getNumKnowledges(AlpsKnowledgeTypeSolution) > 0){
	    AlpsSolution *sol = dynamic_cast<AlpsSolution* >
		(brokerPH.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
	    mibsSol = dynamic_cast<MibSSolution* >(sol);
    }

    if(mibsSol != NULL){
	optSol = new double[uColNum];
	//we know that all UL variables are binary
	for(i = 0; i < uColNum; i++){
	    optSol[i] = floor(mibsSol->getValues()[i] + 0.5);
	}
    }




    delete [] objCoef;
    delete [] rowLower;
    delete [] rowUpper;
    delete [] uColInd;
    delete [] uRowInd;
    delete [] structRowInd;
    delete coefMatrix;
    delete coefRowMatrix;
    delete pHModel;

    return optSol;
}


//#############################################################################
// this function is used when user selects saa approach
void
MibSModel::setupSAA(const CoinPackedMatrix& matrix,
		    const CoinPackedMatrix& rowMatrix,
		    const double* varLB, const double* varUB,
		    const double* objCoef, const double* conLB,
		    const double* conUB, const char * colType,
		    double objSense, int truncNumCols, int truncNumRows,
		    double infinity, const char *rowSense)
{

    //saharSto2: we have one assumption: All SAA problems are feasible
    //saharStoSAA: ask about generating samples
    //saharStoSAA: lcm
    //saharStoSAA: ask about defining symphony then check all COIN_HAS_SYMPHONY

    const int clockType = AlpsPar()->entry(AlpsParams::clockType);
    broker_->timer().setClockType(clockType);
    broker_->subTreeTimer().setClockType(clockType);
    broker_->tempTimer().setClockType(clockType);
    broker_->timer().start();
    
    std::string feasCheckSolver =
	MibSPar_->entry(MibSParams::feasCheckSolver);

    int whichCutsLL =
	MibSPar_->entry(MibSParams::whichCutsLL);

    int maxThreadsLL =
	MibSPar_->entry(MibSParams::maxThreadsLL);
    
    double timeLimit(AlpsPar()->entry(AlpsParams::timeLimit));

    int isSMPSFormat(MibSPar_->entry(MibSParams::isSMPSFormat));

    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));
    
    int replNum(MibSPar_->entry(MibSParams::replNumSAA));//M
    int evalSampleSize(MibSPar_->entry(MibSParams::evalSampSizeSAA));//N'
    int i(0), j(0), m(0);
    int index(0), rowNumElem(0), allScenariosNumSMPS(0);
    double element(0.0);
    double etol(etol_);
    double remainingTime(0.0) ,objU(0.0), objL(0.0);
    double objSAA(0.0);//z_N^i
    double objEval(0.0);//\hat{z}_N'(\hat{x}^i)
    double objBestEvalSol(infinity);//\hat{z}_N'(\hat{x}^*)
    double optGap(0.0);//optimality gap estimator
    double estimatedLowerBound(0.0), varLower(0.0), varUpper(0.0);
    bool isLowerInfeasible(false), isTimeLimReached(false);
    bool allEvalsInfeas(true);
    CoinPackedVector appendedRow;
    CoinPackedVector row;
    int *rowInd = NULL;
    double *rowElem = NULL;
    double *b2Base = NULL;
    //double **fullMatrixA2 = NULL;
    std::vector<std::vector<double> > fullMatrixA2;
    CoinPackedMatrix *matrixA2Base = NULL;
    int truncLColNum(getTruncLowerDim());
    int truncLRowNum(getTruncLowerRowNum());
    int uCols(truncNumCols - truncLColNum);
    int uRows(truncNumRows - truncLRowNum);

    //count the number of all scenarios in smps format
    if(isSMPSFormat == PARAM_ON){
	allScenariosNumSMPS = countScenariosSMPS();
	if(isA2Random != PARAM_OFF){
	    /*fullMatrixA2 = new double*[truncLRowNum];
	    for(i = 0; i < truncLRowNum; i++){
		//saharSto3: check if they are filled with zeros
		fullMatrixA2[i] = new double[uCols]();
		}*/
	    fullMatrixA2.resize(truncLRowNum, std::vector<double>(uCols, 0.0));
	}
	else{
	    matrixA2Base = new CoinPackedMatrix(false, 0, 0);
	    matrixA2Base->setDimensions(0, uCols);
	}
	for(i = uRows; i < truncNumRows; i++){
	    row = rowMatrix.getVector(i);
	    rowInd = row.getIndices();
	    rowElem = row.getElements();
	    rowNumElem = row.getNumElements();
	    for(j = 0; j < rowNumElem; j++){
		index = rowInd[j];
		element = rowElem[j];
		if(index < uCols){
		    appendedRow.insert(index, element);
		    if(isA2Random != PARAM_OFF){
			fullMatrixA2[i - uRows][index] = element;
		    }
		}
	    }
	    if(isA2Random == PARAM_OFF){
		matrixA2Base->appendRow(appendedRow);
	    }
	    appendedRow.clear();
	}

	b2Base = new double[truncLRowNum];
	for(i = 0; i < truncLRowNum; i++){
	    index = i + uRows;
	    if(rowSense[index] == 'L'){
		b2Base[i] = conUB[index];
	    }
	    else{
		b2Base[i] = conLB[index];
	    }
	}
    }
    
    //generating N' samples for evaluation
    int evalLRowNum(truncLRowNum * evalSampleSize);
    double *evalRHS = new double[evalLRowNum];
    CoinZeroN(evalRHS, evalLRowNum);
    CoinPackedMatrix *evalA2Matrix = NULL;
    evalA2Matrix = generateSamples(evalSampleSize, truncNumCols,
				   truncNumRows, allScenariosNumSMPS, rowSense,
				   b2Base, matrixA2Base, fullMatrixA2, evalRHS);
    
    //evaluated upper-level solutions
    std::map<std::vector<double>, double> seenULSolutions;
    OsiSolverInterface * evalLSolver = 0;
    OsiSolverInterface * evalBestSolver = 0;
    //store G2 matrix to avoid extracting it for evaluatuion
    int lcmDenum;
    if(isSMPSFormat == PARAM_ON){
	lcmDenum = 1;
    }
    else{
	int incB2Denum(MibSPar_->entry(MibSParams::incDistB2DenumSAA));
	int incA2Denum(MibSPar_->entry(MibSParams::incDistA2DenumSAA));
	int gcdDenum = greatestCommonDivisor(CoinMin(incB2Denum, incA2Denum), CoinMax(incB2Denum, incA2Denum));
	lcmDenum = (incB2Denum * incA2Denum)/gcdDenum;
    }
    
    CoinPackedMatrix *matrixG2 = new CoinPackedMatrix(false, 0, 0);
    matrixG2->setDimensions(0, truncLColNum);
    for(i = uRows; i < truncNumRows; i++){
	row = rowMatrix.getVector(i);
	rowInd = row.getIndices();
	rowElem = row.getElements();
	rowNumElem = row.getNumElements();
	for(j = 0; j < rowNumElem; j++){
	    index = rowInd[j];
	    if(index >= uCols){
		appendedRow.insert(index - uCols, rowElem[j] * lcmDenum);
	    }
	}
	matrixG2->appendRow(appendedRow);
	appendedRow.clear();
    }

    double *optSolRepl  = NULL;
    std::vector<double> optSolReplVec(uCols);
    double *bestEvalSol = new double[uCols];//\hat{x*}
    CoinZeroN(bestEvalSol, uCols);
    //stores z_N^i for i = 1, ..., m 
    double *objValSAARepls = new double[replNum];
    CoinZeroN(objValSAARepls, replNum);
    double *tmpArr = new double[evalSampleSize + 1];
    CoinZeroN(tmpArr, evalSampleSize + 1);
    //stores Q(\hat{x}, \xi(\omega^n)) for n = 1, ..., N' and c\hat{x*}
    for(m = 0; m < replNum; m++){
	//std::cout << "Replication: " << m + 1 << std::endl;
	optSolRepl = solveSAA(matrix, rowMatrix, varLB, varUB, objCoef, conLB,
			      conUB, colType, lcmDenum, objSense, truncNumCols, truncNumRows,
			      infinity, rowSense, isTimeLimReached, objSAA,
			      allScenariosNumSMPS, b2Base, matrixA2Base, fullMatrixA2);

	if(isTimeLimReached == true){
	    goto TERM_SETUPSAA;
	}

	std::cout << "Optimal solution of replication " << m + 1 << " = " << objSAA << std::endl; 

	if(optSolRepl != NULL){
	    objValSAARepls[m] = objSAA;
	    std::copy(optSolRepl, optSolRepl + uCols, optSolReplVec.begin());
	    if(m == 1){
		memcpy(bestEvalSol, optSolRepl, sizeof(double) * uCols);
	    }
	    if(seenULSolutions.find(optSolReplVec) ==
	       seenULSolutions.end()){
		objU = 0.0;
		objEval = 0.0;
		for(i = 0; i < evalSampleSize; i++){
		    objL = 0.0;
		    isLowerInfeasible = false;
		    evalLSolver = setUpEvalModels(matrixG2, optSolRepl, evalRHS,
						  evalA2Matrix, varLB, varUB, objCoef, rowSense,
						  colType, infinity, i, uCols, uRows);

		    remainingTime = timeLimit - broker_->timer().getTime();
		    remainingTime = CoinMax(remainingTime, 0.00);
		    if(remainingTime <= etol){
			isTimeLimReached = true;
			goto TERM_SETUPSAA;
		    }

		    if (feasCheckSolver == "Cbc"){
			dynamic_cast<OsiCbcSolverInterface *>
			    (evalLSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
		    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
			sym_environment *env = dynamic_cast<OsiSymSolverInterface*>
			    (evalLSolver)->getSymphonyEnvironment();

			sym_set_dbl_param(env, "time_limit", remainingTime);
			sym_set_int_param(env, "do_primal_heuristic", FALSE);
			sym_set_int_param(env, "verbosity", -2);
			sym_set_int_param(env, "prep_level", -1);
			sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
			sym_set_int_param(env, "tighten_root_bounds", FALSE);
			sym_set_int_param(env, "max_sp_size", 100);
			sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
			if (whichCutsLL == 0){
			    sym_set_int_param(env, "generate_cgl_cuts", FALSE);
			}else{
			    sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
			}
#endif
		    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
			evalLSolver->setHintParam(OsiDoReducePrint);
			evalLSolver->messageHandler()->setLogLevel(0);
			CPXENVptr cpxEnv =
			    dynamic_cast<OsiCpxSolverInterface*>(evalLSolver)->getEnvironmentPtr();
			assert(cpxEnv);
			CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
			CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
		    }

		    evalLSolver->branchAndBound();


		    if(feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
			if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface*>
						     (evalLSolver)->getSymphonyEnvironment())){
			    isTimeLimReached = true;
			    goto TERM_SETUPSAA;
			}
#endif
		    }
		       
		    if(!evalLSolver->isProvenOptimal()){
			isLowerInfeasible = true;
			break;
		    }
		    else{
			evalBestSolver = setUpEvalModels(matrixG2, optSolRepl, evalRHS,
							 evalA2Matrix, varLB, varUB, objCoef,
							 rowSense, colType, infinity, i, uCols, uRows,
							 evalLSolver->getObjValue(), false);

			remainingTime = timeLimit - broker_->timer().getTime();
			remainingTime = CoinMax(remainingTime, 0.00);
			if(remainingTime <= etol){
			    isTimeLimReached = true;
			    goto TERM_SETUPSAA;
			}

			if (feasCheckSolver == "Cbc"){
			    dynamic_cast<OsiCbcSolverInterface *>
				(evalBestSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
			}else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
			    sym_environment *env = dynamic_cast<OsiSymSolverInterface*>
				(evalBestSolver)->getSymphonyEnvironment();

			    sym_set_dbl_param(env, "time_limit", remainingTime);
			    sym_set_int_param(env, "do_primal_heuristic", FALSE);
			    sym_set_int_param(env, "verbosity", -2);
			    sym_set_int_param(env, "prep_level", -1);
			    sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
			    sym_set_int_param(env, "tighten_root_bounds", FALSE);
			    sym_set_int_param(env, "max_sp_size", 100);
			    sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
			    if (whichCutsLL == 0){
				sym_set_int_param(env, "generate_cgl_cuts", FALSE);
			    }else{
				sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
			    }
#endif
			}else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
			    evalBestSolver->setHintParam(OsiDoReducePrint);
			    evalBestSolver->messageHandler()->setLogLevel(0);
			    CPXENVptr cpxEnv =
				dynamic_cast<OsiCpxSolverInterface*>(evalBestSolver)->getEnvironmentPtr();
			    assert(cpxEnv);
			    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
			    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
			}

			evalBestSolver->branchAndBound();

			if(feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
			    if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface*>
							 (evalBestSolver)->getSymphonyEnvironment())){
				isTimeLimReached = true;
				goto TERM_SETUPSAA;
			    }
#endif
			}
			
			if(!evalBestSolver->isProvenOptimal()){
			    throw CoinError("evalBestSolver cannot be infeasible",
					    "setupSAA",
					    "MibSModel");
			}
			else{
			    objL = evalBestSolver->getObjValue();
			}
			//const double * values = evalLSolver->getColSolution();
	                //for(j = 0; j < truncLColNum; j++){
			//  objL += objCoef[j + uCols] * values[j];
			//}
			tmpArr[i] = objL/evalSampleSize;
		    }
		}
		if(isLowerInfeasible == true){
		    seenULSolutions[optSolReplVec] = infinity;
		}
		else{
		    for(j = 0; j < uCols; j++){
			objU += objCoef[j] * optSolRepl[j];
		    }
		    tmpArr[evalSampleSize] = objU;
		    allEvalsInfeas = false;
		    for(j = 0; j < evalSampleSize + 1; j++){
			objEval += tmpArr[j]; 
		    }
		    if(objBestEvalSol - objEval > etol){
			objBestEvalSol = objEval;
			memcpy(bestEvalSol, optSolRepl, sizeof(double) * uCols);
		    }
		    seenULSolutions[optSolReplVec] = objEval;
		}
	    }
	    //optSolReplVec.clear();
	    delete [] optSolRepl;
	}
	else{
	    throw CoinError("SAA problem is infeasible",
			    "setupSAA",
			    "MibSModel");
	}
    }//for m

 TERM_SETUPSAA:
    if(evalLSolver){
	delete evalLSolver;
	    }
    if(evalBestSolver){
	delete evalBestSolver;
    }
    delete evalA2Matrix;
    delete [] evalRHS;

    if(isTimeLimReached == true){
	std::cout << "Time limit is reached" << std::endl;
    }
    //else if(allEvalsInfeas == true){
	//std::cout << "All solutions of replicates are infeasible ";
	//std::cout << "with respect to evaluation sample." << std::endl;
	//}
    else{
	double estimatedObj(0.0);
	//we generate a new eval sample to be unbiased
	double *evalRHSNew = new double[evalLRowNum];
	CoinZeroN(evalRHSNew, evalLRowNum);
	CoinPackedMatrix *evalA2MatrixNew = NULL;
	evalA2MatrixNew = generateSamples(evalSampleSize, truncNumCols,
					  truncNumRows, allScenariosNumSMPS, rowSense,
					  b2Base, matrixA2Base, fullMatrixA2, evalRHSNew);

	OsiSolverInterface * evalLSolverNew = 0;
	OsiSolverInterface * evalBestSolverNew = 0;

	CoinZeroN(tmpArr, evalSampleSize + 1);
	objU = 0.0;
	for(i = 0; i < evalSampleSize; i++){
	    objL = 0.0;
	    isLowerInfeasible = false;
	    evalLSolverNew = setUpEvalModels(matrixG2, bestEvalSol, evalRHSNew,
					     evalA2MatrixNew, varLB, varUB, objCoef,
					     rowSense, colType, infinity, i, uCols, uRows);
	    remainingTime = timeLimit - broker_->timer().getTime();
	    remainingTime = CoinMax(remainingTime, 0.00);
	    if(remainingTime <= etol){
		isTimeLimReached = true;
		break;
	    }

	    if (feasCheckSolver == "Cbc"){
		dynamic_cast<OsiCbcSolverInterface *>
		    (evalLSolverNew)->getModelPtr()->messageHandler()->setLogLevel(0);
	    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
		sym_environment *env = dynamic_cast<OsiSymSolverInterface*>
		    (evalLSolverNew)->getSymphonyEnvironment();

		sym_set_dbl_param(env, "time_limit", remainingTime);
		sym_set_int_param(env, "do_primal_heuristic", FALSE);
		sym_set_int_param(env, "verbosity", -2);
		sym_set_int_param(env, "prep_level", -1);
		sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
		sym_set_int_param(env, "tighten_root_bounds", FALSE);
		sym_set_int_param(env, "max_sp_size", 100);
		sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
		if (whichCutsLL == 0){
		    sym_set_int_param(env, "generate_cgl_cuts", FALSE);
		}else{
		    sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
		}
#endif
	    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
		evalLSolverNew->setHintParam(OsiDoReducePrint);
		evalLSolverNew->messageHandler()->setLogLevel(0);
		CPXENVptr cpxEnv =
		    dynamic_cast<OsiCpxSolverInterface*>(evalLSolverNew)->getEnvironmentPtr();
		assert(cpxEnv);
		CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
		CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
	    }

	    evalLSolverNew->branchAndBound();

	    if(feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
		if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface*>
					     (evalLSolverNew)->getSymphonyEnvironment())){
		    isTimeLimReached = true;
		    break;
		}
#endif
	    }
	    if(!evalLSolverNew->isProvenOptimal()){
		isLowerInfeasible = true;
		break;
	    }
	    else{
		evalBestSolverNew = setUpEvalModels(matrixG2, bestEvalSol, evalRHSNew,
						    evalA2MatrixNew, varLB, varUB, objCoef,
						    rowSense, colType, infinity, i, uCols,
						    uRows, evalLSolverNew->getObjValue(), false);

		remainingTime = timeLimit - broker_->timer().getTime();
		remainingTime = CoinMax(remainingTime, 0.00);
		if(remainingTime <= etol){
		    isTimeLimReached = true;
		    break;
		}

		if (feasCheckSolver == "Cbc"){
		    dynamic_cast<OsiCbcSolverInterface *>
			(evalBestSolverNew)->getModelPtr()->messageHandler()->setLogLevel(0);
		}else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
		    sym_environment *env = dynamic_cast<OsiSymSolverInterface*>
			(evalBestSolverNew)->getSymphonyEnvironment();

		    sym_set_dbl_param(env, "time_limit", remainingTime);
		    sym_set_int_param(env, "do_primal_heuristic", FALSE);
		    sym_set_int_param(env, "verbosity", -2);
		    sym_set_int_param(env, "prep_level", -1);
		    sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
		    sym_set_int_param(env, "tighten_root_bounds", FALSE);
		    sym_set_int_param(env, "max_sp_size", 100);
		    sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
		    if (whichCutsLL == 0){
			sym_set_int_param(env, "generate_cgl_cuts", FALSE);
		    }else{
			sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
		    }
#endif
		}else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
		    evalBestSolverNew->setHintParam(OsiDoReducePrint);
		    evalBestSolverNew->messageHandler()->setLogLevel(0);
		    CPXENVptr cpxEnv =
			dynamic_cast<OsiCpxSolverInterface*>(evalBestSolverNew)->getEnvironmentPtr();
		    assert(cpxEnv);
		    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
		    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
		}

		evalBestSolverNew->branchAndBound();

		if(feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
		    if(sym_is_time_limit_reached(dynamic_cast<OsiSymSolverInterface*>
						 (evalBestSolverNew)->getSymphonyEnvironment())){
			isTimeLimReached = true;
			break;
		    }
#endif
		}

		if(!evalBestSolverNew->isProvenOptimal()){
		    throw CoinError("evlaBestSolverNew cannot be infeasible",
				    "setupSAA",
				    "MibSModel");
		}
		else{
		    objL = evalBestSolverNew->getObjValue();
		}
						
		//const double * values = evalLSolverNew->getColSolution();
		//for(j = 0; j < truncLColNum; j++){
		//  objL += objCoef[j + uCols] * values[j];
		//}
		tmpArr[i] = objL;
	    }
	}

	if(isLowerInfeasible == true){
	    std::cout << "Optimality gap is infinite :(" << std::endl;
	}
	else if(isTimeLimReached == true){
	    std::cout << "Time limit is reached" << std::endl;
	}   
	else{
	    for(j = 0; j < uCols; j++){
		objU += objCoef[j] * bestEvalSol[j];
	    }
	    for(j = 0; j < evalSampleSize; j++){
		estimatedObj += tmpArr[j];
	    }
	    estimatedObj = estimatedObj/evalSampleSize;
	    estimatedObj += objU;
	    findOptGapVarSAA(objValSAARepls, estimatedObj,
			     estimatedLowerBound, varLower, varUpper, tmpArr);
	    printSolutionSAA(truncNumCols, estimatedObj, estimatedLowerBound,
			     varLower, varUpper, bestEvalSol);

	    delete evalLSolverNew;
            if(evalBestSolverNew){
		delete evalBestSolverNew;
	    }
	    delete evalA2MatrixNew;
	    delete [] evalRHSNew;
	}
    }
    if(b2Base != NULL){
	delete [] b2Base;
    }
    if(matrixA2Base != NULL){
	delete matrixA2Base;
    }
    delete matrixG2;
    delete [] tmpArr;
    delete [] bestEvalSol;
    delete [] objValSAARepls;

    //saharSto2:exit :(
    exit(0);    

}

//#############################################################################
int
MibSModel::countScenariosSMPS()
{
    int numScenarios(0);

    std::string key;
    
    std::string stoFileName = getStoFile();
    fileCoinReadable(stoFileName);

    std::ifstream stoData_stream(stoFileName.c_str());

    if(!stoData_stream){
	std::cout << "Error opening input stochastic data file. Aborting.\n";
	abort();
    }

    while(stoData_stream >> key){
	if(key == "SC"){
	    numScenarios ++;
	}
    }

    return numScenarios;

}

//############################################################################# 
void
MibSModel::printSolutionSAA(int truncNumCols, double estimatedObj,
			    double estimatedLowerBound, double varLower,
			    double varUpper, double *bestSol)
{

    int i(0);
    int truncLColNum(getTruncLowerDim());
    int uCols(truncNumCols - truncLColNum);
    double nearInt(0.0);
    int sampleSize(MibSPar_->entry(MibSParams::sampSizeSAA));//N
    int evalSampleSize(MibSPar_->entry(MibSParams::evalSampSizeSAA));//N'
    int replNum(MibSPar_->entry(MibSParams::replNumSAA));//M
    
    std::cout <<  "============================================" << std::endl;

    std::cout << "Number of replications : " << replNum << std::endl;
    std::cout << "Number of samples in each replication : " << sampleSize << std::endl;
    std::cout << "Number of samples in evalauation : " << evalSampleSize << std::endl;
    std::cout << "Estimated cost (upper bound) = " << estimatedObj << std::endl;
    for(i = 0; i < uCols; i++){
	if(bestSol[i] > 1.0e-15 || bestSol[i] < -1.0e-15){
	    nearInt = floor(bestSol[i] + 0.5);
	    if(ALPS_FABS(nearInt - bestSol[i]) < 1.0e-6){
		std::cout << "x[" << i << "] = " << nearInt << std::endl;
	    }
	    else{
		std::cout << "x[" << i << "] = " << bestSol[i] << std::endl;
	    }		    
	}
    }
    std::cout << "Estimated lower bound = " << estimatedLowerBound << std::endl;
    std::cout << "Estimated optimality gap = " << estimatedObj - estimatedLowerBound << std::endl;
    std::cout << "Estimated standard deviation of upper bound = " << sqrt(varUpper) << std::endl;
    std::cout << "Estimated standard deviation of lower bound = " << sqrt(varLower) << std::endl;
    std::cout << "Estimated standard deviation of the ptimality gap = ";
    std::cout << sqrt(varUpper + varLower) << std::endl;
    std::cout << "Search CPU time: " << broker_->timer().getCpuTime();
    std::cout << " seconds" << std::endl;
    std::cout << "Search wall-clock time: " << broker_->timer().getWallClock();
    std::cout << " seconds" << std::endl;
    
}

//#############################################################################
void
MibSModel::findOptGapVarSAA(double *objValSAARepls, double estimatedObj,
			    double &estimatedLowerBound, double&varLower,
			    double&varUpper, double *estimatedObjComps)
{
    
    int evalSampleSize(MibSPar_->entry(MibSParams::evalSampSizeSAA));//N'
    int replNum(MibSPar_->entry(MibSParams::replNumSAA));//M 

    int i(0);
    double tmpVal(0.0);
    int truncLRowNum(getTruncLowerRowNum());

    for(i = 0; i < replNum; i++){
	estimatedLowerBound += objValSAARepls[i];
    }
    estimatedLowerBound = estimatedLowerBound/replNum;

    for(i = 0; i < replNum; i++){
	varLower += pow(objValSAARepls[i] - estimatedLowerBound, 2);
    }
    varLower = varLower/((replNum - 1) * replNum);

    tmpVal = estimatedObjComps[evalSampleSize];

    for(i = 0; i < evalSampleSize; i++){
	varUpper += pow(tmpVal + estimatedObjComps[i] -
			   estimatedObj, 2);
    }
    varUpper = varUpper/((evalSampleSize - 1) * evalSampleSize);


}

//#############################################################################
OsiSolverInterface *
MibSModel::setUpEvalModels(CoinPackedMatrix *matrixG2, double *optSol,
			   double *allRHS, CoinPackedMatrix *allA2Matrix,
			   const double *origColLb, const double *origColUb,
			   const double* uObjCoef, const char *rowSense,
			   const char *colType, double infinity, int scenarioIndex,
			   int uCols, int uRows, double optLowerObj,
			   bool isLowerProblem)
{

    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));
    
    std::string feasCheckSolver =
	MibSPar_->entry(MibSParams::feasCheckSolver);

    OsiSolverInterface *nSolver = 0;

    int i(0);
    int cnt(0);
    int lCols = getTruncLowerDim();
    int lRows = getTruncLowerRowNum();
    int numCols = uCols + lCols;
    int numRows(0);
    double objSense(getLowerObjSense());
    double *lObjCoeffs = getLowerObjCoeffs();

    if(isLowerProblem == true){
	numRows = lRows;
    }
    else{
	numRows = lRows + 1;
    }

    //setting col bounds
    double *colUb = new double[lCols];
    double *colLb = new double[lCols];
    CoinDisjointCopyN(origColUb + uCols, lCols, colUb);
    CoinDisjointCopyN(origColLb + uCols, lCols, colLb);
    

    //setting coefficient matrix
    CoinPackedMatrix *coefMatrix = new CoinPackedMatrix(false, 0, 0);
    coefMatrix->setDimensions(0, lCols);
    for(i = 0; i < lRows; i++){
	coefMatrix->appendRow(matrixG2->getVector(i));
    }
    if(isLowerProblem == false){
	CoinPackedVector row;
	for(i = 0; i < lCols; i++){
	    row.insert(i, lObjCoeffs[i] * objSense);
	}
	coefMatrix->appendRow(row);
    }

    //setting row bounds
    double *rowUb = new double[numRows];
    double *rowLb = new double[numRows];
    CoinFillN(rowUb, numRows, infinity);
    CoinFillN(rowLb, numRows, -1 * infinity);
    if(isA2Random != PARAM_OFF){
	cnt = lRows * scenarioIndex;
    }
    CoinPackedMatrix *matrixA2 = new CoinPackedMatrix(false, 0, 0);
    matrixA2->setDimensions(0, uCols);
    for(i = 0; i < lRows; i++){
	matrixA2->appendRow(allA2Matrix->getVector(cnt + i));
    }
    double *multA2XOpt = new double[lRows];
    matrixA2->times(optSol, multA2XOpt);
    for(i = 0; i < lRows; i++){
	if(rowSense[i + uRows] == 'L'){
	    rowUb[i] = allRHS[lRows * scenarioIndex + i] - multA2XOpt[i];
	}
	if(rowSense[i + uRows] == 'G'){
	    rowLb[i] = allRHS[lRows * scenarioIndex + i] - multA2XOpt[i];
	}
    }
    if(isLowerProblem == false){
	rowUb[lRows] = optLowerObj;
    }

    double *objCoeffs = new double[lCols];
    if(isLowerProblem == true){
	memcpy(objCoeffs, lObjCoeffs, sizeof(double) * lCols);
    }
    else{
	CoinDisjointCopyN(uObjCoef + uCols, lCols, objCoeffs);
    }
    

    if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"setUpEvalModels", "MibSModel");
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"setUpEvalModels", "MibSModel");
#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"setUpModel", "MibsBilevel");
    }

    nSolver->loadProblem(*coefMatrix, colLb, colUb,
			 objCoeffs, rowLb, rowUb);

    for(i = uCols; i < numCols; i++){
	if((colType[i] == 'B') || (colType[i] == 'I')){
	    nSolver->setInteger(i - uCols);
	}
    }

    nSolver->setObjSense(objSense); 

    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    delete [] colLb;
    delete [] colUb;
    delete [] objCoeffs;
    delete [] rowLb;
    delete [] rowUb;
    delete [] multA2XOpt;
    delete matrixA2;
    delete coefMatrix;

    return nSolver;

}

//#############################################################################
double *
MibSModel::solveSAA(const CoinPackedMatrix& matrix,
		    const CoinPackedMatrix& rowMatrix,
		    const double* varLB, const double* varUB,
		    const double* objCoef, const double* conLB,
		    const double* conUB, const char * colType, int lcmDenum,
		    double objSense, int truncNumCols, int truncNumRows,
		    double infinity, const char *rowSense, bool &isTimeLimReached,
		    double &objSAA, int allScenariosNumSMPS, double *b2Base,
		    CoinPackedMatrix *matrixA2Base,
		    const std::vector<std::vector<double> > &fullMatrixA2)
{
    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));

    std::string feasCheckSolver
	= MibSPar_->entry(MibSParams::feasCheckSolver);

    double timeLimit(AlpsPar()->entry(AlpsParams::timeLimit));

    int clockType(AlpsPar()->entry(AlpsParams::clockType));

    bool useUBDecompose(MibSPar_->entry(MibSParams::useUBDecompose));

    double etol(etol_);
    int sampleSize(MibSPar_->entry(MibSParams::sampSizeSAA));//N
    int i(0), j(0);
    int index(0), pos(0);
    double optVal(0.0);
    int truncLColNum(getTruncLowerDim());
    int truncLRowNum(getTruncLowerRowNum());
    int uColNum(truncNumCols - truncLColNum);
    int uRowNum(truncNumRows - truncLRowNum);
    int lColNum(sampleSize * truncLColNum);
    int lRowNum(sampleSize * truncLRowNum);
    int numRows(uRowNum + lRowNum);
    
    double remainingTime = timeLimit - broker_->timer().getTime();
    remainingTime = CoinMax(remainingTime, 0.00);

    if(remainingTime < etol){
	isTimeLimReached = true;
	return NULL;
    }
    
    /** Create new MibS model to solve bilevel **/
    MibSModel *modelSAA = new MibSModel();

    /** Set up lp solver **/
    OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);
    modelSAA->setSolver(&lpSolver);
    //set parameters
    //modelSAA->AlpsPar()->setEntry(AlpsParams::msgLevel, -1);
    //modelSAA->AlpsPar()->setEntry(AlpsParams::nodeLimit, boundCutNodeLim);
    modelSAA->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
    modelSAA->BlisPar()->setEntry(BlisParams::heurStrategy, 0);
    modelSAA->MibSPar()->setEntry(MibSParams::feasCheckSolver, feasCheckSolver.c_str());
    modelSAA->MibSPar()->setEntry(MibSParams::printProblemInfo, false);
    modelSAA->MibSPar()->setEntry(MibSParams::useBoundCut, false);
    //saharSto2: think about it
    modelSAA->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyLinking);
    //modelSAA->MibSPar()->setEntry(MibSParams::branchStrategy, MibSBranchingStrategyFractional);
    modelSAA->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
    modelSAA->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useTypeWatermelon, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_ON);
    modelSAA->MibSPar()->setEntry(MibSParams::useTypeTenderIC, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useTypeHybridIC, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::useNewPureIntCut, false);
    modelSAA->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt, PARAM_ON);
    modelSAA->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXVarsInt, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsInt, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed, PARAM_ON);
    modelSAA->MibSPar()->setEntry(MibSParams::computeBestUBWhenXVarsInt, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsInt, PARAM_OFF);
    modelSAA->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed, PARAM_ON);
    modelSAA->MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);

    modelSAA->isInterdict_ = false;
    modelSAA->MibSPar()->setEntry(MibSParams::stochasticityType, "stochasticWithSAA");
    modelSAA->MibSPar()->setEntry(MibSParams::isA2Random, isA2Random);
    modelSAA->MibSPar()->setEntry(MibSParams::useUBDecompose, useUBDecompose);
    modelSAA->AlpsPar()->setEntry(AlpsParams::clockType, clockType);
    

    modelSAA->numScenarios_ = sampleSize;

    double *b2Arr = new double[lRowNum];
    CoinZeroN(b2Arr, lRowNum);
    CoinPackedMatrix *stocA2Matrix = NULL;
    stocA2Matrix = generateSamples(sampleSize, truncNumCols,
				   truncNumRows, allScenariosNumSMPS, rowSense,
				   b2Base, matrixA2Base, fullMatrixA2, b2Arr);

    modelSAA->setTruncLowerDim(truncLColNum);
    modelSAA->setTruncLowerRowNum(truncLRowNum);

    modelSAA->lowerObjCoeffs_ = new double[truncLColNum];
    memcpy(modelSAA->lowerObjCoeffs_, getLowerObjCoeffs(),
	   sizeof(double) * truncLColNum);

    modelSAA->setLowerObjSense(getLowerObjSense());

    modelSAA->lowerColInd_ = new int[truncLColNum];
    modelSAA->lowerRowInd_ = new int[truncLRowNum];
    CoinIotaN(modelSAA->lowerColInd_, truncLColNum, uColNum);
    CoinIotaN(modelSAA->lowerRowInd_, truncLRowNum, uRowNum);

    modelSAA->setStocA2Matrix(stocA2Matrix);

    modelSAA->setLowerDim(lColNum);
    modelSAA->setLowerRowNum(lRowNum);
    
    modelSAA->origRowLb_ = new double[numRows];
    modelSAA->origRowUb_ = new double[numRows];
    CoinFillN(modelSAA->origRowLb_, numRows, -1 * infinity);
    CoinFillN(modelSAA->origRowUb_, numRows, infinity);

    for(i = 0; i < sampleSize; i++){
	pos = i * truncLRowNum;
	for(j = uRowNum; j < truncNumRows; j++){
	    index = pos + j;
	    if(rowSense[j] == 'L'){
		modelSAA->origRowUb_[index] = b2Arr[index - uRowNum];
	    }
	    if(rowSense[j] == 'G'){
		modelSAA->origRowLb_[index] = b2Arr[index - uRowNum];
	    }
	}
    }

    modelSAA->loadProblemData(matrix, rowMatrix, varLB, varUB, objCoef,
			      conLB, conUB, colType, objSense, infinity, rowSense,
			      lcmDenum);


    int argc = 1;
    char** argv = new char* [1];
    argv[0] = (char *) "mibs";

#ifdef  COIN_HAS_MPI
    AlpsKnowledgeBrokerMPI brokerSAA(argc, argv, *modelSAA);
#else
    AlpsKnowledgeBrokerSerial brokerSAA(argc, argv, *modelSAA);
#endif

    brokerSAA.search(modelSAA);

    brokerSAA.printBestSolution();

    MibSSolution *mibsSol = NULL;

    if(brokerSAA.getSolStatus() == AlpsExitStatusOptimal){
	AlpsSolution *sol = dynamic_cast<AlpsSolution* >
	    (brokerSAA.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
	objSAA = brokerSAA.getBestQuality();
	objSAA = objSAA/sampleSize;
	mibsSol = dynamic_cast<MibSSolution* >(sol);
    }

    double *optSol = NULL;
    if(mibsSol != NULL){
	optSol = new double[uColNum];
	for(i = 0; i < uColNum; i++){
	    optVal = mibsSol->getValues()[i];
	    if(modelSAA->solver()->isInteger(i)){
		optSol[i] = floor(optVal + 0.5);
	    }
	    else{
		optSol[i] = optVal;
	    }
	}
    }
    else{
	if(brokerSAA.getSolStatus() == AlpsExitStatusTimeLimit){
	    isTimeLimReached = true;
	}
	else if(brokerSAA.getNumKnowledges(AlpsKnowledgeTypeSolution) <= 0){
	    throw CoinError("SAA problem is infeasible",
			    "solveSAA",
			    "MibSModel");
	}
    }

    delete [] b2Arr;
    delete modelSAA;

    return optSol;
	    
}

//#############################################################################
CoinPackedMatrix *
MibSModel::generateSamples(int size, int truncNumCols, int truncNumRows,
			   int allScenariosNum, const char *rowSense,
			   double *b2Base, CoinPackedMatrix *matrixA2Base,
			   const std::vector<std::vector<double> > &fullMatrixA2,
			   double *rhs)
{
    int isSMPSFormat(MibSPar_->entry(MibSParams::isSMPSFormat));
    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));
    
    int i(0), j(0), n(0), p(0);
    int numGeneratedSamples(0), count(0), posRow(0), posCol(0), index(0);
    int genVal(0), scCount(0);
    bool generatedBefore(false), isFinished(false);
    double tmpVal(0.0), dValue(0.0), element(0.0);
    std::string key;
    CoinPackedVector appendedRow;
    double etol(etol_);
    int truncLowerDim(getTruncLowerDim());
    int truncLowerRowNum(getTruncLowerRowNum());
    int uColNum(truncNumCols - truncLowerDim);
    int uRowNum(truncNumRows - truncLowerRowNum);
    
    CoinPackedVector row;
    CoinPackedMatrix *stocMatrixA2 = 0;
    stocMatrixA2 = new CoinPackedMatrix(false, 0, 0);
    stocMatrixA2->setDimensions(0, uColNum);

#ifdef COIN_HAS_C11
    static unsigned int lastSeed = 123456;
    lastSeed = 1664525 * lastSeed + 1013904223;
    std::default_random_engine generator;
    generator.seed(lastSeed);
#endif
    
    if(isSMPSFormat == PARAM_ON){
	std::string stoFileName = getStoFile();
	std::ifstream stoData_stream(stoFileName.c_str());
#ifdef COIN_HAS_C11
	std::uniform_int_distribution<int> distScenIndex(1, allScenariosNum);
#endif

	int *scenIndex = new int[size];
	while(numGeneratedSamples < size){
#ifdef COIN_HAS_C11
	    index = distScenIndex(generator);
#else
	    index = 1 + int(CoinDrand48() * allScenariosNum);
#endif
	    generatedBefore = false;
	    for(i = 0; i < numGeneratedSamples; i++){
		if(scenIndex[i] == index){
		    generatedBefore = true;
		    break;
		}
	    }
	    if(generatedBefore == false){
		scenIndex[numGeneratedSamples] = index;
		numGeneratedSamples ++;
	    }
	}

	stoData_stream.close();

	std::vector<std::vector<double> > fullMatrixA2Copy;

	if(isA2Random != PARAM_OFF){
	    fullMatrixA2Copy.resize(truncLowerRowNum, std::vector<double>(uColNum, 0.0));
	    for(i = 0; i < truncLowerRowNum; i++){
		for(j = 0; j < uColNum; j++){
		    fullMatrixA2Copy[i][j] = fullMatrixA2[i][j];
		}
	    }
	}
	else{
	    for(i = 0; i < truncLowerRowNum; i++){
		stocMatrixA2->appendRow(matrixA2Base->getVector(i));
	    }
	}

	for(i = 0; i < size; i++){
	    count = 1;
	    index = i * truncLowerRowNum;
	    CoinDisjointCopyN(rhs + index, truncLowerRowNum, b2Base);
	    stoData_stream.open(stoFileName.c_str());
	    while(stoData_stream >> key){
		isFinished = false;
		if(key == "SC"){
		    if(count == scenIndex[i]){
			stoData_stream >> key;
			stoData_stream >> key;
			stoData_stream >> key;
			stoData_stream >> key;
			while(isFinished == false){
			    stoData_stream >> key;
			    if((key == "SC") || (key == "ENDATA")){
				isFinished = true;
			    }
			    else if(key == "RHS"){
				stoData_stream >> key;
				stoData_stream >> dValue;
				posRow = -1;
				for(p = 0; p < truncNumRows; p++){
				    if(rowName_[p] == key){
					posRow = p;
					break;
				    }
				}
				if(posRow < 0){
				    std::cout << key << " does not belong to the list of constraints." << std::endl;
				    throw CoinError("Wrong constriant name",
						    "generateSamples",
						    "MibSModel");
				}
				else if(posRow < uRowNum){
				    throw CoinError("Upper-level RHS cannot be random",
						    "generateSamples",
						    "MibSModel");
				}
				else{
				    index = i * truncLowerRowNum + posRow - uRowNum;
				    if(rowSense[posRow] == 'L'){
					rhs[index] = dValue;
				    }
				    else{
					rhs[index] = dValue;
				    }
				}
				
			    }
			    else{
				posCol = -1;
				for(p = 0; p < truncNumCols; p++){
				    if(columnName_[p] == key){
					posCol = p;
					break;
				    }
				}
				if(posCol < 0){
				    std::cout << key << " does not belong to the list of variables." << std::endl;
				    throw CoinError("Wrong variable name",
						    "generateSamples",
						    "MibSModel");
				}
				else if(posCol >= uColNum){
				    throw CoinError("Lower-level coefficients cannot be random",
						    "generateSamples",
						    "MibSModel");
				}
				else{
				    stoData_stream >> key;
				    stoData_stream >> dValue;
				    posRow = -1;
				    for(p = 0; p < truncNumRows; p++){
					if(rowName_[p] == key){
					    posRow = p;
					    break;
					}
				    }
				    if(posRow < 0){
					std::cout << key << " does not belong to the list of constraints." << std::endl;
					throw CoinError("Wrong constriant name",
							"generateSamples",
							"MibSModel");
				    }
				    else if(posRow < uRowNum){
					throw CoinError("Coefficients of lower-level constraints cannot be random",
							"generateSamples",
							"MibSModel");
				    }
				    else{
					index = i * truncLowerRowNum + posRow - uRowNum;
					fullMatrixA2Copy[index][posCol] = dValue;
				    }
				    
				}
			    }
			}
		    }
		    else{
			count ++;
		    }
		}
		if(isFinished == true){
		    stoData_stream.close();
		    break;
		}
	    }
	    if(isA2Random != PARAM_OFF){
		for(p = 0; p < truncLowerRowNum; p++){
		    for(j = 0; j < uColNum; j++){
			element = fullMatrixA2Copy[p][j];
			if(fabs(element) > etol){
			    appendedRow.insert(j, element);
			}
		    }
		    stocMatrixA2->appendRow(appendedRow);
		    appendedRow.clear();
		}
	    }
	}
	delete [] scenIndex;
    }
    else{
	int intRandVal;
	int lbB2(MibSPar_->entry(MibSParams::lbDistB2SAA));
	int ubB2(MibSPar_->entry(MibSParams::ubDistB2SAA));
	int lbA2(MibSPar_->entry(MibSParams::lbDistA2SAA));
	int ubA2(MibSPar_->entry(MibSParams::ubDistA2SAA));
        //saharSto2: it is assumed that incB2Numer <= incB2Denum
	//and they are relatively prime (the same for A2)
	//strict assumption: numerators are 1
	int incB2Numer(MibSPar_->entry(MibSParams::incDistB2NumerSAA));
        int incB2Denum(MibSPar_->entry(MibSParams::incDistB2DenumSAA));
        int incA2Numer(MibSPar_->entry(MibSParams::incDistA2NumerSAA));
        int incA2Denum(MibSPar_->entry(MibSParams::incDistA2DenumSAA));

	int gcdDenum = greatestCommonDivisor(CoinMin(incB2Denum, incA2Denum), CoinMax(incB2Denum, incA2Denum));
        int lcmDenum = (incB2Denum * incA2Denum)/gcdDenum;

        int tmpB2 = (ubB2 - lbB2) * incB2Denum/incB2Numer + 1;
        int tmpA2 = (ubA2 - lbA2) * incA2Denum/incA2Numer + 1;

#ifdef COIN_HAS_C11
        std::uniform_int_distribution<int> distB2(1, tmpB2);
        std::uniform_int_distribution<int> distA2(1, tmpA2);
#endif
    
        for(n = 0; n < size; n++){
	    //generate A2
	    for(i = 0; i < truncLowerRowNum; i++){
		for(j = 0; j < uColNum; j++){
#ifdef COIN_HAS_C11
		    tmpVal = (distA2(generator) - 1) * incA2Numer * lcmDenum/incA2Denum + lbA2 * lcmDenum;
#else
		    intRandVal = 1 + int(CoinDrand48() * tmpA2);
		    tmpVal = (intRandVal - 1) * incA2Numer * lcmDenum/incA2Denum + lbA2 * lcmDenum;
#endif
		    if(fabs(tmpVal) > etol_){
			genVal = (int)round(tmpVal);
		        row.insert(j, genVal);
		    }
		}
		stocMatrixA2->appendRow(row);
	        row.clear();
#ifdef COIN_HAS_C11
	        tmpVal = (distB2(generator) - 1) * incB2Numer * lcmDenum/incB2Denum + lbB2 * lcmDenum;
#else
		intRandVal = 1 + int(CoinDrand48() * tmpB2);
		tmpVal = (intRandVal - 1) * incB2Numer * lcmDenum/incB2Denum + lbB2 * lcmDenum;
#endif
	        genVal = (int)round(tmpVal);
	        rhs[n * truncLowerRowNum + i] = genVal;
	    }
	}
    }

    return stocMatrixA2;
}

//#############################################################################
int
MibSModel::greatestCommonDivisor(int number1, int number2)
{
    if(number1 == 0){
	return number2;
    }
    return greatestCommonDivisor(number2 % number1, number1);
}

//#############################################################################
void
MibSModel::loadProblemData(const CoinPackedMatrix& matrix,
			   const CoinPackedMatrix& rowMatrix,
			   const double* colLB, const double* colUB,   
			   const double* obj, const double* rowLB,
			   const double* rowUB, const char *types,
			   double objSense, double infinity, const char *rowSense,
			   int lcmDenum)
{
   //FIXME: THIS ISN'T TRUE IF WE LOAD AN INTERDICTION PROBLEM 
   //AS A "GENERAL" PROBLEM.  FOR NOW, IT'S OK SINCE WE ONLY
   //DO THIS FROM KNAP SOLVER, WHICH SHOULD SET THIS ITSELF.

    int numRows(0), numCols(0);
    int numScenarios(getNumScenarios());

    std::string stochasticityType(MibSPar_->entry
				  (MibSParams::stochasticityType));

    int isA2Random(MibSPar_->entry(MibSParams::isA2Random));

    int problemType(MibSPar_->entry(MibSParams::bilevelProblemType));

    int isSMPSFormat(MibSPar_->entry(MibSParams::isSMPSFormat));

    if(stochasticityType == "deterministic"){
	numScenarios_ = 1;
    }

    if(stochasticityType == "deterministic"){

        //int i(0);

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

	numRows = matrix.getNumRows();
	numCols = matrix.getNumCols();
    }
  
   int i(0), j(0);
   int beg(0);

   double *varLB(NULL);  
   double *varUB(NULL);  
   double *conLB(NULL);  
   double *conUB(NULL);  
   double *objCoef(NULL);
   char   *colType(NULL);

   CoinPackedMatrix *newMatrix = NULL;

   //switch (problemType){
      
       //case 0:

   if((stochasticityType != "deterministic") ||
      (problemType == 0) || (interdictCost_ == NULL)){

      if(stochasticityType == "deterministic"){
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
	  int k(0);
	  int colBeg(0), rowBeg(0), index(0);
	  double prob(0.0);
	  std::vector<double> scenarioProb(getScenarioProb());
	  CoinPackedMatrix *stocA2Matrix(getStocA2Matrix());
	  int numLCols(getLowerDim());
	  int numLRows(getLowerRowNum());
	  int truncNumLCols(getTruncLowerDim());
	  int truncNumLRows(getTruncLowerRowNum());
	  int truncNumCols(matrix.getNumCols());
	  int truncNumRows(matrix.getNumRows());
	  int numUCols = truncNumCols - truncNumLCols;
	  int numURows = truncNumRows - truncNumLRows;
	  numCols = numUCols + numLCols;
	  numRows = numURows + numLRows;
	  
	  CoinPackedVector row1;
	  CoinPackedVector row2;
	  int *rowInd = NULL;
	  double *rowElem = NULL;
	  int rowNumElem(0);

	  if (!structRowInd_){
	      structRowInd_ = new int[numRows];
	      CoinIotaN(structRowInd_, numRows, 0);
	      structRowNum_ = numRows;
	  }

	  varLB = new double [numCols];
	  varUB = new double [numCols];
	  conLB = new double [numRows];
	  conUB = new double [numRows];
	  objCoef = new double [numCols];
	  colType = new char [numCols];

	  //set matrix
	  //extract truncG2
	  //int lcmDenum (1);
	  //if(stochasticityType == "stochasticWithSAA"){
	  //  lcmDenum = 1;
	  //  }
	   
	  CoinPackedMatrix *truncMatrixG2 = NULL;
	  truncMatrixG2 = new CoinPackedMatrix(true, 0, 0);
	  truncMatrixG2->setDimensions(truncNumLRows, 0);
	  for(i = numUCols; i < truncNumCols; i++){
	      row1 = matrix.getVector(i);
	      rowInd = row1.getIndices();
	      rowElem = row1.getElements();
	      rowNumElem = row1.getNumElements();
	      for(j = 0; j < rowNumElem; j++){
		  row2.insert(rowInd[j] - numURows, rowElem[j] * lcmDenum);
	      }
	      truncMatrixG2->appendCol(row2);
	      row1.clear();
	      row2.clear();
	  }
	  truncMatrixG2->reverseOrdering();

	  //matrix.reverseOrdering();

	  newMatrix = new CoinPackedMatrix(false, 0, 0);
	  newMatrix->setDimensions(0, numCols);
	  for(i = 0; i < numURows; i++){
	      newMatrix->appendRow(rowMatrix.getVector(i));
	  }

	  for(i = 0; i < numScenarios; i++){
	      for(j = 0; j < truncNumLRows; j++){
		  if(isA2Random != PARAM_OFF){
		      rowBeg = i * truncNumLRows;
		  }
		  row1 = truncMatrixG2->getVector(j);
		  rowInd = row1.getIndices();
		  rowElem = row1.getElements();
		  rowNumElem = row1.getNumElements();
		  row2 = stocA2Matrix->getVector(rowBeg + j);
		  for(k = 0; k < rowNumElem; k++){
		      row2.insert(numUCols + i * truncNumLCols +
				  rowInd[k], rowElem[k]);
		  }
		  newMatrix->appendRow(row2);
	      }
	      row1.clear();
	      row2.clear();
	  }

	  newMatrix->reverseOrdering();

	  delete truncMatrixG2;

	  //set row and col bounds and col type
	  colBeg = truncNumCols;
	  CoinDisjointCopyN(colLB, truncNumCols, varLB);
	  CoinDisjointCopyN(colUB, truncNumCols, varUB);
	  CoinDisjointCopyN(types, truncNumCols, colType);
	  for(i = 0; i < numScenarios - 1; i++){
	      CoinDisjointCopyN(colLB + numUCols, truncNumLCols,
				varLB + colBeg);
	      CoinDisjointCopyN(colUB + numUCols, truncNumLCols,
				varUB + colBeg);
	      CoinDisjointCopyN(types + numUCols, truncNumLCols,
				colType + colBeg);

	      colBeg += truncNumLCols;
	  }

	  //saharSto: can we avoid defining conLB and conUB?
	  CoinDisjointCopyN(rowLB, numURows, origRowLb_);
	  CoinDisjointCopyN(rowUB, numURows, origRowUb_);
	  memcpy(conLB, origRowLb_, sizeof(double) * numRows);
	  memcpy(conUB, origRowUb_, sizeof(double) * numRows);

	  //set objective
	  index = numUCols;
	  CoinDisjointCopyN(obj, numUCols, objCoef);
	  if(stochasticityType == "stochasticWithoutSAA"){
	      double *subObjCoef = new double[truncNumLCols];
	      CoinDisjointCopyN(obj + numUCols, truncNumLCols, subObjCoef);
	      for(i = 0; i < numScenarios; i++){
		  prob = scenarioProb[i];
		  for(j = 0; j < truncNumLCols; j++){
		      objCoef[index] = prob * subObjCoef[j];
		      index ++;
		  }
	      }

	      delete [] subObjCoef;
	  }
	  else{
	      for(i = 0; i < numUCols; i++){
		  objCoef[i] = objCoef[i] * numScenarios;
	      }
	      for(i = 0; i < numScenarios; i++){
		  index = numUCols + i * truncNumLCols;
		  CoinDisjointCopyN(obj + numUCols, truncNumLCols, objCoef + index);
	      }
	  }

	  //saharSto: check it
	  origRowSense_ = new char[truncNumRows];
	  memcpy(origRowSense_, rowSense, truncNumRows);
      }
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
   if((stochasticityType == "deterministic") || (stochasticityType == "stochasticWithoutSAA")){
       instanceStructure(newMatrix, conLB, conUB, rowSense);
   }
}

//#############################################################################
void 
MibSModel::setUpperColData()
{

    std::string stochasticityType(MibSPar_->entry
				  (MibSParams::stochasticityType));
    
   int lowerColNum(lowerDim_);
   
   if (!upperDim_){
      upperDim_ = numVars_ - lowerDim_;

      if(!getUpperColInd())
	  upperColInd_ = new int[upperDim_];

      if(stochasticityType == "deterministic"){
	  int * lowerColInd = getLowerColInd();
      
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
      else{
	  CoinIotaN(upperColInd_, upperDim_, 0);
      }
   }
   //numOrigVars_ = lowerDim_ + upperDim_;
}

//#############################################################################
void 
MibSModel::setUpperRowData()
{

   //FIXME: MAKE THIS MORE SIMPLE

    std::string stochasticityType(MibSPar_->entry
				  (MibSParams::stochasticityType));

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

   int i(0), cnt(0), beg(0), end(0);

   if(stochasticityType == "deterministic"){
       for(i = 0; i < lowerRowNum_ + upperRowNum_; i++){
	   if(!findIndex(i, lowerRowNum, lowerRowInd)){
	       upperRowInd_[cnt] = i;
	       cnt++;
	   }
       }
       assert(cnt == upperRowNum_);
   }
   else{
       beg = lowerRowInd[0];
       end = beg + lowerRowNum_;

       CoinIotaN(upperRowInd_, beg, 0);
       CoinIotaN(upperRowInd_ + beg, numCons_ - end, end);
   }

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
   
   std::string dataFile(AlpsPar_->entry(AlpsParams::instance));

   if((numCols_ == 0) && (dataFile == "NONE")){
       std::cout << "Error: data file is not specefied. Aborting." << std::endl;
       abort();
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

  int numScenarios(getNumScenarios());
  
  if(0)
    solver()->writeLp("userfeasible1");

  int i(0), index(0);
  double upperObj(0.0);
  bool isHeurSolution(true);
  int * upperColInd = getUpperColInd();
  int * lowerColInd = getLowerColInd();
  double * lpSolution = NULL;
  //double * lpSolution = new double[getNumCols()];
  //CoinFillN(lpSolution, getNumCols(), 0.0);
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
      lpSolution = new double[getNumCols()];
      CoinFillN(lpSolution, getNumCols(), 0.0);
      for(i = 0; i < upperDim_; i++){
	  index = upperColInd[i];
	  lpSolution[index] = bS_->optUpperSolutionOrd_[i];
	  upperObj +=
	      bS_->optUpperSolutionOrd_[i] * solver()->getObjCoefficients()[index];
      }
      for(i = 0; i < lowerDim_; i++){
	  if(numScenarios == 1){
	      index = lowerColInd[i];
	  }
	  else{
	      index = upperDim_ + i;
	  }
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
				upperObj, this);
  }
  else if(solType == MibSHeurSol){
      //in stochastic case, we assume that lower-level variables
      //do not participate in the upper-level constraints.
      if((!bS_->isUBSolved_) && (numScenarios == 1)){
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
  if(lpSolution){
      delete [] lpSolution;
  }
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

  double lhs(0.0), value(0.0);
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
      if((getLColLbInLProb() != NULL) || (getLColUbInLProb() != NULL)){
	  int lCols(getLowerDim());
	  int * lColIndices(getLowerColInd());
	  const double * origColLb(getOrigColLb());
	  const double * origColUb(getOrigColUb());
	  for(i = 0; i < lCols; i++){
	      index1 = lColIndices[i];
	      value = solution[index1];
	      if(((origColLb[index1] - value) > etol_)
		 || ((value - origColUb[index1]) > etol_)){
		  upperFeasible = false;
		  break;
	      }
	  }
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

  int numScenarios(getNumScenarios());
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
      if(numScenarios == 1){
	  index = lIndices[i];
      }
      else{
	  index = udim + i;
      }
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
  std::string stochasticityType(MibSPar_->entry
				(MibSParams::stochasticityType));
  double * varlower = varLB();
  double * varupper = varUB();

  if(stochasticityType == "deterministic"){
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
  else{
      int truncNumCols = getUpperDim() + getTruncLowerDim();
      origColLb_ = new double[truncNumCols];
      origColUb_ = new double[truncNumCols];
      CoinDisjointCopyN(varlower, truncNumCols, origColLb_);
      CoinDisjointCopyN(varupper, truncNumCols, origColUb_);
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
    //saharSto: we should change the structure of fixedInd
    std::string stochasticityType(MibSPar_->entry
				  (MibSParams::stochasticityType));

    int uCols(upperDim_);
    int lRows(lowerRowNum_);
    int * upperColInd = getUpperColInd();
    int * lowerRowInd = getLowerRowInd();

    const double * matElements = newMatrix->getElements();
    const int * matIndices = newMatrix->getIndices();
    const int * matStarts = newMatrix->getVectorStarts();

    int index1, rowIndex, posRow, start, end, begPos, endPos;
    int i, j;
    int num(0);

    if(!fixedInd_){
	fixedInd_ = new int[numVars_]();
    }

    if(stochasticityType == "deterministic"){
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
    else{
	begPos = lowerRowInd[0];
	endPos = begPos + lRows - 1;
	for(i = 0; i < numVars_; i++){
	    fixedInd_[i] = 0;
	    if(i < uCols){
		start = matStarts[i];
		end = start + newMatrix->getVectorSize(i);
		for(j = start; j < end; j++){
		    rowIndex = matIndices[j];
		    if((begPos <= rowIndex) && (rowIndex <= endPos)){
			fixedInd_[i] = 1;
			sizeFixedInd_ ++;
			break;
		    }
		}
	    }
	}
    }

}

//#############################################################################
void                                                                                                                                                       
MibSModel::instanceStructure(const CoinPackedMatrix *newMatrix,
                             const double* rowLB, const double* rowUB,
                             const char *rowSense)
{

   bool printProblemInfo(MibSPar_->entry(MibSParams::printProblemInfo));
    
   std::string stochasticityType(MibSPar_->entry
				(MibSParams::stochasticityType));

   int paramValue(0);

   //set the default value of time limit if it is not set by user
   //Alps sets the default value of this parameter to a very large value
   //but it is very large and creats memory issues for setting the time limit
   //when an extra mip should be solved. So, we set it to a smaller large value.
   if (AlpsPar()->entry(AlpsParams::timeLimit) > 1.0e200){
       AlpsPar()->setEntry(AlpsParams::timeLimit, 1.0e200);
   }

   if(stochasticityType == "stochasticWithoutSAA"){
       BlisPar()->setEntry(BlisParams::heurStrategy, 0);

       //Param: "MibS_usePreprocessor"
       paramValue = MibSPar_->entry(MibSParams::usePreprocessor);

       if (paramValue == PARAM_NOTSET){
	   MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
       }else if (paramValue == PARAM_ON){
	   std::cout << "The preprocessor is not currently functional.";
	   std::cout << std::endl;
	   MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
       }

       //Param: MibS heuristics
       if((MibSPar_->entry(MibSParams::useLowerObjHeuristic) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useObjCutHeuristic) == PARAM_ON)||
	  (MibSPar_->entry(MibSParams::useWSHeuristic) == PARAM_ON)||
	  (MibSPar_->entry(MibSParams::useGreedyHeuristic) == PARAM_ON)){
	   std::cout << "The heuristics are  not currently functional.";
	   std::cout << std::endl;
       }
       MibSPar()->setEntry(MibSParams::useLowerObjHeuristic, PARAM_OFF);
       MibSPar()->setEntry(MibSParams::useObjCutHeuristic, PARAM_OFF);
       MibSPar()->setEntry(MibSParams::useWSHeuristic, PARAM_OFF);
       MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);


       if((MibSPar_->entry(MibSParams::useBoundCut) == true) ||
	  (MibSPar_->entry(MibSParams::useBendersCut) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useNoGoodCut) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useGeneralNoGoodCut) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useTypeIC) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useTypeWatermelon) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useTypeTenderIC) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useTypeHybridIC) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::useIncObjCut) == PARAM_ON) ||
	  (MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_ON)){
	   throw CoinError("Only hypercube intersection cut is currently functional for the stochastic problems",
			   "instanceStructure", "MibSModel");
       }

       if((MibSPar_->entry(MibSParams::useBendersCut) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useNoGoodCut) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useGeneralNoGoodCut) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useTypeIC) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useTypeWatermelon) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useTypeTenderIC) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useTypeHybridIC) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::useIncObjCut) == PARAM_NOTSET) &&
	  (MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_NOTSET)){
	   MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_ON);
	   if (printProblemInfo == true){
	       std::cout <<
		   "Default cut is hypercube intersection cut.";
	       std::cout << std::endl;
	   }
       }
       else if(MibSPar_->entry(MibSParams::useTypeHypercubeIC) == PARAM_ON){
	   std::cout << "Hypercub intersection cut generator is on." << std::endl;
       }
   
       //Param: "MibS_branchProcedure"
       MibSBranchingStrategy branchParSto = static_cast<MibSBranchingStrategy>
	   (MibSPar_->entry(MibSParams::branchStrategy));

       if (branchParSto == MibSBranchingStrategyNotSet){
	   MibSPar()->setEntry(MibSParams::branchStrategy,
			       MibSBranchingStrategyLinking);
	   if (printProblemInfo == true){
	       std::cout <<
		   "Default branching strategy is MibSBranchingStrategyLinking.";
	       std::cout << std::endl;
	   }
       }else if (branchParSto == MibSBranchingStrategyLinking){
	   std::cout << "Branching strategy is MibSBranchingStrategyLinking.";
	   std::cout << std::endl;
       }else{
	   std::cout << "Branching procedure is MibSBranchingStrategyFractional.";
           std::cout << std::endl;
       }

       //Setting parameters of solving (VF) and (UB)
       int solveSecondLevelWhenXYVarsIntSto(MibSPar_->entry(MibSParams::solveSecondLevelWhenXYVarsInt));
       int solveSecondLevelWhenXVarsIntSto(MibSPar_->entry(MibSParams::solveSecondLevelWhenXVarsInt));
       int solveSecondLevelWhenLVarsIntSto(MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsInt));
       int solveSecondLevelWhenLVarsFixedSto(MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsFixed));
       int computeBestUBWhenXVarsIntSto(MibSPar_->entry(MibSParams::computeBestUBWhenXVarsInt));
       int computeBestUBWhenLVarsIntSto(MibSPar_->entry(MibSParams::computeBestUBWhenLVarsInt));
       int computeBestUBWhenLVarsFixedSto(MibSPar_->entry(MibSParams::computeBestUBWhenLVarsFixed));
   

       if ((solveSecondLevelWhenXYVarsIntSto == PARAM_NOTSET) &&
	   (solveSecondLevelWhenXVarsIntSto == PARAM_NOTSET) &&
           (solveSecondLevelWhenLVarsIntSto == PARAM_NOTSET) &&
           (solveSecondLevelWhenLVarsFixedSto == PARAM_NOTSET) &&
           (computeBestUBWhenXVarsIntSto == PARAM_NOTSET) &&
           (computeBestUBWhenLVarsIntSto == PARAM_NOTSET) &&
           (computeBestUBWhenLVarsFixedSto == PARAM_NOTSET)){
	   MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt, PARAM_ON);
           MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed, PARAM_ON);
           MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed, PARAM_ON);
       } 

       if (printProblemInfo == true){
	   if (MibSPar_->entry(MibSParams::solveSecondLevelWhenXYVarsInt) == PARAM_ON){
	       std::cout << "solveSecondLevelWhenXYVarsInt is set." << std::endl;
	   }

           if (MibSPar_->entry(MibSParams::solveSecondLevelWhenXVarsInt) == PARAM_ON){
	       std::cout << "solveSecondLevelWhenXVarsInt is set." <<std::endl;
	   }

           if (MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsInt) == PARAM_ON){
	       std::cout << "solveSecondLevelWhenLVarsInt is set." <<std::endl;
	   }

           if (MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsFixed) == PARAM_ON){
	       std::cout << "solveSecondLevelWhenLVarsFixed is set." <<std::endl;
	   }

           if (MibSPar_->entry(MibSParams::computeBestUBWhenXVarsInt) == PARAM_ON){
	       std::cout << "computeBestUBWhenXVarsInt is set." <<std::endl;
	   }

           if (MibSPar_->entry(MibSParams::computeBestUBWhenLVarsInt) == PARAM_ON){
	       std::cout << "computeBestUBWhenLVarsInt is set." <<std::endl;
	   }

           if (MibSPar_->entry(MibSParams::computeBestUBWhenLVarsFixed) == PARAM_ON){
	       std::cout << "computeBestUBWhenLVarsFixed is set." <<std::endl;
	   }
       }

       //Setting "saveSeenLinkingSols" parameter

       if (MibSPar_->entry(MibSParams::useLinkingSolutionPool) == PARAM_NOTSET){
	   MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
       }

       if (printProblemInfo == true){
	   if(MibSPar_->entry(MibSParams::useLinkingSolutionPool) == PARAM_ON){
	       std::cout << "Linking solution pool will be used." << std::endl;
	   }
           else{
	       std::cout << "Linking solution pool will not be used." << std::endl;
	   }
       }
       std::cout << std::endl;

       return;
   }

   
    
   /** Determines the properties of instance **/
   if(printProblemInfo == true){
      std::cout<<"======================================="<<std::endl;
      std::cout<<"Analyzing problem structure            "<<std::endl;
      std::cout<<"======================================="<<std::endl;
      std::cout << std::endl;
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
         if (colType_[index] != 'C'){
            numUpperInt ++;
         }
      }
      if (printProblemInfo == true){
         std::cout <<"Number of integer UL Variables: ";
         std::cout << numUpperInt << std::endl;
      }
      for (i = 0; i < lCols; i++){
         index = lColIndices[i];
         if (colType_[index] != 'C'){
            numLowerInt++;
         }
      }
      if (printProblemInfo == true){
         std::cout <<"Number of integer LL Variables: ";
         std::cout << numLowerInt << std::endl;
      }
   }
   
   
   if (isInterdict_ == false){
      CoinDisjointCopyN(rowSense, numRows, newRowSense);
   }else{
      for (i = 0; i < numRows; i++){
         newRowSense[i] = 'L';
      }
      for (i = 0; i < lRows-lCols; i++){
         newRowSense[i+uRows] = rowSense[i];
      }
   }
   
   //Checks general or interdiction 
   if (isInterdict_ == true){
      if (printProblemInfo == true){
         std::cout << "This instance is an interdiction problem." << std::endl;
      }
   }
   
   //Checks type of variables
   if (isInterdict_ == false){
      for (i = 0; i < numCols; i++){
         if (colType_[i] == 'C'){
            isPureInteger_ = false;
            break;
         }
      }
      
      for (i = 0; i < numCols; i++){
         if (colType_[i] != 'B'){
            if (binarySearch(0, lCols - 1, i, lColIndices) < 0){
               if(fixedInd_[i] == 1){
                  allLinkingBin_ = false;
               }
               allUpperBin_ = false;
            }
            else{
               allLowerBin_ = false;
            }
            if ((!allUpperBin_) && (!allLowerBin_) && (!allLinkingBin_)){
               break;
            }
         }
      }
   }
   
   if (printProblemInfo == true){
      if (isPureInteger_ == true){
         std::cout << "This instance is a pure integer problem." << std::endl;
      }
      else{
         std::cout << "This instance is a mixed integer problem." << std::endl;
      }
      
      if (allUpperBin_ == true){
         std::cout << "All upper level variables are binary." << std::endl;
      }
      
      if (allLowerBin_ == true){
         std::cout << "All lower level variables are binary." << std::endl; 
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
   for (i = 0; i < numCols; i++){
      counterStart = matStarts[i];
      counterEnd = matStarts[i+1];
      for (j = counterStart; j < counterEnd; j++){                          
         rowIndex = matIndices[j];
         if(newRowSense[rowIndex] == 'L'){
            mult = 1;
         }
         else{
            mult = -1;
         }
         posRow = binarySearch(0, lRows - 1, rowIndex, lRowIndices);       
         posCol = binarySearch(0, lCols - 1, i, lColIndices);
         if ((fabs(matElements[j] - floor(matElements[j])) > etol_) &&
            (fabs(matElements[j] - ceil(matElements[j])) > etol_)){
            if (posRow < 0){
               isUpperCoeffInt_ = false;
            }else{
               isLowerCoeffInt_ = false;
            }
         }
         if (mult * matElements[j] < -etol_){
            if (posRow < 0){
               if (posCol < 0){
                  positiveA1_ = false;
               }  
               else{
                  positiveG1_ = false;
               }  
            }else{  
               if (posCol < 0){
                  positiveA2_ = false;
               }else{
                  positiveG2_ = false;
               }
            }
         }
      }
    }

    if ((isUpperCoeffInt_ == true) || (isLowerCoeffInt_ == true)){
       for (i = 0; i < numRows; i++){
          switch(newRowSense[i]){
           case 'L':
             rhs = rowUB[i];
             break;
           case 'G':
             rhs = rowLB[i];
             break;
           case 'E':
             std::cout << "MibS cannot currently handle equality constraints.";
             std::cout << std::endl; 
             abort();
             break;
           case 'R':
             std::cout << "MibS cannot currently handle range constraints.";
             std::cout << std::endl;
             abort();
             break;
          }
          if ((fabs(rhs - floor(rhs)) > etol_) &&
              (fabs(rhs - ceil(rhs)) > etol_)){
             posRow = binarySearch(0, lRows - 1, rowIndex, lRowIndices);
             if (posRow < 0){
                isUpperCoeffInt_ = false;
             }else{
                isLowerCoeffInt_ = false;
             }
          }
       }
    }
    
    for (i = 0; i < numCols; i++){
       if (fixedInd_[i] == 1){
          if (colType_[i] == 'C'){
             std::cout << "All linking variables should be discrete";
             std::cout << std::endl;
             i = -1;
             assert(i > 0);
          }
       }
    }

    if (printProblemInfo == true){
       if (positiveA1_ == true){
          std::cout << "Coefficient matrix of upper level variables in upper level problem is non-negative." << std::endl;
       }
       
       if (positiveG1_ == true){
          std::cout << "Coefficient matrix of lower level variables in upper level problem is non-negative." << std::endl;
       }
       
       if (positiveA2_ == true){
          std::cout << "Coefficient matrix of upper level variables in lower level problem is non-negative." << std::endl;
       }
       
       if (positiveG2_ == true){
          std::cout << "Coefficient matrix of lower level variables in upper level problem is non-negative." << std::endl;
       }
    }

    std::cout << std::endl;
    std::cout << "=======================================" <<std::endl;
    std::cout << "Setting parameters                     " <<std::endl;
    std::cout << "=======================================" <<std::endl;
    std::cout << std::endl;
   
   
    //MibSIntersectionCutType cutType(MibSIntersectionCutTypeNotSet);
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
    
    if (paramValue == PARAM_NOTSET){
       MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
    }else if (paramValue == PARAM_ON){
       std::cout << "The preprocessor is not currently functional.";
       std::cout << std::endl;
       MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
    }
    
    //Param: "MibS_useGreedyHeuristic"
    paramValue = MibSPar_->entry(MibSParams::useGreedyHeuristic);
    
    if (paramValue == PARAM_NOTSET){
       MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);
    }else if (paramValue == PARAM_ON){
       if (isInterdict_ == false){
          MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);
          std::cout << "Greedy heuristic only works for interdiction problems.";
          std::cout << std::endl; 
       }
    }
    
    //Param: "MibS_useIncObjCut"
    if ((turnOffOtherCuts == true) &&
        (MibSPar_->entry(MibSParams::useIncObjCut) == PARAM_NOTSET)){
       MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useIncObjCut);
    
    if (paramValue == PARAM_NOTSET){
	MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
    }else if (paramValue == PARAM_ON){
       if ((allLinkingBin_ == false) || (positiveA2_ == false)){
          std::cout << "The increasing objective cut is not valid for this problem.";
          std::cout << std::endl;
          MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
       }
    }
    
    //Param: "MibS_useInterSectionCut" (MibSIntersectionCutTypeIC)
    paramValue = MibSPar_->entry(MibSParams::useTypeIC);
    //cutType = static_cast<MibSIntersectionCutType>
    //	(MibSPar_->entry(MibSParams::intersectionCutType));
    
    if (paramValue == PARAM_NOTSET){
       MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
    }else if (paramValue == PARAM_ON){
       if ((isPureInteger_ == false) || (isLowerCoeffInt_ == false)){
          std::cout << "The intersection cut IC is not valid for this problem.";
          std::cout << std::endl;
          MibSPar()->setEntry(MibSParams::useTypeIC, PARAM_OFF);
       }
    }
    /*if((paramValue == PARAM_ON) && (cutType == MibSIntersectionCutTypeIC)){
      if((isPureInteger_ == false) || (isLowerCoeffInt_ == false)){
	    std::cout << "The intersection cut IC is not valid problem.";
	    std::cout << std::endl;
	    MibSPar()->setEntry(MibSParams::useIntersectionCut, PARAM_NOTSET);
	}
    }*/    

    //Param: "MibS_useNoGoodCut"
    paramValue = MibSPar_->entry(MibSParams::useNoGoodCut);

    if (paramValue == PARAM_NOTSET){
	MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    } else if ((paramValue == PARAM_ON) && (allUpperBin_ == false)){
       std::cout << "The no-good cut is not valid for this problem.";
       std::cout << std::endl;
       MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
    }

    //Param: "MibS_useBendersCut"
    if ((turnOffOtherCuts == true) &&
        (MibSPar_->entry(MibSParams::useBendersCut) == PARAM_NOTSET)){
       MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useBendersCut);

    if (paramValue == PARAM_NOTSET){
	if ((isInterdict_ == false) || (positiveG2_ == false)){
	    MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	}else{
	    MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_ON);
	    MibSPar()->setEntry(MibSParams::bendersCutType,
                                MibSBendersCutTypeJustOneCut);
	}
    }
    else if (paramValue == PARAM_ON){
	if ((isInterdict_ == false) || (positiveG2_ == false)){
           std::cout << "The Benders cut is not valid for this problem.";
           std::cout << std::endl;
           MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	}
    }
    if (MibSPar_->entry(MibSParams::useBendersCut) == PARAM_ON){
       defaultCutIsOn = true;
    }

    //Param: "MibS_useGeneralNoGoodCut"
    if ((turnOffOtherCuts == true) &&
        (MibSPar_->entry(MibSParams::useNoGoodCut) == PARAM_NOTSET)){
       MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useGeneralNoGoodCut);

    if (paramValue == PARAM_NOTSET){
       if (allLinkingBin_ == false){
          MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
       }else if (defaultCutIsOn == false){
          MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_ON);
       }
    }else if ((paramValue == PARAM_ON) && (allLinkingBin_ == false)){
       std::cout << "Generalized no-good cut is not valid for this problem.";
       std::cout << std::endl;
       MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
    }
    if (MibSPar_->entry(MibSParams::useGeneralNoGoodCut) == PARAM_ON){
       defaultCutIsOn = true;
    }
    
    //Param: "MibS_usePureIntegerCut"
    if ((turnOffOtherCuts == true) &&
        (MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_NOTSET)){
       MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::usePureIntegerCut);

    if (paramValue == PARAM_NOTSET){
       if ((isPureInteger_ == false) || (isUpperCoeffInt_ == false)
           || (isLowerCoeffInt_ == false)){
          MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
       }else if (defaultCutIsOn == false){
          MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_ON);
       }
    }else if((paramValue == PARAM_ON) &&
             ((isPureInteger_ == false) ||
              (isUpperCoeffInt_ == false) ||
              (isLowerCoeffInt_ == false))){
       std::cout << "The pure integer cut is not valid for this problem.";
       std::cout << std::endl;
       MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
    }
    if (MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_ON){
	defaultCutIsOn = true;
    }

    //Param: "MibS_useInterSectionCut" (MibSIntersectionCutTypeHypercubeIC) 
    if ((turnOffOtherCuts == true) &&
        (MibSPar_->entry(MibSParams::useTypeHypercubeIC) == PARAM_NOTSET)){
       MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_OFF);
    }
    
    paramValue = MibSPar_->entry(MibSParams::useTypeHypercubeIC);
    
    if (paramValue == PARAM_NOTSET){
	if(defaultCutIsOn == false){
           MibSPar()->setEntry(MibSParams::useTypeHypercubeIC, PARAM_ON);
           /*MibSPar()->setEntry(MibSParams::intersectionCutType,
             MibSIntersectionCutTypeHypercubeIC);*/
	}
    }
    
    //Param: "MibS_branchProcedure"
    MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
	  (MibSPar_->entry(MibSParams::branchStrategy));
    if (branchPar == MibSBranchingStrategyNotSet){
       if ((isInterdict_ == true) || (numUpperInt <= numLowerInt)){
          MibSPar()->setEntry(MibSParams::branchStrategy,
                              MibSBranchingStrategyLinking);
          if (printProblemInfo == true){
             std::cout <<
                "Default branching strategy is MibSBranchingStrategyLinking.";
             std::cout << std::endl;
          }
       }else{
          MibSPar()->setEntry(MibSParams::branchStrategy,
                              MibSBranchingStrategyFractional);
          if (printProblemInfo == true){
             std::cout <<
                "Default branching strategy is MibSBranchingStrategyFractional.";
             std::cout << std::endl;
          }
       }
    }else if (MibSPar_->entry(MibSParams::branchStrategy) ==
              MibSBranchingStrategyLinking){
       std::cout << "Branching strategy is MibSBranchingStrategyLinking.";
       std::cout << std::endl;
    }else{
       std::cout << "Branching procedure is MibSBranchingStrategyFractional.";
       std::cout << std::endl;
    }

    if (printProblemInfo == true){
	if(MibSPar_->entry(MibSParams::useIncObjCut) == PARAM_ON){
           std::cout << "Increasing objective cut generator is on.";
           std::cout << std::endl;
	}

	if (MibSPar_->entry(MibSParams::useBendersCut) == PARAM_ON){
           if (MibSPar_->entry(MibSParams::bendersCutType) ==
               MibSBendersCutTypeJustOneCut){
              std::cout << "Benders cut generator (just one cut) is on.";
              std::cout << std::endl;
           }else{
              std::cout << "Benders cut generator (multiple cuts) is on.";
              std::cout << std::endl;
           }
	}
        
        if (MibSPar_->entry(MibSParams::usePureIntegerCut) == PARAM_ON){
           std::cout << "Pure integer cut generator is on."<< std::endl;
	}
        
	if (MibSPar_->entry(MibSParams::useNoGoodCut) == PARAM_ON){
           std::cout << "No-good cut generator is on."<< std::endl;
	}

	if (MibSPar_->entry(MibSParams::useGeneralNoGoodCut) == PARAM_ON){
           std::cout << "Generalized no-good cut generator is on."<< std::endl;
	}

        //if(MibSPar_->entry(MibSParams::useIntersectionCut) == PARAM_ON){
	if (MibSPar_->entry(MibSParams::useTypeIC) == PARAM_ON){
           std::cout << "Intersection cut IC generator is on." << std::endl;
	}
	if (MibSPar_->entry(MibSParams::useTypeWatermelon) == PARAM_ON){
           std::cout << "watermelon IC generator is on." << std::endl;
	}
	if (MibSPar_->entry(MibSParams::useTypeHypercubeIC) == PARAM_ON){
           std::cout << "hypercube IC generator is on." << std::endl;
	}
        //}
    }
	
    //Setting parameters of solving (VF) and (UB)
    int solveSecondLevelWhenXYVarsInt(MibSPar_->entry(MibSParams::solveSecondLevelWhenXYVarsInt));
    int solveSecondLevelWhenXVarsInt(MibSPar_->entry(MibSParams::solveSecondLevelWhenXVarsInt));
    int solveSecondLevelWhenLVarsInt(MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsInt));
    int solveSecondLevelWhenLVarsFixed(MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsFixed));
    int computeBestUBWhenXVarsInt(MibSPar_->entry(MibSParams::computeBestUBWhenXVarsInt));
    int computeBestUBWhenLVarsInt(MibSPar_->entry(MibSParams::computeBestUBWhenLVarsInt));
    int computeBestUBWhenLVarsFixed(MibSPar_->entry(MibSParams::computeBestUBWhenLVarsFixed));

    if ((solveSecondLevelWhenXYVarsInt == PARAM_NOTSET) &&
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

    if (printProblemInfo == true){
       if (MibSPar_->entry(MibSParams::solveSecondLevelWhenXYVarsInt) == PARAM_ON){
          std::cout << "solveSecondLevelWhenXYVarsInt is set." << std::endl;
       }
    
       if (MibSPar_->entry(MibSParams::solveSecondLevelWhenXVarsInt) == PARAM_ON){
          std::cout << "solveSecondLevelWhenXVarsInt is set." <<std::endl;
       }

       if (MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsInt) == PARAM_ON){
          std::cout << "solveSecondLevelWhenLVarsInt is set." <<std::endl;
       }
       
       if (MibSPar_->entry(MibSParams::solveSecondLevelWhenLVarsFixed) == PARAM_ON){
          std::cout << "solveSecondLevelWhenLVarsFixed is set." <<std::endl;
       }
       
       if (MibSPar_->entry(MibSParams::computeBestUBWhenXVarsInt) == PARAM_ON){
          std::cout << "computeBestUBWhenXVarsInt is set." <<std::endl;
	}

       if (MibSPar_->entry(MibSParams::computeBestUBWhenLVarsInt) == PARAM_ON){
          std::cout << "computeBestUBWhenLVarsInt is set." <<std::endl;
       }
       
       if (MibSPar_->entry(MibSParams::computeBestUBWhenLVarsFixed) == PARAM_ON){
          std::cout << "computeBestUBWhenLVarsFixed is set." <<std::endl;
       }
    }
    
    //Setting "saveSeenLinkingSols" parameter
    
    if (MibSPar_->entry(MibSParams::useLinkingSolutionPool) == PARAM_NOTSET){
       MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
    }
    
    if (printProblemInfo == true){
       if(MibSPar_->entry(MibSParams::useLinkingSolutionPool) == PARAM_ON){
          std::cout << "Linking solution pool will be used." << std::endl;
       }
       else{
          std::cout << "Linking solution pool will not be used." << std::endl;
	}
    }
    std::cout << std::endl;
    
    delete [] newRowSense;
}
