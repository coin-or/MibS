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

#include "OsiCbcSolverInterface.hpp"
#include "OsiSymSolverInterface.hpp"
#include "symphony.h"

#include "MibSCutGenerator.h"
#include "MibSParams.h"
#include "MibSTreeNode.h"
#include "MibSSolution.h"

#include "BlisConGenerator.h"
#include "BlisConstraint.h"
#include "BlisHelp.h"
#include "BlisVariable.h"

//#############################################################################
MibSCutGenerator::MibSCutGenerator(MibSModel *mibs)
{
  localModel_ = mibs;
  setName("MIBS");
  setModel(localModel_);
  auxCount_ = 0;
  upper_ = 0.0;
  maximalCutCount_ = 0;
}

//#############################################################################
MibSCutGenerator::~MibSCutGenerator()
{

}

//#############################################################################
int
MibSCutGenerator::bilevelFeasCut1(BcpsConstraintPool &conPool)
{

  /** Add bilevel feasibility cut **/

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const CoinPackedMatrix * matrix = solver->getMatrixByCol();
  MibSBilevel *bS = localModel_->bS_;
  int i(0), j(0);
  int numCuts(0);
  double etol(localModel_->etol_);
  const char * rowsense = solver->getRowSense();

  if(!bS->isIntegral_) return 0;

  bS->isIntegral_ = false;

  /** Find binding constraints at current solution **/
  int * binding = getBindingCons();

  if(!binding){
    std::cout << "No binding constraints set." << std::endl;
    abort();
  }

  double upper(getCutUpperBound());
  //double upper(- getCutUpperBound());
  double cutlb(- solver->getInfinity());
  double cutub(upper - 1.0);
  
  //double EPS(0.01);
  //double cutub(upper - EPS);

  std::vector<int> indexList;
  std::vector<double> valsList;
  double * tempVals = new double[numCols];
  CoinFillN(tempVals, numCols, 0.0);
  double value(0.0);  
  int mult(0);

  for(i = 0; i < numCols; i++){
    for(j = 0; j < numRows; j++){
      value = matrix->getCoefficient(j, i);
      switch(rowsense[j]){
      case 'L':
	mult = 1;
	break;
      case 'G':
	mult = -1;
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
      tempVals[i] += binding[j] * value * mult;
    }

    tempVals[i] += binding[numRows + i]; // variable upper bound
    tempVals[i] += - binding[numRows + numCols + i]; // variable lower bound

    if((tempVals[i] > etol) || (tempVals[i] < - etol)){
      indexList.push_back(i);
      valsList.push_back(tempVals[i]);
    }
  }

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);

  delete [] binding;
  delete [] tempVals;

  return numCuts;
  
}

//#############################################################################
int
MibSCutGenerator::bilevelFeasCut2(BcpsConstraintPool &conPool)
{

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  int i(0);
  int numCuts(0);
  double etol(localModel_->etol_);



  double cutlb(- solver->getInfinity());
  double cutub(solver->getInfinity());

  std::vector<int> indexList;
  std::vector<double> valsList;

  /* 
     returns a double with the values [alpha | beta | gamma]
     indexed 0..n1-1, n1..n1+n2-1, n1+n2
  */

  double * cutVals = findDeepestLandPCut_ValFunc();  
  //double * cutVals = findDeepestLandPCut1();  
  double val(0.0);

  /* find the nonzero entries and convert to <= constraint */

  if(0)
     std::cout << "cut coeffs" << std::endl;

  for(i = 0; i < numCols; i++){
    val = cutVals[i];
    if(0)
      std::cout << i << ": " << -val << std::endl;
    if((val > etol) || (val < - etol)){
      // entry is nonzero
      indexList.push_back(i);
      valsList.push_back(- val);
    }
  }

  if(0)
    std::cout << "cutub: " << -cutVals[numCols] << std::endl;

  cutub = - cutVals[numCols]; //gamma

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);

  delete [] cutVals;
  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::feasibilityCuts(BcpsConstraintPool &conPool)
{

  /** Add bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIBS Cuts." << std::endl;

  //bool simpleCutOnly 
  //= localModel_->MibSPar_->entry(MibSParams::simpleCutOnly);

  bool usePureIntegerCut 
    = localModel_->MibSPar_->entry(MibSParams::usePureIntegerCut);

   bool useValFuncCut 
    = localModel_->MibSPar_->entry(MibSParams::useValFuncCut);

  if(usePureIntegerCut && !useValFuncCut){
    return bilevelFeasCut1(conPool) ? true : false;
  }
  else if(!usePureIntegerCut && useValFuncCut){
    return bilevelFeasCut2(conPool) ? true : false;
  }
  else if(usePureIntegerCut && useValFuncCut){
    return (bilevelFeasCut1(conPool) && 
	    bilevelFeasCut2(conPool)) ? true : false;
  }
  else{
    //std::cout << "No MIBS Cuts generated" << std::endl;
    return 0;
  }

}

//#############################################################################
double *
MibSCutGenerator::findDeepestLandPCut_ValFunc()
{

  // VERIFY THAT THIS WORKS WITH > 1 UL CONSTRAINT

  /* Construct the CGLP to find the deepest Lift and Project Cut */

  OsiSolverInterface * cglpSolver = new OsiClpSolverInterface();
  //OsiSolverInterface * cglpSolver = new OsiSymSolverInterface();
  OsiSolverInterface * oSolver = localModel_->solver();

  int numOrigRows = localModel_->numOrigCons_;
  int numLowerRows = localModel_->lowerRowNum_;
  int numUpperRows = numOrigRows - numLowerRows;
  const int numCols = oSolver->getNumCols();
  const int numRows = oSolver->getNumRows();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByCol();
  const double * rightHandSide = oSolver->getRightHandSide();
  MibSBilevel *bS = localModel_->bS_;
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  int * lowerRowInd = localModel_->getLowerRowInd();
  double * lowerObjCoeffs = localModel_->getLowerObjCoeffs();

  int numEntries(numCols + 1);
  double * cutVals = new double[numEntries];
  CoinZeroN(cutVals, numEntries);

  int numULCols = localModel_->getUpperDim();
  int numLLCols = localModel_->getLowerDim();
  double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

  /* the negative and positive slopes of the LL value function */
  double leftSlope = localModel_->leftSlope_;
  double rightSlope = localModel_->rightSlope_;

  /* 
     For each polyhedron:
     Column for each variable in the LR + gamma 
     + 2 columns for each constraint,including the disjunction;
     Note: vars associated with = constraints are split into 2 vars
     Row for each variable in the LR + 2 gamma constraints + normalization
  */

  int numAuxRows(numRows - numOrigRows);
  //int numAuxRows(0);
  int numCGCols(0);
  numCGCols = numCols + 1+ 4 * numOrigRows + 4 + 2 * numAuxRows;
  int numCGRows(2 * numCols + 2 + 1);
  
  double * colUb = new double[numCGCols];
  double * colLb = new double[numCGCols];
  double * rowUb = new double[numCGRows];
  double * rowLb = new double[numCGRows];
  double * objCoeffs = new double[numCGCols];
  
  /* 
     variables are ordered: alpha, beta, gamma, u,v,w,z,etc 
     alpha, beta, and gamma are free, rest nonnegative
  */

  CoinFillN(colLb, numCols + 1, - cglpSolver->getInfinity());
  CoinFillN(colLb + numCols + 1, (numCGCols - (numCols + 1)), 0.0);

  CoinFillN(colUb, numCGCols, cglpSolver->getInfinity());

  CoinFillN(rowLb, numCGRows - 1, - cglpSolver->getInfinity());
  CoinFillN(rowUb, numCGRows - 1, 0.0);

  /* normalization constraint */
  rowLb[numCGRows - 1] = 1.0;
  rowUb[numCGRows - 1] = 1.0;

  int i(0), j(0);

  CoinZeroN(objCoeffs, numCGCols);

  /* alpha */
  for(i = 0; i < numULCols; i++)
    objCoeffs[i] = upperSol[i];

  /* beta */
  for(i = numULCols; i < numULCols + numLLCols; i++)
    objCoeffs[i] = lowerSol[i - numULCols];

  /* gamma */
  objCoeffs[numULCols + numLLCols] = - 1.0;

  CoinPackedMatrix * cglpMat = new CoinPackedMatrix(false, 0, 0);
  cglpMat->setDimensions(0, numCGCols);

  /* Set up the matrix for the CGLP */

  /* first, the rows constraining alpha */

  int rowIndex(0);
  int colIndex(0), insertIndex1(0), insertIndex2(0);
  int cnt(0);
  double value(0.0);  

  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    row1.insert(i, - 1.0); // alpha coefficient
    row2.insert(i, - 1.0); // alpha coefficient

    for(j = 0; j < numOrigRows; j++){

      rowIndex = j;
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex + cnt;
        
      row1.insert(insertIndex1, value);
      row1.insert(insertIndex1 + 1, - value);

      //changed this 8/27 (sd)
      //insertIndex2 = numULCols + numLLCols + 1 + rowIndex + cnt 
      //+ 2 * numUpperRows + 4 * numLowerRows;
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex + cnt 
	+ 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows;
      
      row2.insert(insertIndex2, value);
      row2.insert(insertIndex2 + 1, - value);
      
      cnt++;
    }

    cnt = 0;

    for(j = 0; j < numLowerRows; j++){
      
      //matIndex = lowerRowInd[j];
      rowIndex = lowerRowInd[j];
      
      value = matrix->getCoefficient(rowIndex, colIndex);

      //changed this 8/27 (sd)
      //insertIndex1 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex - numUpperRows;
      insertIndex1 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + j;
      row1.insert(insertIndex1, value);
      
      insertIndex1 += numLowerRows;
      row1.insert(insertIndex1, - leftSlope * value);
      //row1.insert(insertIndex1, leftSlope * value);
     
      //changed this 8/27 
      //insertIndex2 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex
      //+ 2 * numUpperRows + 4 * numLowerRows - numUpperRows;
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + j
	+ 2 * numUpperRows + 4 * numLowerRows + numAuxRows;
      row2.insert(insertIndex2, - value);
      
      insertIndex2 += numLowerRows;
      row2.insert(insertIndex2, - rightSlope * value); 
      
    }
 
    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      //changed this 8/27 (sd)      
      //insertIndex1 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex 
      //- numUpperRows + numLowerRows;
      insertIndex1 = numULCols + numLLCols + 1 + 
      	2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + rowIndex 
	- numOrigRows;
      
      row1.insert(insertIndex1, - value);

      //changed this 8/27 (sd)      
      //insertIndex2 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + rowIndex
      //+ 2 * numUpperRows + 4 * numLowerRows 
      //- numUpperRows + numLowerRows;
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + numLowerRows + 1 
	+ 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows 
	+ rowIndex - numOrigRows; 
      
      row2.insert(insertIndex2, - value);
      	
    }
 
    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }

  /* then, the rows constraining beta */

  for(i = 0; i < numLLCols; i++){
    
    colIndex = lowerColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    insertIndex1 = i + numULCols;
    row1.insert(insertIndex1, - 1.0); // beta coefficient

    insertIndex2 = i + numULCols;
    row2.insert(insertIndex1, - 1.0); // beta coefficient

    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
    
      //changed this 8/27 (sd)
      //insertIndex1 = numULCols + numLLCols + 1 
      //+ 2 * numUpperRows + rowIndex - numUpperRows;
      insertIndex1 = numULCols + numLLCols + 1 
	 + 2 * numUpperRows + j;
      
      row1.insert(insertIndex1, value);
      row1.insert(insertIndex1 + 1, - value);
      
      //changed this 8/27 (sd)
      //insertIndex2 = numULCols + numLLCols + 1 + rowIndex - numUpperRows 
      //+ 2 * numUpperRows + 2 * numUpperRows + 4 * numLowerRows;
      insertIndex2 = numULCols + numLLCols + 1 + j
	+ 4 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows;
      
      row2.insert(insertIndex2, value);
      row2.insert(insertIndex2 + 1, - value);
    }
    
    
    value = lowerObjCoeffs[i];
    
    insertIndex1 = numULCols + numLLCols + 1 + 
      2 * (numUpperRows + numLowerRows) + numLowerRows;
    
    row1.insert(insertIndex1, - value);
    
    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + 
    //4 * (numUpperRows + numLowerRows) + 3 * numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + 
      4 * (numUpperRows + numLowerRows) + 2 * numLowerRows + 1 + numAuxRows;

    row2.insert(insertIndex2, - value);

    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);
    	
      //changed this 8/27 (sd)
      //insertIndex1 = numULCols + numLLCols + 1 + 
      //2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + rowIndex 
      //- numUpperRows;
      insertIndex1 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + rowIndex 
	- numOrigRows;
      
      row1.insert(insertIndex1, - value);
      
      //changed this 8/27 (sd)
      //insertIndex2 = numULCols + numLLCols + 1 + 
      //4 * (numUpperRows + numLowerRows) + 
      //3 * numLowerRows + 1 + rowIndex - numUpperRows;
      insertIndex2 = numULCols + numLLCols + 1 + 
	4 * (numUpperRows + numLowerRows) + 
	2 * numLowerRows + 2 + rowIndex - numOrigRows + numAuxRows;
            
      row2.insert(insertIndex2, - value);
      
    }

    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);
    
  }
  
  
  /* then, the rows constraining gamma */

  CoinPackedVector row1;
  CoinPackedVector row2;

  insertIndex1 = numULCols + numLLCols;   
  row1.insert(insertIndex1, 1.0); // gamma coefficient
  
  insertIndex2 = numULCols + numLLCols;
  row2.insert(insertIndex2, 1.0); // gammaa coefficient
  
  for(j = 0; j < numUpperRows; j++){
    
    value = rightHandSide[j];

    //changed this 8/27 (sd)
    //insertIndex1 = numULCols + numLLCols + 1 + j;
    insertIndex1 = numULCols + numLLCols + 1 + j + cnt;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + j
    //+ 2 * numUpperRows + 4 * numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + j
      + 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows + cnt;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    
    cnt++;
  }

  cnt = 0;
  
  for(j = 0; j < numLowerRows; j++){
    
    value = rightHandSide[numUpperRows + j];
   
    //changed this 8/27 (sd)
    //insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j;
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j + cnt;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
    //+ 2 * numUpperRows + 4 * numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numUpperRows + 3 * numLowerRows + 1 + numAuxRows + cnt;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    cnt++;    
  }

  cnt = 0;

  double mult1(0.0), mult2(0.0);
  
  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];
    
    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
      mult1 += value * upperSol[i];
    }
  }
  
  for(j = 0; j < numLLCols; j++){

    value = lowerObjCoeffs[j];
    mult2 += value * optLowerSol[j];

  }
    
  insertIndex1 = numULCols + numLLCols + 1 + 2 * (numUpperRows + numLowerRows); 
  
  row1.insert(insertIndex1, - mult1);
  
  insertIndex1 += numLowerRows;
  
  //row1.insert(insertIndex1, leftSlope * mult1 - mult2); 
  row1.insert(insertIndex1, leftSlope * mult1 + mult2); 
  
  //changed this 8/27 (sd)  
  //insertIndex2 = numULCols + numLLCols + 1 
  //+ 4 * (numUpperRows + numLowerRows) + 2 * numLowerRows; 
  insertIndex2 = numULCols + numLLCols + 1 
    + 4 * (numUpperRows + numLowerRows) + numLowerRows + 1 + numAuxRows; 
  
  row2.insert(insertIndex2, mult1);
  
  insertIndex2 += numLowerRows;
  
  //row2.insert(insertIndex2, rightSlope * mult1 - mult2); 
  row2.insert(insertIndex2, rightSlope * mult1 + mult2); 
  
  for(j = numOrigRows; j < numRows; j++){
    
    //value = rightHandSide[j];
    value = - rightHandSide[j];
  
    //changed this 8/27 (sd)
    //insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
    //+ 2 * numLowerRows + numUpperRows + numLowerRows; 
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows 
      + 2 * numLowerRows + numUpperRows + 1 + j - numOrigRows; 
    
    row1.insert(insertIndex1, - value);

    //changed this 8/27 (sd)
    //insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
    //+ 2 * numUpperRows + 4 * numLowerRows
    //+ 2 * numLowerRows + numUpperRows + numLowerRows;
    insertIndex2 = numULCols + numLLCols + 1 + 4 * numUpperRows
      + 6 * numLowerRows + 2 + j - numOrigRows + numAuxRows;

    row2.insert(insertIndex2, - value);

  }
  
  cglpMat->appendRow(row1);
  cglpMat->appendRow(row2);

  /* then a normalization constraint */

  CoinPackedVector row;

  for(i = numULCols + numLLCols + 1; i < numCGCols; i++)
    row.insert(i, 1.0);

  cglpMat->appendRow(row);

  cglpSolver->assignProblem(cglpMat, colLb, colUb, objCoeffs, rowLb, rowUb);

  if(0)
    cglpSolver->writeLp("cglp");

  if(0){
    dynamic_cast<OsiClpSolverInterface *> 
      (cglpSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }
  else{
    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("max_active_nodes", 1);
  }
  cglpSolver->initialSolve();
  
  if(cglpSolver->isProvenOptimal()){

    const double * cglpSol = cglpSolver->getColSolution();
    double cutViolation = cglpSolver->getObjValue();

    for(i = 0; i < numULCols; i++){
      cutVals[upperColInd[i]] = cglpSol[i];
      if(0)
	std::cout << upperColInd[i] << ": " << cglpSol[i] << std::endl;
    }

    for(i = 0; i < numLLCols; i++){
      cutVals[lowerColInd[i]] = cglpSol[numULCols + i];
      if(0){
	std::cout << lowerColInd[i] << ": " 
		  << cglpSol[numULCols + i] << std::endl;
      }
    }
    
    cutVals[numULCols + numLLCols] = cglpSol[numULCols + numLLCols];

    std::cout << "Violated inequality found, with violation of " 
	      << cutViolation << "." << std::endl;

  }
  else{

    std::cout << "CGLP failed to find a violated inequality." << std::endl;
    
  }

  delete cglpSolver;

  /* Return the entries of the deepest cut - alpha * x + beta * y >= gamma */
  return cutVals;

}

//#############################################################################
double *
MibSCutGenerator::findDeepestLandPCut_IncObj(double * upperSol, 
					     double * lowerSol,
					     double * optLowerSol)
{

  //NOTE: THIS IS ONLY CASE CASE OF THE CGLP (I.E. A^2 POS AS IN MIPINT)
  //      THIS ASSUMES INEQUALITY CONSTRAINTS
  //      SHOULD ADD A CHECK FOR THAT AND INTRODUCE +,- VARS FOR =-CONSTRAINTS

  //PROBABLY WANT TO ADD AN ARGUMENT OF UPPERSOL
  //THAT WAY CAN CALL THE SAME FUNCTION FOR BOTH THE CURRENT SOL
  //AND THE MAXIMAL VECTORS

  /* Construct the CGLP to find the deepest Lift and Project Cut */

  OsiSolverInterface * cglpSolver = new OsiClpSolverInterface();
  OsiSolverInterface * oSolver = localModel_->solver();

  int numOrigRows = localModel_->numOrigCons_;
  int numLowerRows = localModel_->lowerRowNum_;
  int numUpperRows = numOrigRows - numLowerRows;
  const int numCols = oSolver->getNumCols();
  const int numRows = oSolver->getNumRows();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByCol();
  //const CoinPackedMatrix * matrix = oSolver->getMatrixByRow();
  const double * rightHandSide = oSolver->getRightHandSide();
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  int * lowerRowInd = localModel_->getLowerRowInd();
  int * upperRowInd = localModel_->getUpperRowInd();
  double * lowerObjCoeffs = localModel_->getLowerObjCoeffs();
  double lowerObjSense = localModel_->getLowerObjSense();
  double etol = localModel_->etol_;

  int numEntries(numCols + 1);
  double * cutVals = new double[numEntries];
  CoinZeroN(cutVals, numEntries);

  int numULCols = localModel_->getUpperDim();
  int numLLCols = localModel_->getLowerDim();

  //PROBABLY WANT TO CHANGE THIS TO BE AN ARGUMENT
  //double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  //double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  //double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

  /* 
     Column for each variable in the LR + gamma 
     + 1 column for each inequality constraint,including the disjunction
     + 2 columns for each equality constraint;
     Note: vars associated with = constraints are split into 2 vars
     Row for each variable in the LR + 2 gamma constraints + normalization
  */

  int numCGCols(0);
  int numAuxRows(numRows - numOrigRows);
  numCGCols = numCols + 1 + 2 * numOrigRows + 4 + 2 * numAuxRows;
  int numCGRows(2 * numCols + 2 + 1);
  

  double * colUb = new double[numCGCols];
  double * colLb = new double[numCGCols];
  double * rowUb = new double[numCGRows];
  double * rowLb = new double[numCGRows];
  double * objCoeffs = new double[numCGCols];
  
  /* 
     variables are ordered: alpha, beta, gamma, u,v,w,z,etc 
     alpha, beta, and gamma are free, rest nonnegative
  */

  CoinFillN(colLb, numCols + 1, - cglpSolver->getInfinity());
  CoinFillN(colLb + numCols + 1, (numCGCols - (numCols + 1)), 0.0);

  CoinFillN(colUb, numCGCols, cglpSolver->getInfinity());

  CoinFillN(rowLb, numCGRows - 1, - cglpSolver->getInfinity());
  CoinFillN(rowUb, numCGRows - 1, 0.0);

  /* normalization constraint */
  rowLb[numCGRows - 1] = 1.0;
  rowUb[numCGRows - 1] = 1.0;

  int i(0), j(0);

  CoinZeroN(objCoeffs, numCGCols);

  /* alpha */
  for(i = 0; i < numULCols; i++)
    objCoeffs[i] = upperSol[i];

  /* beta */
  for(i = numULCols; i < numULCols + numLLCols; i++)
    objCoeffs[i] = lowerSol[i - numULCols];

  /* gamma */
  objCoeffs[numULCols + numLLCols] = - 1.0;

  CoinPackedMatrix * cglpMat = new CoinPackedMatrix(false, 0, 0);
  cglpMat->setDimensions(0, numCGCols);

  /* Set up the matrix for the CGLP */

  /* first, the rows constraining alpha */

  int rowIndex(0);
  int colIndex(0), insertIndex1(0), insertIndex2(0);
  double value(0.0);  

  //matrix->dumpMatrix();

  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    row1.insert(i, 1.0); // alpha coefficient
    row2.insert(i, 1.0); // alpha coefficient
    /*
    for(j = 0; j < numOrigRows; j++){

      rowIndex = j;
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex;
      
      row1.insert(insertIndex1, - value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex 
	+ numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
      
    }
    */
    for(j = 0; j < numUpperRows; j++){

      rowIndex = upperRowInd[j];
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex;
      
      row1.insert(insertIndex1, - value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex 
	+ numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
      
    }

    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];
      value = matrix->getCoefficient(rowIndex, colIndex);
      
      insertIndex1 = numULCols + numLLCols + 1 + rowIndex;
      
      row1.insert(insertIndex1, -value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 + rowIndex 
	+ numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, -value);
      //row2.insert(insertIndex2 + 1, value);
      
    }

    //cnt = 0;

    //FIXME: SHOULD BE 1 IF < ETOL!
    //w variable
    if(upperSol[i] > etol)
      value = 0.0;
    else 
      value = 1.0;

    insertIndex1 = numULCols + numLLCols + 1 + 
      numUpperRows + numLowerRows;
    row1.insert(insertIndex1, - value);

    insertIndex2 = numULCols + numLLCols + 1 
      + 2 * (numLowerRows + numUpperRows) + 2 + numAuxRows;

    row2.insert(insertIndex2, value);

    //multipliers associated with previously-added cuts
    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);

      insertIndex1 = numULCols + numLLCols + 1 + 
	numUpperRows + numLowerRows + 2 + (rowIndex - numOrigRows);
      row1.insert(insertIndex1, -value);
      
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows + 2) 
	+ (rowIndex - numOrigRows) + numAuxRows; 
      row2.insert(insertIndex2, -value);
      	
    }
        
    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }


  /* then, the rows constraining beta */

  for(i = 0; i < numLLCols; i++){
    
    colIndex = lowerColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    insertIndex1 = i + numULCols;
    row1.insert(insertIndex1, 1.0); // beta coefficient

    insertIndex2 = i + numULCols;
    row2.insert(insertIndex1, 1.0); // beta coefficient

    for(j = 0; j < numUpperRows; j++){

      rowIndex = upperRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
    
      insertIndex1 = numULCols + numLLCols + 1 
	+ rowIndex;
      
      row1.insert(insertIndex1, - value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 
	+ rowIndex + numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
    }

    for(j = 0; j < numLowerRows; j++){

      rowIndex = lowerRowInd[j];

      value = matrix->getCoefficient(rowIndex, colIndex);
    
      insertIndex1 = numULCols + numLLCols + 1 
	+ rowIndex;
      
      row1.insert(insertIndex1, -value);
      //row1.insert(insertIndex1 + 1, value);
      
      insertIndex2 = numULCols + numLLCols + 1 
	+ rowIndex + numUpperRows + numLowerRows + 2 + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      //row2.insert(insertIndex2 + 1, value);
    }
    
    value = lowerObjCoeffs[i] * lowerObjSense;
    
    insertIndex1 = numULCols + numLLCols + 1 + 
      numUpperRows + numLowerRows + 1;

    row1.insert(insertIndex1, -value);

    //insertIndex2 = numULCols + numLLCols + 1 + 
    //4 * (numUpperRows + numLowerRows) + 3 * numLowerRows;

    //row2.insert(insertIndex2, - value);

    for(j = numOrigRows; j < numRows; j++){

      rowIndex = j;
      
      value = matrix->getCoefficient(rowIndex, colIndex);
	
      insertIndex1 = numULCols + numLLCols + 1 + 
	numUpperRows + numLowerRows + 2 + (rowIndex - numOrigRows);
      
      row1.insert(insertIndex1, - value);
      
      insertIndex2 = numULCols + numLLCols + 1 + 
	2 * (numUpperRows + numLowerRows + 2)
	+ (rowIndex - numOrigRows) + numAuxRows;
      
      row2.insert(insertIndex2, - value);
      
    }
    
    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);
    
  }

  /* then, the rows constraining gamma */

  CoinPackedVector row1;
  CoinPackedVector row2;

  insertIndex1 = numULCols + numLLCols;   
  row1.insert(insertIndex1, - 1.0); // gamma coefficient
  
  insertIndex2 = numULCols + numLLCols;
  row2.insert(insertIndex2, - 1.0); // gammaa coefficient
  
  for(j = 0; j < numUpperRows; j++){
    
    //value = rightHandSide[j];
    value = rightHandSide[upperRowInd[j]];

    insertIndex1 = numULCols + numLLCols + 1 + j;

    row1.insert(insertIndex1, value);
    //row1.insert(insertIndex1 + 1, - value);

    insertIndex2 = numULCols + numLLCols + 1 + j
      + numUpperRows + numLowerRows + 2 + numAuxRows;

    row2.insert(insertIndex2, value);
    //row2.insert(insertIndex2 + 1, - value);
    
  }
  
  //cnt = 0;

  for(j = 0; j < numLowerRows; j++){
    
    //value = rightHandSide[numUpperRows + j];
    value = rightHandSide[lowerRowInd[j]];
   
    insertIndex1 = numULCols + numLLCols + 1 + numUpperRows + j;

    row1.insert(insertIndex1, value);
    //row1.insert(insertIndex1 + 1, - value);

    insertIndex2 = numULCols + numLLCols + 1 
      + numUpperRows + numLowerRows + 2 + numUpperRows
      + numAuxRows + j;

    row2.insert(insertIndex2, value);
    //row2.insert(insertIndex2 + 1, - value);
    
  }

  double mult2(0.0);

  for(j = 0; j < numLLCols; j++){

    value = lowerObjCoeffs[j] * lowerObjSense;
    mult2 += value * optLowerSol[j];

  }

  insertIndex1 = numULCols + numLLCols + 1 + numUpperRows + numLowerRows + 1; 
  
  row1.insert(insertIndex1, mult2); 
  
  insertIndex2 = numULCols + numLLCols + 1 
    + 2 * (numUpperRows + numLowerRows) + 2 + numAuxRows;

  //row2.insert(insertIndex2, value);
  row2.insert(insertIndex2, - 1);

  for(j = numOrigRows; j < numRows; j++){
    
    value = rightHandSide[j];
  
    insertIndex1 = numULCols + numLLCols + 1 
      + numUpperRows + numLowerRows + 2 + (j - numOrigRows);
    
    row1.insert(insertIndex1, value);

    insertIndex2 = numULCols + numLLCols + 1 
      + 2 * (numUpperRows + numLowerRows + 2) 
      + (j - numOrigRows) + numAuxRows;

    row2.insert(insertIndex2, value);

  }

  cglpMat->appendRow(row1);
  cglpMat->appendRow(row2);

  /* then a normalization constraint */

  CoinPackedVector row;

  for(i = numULCols + numLLCols + 1; i < numCGCols; i++)
    row.insert(i, 1.0);

  cglpMat->appendRow(row);

  cglpSolver->assignProblem(cglpMat, colLb, colUb, objCoeffs, rowLb, rowUb);

  cglpSolver->setObjSense(-1); //maximization

  if(0)
    cglpSolver->writeLp("cglp");
  //cglpSolver->writeMps("cglp");

  if(0){
    dynamic_cast<OsiClpSolverInterface *> 
      (cglpSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }    
  else{
    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("max_active_nodes", 1);
  }

  cglpSolver->initialSolve();

  if(cglpSolver->isProvenOptimal()){
    
    const double * cglpSol = cglpSolver->getColSolution();
    double cutViolation = cglpSolver->getObjValue();

    for(i = 0; i < numULCols; i++){
      cutVals[upperColInd[i]] = cglpSol[i];
      if(0)
	std::cout << upperColInd[i] << ": " << cglpSol[i] << std::endl;
    }
      
    for(i = 0; i < numLLCols; i++){
      cutVals[lowerColInd[i]] = cglpSol[numULCols + i];
      if(0){
	std::cout << lowerColInd[i] << ": " 
		  << cglpSol[numULCols + i] << std::endl;
      }
    }

    cutVals[numULCols + numLLCols] = cglpSol[numULCols + numLLCols];

    std::cout << "Violated inequality found, with violation of " 
	      << cutViolation << "." << std::endl;

  }
  else{

    std::cout << "CGLP failed to find a violated inequality." << std::endl;
    
  }

  delete cglpSolver;

  /* Return the entries of the deepest cut - alpha * x + beta * y >= gamma */
  return cutVals;

}
//#############################################################################
double *
MibSCutGenerator::findDeepestLandPCut1()
{
  // OLD FUNCTION.   
 
  /* Construct the CGLP to find the deepest Lift and Project Cut */

  OsiSolverInterface * cglpSolver = new OsiClpSolverInterface();
  //OsiSolverInterface * cglpSolver = new OsiCbcSolverInterface();
  //OsiSolverInterface * cglpSolver = new OsiSymSolverInterface();
  OsiSolverInterface * oSolver = localModel_->solver();

  int numOrigRows = localModel_->numOrigCons_;
  int numLowerRows = localModel_->lowerRowNum_;
  int numUpperRows = numOrigRows - numLowerRows;
  const int numCols = oSolver->getNumCols();
  const int numRows = oSolver->getNumRows();
  const CoinPackedMatrix * matrix = oSolver->getMatrixByCol();
  const double * matElements = matrix->getElements();
  const int * matIndices = matrix->getIndices();
  const int * matStarts = matrix->getVectorStarts();
  const double * rightHandSide = oSolver->getRightHandSide();
  MibSBilevel *bS = localModel_->bS_;
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  double * lowerObjCoeffs = localModel_->getLowerObjCoeffs();

  int numEntries(numCols + 1);
  double * cutVals = new double[numEntries];
  CoinZeroN(cutVals, numEntries);

  int numULCols = localModel_->getUpperDim();
  int numLLCols = localModel_->getLowerDim();
  double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

  /* the negative and positive slopes of the LL value function */
  double leftSlope = localModel_->leftSlope_;
  double rightSlope = localModel_->rightSlope_;

  /* 
     Column for each variable in the LR + gamma 
     + 2 columns for each constraint,including the disjunction;
     Note: vars associated with = constraints are split into 2 vars
     Row for each variable in the LR + 2 gamma constraints + normalization
  */

  int numCGCols(0);
  numCGCols = numCols + 1 + 4 * numOrigRows + 2 * (numRows - numOrigRows + 2);
  int numCGRows(2 * numCols + 2 + 1);
  
  int numAuxRows(numRows - numOrigRows);

  double * colUb = new double[numCGCols];
  double * colLb = new double[numCGCols];
  double * rowUb = new double[numCGRows];
  double * rowLb = new double[numCGRows];
  double * objCoeffs = new double[numCGCols];
  
  /* 
     variables are ordered: alpha, beta, gamma, u,v,w,z,etc 
     alpha, beta, and gamma are free, rest nonnegative
  */

  CoinFillN(colLb, numCols + 1, - cglpSolver->getInfinity());
  CoinFillN(colLb + numCols + 1, (numCGCols - (numCols + 1)), 0.0);

  CoinFillN(colUb, numCGCols, cglpSolver->getInfinity());

  CoinFillN(rowLb, numCGRows - 1, - cglpSolver->getInfinity());
  CoinFillN(rowUb, numCGRows - 1, 0.0);

  /* normalization constraint */
  rowLb[numCGRows - 1] = 1.0;
  rowUb[numCGRows - 1] = 1.0;

  int i(0), j(0), start(0), end(0);

  CoinZeroN(objCoeffs, numCGCols);

  /* alpha */
  for(i = 0; i < numULCols; i++)
    objCoeffs[i] = upperSol[i];

  /* beta */
  for(i = numULCols; i < numULCols + numLLCols; i++)
    objCoeffs[i] = lowerSol[i - numULCols];

  /* gamma */
  objCoeffs[numULCols + numLLCols] = - 1.0;

  CoinPackedMatrix * cglpMat = new CoinPackedMatrix(false, 0, 0);
  cglpMat->setDimensions(0, numCGCols);

  /* Set up the matrix for the CGLP */

  /* first, the rows constraining alpha */

  int colIndex(0), matIndex(0), insertIndex1(0), insertIndex2(0);
  int size(0), cnt(0);
  double value(0.0);  

  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    row1.insert(i, - 1.0); // alpha coefficient
    row2.insert(i, - 1.0); // alpha coefficient

    size = matrix->getVectorSize(colIndex);

    start = matStarts[colIndex];

    end = start + size;

    for(j = start; j < end; j++){

      matIndex = matIndices[j];
      if(matIndex < numOrigRows){
	
	value = matElements[matIndex + start];
	
	insertIndex1 = numULCols + numLLCols + 1 + matIndex + cnt;

	row1.insert(insertIndex1, value);
	row1.insert(insertIndex1 + 1, - value);

	insertIndex2 = numULCols + numLLCols + 1 + matIndex + cnt 
	  + 2 * numUpperRows + 4 * numLowerRows;

	row2.insert(insertIndex2, value);
	row2.insert(insertIndex2 + 1, - value);

	cnt++;
      }
    }
    cnt = 0;

    //FIXME:  THERE IS A BUG HERE
    // NEED TO ADD IN THE NUMBER OF CUTS THAT HAVE BEEN ADDED
    // TO THE INDICES BELOW
    // RIGHT NOW JUST USING NUM ORIGINAL ROWS
    // INDEX IS OFF
    // 8/26: IS THIS STILL TRUE?  SEEMS TO BE WORKING...

    for(j = start; j < end; j++){

       matIndex = matIndices[j];
       if((matIndex > numUpperRows - 1) && (matIndex < numOrigRows)){

	  value = matElements[matIndex + start];

	  insertIndex1 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex - numUpperRows;
	  row1.insert(insertIndex1, value);
	  
	  insertIndex1 += numLowerRows;
	  row1.insert(insertIndex1, - leftSlope * value);

	  insertIndex2 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex
	    + 2 * numUpperRows + 4 * numLowerRows - numUpperRows;
	  row2.insert(insertIndex2, - value);

	  insertIndex2 += numLowerRows;
	  row2.insert(insertIndex2, - rightSlope * value); 
       }
    }

    for(j = start; j < end; j++){

       matIndex = matIndices[j];
       if(matIndex > numOrigRows - 1){
	  value = matElements[matIndex + start];
	  
	  insertIndex1 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex 
	    - numUpperRows + numLowerRows;
	  row1.insert(insertIndex1, - value);

	  // row1.insert(numULCols + numLLCols + 1 +  
	  //     2 * (numUpperRows + numLowerRows) 
	  //      + numLowerRows + index2 - numUpperRows + numLowerRows,
	  //     - leftSlope * value);
 
	  insertIndex2 = numULCols + numLLCols + 1 + 
	    2 * (numUpperRows + numLowerRows) + matIndex
	    + 2 * numUpperRows + 4 * numLowerRows 
	    - numUpperRows + numLowerRows;
	  row2.insert(insertIndex2, - value);

	  //row2.insert(numULCols + numLLCols + 1 +  
	  //     2 * (numUpperRows + numLowerRows) + numLowerRows 
	  //      + index2 + 2 * numUpperRows + 4 * numLowerRows 
	  //      - numUpperRows + numLowerRows, - rightSlope * value); 
       }
    }


    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }

  /* then, the rows constraining beta */

  for(i = 0; i < numLLCols; i++){
    
    colIndex = lowerColInd[i];

    CoinPackedVector row1;
    CoinPackedVector row2;

    insertIndex1 = i + numULCols;
    row1.insert(insertIndex1, - 1.0); // beta coefficient

    insertIndex2 = i + numULCols;
    row2.insert(insertIndex1, - 1.0); // beta coefficient

    size = matrix->getVectorSize(colIndex);

    start = matStarts[colIndex];

    end = start + size;

    for(j = start; j < end; j++){

      matIndex = matIndices[j];

      //should this be searching lowerRowIndices?
      if(matIndex < numOrigRows){
	
	value = matElements[matIndex + start - numUpperRows];
	
	insertIndex1 = numULCols + numLLCols + 1 
	  + 2 * numUpperRows + matIndex - numUpperRows;

	row1.insert(insertIndex1, value);
	row1.insert(insertIndex1 + 1, - value);
	
	insertIndex2 = numULCols + numLLCols + 1 + matIndex - numUpperRows 
	  + 2 * numUpperRows + 2 * numUpperRows + 4 * numLowerRows;

	row2.insert(insertIndex2, value);
	row2.insert(insertIndex2 + 1, - value);
      }
    }

    
    value = lowerObjCoeffs[i];
    
    insertIndex1 = numULCols + numLLCols + 1 + 
      2 * (numUpperRows + numLowerRows) + numLowerRows;

    row1.insert(insertIndex1, - value);

    insertIndex2 = numULCols + numLLCols + 1 + 
      4 * (numUpperRows + numLowerRows) + 3 * numLowerRows;

    row2.insert(insertIndex2, - value);

    for(j = start; j < end; j++){

      matIndex = matIndices[j];

      if(matIndex > numOrigRows - 1){
	
	value = matElements[matIndex + start - numUpperRows];
	//row1.insert(numULCols + numLLCols + 1 + 2 * numUpperRows 
	//	    + index2 - numUpperRows + 1, value);
	
	insertIndex1 = numULCols + numLLCols + 1 + 
	  2 * (numUpperRows + numLowerRows) + numLowerRows + 1 + matIndex 
	  - numUpperRows;
	
	row1.insert(insertIndex1, - value);

	//row1.insert(numULCols + numLLCols + 1 + index2 + 2 * numUpperRows + 1
	//    - numUpperRows + 1, - value);
	//	row2.insert(numULCols + numLLCols + 1 + index2 - numUpperRows 
	//    + 2 * numUpperRows + 2 * numUpperRows 
	//    + 4 * numLowerRows + 1, value);
	
	insertIndex2 = numULCols + numLLCols + 1 + 
	  4 * (numUpperRows + numLowerRows) + 
	  3 * numLowerRows + 1 + matIndex - numUpperRows;

	row2.insert(insertIndex2, - value);
    
	//row2.insert(numULCols + numLLCols + 1 + index2  - numUpperRows 
	//    + 2 * numUpperRows + 1 + 2 * numUpperRows 
	//    + 4 * numLowerRows + 1, - value);
      }
    }

    cglpMat->appendRow(row1);
    cglpMat->appendRow(row2);

  }

  
  /* then, the rows constraining gamma */

  CoinPackedVector row1;
  CoinPackedVector row2;

  insertIndex1 = numULCols + numLLCols;   
  row1.insert(insertIndex1, 1.0); // gamma coefficient
  
  insertIndex2 = numULCols + numLLCols;
  row2.insert(insertIndex2, 1.0); // gammaa coefficient
  
  for(j = 0; j < numUpperRows; j++){
    
    value = rightHandSide[j];

    insertIndex1 = numULCols + numLLCols + 1 + j;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    insertIndex2 = numULCols + numLLCols + 1 + j
      + 2 * numUpperRows + 4 * numLowerRows;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    
  }
  
  for(j = 0; j < numLowerRows; j++){
    
    value = rightHandSide[numUpperRows + j];
   
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j;

    row1.insert(insertIndex1, - value);
    row1.insert(insertIndex1 + 1, value);

    insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numUpperRows + 4 * numLowerRows;

    row2.insert(insertIndex2, - value);
    row2.insert(insertIndex2 + 1, value);
    
  }

  double mult1(0.0), mult2(0.0);
  
  for(i = 0; i < numULCols; i++){
    
    colIndex = upperColInd[i];
    
    size = matrix->getVectorSize(colIndex);
    
    start = matStarts[colIndex];
    
    end = start + size;
    
    for(j = start; j < end; j++){
      
      matIndex = matIndices[j];
      if((matIndex > numUpperRows - 1) && (matIndex < numOrigRows)){
	value = matElements[matIndex + start];
	mult1 += value * upperSol[i];
	
      }
    }
  }
  
  for(j = 0; j < numLLCols; j++){
      
    //index2 = lowerColInd[j];
    value = lowerObjCoeffs[j];
    mult2 += value * optLowerSol[j];
  }
    

  insertIndex1 = numULCols + numLLCols + 1 + 2 * (numUpperRows + numLowerRows); 
  
  row1.insert(insertIndex1, - mult1);
  
  insertIndex1 += numLowerRows;
  
  row1.insert(insertIndex1, leftSlope * mult1 - mult2); 
  
  insertIndex2 = numULCols + numLLCols + 1 
    + 4 * (numUpperRows + numLowerRows) + 2 * numLowerRows; 

  row2.insert(insertIndex2, mult1);

  insertIndex2 += numLowerRows;

  row2.insert(insertIndex2, rightSlope * mult1 - mult2); 
  
  for(j = 0; j < numAuxRows; j++){
    
    value = rightHandSide[numOrigRows + j];
  
    insertIndex1 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numLowerRows + numUpperRows + numLowerRows; 
    
    row1.insert(insertIndex1, - value);

    //row1.insert(numULCols + numLLCols + 1 + 2 * numUpperRows + j + 1, value);
    
    insertIndex2 = numULCols + numLLCols + 1 + 2 * numUpperRows + j
      + 2 * numUpperRows + 4 * numLowerRows
      + 2 * numLowerRows + numUpperRows + numLowerRows;

    row2.insert(insertIndex2, - value);

    //row2.insert(numULCols + numLLCols + 1 + 2 * numUpperRows + j + 1
    //	+ 2 * numUpperRows + 4 * numLowerRows, value);
    
  }


  
  cglpMat->appendRow(row1);
  cglpMat->appendRow(row2);

  /* then a normalization constraint */

  CoinPackedVector row;

  for(i = numULCols + numLLCols + 1; i < numCGCols; i++)
    row.insert(i, 1.0);

  cglpMat->appendRow(row);

  cglpSolver->assignProblem(cglpMat, colLb, colUb, objCoeffs, rowLb, rowUb);

  if(0)
    cglpSolver->writeLp("cglp");
  //cglpSolver->writeMps("cglp");

  if(0){
    dynamic_cast<OsiClpSolverInterface *> 
      (cglpSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }
  else{
    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (cglpSolver)->setSymParam("max_active_nodes", 1);
  }

  cglpSolver->initialSolve();
  
  if(cglpSolver->isProvenOptimal()){

    const double * cglpSol = cglpSolver->getColSolution();
    double cutViolation = cglpSolver->getObjValue();
    /*
    for(i = 0; i < numEntries; i++)
      cutVals[i] = cglpSol[i];
    */

    for(i = 0; i < numULCols; i++)
      cutVals[upperColInd[i]] = cglpSol[i];

    for(i = 0; i < numLLCols; i++)
      cutVals[lowerColInd[i]] = cglpSol[numULCols + i];

    cutVals[numULCols + numLLCols] = cglpSol[numULCols + numLLCols];

    std::cout << "Violated inequality found, with violation of " 
	      << cutViolation << "." << std::endl;

  }
  else{

    std::cout << "CGLP failed to find a violated inequality." << std::endl;
    
  }

  delete cglpSolver;

  /* Return the entries of the deepest cut - alpha * x + beta * y >= gamma */
  return cutVals;

}

//#############################################################################
int
MibSCutGenerator::interdictionCuts(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIPINT Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();
  
  int useBendersCut = 
     localModel_->MibSPar_->entry(MibSParams::useBendersCut);

  int numCuts(0);
  MibSTreeNode * node = 
    dynamic_cast<MibSTreeNode *>(localModel_->activeNode_); 
  double maxLowerObj(node->getLowerUB());
  double etol(localModel_->etol_);
  int lN(localModel_->getLowerDim());
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
 
  int i(0), index(0);
  double cutub(- 1.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  for(i = 0; i < uN; i++){
    index = upperColInd[i];
    indexList.push_back(index);
    if(sol[index] > etol){
      valsList.push_back(1.0);
      cutub += 1.0;
    }
    else{
      valsList.push_back(- 1.0);
    }
  }

#if 0
  indexList.push_back(uN + lN);
  valsList.push_back(- 1.0);
#endif

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, false);
  
  indexList.clear();
  valsList.clear();

#if 0
  double m(solver->getObjValue() - maxLowerObj);
  //double m(maxLowerObj - 
  cutub = m + maxLowerObj;

  for(i = 0; i < lN; i++){
    index = lowerColInd[i];
    if(fabs(lObjCoeffs[i]) != 0.0){
      indexList.push_back(index);
      valsList.push_back(lObjCoeffs[i]);
    }
  }
  
  if(fabs(m) != 0.0){
    indexList.push_back(uN + lN);
    valsList.push_back(m);
  }

  //for(i = 0; i < valsList.size(); i++){
  //  std::cout << "index: " << indexList.at(i) << std::endl;
  //  std::cout << "value: " << valsList.at(i);
  //}

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, false);
#endif

  if (useBendersCut){
     numCuts += bendersInterdictionCuts(conPool);
  }

  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::weakIncObjCutCurrent(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIPINT Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();

  int numCuts(0);
  MibSBilevel * bS = localModel_->bS_;

  //double maxLowerObj(node->getLowerUB());
  double minLowerObj(solver->getInfinity());
  //double maxLowerObj = localModel_->lowerObjectiveBound(true);
  double etol(localModel_->etol_);
  int lN(localModel_->getLowerDim());
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
  double lowerObjSense = localModel_->getLowerObjSense(); 

  double * tmpsol = new double[lN + uN];
  CoinZeroN(tmpsol, lN + uN);
  OsiSolverInterface * lSolver = bS->setUpModel(solver, true, tmpsol);

  if(0){
    dynamic_cast<OsiCbcSolverInterface *> 
      (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }    
  else{
    dynamic_cast<OsiSymSolverInterface *> 
      (lSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (lSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (lSolver)->setSymParam("max_active_nodes", 1);
  }

  if(0)
    lSolver->writeLp("boundfinder");

  lSolver->branchAndBound();
  minLowerObj = lSolver->getObjValue();

  //if it's exactly zero, just add a little bit to comply with blis
  if(minLowerObj == 0.0)
    minLowerObj = minLowerObj + etol;

  if(0)
    std::cout << "minLowerObj " << minLowerObj << std::endl;

  double bigM(10000);
  int i(0), index(0);
  double cutub(0.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  for(i = 0; i < uN; i++){
    index = upperColInd[i];
    if(sol[index] < etol){
      indexList.push_back(index);
      //valsList.push_back(minLowerObj);
      valsList.push_back(-bigM);
    }
  }

  double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 
  
  for(i = 0; i < lN; i++){
    index = lowerColInd[i];
    //cutub += lObjCoeffs[i] * sol[index]; 
    if(fabs(lObjCoeffs[i]) > etol){
      cutub += lObjCoeffs[i] * optLowerSol[i]; //should this be position? 
      indexList.push_back(index);
      valsList.push_back(lowerObjSense * lObjCoeffs[i]);
    }
  }

  assert(indexList.size() == valsList.size());

  //numCuts += addCut(conPool, cutlb, 
  //	    lowerObjSense * cutub, indexList, valsList, false);
  numCuts += addCut(conPool, cutlb, 
		    lowerObjSense * cutub, indexList, valsList, true);

  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::weakIncObjCutMaximal(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating MIPINT Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();

  int numCuts(0);
  MibSBilevel * bS = localModel_->bS_;
  //MibSTreeNode * node = 
  //dynamic_cast<MibSTreeNode *>(localModel_->activeNode_); 
  //double maxLowerObj(node->getLowerUB());
  double maxLowerObj(solver->getInfinity());
  //double maxLowerObj = localModel_->lowerObjectiveBound(true);
  double etol(localModel_->etol_);
  int lN(localModel_->getLowerDim());
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
  //const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
  double lowerObjSense = localModel_->getLowerObjSense(); 

  //if we keep this as our big M, should move (same every time)
  double * tmpsol = new double[lN + uN];
  CoinZeroN(tmpsol, lN + uN);
  OsiSolverInterface * lSolver = bS->setUpModel(solver, tmpsol);
  delete [] tmpsol;
  if(0){
    dynamic_cast<OsiCbcSolverInterface *> 
      (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }
  else{
    dynamic_cast<OsiSymSolverInterface *> 
      (lSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (lSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (lSolver)->setSymParam("max_active_nodes", 1);
  }
  //lSolver->setObjSense(lSolver->getObjSense());
  lSolver->branchAndBound();
  maxLowerObj = lSolver->getObjValue();

  int i(0), index(0);
  double cutub(0.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  if(0)
    std::cout << "maxLowerObj " << maxLowerObj << std::endl;

  const double * sol = findMaximalUpperSol(solver);

  //const double * maximalupper = findMaximalUpperSol(solver);

  if(!sol){

    return 0;
  }
  else{
  
    OsiSolverInterface * lSolver = bS->setUpModel(solver, sol);  

    if(0){
      dynamic_cast<OsiCbcSolverInterface *> 
	(lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }
    else{
      dynamic_cast<OsiSymSolverInterface *> 
	 (lSolver)->setSymParam("prep_level", -1);
      
      dynamic_cast<OsiSymSolverInterface *> 
	 (lSolver)->setSymParam("verbosity", -2);

      dynamic_cast<OsiSymSolverInterface *> 
	 (lSolver)->setSymParam("max_active_nodes", 1);
    }

    lSolver->branchAndBound();

    const double * optLowerSol = lSolver->getColSolution();
    
    for(i = 0; i < uN; i++){
      index = upperColInd[i];
      if(sol[index] < etol){
	indexList.push_back(index);
	valsList.push_back(maxLowerObj);
      }
    }

    for(i = 0; i < lN; i++){
      index = lowerColInd[i];
      //cutub += lObjCoeffs[i] * sol[index]; 
      if(fabs(lObjCoeffs[i]) > etol){
	cutub += lObjCoeffs[i] * optLowerSol[i]; 
	indexList.push_back(index);
	valsList.push_back(lowerObjSense * lObjCoeffs[i]);
      }
    }
    
    assert(indexList.size() == valsList.size());
    
    //numCuts += addCut(conPool, cutlb, 
    //    lowerObjSense * cutub, indexList, valsList, false);
    numCuts += addCut(conPool, cutlb, 
		      lowerObjSense * cutub, indexList, valsList, true);

    maximalCutCount_++;

    delete [] sol;
    return numCuts;
  }

}

//#############################################################################
int
MibSCutGenerator::binaryCuts(BcpsConstraintPool &conPool)
{

  //FIXME: NEED TO CHECK FOR ROW TYPES FOR CGLP.  
  //COMING FROM KNAP SOLVER THEY ARE RANGED.

  bool useNoGoodCut 
    = localModel_->MibSPar_->entry(MibSParams::useNoGoodCut);

   bool useIncObjCut 
    = localModel_->MibSPar_->entry(MibSParams::useIncObjCut);

  if(useNoGoodCut && !useIncObjCut){
    return noGoodCut(conPool) ? true : false;
  }
  else if(!useNoGoodCut && useIncObjCut){
    return incObjCut(conPool) ? true : false;
  }
  else if(useNoGoodCut && useIncObjCut){
    return (noGoodCut(conPool) && 
	    incObjCut(conPool)) ? true : false;
  }
  else{
    //std::cout << "No BINARY Cuts generated" << std::endl;
    return 0;
  }
}

//#############################################################################
int
MibSCutGenerator::incObjCut(BcpsConstraintPool &conPool)
{

  int numCuts(0);

  //numCuts += incObjCutCurrent(conPool);

  //numCuts += incObjCutMaximal(conPool);

  numCuts += weakIncObjCutCurrent(conPool);
  
  //if(!maximalCutCount_)
  //numCuts += weakIncObjCutMaximal(conPool);

  return numCuts;

}

//#############################################################################
int
MibSCutGenerator::incObjCutMaximal(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  MibSBilevel * bS = localModel_->bS_;
  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  int i(0);
  int numCuts(0);
  double etol(localModel_->etol_);
  int uCols = localModel_->getUpperDim();
  int lCols = localModel_->getLowerDim();
 
  //ACTUALLY, WE WANT THE FULL SOLUTION FROM FIND MAXIMAL
  //THIS GIVES US A FULL CURRENT SOLUTION
  //THEN, WANT TO TAKE UPPER PORTION AND FIND OPTIMAL LL SOL
  //SHOULD BE USING setUpModel(solver, sol) (from MibSBilevel?)

  //figure out how these need to be ordered, cglp is expecting ordered sols

  const double * maximalupper = findMaximalUpperSol(solver);

  double * upperSol = new double[uCols];
  double * lowerSol = new double[lCols];
  double * optLowerSol = new double[lCols];
  CoinZeroN(upperSol, uCols);
  CoinZeroN(lowerSol, lCols);
  CoinZeroN(optLowerSol, lCols);

  /*  
  std::cout << "The maximal solution is: " << std::endl;

  for(i = 0; i < uCols; i++){
    upperSol[i] = maximalupper[uColInd[i]];
    std::cout << uColInd[i] << ": " << upperSol[i] << endl;
  }

  for(i = 0; i < uCols; i++){
    lowerSol[i] = maximalupper[lColInd[i]];
    std::cout << lColInd[i] << ": " << lowerSol[i] << endl;
  }
  */
    //double * lSol = new double[localModel_->getLowerDim()];
    //CoinZeroN(lSol, localModel_->getLowerDim());

  if(!maximalupper){

    return 0;
  }
  else{
  
    OsiSolverInterface * lSolver = bS->setUpModel(solver, maximalupper);  

    if(0){
      dynamic_cast<OsiCbcSolverInterface *> 
	(lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }
    else{
      dynamic_cast<OsiSymSolverInterface *> 
	 (lSolver)->setSymParam("prep_level", -1);
      
      dynamic_cast<OsiSymSolverInterface *> 
	 (lSolver)->setSymParam("verbosity", -2);

      dynamic_cast<OsiSymSolverInterface *> 
	 (lSolver)->setSymParam("max_active_nodes", 1);
    }
    lSolver->branchAndBound();

    int lCols = localModel_->getLowerDim();
    
    //CoinCopyN(lSolver->getColSolution(), localModel_->getLowerDim(), lSol);
    for(i = 0; i < lCols; i++)
      optLowerSol[i] = lSolver->getColSolution()[i];

    double cutlb(- solver->getInfinity());
    double cutub(solver->getInfinity());
    
    std::vector<int> indexList;
    std::vector<double> valsList;
    
    //replace these with results from above
    //NOT CORRECT NOW
    //MibSBilevel * bS = localModel_->bS_; //not needed once fixed below
    //double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
    //double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
    //double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

    /* 
       returns a double with the values [alpha | beta | gamma]
       indexed 0..n1-1, n1..n1+n2-1, n1+n2
    */
    
    double * cutVals = findDeepestLandPCut_IncObj(upperSol, 
						  lowerSol, 
						  optLowerSol);  
    double val(0.0);
    
    /* find the nonzero entries */
    if(0)
      std::cout << "maximal solution." << std::endl;

    for(i = 0; i < numCols; i++){
      val = cutVals[i];
      if((val > etol) || (val < - etol)){
	// entry is nonzero
	indexList.push_back(i);
	if(0)
	  std::cout << "i: " << i << " -- "  << val << std::endl;
	valsList.push_back(val);
      }
    }
    
    cutub = cutVals[numCols]; //gamma

    if(0)
      std::cout << "cutub: " << cutub << std::endl;
    
    assert(indexList.size() == valsList.size());
    
    numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);
    
    delete [] cutVals;
    return numCuts;
  }

}

//#############################################################################
const double *
MibSCutGenerator::findMaximalUpperSol(OsiSolverInterface * si)
{

  OsiSolverInterface * maxSolver = new OsiCbcSolverInterface(si);

  int numCols = maxSolver->getNumCols();
  int numULCols = localModel_->getUpperDim();
  int * upperColInd = localModel_->getUpperColInd();

  //  double * maximalsol = new double[numCols];
  //CoinZeroN(maximalsol, numCols);

  int i(0);
  double * objectives = new double[numCols];
  CoinZeroN(objectives, numCols);

  //maximizing the sum of the upper-level variables 
  for(i = 0; i < numULCols; i++)
    objectives[upperColInd[i]] = 1.0;
 
  maxSolver->setObjective(objectives);
  maxSolver->setObjSense(-1); //maximization

  if(0){
    dynamic_cast<OsiCbcSolverInterface *> 
      (maxSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
  }
  else{
    dynamic_cast<OsiSymSolverInterface *> 
      (maxSolver)->setSymParam("prep_level", -1);
    
    dynamic_cast<OsiSymSolverInterface *> 
      (maxSolver)->setSymParam("verbosity", -2);

    dynamic_cast<OsiSymSolverInterface *> 
      (maxSolver)->setSymParam("max_active_nodes", 1);
  }
  maxSolver->branchAndBound();
  //maxSolver->initialSolve();

  if(0)
    maxSolver->writeLp("maxsolver");

  if(maxSolver->isProvenOptimal()){

    double * maxsol = new double[maxSolver->getNumCols()];
    CoinCopyN(maxSolver->getColSolution(), maxSolver->getNumCols(), maxsol);
    
    if(1){
      std::cout << "Maximal UL Solution found; cardinality of " 
		<< maxSolver->getObjValue() << std::endl;
    }

    delete [] objectives;
    //delete maxSolver;

    return (maxsol);

  }
  else{

    std::cout << "Unable to find maximal solution." << std::endl;
    delete [] objectives;
    //delete maxSolver;
    return NULL;

  }


}

//#############################################################################
int
MibSCutGenerator::incObjCutCurrent(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  OsiSolverInterface * solver = localModel_->solver();
  MibSBilevel * bS = localModel_->bS_;

  const int numCols = solver->getNumCols();
  int i(0);
  int numCuts(0);
  double etol(localModel_->etol_);
 
  double cutlb(- solver->getInfinity());
  double cutub(solver->getInfinity());

  std::vector<int> indexList;
  std::vector<double> valsList;

  double * upperSol = bS->upperSolutionOrd_; // UL portion from LR soln
  double * lowerSol = bS->lowerSolutionOrd_; // LL portion from LR soln
  double * optLowerSol = bS->optLowerSolution_; // optimal LL solution 

  /* 
     returns a double with the values [alpha | beta | gamma]
     indexed 0..n1-1, n1..n1+n2-1, n1+n2
  */

  double * cutVals = findDeepestLandPCut_IncObj(upperSol, 
						lowerSol, 
						optLowerSol);  
  double val(0.0);

  /* find the nonzero entries */

  for(i = 0; i < numCols; i++){
    val = cutVals[i];
    if((val > etol) || (val < - etol)){
      // entry is nonzero
      indexList.push_back(i);
      valsList.push_back(val);
    }
  }

  cutub = cutVals[numCols]; //gamma

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, true);

  delete [] cutVals;
  return numCuts;

}

//###########################################################################
int
MibSCutGenerator::noGoodCut(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating No-Good Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();

  int numCuts(0);
  double etol(localModel_->etol_);
  const double * sol = solver->getColSolution();
  int uN(localModel_->upperDim_);
  int * upperColInd = localModel_->getUpperColInd();
 
  int i(0), index(0);
  double cutub(- 1.0);
  double cutlb(- solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  for(i = 0; i < uN; i++){
    index = upperColInd[i];
    indexList.push_back(index);
    if(sol[index] > etol){
      valsList.push_back(1.0);
      cutub += 1.0;
    }
    else{
      valsList.push_back(- 1.0);
    }
  }

  assert(indexList.size() == valsList.size());

  numCuts += addCut(conPool, cutlb, cutub, indexList, valsList, false);
  
  return numCuts;

}

//###########################################################################
int
MibSCutGenerator::bendersInterdictionCuts(BcpsConstraintPool &conPool)
{

  /** Add specialized bilevel feasibility cuts, as appropriate **/

  //std::cout << "Generating No-Good Cuts." << std::endl;

  OsiSolverInterface * solver = localModel_->solver();
  OsiSolverInterface * lSolver = localModel_->getMibSBilevel()->solver_;

  const double *lpSol = solver->getColSolution();
  
  sym_environment *env = 
     dynamic_cast<OsiSymSolverInterface *>(lSolver)->getSymphonyEnvironment();
  
  int sp_size(0);
  sym_get_sp_size(env, &sp_size);
  
  int numCuts(0);
  double etol(localModel_->etol_);
  int uN(localModel_->upperDim_);
  int lN(localModel_->getLowerDim());
  int * upperColInd = localModel_->getUpperColInd();
  int * lowerColInd = localModel_->getLowerColInd();
  double * lObjCoeffs = localModel_->getLowerObjCoeffs();
 
  int i(0), j(0);
  double cutub(solver->getInfinity());
  std::vector<int> indexList;
  std::vector<double> valsList;

  int indexU(0), indexL(0);

#if 0
  double ll_objval;
  double * ll_sol = new double[lN];

  for (i = 0; i < sp_size; i++){
     sym_get_sp_solution(env, i, ll_sol, &ll_objval);
     for(j = 0; j < uN; j++){
	indexU = upperColInd[j];
	indexL = lowerColInd[j];
	indexList.push_back(indexL);
	valsList.push_back(-lObjCoeffs[j]);
	if (ll_sol[j] > etol){
	   indexList.push_back(indexU);
	   valsList.push_back(-lObjCoeffs[j]);
	}
     }
     assert(indexList.size() == valsList.size());
     numCuts += addCut(conPool, -ll_objval, cutub, indexList, valsList, 
		       false);
     indexList.clear();
     valsList.clear();
  }

  delete[] ll_sol;
#endif

  std::vector<std::pair<AlpsKnowledge*, double> > solutionPool;
  localModel_->getKnowledgeBroker()->
     getAllKnowledges(AlpsKnowledgeTypeSolution, solutionPool);
  
  const double * sol; 
  double objval, lhs(0);
  BlisSolution* blisSol;
  std::vector<std::pair<AlpsKnowledge*, double> >::const_iterator si;
  for (si = solutionPool.begin(); si != solutionPool.end(); ++si){
     lhs = 0;
     blisSol = dynamic_cast<BlisSolution*>(si->first);
     sol = blisSol->getValues();
     objval = blisSol->getQuality();
     for(j = 0; j < uN; j++){
	indexU = upperColInd[j];
	indexL = lowerColInd[j];
	indexList.push_back(indexL);
	valsList.push_back(-lObjCoeffs[j]);
	lhs -= lObjCoeffs[j]*lpSol[indexL];
	if (sol[indexL] > etol){
	   indexList.push_back(indexU);
	   valsList.push_back(-lObjCoeffs[j]*sol[indexL]);
	   lhs -= lObjCoeffs[j]*lpSol[indexU];
	}
     }
     assert(indexList.size() == valsList.size());
     if (lhs < objval){
	numCuts += addCut(conPool, objval, cutub, indexList, valsList, 
			  false);
     }
     indexList.clear();
     valsList.clear();
  }

  return numCuts;

}

//###########################################################################
bool 
MibSCutGenerator::generateConstraints(BcpsConstraintPool &conPool)
{

  //FIXME: MAKE THIS MORE SIMPLE
  int numCuts(0);

  int cutStrategy =
    localModel_->MibSPar_->entry(MibSParams::cutStrategy);

  if(cutStrategy > 0){
     //use bilevel feasibility cut
    MibSBilevel *bS = localModel_->bS_;
    //int problemType = 
    //  localModel_->MibSPar_->entry(MibSParams::bilevelProblemType);
    int cutTypes = 
      localModel_->MibSPar_->entry(MibSParams::bilevelCutTypes);
    
    CoinPackedVector *sol = localModel_->getSolution();

    if(localModel_->solIsUpdated_)
      bS = localModel_->bS_;
    else
      bS->createBilevel(sol, localModel_);

    localModel_->solIsUpdated_ = false;

    if(bS->isIntegral_ && cutTypes == 0){
      //general type of problem, no specialized cuts
      delete sol;
      return feasibilityCuts(conPool) ? true : false;
    }
    else if(bS->isIntegral_ && cutTypes == 1){
      //interdiction problem
      delete sol;
      return (feasibilityCuts(conPool) &&
	      interdictionCuts(conPool)) ? true : false;
    }
    else if(bS->isIntegral_ && cutTypes == 2){
      //problem with binary UL variables and integer LL variables
      delete sol;
      return (feasibilityCuts(conPool) &&
	      binaryCuts(conPool)) ? true : false;
    }
    else if(bS->isIntegral_ && cutTypes == 3){
      //problem with binary UL variables and general LL variables
      delete sol;
      return binaryCuts(conPool) ? true : false;
    }
  
    /*
      if(bS->isIntegral_){
      delete sol;
      return feasibilityCuts(conPool) ? true : false;
      }
    */
    
    /*
      May want to add our own cuts later to separate fractional solutions
    */
    
    delete sol;
  }

  return numCuts ? true : false;

}

//#############################################################################
int 
MibSCutGenerator::addCut(BcpsConstraintPool &conPool,
			 double cutlb, double cutub,
			 std::vector<int> & indexList, 
			 std::vector<double> &elementList,
			 bool removable)
{

  int size(indexList.size());
  int * indices = new int[size];
  double * elements = new double[size];
  int capacity(localModel_->solver()->getNumCols());

  int i(0);

  for(i = 0; i < size; i++){
    indices[i] = indexList.at(i);
    elements[i] = elementList.at(i);
  }

  OsiRowCut *cut = new OsiRowCut(cutlb, cutub, capacity, 
				 size, indices, elements);

  BlisConstraint *blisCon = BlisOsiCutToConstraint(cut);

  if(!removable)
    blisCon->setStatus(BCPS_NONREMOVALBE);
  
  conPool.addConstraint(blisCon);

  delete cut;
  return 1;

}

//#############################################################################
int *
MibSCutGenerator::getBindingCons()
{

   //std::string method = model_->bindingMethod_; 
   int method = localModel_->MibSPar_->entry(MibSParams::whichActiveConMethod);
   int * binding = NULL;

  if(method == 0){
     binding = getBindingConsSimple();
  }
  else if(method == 1){
     binding = getBindingConsBasis();
  }
  else{
     std::cout << "No method for binding constraint set." << std::endl;
     return NULL;
  }
  
  return binding;
  
}

//#############################################################################
int *
MibSCutGenerator::getBindingConsSimple()
{

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const char * rowsense = solver->getRowSense();
  const double * rowrange = solver->getRowRange();
  const double * sol = solver->getColSolution();
  const double * rhs = solver->getRightHandSide();
  int i(0), index(0);
  double etol(localModel_->etol_);

  //FIXME: CHECK THIS
  const double * collb;
  const double * colub;

  if(localModel_->isRoot_){
    collb = localModel_->getSolver()->getColLower();
    colub = localModel_->getSolver()->getColUpper();
  }
  else{
    collb = solver->getColLower();
    colub = solver->getColUpper();
  }

  int *colStatus = new int[numCols];
  int *rowStatus = new int[numRows];
  solver->getBasisStatus(colStatus, rowStatus);

  /** Find all binding constraints at current solution **/

  int bindingCons(0);
  double slack(0.0);
  double upper(0.0);
  
  int * binding = new int[numRows + 2 * numCols];
  double * slackVal = new double[numRows + 2 * numCols];
  const double * rowActivity = solver->getRowActivity();
  int MAXCONS(numCols);
  //int MAXCONS(100);  

  CoinZeroN(binding, numRows + 2 * numCols);
  CoinFillN(slackVal, numRows + 2 * numCols, 0.0);
  
  for(i = 0; i < numRows; i++){
     
    if(1){//if(rowStatus[i] != 3) && (rowStatus[i] != 2)){
      if((rowsense[i] =='R') && (rowrange[i] < 0))
	slackVal[i] = - rhs[i] + rowActivity[i];
      else
	slackVal[i] = rhs[i] - rowActivity[i];
      
      slack = slackVal[i];
      
      switch(rowsense[i]){
      case 'L':
	if(slack > etol){
	  if(0) 
	    std::cout << "constraint " << i << " is not tight." << std::endl;
	  binding[i] = 0;
	}
	else if((slack > - etol) && (bindingCons < MAXCONS)){
	  if(0){
	    std::cout << "constraint " << i << " is tight." << std::endl;
	    std::cout << "status is " << rowStatus[i] << std::endl;
	  }
	  binding[i] = 1;
	  upper += rhs[i];
	  //upper += rhs[i] * solver->getRowPrice()[i];
	  bindingCons++;
	}
	break;
      case 'G':
	if(slack < - etol){
	  if(0) 
	    std::cout << "constraint " << i << " is not tight." << std::endl;
	  binding[i] = 0;
	}
	else if((slack < etol) && (bindingCons < MAXCONS)){
	  if(0){
	    std::cout << "constraint " << i << " is tight." << std::endl;
	    std::cout << "status is " << rowStatus[i] << std::endl;
	  }
	  binding[i] = 1;
	  upper -= rhs[i];
	  //upper -= rhs[i] * solver->getRowPrice()[i];
	  bindingCons++;
	}
	break;
      case 'E':
	if((slack < - etol) || (slack > etol)){
	  if(0) 
	    std::cout << "constraint " << i << " is not tight." << std::endl;
	  binding[i] = 0;
	}
	else if(bindingCons < MAXCONS){
	  if(0){
	    std::cout << "constraint " << i << " is tight." << std::endl;
	    std::cout << "status is " << rowStatus[i] << std::endl;
	  }
	  binding[i] = 1;
	  upper += rhs[i];
	  //upper += rhs[i] * solver->getRowPrice()[i];
	  bindingCons++;
	}
	break;
      case 'R':
	if(rowrange[i] < 0){
	  if((slack > rowrange[i] + etol) || (slack < - etol)){
	    //std::cout << "constraint " 
	    //	<< i << " is not tight." << std::endl;
	    binding[i] = 0;
	  }
	  else if(((slack > rowrange[i] - etol) || (slack < etol)) 
		  && (bindingCons < MAXCONS)){
	    if(0){
	      std::cout << "constraint " << i << " is tight." << std::endl;
	      std::cout << "status is " << rowStatus[i] << std::endl;
	    }	   
	    binding[i] = 1;
	    upper += rhs[i];
	    //upper += rhs[i] * solver->getRowPrice()[i];
	    bindingCons++;
	  }
	}
	else{
	  if((slack < rowrange[i] - etol) || (slack > etol)){
	    if(0)
	      std::cout << "constraint " 
			<< i << " is not tight." << std::endl;
	    binding[i] = 0;
	  }
	  else if(((slack < rowrange[i] + etol) || (slack > - etol)) 
		  && (bindingCons < MAXCONS)){
	    if(0){
	      std::cout << "constraint " << i << " is tight." << std::endl;
	      std::cout << "status is " << rowStatus[i] << std::endl;
	    }
	    binding[i] = 1;
	    upper += rhs[i];
	    //upper += rhs[i] * solver->getRowPrice()[i];
	    bindingCons++;
	  }
	}
	break;
      }
    }
  }

  double upSlack(0.0);
  double downSlack(0.0);

  for(i = 0; i < numCols; i++){
    if(0){
      std::cout << "status of column " << i << ": " 
		<< colStatus[i] << std::endl;
    }
    index = numRows + i;
   
    /*   
    if((fabs(colub[i] - collb[i]) < etol) && (colStatus[i] == 3)){
      if(1)
	std::cout << "col " << i << " is fixed with status." << colStatus[i] 
		  << std::endl;
      //binding[index] = 1;
      //upper += colub[i];
      //bindingCons++;
    }
    */
    //else{// 
    if(1){//if(colStatus[i] != 1){ 
      upSlack = colub[i] - sol[i];
      downSlack = sol[i] - collb[i];
      
      if((upSlack > - etol) && (upSlack < etol)){
	if(0){
	  std::cout << "upper bound " << i << " is tight." << std::endl;
	  std::cout << "status is " << colStatus[i] << std::endl;
	}
	binding[index] = 1;
	upper += colub[i];
	bindingCons++;
      }
      else if((downSlack > - etol) && (downSlack < etol)){
	if(0){
	  std::cout << "lower bound " << i << " is tight." << std::endl;
	  std::cout << "status is " << colStatus[i] << std::endl;
	}
	binding[index + numCols] = 1;
	upper -= collb[i];
	bindingCons++;
      }
    }
  }

  
  if(!bindingCons){
     std::cout << "No constraints are binding." << std::endl;
     //return 0;
     abort();
  }
  
  /** Set the upper bound for the current cut **/ 
  setCutUpperBound(upper);
  
  delete [] slackVal;  

  return binding;
  
}

//#############################################################################
int *
MibSCutGenerator::getBindingConsBasis()
{

  OsiSolverInterface * solver = localModel_->solver();

  const int numCols = solver->getNumCols();
  const int numRows = solver->getNumRows();
  const char * rowsense = solver->getRowSense();
  const double * collb = solver->getColLower();
  const double * colub = solver->getColUpper();
  const double * rhs = solver->getRightHandSide();
  int i(0), index(0);

  /** Find binding constraints at current solution using basis information **/

  int bindingCons(0);
  double upper(0.0);

  int * binding = new int[numRows + 2 * numCols];
  CoinZeroN(binding, numRows + 2 * numCols);

  CoinWarmStartBasis::Status rowStatus;
  CoinWarmStartBasis * ws = 
     dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
  
  for(i = 0; i < numRows; i++){
    
    rowStatus = ws->getArtifStatus(i);
    
     switch(rowsense[i]){
      case 'L':
	if(rowStatus == CoinWarmStartBasis::basic){
	  binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper += rhs[i];
	  bindingCons++;
	}
	break;
      case 'G':
	if(rowStatus == CoinWarmStartBasis::basic){
	  binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper -= rhs[i];
	  bindingCons++;
	}
	break;
      case 'E':
	if(rowStatus == CoinWarmStartBasis::basic){
	   binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper += rhs[i];
	  bindingCons++;
	}
	break;
      case 'R':
	if(rowStatus == CoinWarmStartBasis::basic){
	   binding[i] = 0;
	}
	else{
	  //std::cout << "artificial " << i << " is nonbasic." << std::endl;
	  binding[i] = 1;
	  upper += rhs[i];
	  bindingCons++;
	}
	break;
     }
  }
  
  /** Find binding bounds at current solution **/

  const double * sol = solver->getColSolution();
  double etol(localModel_->etol_);

  double upSlack(0.0);
  double downSlack(0.0);
  
  for(i = 0; i < numCols; i++){
    index = numRows + i;
    upSlack = colub[i] - sol[i];
    downSlack = sol[i] - collb[i];
    
    if((upSlack > - etol) && (upSlack < etol)){
      //std::cout << "upper bound " << i << " is tight." << std::endl;
      binding[index] = 1;
      upper += colub[i];
      bindingCons++;
     }
    else if((downSlack > - etol) && (downSlack < etol)){
      binding[index + numCols] = 1;
      upper -= collb[i];
      bindingCons++;
    }
  }

  
  if(!bindingCons){
    std::cout << "Invalid basis." << std::endl;
    //return 0;
    abort();
  }
  
  /** Set the upper bound for the current cut **/ 
  setCutUpperBound(upper);
  
  return binding;

}
