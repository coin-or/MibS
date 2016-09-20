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

#ifndef MibsCutGenerator_h_
#define MibsCutGenerator_h_

#include "MibSModel.hpp"

//#############################################################################

class MibSCutGenerator : public BlisConGenerator {

 private:

   MibSModel * localModel_;
   double upper_;
   int auxCount_;
   int maximalCutCount_;
   
 public:
   
   /** Default Constructor **/
   MibSCutGenerator(MibSModel *mibs=0);
   
   /** Destructor **/
   ~MibSCutGenerator();
   
   bool generateConstraints(BcpsConstraintPool &conPool);

   int addCut(BcpsConstraintPool &conPool,
	      double cutlb, double cutub,
	      std::vector<int> & indexList, 
   	      std::vector<double> &elementList,
	      bool removable);
   
   void setLocalModel(MibSModel *mibs) { localModel_ = mibs; }
   
   /** Add a cut pi_1x + pi_2y <= pi_0 - 1**/
   int bilevelFeasCut1(BcpsConstraintPool &conPool);
   
   /** Add a a lift-and-project cut for bilevel infeasible solution **/
   int bilevelFeasCut2(BcpsConstraintPool &conPool);
   
   /** Add appropriate bilevel feasiblity cuts **/
   int feasibilityCuts(BcpsConstraintPool &conPool);

   /** Add specialized bilevel feasiblity cuts **/
   int interdictionCuts(BcpsConstraintPool &conPool);

   /** Add specialized bilevel feasiblity cuts **/
   int binaryCuts(BcpsConstraintPool &conPool);

   /** Add no-good cuts for binary upper-level variables **/
   int noGoodCut(BcpsConstraintPool &conPool);

   /** Add Benders-type cuts for interdiction problems **/
   int bendersInterdictionCuts(BcpsConstraintPool &conPool);

   /** Add Benders-type cuts for zero sum problems **/
   int bendersZeroSumCuts(BcpsConstraintPool &conPool);

   /** Add intersection cuts for discrete problems **/
   int intersectionCuts(BcpsConstraintPool &conPool, double *optLowerSolution);

   /** Helper function for intersection cut*/
   void getAlphaIntersectionCut(double** extRay, double* lowerSolution,
				int numStruct, int numNonBasic,
				const double* lpSol, std::vector<double> &alphaVec);

   /** Helper function for intersection cut*/
   double solveModelIntersectionCut(const CoinPackedMatrix* matrix,
				    double** extRay, double* rowLb,double* rowUb,
				    int lRows, int numRows, int numNonBasic, int cnt);
		     
   /** Add bound cuts for general problems **/
   int boundCuts(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables **/
   int incObjCut(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables (current sol)**/
   int incObjCutCurrent(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables (maximal sol) **/
   int incObjCutMaximal(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables (current sol)**/
   int weakIncObjCutCurrent(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables (maximal sol) **/
   int weakIncObjCutMaximal(BcpsConstraintPool &conPool);

   /** Use the cut generator LP to find the deepest L and P cut **/
   double * findDeepestLandPCut_ValFunc();

   /** Use the cut generator LP to find the deepest L and P cut **/
   double * findDeepestLandPCut_IncObj(double * uppersol, 
				       double * lowersol,
				       double * optlowersol);
   
   /** Find the maximal UL binary solution over current constraint region **/
   const double * findMaximalUpperSol(OsiSolverInterface * si);

   /** old function, delete eventually **/
   double * findDeepestLandPCut1();//stable (but maybe wrong)
   
 private:

   /** Set the upper bound for the current cut **/
   inline void setCutUpperBound(double val) {upper_ = val;}
   
   /** Get the upper bound for the current cut **/
   double getCutUpperBound() {return upper_;}
   
   /** Get the binding constraints at the current solution **/
   int * getBindingCons();

   /** Get all the binding constraints **/
   int * getBindingConsSimple();

   /** Get the binding constraints using basis information **/
   int * getBindingConsBasis();

};

#endif
