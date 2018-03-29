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
   bool isBigMIncObjSet_;
   double bigMIncObj_;
   OsiSolverInterface * watermelonICSolver_;
    
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

   int generalNoGoodCut(BcpsConstraintPool &conPool);

   /** Add Benders-type cuts for interdiction problems **/
   int bendersInterdictionOneCut(BcpsConstraintPool &conPool,
				   double *lSolution);

   int bendersInterdictionMultipleCuts(BcpsConstraintPool &conPool);

   /** Add Benders-type cuts for zero sum problems **/
   int bendersZeroSumCuts(BcpsConstraintPool &conPool);

    /** Add intersection cuts for general problems (IC: discrete, hypercube, tender: general) **/
    int intersectionCuts(BcpsConstraintPool &conPool,
			 double *optLowerSolution);
    /** Helper function for IC*/
    void findLowerLevelSol(double *uselessIneqs, double *lowerLevelSol, const double *sol,
			   bool &isTimeLimReached);

    /** Helper function for IC*/
    bool getAlphaIC(double** extRay, double *uselessIneqs, double* lowerSolution, int numStruct,
		    int numNonBasic, const double* lpSol, std::vector<double> &alphaVec);

    /** Helper function for IC*/
    double solveModelIC(double *uselessIneqs, double *ray, double *rhs, int numNonBasic);

    /** Helper function for watermelon IC **/
    void findLowerLevelSolWatermelonIC(double *uselessIneqs, double *lowerLevelSol,
				       double* lpSol, bool &isTimeLimReached);

    /** Helper function for watermelon IC*/
    bool getAlphaWatermelonIC(double** extRay, double *uselessIneqs, double* lowerSolution,
			      int numStruct, int numNonBasic, double* lpSol,
			      std::vector<double> &alphaVec);

    /** Helper function for hypercube IC*/
    void storeBestSolHypercubeIC(const double* lpSol, double optLowerObj, bool &isTimeLimReached);

    /** Helper function for hypercube IC*/
    void getAlphaHypercubeIC(double** extRay, int numStruct, int numNonBasic,
			     std::vector<double> &alphaVec);

    /** Helper function for Tender IC*/
    void storeBestSolTenderIC(const double* lpSol,
			      double optLowerObj);

    /** Helper function for Tender IC*/
    void getAlphaTenderIC(double** extRay, int numNonBasic,
			  std::vector<double> &alphaVec);
		     
   /** Add bound cuts for general problems **/
   int boundCuts(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables **/
   int incObjCut(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables (current sol)**/
   int incObjCutCurrent(BcpsConstraintPool &conPool);

   /** Add disjunctive cuts for binary upper-level variables (maximal sol) **/
   int incObjCutMaximal(BcpsConstraintPool &conPool);

   int generalWeakIncObjCutCurrent(BcpsConstraintPool &conPool);

   double findBigMIncObjCut();

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

   /** Store the matrices A2, G2 and lower-level coeffs (all constraints are 'L') **/
   void getLowerMatrices(bool getLowerConstCoefMatrix, bool getA2Matrix, bool getG2Matrix);
   
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
