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
 *===========================================================================*/

#ifndef MibSSolTypes_h_
#define MibSSolTypes_h_

//#############################################################################
class mcSol{

 private:

  /** pair of upper- and lower-level objectives from MC algo **/
  std::pair<double, double> objPair_;

  /** the column solution that corresponds to the obj pair **/
  const double * colSol_;

 public:
  
  mcSol()
    : objPair_(0.0, 0.0), colSol_(0) {
  }

 mcSol(std::pair<double, double> objpair, const double * colsol)
   : objPair_(objpair), colSol_(colsol) {

  }

  ~mcSol()
    {
      //if(colSol_) delete [] colSol_;
    }

  std::pair<double, double> getObjPair() {return objPair_;}
  const double * getColumnSol() {return colSol_;}

  /*
  mcSol * operator=(mcSol * sol)
    {
      mcSol tmp;
      tmp.objPair_ = sol->getObjPair();
      tmp.colSol_ = sol>getColumnSol();

    }
  */
};

//#############################################################################
class bfSol{

 private:

  /** upper-level objective value **/
  double objVal_;

  /** the column solution that corresponds to the obj value **/
  const double * colSol_;

 public:
  
 bfSol()
   : objVal_(0.0), colSol_(0) {
  }

 bfSol(double objval, const double * colsol)
   : objVal_(objval), colSol_(colsol) {
    
  }
  
  ~bfSol()
    {
      //if(colSol_) delete [] colSol_;
    }

  double getObjVal() {return objVal_;}
  const double * getColumnSol() {return colSol_;}

  void setObjVal(double val) {objVal_ = val;}
  void setColumnSol(double * sol) {colSol_ = sol;}


};

//#############################################################################

class SolPair{

 private:

  std::pair<double, double> pair1_;
  std::pair<double, double> pair2_;

 public:

  SolPair()
    : pair1_(0.0, 0.0), pair2_(0.0, 0.0) {
  }

  SolPair(std::pair<double, double> pair1, std::pair<double, double> pair2)
   : pair1_(pair1), pair2_(pair2) {
  }

  ~SolPair() {}

  std::pair<double, double> getFirst() {return pair1_;}
  std::pair<double, double> getSecond() {return pair2_;}

  void setPairs(std::pair<double, double> first, 
		std::pair<double, double> second) 
  {
    pair1_ = first;
    pair2_ = second;
  }
    
};

#endif
