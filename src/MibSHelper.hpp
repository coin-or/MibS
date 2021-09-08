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

#ifndef MibSHelper_h_
#define MibSHelper_h_

//#############################################################################

struct LINKING_SOLUTION{
    int tag;
    std::vector<double> lowerObjValue;
    std::vector<double> brLowerObjValue; // YX: container for SL-MILP with gap results; vector for stochastic case?
    double UBObjValue;
    std::vector<double> lowerSolution;
    std::vector<double> brLowerSolution; // YX: container for SL-MILP with gap results;
    std::vector<double> UBSolution;
};

#endif
