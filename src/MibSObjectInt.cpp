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

#include "MibSObjectInt.hpp"
#include "MibSModel.hpp"
#include "MibSBranchObjectInt.hpp"


//#############################################################################
MibSObjectInt::~MibSObjectInt()
{

}

//#############################################################################

// Creates a branching object from this infeasible object.
BcpsBranchObject * 
MibSObjectInt::createBranchObject(BcpsModel *m, int direction) const
{
    BlisModel *model = dynamic_cast<BlisModel* >(m);
    OsiSolverInterface * solver = model->solver();
    
    //double integerTolerance 
    //= model->BlisPar()->entry(BlisParams::integerTol);
    
    const double * solution = solver->getColSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    
    double value = solution[columnIndex_];
    //std::cout << "COL"<< columnIndex_ << ": x = " << value << std::endl;
    
    // Force value in bounds.
    value = CoinMax(value, lower[columnIndex_]);
    value = CoinMin(value, upper[columnIndex_]);
    
    //FIXME::CHANGE ARGUMENTS
    return new MibSBranchObjectInt(model, objectIndex_, direction, value);
    //return new BlisBranchObjectInt(model, objectIndex_, direction, value);
}
