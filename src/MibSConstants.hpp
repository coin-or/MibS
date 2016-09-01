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

#ifndef MibSConstants_h_
#define MibSConstants_h_

#define BIGCONSTANT 10000

#define PARAM_NOTSET  -1
#define PARAM_OFF      0
#define PARAM_ON       1

/*---------  which_active_con_method choices --------------------------------*/
#define SIMPLE  0
#define BASIS  1

/*---------  bilevel_problem_type choices --------------------------------*/
#define GENERAL  0
#define INTERDICT  1

/*---------  bilevel_cut_typeschoices --------------------------------*/
#define GENERALONLY  0
#define GENERALINTERDICT  1
#define GENERALBINARYUL  2

/*---------  cut_strategy choices ----------------------------------------*/
#define BRANCHONLY 0
#define CUTONLY 1
#define BRANCHANDCUT 2

/*---------  obj_bound_strategy choices -------------------------------------*/
#define LPBOUND 0
#define KNAPBOUND 1

#endif
