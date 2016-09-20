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

/* include the COIN-wide system specific configure header */
#include "configall_system.h"

/* include the public project specific macros */
#include "config_mibs_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to 1 if the CoinUtils package is used */
#define COIN_HAS_COINDEPEND 1

/* Define to 1 if the SYMPHONY package is used */
#define COIN_HAS_SYMPHONY 1

/* Define to 1 if the SYMPHONY package is used */
#define COIN_HAS_CPLEX 1
