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

/* Define to 1 if Alps is available. */
#define MIBS_HAS_ALPS 1

/* Define to 1 if Bcps is available. */
#define MIBS_HAS_BCPS 1

/* Define to 1 if Blis is available. */
#define MIBS_HAS_BLIS 1

/* Define to 1 if Cbc is available. */
#define MIBS_HAS_CBC 1

/* Define to 1 if Cgl is available. */
#define MIBS_HAS_CGL 1

/* Define to 1 if Clp is available. */
#define MIBS_HAS_CLP 1

/* Define to 1 if CoinUtils is available. */
#define MIBS_HAS_COINUTILS 1

/* Define to 1 if Cplex is available. */
/* #undef MIBS_HAS_CPLEX */

/* Define to 1 if Osi is available. */
#define MIBS_HAS_OSI 1

/* Define to 1 if Sample is available. */
#define MIBS_HAS_SAMPLE 1

/* Define to 1 if SYMPHONY is available. */
#define MIBS_HAS_SYMPHONY 1
