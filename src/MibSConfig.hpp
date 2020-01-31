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

/* Include file for the configuration of MibS.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these file and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h files.
 *
 */

#ifndef MibSConfig_hpp_
#define MibSConfig_hpp_

#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__MINGW64__)
/* Different function call in Windows */
#define SRANDOM(seed) srand(seed)
#define RANDOM() rand()
#else
#define SRANDOM(seed) srandom(seed)
#define RANDOM() random()
#endif

#ifdef HAVE_CONFIG_H
#ifdef MIBS_BUILD
#include "config.h"
#else
#include "config_mibs.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef MIBS_BUILD
#include "config_default.h"
#else
#include "config_mibs_default.h"
#endif

#endif /* HAVE_CONFIG_H */

#endif /*__MIBSCONFIG_H__*/

//supporting c++11
#if __cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1900)
#define COIN_HAS_C11
#endif
