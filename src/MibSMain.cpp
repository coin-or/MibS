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

#include <iostream>

#include "CoinError.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

#include "MibSConfig.hpp"
#include "MibSModel.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

//#############################################################################
//#############################################################################

int main(int argc, char* argv[])
{

   std::cout 
      << "\n========================================================"
      << "\n========================================================"
      <<   "\nWelcome to MibS (Mixed Integer Bilevel Solver)"
      <<   "\nCopyright (C) 2007-2019 Lehigh University, Scott DeNegre, Ted Ralphs"
      <<   "\nand Sahar Tahernejad."
      <<   "\nAll Rights Reserved."
      <<   "\nThis software is licensed under the Eclipse Public License. Please see"
      <<   "\naccompanying file for terms."
      <<   "\nVersion: " << MIBS_VERSION
      <<   "\nBuild Date: " << __DATE__
      << "\n========================================================"
      << "\n========================================================"
      << "\n";

    try{
       
      /** Set up lp solver **/
      OsiClpSolverInterface lpSolver;
      lpSolver.getModelPtr()->setDualBound(1.0e10);
      lpSolver.messageHandler()->setLogLevel(0);
      
      /** Create MibS model **/
      MibSModel model;
      model.setSolver(&lpSolver);


#ifdef  COIN_HAS_MPI
        AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
        if(argc == 1){
           std::cerr
              << "================================================" << std::endl
              << "Usage:" << std::endl
              << "1)./mibs -Alps_instance file.mps" << std::endl
              << "       -MibS_auxiliaryInfoFile aux_file.aux" << std::endl
              << "       --BlockFile /FilePath/ABC.block" << std::endl
              << "2)./mibs -param mibs.par" << std::endl
              << "================================================" << std::endl
              << std::endl;
           return 0;
	}
        AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif

	broker.search(&model);
	broker.printBestSolution();

    }
    catch(CoinError& er) {
	std::cerr << "ERROR:" << er.message() << std::endl
		  << " from function " << er.methodName() << std::endl
		  << " from class " << er.className() << std::endl;
    }
    catch(...) {
	std::cerr << "Something went wrong!" << std::endl;
    }

    return 0;
}
