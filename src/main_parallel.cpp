/////////////////////////////////////////////////////////////////
// 24/04/2015 Generate child subproblems of the root node
// 10/10/2014 Read a miplib 0-1 program and builds the IBLP MIP2
// in MODE_IBLP == 1
//      USAGE: iblp.exe <-objvalue>  --- objvalue is the incumbent value in maximization form
//      print probname_iblp.lp
//      solve iblp
//      print to resfile.txt the coversize, cover status (optimal/heuristic) and the list of vars in the cover
//
//      write mps files for all child subproblems of the root node: <probname_rootchild_varname.mps>
//                      varname is the variable index fixed to one
//                      all variables indexing previously printed problems are fixed to zero
//
//      write sol files for all subproblems
//
// parameters: MYNODELIM/MYTILIM is the node/time limit for cover search
//
// in MODE_IBLP == 2
//      USAGE: iblp.exe <-objvalue> [filesol.sol]
//      read a MIP problem inmps format
//      read a solution in .sol format if present in the input
//      solve the problem 
//      print to resfile_MIP  
/////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  

#include <cplex.h>
//#include <cplexx.h>

//#define EXACT
#ifdef EXACT
#include "CoinError.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

#include "MibSModel.h"
#include "MibSSolution.h"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif
#endif

#define M 10E6;
#define EPSILON 0.00000001
#define MYNODELIM 10 // node limit for cover search
#define MYTILIM 3600 // time limit for cover search
#define MYTILIM_SUB 10 // time limit for subMIPs
#define MODE_IBLP 1 // 1 build and solve IBLP - 2 read MIP load sol and solve

static int CPXPUBLIC
cutcallback (CPXCENVptr env,
		 void       *cbdata,
		 int        wherefrom,
		 void       *cbhandle,
		 int        *useraction_p)
{
   int status = 0;
   int nodecount;
   CPXCLPptr lp;
   
   status = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                CPX_CALLBACK_INFO_NODE_COUNT,
                                &nodecount);
   if (nodecount <= 1){
      status = CPXgetcallbacklp (env, cbdata, wherefrom, &lp);
      status =  CPXwriteprob (env, lp, "root.lp", "LP");
   }
}

int  main (int argc, char *argv[])
{

   clock_t start, end;
   double cputime_iblp, cputime_MIP;
     
   char         *inputfilename, *iblpname, *resfilename;
   char         *buffer, *instancefile, *probdir, *solfilename, *solfile;
   CPXENVptr     env = NULL;
   CPXLPptr      lp = NULL, dualp = NULL;
   int           status = 0;
   int           primal_numrows, primal_numcols, dual_numrows, dual_numcols;
   int           numnz;
   //double *x;
   double incumbentval, objval;
   bool loadsol = false;
   int lpstat;
   int GENERATE_INCUMBENT = 0;
     
   iblpname = new char[100];
   resfilename = new char[100];
   buffer = new char[100];
   instancefile = new char[100];
   solfile = new char[100];
     
   /* Check the command line arguments */
     
   switch (argc) {
    case 2:
      inputfilename = argv [1];
      GENERATE_INCUMBENT = 1;
      break;
      
    case 3:
      inputfilename = argv [1];
      incumbentval = atof(argv [2]);
      break;
	
    case 4:
      inputfilename = argv [1];
      incumbentval = atof(argv [2]);
      solfilename = argv [3];
      loadsol = true;
	
      break;
	
    default:
      std::cout << "Usage: <IBLP> <inputfilename> <-objvalue>" 
		<< std::endl;
      exit (1000);
      break;
   }

   probdir = "";

   if (probdir != ""){
      strcpy (instancefile, probdir);
   }
   strcat (instancefile, inputfilename);
     
   if (loadsol){
      strcpy (solfile, probdir);
      strcat (solfile, solfilename);
   }
     
   strcpy (buffer, inputfilename);
   buffer [strlen (buffer) - 4] = '\0';
   strcat (buffer, "_iblp.lp");
   strcpy (iblpname, buffer);
     
   // results file
   std::ofstream resfile_MIP;
   if (MODE_IBLP == 2){
      resfile_MIP.open ("results_MIP.txt", std::ofstream::app);
      resfile_MIP << "Problem" << '\t' << inputfilename << std::endl;
      resfile_MIP << "Time limit" << '\t' << MYTILIM << std::endl;
      resfile_MIP << std::endl;
   }
   
#ifdef COUT_TO_FILE
   std::ofstream resfile;
   if (MODE_IBLP == 1){
      // If we want to put results in a file
      strcpy (buffer, inputfilename);
      buffer [strlen (buffer) - 4] = '\0';
      strcat (buffer, "_results.txt");
      strcpy (resfilename, buffer);

      resfile.open (resfilename, std::ofstream::app);
      std::cout.rdbuf(resfile.rdbuf());
   }
#endif
   
   /* Initialize the CPLEX environment */
     
   env = CPXopenCPLEX (&status);
     
   if ( env == NULL ) {
      char  errmsg[1024];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      exit(0);
   }
     
   status = CPXsetusercutcallbackfunc (env, cutcallback, NULL);
   //status = CPXsetintparam (env, CPXPARAM_Preprocessing_Linear, 0);
   //status = CPXsetintparam (env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);

   /* Turn on output to the screen */
     
   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) {
      fprintf (stderr,
	       "Failure to turn on screen indicator, error %d.\n", status);
      exit(0);
   }
     
   // time and node limits valid for all MODE_IBLP values, that is, all usage modes
     
   // status = CPXsetintparam(env, CPX_PARAM_NODELIM, MYNODELIM);
   status = CPXsetintparam(env, CPX_PARAM_CLOCKTYPE ,  1); 
   status = CPXsetdblparam(env, CPX_PARAM_TILIM, MYTILIM);
     
   /* Create the problem, using the filename as the problem name */
     
   lp = CPXcreateprob (env, &status, instancefile);
     
   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      exit(0);
   }

   /* Now read the file, and copy the data into the created lp */
   //status = CPXsetintparam(env, CPX_PARAM_ADVIND ,  1); 
   status = CPXreadcopyprob (env, lp, instancefile, NULL);
   CPXLPptr clone = CPXcloneprob (env, lp, &status);
   if ( status ) {
      fprintf (stderr, "Failed to clone the problem.\n");
      exit(0);
   }

   if (GENERATE_INCUMBENT == 1){
      status = CPXsetintparam(env, CPX_PARAM_NODELIM, 1);
      status = CPXmipopt(env, clone);
      
      lpstat = CPXgetstat(env, clone);
      if (lpstat == CPXMIP_NODE_LIM_INFEAS){
	 status = CPXsetintparam(env, CPX_PARAM_NODELIM, 10000000000);
	 status = CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
	 status = CPXmipopt(env, clone);
	 lpstat = CPXgetstat(env, clone);
	 status = CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1000000000000);
      }

      if(lpstat == CPXMIP_TIME_LIM_FEAS || lpstat == CPXMIP_SOL_LIM ||
	 lpstat == CPXMIP_OPTIMAL || lpstat == CPXMIP_OPTIMAL_TOL ||
	 lpstat == CPXMIP_NODE_LIM_FEAS){
	 CPXgetobjval(env, clone, &incumbentval);
	 incumbentval *= -1;
      }else{
	 exit(0);
      }
      status = CPXsetintparam(env, CPX_PARAM_NODELIM, 10000000000);
   }

   status = CPXreadcopyprob (env, lp, "root.lp", NULL);

   status = CPXsetusercutcallbackfunc (env, NULL, NULL);
   
   std::cout << "Incumbent value: " << incumbentval << std::endl;

#ifdef EXACT //EXACT Cover

   int i(0), j(0), rc(-1);
   int numInterdictNZ(0);
   
   /** Use MibS to find a cover of minimum size **/

   /** Set up lp solver **/
   OsiClpSolverInterface lpSolver;
   lpSolver.getModelPtr()->setDualBound(1.0e10);
   lpSolver.messageHandler()->setLogLevel(0);
   
   /** Create MibS model **/
   MibSModel model;
   model.setSolver(&lpSolver);
   
   /** read in instance data **/
   CoinMpsIO *mps = new CoinMpsIO;
   rc = mps->readMps(instancefile, "");

   int numCols = mps->getNumCols(); 
   int numRows = mps->getNumRows();
   
   const double *mpsObj =  mps->getObjCoefficients();

   for(i = 0; i < numCols; i++){
      if((mpsObj[i] > model.getTolerance()) || (mpsObj[i] < - model.getTolerance())){
	 numInterdictNZ++;
      }
   }
      
   int numTotalCols = 2 * numCols + 1;
   int numTotalRows = numRows + numCols + 1;
      
#if 0
   //Not sure about this for now
   int structRows(numRows+1);
   model.structRowInd_ = new int[structRows];
   CoinIotaN(structRowInd_, structRows, 0);
   structRowNum_ = structRows;
#endif
      
   double objSense(0.0);
   
   //------------------------------------------------------
   // Set bounds on variables and constraints
   //------------------------------------------------------
   
   double *varLB = new double [numTotalCols];
   double *varUB = new double [numTotalCols];
   
   double *conLB = new double [numTotalRows];
   double *conUB = new double [numTotalRows];
   
   CoinDisjointCopyN(mps->getColLower(), numCols, varLB + numCols);
   CoinDisjointCopyN(mps->getColUpper(), numCols, varUB + numCols);
   
   CoinFillN(varLB, numCols, 0.0); 
   CoinFillN(varUB, numCols, 1.0); 

   /* This is the auxiliary variable we need for adding interdiction cuts */
   CoinFillN(varLB + 2 * numCols, 1, 0.0); 
   CoinFillN(varUB + 2 * numCols, 1, 1.0); 
   
   /* Bound constraint row */
   conLB[0] = -incumbentval;
   conUB[0] = mps->getInfinity();

   /* Lower Rows */
   CoinDisjointCopyN(mps->getRowLower(), numRows, conLB+1);
   CoinDisjointCopyN(mps->getRowUpper(), numRows, conUB+1);
   
   /* Add VUB rows */
   CoinFillN(conLB + (numTotalRows - numCols),
	     numCols, - 1 * mps->getInfinity());
   CoinDisjointCopyN(mps->getColUpper(), numCols, 
		     conUB + (numTotalRows - numCols));
      
   /* Construct objective (minimize size of cover) */ 
   double *objCoef = new double [numTotalCols];
   CoinFillN(objCoef, numCols, 1.0);
   CoinFillN(objCoef+numCols, numTotalCols - numCols, 0.0);
   
   //------------------------------------------------------
   // Set colType_
   //------------------------------------------------------
   
   char *colType = new char [numTotalCols];   
   
   for(j = 0; j < numCols; ++j) {
      colType[j] = 'B';
   }
      
   for(j = 0; j < numCols; ++j) {
      if (mps->isContinuous(j)) {
	 colType[j + numCols] = 'C';
      }
      else {
	 if (varLB[j] == 0 && varUB[j] == 1.0) {
	    colType[j + numCols] = 'B';
	 }
	 else {
	    colType[j = numCols] = 'I';
	 }
      }
   }
   
   //------------------------------------------------------
   // Build coefficient matrix
   //------------------------------------------------------
   
   CoinPackedMatrix matrix = *(mps->getMatrixByCol());
   matrix.reverseOrdering();
   const double * matElements = matrix.getElements();
   const int * matIndices = matrix.getIndices();
   const int * matStarts = matrix.getVectorStarts();
      
   CoinPackedMatrix *newMatrix = new CoinPackedMatrix(false, 0, 0);
   newMatrix->setDimensions(0, numTotalCols);
   int row_start(0), row_end(0), tmp(0), index(0);
      
   /* Add bound constraint row */

   {
      CoinPackedVector row;
      for (j = numCols; j < 2*numCols; j++){
	 row.insert(j, mpsObj[j-numCols]);
      }
      newMatrix->appendRow(row);
   }

   /* lower-level rows */
      
   for(i = 0; i < numRows; i++){
      CoinPackedVector row;
      row_start = matStarts[i];
      row_end = row_start + matrix.getVectorSize(i);
      for(j = row_start; j < row_end; j++){
	 index = matIndices[j] + numCols;
	 row.insert(index, matElements[j]);
      }
      newMatrix->appendRow(row);
   }
      
   /* Add VUB constraints */
   
   for(i = 0; i < numCols; i++){
      CoinPackedVector row;
      row.insert(i, mps->getColUpper()[i]);
      row.insert(i + numCols, 1.0);
      newMatrix->appendRow(row);
   }
   
   newMatrix->reverseOrdering();

   //------------------------------------------------------
   // Construct lower problem data
   //------------------------------------------------------
   
   int upperColNum = numCols;
   int upperRowNum = 1;
   int lowerColNum = numCols;
   int lowerRowNum = numTotalRows-upperRowNum;
   int structRowNum = numTotalRows-numCols;
   int *upperColInd = new int[upperColNum];
   int *upperRowInd = new int[upperRowNum];
   int *lowerColInd = new int[lowerColNum];
   int *lowerRowInd = new int[lowerRowNum];
   int *structRowInd = new int[structRowNum];
   
   CoinIotaN(upperColInd, upperColNum, 0);
   CoinIotaN(upperRowInd, upperRowNum, 0);
   CoinIotaN(lowerColInd, lowerColNum, upperColNum);
   CoinIotaN(lowerRowInd, lowerRowNum, upperRowNum);
   CoinIotaN(structRowInd, structRowNum, 0);
   
   model.loadAuxiliaryData(lowerColNum, lowerRowNum, lowerColInd, lowerRowInd,
			   1.0, mpsObj,
			   upperColNum, upperRowNum, upperColInd, upperRowInd,
			   structRowNum, structRowInd, 0, NULL);

   model.loadProblemData(*newMatrix, varLB, varUB, objCoef, conLB, 
			 conUB, colType, 1, mps->getInfinity());

   model.AlpsPar()->setEntry(AlpsParams::timeLimit, MYTILIM);
   model.MibSPar()->setEntry(MibSParams::bilevelProblemType, 1);
   //model.MibSPar()->setEntry(MibSParams::cutStrategy, 1);
   model.MibSPar()->setEntry(MibSParams::bilevelCutTypes, 1);
   model.MibSPar()->setEntry(MibSParams::useBendersCut, true);
   //model.MibSPar()->setEntry(MibSParams::useLowerObjHeuristic, false);
   //model.MibSPar()->setEntry(MibSParams::useObjCutHeuristic, false);
   //model.MibSPar()->setEntry(MibSParams::useWSHeuristic, false);
   //model.MibSPar()->setEntry(MibSParams::useGreedyHeuristic, false);

   int argc = 1;
   char** argv = new char* [1];
   argv[0] = "mibs";

#ifdef  COIN_HAS_MPI
   AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
   AlpsKnowledgeBrokerSerial broker(argc, argv, model);
#endif
   
   broker.search(&model);

   MibSSolution *solution = dynamic_cast<MibSSolution* >
      (broker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
   
   double *y = new double[upperColNum];
   
   for (j = 0; j < upperColNum; j++){
      y[j] = floor(solution->getValues()[j] + 0.5);
   }

   delete mps;

   primal_numcols = numCols;
   primal_numrows = numRows;

#else //INEXACT cover

   /** Use CPLEX to find a heuristic cover **/
     
   status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);
   if ( status ) {
      fprintf (stderr,
	       "Failed to set the optimization method, error %d.\n", status);
      exit(0);
   }
     
     
   // solve the MIP problem loading the initial solution
   if (MODE_IBLP == 2){
      if (loadsol){
	 status = CPXreadcopysol(env, lp, solfile); 
	   
	 // check if solution has been loaded
	 // does not work!!
	 //  status = CPXgetbestobjval(env, lp, &loadobjval); 
	 //  resfile_MIP << "initial objval: " << loadobjval << std::endl;
	 if (status == 0)
	    resfile_MIP << "MIP start solution loaded successfully " << std::endl;
      }
	
      start = clock();
      status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
      status = CPXmipopt (env, lp);
      end = clock();
      if ( status ) {
	 fprintf (stderr, "Failed to optimize MIP.\n");
	 exit(0);
      }
	
      status = CPXgetobjval (env, lp, &objval);
      if ( status ) {
	 fprintf (stderr, "Failed to obtain objective value.\n");
	 exit(0);
      }
      cputime_MIP = ((double) (end - start)) / CLOCKS_PER_SEC;
	
      lpstat = CPXgetstat(env, lp);
	
      if(lpstat == CPXMIP_TIME_LIM_FEAS || lpstat == CPXMIP_OPTIMAL || lpstat == CPXMIP_OPTIMAL_TOL){
	 resfile_MIP << "Obj value " << objval << std::endl;
	 resfile_MIP << "Number of evaluated subproblems " << CPXgetnodecnt(env, lp ) << std::endl;
	 resfile_MIP << "CPU time for MIP " << cputime_MIP << std::endl;
	 resfile_MIP << std::endl;
      }
      else{
	 resfile_MIP << "Premature exit" << std::endl;
      }
      exit(0);
   }
     
   primal_numcols = CPXgetnumcols (env, lp);
   primal_numrows = CPXgetnumrows (env, lp);
     
   /* Retrieve solution vector */
   /*
     x = (double *) malloc (primal_numcols*sizeof(double));
     if ( x == NULL ) {
     fprintf (stderr, "No memory for solution.\n");
     exit(0);
     }
       
     status = CPXgetx (env, lp, x, 0, cur_numcols-1);
     if ( status ) {
     fprintf (stderr, "Failed to obtain primal solution.\n");
     exit(0);
     }
   */
   // writes the dual problem of the LP relaxation
   status = CPXchgprobtype (env, lp, CPXPROB_LP);
   double objshift_p;
   status =  CPXdualwrite (env, lp, "dualprob.dua", &objshift_p);
     
   // creates and reads the dual problem
   dualp = CPXcreateprob (env, &status, iblpname);
     
   if ( dualp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      exit(0);
   }
     
   status = CPXreadcopyprob (env, dualp, "dualprob.dua", NULL);
   if ( status ) {
      fprintf (stderr, "Failed to read and copy the dual %problem data.\n");
      exit(0);
   }
     
   dual_numcols = CPXgetnumcols (env, dualp);
   dual_numrows = CPXgetnumrows (env, dualp);
   numnz        = CPXgetnumnz(env, dualp);

   /* stuff to get the variables in the cover */
   double *y;
   y = (double *) malloc ((primal_numcols)*sizeof(double));
   if ( y == NULL ) {
      fprintf (stderr, "No memory for iblp solution.\n");
      exit(0);
   }
     
   // build IBLP by modifying the dual
     
   // translate the dual obj into the UB constraint
   int rmatbeg_one = 0;
   char sense = 'L';
   double rhs = incumbentval;
   double *duaobj; // gather the original dual ojective
     
   // the size of this includes possible zeros 
   int *onerowmatind;
   onerowmatind = (int *) malloc (dual_numcols*sizeof(int));
   if ( onerowmatind == NULL ) {
      fprintf (stderr, "No memory for the new constraint.\n");
      exit(0);
   }
     
   // the size of this includes possible zeros 
   double *onerowmatval;
   onerowmatval = (double *) malloc (dual_numcols*sizeof(double));
   if ( onerowmatval == NULL ) {
      fprintf (stderr, "No memory for the new constraint.\n");
      exit(0);
   }
     
   // the size of this includes possible zeros
   duaobj = (double *) malloc (dual_numcols*sizeof(double));
   if ( duaobj == NULL ) {
      fprintf (stderr, "No memory for original dual obj values.\n");
      exit(0);
   }
     
   status = CPXgetobj(env, dualp, duaobj, 0, dual_numcols - 1); 
     
   int objduanzcnt = 0;
   for(int j = 0; j < dual_numcols; j++){
      if (fabs(duaobj[j]) > EPSILON){
	 onerowmatind[objduanzcnt] = j;
	 onerowmatval[objduanzcnt] = duaobj[j];
	 objduanzcnt ++;
      }
   }
     
   status = CPXaddrows(env, dualp, 0, 1, objduanzcnt, &rhs, &sense,
		       &rmatbeg_one, 
		       onerowmatind, onerowmatval, NULL,NULL);
     
   // reset obj coefficient of the dual vars
   int *dualindices;
   dualindices = (int *) malloc (dual_numcols*sizeof(int));
   if ( dualindices == NULL ) {
      fprintf (stderr, "No memory for dual indices.\n");
      exit(0);
   }
     
   double *dualobjvalues;
   dualobjvalues = (double *) malloc (dual_numcols*sizeof(double));
   if ( dualobjvalues == NULL ) {
      fprintf (stderr, "No memory for dual obj values.\n");
      exit(0);
   }
     
   for (int k = 0; k < dual_numcols; k++){ 
      dualindices[k] = k;	  
      dualobjvalues[k] = 0.0;
   }
     
   status = CPXchgobj(env, dualp, dual_numcols, dualindices, dualobjvalues);
     
   //Allocate space for y variables
   double *yobj;
   yobj = (double *) malloc (primal_numcols*sizeof(double));
   if ( yobj == NULL ) {
      fprintf (stderr, "No memory for y vars.\n");
      exit(0);
   }
   for (int j = 0; j < primal_numcols; j++) yobj[j] = 1.0;
     
   // introduce indicator constraints

   int *rmatbeg;
   rmatbeg = (int *) malloc ((dual_numrows+1)*sizeof(int));
   if ( rmatbeg == NULL ) {
      fprintf (stderr, "No memory for rmatbeg.\n");
      exit(0);
   }
     
   int *rmatind;
   rmatind = (int *) malloc (numnz*sizeof(int));
   if ( rmatind == NULL ) {
      fprintf (stderr, "No memory for rmatind.\n");
      exit(0);
   }
     
   double *rmatval;
   rmatval = (double *) malloc (numnz*sizeof(double));
   if ( rmatval == NULL ) {
      fprintf (stderr, "No memory for rmatval.\n");
      exit(0);
   }

   double *rhs_vals;
   rhs_vals = (double *) malloc (dual_numrows*sizeof(double));
   if ( rhs_vals == NULL ) {
      fprintf (stderr, "No memory for rhs_vals.\n");
      exit(0);
   }

   char *senses;
   senses = (char *) malloc (dual_numrows*sizeof(char));
   if ( senses == NULL ) {
      fprintf (stderr, "No memory for senses.\n");
      exit(0);
   }

   // Get the matrix in row order
   int space;
   CPXgetrows(env, dualp, &numnz, rmatbeg, rmatind, rmatval, numnz,
	      &space, 0, dual_numrows-1);
   
   CPXgetrhs(env, dualp, rhs_vals, 0, dual_numrows-1);
   
   CPXgetsense(env, dualp, senses, 0, dual_numrows-1);
   
   // Introduce big-M coefficients

   int *cmatbeg;
   cmatbeg = (int *) malloc (primal_numcols*sizeof(int));
   if ( cmatbeg == NULL ) {
      fprintf (stderr, "No memory for cmatbeg.\n");
      exit(0);
   }
     
   int *cmatind;
   cmatind = (int *) malloc (primal_numcols*sizeof(int));
   if ( cmatind == NULL ) {
      fprintf (stderr, "No memory for cmatind.\n");
      exit(0);
   }
     
   double *cmatval;
   cmatval = (double *) malloc (primal_numcols*sizeof(double));
   if ( cmatval == NULL ) {
      fprintf (stderr, "No memory for cmatval.\n");
      exit(0);
   }
     
   double *lb;
   lb = (double *) malloc (primal_numcols*sizeof(double));
   if ( lb == NULL ) {
      fprintf (stderr, "No memory for lb.\n");
      exit(0);
   }
     
   double *ub;
   ub = (double *) malloc (primal_numcols*sizeof(double));
   if ( ub == NULL ) {
      fprintf (stderr, "No memory for ub.\n");
      exit(0);
   }
     
   char **colname;
   colname = (char **) malloc (primal_numcols*sizeof(char*));
   for (int i = 0; i < primal_numcols; i++)
      colname [i] = (char *) malloc (primal_numrows*sizeof(char));
     
     
   for (int j = 0; j < primal_numcols; j++){
	
      cmatbeg[j] = j;
      cmatind[j] = j;
      cmatval[j] = M;
      lb[j] = 0.0;
      ub[j] = 1.0;
	
      sprintf (colname [j], "y%d", j+1);
   }
     
   status = CPXaddcols(env, dualp, primal_numcols, primal_numcols, yobj, 
		       cmatbeg, cmatind, cmatval, lb, ub, colname);
     
   // Add the indicator constraints

   char **indconstr_name;
   indconstr_name = (char **) malloc (dual_numrows*sizeof(char*));
   for (int i = 0; i < dual_numrows; i++){
      indconstr_name[i] = (char *) malloc (dual_numrows*sizeof(char));
      sprintf (indconstr_name[i], "ind%d", i+1);
   }

   rmatbeg[dual_numrows] = numnz;
   for (int j = 0; j < dual_numrows; j++){
      CPXaddindconstr(env, dualp, dual_numcols+j, 1, rmatbeg[j+1]-rmatbeg[j],
		      rhs_vals[j], senses[j], rmatind+rmatbeg[j],
		      rmatval+rmatbeg[j], indconstr_name[j]);
   }
   
   // convert the problem to a MILP
     
   int *indices;
   indices = (int *) malloc (primal_numcols*sizeof(int));
   if ( indices == NULL ) {
      fprintf (stderr, "No memory for indices.\n");
      exit(0);
   }

   char *ctype;
   ctype = (char *) malloc (primal_numcols*sizeof(char));
   if ( ctype == NULL ) {
      fprintf (stderr, "No memory for ctype.\n");
      exit(0);
   }

   // Assume that the new cols are indexed after the original ones
   for (int j = 0; j < primal_numcols; j++){
      indices[j] = dual_numcols + j;
      ctype[j] = 'B';
   }

   status = CPXchgctype(env, dualp, primal_numcols, indices, ctype);

   status = CPXchgprobtype(env, dualp, CPXPROB_MILP);

   status = CPXchgprobname (env, dualp, iblpname);

   status =  CPXwriteprob (env, dualp, iblpname, "LP");

   // set integer tolerance to 10^{-9}
   //status = CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_Integrality,  1e-09); 
   //	  status = CPXsetintparam (env, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);

   // solve IBLP with TLIM
   start = clock();
   status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
   status = CPXmipopt (env, dualp);
   end = clock();
   if ( status ) {
      fprintf (stderr, "Failed to optimize LP.\n");
      exit(0);
   }

   cputime_iblp = ((double) (end - start)) / CLOCKS_PER_SEC;
   std::cout << "CPU time for iblp " << cputime_iblp << std::endl;
   std::cout << "Number of evaluated subproblems " << CPXgetnodecnt(env, dualp ) << std::endl;

   lpstat = CPXgetstat(env, dualp);

   double coversize;
   //	  if(lpstat == CPXMIP_NODE_LIM_FEAS || lpstat == CPXMIP_OPTIMAL || lpstat == CPXMIP_OPTIMAL_TOL){
   if(lpstat == CPXMIP_TIME_LIM_FEAS || lpstat == CPXMIP_OPTIMAL || lpstat == CPXMIP_OPTIMAL_TOL){

      /* Retrieve last part of solution vector, that is, the variables in the cover */
      status = CPXgetobjval (env, dualp, &coversize);
      if ( status ) {
	 fprintf (stderr, "Failed to obtain coversize value.\n");
	 exit(0);
      }

      status = CPXgetx (env, dualp, y, dual_numcols, dual_numcols + primal_numcols -1);
      if ( status ) {
	 fprintf (stderr, "Failed to obtain iblp solution.\n");
	 exit(0);
      }

      // wites the solution to std::cout
      std::cout << "cover size: " << coversize << std::endl;

      if (coversize == 0) std::cout << "incumbent value is certitied OPTIMAL!" << std::endl;
      else if(lpstat == CPXMIP_OPTIMAL || lpstat == CPXMIP_OPTIMAL_TOL) std::cout << "cover is optimal" << std::endl;
      else std::cout << "cover is not proven to be optimal" << std::endl;

      std::cout << "variables in the cover: " << std::endl;
      for (int j = 0; j < primal_numcols; j++){
	 if (y[j] > 0.5)
	    std::cout << colname[j] << '\t';
      }
      std::cout << std::endl;
   }
   else{
      std::cout << "cover NOT found" << std::endl;
      std::cout << std::endl;
   }

   free (yobj);
   free (cmatbeg);
   free (cmatind);
   free (cmatval);
   free (lb);
   free (ub);
   for (int i = 0; i < primal_numcols; i++) free (colname [i]);
   free (colname);

   free(dualindices);
   free(dualobjvalues);
   free(indices);
   free(ctype);
   free(duaobj);

   free(onerowmatind);
   free(onerowmatval);

#endif

   // At this point we can generate the subproblems of the root node
   // LOOP: at iteration j
   //		- create lp_j
   //		- print	 lp_j

   double	cputime_sub;
   char		*subprobname, *solname;
   CPXLPptr	lpsub = NULL;
   subprobname = new char[100];
   solname = new char[100];
   int index2one;
   int *indices2zero;
   int numfixedzero = 0;

   char lowerbound = 'L';
   double valueone = 1.0;
   char upperbound = 'U';
   double valuezero = 0.0;

   indices2zero = new int [primal_numcols];

   double best_val = -incumbentval, cur_val;
   for (int j = 0, count = 0; j < primal_numcols && count <= 100; j++){
      if (y[j] > 0.5){
	 count++;
	 
	 index2one = j;

	 // create child with var index fixed to 1
	 strcpy (subprobname, inputfilename);
	 subprobname [strlen (subprobname) - 4] = '\0';
	 strcat (subprobname, "_rootchild_");

	 sprintf(buffer,"%d",j+1);
	 strcat (subprobname, buffer);
	 strcpy (solname, subprobname);
	 lpsub = CPXcreateprob (env, &status, subprobname);
			
	 strcat (subprobname, ".mps");
	 strcat (solname, ".sol");

	 // inizializza il problema copiando il problema root e fissando a 1 la var j
	 status = CPXreadcopyprob (env, lpsub, instancefile, NULL);
	 if ( status ) {
	    fprintf (stderr, "Failed to read and copy the problem data.\n");
	    exit(0);
	 }

	 // stupid cplex copies also the name
	 status = CPXchgprobname (env, lpsub, subprobname);

	 // fix the current index to one
	 status = CPXchgbds(env, lpsub, 1, &index2one, &lowerbound, &valueone);

	 // fix the previous indices to zero
	 for (int r = 0; r < numfixedzero; r++){
	    status = CPXchgbds(env, lpsub, 1, &indices2zero[r], &upperbound, &valuezero);
	 }

	 // print subproblem
	 status =  CPXwriteprob (env, lpsub, subprobname, "MPS");

	 // solve subproblem with TLIM
	 start = clock();
	 status = CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_UpperCutoff,
				 best_val - 10e-5);
	 status = CPXsetdblparam(env, CPX_PARAM_TILIM, MYTILIM_SUB);
	 status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	 status = CPXmipopt (env, lpsub);
	 lpstat = CPXgetstat (env, lpsub);
	 if (lpstat == CPXMIP_TIME_LIM_FEAS || lpstat == CPXMIP_OPTIMAL ||
	     lpstat == CPXMIP_OPTIMAL_TOL){
	    status = CPXgetobjval(env, lpsub, &cur_val);
	    if (cur_val < best_val){
	       best_val = cur_val;
#if 0
	       status = CPXsolwrite(env, lpsub, solname);
	       if ( status ) {
		  fprintf (stderr, "Failed to write subproblem solution.\n");
		  exit(0);
	       }
#endif
	    }
	 }
	 
	 end = clock();
	 if ( status ) {
	    fprintf (stderr, "Failed to optimize subproblem.\n");
	    exit(0);
	 }

	 cputime_sub = ((double) (end - start)) / CLOCKS_PER_SEC;

	 std::cout << std::endl;
	 std::cout << "Stats for subproblem " << subprobname << std::endl;
	 std::cout << "CPU time " << cputime_sub << std::endl;
	 std::cout << "Status   " << lpstat << std::endl;
	 if (lpstat == CPXMIP_TIME_LIM_FEAS || lpstat == CPXMIP_OPTIMAL ||
	     lpstat == CPXMIP_OPTIMAL_TOL){
	    std::cout << "Value  " << cur_val << std::endl;
	 }
	 std::cout << "Number of evaluated subproblems " << CPXgetnodecnt(env, lpsub ) << std::endl;

	 // include current index in the fix2zero list
	 indices2zero[numfixedzero] = index2one;
	 numfixedzero ++;

      }
   }

   std::cout << std::endl;
   std::cout << "Best solution found: " << best_val << std::endl;

   //free (x);
   free (y);

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
	 fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */
   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      if ( status ) {
	 char	errmsg[1024];
	 fprintf (stderr, "Could not close CPLEX environment.\n");
	 CPXgeterrorstring (env, status, errmsg);
	 fprintf (stderr, "%s", errmsg);
      }
   }

   return (status);
	  
}  /* END main */
