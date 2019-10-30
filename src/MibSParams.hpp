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

#ifndef MibSParams_h
#define MibSParams_h

#include "AlpsKnowledge.h"
#include "AlpsParameterBase.h"

//#############################################################################

//** Parameters used in MibS. */
class MibSParams : public AlpsParameterSet {
 public:
  /** Character parameters. All of these variable are used as booleans
      (true = 1, false = 0). */
  enum boolParams{
     useValFuncCut,
     useBoundCut,
     boundCutOptimal,
     boundCutRelaxUpper,
     useIpBound,
     isBilevelBranchProb,
     warmStartLL,
     doDualFixing,
     useUBDecompose,
     turnOffOtherCuts,
     printProblemInfo,
     allowRemoveCut,
     useNewPureIntCut,
     useProgresHedg,
     endOfBoolParams
  };
  
  /** Integer paramters. */
  enum intParams{
    //     verbosity,
     whichActiveConMethod,
     maxNumActiveCons,
     bilevelProblemType,
     bilevelCutTypes,
     cutStrategy,
     objBoundStrategy,
     blisCutStrategy,
     blisBranchStrategy,
     branchStrategy,
     upperFileFormat,
     maxThreadsLL,
     whichCutsLL,
     usePreprocessor,
     useLowerObjHeuristic,
     useObjCutHeuristic,
     useWSHeuristic,
     useGreedyHeuristic,
     usePureIntegerCut,
     useNoGoodCut,
     useGeneralNoGoodCut,
     useIncObjCut,
     useBendersCut,
     bendersCutType,
     useIntersectionCut,
     //intersectionCutType,
     useTypeIC,
     useTypeWatermelon,
     useTypeHypercubeIC,
     useTypeTenderIC,
     useTypeHybridIC,
     bilevelFreeSetTypeIC,
     solveSecondLevelWhenXYVarsInt,
     solveSecondLevelWhenXVarsInt,
     solveSecondLevelWhenLVarsInt,
     solveSecondLevelWhenLVarsFixed,
     computeBestUBWhenXVarsInt,
     computeBestUBWhenLVarsInt,
     computeBestUBWhenLVarsFixed,
     useLinkingSolutionPool,
     newPureIntCutDepthLb,
     newPureIntCutDepthUb,
     boundCutOptimalType,
     boundCutDepthLb,
     boundCutDepthUb,
     boundCutFreq,
     boundCutNodeLim,
     relaxTypeParamBoundCut,
     maxActiveNodes,
     isA2Random,
     isSMPSFormat,
     sampSizeSAA,
     evalSampSizeSAA,
     replNumSAA,
     lbDistB2SAA,
     ubDistB2SAA,
     lbDistA2SAA,
     ubDistA2SAA,
     incDistB2NumerSAA,
     incDistB2DenumSAA,
     incDistA2NumerSAA,
     incDistA2DenumSAA,
     iterationLimitPH,
     nodeLimitPHSubprob,
     endOfIntParams
  };

  /** Double parameters. */
  enum dblParams{
      boundCutTimeLim,
      optimalRelGapLimitPHSubprob,
      endOfDblParams
  };

  /** String parameters. */
  enum strParams{
      strDummy,
      auxiliaryInfoFile,
      auxiliaryTimFile,
      auxiliaryStoFile,
      feasCheckSolver,
      inputFormat,
      stochasticityType,
      endOfStrParams
  };

  /** There are no string array parameters. */
  enum strArrayParams{
      strArrayDummy,
      ///
      endOfStrArrayParams
  };

 public:
  /**@name Constructors. */
  /*@{*/
  /** The default constructor creates a parameter set with from the template
      argument structure. The keyword list is created and the defaults are
      set. */
  MibSParams() :
    AlpsParameterSet(
		     static_cast<int>(endOfBoolParams),
		     static_cast<int>(endOfIntParams),
		     static_cast<int>(endOfDblParams),
		     static_cast<int>(endOfStrParams),
		     static_cast<int>(endOfStrArrayParams)
		     )
    {
      createKeywordList();
      setDefaultEntries();
    }
  /*@}*/

  /** Method for creating the list of keyword looked for in the parameter
      file. */
  virtual void createKeywordList();
  /** Method for setting the default values for the parameters. */
  virtual void setDefaultEntries();
  /*@}*/


 public:
  //====================================================
  /** For user's application: 
   *    Copy following code exactly (till the end of this class) and 
   *    do NOT change anything. 
   *
   *  The reason can not put following functions in base class 
   *  <CODE> AlpsParameterSet </CODE> is:
   *
   *    <CODE> boolParams </CODE> and <CODE> endOfBoolParams </CODE> etc. 
   *    can NOT be declared in base class. They are different types for
   *    each derived classes.
   */
  //====================================================
  
  
  /**@name Query methods     
     The members of the parameter set can be queried for using the 
     overloaded entry() method. Using the example in the class
     documentation the user can get a parameter with the
     "<code>param.entry(USER_par::parameter_name)</code>" expression.
  */
  /*@{*/
  ///
  inline bool entry(const boolParams key) const { return bpar_[key]; }
  ///
  inline int entry(const intParams key) const { return ipar_[key]; }
  ///
  inline double entry(const dblParams key) const { return dpar_[key]; }
  ///
  inline const std::string&
    entry(const strParams key) const { return spar_[key]; }
  ///
  inline const std::vector<std::string>&
    entry(const strArrayParams key) const { return sapar_[key]; }
  /*@}*/

  //----------------------------------------------------
  
  /// char* is true(1) or false(0), not used
  void setEntry(const boolParams key, const char * val) {
    bpar_[key] = atoi(val) ? true : false; }
  /// char is true(1) or false(0), not used
  void setEntry(const boolParams key, const char val) {
    bpar_[key] = val ? true : false; }
  /// This method is the one that ever been used.
  void setEntry(const boolParams key, const bool val) {
    bpar_[key] = val; }
  ///
  void setEntry(const intParams key, const char * val) {
    ipar_[key] = atoi(val); }
  ///
  void setEntry(const intParams key, const int val) {
    ipar_[key] = val; }
  ///
  void setEntry(const dblParams key, const char * val) {
    dpar_[key] = atof(val); }
  ///
  void setEntry(const dblParams key, const double val) {
    dpar_[key] = val; }
  ///
  void setEntry(const strParams key, const char * val) {
    spar_[key] = val; }
  ///
  void setEntry(const strArrayParams key, const char *val) {
    sapar_[key].push_back(val); }
  
  //----------------------------------------------------
  
  /**@name Packing/unpacking methods */
  /*@{*/
  /** Pack the parameter set into buf. */
  void pack(AlpsEncoded& buf) {
    buf.writeRep(bpar_, endOfBoolParams)
      .writeRep(ipar_, endOfIntParams)
      .writeRep(dpar_, endOfDblParams);
    for (int i = 0; i < endOfStrParams; ++i)
      buf.writeRep(spar_[i]);
    for (int i = 0; i < endOfStrArrayParams; ++i) {
      buf.writeRep(sapar_[i].size());
      for (size_t j = 0; j < sapar_[i].size(); ++j)
	buf.writeRep(sapar_[i][j]);
    }
  }
  
  /** Unpack the parameter set from buf. */
  void unpack(AlpsEncoded& buf) {
    int dummy;
    // No need to allocate the arrays, they are of fixed length
    dummy = static_cast<int>(endOfBoolParams);
    buf.readRep(bpar_, dummy, false);
    dummy = static_cast<int>(endOfIntParams);
    buf.readRep(ipar_, dummy, false);
    dummy = static_cast<int>(endOfDblParams);
    buf.readRep(dpar_, dummy, false);
    for (int i = 0; i < endOfStrParams; ++i)
      buf.readRep(spar_[i]);
    for (int i = 0; i < endOfStrArrayParams; ++i) {
      size_t str_size;
      buf.readRep(str_size);
      sapar_[i].reserve(str_size);
      for (size_t j = 0; j < str_size; ++j){
	//	sapar_[i].unchecked_push_back(std::string());
	sapar_[i].push_back(std::string());
	buf.readRep(sapar_[i].back());
      }
    }
  }
  /*@}*/
  
};

#endif
