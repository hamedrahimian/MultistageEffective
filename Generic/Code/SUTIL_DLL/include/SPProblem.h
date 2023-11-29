/*  
 *     SUTIL -- A Stochastic Programming Utility Library
 *
 *     VERSION 0.1
 *
 *     Authors:   Joe Czyzyk
 *                Northwestern University
 *
 *		  Jeff Linderoth and Jierui Shen
 *                Lehigh University
 *
 *     (C)opyright 2005 - J. Czyzyk, J. Linderoth and J. Shen
 *
 * $Id: SPProblem.h,v 1.31 2008/11/02 16:28:31 linderot Exp $
 */

#pragma once  

#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 


#ifndef SPPROBLEM_H
#define SPPROBLEM_H


extern "C" {
#include "mps.h"
#include "timefile.h"
#include "stochfile.h"
#include "split.h"
#include "memory.h"
#include "externsymbs.h"
}

#include <iostream>
#include <vector>

#if defined( USING_MW )
class MWRMComm;
#endif

// Some little stuff
const int MINIMIZE = 1;
const int MAXIMIZE = -1;
const double SPPROB_INFINITY = 1.0e+200;
const int MAX_EXPLICIT_SCENARIOS = 100000000;

typedef int ChangeType;
/// The maximum number of rows that will be return from the deterministic equivalent
const int MAX_ROWS_DE = 100000;
const int MAX_COLS_DE = 1000000;
const int MAX_NZ_DE = 6 * MAX_COLS_DE;

/** The different scenario handling modes in SUTIL
 */
enum SPScenarioMode {
   UNDEFINED_SCENARIOMODE = -1,
   EXPLICIT,
   EXPLICIT_SAMPLE,
   IMPLICIT_SAMPLE,
   IMPLICIT
};

enum SPSamplingMode {
   UNDEFINED_SAMPLINGMODE = -1,
   MONTE_CARLO,
   LATIN_HYPERCUBE,
   ANTITHETIC_VARIATES
};


class SPProblem 
{
   public:
	  SUTIL_DLL_API SPProblem();
	  SUTIL_DLL_API virtual ~SPProblem();

      /** Read the stochastic LP specified in SMPS format by the specified files.
          @param coreFile Name of Core File
       */
	  SUTIL_DLL_API void readFiles( char *coreFile, char *timeFile, char *stochFile );

      /// How are the scenarios handled.
	  SUTIL_DLL_API void setScenarioMode( SPScenarioMode m );

      /// How are the scenarios handles
	  SUTIL_DLL_API void setScenarioMode( const char scenario_mode[] );

      /// How should sampling be done
	  SUTIL_DLL_API void setSamplingMode( SPSamplingMode m );

      /// How is sampling being done?
	  SUTIL_DLL_API SPSamplingMode getSamplingMode( void ) const;

      /// How are scenarios being done
	  SUTIL_DLL_API SPScenarioMode getScenarioMode() const { return scenarioMode; }

      /// How many samples 
	  SUTIL_DLL_API void setNumSamples( int ns );

      /// How many aggregation periods
	  SUTIL_DLL_API void setNumAggrePeriods( int n );

      /// Returns the number of samples.
	  SUTIL_DLL_API int getNumSamples( void ) const;

      /// Calculate the total number of scenarios -- not the number in the sampled instance
	  SUTIL_DLL_API double numTotalScenarios() const;

      /// The number of scenarios NODES! in a period (in the problem being solved!)
	  SUTIL_DLL_API int numScenarios( int period ) const;

      /// Returns the number of nodes in the scenario tree in the problem being solved
	  SUTIL_DLL_API int numNodesInstance() const;

      /// #TRUE# if a scenario tree has been created.
	  bool scenarioTreeCreated;

      /** Create the scenario tree.

	  If #sample == TRUE#, then only create a full "sampled" scenario tree
	  with #nSamples# samples.
      */
	  SUTIL_DLL_API int createScenarioTree( void );

      /** 
       * Returns the size of the deterministic equivalent.
       * Parameters are double to ensure no integer overruns
       *
       * Does this actually only return an estimate of the number of rows, columns, and nz, since
       * what if some scenarios have many more nonzeroes?          
       */
	  SUTIL_DLL_API void getDeterministicEquivalentSize( double *nrows, double *ncols, double *nz ) const;

	  SUTIL_DLL_API void getDeterministicEquivalentSizeFromNode( int period, int ix,
                                                   double *nrows, double *ncols, double *nz ) const;



	  SUTIL_DLL_API int getBaseLPRows(int period, int *nrows) const;
	  SUTIL_DLL_API int getBaseLPCols(int period, int *ncols) const;
      /** Return the "base" linear program for a period. (The "recourse"
	  matrix and RHS 
      */
	  SUTIL_DLL_API int getBaseLP( int period, int *ncols, int *nrows,
		     int *nnzero, int *os, double **obj, double **rhs,
		     char **sense, int **matbeg, int **matind,
		     double **matval, double **bdl, double **bdu,
		     int **coltype = NULL) const;

      /** Returns the technology matrix for a period */
	  SUTIL_DLL_API int getTechnologyMatrix( int period, int *ncols, int *nrows,
			       int *nnzero, int **matbeg, int **matind, 
			       double **matval ) const;

      /** Returns the range of scenarios between #bix# and #eix# for #period# */
	  SUTIL_DLL_API int getRangeOfScenarios( int period, int bix, int eix,
			       double prob[], int nChanges[],
			       ChangeType *changeEntity[],
			       int *changeEntry[],
			       double *changeValue[] );

      /** Returns scenario indexed #ix# from period #period# 
          @param period.  The period of the node you wish to retrieve
          @param ix.  The index of the node you wish to retrieve
          @param prob.  The CONDITIONAL(?) probability of this scenario
       */
	  SUTIL_DLL_API int getScenario( int period, int ix, double *prob, int *nChanges,
		       ChangeType **changeEntity, int **changeEntry,
		       double **changeValue ) const;

      /** Returns scenario indexed #ix# regardless of which period.  NOT IMPLEMENTED! */
	  SUTIL_DLL_API int getScenario( int ix, double *prob, int *nChanges,
		       ChangeType **changeEntity, int **changeEntry,
		       double **changeValue ) const;


      /** This just returns a random scenario. */
	  SUTIL_DLL_API int getRandomScenario( double *p, int *nChanges, ChangeType **changeEntity,
			     int **changeEntry, double **changeValue );

      /** Given scenario index ix in period 'period', returns the parent index
       *   (in the previous period).  
       * @return -1 if period or ix is out of range
       */
	  SUTIL_DLL_API int getParentScenarioIndex(int period, int ix, int *parentPeriod, int *parentIx) const;

      /** 
       */
	  SUTIL_DLL_API int getChildScenarioIndices(int period, int ix, int *childPeriod,
				  int *nChildren, int **childIx) const;

      /**
       * Return the number of children for node(period, ix)
       */
	  SUTIL_DLL_API int getGrandParentScenarioIndex(int period, int ix,
                                      int *grandParentPeriod, int *grandParentIx);

      /**
       */
	  SUTIL_DLL_API int getGrandParentIndices(int period, const int *ix,
                                int nIx, int *nGrandParent, int **grandParentIx);

      /**
       * Given node (period,ix) returns the number of descendants.
       * @param dT Descendants time periods
       * @param dIx Descendants index
       */
	  SUTIL_DLL_API int getAllDescendants(int period, int ix, int *numDescendants, int **dT, int **dIx) const;

      /**
       */
	  SUTIL_DLL_API int getNumChildren(int period, int ix) const;

      /**
       * get Number of periods in the problem.
       */
	  SUTIL_DLL_API int getNumPeriods() const;

      /**
       * get Number of periods in the problem after aggregation.
       */
	  SUTIL_DLL_API int getNumAAPeriods() const;

      /**
       * get Number of periods to be aggregated in the problem.
       */
	  SUTIL_DLL_API int getNumAggrePeriods() const;

      /** Return probability of scenario #ix# in period #period# (
          This is the total probability
       */ 
	  SUTIL_DLL_API double getProbability( int period, int ix ) const;

      /** Return probability of scenario #ix# in period #period# (
          This is the conditional probability
       */ 
	  SUTIL_DLL_API double getCondProbability( int period, int ix ) const;
  
      /** Writes the deterministic equivalent in an MPS file */
	  SUTIL_DLL_API void writeDeterministicEquivalent( char *filename ) const;

      /** Returns the determinisic equivalent */
	  SUTIL_DLL_API int getDeterministicEquivalent(int *ncols, int *nrows, int *nnzero, int *os,
                                     double **obj, double **rhs, char **sense, 
                                     int **matbeg, int **matind, double **matval,
                                     double **bdl, double **bdu) const;
      
      /** Returns the determinisic equivalent from node (period, ix) */
	  SUTIL_DLL_API int getDeterministicEquivalentFromNode(int period, int ix,
                                             int *ncols, int *nrows, int *nnzero, int *os, 
                                             double **obj, double **rhs, char **sense, 
                                             int **matbeg, int **matind, double **matval,
                                             double **bdl, double **bdu) const;

      /** Creates and returns the "expected value" linear program */
	  SUTIL_DLL_API int getExpectedValueProblem( int *ncols, int *nrows,
				   int *nnzero, int *os, double **obj, 
				   double **rhs,
				   char **sense, int **matbeg, int **matind,
				   double **matval, double **bdl, 
				   double **bdu );
      

      /** This prints out some information */
	  SUTIL_DLL_API void printStochSummary( void ) const;

      ///
	  SUTIL_DLL_API void printStochInfo( void );

      ///
	  SUTIL_DLL_API void printSampledStochInfo( void );
#if defined( USING_PVM )
      /** Packs the class information into the PVM buffer */
      virtual int packInPVMBuffer( void ) const;

      /** Unpacks the class information from the PVM buffer */
      int unpackFromPVMBuffer( void );
#endif

#if defined( USING_MW )
      int MWPack( MWRMComm * ) const;
      int MWUnpack( MWRMComm * );
#endif

      /** Free the memory in the scenario Tree */
	  SUTIL_DLL_API void freeScenarioTree( void );

      /** Free the memory in the split Problem */
	  SUTIL_DLL_API void freeSplitProblem( void );

      /** Create a two stage equivalent problem.  NOT IMPLEMENTED IN BASE RELEASE. */
	  SUTIL_DLL_API void createTwoStage();

      /// 
	  SUTIL_DLL_API friend std::ostream& operator<<(std::ostream&, const SPProblem &);
  
   private:

      /// Number of periods in the problem.
      int numPeriods;

      /// Number of periods to be aggregated in the problem
      int numAggrePeriods;

      /// A cache value to hold the total number of scenarios
      mutable double numTotalScenarios_;

      /// Type that holds core file information
      MPStype *MPS;

      /// Type that holds time file information
      TimeType *Time;

      /// Type that holds the stoch file information
      StochInfoType *stochasticInfo;

      /** Type that holds the info if you sample.  It will create one
	  big "block" parameter, from this you can create the
	  scenario tree or deterministic equivalent or
	  whatever.
      */
      StochInfoType *sampledStochasticInfo;

      /** This is a "compact" form of the sampling information.
	  For discrete distribution and blocks, it contains the realization 
      */
    
      /// This is an ARRAY of TreeType (unlike the other two, which are pointers)
      TreeType *scenarioTree;

      ///
      SplitProblemType *splitProblem;

      ///
      SPScenarioMode scenarioMode;

      /// 
      SPSamplingMode samplingMode;
  
   public:
  
      /// This sets the Instance seed when you are sampling.
	   SUTIL_DLL_API int setInstanceSeed( long seed );

      /// Returns the current instance seed
	   SUTIL_DLL_API long getInstanceSeed( void ) const;

	   SUTIL_DLL_API bool doBootstrap() const { return doBootstrap_; }
	   bool doDoubleBootstrap() const { doDoubleBootstrap_; }
	   SUTIL_DLL_API void setDoBootstrap(bool b) { doBootstrap_ = b; }
	   SUTIL_DLL_API void setDoDoubleBootstrap(bool b) { doDoubleBootstrap_ = b; }
	   SUTIL_DLL_API void setBootstrapSeed(long int sval) { bootstrapSeed_ = sval; }
      
   private:

      bool doBootstrap_;
      bool doDoubleBootstrap_; 
      long int bootstrapSeed_;

      /** Generate a U(0,1) random value using Knuth's subtractive method.
	  To seed the sequence, enter a negative number.
      */
      double knuthRandom( long *idum ) const;

      /** This just gets the next random number in the sequence.    
       */
      double knuthRandom( void ) const;

      /** This is getting a random number using Linear Congruential Generators
       */
      long LCGRandom( long seed ) const;

      /// This sets the seed to the KnuthRandom() generator
      void setRandomSeed( long seed ) const;

      ///
      long instanceSeed;

      ///
      mutable long currentSeed;

      /// The number of samples (==-1 if whole scenario tree has been created).
      int nSamples;

      // Private helper functionss

      void getScenarioMonteCarlo(int period, int ix, double *prob, int *nChanges,
                                 ChangeType **changeEntity, int **changeEntry,
                                 double **changeValue) const;

      void getScenarioLatinHyperCube(int period, int ix, double *prob, int *nChanges,
                                     ChangeType **changeEntity, int **changeEntry,
                                     double **changeValue) const;

      void getScenarioAntitheticVariates(int period, int ix, double *prob, int *nChanges,
                                         ChangeType **changeEntity, int **changeEntry,
                                         double **changeValue) const;

      /** Returns the range of scenarios between #bix# and #eix# for #period# */
      void getRangeOfScenariosMonteCarlo(int period, int bix, int eix,
                                         double prob[], int nChanges[],
                                         ChangeType *changeEntity[],
                                         int *changeEntry[],
                                         double *changeValue[]);

      void getRangeOfScenariosLatinHyperCube(int period, int bix, int eix,
                                             double prob[], int nChanges[],
                                             ChangeType *changeEntity[],
                                             int *changeEntry[],
                                             double *changeValue[]);

      void getRangeOfScenariosAntitheticVariates(int period, int bix, int eix,
                                                 double prob[], int nChanges[],
                                                 ChangeType *changeEntity[],
                                                 int *changeEntry[],
                                                 double *changeValue[]);

      void createScenarioTreeMonteCarlo();
      void createScenarioTreeMonteCarloBootstrap();
      int createScenarioTreeLatinHypercube();
      int createScenarioTreeAntitheticVariates();
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "ellemtel")
// eval: (setq indent-tabs-mode nil)
