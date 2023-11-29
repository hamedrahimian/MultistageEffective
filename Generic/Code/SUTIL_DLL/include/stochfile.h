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
 * $Id: stochfile.h,v 1.4 2005/03/04 23:38:22 linderot Exp $
 */
#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef STOCHFILE_H
#define STOCHFILE_H

static const int   UNKNOWN_TYPE      =  0;

static const int   INDEPENDENT_TYPE  =  1;
static const int   BLOCKS_TYPE       =  2;
static const int   SCENARIOS_TYPE    =  3;
static const int   EXIT_TYPE         =  4;

static const int   DISCRETE_TYPE     = 30;
static const int   UNIFORM_TYPE      = 31;
static const int   NORMAL_TYPE       = 32;
static const int   BETA_TYPE         = 33;
static const int   GAMMA_TYPE        = 34;
static const int   LOGNORMAL_TYPE    = 35;
static const int   SUBROUTINE_TYPE   = 36;
static const int   LINTR_TYPE        = 37;   /* linear transformation */

static const int      TT =  -2;
static const int      WW =  -3;
static const int    COST =  -4;
static const int     RHS =  -5;
static const int      LB =  -6;
static const int      UB =  -7;
static const int      RANGE = -8;
/* static const int   RANGE   -8 */
static const int   BOUND =  -9;
static const int      II = -10;
static const int      MI = -11;

#ifndef STOCHFILE_INCLUDE_FILE
#define STOCHFILE_INCLUDE_FILE

#include "split.h"

typedef struct {
  int      Period;
  int      NumRealizations;
  int      RealizationsSize;     /* number of spaces allocated */
  int      Matrix, Entry;
  double  *Probability;
  double  *Value;
} IndType;

typedef struct {
  double   Prob;
  int      NumEntries;
  int      EntriesSize;          /* number of spaces allocated */
  int     *Matrix, *Entry;
  double  *Value;
} BlockRealization;

typedef struct {
  int                Period;
  int                NumRealizations;
  int                RealizationsSize;  /* number of spaces allocated */
  BlockRealization  *Realization;
} BlockType;

typedef struct {
  char    *Name;
  double   Prob;
  int      BranchPeriod;
  int      ParentScenario;
  int      NumEntries;
  int      EntriesSize;                /* number of spaces allocated */
  int     *Period, *Matrix, *Entry;
  double  *Value;
} ScenarioType;

typedef struct {

  int    NumIndepType;
  int    NumBlockType;
  int    NumScenarioType;

  int    IndepSize;         /* number of spaces allocated */
  int    BlockSize;         /* number of spaces allocated */
  int    ScenarioSize;      /* number of spaces allocated */

  IndType      *IndParameter;
  BlockType    *BlockParameter;
  ScenarioType *ScenarioParameter;

  int    NumRootScenarios;
  int    RootScenariosSize;  /* number of spaces allocated */
  int   *RootScenario;
  

} StochInfoType;

typedef struct {
  int     ParentNode;
  int     NumEntries;
  int     NumT;
  int     NumChildren;
  int    *Child;
  double  Prob;
  int    *Matrix, *Entry;
  double *Value;
} TreeNodeType;

typedef struct {
  int           NumNodes;
  TreeNodeType *Node;
} TreeType;

#endif

extern "C" SUTIL_DLL_API  void  DeleteTreeTypeArray(int n, TreeType a[]);
extern "C" SUTIL_DLL_API int AllocateNewParameter(StochInfoType * Stoch, int Period, int Matrix, int Entry);
extern "C" SUTIL_DLL_API void AddNewParameterRealization(StochInfoType  *Stoch, int             ParamNum, double          Value, double Probability);
extern "C" SUTIL_DLL_API StochInfoType  *NewStoch();
extern "C" SUTIL_DLL_API int  AllocateNewBlock(StochInfoType * Stoch, int Period);
extern "C" SUTIL_DLL_API int AllocateNewBlockRealization(StochInfoType    *Stoch, int               BlockNum, double            Prob);
extern "C" SUTIL_DLL_API void AddNewBlockEntry(StochInfoType* Stoch, int BlockNum, int RealNum, int Matrix, int Entry, double Value);
extern "C" SUTIL_DLL_API int AllocateNewScenario(StochInfoType* Stoch, int BranchPeriod, int ParentScenario, char* Name, double Prob);
extern "C" SUTIL_DLL_API void AddNewScenEntry(StochInfoType* Stoch, int ScenNum, int Period, int Matrix, int Entry, double Value);
extern "C" SUTIL_DLL_API int GetScenNumberGivenName(StochInfoType* Stoch, char* Name);
extern "C" SUTIL_DLL_API void  PrintStoch(StochInfoType   *Stoch);
extern "C" SUTIL_DLL_API void PrintStochToFile(MPStype *MPS, TimeType* Time, StochInfoType* Stoch, SplitProblemType* SplitProblem, char* filename);
extern "C" SUTIL_DLL_API void ComputeRowColGivenEntry(MPStype* MPS, TimeType* Time, SplitProblemType* SplitProblem, int Period, int Matrix, int Entry, int* Row, int* Col);
extern "C" SUTIL_DLL_API void PrintMatrixName(int num);
extern "C" SUTIL_DLL_API void DeleteStoch(StochInfoType * Stoch, int NumberPeriods);

#endif
