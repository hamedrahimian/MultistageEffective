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
 * $Id: split.h,v 1.4 2005/03/04 23:38:22 linderot Exp $
 */
#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef SPLIT_INCLUDE_FILE
#define SPLIT_INCLUDE_FILE

#include "memory.h"
#include "sparse.h"
#include "mps.h"
#include "timefile.h"

typedef struct {

  int NumPeriods;

  SparseShape   *TShape, *WShape;
  SparseMatrix  *T, *W;

  VectorD       *b, *c, *UpBound, *LowBound, *Ranges;
  VectorI       *BoundType, *VarType;

  char         **RowType;

} SplitProblemType;
  
extern "C" SUTIL_DLL_API SplitProblemType *NewProblem(int np);
extern "C" SUTIL_DLL_API void DeleteProblem(SplitProblemType *ptr);
extern "C" SUTIL_DLL_API void SplitMPSfile(MPStype * MPS, TimeType* Time, SplitProblemType* SplitProblem);
extern "C" SUTIL_DLL_API int GetRowPeriod(int row, TimeType* Time);
extern "C" SUTIL_DLL_API int GetColPeriod(int col, TimeType* Time);
extern "C" SUTIL_DLL_API int FindEntry(MPStype* MPS, TimeType* Time, SplitProblemType* SplitProblem, char* Code, char* ColName, char* RowName, int* Period, int* Matrix, int* Entry);

#endif