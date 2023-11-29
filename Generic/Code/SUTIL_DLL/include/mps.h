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
 * $Id: mps.h,v 1.4 2005/03/04 23:38:22 linderot Exp $
 */

#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef MPS_INCLUDE_FILE
#define MPS_INCLUDE_FILE

#include <stdio.h>
#include "hash.h"

static const int PINFTY      = 0;     /* normal */
static const int FREE        = 1;
static const int UPPER       = 2;
static const int LOWER       = 3;
static const int UPPERLOWER  = 4;
static const int FIX         = 5;
static const int MINFTY      = 6;

static const int SUTIL_VAR_CONTINUOUS = 0;
static const int SUTIL_VAR_INTEGER = 1;

typedef struct {

  int  NumRows;  /* indicates number of entries in data structures */
  int  NumCols;
  int  NumEnts;

  int  RowSize;  /* indicates size of data structures allocated */
  int  ColSize;
  int  EntSize;

  struct {
    int    *pBeginRow;
    int    *pEndRow;
    int    *Row;
    double *Value;
  } A;

  double *b, *c;
  double constant;
  int    *BoundType;
  double *UpBound;
  double *LowBound;
  double *Ranges;
  char   *RowType;
  int    *ColType;  

  char   *ProblemName;
  char   *ObjectiveName;
  char   *RHSName;
  char   *RangeName;
  char   *BoundName;

  char  **RowNames;
  char  **ColNames;

  HashTable *RowTable;
  HashTable *ColTable;

} MPStype;

extern "C" SUTIL_DLL_API MPStype *NewMPS( int rows, int cols, int ents );
extern "C" SUTIL_DLL_API void ResizeMPSMatrix(MPStype* MPS, int rows, int cols, int ents);
extern "C" SUTIL_DLL_API void DeleteMPS(MPStype *MPS);
extern "C" SUTIL_DLL_API void AddMPSProblemName(MPStype* MPS, char *name);
extern "C" SUTIL_DLL_API void AddMPSObjectiveRowName(MPStype * MPS, char *name);
extern "C" SUTIL_DLL_API int NoObjectiveRowName(MPStype * MPS);
extern "C" SUTIL_DLL_API int AddMPSRow(MPStype  *MPS, char* name, char* type);
extern "C" SUTIL_DLL_API int AddMPSColumn(MPStype* MPS, char* name, int VarType);
extern "C" SUTIL_DLL_API int AddMPSCoefficient(MPStype * MPS, char* colname, char* rowname, double value);
extern "C" SUTIL_DLL_API void AddMPSRHSName(MPStype *MPS, char* name);
extern "C" SUTIL_DLL_API int NoRHSName(MPStype * MPS);
extern "C" SUTIL_DLL_API void AddMPSRHS(MPStype *MPS, char* rowname, double value);
extern "C" SUTIL_DLL_API void AddMPSRangeName(MPStype * MPS, char *name);
extern "C" SUTIL_DLL_API int NoRangeName(MPStype * MPS);
extern "C" SUTIL_DLL_API void AddMPSRange(MPStype * MPS, char* rowname, double value);
extern "C" SUTIL_DLL_API void AddMPSBoundName(MPStype* MPS, char* name);
extern "C" SUTIL_DLL_API int NoBoundName(MPStype * MPS);
extern "C" SUTIL_DLL_API void AddMPSBound(MPStype* MPS, char* code, char* colname, double value);
extern "C" SUTIL_DLL_API void PrintMPS(MPStype * MPS);
extern "C" SUTIL_DLL_API void PrintMPSMatrixToFile(MPStype * MPS, char* filename);
extern "C" SUTIL_DLL_API void ParseHeaderLine(char* line, char* entry1, char* entry2);
extern "C" SUTIL_DLL_API void ParseDataLine(char* line, char* code, char* name1, char* name2, double* val1, char* name3, double* val2);
extern "C" SUTIL_DLL_API void string_copy(char* dest, char* string, int max);
extern "C" SUTIL_DLL_API int GetLine(char* line, FILE* file, int length);
extern "C" SUTIL_DLL_API int      GetNumberFromName( char **NameList, char *name, int size );
extern "C" SUTIL_DLL_API char    *GetNameFromNumber( char **NameList, int index );
extern "C" SUTIL_DLL_API int  GetRowNumber(MPStype* MPS, char* name);
extern "C" SUTIL_DLL_API int  GetColNumber(MPStype* MPS, char* name);

#endif
