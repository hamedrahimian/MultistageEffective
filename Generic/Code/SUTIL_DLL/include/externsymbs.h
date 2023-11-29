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
* $Id: time.c,v 1.4 2005/03/04 23:38:22 linderot Exp $
*/

#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef EXTERNSYMBS_INCLUDE_FILE
#define EXTERNSYMBS_INCLUDE_FILE

#include "mps.h"
#include "timefile.h"
#include "split.h"
#include "stochfile.h"


/*  This file contains the declarations of some routines defined in readmps.c, readstoch.c, readtime.c, output.c, and build.c
*
*  Written by Hamed Rahimian, September 2017.
*/


extern "C" SUTIL_DLL_API void ReadMPSFile(FILE *, MPStype *);

extern "C" SUTIL_DLL_API void ReadTimeFile(FILE *, TimeType *, MPStype *);

extern "C" SUTIL_DLL_API void ReadStochFile(FILE *, MPStype *, TimeType *, StochInfoType *,
	SplitProblemType *);

extern "C" SUTIL_DLL_API int BuildScenarios(MPStype *, TimeType *, StochInfoType *, TreeType **);

extern "C" SUTIL_DLL_API void  WriteDeterministicEquivalent(MPStype *, TimeType *, SplitProblemType *,
	TreeType *, char *);

extern "C" SUTIL_DLL_API int ReturnDeterministicEquivalent(TimeType *,
	SplitProblemType *,
	TreeType *,
	int *, int *, int *, int *,
	double *, double *, char *,
	int *, int *, double *,
	double *, double *,
	int, int, int);

extern "C" SUTIL_DLL_API int ReturnDeterministicEquivalentFromNode(int, int,
	TimeType *,
	SplitProblemType *,
	TreeType *,
	int *, int *, int *, int *,
	double *, double *, char *,
	int *, int *, double *,
	double *, double *,
	int, int, int);

extern "C" SUTIL_DLL_API void CreateNewTimeStartingFromPeriod(int , TimeType *,
	TimeType  **);

extern "C" SUTIL_DLL_API void CreateNewSplitProblemStartingFromPeriod(int ,
	TimeType *,
	SplitProblemType *,
	SplitProblemType **);


#endif

