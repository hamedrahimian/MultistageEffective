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
 * $Id: timefile.h,v 1.2 2005/03/04 23:38:22 linderot Exp $
 */
#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 


#ifndef TIMEFILE_INCLUDE_FILE
#define TIMEFILE_INCLUDE_FILE

typedef struct  {
  
  char *ProblemName;

  int PeriodSize;   /* indicates size of data structures */
  int NumPeriods;   /* indicates number of entries       */

  int *PeriodColumn; /* PeriodColumn[i] tells left-most column of period i */
  int *PeriodRow;    /* PeriodRow[i] tells first row of period i            */

  char **PeriodNames;

} TimeType;
  
extern "C" SUTIL_DLL_API TimeType *NewTime(int size);
extern "C" SUTIL_DLL_API void ResizeTimeMatrix(TimeType * Time, int size);
extern "C" SUTIL_DLL_API void DeleteTime(TimeType * Time);
extern "C" SUTIL_DLL_API void AddTimeProblemName(TimeType * Time, char* name);
extern "C" SUTIL_DLL_API int AddTimePeriod(TimeType* Time, char* PeriodName, char* ColName, char* RowName, MPStype * MPS);
extern "C" SUTIL_DLL_API void PrintTime(TimeType* Time);
extern "C" SUTIL_DLL_API int FindPeriod(TimeType * Time, char* Name);
extern "C" SUTIL_DLL_API void PrintTimeToFile(TimeType* Time, char* filename);

#endif
