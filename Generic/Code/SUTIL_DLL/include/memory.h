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
 * $Id: memory.h,v 1.3 2007/05/29 19:04:47 jierui Exp $
 */

#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef MEMORY_INCLUDE_FILE
#define MEMORY_INCLUDE_FILE

#include <stdio.h>

typedef double *VectorD;
typedef int    *VectorI;

extern "C" SUTIL_DLL_API void    *SUTILMalloc(size_t      size, char    *message);
extern "C" SUTIL_DLL_API void    *Calloc(int num, size_t size, char *message);
extern "C" SUTIL_DLL_API char    *StrDup(char *s1, char *message);
extern "C" SUTIL_DLL_API double  *NewDouble(int      dimension, char    *message);
extern "C" SUTIL_DLL_API double  **NewDouble2(int       size1, int size2, char     *message);
extern "C" SUTIL_DLL_API int     *NewInt(int      size, char    *message);
extern "C" SUTIL_DLL_API int    **NewInt2(int  size1, int size2, char *message=NULL);
extern "C" SUTIL_DLL_API char    *NewChar(int   size, char *message);
extern "C" SUTIL_DLL_API void    *SUTILRealloc(void    *oldptr, int      size, char    *message);


extern  "C" SUTIL_DLL_API VectorD   *NewVectorD(int   size, char *message);
extern "C" SUTIL_DLL_API VectorD  **NewVectorDPtr(int   size, char *message);
extern "C" SUTIL_DLL_API VectorI   *NewVectorI(int   size, char *message);
extern "C" SUTIL_DLL_API VectorI  **NewVectorIPtr(int   size, char *message);
extern "C" SUTIL_DLL_API double  **NewDoublePtr(int   size, char *message);
extern "C" SUTIL_DLL_API int     **NewIntPtr(int   size, char *message);
extern "C" SUTIL_DLL_API char    **NewCharPtr(int   size, char *message);

extern "C" SUTIL_DLL_API void fail(char* message);


#endif
