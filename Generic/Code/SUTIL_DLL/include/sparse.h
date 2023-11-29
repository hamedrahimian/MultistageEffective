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
 * $Id: sparse.h,v 1.5 2007/05/29 19:04:48 jierui Exp $
 */
#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef SPARSE_INCLUDE_FILE
#define SPARSE_INCLUDE_FILE

/* Sparse Matrix information structure */

typedef struct sparse_matrix_rec {
  double  *Value;
  int      Shape;
  int      Allocated;
} SparseMatrix;

typedef struct sparse_shape_rec {

  int      NumRows;
  int      NumCols;
  int      NumEnts;

  int     *pBeginRow;
  int     *pEndRow;

  int     *Row;

  int      Allocated;
} SparseShape;
  
extern "C" SUTIL_DLL_API void PrintMatrix(SparseShape   Shape, SparseMatrix  Matrix);
extern "C" SUTIL_DLL_API  int FindEntryInMatrix(SparseShape  *Shape, int  RowNum, int ColNum);
extern "C" SUTIL_DLL_API SparseMatrix   **NewMatrixPtr(size_t size, char *message);
extern "C" SUTIL_DLL_API SparseMatrix    *NewMatrix(size_t size, char *message);
//extern SparseShape    **NewShapePtr();
extern "C" SUTIL_DLL_API SparseShape     *NewShape(size_t size, char *message);
extern "C" SUTIL_DLL_API int DeleteShapeArray(int size, SparseShape * ptr);
extern "C" SUTIL_DLL_API int DeleteMatrixArray(int size, SparseMatrix * ptr);

#endif
