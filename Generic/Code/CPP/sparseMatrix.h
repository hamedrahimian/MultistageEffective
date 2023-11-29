/*
*     Nested L-shaped method for DRSO-- Built on top of SUTIL
*
*     VERSION 0.1
*
*     Authors:   Hamed Rahimian
*                The Ohio State University
*
*        Guzin Bayraksan and Tito Homem-de-Mello
*
*		September 21, 2017
*
*     (C)opyright 2017 - H. Rahimian, G. Bayraksan, and T. Homem-de-Mello
*
*/

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

//TO DO: write it as a class

// Compressed sparse column (CSC or CCS) format to store a matrix in sparse format
//CSC is (val, row_ind, col_ptr), where 
//val is an array of the (top-to-bottom, then left-to-right) non-zero values of the matrix; 
//row_ind is the row indices corresponding to the values; and, 
//col_ptr is the list of val indices where each column starts.

typedef struct sparse_matrix_CSC {
	int NumEntries;
	double  *val;
	int     *col_ptr;
	int     *row_ind; 
	int NumRows;
	int NumCols;
} SparseCSC;

/**  Find entry in a CSC formatted Sparse matrix #mat, given #RowNum and #ColNum*/
double FindEntryValueInCSCMatrix(SparseCSC* mat, const int RowNum, const int ColNum);
/**   Find entry column in a CSC formatted Sparse matrix*/
int FindEntryColInCSCMatrix(const int * ncols, const int * element, int **matbeg);
/** Remove an entry from a CSC formatted Sparse Matrix*/
void RemoveEntryCSCMatrix(SparseCSC* mat, const int ent);
#endif