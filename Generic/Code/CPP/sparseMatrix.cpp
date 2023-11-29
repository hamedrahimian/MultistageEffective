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

#include <iostream>
#include "sparseMatrix.h"

double FindEntryValueInCSCMatrix(SparseCSC* mat, const int RowNum, const  int ColNum) {
	int pos = -1; //position of RowNum among non-zero elements of ColNum
	int element=0; //counter to loop over non-zero elements of ColNum
	for (element = mat->col_ptr[ColNum]; element< mat->col_ptr[ColNum + 1]; element++) {
		if (mat->row_ind[element] == RowNum) {
			return (mat->val[element]);
			pos = element - mat->col_ptr[ColNum];
		}
	}

	if (pos == -1) {
		std::cerr << "Unable to find match in FindEntryInCSCMatrix: error" << std::endl;
		return(-1);
	}
} //end FindEntryValueInCSCMatrix

int FindEntryColInCSCMatrix(const int * ncols, const int * element, int **matbeg) {
	int ColNum;

	for (ColNum = 0; ColNum< *ncols; ColNum++) {
		if (*element < (*matbeg)[ColNum + 1] && *element >= (*matbeg)[ColNum]) {
			return(ColNum);
			break;
		}
	}
} //end FindEntryColInCSCMatrix


//void AddEntriesCSCMatrix(SparseCSC* mat, const int nz, const double * value, const int * col, const int * row)
//{
//	int entry, oldnnzeros, oldrows, oldcols;
//
//	oldrows = mat->NumRows;
//	oldcols = mat->NumCols;
//	oldnnzeros = mat->NumEntries;
//
//
//
//	mat->NumEntries = oldnnzeros + nz;
//
//	for (entry = oldnnzeros; entry < mat->NumEntries; entry++) {
//
//	}
//	mat->NumRows = oldrows + 1;
//
//
//}// end AddEntriesCSCMatrix

void RemoveEntryCSCMatrix(SparseCSC* mat, const int ent)
{

	//NumCols and NumRows remain the same
	//row info remians the same, but we need to remove the entry of #ent
	int i, row, j, col; 
	int NumCols = mat->NumCols;
	int NumRows = mat->NumRows;
	int NumEntries = mat->NumEntries;

	col = FindEntryColInCSCMatrix(&NumCols, &ent, &mat->col_ptr);
	row = mat->row_ind[ent];
	
	int * new_row_ind = new  int[NumEntries - 1];
	double * new_val = new double[NumEntries - 1];
	int * new_col_ptr = new int[NumCols + 1];
	
	//update row index and val 
	for (i = 0; i < ent; i++) {
		new_row_ind[i] = mat->row_ind[i];
		new_val[i] = mat->val[i];
	}
	for (i = ent + 1; i < NumEntries; i++) {
		new_row_ind[i - 1] = mat->row_ind[i];
		new_val[i - 1] = mat->val[i];
	}

	//update col pointer
	for (j = 0; j<col; j++) {
		new_col_ptr[j] = mat->col_ptr[j];
	}

	for (j = col + 1; j<NumCols + 1; j++) {
		new_col_ptr[j] = mat->col_ptr[j] - 1;
	}

	//numbe rof non-zeros in col
	int nnz = mat->col_ptr[col + 1] - mat->col_ptr[col];

	if (nnz == 1)
		new_col_ptr[col] = new_col_ptr[col+1]; //this means col has no non-zero element
	else 
		new_col_ptr[col] = mat->col_ptr[col];
	
	
	mat->row_ind = new_row_ind;
	mat->val = new_val;
	mat->col_ptr = new_col_ptr;
	

	mat->NumEntries=NumEntries-1;
	mat->NumCols = NumCols;
	mat->NumRows = NumRows;


}//end RemoveEntryCSCMatrix
