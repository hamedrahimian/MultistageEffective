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

#include "scenariotree.h"
#include <iostream>
#include <utility> //pair


using namespace std;

int TreeStucture::numTotalNodesInstance(const int period) const
{

	if (period >getNumPeriods() - 1) {
		printf("The scenario tree does not have %d periods\n", period);
		return 0;
	}

	if (period == 0)
		return 1;
	else
		return (numScenarios(period) + numTotalNodesInstance(period - 1));
}//end numTotalNodesInstance

int TreeStucture::getNodeOverallIndex(const int period, const int ix) const
{
	if (period == 0)
		return 0;
	else
		return (numTotalNodesInstance(period - 1) + ix);
}//end GetIndex

int TreeStucture::getChildIndex(const int period, const int ix) const
{
	int ChildIndex = 0;
	if (period == 0) {
		printf("This node is the root of the scenario tree");
		return -1;
	}
	int ParentSibilingsIx, ParentIndex, ParentPeriod;
	int PARENTstat = 0;
	PARENTstat = getParentScenarioIndex(period, ix, &ParentPeriod, &ParentIndex);

	ChildIndex = ix;
	for (ParentSibilingsIx = 0; ParentSibilingsIx < ParentIndex; ParentSibilingsIx++) {
		int numChild = getNumChildren(ParentPeriod, ParentSibilingsIx);
		ChildIndex -= numChild;
	}
	return  (ChildIndex);

}//end getChildIndex

int TreeStucture::getChildRelativeIndex(const int period, const int ix, const int j) const
{
	int pIx;
	int index = 0;
	for (pIx = 0; pIx < ix; pIx++)
		index += getNumChildren(period, pIx);

	return(index + j);

}//end getChildRelativeIndex

int TreeStucture::getChildOverallIndex(const int period, const int ix, const int j) const
{
	return(numTotalNodesInstance(period) + getChildRelativeIndex(period, ix, j));

}//end getChildOverallIndex

int TreeStucture::getDepth(const int period, const int ix, const int cperiod, const int cix) const
{
	int i, s;
	int * T;
	int * Inx;
	int ix_period;
	int numDescendants;

	if (cperiod <= period)
		return -1;
	else {
		int DESstat = 0;
		DESstat = getAllDescendants(period, ix, &numDescendants, &T, &Inx);
		if (DESstat != 1) {
			cerr << "getAllDescendants: error" << endl;
			return -1;
		}

		//vector<array<int, 2>> couple
		vector<pair <int, int>> Descendant;

		for (i = 0; i < numDescendants; i++) {
			//array<int, 2> temp = { T[i], Inx[i] };
			pair<int, int> temp(T[i], Inx[i]);
			Descendant.push_back(temp);
		};

		//array <int, 2>	 search = { cperiod, cix };
		pair<int, int> search(cperiod, cix);
		//vector<array<int, 2>>::iterator match;
		vector<pair<int, int>>::iterator match;

		match = std::find(Descendant.begin(), Descendant.end(), search);
		if (match != Descendant.end())
			return(cperiod - period);
		else
			return -1;

		Descendant.~vector();
	}

}// end getDepth


void TreeStucture::PrintScenariosToFile(char *problem_name)
{

	char resName[100];
	sprintf_s(resName, "EX_Scenarios_%s.txt", problem_name);
	const char* Resfilename = resName;
	ofstream fs(Resfilename);
	fs.precision(10);

	int s, ix, w, j;

	int col, row; //dimension of X variables and constraints counters
	double val; //value of a coefficient
	double bound; //current value of a bound in a range
	int change; //to loop over elements that should be changed in a scenario node
	int ent;


	const int numStage = getNumPeriods();

	for (s = 0; s < numStage - 1; s++) {

		IloEnv env;

		int os = 0;
		SparseCSC recourseMat;
		double *cost = NULL;
		double *rhs = NULL;
		char *sense = NULL;
		double *bdl = NULL;
		double *bdu = NULL;
		int * coltype = NULL;


		int LPstat = 0;
		//get subproblem technology matrix (W), coefficients, etc
		LPstat = getBaseLP(s + 1, &recourseMat.NumCols, &recourseMat.NumRows,
			&recourseMat.NumEntries, &os, &cost, &rhs,
			&sense, &recourseMat.col_ptr, &recourseMat.row_ind,
			&recourseMat.val, &bdl, &bdu,
			&coltype);


		if (LPstat != 0) {
			cerr << "getBaseLP: error" << endl;
			return;
		}


		//To hold Technology Matrix info
		SparseCSC techMat;

		int TECHstat = 0;
		//get subproblem Technology matrix
		TECHstat = getTechnologyMatrix(s + 1, &techMat.NumCols, &techMat.NumRows,
			&techMat.NumEntries, &techMat.col_ptr, &techMat.row_ind,
			&techMat.val);

		if (TECHstat != 0) {
			cerr << "getTechnologyMatrix: error" << endl;
			return;
		}


		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);

			int ChildPeriod = 0;
			int numChild = 0;
			int * ChildIndices = NULL;
			int CHILDstat = 0;

			//Gets children scenarios of scenario node #ix in period #s
			CHILDstat = getChildScenarioIndices(s, ix, &ChildPeriod, &numChild, &ChildIndices);

			fs << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;


			for (j = 0; j < numChild; j++) {

				double q = 0;
				int numChanges = 0;
				ChangeType *entityChange = NULL;
				int *entryChange = NULL;
				double *valueChange = NULL;

				int SCENstat = 0;
				//get scenario information
				SCENstat = getScenario(ChildPeriod, ChildIndices[j], &q, &numChanges, &entityChange, &entryChange, &valueChange);

				if (SCENstat != 0) {
					cerr << "getScenario: error" << endl;
					return;
				}

				int w_Child = getNodeOverallIndex(ChildPeriod, ChildIndices[j]);
				double prob = getCondProbability(ChildPeriod, ChildIndices[j]);
				fs << "******* node = " << w_Child << " ** q = " << prob << " *******" << endl;

				for (change = 0; change < numChanges; change++) {
					if (entityChange[change] == COST) {
						col = entryChange[change];
						val = valueChange[change];
						fs << "c." << col << '\t' << val << endl;
					}//COST

					if (entityChange[change] == RHS) {
						row = entryChange[change];
						val = valueChange[change];
						fs << "b." << row + 1 << '\t' << val << endl;
					} //RHS


					  //here sense should be outside varibles sense[period][row]
					if (entityChange[change] == RANGE) {
						row = entryChange[change];
						val = valueChange[change];
						fs << "range." << row + 1 << '\t' << val << endl;
					} //RANGE

					if (entityChange[change] == LB) {
						col = entryChange[change];
						val = valueChange[change];
						fs << "xL." << col << '\t' << val << endl;
					}//LB


					if (entityChange[change] == UB) {
						col = entryChange[change];
						val = valueChange[change];
						fs << "xU." << col << '\t' << val << endl;
					}//UB

					if (entityChange[change] == TT) {
						ent = entryChange[change];
						val = valueChange[change];
						col = FindEntryColInCSCMatrix(&techMat.NumCols, &ent, &techMat.col_ptr);
						row = techMat.row_ind[ent];
						fs << "T." << row + 1 << "." << col << '\t' << val << endl;
					}//TT technology matrix


					if (entityChange[change] == WW) {
						ent = entryChange[change];
						val = valueChange[change];
						col = FindEntryColInCSCMatrix(&recourseMat.NumCols, &ent, &recourseMat.col_ptr);
						row = recourseMat.row_ind[ent];
						fs << "W." << row + 1 << "." << col << '\t' << val << endl;
					}//WW recourse matrix
				}//change

			}//child
			ix++;
		}//ix
	}//stage

	fs.close();

}//end PrintScenariosToFile
