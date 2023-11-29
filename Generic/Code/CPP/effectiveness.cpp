/*
*     Nested L-shaped method for DRSO and Risk-Averse Optimization-- Built on top of SUTIL
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

#include "effectiveness.h"
#include <algorithm>
#include <iostream>

using namespace std;

#define TOLER		0.000001

void Effectiveness::CalculateCDF(const int numChild, const double* nominal_prob, double* cdf)
{
	int j;

	cdf[0] = nominal_prob[0];
	for (j = 1; j < numChild; j++)
		cdf[j] = cdf[j-1] + nominal_prob[j];

}// end CalculateCDF

void Effectiveness::CalculateVaR(const int numChild, const double* cdf, const double * cost, int * indexVaR, double * VaR, const double beta)
{
	int j;

	if (cdf[0] >= beta) {
		*indexVaR = 0;
	}
	else {
		for (j = 1; j < numChild; j++) {
			if (cdf[j] >= beta - TOLER && cdf[j-1] < beta) {
				*indexVaR = j;
				break;
			}
		}
	}

	*VaR = cost[*indexVaR];


}//end CalculateVaR

void Effectiveness::StoreEff(const int it, ScenChildResult * scen, vector<int>* eff, IloInt* num_effect_category, const effect_category status, const bool IsCondEff)
{
	eff->push_back(it);
	(*num_effect_category)++;
	if (IsCondEff) {
		scen->CondStatus = status;
	}
	else {
		scen->PathStatus = status;
	}
}//end StoreEff

void Effectiveness::SwapEff(const int it, ScenChildResult * scen, vector<int>* eff_old, IloInt* num_effect_category_old, vector<int>* eff_new, IloInt *num_effect_category_new, const effect_category status_new)
{

	vector<int>::iterator match;
	match = std::find(eff_old->begin(), eff_old->end(), it);
	eff_old->erase(match);
	eff_new->push_back(it);
	(*num_effect_category_old)--;
	(*num_effect_category_new)++;
	scen->PathStatus = status_new;

}//end SwapEff

void Effectiveness::PrintCondEff(ostream& fe, const int w, const int period,
	OutputType* wOutput, const IloIntArray* num_effect_category)
{
	int j; 

	int numChild = wOutput->numChild;

	fe << "******* ancestor = " << w << " ** stage = " << period << " *******" << endl;
	fe << "w" << '\t' << '\t' << "h" << '\t' << "q" << '\t' << "p" << '\t' << "cdf" << '\t' << "Primal Category" << '\t' << "Status" << endl;
	for (j = 0; j<numChild; j++) {
		fe << wOutput->Child[j].No << '\t' << '\t' << wOutput->Child[j].Cost << '\t' << wOutput->Child[j].Nominal_Prob << '\t';
		fe << wOutput->Child[j].Worst_Prob << '\t' << wOutput->Child[j].cdf << '\t' << wOutput->Child[j].Primal_Category << '\t' << wOutput->Child[j].CondStatus << endl;
	}
	fe << *num_effect_category << endl;
	fe << endl;

}//end PrintCondEff

void Effectiveness::CheckNestednessPath(char *problem_name, vector<int>*** EffPathScen, int numInstance, const double * rho)
{

	char resName[100];
	const char* Resfilename = resName;
	sprintf_s(resName, "EX_Nestedness_Path_%s.txt", problem_name);
	Resfilename = resName;
	ofstream NestednessPath(Resfilename);
	NestednessPath.precision(10);

	vector<int>::iterator it;
	vector<int>::iterator match; 
	int nestedness = 1;

	const int numStage = getNumPeriods();
	
	int s, ix, w, ii;
	ii = 0;


	NestednessPath << "Nestedness of Effective Scenario Paths" << endl;
	NestednessPath << "**************************************" << endl;
	NestednessPath << endl;



	for (ii = numInstance-1; ii > 0; ii--) {
		nestedness = 1;
		NestednessPath << "******* gamma = " << rho[ii] << " *******" << endl;
		NestednessPath << "*************************************************" << endl;
		s = numStage - 2;
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			NestednessPath << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
			if ((*EffPathScen)[ii][w].size()>(*EffPathScen)[ii - 1][w].size()) {
				NestednessPath << "***** Effective scenario paths for node "<< w<< " in gamma " << rho[ii]<< " is larger than that in " << rho[ii - 1] << endl;
			}
			if ((*EffPathScen)[ii][w].size() > 0) {
				for (it = (*EffPathScen)[ii][w].begin(); it != (*EffPathScen)[ii][w].end(); it++) {
					match = std::find((*EffPathScen)[ii-1][w].begin(), (*EffPathScen)[ii-1][w].end(), *it);
					if (match == (*EffPathScen)[ii-1][w].end()) {
						NestednessPath << "Found " << getChildOverallIndex(s, ix, *it) <<  " in gamma " << rho[ii] << " does not exit in gamma " << rho[ii - 1] << endl;
						if (nestedness == 1)
							nestedness = 0;
					}
				}
			}
			NestednessPath << endl;

			ix++;
		}//end node

		if (nestedness == 1) {
			NestednessPath << "***** gamma " << rho[ii] << " is nested in gamma " << rho[ii-1] << endl;
		}
		NestednessPath << "*************************************************" << endl;


	}//end instance

	cout << "*****DONE with Checking Nestedness for Conditionally Effective Scenario Nodes*****" << endl;

}//end CheckNestednessPath

void Effectiveness::CheckNestednessCond(char *problem_name, vector<int>*** EffCondScen, int numInstance, const double * rho)
{

	char resName[100];
	const char* Resfilename = resName;
	sprintf_s(resName, "EX_Nestedness_Cond_%s.txt", problem_name);
	Resfilename = resName;
	ofstream NestednessCond(Resfilename);
	NestednessCond.precision(10);


	vector<int>::iterator it;
	vector<int>::iterator match;
	int nestedness = 1;

	const int numStage = getNumPeriods();

	int s, ix, w, ii;

	ii = 0;

	NestednessCond << "Nestedness of Conditionally Effective Scenario Nodes" << endl;
	NestednessCond << "***************************************************" << endl;
	NestednessCond << endl;


	for (ii = numInstance-1; ii > 0; ii--) {
		nestedness = 1;
		NestednessCond << "******* gamma = " << rho[ii] << " *******" << endl;
		NestednessCond << "*************************************************" << endl;
		for (s = 0; s < numStage - 1; s++) {
			ix = 0;
			while (ix < numScenarios(s)) {
				w = getNodeOverallIndex(s, ix);
				NestednessCond << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
				if ((*EffCondScen)[ii][w].size() > (*EffCondScen)[ii - 1][w].size()) {
					NestednessCond << "***** Condtionally effective scenario nodes for node " << w << " in gamma " << rho[ii]<< " is larger than that in " << rho[ii-1] << endl;
				}
				if ((*EffCondScen)[ii][w].size() > 0) {
					for (it = (*EffCondScen)[ii][w].begin(); it != (*EffCondScen)[ii][w].end(); it++) {
						match = std::find((*EffCondScen)[ii - 1][w].begin(), (*EffCondScen)[ii - 1][w].end(), *it);
						if (match == (*EffCondScen)[ii-1][w].end()) {
							NestednessCond << "Found " << getChildOverallIndex(s, ix, *it) << " in gamma " << rho[ii]<< " does not exit in gamma " << rho[ii - 1] << endl;
							if (nestedness == 1)
								nestedness = 0;
						}
					}
				}
				NestednessCond<< endl;
				ix++;

			}//end node
		}//end stage

		if (nestedness == 1) {
			NestednessCond << "***** gamma " << rho[ii] << " is nested in gamma " << rho[ii-1]<< endl;
		}
		NestednessCond << "*************************************************" << endl;


	}//end instance

	cout << "*****DONE with Checking Nestedness for Effective Scenario Paths*****" << endl;

}//end CheckNestednessPath

