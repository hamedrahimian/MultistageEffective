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

#include "effectivenessTV.h"
#include <algorithm>
#include <iostream>

using namespace std;

#define TOLER		0.000001
#define TOLERE		0.00000001


void EffectivenessTV::FormPrimalCategories(const int numChild, const double * cost, const int VaR_index, const double VaR_value,
	vector<int>** categories, cost_category* primal_category, IloIntArray* num_cost_category)
{
	int j, jj;

	vector<int> b_VaR;
	vector<int> VaR;
	vector<int> a_VaR;
	vector<int> sup;

	double sup_value = cost[numChild-1];
	double lambda = abs(sup_value - VaR_value);

	//form sup category
	sup.push_back(numChild - 1);
	j = numChild - 2;
	while (abs(cost[j] - sup_value) < TOLER && j >= 0) {
		sup.push_back(j);
		j--;
	}

	for (jj = j + 1; jj < numChild; jj++) {
		primal_category[jj] = cost_category::sup;
		(*num_cost_category)[3]++;
	}

	//form VaR category
	j = VaR_index - 1;
	if (VaR_index != 0) {
		while (abs(VaR_value - cost[j]) < TOLER && j >= 0) {
			VaR.push_back(j);
			primal_category[j] = cost_category::VaR;
			(*num_cost_category)[1]++;
			j--;
		}
	}


	//form b_VaR category

	for (jj = 0; jj < j + 1; jj++) {
		b_VaR.push_back(jj);
		primal_category[jj] = cost_category::b_VaR;
		(*num_cost_category)[0]++;
	}

	//continue forming VaR category
	j = VaR_index;
	while (abs(VaR_value - cost[j]) < TOLER && j<numChild) {
		VaR.push_back(j);
		primal_category[j] = cost_category::VaR;
		(*num_cost_category)[1]++;
		j++;
	}

	//form a_VaR category
	if (lambda > TOLER) {
		while (abs(cost[j] - sup_value) > TOLER && j < numChild) {
			a_VaR.push_back(j);
			primal_category[j] = cost_category::a_VaR;
			(*num_cost_category)[2]++;
			j++;
		}
	}

	if (lambda< TOLER) {
		(*num_cost_category)[1] = 0;
		(*num_cost_category)[2] = 0;
	}

	(*categories)[0] = b_VaR;
	(*categories)[1] = VaR;
	(*categories)[2] = a_VaR;
	(*categories)[3] = sup;

	b_VaR.~vector();
	VaR.~vector();
	a_VaR.~vector();
	sup.~vector();


} //end FormPrimalCategories

bool EffectivenessTV::ETCCondEffVaR(const int numChild, const int test_scenario, const double * sorted_cost, const double * sorted_nominal_prob,
	const int old_VaR_index, const double old_VaR_value, const double old_beta)
{
	int j;
	bool IsEffective = false;
	double * new_nominal_prob = new double[numChild];

	double remaining_prob = 1- sorted_nominal_prob[test_scenario];
	//adjsuted level
	double new_beta = (double)((old_beta - sorted_nominal_prob[test_scenario]) / remaining_prob);

	//obtain conditional probability
	for (j = 0; j < numChild; j++) {
		if (j != test_scenario)
			new_nominal_prob[j] = sorted_nominal_prob[j] / remaining_prob;
		else
			new_nominal_prob[j] = 0;
	}

	//calaculate new cdf
	double * cdf = new double[numChild];
	CalculateCDF(numChild, new_nominal_prob, cdf);

	//obtain new VaR and its index
	int new_VaR_index = 0;
	double new_VaR_value = 0;
	CalculateVaR(numChild, cdf, sorted_cost, &new_VaR_index, &new_VaR_value, new_beta);
	

	if (new_VaR_value >= old_VaR_value - TOLER)
		IsEffective=false;
	else if (cdf[new_VaR_index] > new_beta)
		IsEffective = true;
	else {
		//check where there exists a scenario between new_VaR and old_VaR with q>0
		bool ExistScenario = false;
		for (j = 0; j < numChild; j++) {
			if (j != test_scenario) {
				if (sorted_cost[j] > new_VaR_value  && sorted_cost[j] < old_VaR_value  && sorted_nominal_prob[j]>0) {
					ExistScenario = true;
					break;
				}
			}	
		}
		if (ExistScenario)
			IsEffective = true;
	}

	delete[] new_nominal_prob;
	delete[] cdf;

	return (IsEffective);
	
}//end ETCVaRCondEff

void EffectivenessTV::ETCAllCondEff(ostream& fe, const double rho, const IloNumArray* SubObj, const IloNumArray* Worst_p, OutputType** Output, vector<int>** IneffCondScen, vector<int>** EffCondScen, vector<int>** UnknownCondScen,
	IloIntArray2* num_cond_effect_category)
{

	int j, s, w, ix, jj;
	IloEnv env;

	const int numStage = getNumPeriods();
	const int numScenNode = numNodesInstance();
	const int numParents = numTotalNodesInstance(numStage - 2);

	IloIntArray2 num_cost_category(env, numParents);
	IloIntArray2 num_cond_effect_categoryTemp(env, numParents);
	num_cond_effect_category->clear();
	num_cond_effect_category->add(num_cond_effect_categoryTemp);

	for (w = 0; w < numParents; w++) {
		num_cost_category[w] = IloIntArray(env, 4); //categories 0= b_VaR; 1= VaR; 2= a_VaR; 3= sup
		num_cond_effect_categoryTemp[w] = IloIntArray(env, 3); //0= ineffective; 1= effective; 2=unknown
		(*num_cond_effect_category)[w] = num_cond_effect_categoryTemp[w];
	}
	num_cond_effect_categoryTemp.clear();
	num_cond_effect_categoryTemp.end();

	fe << "EASY-to-CHECK CONDITIONS: Conditional Effectiveness of Scenario Nodes.." << endl;
	fe << "******************************************************************" << endl;
	fe << endl;


	for (s = 0; s < numStage-1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);

			int ChildPeriod = 0;
			int numChild = 0;
			int * ChildIndices = NULL;
			int CHILDstat = 0;

			//Gets children scenarios of scenario node #ix in period #s
			CHILDstat = getChildScenarioIndices(s, ix, &ChildPeriod, &numChild, &ChildIndices);

			//overall index
			int * cIx = new int[numChild];
			double *q = new double[numChild];
			double * cost = new double [numChild];
			double * cdf= new double[numChild];
			cost_category * category_cost= new cost_category[numChild];

			OutputType wOutput;
			wOutput.numChild = numChild;
			wOutput.Child = vector<ScenChildResult> (numChild);
			for (j = 0; j < numChild; j++) {
				cIx[j] = getNodeOverallIndex(ChildPeriod, ChildIndices[j]);
				q[j] = getCondProbability(ChildPeriod, ChildIndices[j]);
				cost[j] = (*SubObj)[cIx[j]];
				wOutput.Child[j].ChildIndex = j;
				wOutput.Child[j].No = cIx[j];
				wOutput.Child[j].Cost = cost[j];
				wOutput.Child[j].Nominal_Prob = q[j];
				wOutput.Child[j].Worst_Prob = (*Worst_p)[cIx[j]];
			}


			//sort Output of scenario node #ix in period #s (overall index w) 
			sort(wOutput.Child.begin(), wOutput.Child.end(), CompareScenariosCostIndex());
			double sup_value = wOutput.Child[numChild-1].Cost;

			//obtain cdf

			double *sorted_q = new double[numChild];
			double * sorted_cost = new double[numChild];
			for (j = 0; j < numChild; j++) {
				sorted_q[j] = wOutput.Child[j].Nominal_Prob;
				sorted_cost[j] = wOutput.Child[j].Cost;
			}
			CalculateCDF(wOutput.numChild, sorted_q, cdf);
			/*cout << cdf[3] << endl;*/
			/*wOutput.Child[0].cdf = wOutput.Child[0].Nominal_Prob;
			for (j = 1; j < numChild; j++) 
				wOutput.Child[j].cdf = wOutput.Child[j-1].cdf + wOutput.Child[j].Nominal_Prob;*/

			//obtain VaR and its index
			int VaR_index=0;
			double VaR_value=0;


			CalculateVaR(wOutput.numChild, cdf, sorted_cost, &VaR_index, &VaR_value, rho/2);

			//form primal categories
			vector<int>* categories = new vector<int>[4];
			FormPrimalCategories(numChild, sorted_cost, VaR_index, VaR_value, &categories, 
				category_cost, &(num_cost_category[w]));

			
			/////////////NOW CHECK EFFECTIVENESS 
			vector<int>::iterator it;
			double lambda = abs(sup_value - VaR_value);

			//b_VaR category
			for (it = categories[0].begin(); it != categories[0].end(); it++)
				//cout << *it << endl;
				StoreEff(*it, &wOutput.Child[*it], &(*IneffCondScen)[w], &(*num_cond_effect_category)[w][0], effect_category::ineffective);
				//cout << (*num_cond_effect_category)[w][0] << endl;
				
			
			//sup category
			if (lambda > TOLER) {
				for (it = categories[3].begin(); it != categories[3].end(); it++) {
					/*cout << *it << endl;*/
					if (wOutput.Child[*it].Nominal_Prob <= rho / 2) {
						if (wOutput.Child[*it].Nominal_Prob > 0 || categories[3].size() == 1)
							StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
						else
							StoreEff(*it, &wOutput.Child[*it], &(*UnknownCondScen)[w], &(*num_cond_effect_category)[w][2], effect_category::unknown);
					}
					else
						StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
				}
			}

			//VaR category
			if (w == 64) {
				cout << "H" << endl;
			}
			if (lambda > TOLER) {
				for (it = categories[1].begin(); it != categories[1].end(); it++) {
					//cout << *it << endl;
					if (wOutput.Child[*it].Nominal_Prob <= rho / 2) {
						if (wOutput.Child[*it].Nominal_Prob > 0 && categories[1].size() > 1) {
							bool IsEffective = false;
							IsEffective = ETCCondEffVaR(numChild, *it, sorted_cost, sorted_q, VaR_index, VaR_value, rho/2);
							if (IsEffective)
								StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
							else
								StoreEff(*it, &wOutput.Child[*it], &(*UnknownCondScen)[w], &(*num_cond_effect_category)[w][2], effect_category::unknown);
						}
						else if (categories[1].size() == 1) {
							if (wOutput.Child[*it].Worst_Prob == 0)
								StoreEff(*it, &wOutput.Child[*it], &(*IneffCondScen)[w], &(*num_cond_effect_category)[w][0], effect_category::ineffective);
							else
								StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
						}
						else
							StoreEff(*it, &wOutput.Child[*it], &(*IneffCondScen)[w], &(*num_cond_effect_category)[w][0], effect_category::ineffective);
					}
					else
						StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
				}
			}

			if (lambda < TOLER) {
				for (it = categories[1].begin(); it != categories[1].end(); it++) {
					wOutput.Child[*it].Primal_Category = cost_category::sup;
					if (wOutput.Child[*it].Nominal_Prob <= rho / 2) {
						if (categories[3].size() == 1)
							StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
						else if (wOutput.Child[*it].Nominal_Prob>0){
							bool IsEffective = ETCCondEffVaR(numChild, *it, sorted_cost, sorted_q, VaR_index, VaR_value, rho);
							if (IsEffective)
								StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
							else
								StoreEff(*it, &wOutput.Child[*it], &(*UnknownCondScen)[w], &(*num_cond_effect_category)[w][2], effect_category::unknown);
						}
						else {
							StoreEff(*it, &wOutput.Child[*it], &(*UnknownCondScen)[w], &(*num_cond_effect_category)[w][2], effect_category::unknown);
						}
					}
					else 
						StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
				}
			}

			//a_VaR category
			if (lambda > TOLER) {
				for (it = categories[2].begin(); it != categories[2].end(); it++) {
					//cout << *it << endl;
					if (wOutput.Child[*it].Nominal_Prob <= rho / 2) {
						if (wOutput.Child[*it].Nominal_Prob > 0)
							StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
						else
							StoreEff(*it, &wOutput.Child[*it], &(*IneffCondScen)[w], &(*num_cond_effect_category)[w][0], effect_category::ineffective);
					}
					else
						StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
				}
			}

			for (j = 0; j < numChild; j++) {
				wOutput.Child[j].cdf = cdf[j];
				wOutput.Child[j].Primal_Category = category_cost[j];
			}
			//print results
			PrintCondEff(fe, w, s, &wOutput, &(*num_cond_effect_category)[w]);


			(*Output)[w] = wOutput;

			ix++;

			delete[] q;
			delete[] cost;
			delete[] sorted_cost;
			delete[] sorted_q;
			delete[] cdf;
			delete[] cIx;
			delete[] category_cost;
			delete[] categories;
		}//end ix
	}//end stage

	num_cost_category.clear();
	num_cost_category.end();
	

	cout << "*****DONE with Easy-to-Check Conditions for the Effectiveness of Scenario Nodes*****" << endl;

}//end ETCAllCondEff

void EffectivenessTV::ETCAllPathEff(ostream& fs, ostream& fe, const double rho, OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen,
	IloIntArray* num_path_effect_category)
{

	IloEnv env;
	IloIntArray num_path_effect_categoryTemp(env, 3);
	num_path_effect_category->clear();
	num_path_effect_category->add(num_path_effect_categoryTemp);
	num_path_effect_categoryTemp.clear();
	num_path_effect_categoryTemp.end();

	int ix, s, w, j;
	

	const int numStage = getNumPeriods();
	const int numScenNode = numNodesInstance();
	const int numParents = numTotalNodesInstance(numStage - 2);

	fe << endl;
	fe << "EASY-to-CHECK CONDITIONS: Effectiveness of Scenario Paths.." << endl;
	fe << "***********************************************************" << endl;

	//sort (*Output)[w] based on No for simplicity
	for (s = 0; s < numStage - 1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			int numChild = getNumChildren(s, ix);
			sort((*Output)[w].Child.begin(), (*Output)[w].Child.end(), CompareScenariosIndex());
			ix++;
		}
	}

	ix = 0;
	s = numStage - 2;
	while (ix < numScenarios(s)) {
	
		w = getNodeOverallIndex(s, ix);
		fe << endl;
		fe << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
		fe << "w" << '\t' << "Status" << endl;

		int ChildPeriod = 0;
		int numChild = 0;
		int * ChildIndices = NULL;
		int CHILDstat = 0;

		//Gets children scenarios of scenario node #ix in period #s
		CHILDstat = getChildScenarioIndices(s, ix, &ChildPeriod, &numChild, &ChildIndices);
		
		int *cIx = new int[numChild];
		for (j = 0; j < numChild; j++) {
			cIx[j] = getNodeOverallIndex(ChildPeriod, ChildIndices[j]);
			if ((*Output)[w].Child[j].Nominal_Prob <= rho / 2) {
				bool IsIneffFound = false;
				effect_category breaker;
				int effective = 1;

				if ((*Output)[w].Child[j].CondStatus == effect_category::ineffective) {
					IsIneffFound = true;
					breaker = (*Output)[w].Child[j].CondStatus;
				}
				if ((*Output)[w].Child[j].CondStatus == effect_category::effective) {
					effective *= 1;
				}
				else {
					effective = 0;
				}
				//set #ix in period #s be the current parent 
				int ParentIndex = ix;
				int ParentPeriod = s;
				while (ParentPeriod > 0 && !IsIneffFound) {
					//find the parent of the current parent
					int GrandParentIndex, GrandParentPeriod;
					int GrandPARENTstat = 0;
					GrandPARENTstat = getParentScenarioIndex(ParentPeriod, ParentIndex, &GrandParentPeriod, &GrandParentIndex);
					int w_GrandParent = getNodeOverallIndex(GrandParentPeriod, GrandParentIndex);
					//get the index of the current parent among its parent's children
					int j_Parent = getChildIndex(ParentPeriod, ParentIndex);
					if ((*Output)[w_GrandParent].Child[j_Parent].CondStatus == effect_category::ineffective) {
						IsIneffFound = true;
						breaker = (*Output)[w_GrandParent].Child[j_Parent].CondStatus;
					}
					else {
						if ((*Output)[w_GrandParent].Child[j_Parent].CondStatus == effect_category::effective) {
							effective *= 1;
						}
						else {
							effective = 0;
						}
					}
					//update the current parent
					ParentPeriod = GrandParentPeriod;
					ParentIndex = GrandParentIndex;
				}

				if (IsIneffFound) {
					StoreEff(j, &(*Output)[w].Child[j], &(*IneffPathScen)[w], &(*num_path_effect_category)[0], effect_category::ineffective, false);
				}
				else {
					if (effective == 0) {
						StoreEff(j, &(*Output)[w].Child[j], &(*UnknownPathScen)[w], &(*num_path_effect_category)[2], effect_category::unknown, false);
					}
					else {
						StoreEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, false);
					}
				}
			}//end if prob<=rho/2
			else
			{
				StoreEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, false);
			}
			fe << cIx[j] << '\t' << (*Output)[w].Child[j].PathStatus << endl;
		}//end Children of ix
		delete[] cIx;
		ix++;
	}//end while
		
	
	fe << *num_path_effect_category << endl;
	fs << rho << '\t' << (*num_path_effect_category)[0] << '\t' << (*num_path_effect_category)[1] << '\t' << (*num_path_effect_category)[2] << '\t';

	//sort Output as they were because the next task to run resolving unknown scenarios
	for (s = 0; s < numStage - 1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			int numChild = getNumChildren(s, ix);
			OutputType wOutPut = (*Output)[w];
			sort((*Output)[w].Child.begin(), (*Output)[w].Child.end(), CompareScenariosCostIndex());
			ix++;
		}
	}

	cout << "*****DONE with Easy-to-Check Conditions for the Effectiveness of Scenario Paths*****" << endl;

	env.end();

}//end ETCAllPathEff

void EffectivenessTV::ETCUnknownPathEff(ostream& fs, ostream& fe, const double rho, OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen,
	IloIntArray* num_path_effect_category)
{

	IloEnv env;

	int ix, s, w, j;

	vector<int>::iterator it;

	const int numStage = getNumPeriods();
	const int numScenNode = numNodesInstance();
	const int numParents = numTotalNodesInstance(numStage - 2);

	fe << endl;
	fe << "EASY-to-CHECK CONDITIONS: Effectiveness of Unknown Scenario Paths.." << endl;
	fe << "***********************************************************" << endl;

	//sort (*Output)[w] based on No for simplicity
	for (s = 0; s < numStage - 1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			int numChild = getNumChildren(s, ix);
			sort((*Output)[w].Child.begin(), (*Output)[w].Child.end(), CompareScenariosIndex());
			ix++;
		}
	}

	ix = 0;
	s = numStage - 2;
	while (ix < numScenarios(s)) {

		w = getNodeOverallIndex(s, ix);
		if ((*UnknownPathScen)[w].size() > 0) {
			fe << endl;
			fe << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
			fe << "w" << '\t' << "f(x*)" << '\t' << "f(x^)" << '\t' << "Status" << endl;
		}

		int ChildPeriod = 0;
		int numChild = 0;
		int * ChildIndices = NULL;
		int CHILDstat = 0;

		//Gets children scenarios of scenario node #ix in period #s
		CHILDstat = getChildScenarioIndices(s, ix, &ChildPeriod, &numChild, &ChildIndices);


		for (it = (*UnknownPathScen)[w].begin(); it != (*UnknownPathScen)[w].end(); it++) {
			j = (*Output)[w].Child[*it].ChildIndex;
			
			bool IsIneffFound = false;
			effect_category breaker;
			int effective = 1;

			if ((*Output)[w].Child[j].CondStatus == effect_category::ineffective) {
				IsIneffFound = true;
				breaker = (*Output)[w].Child[j].CondStatus;
			}
			if ((*Output)[w].Child[j].CondStatus == effect_category::effective) {
				effective *= 1;
			}
			else {
				effective = 0;
			}
			//set #ix in period #s be the current parent 
			int ParentIndex = ix;
			int ParentPeriod = s;
			while (ParentPeriod > 0 && !IsIneffFound) {
				//find the parent of the current parent
				int GrandParentIndex, GrandParentPeriod;
				int GrandPARENTstat = 0;
				GrandPARENTstat = getParentScenarioIndex(ParentPeriod, ParentIndex, &GrandParentPeriod, &GrandParentIndex);
				int w_GrandParent = getNodeOverallIndex(GrandParentPeriod, GrandParentIndex);
				//get the index of the current parent among its parent's children
				int j_Parent = getChildIndex(ParentPeriod, ParentIndex);
				if ((*Output)[w_GrandParent].Child[j_Parent].CondStatus == effect_category::ineffective) {
					IsIneffFound = true;
					breaker = (*Output)[w_GrandParent].Child[j_Parent].CondStatus;
				}
				else if ((*Output)[w_GrandParent].Child[j_Parent].CondStatus == effect_category::effective) {
						effective *= 1;
				}
				//update the current parent
				ParentPeriod = GrandParentPeriod;
				ParentIndex = GrandParentIndex;
			}

			if (IsIneffFound) {
				StoreEff(j, &(*Output)[w].Child[j], &(*IneffPathScen)[w], &(*num_path_effect_category)[0], effect_category::ineffective, false);
				(*num_path_effect_category)[2]--;
			}
			else if (effective == 1) {
					StoreEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, false);
					(*num_path_effect_category)[2]--;
			}
			
			fe << (*Output)[w].Child[j].No << '\t' << (*Output)[w].Child[j].PathStatus << endl;
		}//end Children of ix
		ix++;
	}//end while

	fe << *num_path_effect_category << endl;

	fs << (*num_path_effect_category)[0] << '\t' << (*num_path_effect_category)[1] << endl;

	//sort Output as they were because the next task to run resolving unknown scenarios
	for (s = 0; s < numStage - 1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			int numChild = getNumChildren(s, ix);
			OutputType wOutPut = (*Output)[w];
			sort((*Output)[w].Child.begin(), (*Output)[w].Child.end(), CompareScenariosCostIndex());
			ix++;
		}
	}

	cout << "*****DONE with Easy-to-Check Conditions for the Effectiveness of UNKNOWN Scenario Paths*****" << endl;

	env.end();

}//end ETCUnknownPathEff


void EffectivenessTV::ResolveUnknownCondEff(ostream& fe, const double rho, const IloNumArray2* xOpt, const IloNumArray* optSubObj, OutputType** Output, vector<int>** IneffCondScen, vector<int>** EffCondScen,
	vector<int>** UnknownCondScen, IloIntArray2* num_cond_effect_category)
{
	
	IloEnv env1;

	const int numStage = getNumPeriods();
	const int numScenNode = numNodesInstance();
	const int numParents = numTotalNodesInstance(numStage - 2);

	int s, ix, w, k, ww;

	//store the optimal cost function values and optimal solution from the original problem
	IloNumArray Best_SubObj(env1);
	Best_SubObj.add(*optSubObj);
	IloNumArray2 Best_x(env1);
	Best_x = *xOpt;
	for (w = 0; w < numScenNode; w++) {
		Best_x[w] = (*xOpt)[w];
	}
	

	vector<int>::iterator it;
	

	bool IsCondEff = true;
	bool IsPathEff = false;

	fe << endl;
	fe << "RESOLVING: Conditional Effectiveness of (Primarily) UNKNOWN Scenario Nodes.." << endl;
	fe << "*********************************************************************" << endl;

	for (s = 0; s < numStage-1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			if ((*UnknownCondScen)[w].size() > 0) {
				fe << endl;
				fe << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
				fe << "w" << '\t' << "f_" << s << "(x*)" << '\t' << "f_" << s << "(x^)" << '\t' << "Status" << endl;
			}
			//get the parent of node #ix in period #s
			int GrandParentPeriod, GrandParentIndex;
			int GrandPARENTstat = 0;
			if (s>0)
				GrandPARENTstat = getParentScenarioIndex(s, ix, &GrandParentPeriod, &GrandParentIndex);
			int w_GrandParent=0;
			if (s > 0)
				w_GrandParent = getNodeOverallIndex(GrandParentPeriod, GrandParentIndex);
			
			for (it = (*UnknownCondScen)[w].begin(); it != (*UnknownCondScen)[w].end(); it++) {
				int j_Child = (*Output)[w].Child[*it].ChildIndex;
				int jIx = getChildRelativeIndex(s, ix, j_Child);
				IloEnv env;

				double rel_Gap = IloInfinity;
				double Appx = -IloInfinity;
				double FuncVal = IloInfinity;
				double z_hat;
				int nCut = 0;

				///////////////// PARAMETERS ///////////////// 
				IloNumArray2 sol(env, numScenNode);
				IloNumArray2 pi(env, numScenNode);
				vector<gNode> gs(numParents);
				vector<GNode> Gs(numParents);
				vector<PNode> ps(numParents);
				IloNumArray SubObj(env, numScenNode);
				int objSense;

				////////////////////////////////////////////////

				bool contForwardPass = true;
				while (rel_Gap > TOLER && contForwardPass) {
					ForwardPassX(&sol, &pi, &Gs, &gs, &ps, &objSense, &Appx, &SubObj, nCut, DRSO, MULTI, IsCondEff, s, ix, jIx, &(*xOpt)[w_GrandParent]);
					IloNumArray pArray(env);
					vector<cut_P> p(numParents);
					//store the solution to the inner problem for the current iteration
					ForwardPassP(&p, &pArray, &SubObj, &objSense, nCut, &contForwardPass, DRSO, MULTI, rho, IsCondEff, IsPathEff, s, ix, jIx, j_Child);

					for (ww = 0; ww < numParents; ww++)
						ps[ww].push_back(p[ww]);

					p.~vector();
					pArray.clear();
					pArray.end();

					z_hat = SubObj[w];
					if (objSense == MINIMIZE) {
						if (z_hat < FuncVal) {
							FuncVal = z_hat;
						}
					}
					else {
						if (z_hat > FuncVal) {
							FuncVal = z_hat;
						}
					}
					rel_Gap = abs(FuncVal - Appx) / abs(Appx);
					//fe << nCut + 1 << '\t' << z_hat << '\t' << Appx << '\t' << FuncVal << '\t' << rel_Gap * 100 << endl;
					if (rel_Gap > TOLER && contForwardPass) {
						BackwardPass(&sol, &pi, &Gs, &gs, &ps, nCut, DRSO, MULTI, IsCondEff, s, ix, jIx, j_Child);
					}
					nCut++;
					
				}
				gs.~vector();
				Gs.~vector();
				ps.~vector();
				sol.clear();
				sol.end();
				pi.clear();
				pi.end();
				SubObj.clear();
				SubObj.end();
				env.end();
				

				if (contForwardPass) {
					(*Output)[w].Child[*it].CondObjCost = Appx; //the objective function of w after removing j_Child
					if (abs(Appx - Best_SubObj[w])/abs(Best_SubObj[w]) > TOLERE) {
						StoreEff(*it, &(*Output)[w].Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
						(*num_cond_effect_category)[w][2]--;
					}
					else {
						StoreEff(*it, &(*Output)[w].Child[*it], &(*IneffCondScen)[w], &(*num_cond_effect_category)[w][0], effect_category::ineffective);
						(*num_cond_effect_category)[w][2]--;
					}
				}
				else {
					StoreEff(*it, &(*Output)[w].Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
					(*num_cond_effect_category)[w][2]--;
				}

				fe << getChildOverallIndex(s, ix, j_Child) << '\t' << Best_SubObj[w] << '\t' << (*Output)[w].Child[*it].CondObjCost << '\t' << (*Output)[w].Child[*it].CondStatus << endl;
			}
			if ((*UnknownCondScen)[w].size() > 0) 
				fe << (*num_cond_effect_category)[w] << endl;
			ix++;
		}

	}

	Best_SubObj.clear();
	Best_SubObj.end();
	Best_x.clear();
	Best_x.end();
	env1.end();
	cout << "*****DONE with Resolving UNKNOWN Scenario Nodes*****" << endl;


}//end ResolveUnknownCondEff

void EffectivenessTV::ResolveAllCondEff(ostream& fe, const double rho, const IloNumArray2* xOpt, const IloNumArray* optSubObj, OutputType** Output, vector<int>** IneffCondScen, vector<int>** EffCondScen,
	vector<int>** UnknownCondScen, IloIntArray2* num_cond_effect_category)
{

}//end ResolveAllCondEff

void EffectivenessTV::ResolveUnknownPathEff(ostream& fs, ostream& fe, const double rho,
	OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen, IloIntArray* num_path_effect_category, const double objVal)
{

	int  s, j, ix, w, k, ww;

	bool IsCondEff = false;
	bool IsPathEff = true;


	vector<int>::iterator it;


	const int numStage = getNumPeriods();
	const int numScenNode = numNodesInstance();
	const int numParents = numTotalNodesInstance(numStage - 2);

	//sort (*Output)[w] based on No for simplicity
	for (s = 0; s < numStage - 1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			int numChild = getNumChildren(s, ix);
			sort((*Output)[w].Child.begin(), (*Output)[w].Child.end(), CompareScenariosIndex());
			ix++;
		}
	}

	fe << endl;
	fe << "RESOLVING: Effectiveness of UNKNOWN Scenario Paths.." << endl;
	fe << "********************************************" << endl;

	ix = 0;
	s = numStage - 2;
	while (ix < numScenarios(s)) {

		w = getNodeOverallIndex(s, ix);

		if ((*UnknownPathScen)[w].size() > 0) {
			fe << endl;
			fe << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
			fe << "w" << '\t' << "f(x*)" << '\t' << "f(x^)" << '\t' << "Status" << endl;
		}

		for (it = (*UnknownPathScen)[w].begin(); it != (*UnknownPathScen)[w].end(); it++) {
			int j_Child = (*Output)[w].Child[*it].ChildIndex;
			int jIx = getChildRelativeIndex(s, ix, j_Child);
			
			double rel_Gap = IloInfinity;
			double Appx = -IloInfinity;
			double FuncVal = IloInfinity;
			double z_hat;
			int nCut = 0;
			
			///////////////// PARAMETERS ///////////////// 
			IloEnv env;
			IloNumArray2 sol(env, numScenNode);
			IloNumArray2 pi(env, numScenNode);
			vector<gNode> gs(numParents);
			vector<GNode> Gs(numParents);
			vector<PNode> ps(numParents);
			IloNumArray SubObj(env, numScenNode);
			int objSense;

			////////////////////////////////////////////////

			bool contForwardPass = true;
			while (rel_Gap > TOLER && contForwardPass) {
				ForwardPassX(&sol, &pi, &Gs, &gs, &ps, &objSense, &Appx, &SubObj, nCut, DRSO, MULTI);
				IloNumArray pArray(env);
				vector<cut_P> p(numParents);
				ForwardPassP(&p, &pArray, &SubObj, &objSense, nCut, &contForwardPass, DRSO, MULTI, rho, IsCondEff, IsPathEff, s, ix, jIx, j_Child);


				for (ww = 0; ww < numParents; ww++)
					ps[ww].push_back(p[ww]);

				p.~vector();
				pArray.clear();
				pArray.end();

				z_hat = SubObj[0];
				if (objSense == MINIMIZE) {
					if (z_hat < FuncVal) {
						FuncVal = z_hat;
					}
				}
				else {
					if (z_hat > FuncVal) {
						FuncVal = z_hat;
					}
				}
				rel_Gap = abs(FuncVal - Appx) / abs(Appx);
				//fe << nCut + 1 << '\t' << z_hat << '\t' << Appx << '\t' << FuncVal << '\t' << rel_Gap * 100 << endl;
				if (rel_Gap > TOLER && contForwardPass) {
					BackwardPass(&sol, &pi, &Gs, &gs, &ps, nCut, DRSO, MULTI);
				}
				nCut++;
			}
			gs.~vector();
			Gs.~vector();
			ps.~vector();
			sol.clear();
			sol.end();
			pi.clear();
			pi.end();
			SubObj.clear();
			SubObj.end();
			env.end();

			
			if (contForwardPass) {
				(*Output)[w].Child[*it].PathObjCost = Appx; //the objective function of node 0 after removing j_child of w
				if (abs(Appx - objVal) / abs(objVal) > TOLERE) {
					if ((*Output)[w].Child[*it].PathStatus == effect_category::unknown) {
						StoreEff(*it, &(*Output)[w].Child[*it], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, IsCondEff);
						(*num_path_effect_category)[2]--;
					}
				}
				else {
					if ((*Output)[w].Child[*it].PathStatus == effect_category::unknown) {
						StoreEff(*it, &(*Output)[w].Child[*it], &(*IneffPathScen)[w], &(*num_path_effect_category)[0], effect_category::ineffective, IsCondEff);
						(*num_path_effect_category)[2]--;
					}
				}
			}
			else {
				(*Output)[w].Child[*it].PathObjCost = -1;
				if ((*Output)[w].Child[*it].PathStatus == effect_category::unknown) {
					StoreEff(*it, &(*Output)[w].Child[*it], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, IsCondEff);
					(*num_path_effect_category)[2]--;
				}
			}

			fe << getChildOverallIndex(s, ix, j_Child) << '\t' << objVal << '\t' << (*Output)[w].Child[*it].PathObjCost << '\t' << (*Output)[w].Child[*it].PathStatus << endl;
		
		}
		ix++;

	}


	fe << *num_path_effect_category << endl;
	
	fs << (*num_path_effect_category)[0] << '\t' << (*num_path_effect_category)[1] << endl;


	cout << "*****DONE with Resolving ALL Scenario Paths*****" << endl;


}//end ResolveUnknownPathEff

void EffectivenessTV::ResolveAllPathEff(ostream& fs, ostream& fe, const double rho,
	OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen, IloIntArray* num_path_effect_category, const double objVal)
{	

	int  s, j, ix, w, k, ww;

	bool IsCondEff = false;
	bool IsPathEff = true;
	bool IsEffFailed = false;
	bool IsIneffFailed = false;


	const int numStage = getNumPeriods();
	const int numScenNode = numNodesInstance();
	const int numParents = numTotalNodesInstance(numStage - 2);

	//sort (*Output)[w] based on No for simplicity
	for (s = 0; s < numStage - 1; s++) {
		ix = 0;
		while (ix < numScenarios(s)) {
			w = getNodeOverallIndex(s, ix);
			int numChild = getNumChildren(s, ix);
			sort((*Output)[w].Child.begin(), (*Output)[w].Child.end(), CompareScenariosIndex());
			ix++;
		}
	}

	fe << endl;
	fe << "RESOLVING: Effectiveness of Scenario Paths.." << endl;
	fe << "********************************************" << endl;

	ix = 0;
	s = numStage - 2;
	while (ix < numScenarios(s)) {

		w = getNodeOverallIndex(s, ix);
		fe << endl;
		fe << "******* ancestor = " << w << " ** stage = " << s << " *******" << endl;
		fe << "w" << '\t' << "f(x*)" << '\t' << "f(x^)" << '\t' << "Status" << '\t' << "Match" << endl;

		int ChildPeriod = 0;
		int numChild = 0;
		int * ChildIndices = NULL;
		int CHILDstat = 0;

		//Gets children scenarios of scenario node #ix in period #s
		CHILDstat = getChildScenarioIndices(s, ix, &ChildPeriod, &numChild, &ChildIndices);

		int *cIx = new int[numChild];
		for (j = 0; j < numChild; j++) {
			cIx[j] = getNodeOverallIndex(ChildPeriod, ChildIndices[j]);
			int jIx = getChildRelativeIndex(s, ix, j);

			string match;
			if ((*Output)[w].Child[j].Nominal_Prob <= rho / 2) {
				double rel_Gap = IloInfinity;
				double Appx = -IloInfinity;
				double FuncVal = IloInfinity;
				double z_hat;
				int nCut = 0;
				

				///////////////// PARAMETERS ///////////////// 
				IloEnv env;
				IloNumArray2 sol(env, numScenNode);
				IloNumArray2 pi(env, numScenNode);
				vector<gNode> gs(numParents);
				vector<GNode> Gs(numParents);
				vector<PNode> ps(numParents);
				IloNumArray SubObj(env, numScenNode);
				int objSense;

				////////////////////////////////////////////////

				bool contForwardPass = true;
				while (rel_Gap > TOLER && contForwardPass) {
					ForwardPassX(&sol, &pi, &Gs, &gs, &ps, &objSense, &Appx, &SubObj, nCut, DRSO, MULTI);
					IloNumArray pArray(env);
					vector<cut_P> p(numParents);
					ForwardPassP(&p, &pArray, &SubObj, &objSense, nCut, &contForwardPass, DRSO, MULTI, rho, IsCondEff, IsPathEff, s, ix, jIx, j);


					for (ww = 0; ww < numParents; ww++)
						ps[ww].push_back(p[ww]);

					p.~vector();
					pArray.clear();
					pArray.end();

					z_hat = SubObj[0];
					if (objSense == MINIMIZE) {
						if (z_hat < FuncVal) {
							FuncVal = z_hat;
						}
					}
					else {
						if (z_hat > FuncVal) {
							FuncVal = z_hat;
						}
					}
					rel_Gap = abs(FuncVal - Appx) / abs(Appx);
					//fe << nCut + 1 << '\t' << z_hat << '\t' << Appx << '\t' << FuncVal << '\t' << rel_Gap * 100 << endl;
					if (rel_Gap > TOLER && contForwardPass) {
						BackwardPass(&sol, &pi, &Gs, &gs, &ps, nCut, DRSO, MULTI);
					}
					nCut++;
				}
				gs.~vector();
				Gs.~vector();
				ps.~vector();
				sol.clear();
				sol.end();
				pi.clear();
				pi.end();
				SubObj.clear();
				SubObj.end();
				env.end();


				if (contForwardPass) {
					(*Output)[w].Child[j].PathObjCost = Appx; //the objective function of node 0 after removing j of w
					if (abs(Appx - objVal) / abs(objVal) > TOLERE) {
						if ((*Output)[w].Child[j].PathStatus == effect_category::unknown) {
							StoreEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, IsCondEff);
							(*num_path_effect_category)[2]--;
						}
						else {
							if ((*Output)[w].Child[j].PathStatus != effect_category::effective) {
								SwapEff(j, &(*Output)[w].Child[j], &(*IneffPathScen)[w], &(*num_path_effect_category)[0],
									&(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective);
								match = "NO";
								IsIneffFailed = true;
							}
						}

					}
					else {
						if ((*Output)[w].Child[j].PathStatus == effect_category::unknown) {
							StoreEff(j, &(*Output)[w].Child[j], &(*IneffPathScen)[w], &(*num_path_effect_category)[0], effect_category::ineffective, IsCondEff);
							(*num_path_effect_category)[2]--;
						}
						else {
							if ((*Output)[w].Child[j].PathStatus != effect_category::ineffective) {
								SwapEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1],
									&(*IneffPathScen)[w], &(*num_path_effect_category)[0], effect_category::ineffective);
								match = "NO";
								IsEffFailed = true;
							}
						}
					}
				}
				else {
					(*Output)[w].Child[j].PathObjCost = -1;
					if ((*Output)[w].Child[j].PathStatus == effect_category::unknown) {
						StoreEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, IsCondEff);
						(*num_path_effect_category)[2]--;
					}
					else {
						if ((*Output)[w].Child[j].PathStatus != effect_category::effective) {
							SwapEff(j, &(*Output)[w].Child[j], &(*IneffPathScen)[w], &(*num_path_effect_category)[0],
								&(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective);
							match = "NO";
							IsIneffFailed = true;
						}
					}
				}
			}//end if prob<= rho/2
			else {
				(*Output)[w].Child[j].PathObjCost = -1;
				if ((*Output)[w].Child[j].PathStatus == effect_category::unknown) {
					StoreEff(j, &(*Output)[w].Child[j], &(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective, IsCondEff);
					(*num_path_effect_category)[2]--;
				}
				else {
					if ((*Output)[w].Child[j].PathStatus != effect_category::effective) {
						SwapEff(j, &(*Output)[w].Child[j], &(*IneffPathScen)[w], &(*num_path_effect_category)[0],
							&(*EffPathScen)[w], &(*num_path_effect_category)[1], effect_category::effective);
						match = "NO";
						IsIneffFailed = true;
					}
				}
			}
			fe << cIx[j] << '\t' << objVal << '\t' << (*Output)[w].Child[j].PathObjCost << '\t' << (*Output)[w].Child[j].PathStatus << '\t' << match << endl;

		}//end children of ix

		ix++;
		delete[] cIx;

	}


	fe << *num_path_effect_category << endl;
	if (IsEffFailed) {
		fe << "*****Easy-to-Check Conditions for Effective Scenario Paths does NOT match with Resolving*****" << endl;
	}
	if (IsIneffFailed) {
		fe << "*****Easy-to-Check Conditions for Ineffective Scenario Paths does NOT match with Resolving*****" << endl;
	}
	fs << (*num_path_effect_category)[0] << '\t' << (*num_path_effect_category)[1] << endl;


	cout << "*****DONE with Resolving Scenario Paths*****" << endl;


}//end ResolveAllPathEff
