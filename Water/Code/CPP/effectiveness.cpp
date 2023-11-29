/*
*     Nested L-shaped method for DRSO and Risk-Averse Optimization
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
#include "vars.h"
#include "tree.h"
#include "decomposition.h"
#include <algorithm> 
#include <stack>

using namespace std;

#define toler		0.00000001

void FormPrimalCategories(const double * cost, vector<int>** categories)
{
	int j, jj;

	vector<int> b_VaR;
	vector<int> VaR;
	vector<int> a_VaR;
	vector<int> sup;

	double sup_value = cost[numScen - 1];
	double var_value = cost[CVaRno];
	double lambda = abs(sup_value - var_value);

	//form sup category
	sup.push_back(numScen - 1);
	j = numScen - 2;
	while (abs(cost[j] - sup_value) < toler && j >= 0) {
		sup.push_back(j);
		j--;
	}


	//form VaR category
	j = CVaRno - 1;
	if (CVaRno != 0) {
		while (abs(var_value - cost[j]) < toler && j >= 0) {
			VaR.push_back(j);
			j--;
		}
	}


	//form b_VaR category
	for (jj = 0; jj < j + 1; jj++) {
		b_VaR.push_back(jj);
	}

	//continue forming VaR category
	j = CVaRno;
	while (abs(var_value - cost[j]) < toler && j<numScen) {
		VaR.push_back(j);
		j++;
	}

	//form a_VaR category
	if (lambda > toler) {
		while (abs(cost[j] - sup_value) > toler && j < numScen) {
			a_VaR.push_back(j);
			j++;
		}
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

void StoreEff(const int it, Scenario& scen, vector<int>& eff, IloInt& num_effect_category, const effect_category& status, const IloBool IsCondEff)
{
	eff.push_back(it);
	num_effect_category++;
	if (IsCondEff) {
		scen.CondStatus = status;
	}
	else {
		scen.PathStatus = status;
	}
}

void SwapEff(const int it, Scenario & scen, vector<int>& eff_old, IloInt & num_effect_category_old, vector<int>& eff_new, IloInt & num_effect_category_new, const effect_category & status_new)
{
	vector<int>::iterator match;
	match = std::find(eff_old.begin(), eff_old.end(), it);
	eff_old.erase(match);
	eff_new.push_back(it);
	num_effect_category_old--;
	num_effect_category_new++;
	scen.PathStatus = status_new;
}

void EasyCheckCondEff(ostream& Effectiveness, Scenario * Output, vector<int>& IneffCondScen, vector<int>& EffCondScen, vector<int>& UnknownCondScen, const pair<double, double> param, IloIntArray2& num_cost_category, IloIntArray2& num_cond_effect_category)
{
	IloInt w, w2, j, i;
	IloEnv env;
	double rho = param.first;

	IloIntArray2 num_cost_categoryTemp(env, numScenNode_stage_sum[numStage - 2]);
	IloIntArray2 num_cond_effect_categoryTemp(env, numScenNode_stage_sum[numStage - 2]);
	num_cost_category.clear();
	num_cond_effect_category.clear();
	num_cost_category.add(num_cost_categoryTemp);
	num_cond_effect_category.add(num_cond_effect_categoryTemp);
	for (w = 0; w < numScenNode_stage_sum[numStage - 2]; w++) {
		num_cost_categoryTemp[w] = IloIntArray(env, 4); //categories 0= b_VaR; 1= VaR; 2= a_VaR; 3= sup
		num_cond_effect_categoryTemp[w] = IloIntArray(env, 3); //0= ineffective; 1= effective; 2=unknown
		num_cost_category[w] = num_cost_categoryTemp[w];
		num_cond_effect_category[w] = num_cond_effect_categoryTemp[w];
	}
	num_cost_categoryTemp.end();
	num_cond_effect_categoryTemp.end();

	Effectiveness << "EASY-to-CHECK CONDITIONS: Conditional Effectiveness of Scenarios.." << endl;
	Effectiveness << "******************************************************************" << endl;
	Effectiveness << endl;

	for (w = 0; w < numScenNode_stage_sum[numStage - 2]; w++) {
		std::sort(Output + 1 + numScen*w, Output + 1 + numScen*(w + 1), CompareScenarios());

		Output[1 + numScen*w].partial_sum = Output[1 + numScen*w].nom_Prob;
		for (j = 1; j < numScen; j++) {
			Output[descendant[w][j]].partial_sum = Output[descendant[w][j] - 1].partial_sum + Output[descendant[w][j]].nom_Prob;
		}

		vector<int>::iterator it;
		vector<int> Sup_index;

		Sup_index.push_back(descendant[w][numScen - 1]);
		w2 = descendant[w][numScen - 1] - 1;
		while (abs(Output[w2].Cost - Output[Sup_index[0]].Cost) < toler && w2 > descendant[w][0] - 1) {
			Sup_index.push_back(w2);
			w2--;
		}

		for (int ww = w2 + 1; ww <= descendant[w][numScen - 1]; ww++) {
			Output[ww].primal_category = cost_category::sup;
			num_cost_category[w][3]++;
		}

		vector<int> VaR_index;
		if (Output[1 + numScen*w].partial_sum >= rho)
			VaR_index.push_back(1 + numScen*w);

		for (j = 1; j < numScen; j++) {
			if (Output[descendant[w][j]].partial_sum >= rho - toler && Output[descendant[w][j] - 1].partial_sum < rho) {
				VaR_index.push_back(descendant[w][j]);
				break;
			}
		}

		vector<int> VaR;
		if (VaR_index[0] != 1 + numScen*w) {
			w2 = VaR_index[0] - 1;
			while (abs(Output[VaR_index[0]].Cost - Output[w2].Cost) < toler && w2 > descendant[w][0] - 1) {
				VaR.push_back(w2);
				Output[w2].primal_category = cost_category::VaR;
				num_cost_category[w][1]++;
				w2--;
			}
		}
		if (VaR_index[0] != 1 + numScen*w)
			num_cost_category[w][0] = w2 - descendant[w][0] + 1;
		else
			num_cost_category[w][0] = 0;

		if (num_cost_category[w][0] != 0) {
			for (int ww = descendant[w][0]; ww <= w2; ww++) {
				Output[ww].primal_category = cost_category::b_VaR;
				StoreEff(ww, Output[ww], IneffCondScen, num_cond_effect_category[w][0], effect_category::ineffective);
			}
		}

		w2 = VaR_index[0];
		while (abs(Output[VaR_index[0]].Cost - Output[w2].Cost) < toler && w2 < descendant[w][numScen - 1] + 1) {
			VaR.push_back(w2);
			Output[w2].primal_category = cost_category::VaR;
			num_cost_category[w][1]++;
			w2++;
		}

		vector<int> third;
		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost) > toler) {
			while (abs(Output[w2].Cost - Output[Sup_index[0]].Cost) > toler && w2 < descendant[w][numScen - 1]) {
				third.push_back(w2);
				Output[w2].primal_category = cost_category::a_VaR;
				num_cost_category[w][2]++;
				w2++;
			}
		}

		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost)< toler) {
			num_cost_category[w][1] = 0;
			num_cost_category[w][2] = 0;
		}
		/////////////NOW CHECK EFFECTIVENESS 
		//Sup Category
		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost)> toler) {
			for (it = Sup_index.begin(); it != Sup_index.end(); it++) {
				if (Output[*it].nom_Prob <= rho) {
					if (Output[*it].nom_Prob>0 && num_cost_category[w][3]>1) {
						StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
					}
					else if (num_cost_category[w][3] == 1) {
						StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
					}
					else {
						StoreEff(*it, Output[*it], UnknownCondScen, num_cond_effect_category[w][2], effect_category::unknown);
					}
				}
				else {
					StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
				}
			}
		}

		//VaR category
		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost)> toler) {
			for (it = VaR.begin(); it != VaR.end(); it++) {
				if (Output[*it].nom_Prob <= rho ) {
					if (Output[*it].nom_Prob>0 && num_cost_category[w][1]>1) {
						StoreEff(*it, Output[*it], UnknownCondScen, num_cond_effect_category[w][2], effect_category::unknown);
					}
					else if (num_cost_category[w][1] == 1 && Output[*it].worst_Prob == 0) {
						StoreEff(*it, Output[*it], IneffCondScen, num_cond_effect_category[w][0], effect_category::ineffective);
					}
					else if (num_cost_category[w][1] == 1 && Output[*it].worst_Prob>0) {
						StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
					}
					else if (Output[*it].nom_Prob == 0 && num_cost_category[w][1]>1) {
						StoreEff(*it, Output[*it], IneffCondScen, num_cond_effect_category[w][0], effect_category::ineffective);
					}
				}
				else {
					StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
				}
			}
		}

		//lambda=0
		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost)< toler) {
			for (it = VaR.begin(); it != VaR.end(); it++) {
				Output[*it].primal_category = cost_category::sup;
				if (Output[*it].nom_Prob <= rho) {
					if (num_cost_category[w][3] == 1) {
						StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
					}
					else {
						StoreEff(*it, Output[*it], UnknownCondScen, num_cond_effect_category[w][2], effect_category::unknown);
					}
				}
				else {
					StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
				}
			}
		}

		if (abs(Output[VaR_index[0]].Cost - Output[Sup_index[0]].Cost)> toler) {
			for (it = third.begin(); it != third.end(); it++) {
				if (Output[*it].nom_Prob <= rho ) {
					if (Output[*it].nom_Prob>0) {
						StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
					}
					else {
						StoreEff(*it, Output[*it], IneffCondScen, num_cond_effect_category[w][0], effect_category::ineffective);
					}
				}
				else {
					StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w][1], effect_category::effective);
				}
			}
		}

		Effectiveness << "******* ancestor = " << w << " ** stage = " << stageMap[w] << " *******" << endl;
		Effectiveness << "w" << '\t' << '\t' << "h" << '\t' << "q" << '\t' << "p" << '\t' << "cdf" << '\t' << "Primal Category" << '\t' << "Status" << endl;
		for (j = 0; j<numScen; j++) {
			w2 = descendant[w][j];
			Effectiveness << Output[w2].No << '\t' << '\t' << Output[w2].Cost << '\t' << Output[w2].nom_Prob << '\t';
			Effectiveness << Output[w2].worst_Prob << '\t' << Output[w2].partial_sum << '\t' << Output[w2].primal_category << '\t' << Output[w2].CondStatus << endl;
		}
		Effectiveness << num_cond_effect_category[w] << endl;
		Effectiveness << endl;
	}
	cout << "*****DONE with Checking Conditional Effectiveness of Scenarios*****" << endl;


}

void EasyCheckCondEffReal(const IloInt& scenario, const IloInt& nCut, const IloNumArray& subObj_hat, const IloNumArray& p,
	effNodeType* EffCondScen, effNodeType* IneffCondScen, effNodeType* UnknownCondScen, const pair<double, double> param)
{
	IloInt w, w2, j, i;
	IloEnv env;
	double rho = param.first;
	double prob = (double)1 / numScen;

	NodeInfo * wOutput = new NodeInfo[numScen];
	for (j = 0; j < numScen; j++) {
		wOutput[j].No = j;
		wOutput[j].Cost = subObj_hat[descendant[scenario][j]];
		wOutput[j].worst_Prob = p[nCut*numScenNode+descendant[scenario][j]];
	}
	sort(wOutput, wOutput + numScen, CompareScenariosCostIndex());

	double * sorted_cost = new double[numScen];
	for (j = 0; j < numScen; j++) 
		sorted_cost[j] = wOutput[j].Cost;
	
	vector<int>* categories = new vector<int>[4];
	FormPrimalCategories(sorted_cost, &categories);
	
	int var_value = wOutput[CVaRno].Cost;
	int sup_value = wOutput[numScen - 1].Cost;
	
	/////////////NOW CHECK EFFECTIVENESS 
	vector<int>::iterator it;
	double lambda = abs(sup_value - var_value);

	//b_VaR category
	for (it = categories[0].begin(); it != categories[0].end(); it++) {
		int t = wOutput[*it].No;
		IneffCondScen->push_back(t);
	}
		
	//a_VaR category
	if (lambda > toler) {
		for (it = categories[3].begin(); it != categories[3].end(); it++) {
			int t = wOutput[*it].No;
			if (prob <= rho) {
				if (prob > 0 || categories[3].size() == 1) 		
					EffCondScen->push_back(t);
				else 
					UnknownCondScen->push_back(t);
			}
			else
				EffCondScen->push_back(t);
		}
	}
	
	//VaR category
	if (lambda > toler) {
		for (it = categories[1].begin(); it != categories[1].end(); it++) {
			int t = wOutput[*it].No;
			if (prob <= rho ) {
				if (prob > 0 && categories[1].size() > 1) {
					//bool IsEffective = false;
					//IsEffective = ETCCondEffVaR(numChild, *it, sorted_cost, sorted_q, VaR_index, VaR_value, rho / 2);
					//if (IsEffective)
					//	StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
					//else
					UnknownCondScen->push_back(t);
				}
				else if (categories[1].size() == 1) {
					if (wOutput[*it].worst_Prob == 0)
						IneffCondScen->push_back(t);
					else
						EffCondScen->push_back(t);
				}
				else
					IneffCondScen->push_back(t);
			}
			else
				EffCondScen->push_back(t);
		}
	}

	if (lambda < toler) {
		for (it = categories[1].begin(); it != categories[1].end(); it++) {
			int t = wOutput[*it].No;
			if (prob <= rho) {
				if (categories[3].size() == 1)
					EffCondScen->push_back(t);
				else if (prob>0) {
					//bool IsEffective = ETCCondEffVaR(numChild, *it, sorted_cost, sorted_q, VaR_index, VaR_value, rho);
					//if (IsEffective)
					//	StoreEff(*it, &wOutput.Child[*it], &(*EffCondScen)[w], &(*num_cond_effect_category)[w][1], effect_category::effective);
					//else
					  UnknownCondScen->push_back(t);
				}
				else {
					UnknownCondScen->push_back(t);
				}
			}
			else
				EffCondScen->push_back(t);
		}
	}

	//a_VaR category
	if (lambda > toler) {
		for (it = categories[2].begin(); it != categories[2].end(); it++) {
			int t = wOutput[*it].No;
			if (prob <= rho) {
				if (prob > 0)
					EffCondScen->push_back(t);
				else
					IneffCondScen->push_back(t);
			}
			else
				EffCondScen->push_back(t);
		}
	}


	delete[] sorted_cost;
	delete[] categories;
	delete[] wOutput;

}//end EasyCheckCondEffReal

void EasyCheckPathEff(ostream& Summary, ostream& Effectiveness, Scenario * Output, vector<int>& IneffPathScen, vector<int>& EffPathScen, vector<int>& UnknownPathScen, const pair<double, double> param, IloIntArray & num_path_effect_category)
{
	IloEnv env;
	IloIntArray num_path_effect_categoryTemp(env, 3);
	num_path_effect_category.clear();
	num_path_effect_category.add(num_path_effect_categoryTemp);
	num_path_effect_categoryTemp.end();

	IloInt s, w, w2, w_child, w_ancestor, w_grand_ancestor;
	double rho = param.first;
	IloInt start = numScenNode_stage_sum[numStage - 3];

	Effectiveness << endl;
	Effectiveness << "EASY-to-CHECK CONDITIONS: Effectiveness of Scenario Paths.." << endl;
	Effectiveness << "***********************************************************" << endl;

	//sort Output based on No for simplicity
	std::sort(Output, Output + numScenNode, CompareScenarios2());

	for (w_ancestor = start; w_ancestor < numScenNode_stage_sum[numStage - 2]; w_ancestor++) {
		Effectiveness << endl;
		Effectiveness << "******* ancestor = " << w_ancestor << " ** stage = " << stageMap[w_ancestor] << " *******" << endl;
		Effectiveness << "w" << '\t' << "Status" << endl;
		for (w_child = 0; w_child < numScen; w_child++) {
			IloBool IsIneffFound = IloFalse;
			effect_category breaker;
			int effective = 1;
			w = descendant[w_ancestor][w_child];
			if (Output[w].CondStatus == effect_category::ineffective) {
				IsIneffFound = IloTrue;
				breaker = Output[w].CondStatus;
			}
			else {
				if (Output[w].CondStatus == effect_category::effective) {
					effective = effective * 1;
				}
				else {
					effective = 0;
				}
				w2 = w;
				s = numStage - 1;
				while (s > 0 && !IsIneffFound) {
					w_grand_ancestor = FindAncestor(s,w2);
					if (w_grand_ancestor != 0) {
						if (Output[w_grand_ancestor].CondStatus == effect_category::ineffective) {
							IsIneffFound = IloTrue;
							breaker = Output[w_grand_ancestor].CondStatus;
						}
						else {
							if (Output[w_grand_ancestor].CondStatus == effect_category::effective) {
								effective = effective * 1;
							}
							else {
								effective = 0;
							}
						}
					}
					s--;
					w2 = w_grand_ancestor;
				}
			}
			if (IsIneffFound) {
				StoreEff(w, Output[w], IneffPathScen, num_path_effect_category[0], effect_category::ineffective, IloFalse);
			}
			else {
				if (effective == 0) {
					StoreEff(w, Output[w], UnknownPathScen, num_path_effect_category[2], effect_category::unknown, IloFalse);
				}
				else {
					StoreEff(w, Output[w], EffPathScen, num_path_effect_category[1], effect_category::effective, IloFalse);
				}
			}
			Effectiveness << w << '\t' << Output[w].PathStatus << endl;
		}
	}
	Effectiveness << num_path_effect_category << endl;
	Summary << rho << '\t' << num_path_effect_category[0] << '\t' << num_path_effect_category[1] << '\t' << num_path_effect_category[2] << '\t' << endl;

	//sort Output as they were because the next task to run resolving unknown scenarios
	for (w = 0; w < numScenNode_stage_sum[numStage - 2]; w++) {
		std::sort(Output + 1 + numScen*w, Output + 1 + numScen*(w + 1), CompareScenarios());
	}

	cout << "*****DONE with Checking Effectiveness of Scenario Paths*****" << endl;
}

void EasyCheckEffPath(const IloInt& nCut, const IloNumArray& p, vector<effpathNode>* patheffs)
{
	
	int s, j, w;
	
	int tmpT = 0;
	int tmpIx = 0;


	int numParents = numScenNode_stage_sum[numStage - 2];

	effect_category status=effect_category::effective;
	(*patheffs)[0].push_back(status);

	for (w = 1; w < numParents; w++)
		(*patheffs)[w].push_back(effect_category::unknown);

	w = 1;
	for (s = 1; s < numStage-1; s++) {
		while (w < numScenNode_stage_sum[s]) {
			
			if ((*patheffs)[w][nCut]== effect_category::unknown ) {
				if (p[nCut*numScenNode + w] == 0) {
					//w is ineffective; and the whole subtree of w should be ineffective
					stack<int> sT; //period
					stack<int> sIx; //node index
					sT.push(s);
					sIx.push(w);
					while (!sT.empty()) {
						tmpT = sT.top();
						tmpIx = sIx.top();
						sT.pop();
						sIx.pop();
						(*patheffs)[tmpIx][nCut] = effect_category::ineffective;
						if (tmpT < numStage - 2) {
							int childPeriod = tmpT + 1;
							vector<int> childIx;
							for (j = numScen - 1; j >= 0; j--) {
								childIx.push_back(descendant[tmpIx][j]);
								sT.push(childPeriod);
								sIx.push(childIx[numScen - 1 - j]);
							}
						}
					}
				}
				else if (p[nCut*numScenNode + w] >0) {
					(*patheffs)[w][nCut] = effect_category::effective;
					//Cutsforms << "nCut=" << nCut << ", w=" << w << endl;
					//(*lastpatheff)[w] = nCut;

				}
			}	
			w++;
			
		}

	}
	

	
	

}//end EasyCheckEffPath

void ResolveCondEff(ostream& Effectiveness, ModelArray MODEL, CplexArray CPX, Formulation formulation,
	IloCplex CplexSub, IloModel modSub, DistGenFormulation distgen,
	const pair<double, double> param,
	Scenario * Output, vector<int>& IneffCondScen, vector<int>& EffCondScen, vector<int>& UnknownCondScen, IloIntArray2& num_cost_category, IloIntArray2& num_cond_effect_category)
{
	vector<int>::iterator it;
	IloEnv env = MODEL[0].getEnv();
	//store the costs from the original problem
	IloNumArray Best_subObj_hat(env);
	Best_subObj_hat.add(formulation.costsix.subObj_hat);
	IloNumArray Best_y_grand_ancestor(env);
	Best_y_grand_ancestor.add(formulation.optsols.Best_y); //check later
	IloInt i, w;
	IloInt old_ancestor = 0;

	IloBool IsCondEff = IloTrue;
	IloBool IsPathEff = IloFalse;

	Effectiveness << endl;
	Effectiveness << "RESOLVING: Conditional Effectiveness of (Primarily) Unknown Scenarios.." << endl;
	Effectiveness << "*********************************************************************" << endl;
	if (UnknownCondScen.size() > 0) {
		for (it = UnknownCondScen.begin(); it != UnknownCondScen.end(); it++) {
			IloInt w_child = Output[*it].No;
			IloInt w_ancestor = FindAncestor(w_child, stageMap[w_child]);
			IloInt w_grand_ancestor = FindAncestor(w_ancestor, stageMap[w_ancestor]);
			IloInt j = FindChild(stageMap[w_child], w_child);
			if (w_ancestor != old_ancestor && it != UnknownCondScen.begin()) {
				Effectiveness << num_cond_effect_category[old_ancestor] << endl;
			}
			if (w_ancestor != old_ancestor || it == UnknownCondScen.begin()) {
				Effectiveness << endl;
				Effectiveness << "******* ancestor = " << w_ancestor << " ** stage = " << stageMap[w_ancestor] << " *******" << endl;
				Effectiveness << "w" << '\t' << "f_" << stageMap[w_ancestor] << "(x*)" << '\t' << "f_" << stageMap[w_ancestor] << "(x^)" << '\t' << "Status" << endl;
			}
			Termination(MODEL, formulation.cuts, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.optsols);
			Initialization(env, formulation.cuts, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.optsols, formulation.costsix);


			IloNum rel_Gap = IloInfinity;
			IloNum LB = -IloInfinity;
			IloNum UB = IloInfinity;
			IloNum z_hat;
			IloInt nCut = 0;

			IloBool contForwardPass = IloTrue;
			while (rel_Gap > toler && contForwardPass) {
				ForwardPassX(nCut, CPX, formulation.constraints, formulation.cuts, formulation.mainvariables, formulation.appxvariables,
					formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.costsix, LB, param, IsCondEff, stageMap[w_ancestor], w_ancestor, Best_y_grand_ancestor);

				ForwardPassP(nCut, CplexSub, modSub, distgen, formulation.incsols.p, formulation.costsix.subObj_hat, contForwardPass, param, IsCondEff, IsPathEff, stageMap[w_ancestor], w_ancestor, j);

				if (contForwardPass) {
					z_hat = formulation.costsix.subObj_hat[w_ancestor];
					if (z_hat < UB) {
						UB = z_hat;
						for (w = 0; w < numScenNode; w++) {
							for (i = 0; i < numRecharge; i++)
								formulation.optsols.Best_y[i*numScenNode + w] = formulation.incsols.y[i*numScenNode + w];
							formulation.optsols.Worst_p[w] = formulation.incsols.p[nCut*numScenNode + w];
						}
					}
					rel_Gap = abs(UB - LB) / abs(LB);
					if (rel_Gap > toler) {
						BackwardPass(nCut, MODEL, CPX, formulation.constraints, formulation.cuts, formulation.mainvariables.Y,
							formulation.appxvariables, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.costsix, param, stageMap[w_ancestor], w_ancestor);

					}
				}
				nCut++;
			}
			if (contForwardPass) {
				Output[*it].CondObjCost = LB; //the objective function of w_ancestor after removing w_child
				if (abs(Output[*it].CondObjCost - Best_subObj_hat[w_ancestor]) > toler) {
					StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w_ancestor][1], effect_category::effective);
					num_cond_effect_category[w_ancestor][2]--;
				}
				else {
					StoreEff(*it, Output[*it], IneffCondScen, num_cond_effect_category[w_ancestor][0], effect_category::ineffective);
					num_cond_effect_category[w_ancestor][2]--;
				}
			}
			else {
				StoreEff(*it, Output[*it], EffCondScen, num_cond_effect_category[w_ancestor][1], effect_category::effective);
				num_cond_effect_category[w_ancestor][2]--;
			}
			Effectiveness << w_child << '\t' << Best_subObj_hat[w_ancestor] << '\t' << Output[*it].CondObjCost << '\t' << Output[*it].CondStatus << endl;
			old_ancestor = w_ancestor;
		}
		Effectiveness << num_cond_effect_category[old_ancestor] << endl;
	}
	cout << "*****DONE with Resolving Unknown Scenarios*****" << endl;
}

void ResolvePathEff(ostream& Summary, ostream& Effectiveness,
	ModelArray MODEL, CplexArray CPX, Formulation formulation,
	IloCplex CplexSub, IloModel modSub, DistGenFormulation distgen,
	const IloNum objVal,
	const pair<double, double> param,
	Scenario * Output, vector<int>& IneffPathScen, vector<int>& EffPathScen, vector<int>& UnknownPathScen, IloIntArray& num_path_effect_category)
{
	IloEnv env = MODEL[0].getEnv();
	IloInt i, w, w2;
	IloBool IsCondEff = IloFalse;
	IloBool IsPathEff = IloTrue;
	IloBool IsEffFailed = IloFalse;
	IloBool IsIneffFailed = IloFalse;

	IloInt w_child, w_ancestor;
	IloInt start = numScenNode_stage_sum[numStage - 3];

	Effectiveness << endl;
	Effectiveness << "RESOLVING: Effectiveness of Scenario Paths.." << endl;
	Effectiveness << "********************************************" << endl;

	std::sort(Output, Output + numScenNode, CompareScenarios2());

	for (w_ancestor = start; w_ancestor < numScenNode_stage_sum[numStage - 2]; w_ancestor++) {
		Effectiveness << endl;
		Effectiveness << "******* ancestor = " << w_ancestor << " ** stage = " << stageMap[w_ancestor] << " *******" << endl;
		Effectiveness << "w" << '\t' << "f(x*)" << '\t' << "f(x^)" << '\t' << "Status" << '\t' << "Match" << endl;
		for (w_child = 0; w_child < numScen; w_child++) {
			Termination(MODEL, formulation.cuts, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.optsols);
			Initialization(env, formulation.cuts, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.optsols, formulation.costsix);

			IloNum rel_Gap = IloInfinity;
			IloNum LB = -IloInfinity;
			IloNum UB = IloInfinity;
			IloNum z_hat;
			IloInt nCut = 0;
			string match;

			IloBool contForwardPass = IloTrue;
			while (rel_Gap > toler && contForwardPass) {
				ForwardPassX(nCut, CPX, formulation.constraints, formulation.cuts, formulation.mainvariables, formulation.appxvariables,
					formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.costsix, LB, param);

				ForwardPassP(nCut, CplexSub, modSub, distgen, formulation.incsols.p, formulation.costsix.subObj_hat, contForwardPass, param, IsCondEff, IsPathEff, 0, w_ancestor, w_child);


				if (contForwardPass) {
					z_hat = formulation.costsix.subObj_hat[0];
					if (z_hat < UB) {
						UB = z_hat;
						for (w = 0; w < numScenNode; w++) {
							for (i = 0; i < numRecharge; i++)
								formulation.optsols.Best_y[i*numScenNode + w] = formulation.incsols.y[i*numScenNode + w];
							formulation.optsols.Worst_p[w] = formulation.incsols.p[nCut*numScenNode + w];
						}

					}
					rel_Gap = abs(UB - LB) / abs(LB);
					//Effectiveness << nCut + 1 << '\t' << z_hat << '\t' << LB << '\t' << UB << '\t' << rel_Gap * 100 << endl;
					if (rel_Gap > toler) {
						BackwardPass(nCut, MODEL, CPX, formulation.constraints, formulation.cuts, formulation.mainvariables.Y,
							formulation.appxvariables, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.costsix, param);
					}
				}
				nCut++;
			}

			w2 = descendant[w_ancestor][w_child];
			if (contForwardPass) {
				Output[w2].PathObjCost = LB; //the objective function of node 0 after removing w_child of w_ancestor
				if (abs(LB - objVal) / abs(objVal) > toler) {
					if (Output[w2].PathStatus == effect_category::unknown) {
						StoreEff(w2, Output[w2], EffPathScen, num_path_effect_category[1], effect_category::effective, IsCondEff);
						num_path_effect_category[2]--;
					}
					else {
						if (Output[w2].PathStatus != effect_category::effective) {
							SwapEff(w2, Output[w2], IneffPathScen, num_path_effect_category[0], EffPathScen, num_path_effect_category[1], effect_category::effective);
							if (LB > objVal)
								match = "NO**";
							else
								match = "NO";

							IsIneffFailed = IloTrue;
						}
					}

				}
				else {
					if (Output[w2].PathStatus == effect_category::unknown) {
						StoreEff(w2, Output[w2], IneffPathScen, num_path_effect_category[0], effect_category::ineffective, IsCondEff);
						num_path_effect_category[2]--;
					}
					else {
						if (Output[w2].PathStatus != effect_category::ineffective) {
							SwapEff(w2, Output[w2], EffPathScen, num_path_effect_category[1], IneffPathScen, num_path_effect_category[0], effect_category::ineffective);
							match = "NO";
							IsEffFailed = IloTrue;
						}
					}
				}
			}
			else {
				if (Output[w2].PathStatus == effect_category::unknown) {
					StoreEff(w2, Output[w2], EffPathScen, num_path_effect_category[1], effect_category::effective, IsCondEff);
					num_path_effect_category[2]--;
				}
				else {
					if (Output[w2].PathStatus != effect_category::effective) {
						SwapEff(w2, Output[w2], IneffPathScen, num_path_effect_category[0], EffPathScen, num_path_effect_category[1], effect_category::effective);
						match = "NO";
						IsIneffFailed = IloTrue;
					}
				}
			}

			Effectiveness << w2 << '\t' << objVal << '\t' << Output[w2].PathObjCost << '\t' << Output[w2].PathStatus << '\t' << match << endl;
		}

	}
	Effectiveness << num_path_effect_category << endl;
	if (IsEffFailed) {
		Effectiveness << "*****Easy-to-Check Conditions for Effective Scenario Paths does NOT match with Resolving*****" << endl;
	}
	if (IsIneffFailed) {
		Effectiveness << "*****Easy-to-Check Conditions for Ineffective Scenario Paths does NOT match with Resolving*****" << endl;
	}
	Summary << num_path_effect_category[0] << '\t' << num_path_effect_category[1] << endl;
	cout << "*****DONE with Resolving Scenario Paths*****" << endl;
}


void Nested(ostream& Nestedness, vector<int> * EffPathScen)
{
	int ii = 0;
	vector<int>::iterator it;
	vector<int>::iterator match; //if found the item in the next ii
	int nestedness = 1;


	Nestedness << "Nestedness of Effective Scenario Paths" << endl;
	for (ii = 20; ii>0; ii--) {
		if (EffPathScen[ii].size() > 0) {
			nestedness = 1;
			for (it = EffPathScen[ii].begin(); it != EffPathScen[ii].end(); it++) {
				match = std::find(EffPathScen[ii - 1].begin(), EffPathScen[ii - 1].end(), *it);
				if (match == EffPathScen[ii - 1].end()) {
					Nestedness << "Found " << *it << " in gamma " << (double)ii / 10 << " not exits in gamma " << (double)(ii - 1) / 10 << endl;
					if (nestedness == 1)
						nestedness = 0;
				}
			}
			if (nestedness == 1) {
				Nestedness << "gamma " << (double)ii / 10 << " is nested in gamma " << (double)(ii - 1) / 10 << endl;
			}

			if (EffPathScen[ii].size()>EffPathScen[ii - 1].size()) {
				Nestedness << "********** gamma " << (double)ii / 10 << " is larger than " << (double)(ii - 1) / 10 << endl;
			}

		}
	}

	cout << "*****DONE with Checking Nestedness of Scenario Paths*****" << endl;
}