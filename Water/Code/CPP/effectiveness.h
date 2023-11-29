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

#ifndef EFFECTIVENESS_H
#define EFFECTIVENESS_H

#include "types.h"

using namespace std;

void FormPrimalCategories(const double * cost, vector<int>** categories);

void EasyCheckCondEff(ostream& Effectiveness, Scenario * Output, vector<int>& IneffCondScen, vector<int>& EffCondScen, vector<int>& UnknownCondScen, const pair<double, double> param, IloIntArray2& num_cost_category, IloIntArray2& num_cond_effect_category);

void EasyCheckCondEffReal(const IloInt& scenario, const IloInt& nCut, const IloNumArray& subObj_hat, const IloNumArray& p,
	effNodeType* EffCondScen, effNodeType* IneffCondScen, effNodeType* UnknownCondScen, const pair<double, double> param);

void EasyCheckPathEff(ostream& Summary, ostream& Effectiveness, Scenario * Output, vector<int>& IneffPathScen, vector<int>& EffPathScen, vector<int>& UnknownPathScen, const pair<double, double> param, IloIntArray & num_path_effect_category);

void EasyCheckEffPath(const IloInt& nCut, const IloNumArray& p, vector<effpathNode>* patheffs);

void ResolveCondEff(ostream& Effectiveness, ModelArray MODEL, CplexArray CPX, Formulation formulation,
	IloCplex CplexSub, IloModel modSub, DistGenFormulation distgen,
	const pair<double, double> param,
	Scenario * Output, vector<int>& IneffCondScen, vector<int>& EffCondScen, vector<int>& UnknownCondScen, IloIntArray2& num_cost_category, IloIntArray2& num_cond_effect_category);

void ResolvePathEff(ostream& Summary, ostream& Effectiveness,
	ModelArray MODEL, CplexArray CPX, Formulation formulation,
	IloCplex CplexSub, IloModel modSub, DistGenFormulation distgen,
	const IloNum objVal,
	const pair<double, double> param,
	Scenario * Output, vector<int>& IneffPathScen, vector<int>& EffPathScen, vector<int>& UnknownPathScen, IloIntArray& num_path_effect_category);

void Nested(ostream& Nestedness, vector<int> * EffPathScen);

void StoreEff(const int it, Scenario& scen, vector<int>& eff, IloInt& num_effect_category, const effect_category& status, const IloBool IsCondEff = IloTrue);

void SwapEff(const int it, Scenario & scen, vector<int>& eff_old, IloInt & num_effect_category_old, vector<int>& eff_new, IloInt & num_effect_category_new, const effect_category & status_new);

#endif
