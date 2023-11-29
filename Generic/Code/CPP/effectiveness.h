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


#ifndef EFFECTIVENESS_H
#define EFFECTIVENESS_H

#include <ilcplex/ilocplex.h>
#include "decomposition.h"
#include "scenariotree.h"
#include <stdlib.h>   
#include <vector>

/** Enumaration for primal category of a scenario*/
enum   cost_category {
	b_VaR = 1,
	VaR,
	a_VaR,
	sup
};

/** Enumaration for the effectiveness of a scenario*/
enum  effect_category {
	ineffective,
	effective,
	unknown
};


/**Type to store the informatio about a scenario node in order to be used in easy-to-check conditions*/
typedef struct ScenarioChildInfo {
	int ChildIndex; //index over its sibilings from the same parent
	int No; //overall index
	double Cost; //cost function evaluation of a scenario at optimal policy
	double Nominal_Prob; //nominal probabiliy
	double Worst_Prob; //worst-case probability
	double cdf; //conditional cdf among the siblings of this scenario
	double CondObjCost; //objective function value after this scenario is removed conditionally
	double PathObjCost; //objective function value after this scenario path is removed (only for last stage scenario nodes)
	effect_category CondStatus; //status of the conditional effectiveness
	effect_category PathStatus; //status of the effectiveness of scenario path (only for last stage scenario nodes)
	cost_category Primal_Category; //(conditional) primal categroy of this scenario 
} ScenChildResult;

/**Type to store the informatio about children of a node in order to be used in easy-to-check conditions*/
typedef struct OutputInfo {
	int numChild;
	vector<ScenChildResult> Child;
	//ScenChildResult* Child;
} OutputType;

/**Type to compare scenarios base on cost, and using their indices as a tie breaker*/
struct CompareScenariosCostIndex
{
	bool operator()(const ScenChildResult& lhs, const ScenChildResult& rhs)  const
	{
		if (lhs.Cost == rhs.Cost)
			return lhs.No < rhs.No;
		else
			return lhs.Cost < rhs.Cost;
	}
};

/**Type to compare scenarios base on their indices*/
struct CompareScenariosIndex
{
	bool operator()(const ScenChildResult& lhs, const ScenChildResult& rhs)  const
	{
		return lhs.No < rhs.No;
	}
};

class Effectiveness : public Decomposition {
	public:

		/** Obtains cdf*/
		void CalculateCDF(const int numChild, const double* nominal_prob, double* cdf);

		/** Obtains VaR and its index*/
		void CalculateVaR(const int numChild, const double* cdf, const double * cost, int * indexVaR, double * VaR, const double beta);

		/** Stores the effectiveness of a scenario*/
		void StoreEff(const int it, ScenChildResult * scen, vector<int>* eff, IloInt* num_effect_category, const effect_category status, const bool IsCondEff=true);

		/** Swap the effectiveness of a scenarip path in case that resolvin does not match with easy to check conditions */
		void SwapEff(const int it, ScenChildResult * scen, vector<int>* eff_old, IloInt* num_effect_category_old, vector<int>* eff_new, IloInt *num_effect_category_new, const effect_category status_new);

		/** Prints information about the conditional effectiveness of scenarios*/
		void PrintCondEff(ostream& fe, const int w, const int period,
			OutputType* wOutput, const IloIntArray* num_effect_category);

		/** Checks nestedness of effective scenario paths*/
		void CheckNestednessPath(char *problem_name, vector<int>*** EffPathScen, int numInstance, const double * rho);

		/** Checks nestedness of condtionally effective scenario nodes*/
		void CheckNestednessCond(char *problem_name, vector<int>*** EffCondScen, int numInstance, const double * rho);

};


#endif