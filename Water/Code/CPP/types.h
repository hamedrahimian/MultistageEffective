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

#ifndef TYPES_H
#define TYPES_H


#include <ilcplex/ilocplex.h>   
#include <stdlib.h>
#include <vector>

enum  cost_category { b_VaR = 1, VaR, a_VaR, sup };

enum  effect_category { ineffective, effective, unknown };

enum ProblemType {
	RISK_NEUTRAL = 0,
	PRIMAL, //this means DRSO problem
	DUAL //this means risk-averse problem
};

enum CutType {
	SINGLE = 0,
	MULTI
};

enum DecompType {
	D1 = 0,
	D2,
	D3,
	D4,
	D5
};

enum SeparationType {
	COMBINED = 0,
	SEPARATED
};

enum AmbiguityType {
	TV = 0,
	EC
};

enum struct VariantType : char {
	EC_P_M_C,
	EC_D_S_C,
	EC_D_M_C,
	EC_D_S_S_S,
	EC_D_S_M_S,
	EC_D_M_S_S,
	EC_D_M_M_S,
	TV_P_M_C,
	TV_D_S_C,
	TV_D_M_C,
	TV_D_S_S_S,
	TV_D_S_M_S,
	TV_D_M_S_S,
	TV_D_M_M_S
};

struct Variant {
	AmbiguityType ambiguity_type;
	ProblemType problem_type;
	CutType cut1_type, cut2_type;
	DecompType decomp_type;
	SeparationType separation_type;
};

struct Scenario {
	IloNum partial_sum;
	IloInt No;
	IloNum Cost;
	IloNum nom_Prob;
	IloNum worst_Prob;
	IloNum CondObjCost;
	IloNum PathObjCost;
	effect_category CondStatus;
	effect_category PathStatus;
	cost_category primal_category;
};

struct NodeInfo {
	IloInt No; //index among sibilings
	IloNum Cost;
	IloNum worst_Prob;
};

/** Types to store conditional effectiveness */
typedef std::vector<int> effNodeType; //vetor for elements 
typedef std::vector<effNodeType> effNode; //vector for cut

/** Types to store path effectiveness */
typedef std::vector<effect_category> effpathNode; //vetor for ncut

struct CompareScenarios
{
	bool operator()(const Scenario& lhs, const Scenario& rhs)  const
	{
		if (lhs.Cost == rhs.Cost)
			return lhs.No < rhs.No;
		else
			return lhs.Cost < rhs.Cost;
	}
};

struct CompareScenariosCostIndex
{
	bool operator()(const NodeInfo& lhs, const NodeInfo& rhs)  const
	{
		if (lhs.Cost == rhs.Cost)
			return lhs.No < rhs.No;
		else
			return lhs.Cost < rhs.Cost;
	}
};

struct CompareScenarios2
{
	bool operator()(const Scenario& lhs, const Scenario& rhs)  const
	{
		return lhs.No < rhs.No;
	}
};

/** Type to create an array of IloNumVarArray objects, array for periods */
typedef IloArray<IloNumVarArray> VarArray2;
/** Type to create an array of IloRangeArray objects, array for periods */
typedef IloArray<IloRangeArray>  EquationArray;
/** Type to create an array of IloCplex objects, array for periods */
typedef IloArray<IloCplex>       CplexArray;
/** Type to create an array of IloModel objects, array for periods */
typedef IloArray<IloModel>       ModelArray;
/** Type to create an array of IloExpr objects, array for periods */
typedef IloArray<IloExpr>   ObjArray;

struct Constraints {
	IloRangeArray FlowBalance;
	IloRangeArray CapBound;
	IloRangeArray MeetDemand;
	IloRangeArray ReturnFlow;
	IloRangeArray SafeYield;
	IloRangeArray StorageBalance;
	IloRangeArray RFoutflowBound;
};

struct Cuts {
	IloRangeArray optCut;
	IloRangeArray maxCut;
};

struct mainVariables {
	IloNumVarArray Q;
	IloNumVarArray Y;
};

struct appxVariables {
	IloNumVarArray Theta;
	IloNumVarArray Alpha;
};

struct CutCoeffs {
	IloNumArray g;
	IloNumArray G;
	IloNumArray psi;
	IloNumArray psi_x;
	IloNumArray psi_eta;
};

struct Duals {
	IloNumArray pi_balance;
	IloNumArray pi_capacity;
	IloNumArray pi_demand;
	IloNumArray pi_return;
	IloNumArray pi_safe;
	IloNumArray pi_storage;
	IloNumArray pi_flowbound;
	IloNumArray pi_optCut;
};

struct incSols {
	IloNumArray y;
	IloNumArray p;
};

struct optSols {
	IloNumArray Best_y;
	IloNumArray Worst_p;
};

struct CostsIx {
	IloNumArray subObj_hat;
	IloIntArray varIndex;
	IloIntArray supIndex;
};

struct Formulation {
	Constraints constraints;
	Cuts cuts;
	mainVariables mainvariables;
	appxVariables appxvariables;
	CutCoeffs cutcoeffs;
	Duals duals;
	incSols incsols;
	optSols optsols;
	CostsIx costsix;
};

struct DistGenFormulation {
	IloNumVarArray P;
	IloNumVarArray Z;
	IloObjective objSub;
	IloRangeArray conDistancePos;
	IloRangeArray conDistanceNeg;
};




#endif