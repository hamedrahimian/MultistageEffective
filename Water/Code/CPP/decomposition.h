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

#ifndef SOURCE_H
#define SOURCE_H

#include "types.h"


ILOSTLBEGIN


/** Obtains cdf*/
int CalculateCVaRno(const double beta);

/** Obtains VaR and its index*/
double CalculateCVaR(const double * cost, const double beta);

double CalculateExp(const double * cost);

void FindIndex(const IloInt& stage, const IloInt& nCut, CostsIx& costsix);

void FindWorstProb(const IloInt& stage, const IloInt& nCut, IloNumArray& p, IloNumArray& subObj_hat,
	const pair <double, double>& param);

void getnumTheta(void);

void setThetaName(IloNumVarArray& Theta);

//void SetupFiles(ostream& result, ostream& effectiveness, const IloNum rho);

void Initialization(const IloEnv& env, Cuts& cuts, Duals& duals, CutCoeffs& cutcoeffs, incSols& incsols, optSols& optsols, CostsIx& costsix);

void Termination(ModelArray& MODEL, Cuts& cuts, Duals& duals, CutCoeffs& cutcoeffs, incSols& incsols, optSols& optsols);

void addThetaAlphatoOBJ(const IloInt& stage, 
	IloExpr& OBJ, appxVariables& appxvariable, 
	const pair <double, double>& param);

void subtractThetaAlphafromsubObj(const IloInt& stage, IloCplex& CPX, const appxVariables& appxvariable, 
	IloNum& subObj,
	const pair <double, double>& param);

void createMaster(ModelArray& mod, CplexArray& CPX, Constraints& constraints,
	mainVariables& mainvariables, appxVariables& appxvariables, 
	const pair <double, double>& param);

void createSub(IloCplex& CplexSub, IloModel& modSub, DistGenFormulation& distgen, const pair <double, double> param);

void ForwardPassX(const IloInt& nCut, CplexArray& CPX, Constraints& constraints, Cuts& cuts, 
	const mainVariables& mainvariables, const appxVariables& appxvariables,
	Duals& duals, const CutCoeffs& cutscoeffs, incSols& incsols, CostsIx& costsix,
	IloNum& LB,
	const pair <double, double>& param, 
	const IloBool IsCondEff = IloFalse, const IloInt begin_stage = 0, const IloInt w_ancestor = 0, const IloNumArray Best_y_grand_ancestor = 0);

void ForwardPassP(const IloInt& nCut, IloCplex& CplexSub, IloModel& modSub, DistGenFormulation& distgen,
	IloNumArray& p, IloNumArray & subObj_hat, 
	IloBool& contForwardPass, 
	const pair <double, double> param, 
	const IloBool IsCondEff = IloFalse, const IloBool IsPathEff = IloFalse, const IloInt begin_stage = 0, const IloInt w_ancestor = 0, const IloInt child = 0);

void ForwardPassP2(const IloInt& nCut, 
	IloNumArray& p, IloNumArray & subObj_hat,
	const pair <double, double> param,
	const IloBool IsCondEff = IloFalse, const IloBool IsPathEff = IloFalse, const IloInt begin_stage = 0, const IloInt w_ancestor = 0, const IloInt child = 0);

void CalculateFuncVal(const IloInt& nCut, CostsIx& costsix,
	const pair<double, double> param);

void BackwardPass(const IloInt& nCut, ModelArray& MODEL, CplexArray& CPX,
	Constraints& constraints, Cuts& cuts, const IloNumVarArray& Y, const appxVariables& appxvariables,
	Duals& duals, CutCoeffs& cutscoeffs, const incSols& incsols, CostsIx& costsix,
	const pair <double, double>& param, 
	const IloInt begin_stage = 0, const IloInt w_ancestor = 0);

void InitializeCutsnadCoeffs(Cuts& cuts, IloNumArray&  pi_optCut, CutCoeffs& cutcoeffs, 
	const IloInt begin_stage = 0);

void InitializeIndex(CostsIx& costsix);

void AddEmptyCuts(const IloInt& stage, const IloInt& nCut, IloModel& MODEL, const Cuts& cuts, const IloInt begin_stage = 0);

void getConstraintsDuals(const IloInt& stage, const IloInt& scen, const IloInt& nCut, 
	IloCplex& CPX, const Constraints& constraints, Duals& duals);

void getoptCutsDuals(const IloInt& stage, const IloInt& scen, const IloInt& nCut, 
	IloCplex& CPX, const Cuts& cuts, IloNumArray& pi_optCut,
	const IloInt begin_stage = 0);

void GenerateCutCoeffs(const IloInt& stage, const IloInt& scenario, const IloInt& nCut,
	const Duals& duals, CutCoeffs& cutcoeffs, const incSols& incsols, const CostsIx& costsix,
	const pair<double, double>& param, 
	const IloInt begin_stage = 0);

void GenerateCuts(const IloInt& stage, const IloInt& scenario, const IloInt& nCut,
	Cuts& cuts,
	const IloNumVarArray& Y, const IloNumVar& Alpha, const IloNumVarArray& Theta,
	const CutCoeffs& cutcoeffs, const incSols& incsols, const IloInt& varIndex, const IloInt& supIndex,
	const pair<double, double>& param,
	const IloInt begin_stage = 0);

void UpdateOldCuts(const IloInt& stage, const IloInt& scenario, const IloInt& nCut, 
	Cuts& cuts,	
	const IloNumVarArray& Y, const IloNumVar& Alpha, const IloNumVarArray& Theta,
	const CutCoeffs& cutcoeffs, const incSols& incsols, const CostsIx& costsix,
	const pair<double, double> param, 
	const IloInt begin_stage = 0);

void UpdateRhs(const IloInt& stage, const IloInt& scenario, const IloInt& nCut, Constraints& constraints, const IloNumArray& y);

void SolveOriginal(ostream& Result, ModelArray MODEL, CplexArray CPX, Formulation formulation, 
	IloCplex CplexSub, IloModel modSub, DistGenFormulation distgen,
	IloNum& objVal,
	const pair<double, double> param, 
	Scenario * Output);

void usage(char *progname);

#endif