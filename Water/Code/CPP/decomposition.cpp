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

#include "decomposition.h"
#include "effectiveness.h"
#include "vars.h"
#include "tree.h"
#include <algorithm> 


using namespace std;

#define toler		0.00000001


int CalculateCVaRno(const double beta){
	int j;
	double q = prob[1];
	if (q >= beta)
		CVaRno = 0;
	else {
		for (j = 1; j < numScen; j++) {
			if ( (j+1)*q >= beta - toler && j*q < beta) {
				CVaRno = j;
				break;
			}
		}
	}
	return CVaRno;

	
}// endCalculateCVaRno

double CalculateCVaR(const double * cost, const double beta){
	int j;
	double CVaR=0;
	double prob = (double)1 / numScen;

	CVaR = cost[CVaRno];
	for (j = CVaRno + 1; j < numScen; j++) {				
		CVaR += (cost[j] - cost[CVaRno])*prob / (1 - beta);
	}

	return CVaR;
}//end CalculateCVaR

double CalculateExp(const double * cost){
	int j;
	double Exp = 0;
	double prob = (double)1 / numScen;
	for (j = 0; j < numScen; j++) 
		Exp += cost[j] *prob ;
	
	return	 Exp;
}//end CalculateExp

void FindIndex(const IloInt& stage, const IloInt& nCut, CostsIx& costsix) {
	int w, j;

	int numParents = numScenNode_stage_sum[numStage - 2];

	if (stage-1 == 0)
		w = 0;
	else
		w = numScenNode_stage_sum[stage - 2];
	
	while (w < numScenNode_stage_sum[stage-1]) {
		NodeInfo * wOutput = new NodeInfo[numScen];
		for (j = 0; j < numScen; j++) {
			wOutput[j].No = j;
			wOutput[j].Cost = costsix.subObj_hat[descendant[w][j]];
		}
		sort(wOutput, wOutput + numScen, CompareScenariosCostIndex());
		int varix = wOutput[CVaRno].No;
		int supix = wOutput[numScen - 1].No;
		costsix.varIndex[nCut*numParents + w] = varix;
		costsix.supIndex[nCut*numParents + w] = supix;

		delete[] wOutput;
		w++;
	}
			
}//end FindIndex

void FindWorstProb(const IloInt& stage, const IloInt& nCut, IloNumArray& p, IloNumArray& subObj_hat,
	const pair <double, double>& param) {
	int w, j;

	int numParents = numScenNode_stage_sum[numStage - 2];
	
	/*IloEnv env;
	IloNumArray pp(env, numScen);*/

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;

	if (stage - 1 == 0)
		w = 0;
	else
		w = numScenNode_stage_sum[stage - 2];

	while (w < numScenNode_stage_sum[stage - 1]) {
		NodeInfo * wOutput = new NodeInfo[numScen];
		for (j = 0; j < numScen; j++) {
			wOutput[j].No = j;
			wOutput[j].Cost = subObj_hat[descendant[w][j]];
		}
		sort(wOutput, wOutput + numScen, CompareScenariosCostIndex());
		
		int varix = wOutput[CVaRno].No;
		double var_value = wOutput[CVaRno].Cost;
		double sup_value = wOutput[numScen-1].Cost;
		double lambda = abs(var_value - sup_value);
		
		if (variant.ambiguity_type == TV) {
			/*for (j = 0; j < numScen; j++)
				distgen.objSub.setLinearCoef(distgen.P[j], subObj_hat[descendant[w][j]]);

			for (j = 0; j < numScen; j++) {
				distgen.conDistancePos[j].setUb(prob[descendant[w][j]]);
				distgen.conDistanceNeg[j].setUb(-prob[descendant[w][j]]);
			}
			
			CplexSub.extract(modSub);
			CplexSub.solve();
			if (CplexSub.solve()) {
				for (j = 0; j < numScen; j++) {
					pp[j] = CplexSub.getValue(distgen.P[j]);
					cout << pp[j] << endl;
				}
				
			}*/
			if (lambda > toler) {
				for (j = 0; j < CVaRno; j++)
					p[nCut*numScenNode + descendant[w][wOutput[j].No]] = 0;

				p[nCut*numScenNode + descendant[w][wOutput[CVaRno].No]] = (CVaRno + 1)*prob[descendant[w][wOutput[CVaRno].No]] - rho;

				for (j = CVaRno + 1; j < numScen - 1; j++)
					p[nCut*numScenNode + descendant[w][wOutput[j].No]] = prob[descendant[w][wOutput[j].No]];

				p[nCut*numScenNode + descendant[w][wOutput[numScen - 1].No]] = prob[descendant[w][wOutput[numScen - 1].No]] + rho;
			}
			else {
				for (j = 0; j < CVaRno; j++)
					p[nCut*numScenNode + descendant[w][wOutput[j].No]] = 0;

				p[nCut*numScenNode + descendant[w][wOutput[CVaRno].No]] = 1;
			}
			/*for (j = 0; j < numScen; j++)
				if (abs(p[nCut*numScenNode + descendant[w][j]] - pp[j])>toler)
					cout << "not matching" << endl;*/
			

		}
		else if (variant.ambiguity_type == EC)
		{
			for (j = 0; j < numScen; j++) {
				p[nCut*numScenNode + descendant[w][j]] = (1 - lambda)*prob[descendant[w][j]];
			}
			for (j = CVaRno + 1; j < numScen; j++) {
				p[nCut*numScenNode + descendant[w][wOutput[j].No]] += lambda*prob[descendant[w][wOutput[j].No]] / (1 - alpha);
				
			}
			double remained_prob = 1;
			for (j = 0; j < numScen; j++) {
				remained_prob -= p[nCut*numScenNode + descendant[w][j]];
			}
			p[nCut*numScenNode + descendant[w][varix]] += remained_prob;

		}

		delete[] wOutput;
		w++;
	}

}//end FindWorstProb

void getnumTheta(void) {
	switch (types) {
		case  VariantType::EC_P_M_C:
			total_theta = (numStage - 1)*numScen; break;
		case  VariantType::EC_D_S_C:
			total_theta = numStage - 1; break;
		case  VariantType::EC_D_M_C:
			total_theta = (numStage - 1)*numScen; break;
		case  VariantType::EC_D_S_S_S:
			total_theta = (numStage - 1) * 2; break;
		case  VariantType::EC_D_S_M_S:
			total_theta = (numStage - 1) + (numStage - 1)*numScen; break;
		case  VariantType::EC_D_M_S_S:
			total_theta = (numStage - 1)*numScen + numStage - 1; break;
		case  VariantType::EC_D_M_M_S:
			total_theta = (numStage - 1)*numScen * 2; break;
		case  VariantType::TV_P_M_C:
			total_theta = (numStage - 1)*numScen; break;
		case  VariantType::TV_D_S_C:
			total_theta = numStage - 1; break;
		case  VariantType::TV_D_M_C:
			total_theta = (numStage - 1)*numScen; break;
		case  VariantType::TV_D_S_S_S:
			total_theta = (numStage - 1) * 2; break;
		case  VariantType::TV_D_S_M_S:
			total_theta = (numStage - 1) + (numStage - 1)*numScen; break;
		case  VariantType::TV_D_M_S_S:
			total_theta = (numStage - 1) * 2; break;
		case  VariantType::TV_D_M_M_S:
			total_theta = (numStage - 1) + (numStage - 1)*numScen; break;
	}
}//getnumTheta

void setThetaName(IloNumVarArray& Theta) {
	int s, w;
	char varName[100];

	switch (types) {
	case  VariantType::EC_P_M_C:
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "Theta.%d.%d", (int)s, (int)w);
				Theta[s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::EC_D_S_C:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "Theta.%d", (int)s);
			Theta[s].setName(varName);
		}
		break;
	case  VariantType::EC_D_M_C:
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "Theta.%d.%d", (int)s, (int)w);
				Theta[s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::EC_D_S_S_S:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaExp.%d", (int)s);
			Theta[s].setName(varName);
		}
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaCVaR.%d", (int)s);
			Theta[numStage-1 + s].setName(varName);
		}
		break;
	case  VariantType::EC_D_S_M_S:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaExp.%d", (int)s);
			Theta[s].setName(varName);
		}
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "ThetaCVaR.%d.%d", (int)s, (int)w);
				Theta[numStage - 1 + s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::EC_D_M_S_S:
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "ThetaExp.%d.%d", (int)s, (int)w);
				Theta[s*numScen + w].setName(varName);
			}
		}
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaCVaR.%d", (int)s);
			Theta[(numStage - 1)*numScen + s].setName(varName);
		}
		break;
	case  VariantType::EC_D_M_M_S:
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "ThetaExp.%d.%d", (int)s, (int)w);
				Theta[s*numScen + w].setName(varName);
			}
		}
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "ThetaCVaR.%d.%d", (int)s, (int)w);
				Theta[(numStage - 1)*numScen+ s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::TV_P_M_C:
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "Theta.%d.%d", (int)s, (int)w);
				Theta[s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::TV_D_S_C:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "Theta.%d", (int)s);
			Theta[s].setName(varName);
		}
		break;
	case  VariantType::TV_D_M_C:
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "Theta.%d.%d", (int)s, (int)w);
				Theta[s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::TV_D_S_S_S:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaSup.%d", (int)s);
			Theta[s].setName(varName);
		}
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaCVaR.%d", (int)s);
			Theta[numStage - 1 + s].setName(varName);
		}
		break;
	case  VariantType::TV_D_S_M_S:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaSup.%d", (int)s);
			Theta[s].setName(varName);
		}
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "ThetaCVaR.%d.%d", (int)s, (int)w);
				Theta[numStage - 1 + s*numScen + w].setName(varName);
			}
		}
		break;
	case  VariantType::TV_D_M_S_S:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaSup.%d", (int)s);
			Theta[s].setName(varName);
		}
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaCVaR.%d", (int)s);
			Theta[numStage - 1 + s].setName(varName);
		}
		break;
	case  VariantType::TV_D_M_M_S:
		for (s = 0; s<numStage - 1; s++) {
			sprintf_s(varName, "ThetaSup.%d", (int)s);
			Theta[s].setName(varName);
		}
		for (s = 0; s<numStage - 1; s++) {
			for (w = 0; w<numScen; w++) {
				sprintf_s(varName, "ThetaCVaR.%d.%d", (int)s, (int)w);
				Theta[(numStage - 1) + s*numScen + w].setName(varName);
			}
		}
		break;
	}
}//setThetaName

void Initialization(const IloEnv& env, Cuts& cuts, Duals& duals, CutCoeffs& cutcoeffs, incSols& incsols, optSols& optsols, CostsIx& costsix)
{
	IloInt w;
	//IloEnv env;

	cuts.maxCut = IloRangeArray(env);
	cuts.optCut = IloRangeArray(env);
		
	//duals.pi_balance=IloNumArray (env, numBalance*total_pi);
	duals.pi_capacity=IloNumArray (env, numCapacity*total_pi);
	duals.pi_demand=IloNumArray (env, numUser* total_pi);
	duals.pi_return=IloNumArray (env, numPotUser*total_pi);
	duals.pi_safe=IloNumArray (env, numRecharge*total_pi);
	duals.pi_storage=IloNumArray (env, numRecharge*numScenNode);
	duals.pi_flowbound=IloNumArray (env, numRecharge*numScenNode);
	duals.pi_optCut = IloNumArray(env);

	cutcoeffs.g= IloNumArray(env);
	cutcoeffs.G = IloNumArray(env);
	cutcoeffs.psi = IloNumArray(env);
	cutcoeffs.psi_x = IloNumArray(env);
	cutcoeffs.psi_eta = IloNumArray(env);

	incsols.y = IloNumArray(env);
	incsols.p = IloNumArray(env);
	optsols.Best_y = IloNumArray(env, numRecharge*numScenNode);
	optsols.Worst_p = IloNumArray(env, numScenNode);
	
	costsix.subObj_hat=IloNumArray (env, numScenNode);
	costsix.varIndex = IloIntArray(env);
	costsix.supIndex = IloIntArray(env);
	
}

void Termination(ModelArray& MODEL, Cuts& cuts, Duals& duals, CutCoeffs& cutcoeffs, incSols& incsols, optSols& optsols)
{
	int s;
	incsols.y.clear(); incsols.p.clear();
	optsols.Best_y.clear(); optsols.Worst_p.clear();
	cutcoeffs.g.clear(); cutcoeffs.G.clear(); cutcoeffs.psi.clear(); cutcoeffs.psi_x.clear(); cutcoeffs.psi_eta.clear();
	//duals.pi_balance.clear(); 
	duals.pi_return.clear(); duals.pi_safe.clear(); duals.pi_storage.clear(); duals.pi_flowbound.clear(); 
	duals.pi_capacity.clear(); duals.pi_demand.clear(); duals.pi_optCut.clear();
	for (s = 0; s < numStage - 1; s++) {
		MODEL[s].remove(cuts.optCut);
		MODEL[s].remove(cuts.maxCut);
	}
	cuts.optCut.clear();
	cuts.maxCut.clear();
}

void addThetaAlphatoOBJ(const IloInt& stage,
	IloExpr& OBJ, appxVariables& appxvariable,
	const pair <double, double>& param) {

	int w;
	
	double prob = (double)1 / numScen;

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;

	switch (types) {
	case  VariantType::EC_P_M_C:
		OBJ += appxvariable.Alpha[stage]; break;
	case  VariantType::EC_D_S_C:
		OBJ += appxvariable.Theta[stage]; break;
	case  VariantType::EC_D_M_C:
		for (w = 0; w<numScen; w++) {
			OBJ += appxvariable.Theta[stage*numScen+w]*prob;
		}
		break;
	case  VariantType::EC_D_S_S_S:
		OBJ += (1-lambda)*appxvariable.Theta[stage];
		OBJ += lambda*appxvariable.Theta[numStage - 1 + stage];
		break;
	case  VariantType::EC_D_S_M_S:
		OBJ += (1 - lambda)*appxvariable.Theta[stage];
		for (w = 0; w<numScen; w++) {
			OBJ += lambda*appxvariable.Theta[numStage - 1 + stage*numScen + w]*prob;
		}
		break;
	case  VariantType::EC_D_M_S_S:
		for (w = 0; w<numScen; w++) {
			OBJ += (1 - lambda)*appxvariable.Theta[stage*numScen + w] * prob;
		}
		OBJ += lambda*appxvariable.Theta[(numStage - 1)*numScen + stage];
		break;
	case  VariantType::EC_D_M_M_S:
		for (w = 0; w<numScen; w++) {
			OBJ += (1 - lambda)*appxvariable.Theta[stage*numScen + w] * prob;
		}
		for (w = 0; w<numScen; w++) {
			OBJ+=lambda*appxvariable.Theta[(numStage - 1)*numScen + stage*numScen + w];
		}
		break;
	case  VariantType::TV_P_M_C:
		OBJ += appxvariable.Alpha[stage]; break;
	case  VariantType::TV_D_S_C:
		OBJ += appxvariable.Theta[stage]; break;
	case  VariantType::TV_D_M_C:
		for (w = 0; w<numScen; w++) {
			OBJ += appxvariable.Theta[stage*numScen + w] * prob;
		}
		break;
	case  VariantType::TV_D_S_S_S:
		OBJ += rho*appxvariable.Theta[stage];
		OBJ += (1-rho)*appxvariable.Theta[numStage - 1 + stage];
		break;
	case  VariantType::TV_D_S_M_S:
		OBJ += rho*appxvariable.Theta[stage];
		for (w = 0; w<numScen; w++) {
			OBJ += (1 - rho)*appxvariable.Theta[numStage - 1 + stage*numScen + w] * prob;
		}
		break;
	case  VariantType::TV_D_M_S_S:
		OBJ += rho*appxvariable.Theta[stage];
		OBJ += (1 - rho)*appxvariable.Theta[numStage - 1 + stage];
		break;
	case  VariantType::TV_D_M_M_S:
		OBJ += rho*appxvariable.Theta[stage];
		for (w = 0; w<numScen; w++) {
			OBJ += (1 - rho)*appxvariable.Theta[numStage - 1 + stage*numScen + w]*prob;
		}
		break;
	}

}//addThetaAlphatoOBJ

void subtractThetaAlphafromsubObj(const IloInt& stage, IloCplex& CPX, const appxVariables& appxvariable,
	IloNum& subObj,
	const pair <double, double>& param) {

	int w;
	
	double prob = (double)1 / numScen;

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;

	switch (types) {
	case  VariantType::EC_P_M_C:
		subObj -= CPX.getValue(appxvariable.Alpha[stage]); break;
	case  VariantType::EC_D_S_C:
		subObj -= CPX.getValue(appxvariable.Theta[stage]); break;
	case  VariantType::EC_D_M_C:
		for (w = 0; w<numScen; w++) {
			subObj -= CPX.getValue(appxvariable.Theta[stage*numScen + w]) * prob;
		}
		break;
	case  VariantType::EC_D_S_S_S:
		subObj -= (1-lambda)*CPX.getValue(appxvariable.Theta[stage]);
		subObj -= lambda*CPX.getValue(appxvariable.Theta[numStage - 1 + stage]);
		break;
	case  VariantType::EC_D_S_M_S:
		subObj -= (1 - lambda)*CPX.getValue(appxvariable.Theta[stage]);
		for (w = 0; w<numScen; w++) {
			subObj -= lambda*CPX.getValue(appxvariable.Theta[numStage - 1 + stage*numScen + w]) * prob;
		}
		break;
	case  VariantType::EC_D_M_S_S:
		for (w = 0; w<numScen; w++) {
			subObj -= (1 - lambda)*CPX.getValue(appxvariable.Theta[stage*numScen + w]) * prob;
		}
		subObj -= lambda*CPX.getValue(appxvariable.Theta[(numStage - 1)*numScen + stage]);
		break;
	case  VariantType::EC_D_M_M_S:
		for (w = 0; w<numScen; w++) {
			subObj -= (1 - lambda)*CPX.getValue(appxvariable.Theta[stage*numScen + w]) * prob;
		}
		for (w = 0; w<numScen; w++) {
			subObj -= lambda*CPX.getValue(appxvariable.Theta[(numStage - 1)*numScen + stage*numScen + w]);
		}
		break;
	case  VariantType::TV_P_M_C:
		subObj -= CPX.getValue(appxvariable.Alpha[stage]); break;
	case  VariantType::TV_D_S_C:
		subObj -= CPX.getValue(appxvariable.Theta[stage]); break;
	case  VariantType::TV_D_M_C:
		for (w = 0; w<numScen; w++) {
			subObj -= CPX.getValue(appxvariable.Theta[stage*numScen + w]) * prob;
		}
		break;
	case  VariantType::TV_D_S_S_S:
		subObj -= rho*CPX.getValue(appxvariable.Theta[stage]);
		subObj -= (1-rho)*CPX.getValue(appxvariable.Theta[numStage - 1 + stage]);
		break;
	case  VariantType::TV_D_S_M_S:
		subObj -= rho*CPX.getValue(appxvariable.Theta[stage]);
		for (w = 0; w<numScen; w++) {
			subObj -= (1 - rho)*CPX.getValue(appxvariable.Theta[numStage - 1 + stage*numScen + w]) * prob;
		}
		break;
	case  VariantType::TV_D_M_S_S:
		subObj -= rho*CPX.getValue(appxvariable.Theta[stage]);
		subObj -= (1 - rho)*CPX.getValue(appxvariable.Theta[numStage - 1 + stage]);
		break;
	case  VariantType::TV_D_M_M_S:
		subObj -= rho*CPX.getValue(appxvariable.Theta[stage]);
		for (w = 0; w<numScen; w++) {
			subObj -= (1 - rho)*CPX.getValue(appxvariable.Theta[numStage - 1 + stage*numScen + w])*prob;
		}
		break;
	}

}//subtractThetaAlphafromOBJ

void InitializeCutsnadCoeffs(Cuts& cuts, IloNumArray&  pi_optCut, CutCoeffs& cutcoeffs,
	const IloInt begin_stage ) {
	IloEnv env= cutcoeffs.g.getEnv();
	int numMidNodes=numScenNode - 1;
	int numParents = numScenNode_stage_sum[numStage - 2];

	IloNumArray GTemp(env, numRecharge*numMidNodes);
	cutcoeffs.G.add(GTemp);
	GTemp.end();
	IloNumArray gTemp(env, numMidNodes);
	cutcoeffs.g.add(gTemp);
	gTemp.end();
	if (variant.problem_type == DUAL) {
		IloNumArray psi_xTemp(env, numRecharge*numMidNodes);
		cutcoeffs.psi_x.add(psi_xTemp);
		psi_xTemp.end();
		IloNumArray psi_etaTemp(env, numMidNodes);
		cutcoeffs.psi_eta.add(psi_etaTemp);
		psi_etaTemp.end();
		IloNumArray psiTemp(env, numMidNodes);
		cutcoeffs.psi.add(psiTemp);
		psiTemp.end();
	}
	


	switch (types) {
		case  VariantType::EC_P_M_C:
		{
			IloNumArray pi_optCutTemp(env, numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numScen*(numStage - begin_stage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			IloRangeArray maxCutTemp(env, numStage - begin_stage - 1, 0, IloInfinity);
			cuts.maxCut.add(maxCutTemp);
			maxCutTemp.end();
			break;
		}
		
		case  VariantType::EC_D_S_C:
		{
			IloNumArray pi_optCutTemp(env, numParents);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numStage - 1, 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		
		case  VariantType::EC_D_M_C:
		{
			IloNumArray pi_optCutTemp(env, numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numScen*(numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::EC_D_S_S_S:
		{
			IloNumArray pi_optCutTemp(env, 2 * numParents);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, 2 * (numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::EC_D_S_M_S:
		{
			IloNumArray pi_optCutTemp(env, numParents + numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, (numStage - 1) + numScen*(numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::EC_D_M_S_S:
		{
			IloNumArray pi_optCutTemp(env, numMidNodes + numParents);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numScen*(numStage - 1) + (numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::EC_D_M_M_S:
		{
			IloNumArray pi_optCutTemp(env, 2 * numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, 2 * numScen*(numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::TV_P_M_C:
		{
			IloNumArray pi_optCutTemp(env, numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numScen*(numStage - begin_stage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			IloRangeArray maxCutTemp(env, numStage - begin_stage - 1, 0, IloInfinity);
			cuts.maxCut.add(maxCutTemp);
			maxCutTemp.end();
			break;
		}
		case  VariantType::TV_D_S_C:
		{
			IloNumArray pi_optCutTemp(env, numParents);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numStage - 1, 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::TV_D_M_C:
		{
			IloNumArray pi_optCutTemp(env, numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numScen*(numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::TV_D_S_S_S:
		{
			IloNumArray pi_optCutTemp(env, 2 * numParents);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, 2 * (numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::TV_D_S_M_S:
		{
			IloNumArray pi_optCutTemp(env, numParents + numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, (numStage - 1) + numScen*(numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::TV_D_M_S_S:
		{
			IloNumArray pi_optCutTemp(env, numMidNodes + numParents);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, numScen*(numStage - 1) + (numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
		case  VariantType::TV_D_M_M_S:
		{
			IloNumArray pi_optCutTemp(env, 2 * numMidNodes);
			pi_optCut.add(pi_optCutTemp);
			pi_optCutTemp.end();
			IloRangeArray optCutTemp(env, 2 * numScen*(numStage - 1), 0, IloInfinity);
			cuts.optCut.add(optCutTemp);
			optCutTemp.end();
			break;
		}
	}
}//InitializeCutCoeffs

void InitializeIndex(CostsIx& costsix) {
	IloEnv env = costsix.supIndex.getEnv();

	IloIntArray varIndexTemp(env, numScenNode_stage_sum[numStage - 2]);
	costsix.varIndex.add(varIndexTemp);
	varIndexTemp.end();
	IloIntArray supIndexTemp(env, numScenNode_stage_sum[numStage - 2]);
	costsix.supIndex.add(supIndexTemp);
	supIndexTemp.end();
	
}//InitializeIndex

void AddEmptyCuts(const IloInt& stage, const IloInt& nCut, IloModel& MODEL, const Cuts& cuts, const IloInt begin_stage) {
	int w;


	switch (types) {
		case  VariantType::EC_P_M_C:
		{
			MODEL.add(cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)]);
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[(numStage - begin_stage - 1)*nCut*numScen + (stage - begin_stage)*numScen + w]);
			}
			break;
		}
		case  VariantType::EC_D_S_C:
		{
			MODEL.add(cuts.optCut[(numStage - 1)*nCut + stage]);
			break;
		}
		case  VariantType::EC_D_M_C:
		{
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + w]);
			}
			break;
		}
		case  VariantType::EC_D_S_S_S:
		{
			MODEL.add(cuts.optCut[(numStage - 1) * 2 * nCut + stage]);
			MODEL.add(cuts.optCut[(numStage - 1) * 2 * nCut + (numStage - 1) + stage]);
			break;
		}
		case  VariantType::EC_D_S_M_S:
		{
			MODEL.add(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen) * nCut + stage]);
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen) * nCut + (numStage - 1) + stage*numScen + w]);
			}
			break;
		}
		case  VariantType::EC_D_M_S_S:
		{
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[((numStage - 1)*numScen + (numStage - 1)) * nCut + stage*numScen + w]);
			}
			MODEL.add(cuts.optCut[((numStage - 1)*numScen + (numStage - 1)) * nCut + (numStage - 1)*numScen + stage]);
			break;
		}
		case  VariantType::EC_D_M_M_S:
		{
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[2 * (numStage - 1)*numScen * nCut + stage*numScen + w]);
			}
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[2 * (numStage - 1)*numScen * nCut + (numStage - 1)*numScen + stage*numScen + w]);
			}
			break;
		}
		case  VariantType::TV_P_M_C:
		{
			MODEL.add(cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)]);
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[(numStage - begin_stage - 1)*nCut*numScen + (stage - begin_stage)*numScen + w]);
			}
			break;
		}
		case  VariantType::TV_D_S_C:
		{
			MODEL.add(cuts.optCut[(numStage - 1)*nCut + stage]);
			break;
		}
		case  VariantType::TV_D_M_C:
		{
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + w]);
			}
			break;
		}
		case  VariantType::TV_D_S_S_S:
		{
			MODEL.add(cuts.optCut[(numStage - 1) * 2 * nCut + stage]);
			MODEL.add(cuts.optCut[(numStage - 1) * 2 * nCut + (numStage - 1) + stage]);
			break;
		}
		case  VariantType::TV_D_S_M_S:
		{
			MODEL.add(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen) * nCut + stage]);
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen) * nCut + (numStage - 1) + stage*numScen + w]);
			}
			break;
		}
		case  VariantType::TV_D_M_S_S:
		{
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[((numStage - 1)*numScen + (numStage - 1)) * nCut + stage*numScen + w]);
			}
			MODEL.add(cuts.optCut[((numStage - 1)*numScen + (numStage - 1)) * nCut + (numStage - 1)*numScen + stage]);
			break;
		}
		case  VariantType::TV_D_M_M_S:
		{
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[2 * (numStage - 1)*numScen * nCut + stage*numScen + w]);
			}
			for (w = 0; w < numScen; w++) {
				MODEL.add(cuts.optCut[2 * (numStage - 1)*numScen * nCut + (numStage - 1)*numScen + stage*numScen + w]);
			}
			break;
		}
	}
}//InitializeCuts

void getConstraintsDuals(const IloInt& stage, const IloInt& scen, const IloInt& nCut,
	IloCplex& CPX, const Constraints& constraints, Duals& duals) {

	int i, t;

	for (t = 0; t<numYear[stage]; t++) {
		for (i = 0; i<numCapacity; i++)
			duals.pi_capacity[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage] + (scen - numScenNode_stage_sum[stage - 1])*numYear[stage] + t] = CPX.getDual(constraints.CapBound[numYear_stage_sum[stage] * numCapacity + t*numCapacity + i]);

		for (i = 0; i<numUser; i++)
			duals.pi_demand[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage] + (scen - numScenNode_stage_sum[stage - 1])*numYear[stage] + t] = CPX.getDual(constraints.MeetDemand[numYear_stage_sum[stage] * numUser + t*numUser + i]);

		for (i = 0; i<numPotUser; i++)
			duals.pi_return[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage] + (scen - numScenNode_stage_sum[stage - 1])*numYear[stage] + t] = CPX.getDual(constraints.ReturnFlow[numYear_stage_sum[stage] * numPotUser + t*numPotUser + i]);

		for (i = 0; i<numRecharge; i++)
			duals.pi_safe[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage] + (scen - numScenNode_stage_sum[stage - 1])*numYear[stage] + t] = CPX.getDual(constraints.SafeYield[numYear_stage_sum[stage] * numRecharge + t*numRecharge + i]);
	}
	for (i = 0; i<numRecharge; i++) {
		duals.pi_storage[i*numScenNode + scen] = CPX.getDual(constraints.StorageBalance[numYear_stage_sum[stage] * numRecharge + i]);
		duals.pi_flowbound[i*numScenNode + scen] = CPX.getDual(constraints.RFoutflowBound[numYear_stage_sum[stage] * numRecharge + i]);
	}

}//end getConstraintsDual

void getoptCutsDuals(const IloInt& stage, const IloInt& scen, const IloInt& nCut,
	IloCplex& CPX, const Cuts& cuts, IloNumArray& pi_optCut,
	const IloInt begin_stage) {

	int k, w;

	int numMidNodes = numScenNode - 1;
	int numParents = numScenNode_stage_sum[numStage - 2];

	switch (types) {
		case  VariantType::EC_P_M_C:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k*numMidNodes + descendant[scen][w] -1] = CPX.getDual(cuts.optCut[(numStage - begin_stage - 1)*k*numScen + (stage - begin_stage)*numScen + w]);
				}
			}
			break;
		case  VariantType::EC_D_S_C:
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k*numParents + scen] = CPX.getDual(cuts.optCut[(numStage - 1)*k + stage]);
			}
			break;
		case  VariantType::EC_D_M_C:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k*numMidNodes + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[(numStage - 1)*k*numScen + stage*numScen + w]);
				}
			}
			break;
		case  VariantType::EC_D_S_S_S:
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k*2*numParents + scen] = CPX.getDual(cuts.optCut[2*(numStage - 1)*k + stage]);
			}
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k * 2 * numParents + numParents + scen] = CPX.getDual(cuts.optCut[2 * (numStage - 1)*k + (numStage-1)+ stage]);
			}
			break;
		case  VariantType::EC_D_S_M_S:
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k * (numParents+ numMidNodes) + scen] = CPX.getDual(cuts.optCut[ ( (numStage - 1)+(numStage-1)*numScen )*k + stage]);
			}
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k* (numParents + numMidNodes)+ numParents + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen)*k + (numStage-1) +stage*numScen + w]);
				}
			}
			break;
		case  VariantType::EC_D_M_S_S:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k* (numMidNodes + numParents) + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[((numStage - 1)*numScen + (numStage - 1) )*k + stage*numScen + w]);
				}
			}
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k *  (numMidNodes + numParents) + numMidNodes +scen] = CPX.getDual(cuts.optCut[((numStage - 1)*numScen + (numStage - 1))*k + (numStage-1)*numScen + stage]);
			}
			break;
		case  VariantType::EC_D_M_M_S:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[2*k* numMidNodes + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[2*(numStage - 1)*numScen*k + stage*numScen + w]);
				}
			}
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[2 * k* numMidNodes + numMidNodes+ descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[2 * (numStage - 1)*numScen*k + (numStage-1)*numScen+ stage*numScen + w]);
				}
			}
			break;
		case  VariantType::TV_P_M_C:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k*numMidNodes + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[(numStage - begin_stage - 1)*k*numScen + (stage - begin_stage)*numScen + w]);
				}
			}
			break;
		case  VariantType::TV_D_S_C:
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k*numParents + scen] = CPX.getDual(cuts.optCut[(numStage - 1)*k + stage]);
			}
			break;
		case  VariantType::TV_D_M_C:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k*numMidNodes + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[(numStage - 1)*k*numScen + stage*numScen + w]);
				}
			}
			break;
		case  VariantType::TV_D_S_S_S:
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k * 2 * numParents + scen] = CPX.getDual(cuts.optCut[2 * (numStage - 1)*k + stage]);
			}
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k * 2 * numParents + numParents + scen] = CPX.getDual(cuts.optCut[2 * (numStage - 1)*k + (numStage - 1) + stage]);
			}
			break;
		case  VariantType::TV_D_S_M_S:
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k * (numParents + numMidNodes) + scen] = CPX.getDual(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen)*k + stage]);
			}
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k* (numParents + numMidNodes) + numParents + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[((numStage - 1) + (numStage - 1)*numScen)*k + (numStage - 1) + stage*numScen + w]);
				}
			}
			break;
		case  VariantType::TV_D_M_S_S:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[k* (numMidNodes + numParents) + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[((numStage - 1)*numScen + (numStage - 1))*k + stage*numScen + w]);
				}
			}
			for (k = 0; k < nCut + 1; k++) {
				pi_optCut[k *  (numMidNodes + numParents) + numMidNodes + scen] = CPX.getDual(cuts.optCut[((numStage - 1)*numScen + (numStage - 1))*k + (numStage - 1)*numScen + stage]);
			}
			break;
		case  VariantType::TV_D_M_M_S:
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[2 * k* numMidNodes + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[2 * (numStage - 1)*numScen*k + stage*numScen + w]);
				}
			}
			for (k = 0; k<nCut + 1; k++) {
				for (w = 0; w<numScen; w++) {
					pi_optCut[2 * k* numMidNodes + numMidNodes + descendant[scen][w] - 1] = CPX.getDual(cuts.optCut[2 * (numStage - 1)*numScen*k + (numStage - 1)*numScen + stage*numScen + w]);
				}
			}
			break;
	}

}//getoptCutsDuals

void createMaster(ModelArray& mod, CplexArray& CPX, Constraints& constraints,
	mainVariables& mainvariables, appxVariables& appxvariables,
	const pair <double, double>& param)
{
	IloInt i, t, s, w, j;
	IloEnv env = mod.getEnv();

	ModelArray modelTemp(env, numStage);
	mod.clear();
	mod.add(modelTemp);
	modelTemp.end();

	CplexArray cplexTemp(env, numStage);
	CPX.clear();
	CPX.add(cplexTemp);
	cplexTemp.end();


	for (s = 0; s < numStage; s++) {
		mod[s] = IloModel(env);
		CPX[s] = IloCplex(env);
	}

	mainvariables.Q=IloNumVarArray (env, numLink*numYear_stage_sum[numStage], 0, IloInfinity);
	mainvariables.Y=IloNumVarArray (env, numRecharge*numYear_stage_sum[numStage], 0, IloInfinity);
	
	for (s = 0; s<numStage; s++) {
		for (t = 0; t<numYear[s]; t++) {
			for (i = 0; i<numLink; i++) {
				char varName[100];
				sprintf_s(varName, "Q.%d.%d.%d", (int)s, (int)t, (int)i);
				mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + i].setName(varName);
			}

			for (i = 0; i<numRecharge; i++) {
				char varName[100];
				sprintf_s(varName, "Y.%d.%d.%d", (int)s, (int)t, (int)i);
				mainvariables.Y[numYear_stage_sum[s] * numRecharge + t*numRecharge + i].setName(varName);
			}

		} //end year
	} //end stage

	


	appxvariables.Theta=IloNumVarArray (env, total_theta, 0, IloInfinity); //first theta s for mean or sup, then theta s for CVaR: if applicable
	setThetaName(appxvariables.Theta);
	

	appxvariables.Alpha=IloNumVarArray (env, numStage - 1, 0, IloInfinity);


	for (s = 0; s<numStage - 1; s++) {
		char varName[100];
		sprintf_s(varName, "Alpha.%d", (int)s);
		appxvariables.Alpha[s].setName(varName);
	}

	//25a
	ObjArray OBJ(env, numStage);

	//25b
	constraints.FlowBalance=IloRangeArray (env, numBalance*numYear_stage_sum[numStage]);

	//25c
	constraints.MeetDemand=IloRangeArray(env, numUser*numYear_stage_sum[numStage]);

	//25d
	constraints.ReturnFlow=IloRangeArray (env, numPotUser*numYear_stage_sum[numStage]);

	//25e,f
	constraints.StorageBalance=IloRangeArray (env, numRecharge*numYear_stage_sum[numStage]);

	//25g
	constraints.SafeYield=IloRangeArray (env, numRecharge*numYear_stage_sum[numStage]);

	//25h,i
	constraints.RFoutflowBound=IloRangeArray (env, numRecharge*numYear_stage_sum[numStage]);

	//25j,k,l
	constraints.CapBound=IloRangeArray (env, numCapacity*numYear_stage_sum[numStage]);


	for (s = 0; s < numStage; s++) {
		OBJ[s] = IloExpr(env);
		for (t = 0; t < numYear[s]; t++) {
			IloNum val = pow(1 + disRate, -double(numYear_stage_sum[s] + t));
			//25a
			for (i = 0; i < numLink; i++) {
				OBJ[s] +=  val* cost[i]*mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + i];
			}

			//25b
			for (i = 0; i < numBalance; i++) {
				IloExpr LHS(env), RHS(env);
				for (j = 0; j < numLink; j++) {
					if (endNode[j] == balanceID[i])
						LHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j] * (1 - loss[j]);      //inflow
				}
				for (j = 0; j < numLink; j++) {
					if (startNode[j] == balanceID[i])
						RHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j];                 //outflow
				}
				constraints.FlowBalance[numYear_stage_sum[s] * numBalance + t*numBalance + i] = LHS - RHS == 0;
				char conName[100];
				sprintf_s(conName, "FB.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.FlowBalance[numYear_stage_sum[s] * numBalance + t*numBalance + i]);
				constraints.FlowBalance[numYear_stage_sum[s] * numBalance + t*numBalance + i].setName(conName);
				LHS.end(); RHS.end();
			}

			//25c
			for (i = 0; i < numUser; i++) {
				IloExpr LHS(env);
				IloNum RHS; //right hand side is the demand; inflow = demand
				if (s == 0)
					RHS = demand[numYear_stage_sum[s] * numUser + t*numUser + i] * scenFac[numScenNode_stage_sum[s] - 1][t + 1];
				else
					RHS = demand[numYear_stage_sum[s] * numUser + t*numUser + i] * scenFac[numScenNode_stage_sum[s - 1]][t + 1];

				for (j = 0; j < numLink; j++) {
					if (endNode[j] == userID[i])
						LHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j] * (1 - loss[j]);
				}
				constraints.MeetDemand[numYear_stage_sum[s] * numUser + t*numUser + i] = LHS == RHS;
				char conName[100];
				sprintf_s(conName, "DE.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.MeetDemand[numYear_stage_sum[s] * numUser + t*numUser + i]);
				constraints.MeetDemand[numYear_stage_sum[s] * numUser + t*numUser + i].setName(conName);
				LHS.end();
			}
			//25d
			for (i = 0; i < numPotUser; i++) {
				IloExpr LHS(env);
				IloNum RHS; // 98% of the potable uses return to the facility

				if (s == 0)
					RHS = returnRate*demand[numYear_stage_sum[s] * numUser + t*numUser + i] * scenFac[numScenNode_stage_sum[s] - 1][t + 1];
				else
					RHS = returnRate*demand[numYear_stage_sum[s] * numUser + t*numUser + i] * scenFac[numScenNode_stage_sum[s - 1]][t + 1];

				for (j = 0; j < numLink; j++) {
					if (startNode[j] == potUserID[i])
						LHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j];
				}
				constraints.ReturnFlow[numYear_stage_sum[s] * numPotUser + t*numPotUser + i] = LHS == RHS;

				char conName[100];
				sprintf_s(conName, "RE.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.ReturnFlow[numYear_stage_sum[s] * numPotUser + t*numPotUser + i]);
				constraints.ReturnFlow[numYear_stage_sum[s] * numPotUser + t*numPotUser + i].setName(conName);
				LHS.end();
			}

			//25e,f
			for (i = 0; i < numRecharge; i++) {
				IloExpr LHS(env), RHS(env);
				for (j = 0; j < numLink; j++) {
					if (endNode[j] == rechargeID[i])
						LHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j] * (1 - loss[j]);      //inflow
				}
				for (j = 0; j < numLink; j++) {
					if (startNode[j] == rechargeID[i])
						RHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j];                 //outflow
				}

				LHS -= mainvariables.Y[numYear_stage_sum[s] * numRecharge + t*numRecharge + i];
				if (s == 0 && t == 0)
					//inflow -Storage[i][t] - outflow  =  -Storage0[i]
					constraints.StorageBalance[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = LHS - RHS == -Storage0[i];
				else if (s>0 && t == 0) 
					constraints.StorageBalance[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = LHS - RHS == 0;
				else{
					//-Storage[i][t] + Storage[i][t-1] + inflow -outflow=  0
					LHS += mainvariables.Y[numYear_stage_sum[s] * numRecharge + (t - 1)*numRecharge + i];
					constraints.StorageBalance[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = LHS - RHS == 0;
				}
				char conName[100];
				sprintf_s(conName, "SB.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.StorageBalance[numYear_stage_sum[s] * numRecharge + t*numRecharge + i]);
				constraints.StorageBalance[numYear_stage_sum[s] * numRecharge + t*numRecharge + i].setName(conName);
				LHS.end(); RHS.end();
			}

			//25g
			for (i = 0; i < numRecharge; i++) { // recharge facility storage safe yield
				constraints.SafeYield[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = mainvariables.Y[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] <= storageUB[i];
				char conName[100];
				sprintf_s(conName, "SY.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.SafeYield[numYear_stage_sum[s] * numRecharge + t*numRecharge + i]);
				constraints.SafeYield[numYear_stage_sum[s] * numRecharge + t*numRecharge + i].setName(conName);
			}

			//25h,i
			for (i = 0; i < numRecharge; i++) {
				IloExpr LHS(env);
				for (j = 0; j < numLink; j++) {
					if (startNode[j] == rechargeID[i])
						LHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j];
				}
				if (s == 0 && t == 0)
					constraints.RFoutflowBound[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = LHS <= Storage0[i];
				else if (s > 0 && t == 0)
					constraints.RFoutflowBound[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = LHS  <= 0;
				else
					constraints.RFoutflowBound[numYear_stage_sum[s] * numRecharge + t*numRecharge + i] = LHS - mainvariables.Y[numYear_stage_sum[s] * numRecharge + (t - 1)*numRecharge + i] <= 0;

				char conName[100];
				sprintf_s(conName, "RF.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.RFoutflowBound[numYear_stage_sum[s] * numRecharge + t*numRecharge + i]);
				constraints.RFoutflowBound[numYear_stage_sum[s] * numRecharge + t*numRecharge + i].setName(conName);
				LHS.end();
			}

			//25j,k,l
			for (i = 0; i < numCapacity; i++) {
				IloExpr LHS(env);
				IloNum RHS;
				for (j = 0; j < numLink; j++) {
					if (startNode[j] == capacityID[i])
						LHS += mainvariables.Q[numYear_stage_sum[s] * numLink + t*numLink + j];      //outflow
				}
				if (i == 0 && s == 0)
					RHS = CAPamt[numYear_stage_sum[s] + t] * scenFac[numScenNode_stage_sum[s] - 1][0];
				else if (i == 0 && s>0)
					RHS = CAPamt[numYear_stage_sum[s] + t] * scenFac[numScenNode_stage_sum[s - 1]][0];
				else
					RHS = capacity[i];

				constraints.CapBound[numYear_stage_sum[s] * numCapacity + t*numCapacity + i] = LHS <= RHS;
				char conName[100];
				sprintf_s(conName, "CB.%d.%d.%d", (int)s, (int)t, (int)i);
				mod[s].add(constraints.CapBound[numYear_stage_sum[s] * numCapacity + t*numCapacity + i]);
				constraints.CapBound[numYear_stage_sum[s] * numCapacity + t*numCapacity + i].setName(conName);
				LHS.end();
			}

		}//end year
		//add approximation variables to OBJ[s]
		if (s < numStage - 1) {
			addThetaAlphatoOBJ(s, OBJ[s], appxvariables, param);
		}
		mod[s].add(IloMinimize(env, OBJ[s]));
		OBJ[s].end();
		CPX[s].extract(mod[s]);
		
	}//end stage
}//end createMaster

void createSub(IloCplex& CplexSub, IloModel& modSub, DistGenFormulation& distgen, const pair <double, double> param)
{
	double prob = (double)1 / numScen;
	IloInt w;
	IloEnv env = CplexSub.getEnv();
	CplexSub.extract(modSub);

	distgen.P=IloNumVarArray (env, numScen, 0, 1);
	modSub.add(distgen.P);

	for (w = 0; w<numScen; w++) {
		char varName[100];
		sprintf_s(varName, "P.%d", (int)w);
		distgen.P[w].setName(varName);
	}

	distgen.Z=IloNumVarArray (env, numScen, 0, 1);
	modSub.add(distgen.Z);

	for (w = 0; w<numScen; w++) {
		char varName[100];
		sprintf_s(varName, "Z.%d", (int)w);
		distgen.Z[w].setName(varName);
	}

	IloExpr conDistance(env);
	IloRange RngDistance;
	IloExpr conProb(env);
	IloRange RngProb;

	distgen.objSub = IloObjective(env);
	distgen.objSub.setSense(IloObjective::Maximize);
	modSub.add(distgen.objSub);

	distgen.conDistancePos = IloRangeArray(env, numScen);
	distgen.conDistanceNeg = IloRangeArray(env, numScen);

	if (variant.ambiguity_type == TV) {
		for (w = 0; w < numScen; w++) {
			distgen.conDistancePos[w]= distgen.P[w] - distgen.Z[w] <= prob;
			distgen.conDistanceNeg[w]= -distgen.P[w] - distgen.Z[w] <= -prob;
			conDistance += distgen.Z[w];
		}

		double rho = param.first;
		RngDistance = conDistance <= 2*rho;
		modSub.add(RngDistance);
	}
	else {
		double lambda = param.first;
		double alpha = param.second;
		for (w = 0; w < numScen; w++) {
			distgen.conDistancePos[w]=distgen.P[w] <= prob * (1 + lambda*alpha / (1 - alpha));
			distgen.conDistanceNeg[w]=-distgen.P[w] <= -prob * (1 - lambda);
		}
	}

	modSub.add(distgen.conDistancePos);
	modSub.add(distgen.conDistanceNeg);

	
	for (w = 0; w < numScen; w++) {
		conProb += distgen.P[w];
	}
	RngProb = conProb == 1;
	modSub.add(RngProb);
	
}//endcreateSub

void ForwardPassX(const IloInt& nCut, CplexArray& CPX, Constraints& constraints, Cuts& cuts,
	const mainVariables& mainvariables, const appxVariables& appxvariables,
	Duals& duals, const CutCoeffs& cutscoeffs, incSols& incsols, CostsIx& costsix,
	IloNum& LB,
	const pair <double, double>& param,
	const IloBool IsCondEff, const IloInt begin_stage, const IloInt w_ancestor, const IloNumArray Best_y_grand_ancestor)
{
	IloInt s, w, t, i;
	IloEnv env = constraints.FlowBalance.getEnv();

	IloNumArray yTemp(env, numRecharge*numScenNode);
	if (nCut == 0 || variant.problem_type == DUAL)
		incsols.y.add(yTemp);
	yTemp.end();

	if (variant.problem_type == DUAL)
		InitializeIndex(costsix);

	for (w = 0; w<numScenNode; w++)
		costsix.subObj_hat[w] = 0;

	if (begin_stage == 0)
		w = 0;
	else
		w = numScenNode_stage_sum[begin_stage - 1];

	if (IsCondEff && begin_stage>0) {
		IloInt w_grand_ancestor = FindAncestor(begin_stage, w_ancestor);
		for (i = 0; i<numRecharge; i++) {
			//it was w before, I changed it w_grand_ancestor! check results!
			incsols.y[i*numScenNode + w_grand_ancestor] = Best_y_grand_ancestor[i*numScenNode + w_grand_ancestor];
		}
	}

	for (s = begin_stage; s < numStage; s++) {
		while (w < numScenNode_stage_sum[s]) {
			if (w == w_ancestor || ancestor[w_ancestor][w] > 0) {// do the forward pass on the subtree of w_ancestor	
				UpdateRhs(s, w, nCut, constraints, incsols.y);

				if (s < numStage - 1 && nCut>0)
					UpdateOldCuts(s, w, nCut, cuts, mainvariables.Y, appxvariables.Alpha[s], appxvariables.Theta, cutscoeffs, incsols, costsix,
						param, begin_stage);

				CPX[s].solve();
				/*char resName[100];
				sprintf_s(resName, "%d.lp", w);
				CPX[s].exportModel(resName);*/

				if (s == begin_stage) {
					if (!CPX[begin_stage].solve()) {
						env.error() << "Master Problem Infeasible" << endl;
						throw(-1);
					}
				}


				if (s == begin_stage)
					LB = CPX[s].getObjValue();


				if (variant.problem_type == PRIMAL)
					for (i = 0; i<numRecharge; i++)
						incsols.y[i*numScenNode + w] = CPX[s].getValue(mainvariables.Y[numYear_stage_sum[s] * numRecharge + (numYear[s] - 1)*numRecharge + i]);
				else
					for (i = 0; i<numRecharge; i++)
						incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + w] = CPX[s].getValue(mainvariables.Y[numYear_stage_sum[s] * numRecharge + (numYear[s] - 1)*numRecharge + i]);


				//store dual variables at last stage
				costsix.subObj_hat[w] = CPX[s].getObjValue();
				if (s == numStage - 1)
					getConstraintsDuals(s, w, nCut, CPX[s], constraints, duals);
				else
					subtractThetaAlphafromsubObj(s, CPX[s], appxvariables, costsix.subObj_hat[w], param);

			}
			w++;
		}
	}
}

void ForwardPassP(const IloInt& nCut, IloCplex& CplexSub, IloModel& modSub, DistGenFormulation& distgen,
	IloNumArray& p, IloNumArray & subObj_hat,
	IloBool& contForwardPass,
	const pair <double, double> param,
	const IloBool IsCondEff, const IloBool IsPathEff, const IloInt begin_stage, const IloInt w_ancestor, const IloInt child)
{
	IloInt s, w, j;
	IloEnv env = distgen.objSub.getEnv();

	IloNumArray pTemp(env, numScenNode);
	p.add(pTemp);
	pTemp.end();


	IloRange RngZero;
	IloBool addP = IloFalse;
	IloBool IsAssess = IsCondEff || IsPathEff;
	if (IsAssess) {
		RngZero = distgen.P[child] == 0;
	}
	if (IsCondEff) {
		p[nCut*numScenNode + w_ancestor] = 1;
	}
	else {
		p[nCut*numScenNode] = 1;
	}

	for (s = numStage - 1; s > begin_stage && contForwardPass; s--) {
		if (s == 1)
			w = 0;
		else
			w = numScenNode_stage_sum[s - 2];
		
		while (w < numScenNode_stage_sum[s - 1] && contForwardPass) {
			if (IsAssess && w == w_ancestor) {
				addP = IloTrue;
				modSub.add(RngZero);
			}
			if (IsPathEff || w == w_ancestor || ancestor[w_ancestor][w] > 0) {
				for (j = 0; j < numScen; j++) 
					distgen.objSub.setLinearCoef(distgen.P[j], subObj_hat[descendant[w][j]]);
				
				if (variant.ambiguity_type == TV) {
					for (j = 0; j < numScen; j++) {
						distgen.conDistancePos[j].setUb(prob[descendant[w][j]]);
						distgen.conDistanceNeg[j].setUb(-prob[descendant[w][j]]);
					}
				}
				else {
					for (j = 0; j < numScen; j++) {
						double lambda = param.first;
						double alpha = param.second;
						distgen.conDistancePos[j].setUb(prob[descendant[w][j]] * (1 +lambda *alpha / (1 - alpha)));
						distgen.conDistanceNeg[j].setUb(-prob[descendant[w][j]] * (1 - lambda));
					}
				}
				CplexSub.extract(modSub);
				CplexSub.solve();
				if (CplexSub.solve()) {
					for (j = 0; j < numScen; j++) {
						p[nCut*numScenNode + descendant[w][j]] = CplexSub.getValue(distgen.P[j]);
						subObj_hat[w] += p[nCut*numScenNode + descendant[w][j]] * subObj_hat[descendant[w][j]];
					}
				}
				else {
					contForwardPass = IloFalse;
				}
				if (w == w_ancestor && addP) {
					modSub.remove(RngZero);
				}
			}
			w++;
		}
	}
}

void ForwardPassP2(const IloInt& nCut, 
	IloNumArray& p, IloNumArray & subObj_hat,
	const pair <double, double> param,
	const IloBool IsCondEff, const IloBool IsPathEff, const IloInt begin_stage, const IloInt w_ancestor, const IloInt child)
{
	IloInt s, w, j;
	IloEnv env = subObj_hat.getEnv();

	int numParents = numScenNode_stage_sum[numStage - 2];

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;


	IloNumArray pTemp(env, numScenNode);
	p.add(pTemp);
	pTemp.end();

	for (s = numStage - 1; s > begin_stage; s--) {
		if (s == 1)
			w = 0;
		else
			w = numScenNode_stage_sum[s - 2];

		while (w < numScenNode_stage_sum[s - 1] ) {
			NodeInfo * wOutput = new NodeInfo[numScen];
			for (j = 0; j < numScen; j++) {
				wOutput[j].No = j;
				wOutput[j].Cost = subObj_hat[descendant[w][j]];
			}
			sort(wOutput, wOutput + numScen, CompareScenariosCostIndex());
			int varix = wOutput[CVaRno].No;
			double var_value = wOutput[CVaRno].Cost;
			double sup_value = wOutput[numScen - 1].Cost;
			double lambda = abs(var_value - sup_value);
		
			
			if (variant.ambiguity_type == TV) {
				if (lambda > toler) {
					for (j = 0; j < CVaRno - 1; j++)
						p[nCut*numScenNode + descendant[w][wOutput[j].No]] = 0;

					p[nCut*numScenNode + descendant[w][wOutput[CVaRno].No]] = (CVaRno + 1)*prob[descendant[w][wOutput[CVaRno].No]] - rho;

					for (j = CVaRno + 1; j < numScen - 1; j++)
						p[nCut*numScenNode + descendant[w][wOutput[j].No]] = prob[descendant[w][wOutput[j].No]];

					p[nCut*numScenNode + descendant[w][wOutput[numScen - 1].No]] = prob[descendant[w][wOutput[numScen - 1].No]] + rho;
				}
				else {
					for (j = 0; j < CVaRno - 1; j++)
						p[nCut*numScenNode + descendant[w][wOutput[j].No]] = 0;

					p[nCut*numScenNode + descendant[w][wOutput[CVaRno].No]] = 1;
				}
				
				for (j = 0; j < numScen; j++)
					subObj_hat[w] += p[nCut*numScenNode + descendant[w][j]] * subObj_hat[descendant[w][j]];
			}
			else if (variant.ambiguity_type == EC)
			{
				for (j = 0; j < numScen; j++) {
					p[nCut*numScenNode + descendant[w][j]] = (1 - lambda)*prob[descendant[w][j]];
					//subObj_hat[w] += p[nCut*numScenNode + descendant[w][j]] * subObj_hat[descendant[w][j]];//This takes care of Expectation
				}
				for (j = CVaRno+1; j < numScen; j++) {
					p[nCut*numScenNode + descendant[w][wOutput[j].No]] += lambda*prob[descendant[w][wOutput[j].No]]/(1-alpha);
					//subObj_hat[w] += p[nCut*numScenNode + descendant[w][wOutput[j].No]] * subObj_hat[descendant[w][wOutput[j].No]];
				}
				double remained_prob = 1;
				for (j = 0; j < numScen; j++) {
					remained_prob -= p[nCut*numScenNode + descendant[w][j]];
				}
				p[nCut*numScenNode + descendant[w][varix]] += remained_prob;

				for (j = 0; j < numScen; j++) 
					subObj_hat[w] += p[nCut*numScenNode + descendant[w][j]] * subObj_hat[descendant[w][j]];
			
			}
			delete[] wOutput;			
			w++;
		}
	}
}

void CalculateFuncVal(const IloInt& nCut, CostsIx& costsix,
	const pair<double, double> param) {

	int w, s, j; 

	int numParents= numScenNode_stage_sum[numStage- 2];

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;

	for (s = numStage - 1; s > 0; s--) {
		if (s == 1)
			w = 0;
		else
			w = numScenNode_stage_sum[s - 2];

		while (w < numScenNode_stage_sum[s - 1]) {
			NodeInfo * wOutput = new NodeInfo[numScen];
			for (j = 0; j < numScen; j++) {
				wOutput[j].No = j;
				wOutput[j].Cost = costsix.subObj_hat[descendant[w][j]];
			}
			sort(wOutput, wOutput+ numScen, CompareScenariosCostIndex());
			int varix= wOutput[CVaRno].No;
			int supix = wOutput[numScen - 1].No;
			if (s-1 == numStage - 2) {
				costsix.varIndex[nCut*numParents + w]=varix;
				costsix.supIndex[nCut*numParents + w]=supix;
			}

			double * sorted_cost = new double[numScen];
			for (j = 0; j < numScen; j++) 
				sorted_cost[j] = wOutput[j].Cost;
			if (variant.ambiguity_type == TV) {
				double CVaR=CalculateCVaR(sorted_cost, rho);
				double Sup = wOutput[numScen - 1].Cost;
				costsix.subObj_hat[w] += rho*Sup + (1- rho)*CVaR;
			}
			else if (variant.ambiguity_type == EC) 
			{
				double Exp= CalculateExp(sorted_cost);
				double CVaR = CalculateCVaR(sorted_cost, alpha);
				costsix.subObj_hat[w] += (1-lambda)*Exp + lambda*CVaR;
			}	
			
			delete[] sorted_cost;
			delete[] wOutput;
			w++;
		}
		
	}
}//end CalculateFuncVal

void BackwardPass(const IloInt& nCut, ModelArray& MODEL, CplexArray& CPX,
	Constraints& constraints, Cuts& cuts, const IloNumVarArray& Y, const appxVariables& appxvariables,
	Duals& duals, CutCoeffs& cutscoeffs, const incSols& incsols, CostsIx& costsix,
	const pair <double, double>& param,
	const IloInt begin_stage, const IloInt w_ancestor)
{
	IloInt w, i, j, t, s, w2, k;
	IloEnv env = MODEL[begin_stage].getEnv();
	int numParents = numScenNode_stage_sum[numStage - 2];

	InitializeCutsnadCoeffs(cuts, duals.pi_optCut, cutscoeffs, begin_stage);
	
	
	for (s = numStage - 1; s>begin_stage; s--) {
		AddEmptyCuts(s-1, nCut, MODEL[s - 1], cuts, begin_stage);
		
		if (s == 1)
			w = 0;
		else
			w = numScenNode_stage_sum[s - 2];
	
		while (w < numScenNode_stage_sum[s - 1]) {
			if (w == w_ancestor || ancestor[w_ancestor][w] > 0) { //do the backward pass only on the subtree of w_ancestor
				
				UpdateRhs(s-1, w, nCut, constraints, incsols.y);

				UpdateOldCuts(s - 1, w, nCut, cuts, Y, appxvariables.Alpha[s - 1], appxvariables.Theta, cutscoeffs, incsols, costsix,
					param, begin_stage);

				GenerateCutCoeffs(s - 1, w, nCut, duals, cutscoeffs, incsols, costsix, param, begin_stage);
				

				int varix = 0;
				int supix = 0;
				if (variant.problem_type == DUAL) {
					varix = costsix.varIndex[nCut*numParents + w];
					supix = costsix.supIndex[nCut*numParents + w];
				}

				
				GenerateCuts(s - 1, w, nCut, cuts, Y, appxvariables.Alpha[s - 1], appxvariables.Theta, cutscoeffs, incsols, varix, supix,
					param, begin_stage);
				

				if (s>begin_stage + 1) {
					CPX[s - 1].solve();
									

					getConstraintsDuals(s - 1, w, nCut, CPX[s - 1], constraints, duals);		

					getoptCutsDuals(s - 1, w, nCut, CPX[s - 1], cuts, duals.pi_optCut, begin_stage);
					
					if (variant.problem_type == DUAL)
						costsix.subObj_hat[w] = CPX[s - 1].getObjValue();

				}
			}
			w++;
		}
		if (s>1 && variant.problem_type == DUAL)
			FindIndex(s - 1, nCut, costsix);
	}//end backward  
}


void GenerateCutCoeffs(const IloInt& stage, const IloInt& scenario, const IloInt& nCut,
	const Duals& duals, CutCoeffs& cutcoeffs, const incSols& incsols, const CostsIx& costsix,
	const pair<double, double>& param,
	const IloInt begin_stage){

	IloInt i, j, t, k, w2;

	int numParents = numScenNode_stage_sum[numStage - 2];
	int numMidNodes = numScenNode - 1;

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;

	//take care of G
	for (i = 0; i<numRecharge; i++)
		for (j=0; j<numScen; j++)
			cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] = -duals.pi_storage[i*numScenNode + descendant[scenario][j]] + duals.pi_flowbound[i*numScenNode + descendant[scenario][j]];
	

	//tage care of \pi*b part of g
	for (j = 0; j < numScen; j++) {
		for (t = 0; t < numYear[stage + 1]; t++) {
			for (i = 0; i < numUser; i++) {
				cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_demand[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage + 1] + (descendant[scenario][j] - numScenNode_stage_sum[stage])*numYear[stage + 1] + t] * scenFac[descendant[scenario][j]][t + 1] * demand[numYear_stage_sum[stage + 1] * numUser + t*numUser + i];
			}
			for (i = 0; i < numPotUser; i++) {
				cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_return[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage + 1] + (descendant[scenario][j] - numScenNode_stage_sum[stage])*numYear[stage + 1] + t] * scenFac[descendant[scenario][j]][t + 1] * returnRate* demand[numYear_stage_sum[stage + 1] * numUser + t*numUser + i];
			}
			for (i = 1; i < numCapacity; i++) {
				cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_capacity[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage + 1] + (descendant[scenario][j] - numScenNode_stage_sum[stage])*numYear[stage + 1] + t] * capacity[i];
			}
			cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_capacity[total_pi_stage_sum[stage + 1] + (descendant[scenario][j] - numScenNode_stage_sum[stage])*numYear[stage + 1] + t] * CAPamt[numYear_stage_sum[stage + 1] + t] * scenFac[descendant[scenario][j]][0];
			for (i = 0; i < numRecharge; i++) {
				cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_safe[i*total_pi_stage_sum[numStage] + total_pi_stage_sum[stage + 1] + (descendant[scenario][j] - numScenNode_stage_sum[stage])*numYear[stage + 1] + t] * storageUB[i];
			}
		}
	}
	

	switch (types) {

		case  VariantType::EC_P_M_C:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1];
						}
					}
				}
			}
			/*IloEnv env2;
			IloExpr aggcut(env2);
			char varName[100];			
			IloNumVarArray XX(env2, numRecharge, 0, IloInfinity);
			for (i = 0; i < numRecharge; i++) {
				sprintf_s(varName, "Y.%d", (int)i);
				XX[i].setName(varName);
			}
			for (j = 0; j < numScen; j++) {
				aggcut += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] * incsols.p[nCut*numScenNode + descendant[scenario][j]];
				for (i = 0; i < numRecharge; i++) {
					aggcut += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.p[nCut*numScenNode + descendant[scenario][j]] * XX[i];
				}
			}
			Cutsforms << "nCut=" << nCut << ", w=" << scenario << ":" << aggcut << endl;*/

			break;
		}
		case  VariantType::EC_D_S_C:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1]  * (1 - lambda);
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.psi[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * lambda / (1 - alpha);
						}
						int varix = costsix.varIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][varix] - 1] * lambda;
						for (i = 0; i < numRecharge; i++) {
							for (w2 = 0; w2 < numScen; w2++) {
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.psi_x[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][w2] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] * lambda / (1 - alpha);
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.psi_eta[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * cutcoeffs.G[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][varix] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] * lambda / (1 - alpha);
							}
						}
					}
				}
			}
			break;
		}
		case  VariantType::EC_D_M_C:
		{

			break;
		}
		case  VariantType::EC_D_S_S_S:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k * 2 * numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]];
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * cutcoeffs.psi[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] / (1 - alpha);
						}
						int varix = costsix.varIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][varix] - 1];
						for (i = 0; i < numRecharge; i++) {
							for (w2 = 0; w2 < numScen; w2++) {
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_x[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][w2] - 1) + i] * duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * prob[descendant[descendant[scenario][j]][w2]] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - alpha);
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_eta[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * cutcoeffs.G[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][varix] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - alpha);
							}
						}
					}
				}
			}
			break;
		}
		case  VariantType::EC_D_S_M_S:
		{

			break;
		}
		case  VariantType::EC_D_M_S_S:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*(numMidNodes + numParents) + descendant[descendant[scenario][j]][w2] - 1] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1];
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * cutcoeffs.psi[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] / (1 - alpha);
						}
						int varix = costsix.varIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][varix] - 1];
						for (i = 0; i < numRecharge; i++) {
							for (w2 = 0; w2 < numScen; w2++) {
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_x[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][w2] - 1) + i] * duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * prob[descendant[descendant[scenario][j]][w2]] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - alpha);
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_eta[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * cutcoeffs.G[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][varix] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - alpha);
							}
						}
					}
				}
			}
			break;
		}
		case  VariantType::EC_D_M_M_S:
		{

			break;
		}
		case  VariantType::TV_P_M_C:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1];
						}
					}
				}
			}
			break;
		}
		case  VariantType::TV_D_S_C:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						int supix = costsix.supIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][supix] - 1] * rho;
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.psi[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] * (1 - rho) / (1 - rho);
						}
						int varix = costsix.varIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][varix] - 1] * (1 - rho);
						for (i = 0; i < numRecharge; i++) {
							for (w2 = 0; w2 < numScen; w2++) {
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_x[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][w2] - 1) + i] * duals.pi_optCut[k*numParents + descendant[scenario][j]] * prob[descendant[descendant[scenario][j]][w2]] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] * (1 - rho) / (1 - rho);
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_eta[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*numParents + descendant[scenario][j]] * cutcoeffs.G[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][varix] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] * (1 - rho) / (1 - rho);
							}
						}
					}
				}
			}
			break;
		}
		case  VariantType::TV_D_M_C:
		{

			break;
		}
		case  VariantType::TV_D_S_S_S:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						int supix = costsix.supIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k * 2 * numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][supix] - 1];
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * cutcoeffs.psi[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] / (1 - rho);
						}
						int varix = costsix.varIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][varix] - 1];
						for (i = 0; i < numRecharge; i++) {
							for (w2 = 0; w2 < numScen; w2++) {
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_x[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][w2] - 1) + i] * duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * prob[descendant[descendant[scenario][j]][w2]] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - rho);
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_eta[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k * 2 * numParents + numParents + descendant[scenario][j]] * cutcoeffs.G[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][varix] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - rho);
							}
						}
					}
				}
			}
			break;
		}
		case  VariantType::TV_D_S_M_S:
		{

			break;
		}
		case  VariantType::TV_D_M_S_S:
		{
			if (stage + 1 < numStage - 1) {
				for (j = 0; j < numScen; j++) {
					for (k = 0; k < nCut + 1; k++) {
						for (w2 = 0; w2 < numScen; w2++) {
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*(numMidNodes + numParents) + descendant[descendant[scenario][j]][w2] - 1] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1];
							cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * cutcoeffs.psi[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] / (1 - rho);
						}
						int varix = costsix.varIndex[k*numParents + descendant[scenario][j]];
						cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] += duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * cutcoeffs.g[k*numMidNodes + descendant[descendant[scenario][j]][varix] - 1];
						for (i = 0; i < numRecharge; i++) {
							for (w2 = 0; w2 < numScen; w2++) {
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_x[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][w2] - 1) + i] * duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * prob[descendant[descendant[scenario][j]][w2]] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - rho);
								cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] -= cutcoeffs.psi_eta[k*numMidNodes + descendant[descendant[scenario][j]][w2] - 1] * prob[descendant[descendant[scenario][j]][w2]] * duals.pi_optCut[k*(numMidNodes + numParents) + numMidNodes + descendant[scenario][j]] * cutcoeffs.G[k*numRecharge*numMidNodes + numRecharge*(descendant[descendant[scenario][j]][varix] - 1) + i] * incsols.y[k*numRecharge*numScenNode + i*numScenNode + descendant[scenario][j]] / (1 - rho);
							}
						}
					}
				}
			}
			break;
		}
		case  VariantType::TV_D_M_M_S:
		{
			break;
		}

	}

	/*if (variant.problem_type == DUAL) {
		int varix = costsix.varIndex[nCut*numParents + scenario];
		IloEnv env = cutcoeffs.g.getEnv();
		IloNumArray cost(env, numScen);
		for (j = 0; j < numScen; j++) {
			cost[j] = 0;
			cost[j] += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1];
			for (i = 0; i < numRecharge; i++)
				cost[j] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode +i*numScenNode + scenario];
		}

		for (j = 0; j < numScen; j++) {
			if (cost[j] - cost[varix] > 0) {
				cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] = -1;
				cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] = cost[j] - cost[varix];
				for (i = 0; i < numRecharge; i++) {
					cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
				}
			}
		}
	}*/


	if (variant.problem_type == DUAL) {
		int varix = costsix.varIndex[nCut*numParents + scenario];
		for (j = 0; j < numScen; j++) {
			if (costsix.subObj_hat[descendant[scenario][j]] - costsix.subObj_hat[descendant[scenario][varix]] > 0) {
				cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] = -1;
				cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] = costsix.subObj_hat[descendant[scenario][j]] - costsix.subObj_hat[descendant[scenario][varix]];
				for (i = 0; i < numRecharge; i++) {
					cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
				}
			}
		}
	}
}//GenerateCutCoeffs


void GenerateCuts(const IloInt& stage, const IloInt& scenario, const IloInt& nCut,
	Cuts& cuts, 
	const IloNumVarArray& Y, const IloNumVar& Alpha, const IloNumVarArray& Theta,
	const CutCoeffs& cutcoeffs, const incSols& incsols, const IloInt& varIndex, const IloInt& supIndex,
	const pair<double, double>& param,
	const IloInt begin_stage)
{
	IloEnv env = cuts.optCut.getEnv();
	IloInt i, j, t, k, w2;

	int numParents = numScenNode_stage_sum[numStage - 2];
	int numMidNodes = numScenNode - 1;

	double lambda = param.first;
	double alpha = param.second;
	double rho = param.first;
	
	switch (types) {

		case  VariantType::EC_P_M_C:
		{


			IloNumArray RHS(env, numScen);
			IloNumArray LHS(env, numScen*numRecharge);
			for (j = 0; j < numScen; j++) {
				RHS[j] = cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1];
				cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + j].setBounds(RHS[j], IloInfinity);
				cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + j].setLinearCoef(Theta[stage*numScen + j], 1);
				for (i = 0; i < numRecharge; i++) {
					LHS[numRecharge*j + i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
					cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + j].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -LHS[numRecharge*j + i]);
				}
				
			}
			cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)].setBounds(0, IloInfinity);
			cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)].setLinearCoef(Alpha, 1);
			for (j = 0; j < numScen; j++) {
				cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)].setLinearCoef(Theta[stage*numScen + j], -incsols.p[nCut*numScenNode + descendant[scenario][j]]);
			}

			RHS.end();
			LHS.end();

			break;
		}
		case  VariantType::EC_D_S_C:
		{

			IloNum RHS = 0;
			IloNumArray LHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				RHS += prob[descendant[scenario][j]] * cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1] * (1 - lambda);
				RHS += prob[descendant[scenario][j]] * cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] * lambda / (1 - alpha);
				for (i = 0; i < numRecharge; i++) {
					RHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] * lambda / (1 - alpha);
					RHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] * lambda / (1 - alpha);
				}
			}
			RHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][varIndex] - 1] * lambda;
			cuts.optCut[(numStage - 1)*nCut + stage].setBounds(RHS, IloInfinity);
			cuts.optCut[(numStage - 1)*nCut + stage].setLinearCoef(Theta[stage], 1);

			for (i = 0; i < numRecharge; i++) {
				LHS[i] = 0;
				for (j = 0; j < numScen; j++) {
					LHS[i] += prob[descendant[scenario][j]] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * (1 - lambda);
					LHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * lambda / (1 - alpha);
					LHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * lambda / (1 - alpha);
				}
				LHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * lambda;

				cuts.optCut[(numStage - 1)*nCut + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -LHS[i]);
			}
			LHS.end();

			

			break;
		}
		case  VariantType::EC_D_M_C:
		{

			break;
		}
		case  VariantType::EC_D_S_S_S:
		{
			IloNum ExpRHS = 0;
			IloNumArray ExpLHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				ExpRHS += prob[descendant[scenario][j]] * cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1];
			}
			cuts.optCut[2 * (numStage - 1)*nCut + stage].setBounds(ExpRHS, IloInfinity);
			cuts.optCut[2 * (numStage - 1)*nCut + stage].setLinearCoef(Theta[stage], 1);
			for (i = 0; i < numRecharge; i++) {
				ExpLHS[i] = 0;
				for (j = 0; j < numScen; j++) {
					ExpLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
				}
				cuts.optCut[2 * (numStage - 1)*nCut + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -ExpLHS[i]);
			}

			IloNum CVaRRHS = 0;
			IloNumArray CVaRLHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				CVaRRHS += prob[descendant[scenario][j]] * cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] / (1 - alpha);
				for (i = 0; i < numRecharge; i++) {
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - alpha);
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - alpha);
				}
			}
			CVaRRHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][varIndex] - 1];
			cuts.optCut[2 * (numStage - 1)*nCut + numStage - 1 + stage].setBounds(CVaRRHS, IloInfinity);
			cuts.optCut[2 * (numStage - 1)*nCut + numStage - 1 + stage].setLinearCoef(Theta[numStage - 1 + stage], 1);

			for (i = 0; i < numRecharge; i++) {
				CVaRLHS[i] = 0;
				for (j = 0; j < numScen; j++) {
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] / (1 - alpha);
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] / (1 - alpha);
				}
				CVaRLHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i];

				cuts.optCut[2 * (numStage - 1)*nCut + numStage - 1 + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -CVaRLHS[i]);
			}
			ExpLHS.end();
			CVaRLHS.end();
			break;
		}
		case  VariantType::EC_D_S_M_S:
		{

			break;
		}
		case  VariantType::EC_D_M_S_S:
		{
			IloNumArray ExpRHS(env, numScen);
			IloNumArray ExpLHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				ExpRHS[j] = cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1];
				cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + stage*numScen + j].setBounds(ExpRHS[j], IloInfinity);
				cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + stage*numScen + j].setLinearCoef(Theta[stage*numScen + j], 1);
			}

			for (j = 0; j < numScen; j++) {
				for (i = 0; i < numRecharge; i++) {
					ExpLHS[i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
					cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + stage*numScen + j].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -ExpLHS[i]);
				}
			}

			IloNum CVaRRHS = 0;
			IloNumArray CVaRLHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				CVaRRHS += prob[descendant[scenario][j]] * cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] / (1 - alpha);
				for (i = 0; i < numRecharge; i++) {
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - alpha);
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - alpha);
				}
			}
			CVaRRHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][varIndex] - 1];
			cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + (numStage - 1)*numScen + stage].setBounds(CVaRRHS, IloInfinity);
			cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + (numStage - 1)*numScen + stage].setLinearCoef(Theta[numScen*(numStage - 1) + stage], 1);

			for (i = 0; i < numRecharge; i++) {
				CVaRLHS[i] = 0;
				for (j = 0; j < numScen; j++) {
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] / (1 - alpha);
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] / (1 - alpha);
				}
				CVaRLHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i];

				cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + (numStage - 1)*numScen + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -CVaRLHS[i]);
			}
			ExpRHS.end();
			ExpLHS.end();
			CVaRLHS.end();
			break;
		}
		case  VariantType::EC_D_M_M_S:
		{

			break;
		}
		case  VariantType::TV_P_M_C:
		{
			//Cutsforms << "nCut=" << nCut << ", w=" << scenario << ": ";
			//vector<int>::iterator match; //if found the item in the effective scenarios

			IloNumArray RHS(env, numScen);
			IloNumArray LHS(env, numScen*numRecharge);
			
			for (j = 0; j < numScen; j++) {
				//match = std::find(eff->begin(), eff->end(), j);
				//if (match != eff->end()) {
				//if (incsols.p[nCut*numScenNode + descendant[scenario][j]]>0){ //these are effective
					//Cutsforms << j << ", ";
					RHS[j] = cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1];
					cuts.optCut[(numStage - begin_stage - 1)*nCut*numScen + (stage - begin_stage)*numScen + j].setBounds(RHS[j], IloInfinity);
					cuts.optCut[(numStage - begin_stage - 1)*nCut*numScen + (stage - begin_stage)*numScen + j].setLinearCoef(Theta[stage*numScen + j], 1);
					for (i = 0; i < numRecharge; i++) {
						LHS[numRecharge*j + i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
						cuts.optCut[(numStage - begin_stage - 1)*nCut*numScen + (stage - begin_stage)*numScen + j].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -LHS[numRecharge*j + i]);
					}					
				//}
				//else {
				//	cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + j].setBounds(0, IloInfinity);
				//	cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + j].setLinearCoef(Theta[stage*numScen + j], 0);
				//	for (i = 0; i < numRecharge; i++) {
				//		cuts.optCut[(numStage - 1)*nCut*numScen + stage*numScen + j].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], 0);
				//	}
				//}
			}
			cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)].setBounds(0, IloInfinity);
			cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)].setLinearCoef(Alpha, 1);
			for (j = 0; j < numScen; j++) {
				cuts.maxCut[(numStage - begin_stage - 1)*nCut + (stage - begin_stage)].setLinearCoef(Theta[stage*numScen + j], -incsols.p[nCut*numScenNode + descendant[scenario][j]]);
			}
			
			RHS.end();
			LHS.end();
			
			//Cutsforms << endl;
			break;
		}
		case  VariantType::TV_D_S_C:
		{
			IloNum RHS = 0;
			IloNumArray LHS(env, numRecharge);
			RHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][supIndex] - 1] * rho;
			for (j = 0; j < numScen; j++) {
				RHS += prob[descendant[scenario][j]] * cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] * (1 - rho) / (1 - rho);
				for (i = 0; i < numRecharge; i++) {
					RHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] * (1 - rho) / (1 - rho);
					RHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] * (1 - rho) / (1 - rho);
				}
			}
			RHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][varIndex] - 1] * (1 - rho);
			cuts.optCut[(numStage - 1)*nCut + stage].setBounds(RHS, IloInfinity);
			cuts.optCut[(numStage - 1)*nCut + stage].setLinearCoef(Theta[stage], 1);

			for (i = 0; i < numRecharge; i++) {
				LHS[i] = 0;
				LHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][supIndex] - 1) + i] * rho;
				for (j = 0; j < numScen; j++) {
					LHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * (1 - rho) / (1 - rho);
					LHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * (1 - rho) / (1 - rho);
				}
				LHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * (1 - rho);

				cuts.optCut[(numStage - 1)*nCut + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -LHS[i]);
			}
			LHS.end();
			break;
		}
		case  VariantType::TV_D_M_C:
		{

			break;
		}
		case  VariantType::TV_D_S_S_S:
		{
			IloNum SupRHS = 0;
			IloNumArray SupLHS(env, numRecharge);
			SupRHS = cutcoeffs.g[nCut*numMidNodes + descendant[scenario][supIndex] - 1];
			cuts.optCut[2 * (numStage - 1)*nCut + stage].setBounds(SupRHS, IloInfinity);
			cuts.optCut[2 * (numStage - 1)*nCut + stage].setLinearCoef(Theta[stage], 1);
			for (i = 0; i < numRecharge; i++) {
				SupLHS[i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][supIndex] - 1) + i];
				cuts.optCut[2 * (numStage - 1)*nCut + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -SupLHS[i]);
			}

			IloNum CVaRRHS = 0;
			IloNumArray CVaRLHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				CVaRRHS += prob[descendant[scenario][j]] * cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] / (1 - rho);
				for (i = 0; i < numRecharge; i++) {
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - rho);
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - rho);
				}
			}
			CVaRRHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][varIndex] - 1];
			cuts.optCut[2 * (numStage - 1)*nCut + numStage - 1 + stage].setBounds(CVaRRHS, IloInfinity);
			cuts.optCut[2 * (numStage - 1)*nCut + numStage - 1 + stage].setLinearCoef(Theta[numStage - 1 + stage], 1);

			for (i = 0; i < numRecharge; i++) {
				CVaRLHS[i] = 0;
				for (j = 0; j < numScen; j++) {
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] / (1 - rho);
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] / (1 - rho);
				}
				CVaRLHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i];

				cuts.optCut[2 * (numStage - 1)*nCut + numStage - 1 + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -CVaRLHS[i]);
			}
			SupLHS.end();
			CVaRLHS.end();
			break;
		}
		case  VariantType::TV_D_S_M_S:
		{

			break;
		}
		case  VariantType::TV_D_M_S_S:
		{
			IloNumArray SupRHS(env, numScen);
			IloNumArray SupLHS(env, numScen*numRecharge);
			for (j = 0; j < numScen; j++) {
				SupRHS[j] = cutcoeffs.g[nCut*numMidNodes + descendant[scenario][j] - 1];
				cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + stage*numScen + j].setBounds(SupRHS[j], IloInfinity);
				cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + stage*numScen + j].setLinearCoef(Theta[stage], 1);
			}

			for (j = 0; j < numScen; j++) {
				for (i = 0; i < numRecharge; i++) {
					SupLHS[j*numRecharge + i] = cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i];
					cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + stage*numScen + j].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -SupLHS[j*numRecharge + i]);
				}
			}

			IloNum CVaRRHS = 0;
			IloNumArray CVaRLHS(env, numRecharge);
			for (j = 0; j < numScen; j++) {
				CVaRRHS += prob[descendant[scenario][j]] * cutcoeffs.psi[nCut*numMidNodes + descendant[scenario][j] - 1] / (1 - rho);
				for (i = 0; i < numRecharge; i++) {
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - rho);
					CVaRRHS -= prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] * incsols.y[nCut*numRecharge*numScenNode + i*numScenNode + scenario] / (1 - rho);
				}
			}
			CVaRRHS += cutcoeffs.g[nCut*numMidNodes + descendant[scenario][varIndex] - 1];
			cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + (numStage - 1)*numScen + stage].setBounds(CVaRRHS, IloInfinity);
			cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + (numStage - 1)*numScen + stage].setLinearCoef(Theta[numStage - 1 + stage], 1);

			for (i = 0; i < numRecharge; i++) {
				CVaRLHS[i] = 0;
				for (j = 0; j < numScen; j++) {
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_x[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][j] - 1) + i] / (1 - rho);
					CVaRLHS[i] += prob[descendant[scenario][j]] * cutcoeffs.psi_eta[nCut*numMidNodes + descendant[scenario][j] - 1] * cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i] / (1 - rho);
				}
				CVaRLHS[i] += cutcoeffs.G[nCut*numRecharge*numMidNodes + numRecharge*(descendant[scenario][varIndex] - 1) + i];

				cuts.optCut[(numScen* (numStage - 1) + numStage - 1)*nCut + (numStage - 1)*numScen + stage].setLinearCoef(Y[numYear_stage_sum[stage] * numRecharge + (numYear[stage] - 1)*numRecharge + i], -CVaRLHS[i]);
			}
			SupRHS.end();
			SupLHS.end();
			CVaRLHS.end();
			break;
		}
		case  VariantType::TV_D_M_M_S:
		{

			break;
		}
	}
	
}


void UpdateOldCuts(const IloInt& stage, const IloInt& scenario, const IloInt& nCut,
	Cuts& cuts,
	const IloNumVarArray& Y, const IloNumVar& Alpha, const IloNumVarArray& Theta,
	const CutCoeffs& cutcoeffs, const incSols& incsols, const CostsIx& costsix,
	const pair<double, double> param,
	const IloInt begin_stage)
{
	int k;
	int numParents = numScenNode_stage_sum[numStage - 2];

	
	for (k = 0; k < nCut; k++) {
		int varix = 0;
		int supix = 0;
		if (variant.problem_type == DUAL) {
			varix = costsix.varIndex[k*numParents + scenario];
			supix = costsix.supIndex[k*numParents + scenario];
		}
		GenerateCuts(stage, scenario, k, cuts, Y, Alpha, Theta, cutcoeffs, incsols,
				varix, supix, param, begin_stage);
		
	}
	

}//end UpdateOldCuts


void UpdateRhs(const IloInt& stage, const IloInt& scenario, const IloInt& nCut, Constraints& constraints, const IloNumArray& y)
{
	IloInt i, t; 
	for (t = 0; t < numYear[stage]; t++) {
		constraints.CapBound[numYear_stage_sum[stage]*numCapacity + t*numCapacity].setBounds(-IloInfinity, CAPamt[numYear_stage_sum[stage] + t] * scenFac[scenario][0]);
		for (i = 0; i < numUser; i++)
			constraints.MeetDemand[numYear_stage_sum[stage] * numUser + t*numUser + i].setBounds(scenFac[scenario][t + 1] * demand[numYear_stage_sum[stage] * numUser + t*numUser + i], scenFac[scenario][t + 1] * demand[numYear_stage_sum[stage] * numUser + t*numUser + i]);
		for (i = 0; i < numPotUser; i++)
			constraints.ReturnFlow[numYear_stage_sum[stage] * numPotUser + t*numPotUser + i].setBounds(scenFac[scenario][t + 1] * returnRate*demand[numYear_stage_sum[stage] * numUser + t*numUser + i], scenFac[scenario][t + 1] * returnRate*demand[numYear_stage_sum[stage] * numUser + t*numUser + i]);
	}
	IloInt w2 = FindAncestor(stage, scenario);
	if (stage > 0) {
		if (variant.problem_type == PRIMAL) {
			for (i = 0; i<numRecharge; i++) {
				constraints.StorageBalance[numYear_stage_sum[stage] * numRecharge + i].setBounds(-y[i*numScenNode + w2], -y[i*numScenNode + w2]);
				constraints.RFoutflowBound[numYear_stage_sum[stage] * numRecharge + i].setBounds(-IloInfinity, y[i*numScenNode + w2]);
			}
		}
		else {
			for (i = 0; i < numRecharge; i++) {
				constraints.StorageBalance[numYear_stage_sum[stage] * numRecharge + i].setBounds(-y[nCut*numRecharge*numScenNode + i*numScenNode + w2], -y[nCut*numRecharge*numScenNode + i*numScenNode + w2]);
				constraints.RFoutflowBound[numYear_stage_sum[stage] * numRecharge + i].setBounds(-IloInfinity, y[nCut*numRecharge*numScenNode + i*numScenNode + w2]);
			}
		}
	}
	
}

void SolveOriginal(ostream& Result, ModelArray MODEL, CplexArray CPX, Formulation formulation,
	IloCplex CplexSub, IloModel modSub, DistGenFormulation distgen,
	IloNum& objVal,
	const pair<double, double> param,
	Scenario * Output)
{
	IloEnv env = MODEL[0].getEnv();
	IloInt i, w;
	Result << "Iter  \t    Z_hat    \t     LB      \t  UB  \t Gap%" << endl;

	IloNum rel_Gap = IloInfinity;
	IloNum LB = -IloInfinity;
	IloNum UB = IloInfinity;
	IloNum z_hat;
	IloInt nCut = 0;

	const int numParents = numScenNode_stage_sum[numStage-2];
	
	const clock_t begin_time = clock();
	while (rel_Gap > toler  && float(clock() - begin_time) / CLOCKS_PER_SEC <= 7200) {
		ForwardPassX(nCut, CPX, formulation.constraints, formulation.cuts, formulation.mainvariables, formulation.appxvariables,
			formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.costsix, LB, param);
		
		IloBool contForwardPass = IloTrue;
		if (variant.problem_type == DUAL) {
			CalculateFuncVal(nCut, formulation.costsix, param);
		}
		else {
			ForwardPassP(nCut, CplexSub, modSub, distgen, formulation.incsols.p, formulation.costsix.subObj_hat, contForwardPass, param);
			//ForwardPassP2(nCut, formulation.incsols.p, formulation.costsix.subObj_hat, param);
		}

	

		z_hat = formulation.costsix.subObj_hat[0];
		if (z_hat < UB) {
			UB = z_hat;
			for (w = 0; w<numScenNode; w++) {
				if (variant.problem_type == PRIMAL) {
					for (i = 0; i<numRecharge; i++)
						formulation.optsols.Best_y[i*numScenNode + w] = formulation.incsols.y[i*numScenNode + w];
					formulation.optsols.Worst_p[w] = formulation.incsols.p[nCut*numScenNode + w];
					Output[w].worst_Prob = formulation.optsols.Worst_p[w];
				}
				else {
					for (i = 0; i<numRecharge; i++)
						formulation.optsols.Best_y[i*numScenNode + w] = formulation.incsols.y[nCut*numRecharge*numScenNode+i*numScenNode + w];
				}
				Output[w].Cost = formulation.costsix.subObj_hat[w];
			}
		}
		
		rel_Gap = abs(UB - LB) / abs(LB);
		Result << nCut + 1 << '\t' << z_hat << '\t' << LB << '\t' << UB << '\t' << rel_Gap*100 << endl;
		//Result << "time = " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
		if (rel_Gap > toler  && float(clock() - begin_time) / CLOCKS_PER_SEC <= 7200) {
			BackwardPass(nCut, MODEL, CPX, formulation.constraints, formulation.cuts, formulation.mainvariables.Y,
					formulation.appxvariables, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.costsix, param);
		}
		
	
		nCut++;
	}
	objVal = LB;
	Result << "time = " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
	Result << "Y[0] = [";
	for (i = 0; i<numRecharge-1; i++)
		Result << formulation.optsols.Best_y[i]<< ", ";
	Result << formulation.optsols.Best_y[numRecharge - 1] << "]" << endl;
	Result << "cost = " << formulation.costsix.subObj_hat << endl;
	if (variant.problem_type==PRIMAL)
		Result << "worst prob = " << formulation.optsols.Worst_p << endl;
	

	cout << "*****DONE with Original Problem*****" << endl;
}

