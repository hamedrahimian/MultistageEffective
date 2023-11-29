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

#include "tree.h"
#include "vars.h"

void TreeStructure(IloIntArray2& ancestor, IloIntArray2& descendant, IloIntArray& stageMap)
{
	IloInt t, s, w, w2, i, j;
	IloEnv env;

	numScenNode_stage_sum.push_back(1);
	for (s = 1; s < numStage; s++) {
		numScenNode_stage_sum.push_back(numScenNode_stage_sum[s - 1] + pow(numScen, s));
	}
	//determine stage of each node in the tree
	stageMap = IloIntArray(env, numScenNode);
	w = 0;
	for (s = 0; s < numStage; s++) {
		while (w < numScenNode_stage_sum[s]) {
			stageMap[w] = s;
			w++;
		}
	}

	ancestor = IloIntArray2(env, numScenNode_stage_sum[numStage - 2]);
	for (w = 0; w < numScenNode_stage_sum[numStage - 2]; w++) {
		ancestor[w] = IloIntArray(env, numScenNode);
		for (w2 = 0; w2 < numScenNode; w2++) {
			ancestor[w][w2] = 0;
		}
	}

	//ancestor
	for (t = 1; t < numStage; t++) {
		for (i = 0; i < numScenNode_stage_sum[numStage - t - 1]; i++) {
			for (j = numScenNode_stage_sum[t] + pow(numScen, t)*(i - 1); j < numScenNode_stage_sum[t] + pow(numScen, t)*i; j++) {
				ancestor[i][j] = t; //i is ancestor j with "t" step // t=0 means i is not ancestor j
			}
		}
	}


	//denotes decendant of each node in the scenario tree by an array, whose elements are the child node number
	descendant = IloIntArray2(env, numScenNode_stage_sum[numStage - 2]);
	for (w = 0; w < numScenNode_stage_sum[numStage - 2]; w++) {
		descendant[w] = IloIntArray(env, numScen);
		for (j = 0; j < numScen; j++) {
			descendant[w][j] = 0;
		}
	}

	for (i = 0; i < numScenNode_stage_sum[numStage - 2]; i++) {
		for (j = 0; j < numScen; j++) {
			descendant[i][j] = 1 + numScen + numScen*(i - 1) + j;
		}
	}


}

IloInt FindAncestor(const IloInt& stage, const IloInt& scenario, const IloInt t)
{
	IloInt w2;

	if (stage>0 && stage >= t) {
		IloBool IsAncestor = IloFalse;
		if (stage == 1) {
			w2 = 0;
		}
		else
		{
			w2 = numScenNode_stage_sum[stage - t - 1];
		}
		while (w2<numScenNode_stage_sum[stage - t] && !IsAncestor) {
			if (ancestor[w2][scenario] == t) {
				IsAncestor = IloTrue;
				return w2;
			}
			w2++;
		}
	}
	else
	{
		return 0;
	}
}

IloInt FindChild(const IloInt& stage, const IloInt& scenario)
{
	IloInt j = 0;
	IloBool IsChild = IloFalse;
	IloInt w_ancestor = FindAncestor(stage, scenario);
	while (j<numScen && !IsChild) {
		if (descendant[w_ancestor][j] == scenario) {
			IsChild = IloTrue;
			return j;
		}
		j++;
	}
}