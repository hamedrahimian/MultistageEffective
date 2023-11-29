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

#ifndef TREE_H
#define TREE_H

#include <ilcplex/ilocplex.h>

void TreeStructure(IloIntArray2& ancestor, IloIntArray2& descendant, IloIntArray& stageMap);

IloInt FindAncestor(const IloInt& stage, const IloInt& scenario, const IloInt t = 1);

/* Given #scenario and its #stage, finds #scenario corresponds to which child number of #scenario's ancestor */
IloInt FindChild(const IloInt& stage, const IloInt& scenario);

#endif

