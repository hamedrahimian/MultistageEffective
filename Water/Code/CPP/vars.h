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

#ifndef VARS_H
#define VARS_H

#include "types.h"

extern int numStage, numScen, numScenNode, numLink;
extern int numUser, numPotUser, numBalance, numCapacity, numRecharge;
extern int total_pi;
extern int total_theta;
extern std::vector<int> numYear_stage_sum;
extern std::vector<int> numScenNode_stage_sum, total_pi_stage_sum;
extern const double  disRate;
extern std::vector<int> userID, potUserID, rechargeID, balanceID, capacityID;
extern IloIntArray startNode, endNode, numYear;
extern IloNumArray2 scenFac;
extern IloNumArray   loss, cost, capacity, Storage0, storageUB, demand, CAPamt;
extern IloNum   returnRate;
extern IloNumArray prob;
extern IloIntArray stageMap;
extern IloIntArray2 ancestor, descendant;
extern int CVaRno;

extern Variant variant;
extern VariantType types;


#endif