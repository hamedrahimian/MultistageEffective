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

#include "vars.h"

int numScen = 8;
int total_pi = 0;
int total_theta = 0;
const double  disRate = 0.02;

int numStage, numScenNode, numLink;
int numUser, numPotUser, numBalance, numCapacity, numRecharge;
std::vector<int> numYear_stage_sum;
std::vector<int> numScenNode_stage_sum, total_pi_stage_sum;
std::vector<int> userID, potUserID, rechargeID, balanceID, capacityID;
IloIntArray startNode, endNode, numYear;
IloNumArray2 scenFac;
IloNumArray   loss, cost, capacity, Storage0, storageUB, demand, CAPamt;
IloNum   returnRate;
IloNumArray prob;
IloIntArray stageMap;
IloIntArray2 ancestor, descendant;
int CVaRno;

Variant variant;
VariantType types;



