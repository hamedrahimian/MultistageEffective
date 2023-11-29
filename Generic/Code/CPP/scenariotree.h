/*
*     Nested L-shaped method for DRSO-- Built on top of SUTIL
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

#ifndef SCENARIOTREE_H
#define SCENARIOTREE_H

#include <ilcplex/ilocplex.h>
#include <SPProblem.h>
#include <vector>
#include <array>
#include "sparseMatrix.h"

class TreeStucture : public SPProblem {
	//friend class Decomposition;
	//friend class Effectiveness;
	public:

		/** Gets the total number of scenarios nodes by period #period*/
		int numTotalNodesInstance(const int period) const;

		/** Gets the overall node index of node #ix in period #period*/
		int getNodeOverallIndex(const int period, const int ix) const;

		/** Gets the index of node #ix in period #period, among its sibilings*/
		int getChildIndex(const int period, const int ix) const;

		/** Gets the overall node index of #j-th child of node #ix in period #period*/
		int getChildOverallIndex(const int period, const int ix, const int j) const;

		/** Gets the node index of #j-th child of node #ix in period #period, in period #period+1*/
		int getChildRelativeIndex(const int period, const int ix, const int j) const;

		/** Returns the depth (the number of edges) from scenrio node #ix in period #period to #ix with a positive number */
		// if the output is zero, #ix is not in the subtree formed fromed by #w_ancestor. 
		int getDepth(const int period, const int ix, const int cperiod, const int cix) const;

		/** Print all scenarios to output file*/
		void PrintScenariosToFile(char *problem_name);

};




#endif

