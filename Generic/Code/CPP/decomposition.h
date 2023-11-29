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


#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <ilcplex/ilocplex.h>
#include "scenariotree.h"

ILOSTLBEGIN

enum ProblemType {
	RISK_NEUTRAL = 0,
	DRSO
};

enum CutType {
	SINGLE = 0,
	MULTI
};

/** Type to create an array of IloNumVarArray objects, array for periods */
typedef IloArray<IloNumVarArray> VarArray2;
/** Type to create an array of IloRangeArray objects, array for periods */
typedef IloArray<IloRangeArray>  EquationArray;
/** Type to create an array of IloCplex objects, array for periods */
typedef IloArray<IloCplex>       CplexArray;
/** Type to create an array of IloModel objects, array for periods */
typedef IloArray<IloModel>       ModelArray;
/** Type to create an array of IloObjective objects, array for periods */
typedef IloArray<IloObjective>   ObjectiveArray;


/** Types to store optimality cut intercepts */
typedef double gNodeType; 
typedef vector<gNodeType>  cut_g; //vector for children
typedef vector<cut_g> gNode; //vector for cut

/** Types to store optimality cut slopes */
typedef vector<double> GNodeType; //vector for dimension of Tx portion
typedef  vector<GNodeType>  cut_G; //vector for children
typedef vector<cut_G> GNode; //vector for cut


/** Type to store max cut slopes */
typedef double PNodeType;
typedef  vector<PNodeType>  cut_P; //vector for children
typedef vector<cut_P> PNode; //vector for cut



class Decomposition : public TreeStucture {
	public: 
		
		/** Builds base subproblems */
		void BuildBaseSubModel(IloModel* model, IloNumVarArray* X, IloNumVar* Alpha, IloObjective* obj, IloRangeArray* rng, 
			SparseCSC * Wmat, char ** rngSense, double ** b, double ** c, double ** xl, double ** xu, int * objSense, const int period, const ProblemType problem_type=DRSO);

		/** Retrieve base subproblems after updates is made according to a scenario*/
		void RetrieveBaseSubModel(IloNumVarArray* X, IloObjective* obj, IloRangeArray* rng, 
			SparseCSC * Wmat, const char * rngSense, double ** b, double ** c, double ** bdl, double ** bdu);

		/** Print base lp models to output file*/
		void PrintBaseModelToFile(char *problem_name, const int period);

		/** Builds the inner (row generation) model for node #ix in period #period*/
		void BuildRowModel(IloModel * model, IloNumVarArray * P, const IloNumArray* subObj_hat, const int* objSense, 
			const int period, const int ix, const double rho=0);

		/** Performs the forward pass on X, at iteration #nCut */
		//IsCondEff determines whether this is the forward pass for checking condoitional effectiveness of #cIx
		//If IsCondEff, #pPeriod determines the stage of parent scenario #pIx. Also, #Best_x is the fixed value of #x in the current iteration for the ancestor of #pIx
		void ForwardPassX(IloNumArray2* x, IloNumArray2* rng_pi, 
			vector<GNode>* Gs, vector<gNode>* gs, vector<PNode>* ps,
			int* objSense, IloNum* LB, IloNumArray* subObj, const int nCut, const ProblemType problem_type= DRSO, const CutType cut_type= MULTI,
			const bool IsCondEff = false, const int pPeriod = 0, const int pIx = 0, const int cIx = 0, const IloNumArray* Best_x = 0);

		/** Performs the forward pass on P given current #subObj_hat (given x), at iteration #nCut */
		//IsCondEff and IsPathEff determines whether this is the forward pass for checking condoitional effectiveness or path effectiveness
		//If IsCondEff, #pPeriod determines the stage of parent scenario #pIx, where scenario #child (with relative index among its siblings) with index #cIx is tested. Otherwise (if IsPathEFf), #pPeriod=0 is root node
		void ForwardPassP(vector<cut_P>* ps, IloNumArray * pArray, IloNumArray* subObj_hat, const int* objSense,
			const int nCut, bool* contForwardPass, const ProblemType problem_type = DRSO, const CutType cut_type = MULTI, const double rho = 0,
			const bool IsCondEff = false, const bool IsPathEff = false,
			const int pPeriod = 0, const int pIx = 0, const int cIx = 0, const int child = 0);

		/** Calculate Upper Bounds (min) or Lower Bound (max) for Risk-Neutal model */
		void CalculateFuncVal(IloNumArray* subObj_hat);

		/** Performs the backward pass given current value of x, etc, at iteration #nCut */
		//If IsCondEff, #pPeriod determines the stage of parent scenario #pIx, where scenario #child (with relative index among its siblings) with index #cIx is tested. 	
		void BackwardPass(IloNumArray2* x, IloNumArray2* rng_pi, 
			vector<GNode>* Gs, vector<gNode>* gs, vector<PNode>* ps,
			const int nCut, const ProblemType problem_type=DRSO, const CutType cut_type=MULTI,
			const bool IsCondEff = false, const int pPeriod = 0, const int pIx = 0, const int cIx = 0, const int child = 0);

		/** Gets RHS for node #ix in period #period*/
		void Decomposition::GetScenarioRHSTechnologyMatrix(double ** b, double ** blo, double ** bup, SparseCSC* Tmat,
			const int period, const int ix);

		/** Updates the scenario information for node #ix in period #period*/
		void UpdateScenario(IloNumVarArray* X, IloObjective* obj, IloRangeArray* rng, 
			SparseCSC* Wmat, const char * rngSense, 
			double ** b, const IloNumArray* xhat, const int period, const int ix);

		/** Updates the Theta variables for node #ix in period #period*/
		void UpdateThetaVariables(IloNumVarArray* Theta, IloObjective* obj, const int period, const int ix,
			const ProblemType problem_type, const CutType cut_type);

		/** Generates the optimality cut coefficients for node #ix in last period , which later will be used to generate a cut for its ancestor in period #T-1 */
		void GenerateOptCutCoefficients(GNodeType * G, gNodeType *  g, const SparseCSC * Tmat, const double * b, const double * bl, const double * bu,
			const IloNumArray* rng_pi, gNode* gChild, const IloNumArray* pi_optCut, const int period, int ix, const int nCut,
			const ProblemType problem_type=DRSO, const CutType cut_type=MULTI);
				
		
		/** Generate optimility and max cuts for node #ix in period #period, given corresponding #G, #g, and #p from its children */
		void GenerateCuts(const cut_G * G, const cut_g *  g, IloRangeArray * optCut, const cut_P* p, IloRange * maxCut,
			const int * ncols, const int *objSense,
			const IloNumVarArray*  X, const IloNumVar*  Alpha, const IloNumVarArray*  Theta,
			const int nCut, const int ix, const int period, const ProblemType problem_type=DRSO, const CutType cut_type=MULTI,
			const bool IsCondEff = false, const int child = 0);


		/** Updates the previous cuts for node #ix in period #period, given corresponding #Gs, #gs, and #ps from previous iterations*/
		void UpdateOldCuts(GNode* Gs, gNode* gs, IloRangeArray* optCut, PNode* ps, IloRangeArray* maxCut,
			const int * ncols, const int *objSense,
			const IloNumVarArray* X, const IloNumVar* Alpha, const IloNumVarArray* Theta,
			const int period, const int ix, const int nCut, 
			const ProblemType problem_type=DRSO, const CutType cut_type=MULTI, const bool IsCondEff = false, const int child = 0);
		
		/** Solves original problem*/
		void SolveOriginal(char *problem_name, const double rho, double* objVal, IloNumArray2* xOpt, IloNumArray* SubObjOpt, IloNumArray* Worst_pOpt, const ProblemType problem_type=DRSO, const CutType cut_type=MULTI);
			

};


#endif
