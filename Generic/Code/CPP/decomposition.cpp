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

#include "decomposition.h"
#include "scenariotree.h"
#include <algorithm>



using namespace std;

#define TOLER		0.000001
#define TIME_LIMIT  15000
#define MYINFINITY  1.0e+10



void Decomposition::BuildBaseSubModel(IloModel* model, IloNumVarArray* X, IloNumVar* Alpha, IloObjective* obj, IloRangeArray* rng, 
	SparseCSC * Wmat, char ** rngSense, double ** b, double ** c, double ** xl, double ** xu, int * objSense, const int period, const ProblemType problem_type)
{

	int col, row; //dimension of X variables and constraints counters
	double val; //value of a coefficient
	int element;	//counter for loop over entries of sparse matrix Wmat


	char varName[100];
	char rngName[100];

	IloEnv env = model->getEnv();
	const int numStage = getNumPeriods();

	//To hold LP info
	int os = 0;
	SparseCSC recourseMat;
	double *cost = NULL;
	double *rhs = NULL;
	char *sense = NULL;
	double *bdl = NULL;
	double *bdu = NULL;
	int * coltype = NULL;


	int LPstat = 0;
	//get subproblem technology matrix (W), coefficients, etc
	//it does not have scenario node information
	LPstat = getBaseLP(period, &recourseMat.NumCols, &recourseMat.NumRows,
		&recourseMat.NumEntries, &os, &cost, &rhs,
		&sense, &recourseMat.col_ptr, &recourseMat.row_ind,
		&recourseMat.val, &bdl, &bdu,
		&coltype);

	//rngSense, b, and Wmat  are output of this function, used in ChangeSubScenario function
	*rngSense = sense;
	*b = rhs;
	*Wmat = recourseMat;
	*objSense = os;
	*xu = bdu;
	*xl = bdl;
	*c = cost;


	if (LPstat != 0) {
		cerr << "getBaseLP: error" << endl;
		return;
	}


	//create X variables 
	X->clear();
	for (col = 0; col<Wmat->NumCols; col++) {
		sprintf_s(varName, "X.%d", (int)col);
		if (coltype[col]== SUTIL_VAR_CONTINUOUS)
			X->add(IloNumVar(env, bdl[col], bdu[col]));
		else
			X->add(IloIntVar(env, bdl[col], bdu[col]));
		(*X)[col].setName(varName);
	}
	model->add(*X);


	//create Alpha variable
	if (period< numStage-1 && problem_type == DRSO) {
		*Alpha = IloNumVar(env, -MYINFINITY, MYINFINITY);
		sprintf_s(varName, "Alpha.%d", (int)period);
		Alpha->setName(varName);
		model->add(*Alpha);
		obj->setLinearCoef(*Alpha, 1);
	}
	

	//Set objective function

	if (os == MAXIMIZE) 
		obj->setSense(IloObjective::Maximize);
	else 
		obj->setSense(IloObjective::Minimize);
	

	for (col = 0; col<Wmat->NumCols; col++) {
		obj->setLinearCoef((*X)[col], cost[col]);
	}

	model->add(*obj);

	
	//set constraints
	for (row = 0; row<Wmat->NumRows; row++) {
		switch ((*rngSense)[row]) {
		case 'E':
		case 'e':
			rng->add(IloRange(env, (*b)[row], (*b)[row]));
			break;
		case 'L':
		case 'l':
			rng->add(IloRange(env, -IloInfinity, (*b)[row]));
			break;
		case 'G':
		case 'g':
			rng->add(IloRange(env, (*b)[row], IloInfinity));
			break;
		}
	}


	for (element = 0; element<Wmat->NumEntries; element++) {
		col = FindEntryColInCSCMatrix(&Wmat->NumCols, &element, &Wmat->col_ptr);
		row = Wmat->row_ind[element];
		val = Wmat->val[element];
		(*rng)[row].setLinearCoef((*X)[col], val);
	}

	model->add(*rng);

	
	
}//end BuildSubModel


void Decomposition::RetrieveBaseSubModel(IloNumVarArray* X, IloObjective* obj, IloRangeArray* rng, 
	SparseCSC * Wmat, const char * rngSense, double ** b, double ** c, double ** bdl, double ** bdu)
{
	IloEnv env = X->getEnv();
	int row, col, element;
	double val;

	for (col = 0; col<Wmat->NumCols; col++) {
		obj->setLinearCoef((*X)[col], (*c)[col]);
		(*X)[col].setBounds((*bdl)[col], (*bdu)[col]);
	}
	
	for (row = 0; row<Wmat->NumRows; row++) {
		switch (rngSense[row]) {
		case 'E':
		case 'e':
			(*rng)[row].setBounds((*b)[row], (*b)[row]);
			break;
		case 'L':
		case 'l':
			(*rng)[row].setBounds(-IloInfinity, (*b)[row]);
			break;
		case 'G':
		case 'g':
			(*rng)[row].setBounds((*b)[row], IloInfinity);
			break;
		}
	}

	for (element = 0; element<Wmat->NumEntries; element++) {
		col = FindEntryColInCSCMatrix(&Wmat->NumCols, &element, &Wmat->col_ptr);
		row = Wmat->row_ind[element];
		val = Wmat->val[element];
		(*rng)[row].setLinearCoef((*X)[col], val);
	}

}//end RetrieveBaseSubModel

void Decomposition::PrintBaseModelToFile(char *problem_name, const int period)
{
	IloEnv env;

	SparseCSC  Wmat;
	char * rngSense;
	double * b, *c, *bdl, *bdu;
	int os = 0;

	IloModel MODEL(env);
	IloCplex CPX(env);
	IloNumVarArray X(env);
	IloNumVar Alpha(env);
	IloObjective OBJ(env);
	IloRangeArray RNG(env);
	

	BuildBaseSubModel(&MODEL, &X, &Alpha, &OBJ, &RNG, &Wmat, &rngSense, &b, &c, &bdl, &bdu, &os, period);

	CPX.extract(MODEL);

	char resName[100];
	sprintf_s(resName, "EX_%s_%d.lp", problem_name, period);
	CPX.exportModel(resName);

	X.clear();
	X.end();
	Alpha.end();
	OBJ.end();
	RNG.clear();
	RNG.end();
	MODEL.end();
	CPX.end();
	env.end();

}//end PrintBaseModelToFile

 void Decomposition::BuildRowModel(IloModel * model, IloNumVarArray * P, const IloNumArray* subObj_hat, const int* objSense, 
	 const int period, const int ix, const double rho) {
	int child; // child scenario nodes counter

	char varName[100];
	IloEnv env = model->getEnv();

	int ChildPeriod = 0;
	int numChild = 0;
	int * ChildIndices = NULL;
	int CHILDstat = 0;
	//Gets children scenarios of scenario node #ix in period #period
	CHILDstat=getChildScenarioIndices(period, ix, &ChildPeriod, &numChild, &ChildIndices);

	if (CHILDstat != 0) {
		cerr << "getChildScenarioIndices: error" << endl;
		return;
	}

	//(conditional) nominal probabilities 
	double *  q= new double [numChild];
	//Gets condtional probabilities of children scenarios
	for (child = 0; child<numChild; child++) {
		q[child] = getCondProbability(ChildPeriod, ChildIndices[child]);
	}

	IloNumVarArray Z(env);
	P->clear();
	//Creates P and Z variables to lineraize model
	for (child = 0; child<numChild; child++) {
		P->add(IloNumVar(env, 0, 1));
		sprintf_s(varName, "P_%d", (int)child);
		(*P)[child].setName(varName);

		Z.add(IloNumVar(env, 0, 1));
		sprintf_s(varName, "Z_%d", (int)child);
		Z[child].setName(varName);
	}
	model->add(*P);
	model->add(Z);

	IloObjective obj(env);
	if (*objSense==MAXIMIZE)
		obj.setSense(IloObjective::Minimize);	
	else
		obj.setSense(IloObjective::Maximize);

	for (child = 0; child<numChild; child++) 
		obj.setLinearCoef((*P)[child], (*subObj_hat)[child]);
	

	model->add(obj);

	IloRangeArray conDistancePos(env);
	IloRangeArray conDistanceNeg(env);

	for (child = 0; child<numChild; child++) {
		conDistancePos.add((*P)[child] - Z[child] <= q[child]);
		conDistanceNeg.add(-(*P)[child] - Z[child] <= -q[child]);
	}

	model->add(conDistancePos);
	model->add(conDistanceNeg);

	IloExpr conDistance(env);
	IloRange RngDistance;
	IloExpr conProb(env);
	IloRange RngProb;

	for (child = 0; child<numChild; child++) {
		conDistance += Z[child];
		conProb += (*P)[child];
	}

	RngDistance = conDistance <= 2*rho;
	model->add(RngDistance);
	RngProb = conProb == 1;
	model->add(RngProb);


	delete[] q;


} //end BuildRowModel

 void Decomposition::ForwardPassX(IloNumArray2* x, IloNumArray2* rng_pi, 
	 vector<GNode>* Gs, vector<gNode>* gs, vector<PNode>* ps,
	 int* objSense, IloNum* LB, IloNumArray* subObj, const int nCut, const ProblemType problem_type, const CutType cut_type,
	 const bool IsCondEff, const int pPeriod, const int pIx, const int cIx, const IloNumArray* Best_x)
 {
	 int j, s, w, ix;

	 IloEnv env=x->getEnv();

	 const int numStage = getNumPeriods();
	 const int numScenNode = numNodesInstance();

	 
	 int w_Parent = getNodeOverallIndex(pPeriod, pIx);
	 int w_Child = getNodeOverallIndex(pPeriod + 1, cIx);
	 int j_Child = getChildIndex(pPeriod + 1, cIx);

	 if (IsCondEff && pPeriod>0) {
		 int GrandParentIndex, GrandParentPeriod;
		 int PARENTstat = 0;
		 PARENTstat = getParentScenarioIndex(pPeriod, pIx, &GrandParentPeriod, &GrandParentIndex);
		 int w_GrandParent = getNodeOverallIndex(GrandParentPeriod, GrandParentIndex);
		 (*x)[w_GrandParent] = *Best_x;
	 }

	 int os=0;
	 
	 for (s = pPeriod; s < numStage; s++) {

		 //Build  period #s base model
		 SparseCSC  Wmat;
		 char * rngSense;
		 double * b, *c, *bdl, *bdu;
		 

		 IloModel MODEL(env);
		 IloCplex CPX(env);
		 IloNumVarArray X(env);
		 IloNumVar Alpha(env);
		 IloObjective OBJ(env);
		 IloRangeArray RNG(env);
		 
		 BuildBaseSubModel(&MODEL, &X, &Alpha, &OBJ, &RNG, &Wmat, &rngSense, &b, &c, &bdl, &bdu, &os, s, problem_type);
		 

		 *objSense = os;
		 ix = 0;
		 while (ix<numScenarios(s)){		

			 w = getNodeOverallIndex(s, ix);

			 bool IsInParentSubtree = true;
			 bool IsInChildSubtree = false;
			 if (IsCondEff) {
				 IsInParentSubtree = (w_Parent == w) || (getDepth(pPeriod, pIx, s, ix) > 0);
				 IsInChildSubtree = (w_Child == w) || (getDepth(pPeriod + 1, cIx, s, ix) > 0);
			 }

			 if (IsInParentSubtree && !IsInChildSubtree) {// do the forward pass on the subtree of #pIx and not in the subtree of #cIx	 

			//Add Theta variables for node #w
				 IloNumVarArray Theta(env);
				 UpdateThetaVariables(&Theta, &OBJ, s, ix, problem_type, cut_type);

				 //Update scenario information for node #ix (stochasticity on c, b, W, T), given the passed x #x_hat from its ancestor
				 //x from the ancestor of node #ix
				 IloNumArray x_hat(env);
				 int parentperiod, parentIx;
				 int PARENTstat = 0;

				 if (s > 0) {
					 PARENTstat = getParentScenarioIndex(s, ix, &parentperiod, &parentIx);
					 int oIx = getNodeOverallIndex(parentperiod, parentIx);
					 x_hat = (*x)[oIx];
					 UpdateScenario(&X, &OBJ, &RNG, &Wmat, rngSense, &b, &x_hat, s, ix);
				 }

				 //Update previous cuts and add them to model #MODEL
				 IloRangeArray optCut(env);
				 IloRangeArray maxCut(env);
				 if (s < numStage - 1 && nCut>0) {
					 if (w == w_Parent && IsCondEff)
						 UpdateOldCuts(&(*Gs)[w], &(*gs)[w], &optCut, &(*ps)[w], &maxCut, &Wmat.NumCols, &os, &X, &Alpha, &Theta, s, ix, nCut,
							 problem_type, cut_type, IsCondEff, j_Child);
					 else
						 UpdateOldCuts(&(*Gs)[w], &(*gs)[w], &optCut, &(*ps)[w], &maxCut, &Wmat.NumCols, &os, &X, &Alpha, &Theta, s, ix, nCut,
							 problem_type, cut_type);

					 MODEL.add(optCut);
					 if (problem_type == DRSO)
						 MODEL.add(maxCut);
				 }

				 //Create Cplex object #CPLEX and solve model #MODEL
				 CPX.extract(MODEL);
				 CPX.solve();



				 if (!CPX.solve()) {
					 env.error() << "Problem in Period " << s << " is Infeasible" << endl;
					 throw(-1);
				 }


				 //Store lower bound and solution #x
				 if (s == 0) {
					 *LB = CPX.getObjValue();
				 }


				 IloNumArray xsol(env);
				 CPX.getValues(xsol, X);
				 (*x)[w] = xsol;


				 //Store cost and dual variables at last stage	
				 if (s == numStage - 1) {

					 (*subObj)[w] = CPX.getObjValue();

					 IloNumArray RNG_pisol(env);
					 CPX.getDuals(RNG_pisol, RNG);
					 (*rng_pi)[w] = RNG_pisol;



				 }
				 else {
					 if (problem_type == DRSO) {
						 (*subObj)[w] = CPX.getObjValue() - CPX.getValue(Alpha);
					 }
					 else if (cut_type == SINGLE) {
						 (*subObj)[w] = CPX.getObjValue() - CPX.getValue(Theta[0]);
					 }
					 else {
						 (*subObj)[w] = CPX.getObjValue();

						 int ChildPeriod = 0;
						 int numChild = 0;
						 int * ChildIndices = NULL;
						 int CHILDstat = 0;

						 //Gets children scenarios of scenario node #ix in period #s
						 CHILDstat = getChildScenarioIndices(s, ix, &ChildPeriod, &numChild, &ChildIndices);
						 for (j = 0; j < numChild; j++)
							 (*subObj)[w] -= getCondProbability(ChildPeriod, ChildIndices[j])* CPX.getValue(Theta[j]);


					 }
					 //now remove added cuts so that MODEL can be used for other nodes in period #s-1

					 if (nCut > 0) {
						 MODEL.remove(optCut);
						 optCut.clear();
						 optCut.end();
						 if (problem_type == DRSO) {
							 MODEL.remove(maxCut);
							 maxCut.clear();
							 maxCut.end();
						 }
					 }

				 }


				 if (problem_type == RISK_NEUTRAL) {
					 for (j = 0; j < Theta.getSize(); j++)
						 OBJ.setLinearCoef(Theta[j], 0);
				 }
			}//end subtree
			ix++;
			RetrieveBaseSubModel(&X, &OBJ, &RNG, &Wmat, rngSense, &b, &c, &bdl, &bdu);
		 }//end node #w
		 X.clear();
		 X.end();
		 Alpha.end();
		 OBJ.end();
		 RNG.clear();
		 RNG.end();
		 MODEL.end();
		 CPX.end();

	 }//end stage #s


 }//end ForwardPassX


 void Decomposition::CalculateFuncVal(IloNumArray* subObj_hat)
 {
	 int ix, s, w, j;
	 IloEnv env= subObj_hat->getEnv();

	 const int numStage = getNumPeriods();

	 for (s = numStage - 1; s > 0; s--) {
		 ix = 0;

		 while (ix < numScenarios(s - 1)) {

			 w = getNodeOverallIndex(s - 1, ix);

			 int ChildPeriod = 0;
			 int numChild = 0;
			 int * ChildIndices = NULL;
			 int CHILDstat = 0;
			 
			 //Gets children scenarios of scenario #w in period #s-1
			 CHILDstat = getChildScenarioIndices(s - 1, ix, &ChildPeriod, &numChild, &ChildIndices);
			 if (CHILDstat != 0) {
				 cerr << "getChildScenarioIndices: error" << endl;
				 return;
			 }
			 IloNumArray  h(env, numChild);
			 int * oIx = new int[numChild];
			 double * q = new double[numChild];
			 
			 for (j = 0; j < numChild; j++) {
				 oIx[j] = getNodeOverallIndex(ChildPeriod, ChildIndices[j]);
				 h[j] = (*subObj_hat)[oIx[j]];
				 q[j] = getCondProbability(ChildPeriod, ChildIndices[j]);
				 (*subObj_hat)[w] += q[j] * h[j];
			 }
			 
			 delete[] oIx;
			 delete[] q;
			 ix++;
		 }
	 }
	

 }//end CalculateFuncVal

 void Decomposition::ForwardPassP(vector<cut_P>* ps, IloNumArray * pArray, IloNumArray* subObj_hat, const int* objSense,
	 const int nCut, bool* contForwardPass, const ProblemType problem_type, const CutType cut_type, const double rho,
	 const bool IsCondEff, const bool IsPathEff,
	 const int pPeriod, const int pIx, const int cIx, const int child)
 {
	 int ix, s, w, j;
	 IloEnv env = subObj_hat->getEnv();


	 const int numStage = getNumPeriods();
	 const int numScenNode = numNodesInstance();
	 const int numParents = numTotalNodesInstance(numStage - 2);
	 IloNumArray p(env, numScenNode);


	 bool IsAssess = IsCondEff || IsPathEff;

	 int w_Parent = getNodeOverallIndex(pPeriod, pIx);
	 int w_Child = getChildOverallIndex(pPeriod, pIx, child);

	 if (IsCondEff)
		 p[w_Parent] = 1;
	 else
		 p[0] = 1;

	 int begin_stage = 0;
	 if (IsCondEff)
		 begin_stage = pPeriod;


	 for (s = numStage - 1; s > begin_stage && *contForwardPass; s--) {
		 ix = 0;


		 while (ix < numScenarios(s - 1) && *contForwardPass) {
			 w = getNodeOverallIndex(s - 1, ix);
			 bool IsInParentSubtree = true;
			 bool IsInChildSubtree = false;
			 if (IsCondEff) {
				 IsInParentSubtree = (w_Parent == w) || (getDepth(pPeriod, pIx, s - 1, ix) > 0);
				 IsInChildSubtree = (w_Child == w) || (getDepth(pPeriod + 1, cIx, s - 1, ix) > 0);
			 }

			 if (IsInParentSubtree && !IsInChildSubtree) {// do the pass on the subtree of #pIx and not in the subtree of #cIx	 
				 int ChildPeriod = 0;
				 int numChild = 0;
				 int * ChildIndices = NULL;
				 int CHILDstat = 0;

				 //Gets children scenarios of scenario #w in period #s-1
				 CHILDstat = getChildScenarioIndices(s - 1, ix, &ChildPeriod, &numChild, &ChildIndices);
				 if (CHILDstat != 0) {
					 cerr << "getChildScenarioIndices: error" << endl;
					 return;
				 }

				 IloNumArray  h(env, numChild);
				 int * oIx = new int[numChild];
				 for (j = 0; j < numChild; j++) {
					 oIx[j] = getNodeOverallIndex(ChildPeriod, ChildIndices[j]);
					 h[j] = (*subObj_hat)[oIx[j]];
				 }

				 IloNumVarArray P(env);
				 IloModel model(env);
				 IloCplex cplex(env);

				 BuildRowModel(&model, &P, &h, objSense, s - 1, ix, rho);

				 if (IsAssess && w_Parent == w) {
					 IloRange RngZero;
					 RngZero = P[child] == 0;
					 model.add(RngZero);
				 }

				 cplex.extract(model);
				 cplex.solve();
				 if (cplex.solve()) {
					 for (j = 0; j < numChild; j++) {
						 p[oIx[j]] = cplex.getValue(P[j]);
						 (*ps)[w].push_back(p[oIx[j]]);
						 (*subObj_hat)[w] += p[oIx[j]] * h[j];
					 }

				 }
				 else {
					 *contForwardPass = false;
				 }
				 P.end();
				 model.end();
				 cplex.end();
				 delete[] oIx;
				 h.clear();
				 h.end();
			 }
			 ix++;
		 }
	 }
	 pArray->add(p);
	 p.clear();
	 p.end();
	 
 } //endForwardPassP


 void Decomposition::BackwardPass(IloNumArray2* x, IloNumArray2* rng_pi, 
	 vector<GNode>* Gs, vector<gNode>* gs, vector<PNode>* ps,
	 const int nCut, const ProblemType problem_type, const CutType cut_type,
	 const bool IsCondEff, const int pPeriod, const int pIx, const int cIx, const int child)
 {
	 int ix, w, j, s, k, col;
	
	 IloEnv env=x->getEnv(); 

	 const int numStage = getNumPeriods();
	 const int numParents = numTotalNodesInstance(numStage - 2);

	 //pi_optCut[w] is the array of all cuts so far for node #w: 
	 //e.g. k*numChild +j ::  gives the (k+1)-th cut for from j-th child of node #w in case of DRSO or multi-cut
	 //e.g. k : gives the (k+1)-th cut for node #w in case of single-cut
	 IloNumArray2  pi_optCut(env, numParents);
	 for (w = 0; w < numParents; w++)
		 pi_optCut[w] = IloNumArray(env);
	 	 
	 
	 int w_Parent = getNodeOverallIndex(pPeriod, pIx);
	 int w_Child = getChildOverallIndex(pPeriod, pIx, child);

	 for (s = numStage - 1; s>pPeriod; s--) {
		 
		 //Build  period #s-1 base model
		 
		 SparseCSC  Wmat;
		 char * rngSense;
		 double * b, *c, *bdl, *bdu;
		 int os;
		 IloModel MODEL(env);
		 IloCplex CPX(env);
		 IloNumVarArray X(env);
		 IloNumVar Alpha(env);
		 IloObjective OBJ(env);
		 IloRangeArray RNG(env);
		 

		 BuildBaseSubModel(&MODEL, &X, &Alpha, &OBJ, &RNG, &Wmat, &rngSense, &b, &c, &bdl, &bdu, &os, s-1, problem_type);
		 ix = 0; 
		 while (ix < numScenarios(s-1)) {
			 w = getNodeOverallIndex(s - 1, ix);
			 
			 bool IsInParentSubtree = true;
			 bool IsInChildSubtree = false;
			 if (IsCondEff) {
				 IsInParentSubtree = (w_Parent == w) || (getDepth(pPeriod, pIx, s - 1, ix) > 0);
				 IsInChildSubtree = (w_Child == w) || (getDepth(pPeriod + 1, cIx, s - 1, ix) > 0);
			 }

			 if (IsInParentSubtree && !IsInChildSubtree) {// do the pass on the subtree of #pIx and not in the subtree of #cIx	 

			//Add Theta variables for node #w
				 IloNumVarArray Theta(env);
				 UpdateThetaVariables(&Theta, &OBJ, s - 1, ix, problem_type, cut_type);

				 //Update scenario information for node #w (stochasticity on c, b, W, T), given the passed x #x_hat from its ancestor
				 //x from the ancestor of node #w
				 IloNumArray x_hat(env);
				 int parentperiod, parentIx = 0;
				 int PARENTstat = 0;
				 if (s > pPeriod + 1) {
					 PARENTstat = getParentScenarioIndex(s - 1, ix, &parentperiod, &parentIx);
					 int oIx = getNodeOverallIndex(parentperiod, parentIx);
					 x_hat = (*x)[oIx];
					 UpdateScenario(&X, &OBJ, &RNG, &Wmat, rngSense, &b, &x_hat, s - 1, ix);
				 }
				 /*x_hat.clear();
				 x_hat.end();*/

				 //Update previous cuts and add them to model #MODEL
				 IloRangeArray OldoptCut(env);
				 IloRangeArray OldmaxCut(env);
				 if (s > pPeriod + 1 && nCut > 0) {
					 if (w == w_Parent && IsCondEff)
						 UpdateOldCuts(&(*Gs)[w], &(*gs)[w], &OldoptCut, &(*ps)[w], &OldmaxCut, &Wmat.NumCols, &os, &X, &Alpha, &Theta, s - 1, ix, nCut,
							 problem_type, cut_type, IsCondEff, child);
					 else
						 UpdateOldCuts(&(*Gs)[w], &(*gs)[w], &OldoptCut, &(*ps)[w], &OldmaxCut, &Wmat.NumCols, &os, &X, &Alpha, &Theta, s - 1, ix, nCut,
							 problem_type, cut_type);

					 MODEL.add(OldoptCut);
					 if (problem_type == DRSO)
						 MODEL.add(OldmaxCut);
				 }

				 CPX.extract(MODEL);

				 //generate cut coefficients for children of node #w {g(child), G(child), and P(child)}
				 int childperiod = 0;
				 int numChild = 0;
				 int * ChildIndices = NULL;
				 int CHILDstat = 0;
				 //get children scenarios of scenario node #ix in period #s-1
				 CHILDstat = getChildScenarioIndices(s - 1, ix, &childperiod, &numChild, &ChildIndices);

				 if (CHILDstat != 0) {
					 cerr << "getChildScenarioIndices: error" << endl;
					 return;
				 }

				 cut_g gw;
				 cut_G Gw;
				 cut_P Pw;

				 //overall index 
				 int * cIx = new int[numChild];

				 for (j = 0; j < numChild; j++) {
					 GNodeType G;
					 gNodeType g;
					 double * bChild;
					 double * loChild;
					 double * upChild;
					 SparseCSC TChild;
					 //Get RHS and T of children of node #w
					 GetScenarioRHSTechnologyMatrix(&bChild, &loChild, &upChild, &TChild, childperiod, ChildIndices[j]);

					 if (w != w_Parent || !IsCondEff || j != child) {
						 cIx[j] = getNodeOverallIndex(childperiod, ChildIndices[j]);

						 IloNumArray rng_piChild(env);
						 rng_piChild = (*rng_pi)[cIx[j]];



						 gNode gChild;
						 IloNumArray piOptCutChild(env);//array for cuts

						 if (childperiod < numStage - 1)
						 {
							 gChild = (*gs)[cIx[j]];
							 piOptCutChild.add(pi_optCut[cIx[j]]);
						 }


						 GenerateOptCutCoefficients(&G, &g, &TChild, bChild, loChild, upChild, &rng_piChild, &gChild, &piOptCutChild,
							 childperiod, ChildIndices[j], nCut, problem_type, cut_type);

					 }
					 gw.push_back(g);
					 Gw.push_back(G);

				 }
				 delete[] cIx;

				 (*gs)[w].push_back(gw);
				 (*Gs)[w].push_back(Gw);



				 if (problem_type == DRSO)
					 Pw = (*ps)[w][nCut];

				 //use cut coefficients to form the cut for node #ix in current iteration
				 IloRangeArray NewoptCut(env);
				 IloRange NewmaxCut;
				 if (w == w_Parent && IsCondEff)
					 GenerateCuts(&Gw, &gw, &NewoptCut, &Pw, &NewmaxCut, &Wmat.NumCols, &os, &X, &Alpha, &Theta,
						 nCut, ix, s - 1, problem_type, cut_type, IsCondEff, child);
				 else
					 GenerateCuts(&Gw, &gw, &NewoptCut, &Pw, &NewmaxCut, &Wmat.NumCols, &os, &X, &Alpha, &Theta,
						 nCut, ix, s - 1, problem_type, cut_type);


				 if (s > pPeriod + 1) {
					 MODEL.add(NewoptCut);
					 if (problem_type == DRSO)
						 MODEL.add(NewmaxCut);
				 }


				 //Create Cplex object #CPLEX and solve model #MODEL 
				 if (s > pPeriod + 1) {
					 CPX.extract(MODEL);
					 CPX.solve();
					 //Dual values for constraints
					 IloNumArray rng_pisol(env);
					 CPX.getDuals(rng_pisol, RNG);
					 (*rng_pi)[w] = rng_pisol;



					 //Dual values for optimality old and new cuts
					 if (nCut > 0) {
						 IloNumArray pi_OldOptCut(env);
						 CPX.getDuals(pi_OldOptCut, OldoptCut);
						 pi_optCut[w].add(pi_OldOptCut);
						 pi_OldOptCut.clear();
						 pi_OldOptCut.end();
					 }

					 IloNumArray pi_NewOptCut(env);
					 CPX.getDuals(pi_NewOptCut, NewoptCut);
					 pi_optCut[w].add(pi_NewOptCut);
					 pi_NewOptCut.clear();
					 pi_NewOptCut.end();
					 //now remove added cuts so that MODEL can be used for other nodes in period #s-1

					 MODEL.remove(NewoptCut);
					 NewoptCut.clear();
					 NewoptCut.end();
					 if (problem_type == DRSO) {
						 MODEL.remove(NewmaxCut);
						 NewmaxCut.end();
					 }

					 if (nCut > 0) {
						 MODEL.remove(OldoptCut);
						 OldoptCut.clear();
						 OldoptCut.end();
						 if (problem_type == DRSO) {
							 MODEL.remove(OldmaxCut);
							 OldmaxCut.clear();
							 OldmaxCut.end();
						 }
					 }

				 }
				 for (j = 0; j < Theta.getSize(); j++)
					 OBJ.setLinearCoef(Theta[j], 0);

			 }
			 ix++;
			 RetrieveBaseSubModel(&X, &OBJ, &RNG, &Wmat, rngSense, &b, &c, &bdl, &bdu);
		 }
		 X.clear();
		 X.end();
		 Alpha.end();
		 OBJ.end();
		 RNG.clear();
		 RNG.end();
		 MODEL.end();
		 CPX.end();
		 
	 }
	 pi_optCut.clear();
	 pi_optCut.end();
	 
 
 }// end BackwardPass


 void Decomposition::UpdateScenario(IloNumVarArray* X, IloObjective* obj, IloRangeArray* rng, 
	 SparseCSC* Wmat, const char * rngSense, 
	 double ** b, const IloNumArray* xhat, const int period, const int ix)
 {

	 IloEnv env = X->getEnv();

	 int col, row; //dimension of X variables and constraints counters
	 double val; //value of a coefficient
	 double bound; //current value of a bound in a range
	 int change; //to loop over elements that should be changed in a scenario node
	 int ent;

	 //To hold Technology Matrix info
	 SparseCSC techMat;

	 int TECHstat = 0;
	 //get subproblem Technology matrix
	 TECHstat = getTechnologyMatrix(period, &techMat.NumCols, &techMat.NumRows,
		 &techMat.NumEntries, &techMat.col_ptr, &techMat.row_ind,
		 &techMat.val);

	 if (TECHstat != 0) {
		 cerr << "getTechnologyMatrix: error" << endl;
		 return;
	 }

	 //take care of b portion
	 IloNumArray rhs(env, techMat.NumRows);
	 for (row = 0; row < techMat.NumRows; row++)
		 rhs[row] = (*b)[row];

	 //take care of Tx portion
	 IloNumArray Tx_temp(env, techMat.NumRows);
	 for (row = 0; row < techMat.NumRows; row++)
		 Tx_temp[row] = 0;


	 if (techMat.NumCols != xhat->getSize()) {
		 cerr << "UpdateScenario- Wrong size xhat passed: error" << endl;
		 return;
	 }

	 double q = 0;
	 int numChanges = 0;
	 ChangeType *entityChange = NULL;
	 int *entryChange = NULL;
	 double *valueChange = NULL;

	 int SCENstat = 0;
	 //get scenario information
	 SCENstat = getScenario(period, ix, &q, &numChanges, &entityChange, &entryChange, &valueChange);

	 if (SCENstat != 0) {
		 cerr << "getScenario: error" << endl;
		 return;
	 }

	 //cout << *xhat << endl;
	 for (change = 0; change < numChanges; change++) {
		 if (entityChange[change] == COST) {
			 col = entryChange[change];
			 val = valueChange[change];
			 obj->setLinearCoef((*X)[col], val);
		 }//COST

		 if (entityChange[change] == RHS) {
			 row = entryChange[change];
			 val = valueChange[change];
			 //rhs_changed[row] = 1;
			 rhs[row] = val;
		 } //RHS


		   //here sense should be outside varibles sense[period][row]
		 if (entityChange[change] == RANGE) {
			 row = entryChange[change];
			 val = valueChange[change];
			 switch (rngSense[row]) {
			 case 'E':
			 case 'e':
				 if (val > 0) {
					 bound = (*rng)[row].getLB();
					 //?
					 (*rng)[row].setBounds(bound, bound + abs(val));
				 }
				 else {
					 bound = (*rng)[row].getUB();
					 (*rng)[row].setBounds(bound - abs(val), bound);
				 }
				 break;
			 case 'L':
			 case 'l':
				 bound = (*rng)[row].getUB();
				 (*rng)[row].setBounds(bound - abs(val), bound);
				 break;
			 case 'G':
			 case 'g':
				 bound = (*rng)[row].getLB();
				 (*rng)[row].setBounds(bound, bound + abs(val));
				 break;
			 }
		 } //RANGE

		 if (entityChange[change] == LB) {
			 col = entryChange[change];
			 val = valueChange[change];
			 (*X)[col].setLB(val);
		 }//LB


		 if (entityChange[change] == UB) {
			 col = entryChange[change];
			 val = valueChange[change];
			 (*X)[col].setUB(val);
		 }//UB

		 if (entityChange[change] == TT) {
			 ent = entryChange[change];
			 val = valueChange[change];
			 techMat.val[ent] = val;
		 }//TT technology matrix


		 if (entityChange[change] == WW) {
			 ent = entryChange[change];
			 //if (ent != -1) {
			col = FindEntryColInCSCMatrix(&Wmat->NumCols, &ent, &Wmat->col_ptr);
			row = Wmat->row_ind[ent];
			val = valueChange[change];
			(*rng)[row].setLinearCoef((*X)[col], val);
			 //}
		 }//WW recourse matrix
	 }

	 //set rhs
	 for (row = 0; row < techMat.NumRows; row++){
		 for (ent = 0; ent < techMat.NumEntries; ent++){
			 if (techMat.row_ind[ent] == row) {
				 col = FindEntryColInCSCMatrix(&techMat.NumCols, &ent, &techMat.col_ptr);
				 Tx_temp[row] += (*xhat)[col] * techMat.val[ent];
			 }
		 }
	 }
		 


	 for (row = 0; row<techMat.NumRows; row++) {
		 
		 switch (rngSense[row]) {
		 case 'E':
		 case 'e':
			 (*rng)[row].setBounds(rhs[row] - Tx_temp[row], rhs[row] -Tx_temp[row]);
			break;
		 case 'L':
		 case 'l':
			 (*rng)[row].setUB(rhs[row] - Tx_temp[row]);
			 break;
		 case 'G':
		 case 'g':
			 (*rng)[row].setLB(rhs[row] - Tx_temp[row]);
			 break;
		 }

	 }
	 Tx_temp.clear();
	 Tx_temp.end();


 }//end UpdateScenario


 void Decomposition::GetScenarioRHSTechnologyMatrix(double ** b, double ** blo, double ** bup, SparseCSC* Tmat, 
	 const int period, const int ix)
 {

	 int col, row; //dimension of X variables and constraints counters
	 double val; //value of a coefficient
	 int change; //to loop over elements that should be changed in a scenario node
	 int i, ent;


	 //To hold LP info
	 int os = 0;
	 SparseCSC recourseMat;
	 double *cost = NULL;
	 double *rhs = NULL;
	 char *sense = NULL;
	 double *bdl = NULL;
	 double *bdu = NULL;
	 int * coltype = NULL;


	 int LPstat = 0;
	 //get b for Base LP 
	 LPstat = getBaseLP(period, &recourseMat.NumCols, &recourseMat.NumRows,
		 &recourseMat.NumEntries, &os, &cost, &rhs,
		 &sense, &recourseMat.col_ptr, &recourseMat.row_ind,
		 &recourseMat.val, &bdl, &bdu,
		 &coltype);

	 //take care of b portion
	 //(*b)[row]
	 *b = rhs;
	 *blo = bdl;
	 *bup = bdu;


	 
	 if (LPstat != 0) {
		 cerr << "getBaseLP: error" << endl;
		 return;
	 }

	 //get T for Base LP
	 SparseCSC TmatBase;

	 int TECHstat = 0;
	 //get subproblem Technology matrix
	 TECHstat = getTechnologyMatrix(period, &TmatBase.NumCols, &TmatBase.NumRows,
		 &TmatBase.NumEntries, &TmatBase.col_ptr, &TmatBase.row_ind,
		 &TmatBase.val);

	 if (TECHstat != 0) {
		 cerr << "getTechnologyMatrix: error" << endl;
		 return;
	 }

	 *Tmat = TmatBase;

	 double q = 0;
	 int numChanges = 0;
	 ChangeType *entityChange = NULL;
	 int *entryChange = NULL;
	 double *valueChange = NULL;

	 int SCENstat = 0;
	 //get scenario information
	 SCENstat = getScenario(period, ix, &q, &numChanges, &entityChange, &entryChange, &valueChange);

	 if (SCENstat != 0) {
		 cerr << "getScenario: error" << endl;
		 return;
	 }

	 for (change = 0; change < numChanges; change++) {

		 if (entityChange[change] == RHS) {
			 row = entryChange[change];
			 val = valueChange[change];
			 (*b)[row] = val;
		 } //RHS

		 if (entityChange[change] == TT) {
			 ent = entryChange[change];
			 val = valueChange[change];
			 Tmat->val[ent] = val;
			/* if (val != 0)
				 Tmat->val[ent] = val;
			 else
				 RemoveEntryCSCMatrix(Tmat, ent);*/
			 
		 }//TT technology matrix


		 if (entityChange[change] == LB) {
			 col = entryChange[change];
			 val = valueChange[change];
			 (*blo)[col] = val;
		 }//LB


		 if (entityChange[change] == UB) {
			 col = entryChange[change];
			 val = valueChange[change];
			 (*bup)[col] = val;
		 }//UB



	 } 
 }//end GetScenarioInfo


 void Decomposition::UpdateThetaVariables(IloNumVarArray* Theta, IloObjective* obj, const int period, const int ix,
	 const ProblemType problem_type, const CutType cut_type)
 {
	 IloEnv env = Theta->getEnv();
	 
	 char varName[100];

	 const int numStage = getNumPeriods();
	 int child;

	 //add Theta's to objective function if risk-neutral
	 if (period < numStage - 1) {

		 int childperiod = 0;
		 int numChild = 0;
		 int * ChildIndices = NULL;
		 int CHILDstat = 0;
		 //get children scenarios of scenario node #ix in period #period
		 CHILDstat = getChildScenarioIndices(period, ix, &childperiod, &numChild, &ChildIndices);

		 if (CHILDstat != 0) {
			 cerr << "getChildScenarioIndices: error" << endl;
			 return;
		 }

		 if (problem_type == DRSO || cut_type == MULTI) {

			IloNumVarArray ThetaTemp(env, numChild, -MYINFINITY, MYINFINITY);
			*Theta = ThetaTemp;
			double * q=new double[numChild];
			//get condtional probabilities of children scenario
			for (child = 0; child < numChild; child++) {
				sprintf_s(varName, "Theta.%d", (int)child);
				(*Theta)[child].setName(varName);
				if (problem_type != DRSO) {
					q[child] = getCondProbability(childperiod, ChildIndices[child]);
					obj->setLinearCoef((*Theta)[child], q[child]);
				}
			}
			delete[] q;
			/*ThetaTemp.clear();
			ThetaTemp.end();*/
		 }
		 if  (cut_type == SINGLE){
			 IloNumVarArray ThetaTemp(env, 1, -MYINFINITY, MYINFINITY);
			 *Theta = ThetaTemp;
			 sprintf_s(varName, "Theta.%d", 0);
			 (*Theta)[0].setName(varName);
			 if (problem_type != DRSO)
				obj->setLinearCoef((*Theta)[0], 1);

			 /*ThetaTemp.clear();
			 ThetaTemp.end();*/
		}
	 }//end period


 } //end UpdateThetaVariables

 
 void Decomposition::UpdateOldCuts(GNode* Gs, gNode* gs, IloRangeArray* optCut, PNode* ps, IloRangeArray* maxCut,
	 const int * ncols, const int *objSense,
	 const IloNumVarArray* X, const IloNumVar* Alpha, const IloNumVarArray* Theta,
	 const int period, const int ix, const int nCut, 
	 const ProblemType problem_type, const CutType cut_type, const bool IsCondEff, const int child)
 {
	 
	 IloEnv env = maxCut->getEnv();

	 const int numStage = getNumPeriods();
	 
	 int childperiod = 0;
	 int numChild = 0;
	 int * ChildIndices = NULL;
	 int CHILDstat = 0;
	 //get children scenarios of scenario node #ix in period #period
	 CHILDstat = getChildScenarioIndices(period, ix, &childperiod, &numChild, &ChildIndices);

	 if (CHILDstat != 0) {
		 cerr << "getChildScenarioIndices: error" << endl;
		 return;
	 }

	 int k, j, col;

	 //initialize optCut and maxCut
	 if (problem_type == DRSO) {
		 IloRangeArray maxCut_hat(env, nCut, -IloInfinity, IloInfinity);
		 *maxCut = maxCut_hat;
		 /*maxCut_hat.clear();
		 maxCut_hat.end();*/
	 }
	 if (problem_type == DRSO || cut_type == MULTI) {
		 IloRangeArray optCut_hat(env, nCut*numChild, -IloInfinity, IloInfinity);
		 *optCut = optCut_hat;
		 /*optCut_hat.clear();
		 optCut_hat.end();*/
	 }
	

	 if (cut_type == SINGLE) {
		 IloRangeArray optCut_hat(env, nCut, -IloInfinity, IloInfinity);
		 *optCut = optCut_hat;
		 /*optCut_hat.clear();
		 optCut_hat.end();*/
	 }

	 

	if (problem_type == DRSO) {

		for (k = 0; k < nCut; k++) {

			(*maxCut)[k].setLinearCoef(*Alpha, 1);
			if (*objSense == MINIMIZE)
				(*maxCut)[k].setBounds(0, IloInfinity);
			else
				(*maxCut)[k].setBounds(-IloInfinity, 0);

			for (j = 0; j < numChild; j++) {
				if (!IsCondEff || j != child)
					(*maxCut)[k].setLinearCoef((*Theta)[j], -(*ps)[k][j]);
			}
				

		}
	}//end if 

	if (problem_type == DRSO || cut_type == MULTI) {

		for (k = 0; k < nCut; k++) {
			 
			for (j = 0; j < numChild; j++) {
				
				if (*objSense == MINIMIZE)
					(*optCut)[k*numChild + j].setBounds((*gs)[k][j], IloInfinity);
				else
					(*optCut)[k*numChild + j].setBounds(-IloInfinity, (*gs)[k][j]);

				(*optCut)[k*numChild + j].setLinearCoef((*Theta)[j], 1);


				//take care of G portion 
				for (col = 0; col < *ncols; col++)
					(*optCut)[k*numChild + j].setLinearCoef((*X)[col], (*Gs)[k][j][col]);
			}
		}

	}//end if

	if (cut_type == SINGLE) {
		double * q = new double[numChild];
		for (j = 0; j < numChild; j++)
			q[j] = getCondProbability(childperiod, ChildIndices[j]);

		for (k = 0; k < nCut; k++) {

			//take care of g portion
			double gAgg=0;
			for (j = 0; j < numChild; j++)
				gAgg += q[j] * (*gs)[k][j];
			 
			
			if (*objSense == MINIMIZE) 
				(*optCut)[k].setBounds(gAgg, IloInfinity);
			else
				(*optCut)[k].setBounds(-IloInfinity, gAgg);

			(*optCut)[k].setLinearCoef((*Theta)[0], 1);

			//take care of G portion
			double *  GAgg=new double[*ncols];
			for (col = 0; col < *ncols; col++) {
				GAgg[col] = 0;
				for (j = 0; j < numChild; j++)
					GAgg[col] += q[j] * (*Gs)[k][j][col];

				(*optCut)[k].setLinearCoef((*X)[col], GAgg[col]);
			}
			delete[] GAgg;
		}
	delete[] q;
	
	}//end if

 } //end UpdateOldCuts
 
 void Decomposition::GenerateOptCutCoefficients(GNodeType * G, gNodeType *  g, const SparseCSC * Tmat, const double * b, const double * bl, const double * bu,
	 const IloNumArray* rng_pi, gNode* gChild, const IloNumArray* pi_optCut, const int period, const int ix, const int nCut,
	 const ProblemType problem_type, const CutType cut_type)
 {

 
	 int row, col; //row and col counters
	 double val=0; //entry value
	 int entry; //entry counter
	 int child, k; 
	

	//store cut intercept
	 //take care of b
	 gNodeType  gTemp=0;

	 for (row = 0; row< Tmat->NumRows; row++) 
		 gTemp += (*rng_pi)[row] * b[row];

	 

	//now, take care of the rest of g
	 const int numStage = getNumPeriods();

	 int childperiod = 0;
	 int numChild = 0;
	 int * ChildIndices = NULL;
	 int CHILDstat = 0;
	 //get children scenarios of scenario node #ix in period #period
	 if (period  < numStage - 1) {
		 CHILDstat = getChildScenarioIndices(period, ix, &childperiod, &numChild, &ChildIndices);

		 if (CHILDstat != 0) {
			 cerr << "getChildScenarioIndices: error" << endl;
			 return;
		 }
	 }

	 if (problem_type == DRSO || cut_type == MULTI) {
		 if (period < numStage - 1) {
			 for (k = 0; k < nCut + 1; k++) {
				 for (child = 0; child < numChild; child++) 
					 gTemp += (*pi_optCut)[k*numChild + child] * (*gChild)[k][child];
			 }
		 }
	 }


	 if (cut_type == SINGLE) {
		 
		 if (period  < numStage - 1) {
			 
			 double * q = new double[numChild];
			 
			 for (child = 0; child < numChild; child++)
				 q[child] = getCondProbability(childperiod, ChildIndices[child]);

			 for (k = 0; k < nCut + 1; k++) {
				 double gAgg = 0;
				 
				 for (child = 0; child < numChild; child++) 
					 gAgg += q[child] * (*gChild)[k][child];
				 
				 gTemp += (*pi_optCut)[k] * gAgg;
			 }
			 delete[] q;
		 }
	 }

	 *g = gTemp;
	 
	 //store cut slope
	
	 GNodeType GTemp (Tmat->NumCols);
	 for (col = 0; col < Tmat->NumCols; col++) {
		 GTemp[col] = 0;
		 for (entry = Tmat->col_ptr[col]; entry < Tmat->col_ptr[col + 1]; entry++) {
			 row = Tmat->row_ind[entry];
			 val = Tmat->val[entry];
			 GTemp[col] += (*rng_pi)[row] * val;
		 }
	 }
	 *G = GTemp;
 }//end GenerateOptCutCoefficients


 void Decomposition::GenerateCuts(const cut_G * G, const cut_g *  g, IloRangeArray * optCut, const cut_P* p, IloRange * maxCut, 
	 const int * ncols, const int *objSense,
	const IloNumVarArray*  X, const IloNumVar*  Alpha, const IloNumVarArray*  Theta,
	 const int nCut, const int ix, const int period, const ProblemType problem_type, const CutType cut_type,
	 const bool IsCondEff, const int child)
 {
	 
	 int col, j;

	 IloEnv env= optCut->getEnv();

	 const int numScenNode = numNodesInstance();

	 int childperiod = 0;
	 int numChild = 0;
	 int * ChildIndices = NULL;
	 int CHILDstat = 0;
	 //get children scenarios of scenario node #ix in period #period
	 CHILDstat = getChildScenarioIndices(period, ix, &childperiod, &numChild, &ChildIndices);

	 if (CHILDstat != 0) {
		 cerr << "getChildScenarioIndices: error" << endl;
		 return;
	 }


	 //initialize optCut and maxCut
	 if (problem_type == DRSO) {
		 IloRange maxCut_hat(env, -IloInfinity, IloInfinity);
		 *maxCut = maxCut_hat;
	 }
	 if (problem_type == DRSO || cut_type == MULTI) {
		 IloRangeArray optCut_hat(env, numChild, -IloInfinity, IloInfinity);
		 *optCut = optCut_hat;
		 /*optCut_hat.clear();
		 optCut_hat.end();*/
	 }
	 if (cut_type == SINGLE) {
		 IloRangeArray optCut_hat(env, 1, -IloInfinity, IloInfinity);
		 *optCut = optCut_hat;
		 /*optCut_hat.clear();
		 optCut_hat.end();*/
	 }
	 //take care of max cut
	 if (problem_type == DRSO) {
		 maxCut->setLinearCoef(*Alpha, 1);
		 if (*objSense == MINIMIZE)
			 maxCut->setBounds(0, IloInfinity);
		 else
			 maxCut->setBounds(-IloInfinity, 0);

		 for (j = 0; j < numChild; j++) {
			 if (!IsCondEff || j != child)
				 maxCut->setLinearCoef((*Theta)[j], -(*p)[j]);
		 }
	 }

	 if (problem_type == DRSO || cut_type == MULTI) {
		 //take care of g portion 
		 for (j = 0; j < numChild; j++) {
			if (*objSense == MINIMIZE)
				(*optCut)[j].setBounds((*g)[j], IloInfinity);
			else
				(*optCut)[j].setBounds(-IloInfinity, (*g)[j]);

			(*optCut)[j].setLinearCoef((*Theta)[j], 1);
			 	 
		 }
		 
		 //take care of G portion
		 for (j = 0; j < numChild; j++) { 
			for (col = 0; col < *ncols; col++)
				(*optCut)[j].setLinearCoef((*X)[col], (*G)[j][col]);
		 }
		 
	 }//end if 

	 
	 if (cut_type == SINGLE) {
		 double * q=new double [numChild];
		//get condtional probabilities of children scenario
		 for (j = 0; j < numChild; j++)
			 q[j] = getCondProbability(childperiod, ChildIndices[j]);

		 //take care of g portion 
		 double gAgg=0;
		 for (j = 0; j < numChild; j++) 
			 gAgg += q[j] * (*g)[j];

		 if (*objSense == MINIMIZE)
			 (*optCut)[0].setBounds(gAgg, IloInfinity);
		 else
			 (*optCut)[0].setBounds(-IloInfinity, gAgg);

		 (*optCut)[0].setLinearCoef((*Theta)[0], 1);

		 
		 //take care of G portion
		 double * GAgg=new double [*ncols];
		 for (col = 0; col < *ncols; col++) {
			 GAgg[col] = 0;
			 for (j = 0; j < numChild; j++)
				 GAgg[col] += q[j] * (*G)[j][col];
			 (*optCut)[0].setLinearCoef((*X)[col], GAgg[col]);
		 }
		delete[] GAgg;
		delete[] q;
	 }//end if


 }// end GenerateCuts


 void Decomposition::SolveOriginal(char *problem_name, const double rho, double* objVal, IloNumArray2* xOpt, IloNumArray* SubObjOpt, IloNumArray* Worst_pOpt, const ProblemType problem_type, const CutType cut_type)
 {
	 char * problem_t=NULL;
	 if (problem_type == DRSO)
		 problem_t = "DRSO";
	 else
		 problem_t = "RISK_NEUTRAL";

	 char * cut_t=NULL;
	 if (cut_type == SINGLE)
		 cut_t = "SINGLE";
	 else
		 cut_t = "MULTI";

	 char resName[100];
	 sprintf_s(resName, "EX_Results_%s_%s_%s_%0.2f.txt", problem_name, problem_t, cut_t, rho);
	 const char* Resfilename = resName;
	 ofstream Result(Resfilename);
	 Result.precision(15);

	 int k, w;
	 Result << "Iter  \t    Z_hat    \t     LB      \t  UB  \t Gap%" << endl;

	 double rel_Gap = IloInfinity;
	 double Appx = -IloInfinity;
	 double FuncVal = IloInfinity;
	 double z_hat;
	 int nCut = 0;

	 const int numScenNode = numNodesInstance();
	 const int numStage = getNumPeriods();
	 IloEnv env;

	 ///////////////// PARAMETERS ///////////////// 
	 IloNumArray2 rng_pi(env, numScenNode);
	 const int numParents = numTotalNodesInstance(numStage - 2);
	 vector<gNode> gs (numParents);
	 vector<GNode> Gs(numParents);
	 vector<PNode> ps(numParents);
	 IloNumArray SubObj(env, numScenNode);
	 IloNumArray Worst_p(env, numScenNode);
	 IloNumArray2 sol(env, numScenNode);
	 *xOpt = sol;
	 int objSense;

	 ////////////////////////////////////////////////

	 const clock_t begin_time = clock();
	 while (rel_Gap > TOLER  && float(clock() - begin_time) / CLOCKS_PER_SEC <= TIME_LIMIT) {

		 ForwardPassX(&sol, &rng_pi, &Gs, &gs, &ps, &objSense, &Appx, &SubObj, nCut, problem_type, cut_type);
		 bool contForwardPass = true;
		 if (nCut == 0 && objSense == MAXIMIZE)
			 FuncVal = -IloInfinity;
		 
		 IloNumArray pArray(env);
		 if (problem_type == DRSO) {
			 vector<cut_P> p (numParents);
			 ForwardPassP(&p, &pArray, &SubObj, &objSense, nCut, &contForwardPass, problem_type, cut_type, rho);
			 
			 //store the solution to the inner problem for the current iteration
		

			 for (w = 0; w < numParents; w++)
				 ps[w].push_back(p[w]);

			 p.~vector();
			 
		 }
		 else {
			 CalculateFuncVal(&SubObj);
		 }


		 z_hat = SubObj[0];
		 if (objSense == MINIMIZE) {
			 if (z_hat < FuncVal) {
				 FuncVal = z_hat;
				 *Worst_pOpt=pArray;
				 *SubObjOpt= SubObj;
				 for (w = 0; w < numScenNode; w++) {
					 (*xOpt)[w]=sol[w];
				 }			 
			 }
		 }
		 else {
			 if (z_hat > FuncVal) {
				 FuncVal = z_hat;
				 *Worst_pOpt = pArray;
				 *SubObjOpt = SubObj;
				 for (w = 0; w < numScenNode; w++) {
					 (*xOpt)[w] = sol[w];
				 }
			 }
		 }
		 rel_Gap = abs(FuncVal - Appx) / abs(Appx);
		 Result << nCut + 1 << '\t' << z_hat << '\t' << Appx << '\t' << FuncVal << '\t' << rel_Gap * 100 << endl;
		 if (rel_Gap > TOLER  && float(clock() - begin_time) / CLOCKS_PER_SEC <= TIME_LIMIT) 
			 BackwardPass(&sol, &rng_pi, &Gs, &gs, &ps, nCut, problem_type, cut_type);
		 
		 nCut++;
	 }
	 gs.~vector();
	 Gs.~vector();
	 ps.~vector();
	 *objVal = Appx;
	 Result << "time = " << float(clock() - begin_time) / CLOCKS_PER_SEC << endl;
	 Result << "cost = " << *SubObjOpt<< endl;
	 Result << endl;
	 Result << "worst prob = "<< *Worst_pOpt<<endl;
	 Result << endl;
	 Result << "x[0] = " << (*xOpt)[0] << endl;
	 Result << endl;

	 rng_pi.clear();
	 rng_pi.end();
	 
	 
	 Result.close();

	 cout << "*****DONE with Original Problem*****" << endl;

 }//end SolveOriginal

