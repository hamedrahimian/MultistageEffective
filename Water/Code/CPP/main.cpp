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
#include "input.h"
#include "tree.h"
#include "decomposition.h"
#include "effectiveness.h"


int main(int argc, char **argv) {
	IloEnv env;
	try {

		IloInt w;
		int numInstance;
		vector < pair<double, double> > config;

		if (argc != 6) {
			cout << "WRONG!  Usage: testlib <ambiguity_type> <problem_type> <cut1_type> <cut2_type> <separation_type>" << endl;
		}
		setTypes(argv);
		types = getVariantType();
		char* variant_name = setFilename();

		ReadCoeff(argv[1], numInstance, config);


		///////////////// READ DATA ////////////////////////////////
		ReadData(numScen, argv[0]);

		getnumTheta();

		prob = IloNumArray(env, numScenNode);


		prob[0] = 1;
		for (w = 1; w < numScenNode; w++)
			prob[w] = (double)1 / numScen;


		///////////////// DETERMINE THE STRUCTURE OF THE TREE /////////////////
		TreeStructure(ancestor, descendant, stageMap);

		///////////////// STORE INEFFECTIVE, EFFECTIVE, and UNKNOWN SCENARIOS /////////////////
		
		//vector to stores each gamma ineffective scenarios
		vector<int> * IneffCondscen = new vector<int>[21];
		//vector to stores each gamma effective scenarios
		vector<int> * EffCondscen = new vector<int>[21];
		//vector to stores each gamma undetermined scenarios
		vector<int> * UnknownCondscen = new vector<int>[21];

		//vector to stores each gamma ineffective scenarios
		vector<int> * IneffPathscen = new vector<int>[21];
		//vector to stores each gamma effective scenarios
		vector<int> * EffPathscen = new vector<int>[21];
		//vector to stores each gamma undetermined scenarios
		vector<int> * UnknownPathscen = new vector<int>[21];
		

		char resName[100];
		const char* Resfilename = resName;
		sprintf_s(resName, "EX_%d_Nestedness.txt", numScen);
		Resfilename = resName;
		ofstream Nestedness(Resfilename);
		Nestedness.precision(10);

		sprintf_s(resName, "EX_%d_Summary.txt", numScen);
		Resfilename = resName;
		ofstream Summary(Resfilename);
		Summary.precision(10);
		Summary << "Number of Scenario Paths.." << endl;
		Summary << '\t';
		Summary << "I" << '\t' << "E" << '\t' << "U" << '\t' << "I" << '\t' << "E" << endl;

		for (int ii = 0; ii < numInstance; ii++) {
			CVaRno = CalculateCVaRno(config[ii].second);

			/////////////////// MODELS /////////////////////////////
			Formulation formulation;
			ModelArray MODEL(env);
			CplexArray CPX(env);


			createMaster(MODEL, CPX, formulation.constraints, formulation.mainvariables, formulation.appxvariables, config[ii]);

			DistGenFormulation distgen;
			IloCplex CplexSub(env);
			IloModel modSub(env);

			if (variant.problem_type == PRIMAL)
				createSub(CplexSub, modSub, distgen, config[ii]);

			///////////////// PARAMETERS ///////////////// 
			IloNum objVal;
			IloIntArray2 num_cost_category(env);
			IloIntArray2 num_cond_effect_category(env);
			IloIntArray num_path_effect_category(env);

			Scenario * Output = new Scenario[numScenNode];
			for (w = 0; w < numScenNode; w++) {
				Output[w].partial_sum = 0;
				Output[w].No = w;
				Output[w].nom_Prob = prob[w];
				Output[w].CondStatus = effect_category::unknown;
				Output[w].PathStatus = effect_category::unknown;
				Output[w].primal_category = cost_category::a_VaR;
				Output[w].CondObjCost = -1;
				Output[w].PathObjCost = -1;
			}


			char resName[100];
			sprintf_s(resName, "EX_%d_%s_%0.1f_%0.1f.txt", numScen, variant_name, config[ii].first, config[ii].second);
			Resfilename = resName;
			ofstream Result(Resfilename);
			Result.precision(10);

			sprintf_s(resName, "EX_%d_%0.2f_Effectiveness.txt", numScen, config[ii].first);
			Resfilename = resName;
			ofstream Effectiveness(Resfilename);
			Effectiveness.precision(15);


			Initialization(env, formulation.cuts, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.optsols, formulation.costsix);

			SolveOriginal(Result, MODEL, CPX, formulation, CplexSub, modSub, distgen, objVal, config[ii], Output);

			

			//Termination(MODEL, formulation.cuts, formulation.duals, formulation.cutcoeffs, formulation.incsols, formulation.optsols);

			if (variant.ambiguity_type == TV && variant.problem_type == PRIMAL) {
				EasyCheckCondEff(Effectiveness, Output, IneffCondscen[ii], EffCondscen[ii], UnknownCondscen[ii], config[ii], num_cost_category, num_cond_effect_category);

				if (numStage > 2)
					EasyCheckPathEff(Summary, Effectiveness, Output, IneffPathscen[ii], EffPathscen[ii], UnknownPathscen[ii], config[ii], num_path_effect_category);

				/*ResolveCondEff(Effectiveness, MODEL, CPX, formulation, CplexSub, modSub, distgen, config[ii], Output,
					IneffCondscen[ii], EffCondscen[ii], UnknownCondscen[ii], num_cost_category, num_cond_effect_category);

				if (numStage > 2)
					ResolvePathEff(Summary, Effectiveness, MODEL, CPX, formulation, CplexSub, modSub, distgen, objVal, config[ii], Output, IneffPathscen[ii], EffPathscen[ii], UnknownPathscen[ii], num_path_effect_category); */

			}

			delete[] Output;
			MODEL.end();
			CPX.end();
			modSub.end();
			CplexSub.end();
		}

		//Nested(Nestedness, EffPathscen);
	}

	catch (IloException &e) {
		env.out() << "ERROR: " << e << endl;
	}
	catch (...) {
		env.out() << "Unknown exception" << endl;
	}

	env.end();

	getchar();
	return 0;
}

