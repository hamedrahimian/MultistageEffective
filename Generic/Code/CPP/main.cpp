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

#include "effectivenessTV.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <io.h>
#include <iostream>


int main(int argc, char *argv[])
{
	IloEnv env2;
	try {

		if (argc < 4){
			cerr << "WRONG!  Usage: testlib <core file> <time file> <sto file>" << endl;
			exit(1);
		}

		if (argc != 6) {
			cout << "WRONG!  Usage: testlib <core file> <time file> <sto file> <problem_type> <cut_type>" << endl;
		}


		ProblemType problem_type;
		if (strncmp(argv[4], "DRSO", 4) == 0)
			problem_type = DRSO;
		else
			problem_type = RISK_NEUTRAL;

		CutType cut_type;
		if (strncmp(argv[5], "MULTI", 5) == 0)
			cut_type = MULTI;
		else
			cut_type = SINGLE;

		if (argv[4] == NULL && argv[5] == NULL) {
			problem_type = DRSO;
			cut_type = MULTI;
		}

		string problem_name=argv[3];
		//remove "data/"
		problem_name.erase(problem_name.begin(), problem_name.begin()+5);
		//remove ".sto"
		problem_name.erase(problem_name.end()-4, problem_name.end());


		EffectivenessTV prob;
		prob.readFiles(argv[1], argv[2], argv[3]);
		prob.printStochSummary();
		prob.setScenarioMode(EXPLICIT);
		prob.createScenarioTree();

		prob.writeDeterministicEquivalent("def.mps");

		//read differenet values of rho
		int numInstance;
		vector<double> config;
		if (strncmp(argv[4], "DRSO", 4) == 0) {
			string line;
			ifstream riskpar("Data\\config.txt");
			vector<double> vect;
			int kk = 0;
			if (riskpar.is_open()) {
				while (getline(riskpar, line)) {
					kk++;
					stringstream ss(line);
					double i;
					while (ss >> i) {
						vect.push_back(i);
						if (ss.peek() == ' ')
							ss.ignore();
					}
				}
			}
			riskpar.close();
			numInstance = kk;
			for (int i = 0; i < vect.size(); i++) {
				config.push_back(vect[i]);
			}
		}
		
		
		char *problem_c = new char[problem_name.size() + 1];
		problem_c[problem_name.size()] = 0;
		memcpy(problem_c, problem_name.c_str(), problem_name.size());


		/////////////////// OUTPUT /////////////////////////////
		char resName[100];
		const char* Resfilename = resName;
		sprintf_s(resName, "EX_Summary_%s.txt", problem_c);
		Resfilename = resName;
		ofstream Summary(Resfilename);
		Summary.precision(10);
		Summary << "Number of Scenario Paths.." << endl;
		Summary << '\t';
		Summary << "I" << '\t' << "E" << '\t' << "U" << '\t' << "I" << '\t' << "E" << endl;


		const int numStage = prob.getNumPeriods();

		int period;
		//print base lp models
		for (period = 0; period < numStage; period++)
			prob.PrintBaseModelToFile(problem_c, period);

		//print scenarios
		prob.PrintScenariosToFile(problem_c);

		if (problem_type != DRSO)
			numInstance = 1;
		double * rho = new double[numInstance];


		if (problem_type != DRSO)
			rho[0] = 0;

		for (int ii = 0; ii < numInstance; ii++) {
			if (problem_type == DRSO)
				rho[ii] = config[ii];

			IloEnv env;

			/////////////////// OUTPUT for CURRENT INSTANCE/////////////////////////////
			sprintf_s(resName, "EX_Effectiveness_%s_%0.1f.txt", problem_c, config[ii]);
			Resfilename = resName;
			ofstream Effectiveness(Resfilename);
			Effectiveness.precision(10);

			/////////////////// ORIGINAL PROBLEM /////////////////////////////
			IloNum objVal;
			IloNumArray2 x(env);
			IloNumArray SubObj(env);
			IloNumArray Worst_p(env);
			prob.SolveOriginal(problem_c, rho[ii], &objVal, &x, &SubObj, &Worst_p, problem_type, cut_type);

			/////////////////// EFFECTIVENESS RESULTS /////////////////////////////

			

			int s, w, j;

			IloIntArray2 num_cond_effect_category(env);
			IloIntArray num_path_effect_category(env);

			const int numParents = prob.numTotalNodesInstance(numStage - 2);

			OutputType * Output = new OutputType[numParents];

			vector<int>* IneffCondscen = new vector<int>[numParents];
			vector<int>* EffCondscen = new vector<int>[numParents];
			vector<int>* UnknownCondscen = new vector<int>[numParents];
			vector<int>* IneffPathscen = new vector<int>[numParents];
			vector<int>* EffPathscen = new vector<int>[numParents];
			vector<int>* UnknownPathscen = new vector<int>[numParents];

					
			prob.ETCAllCondEff(Effectiveness, config[ii], &SubObj, &Worst_p, &Output, &IneffCondscen, &EffCondscen, &UnknownCondscen, &num_cond_effect_category);
			if (numStage > 2)
				prob.ETCAllPathEff(Summary, Effectiveness, config[ii], &Output, &IneffPathscen, &EffPathscen, &UnknownPathscen, &num_path_effect_category);
			

			//prob.ResolveUnknownCondEff(Effectiveness, config[ii], &x, &SubObj, &Output, &IneffCondscen, &EffCondscen, &UnknownCondscen, &num_cond_effect_category);
			
			if (numStage > 2) {
				//Do this for larger problems
					prob.ETCUnknownPathEff(Summary, Effectiveness, config[ii], &Output, &IneffPathscen, &EffPathscen, &UnknownPathscen, &num_path_effect_category);
				//Do this for smaller problems
				//prob.ResolveAllPathEff(Summary, Effectiveness, config[ii], &Output, &IneffPathscen, &EffPathscen, &UnknownPathscen, &num_path_effect_category, objVal);
			}
			

			delete[] Output;
			delete[] IneffCondscen;
			delete[] EffCondscen;
			delete[] UnknownCondscen;
			delete[] IneffPathscen;
			delete[] EffPathscen;
			delete[] UnknownPathscen;
			num_cond_effect_category.clear();
			num_cond_effect_category.end();
			num_path_effect_category.clear();
			num_path_effect_category.end();

			

			x.clear();
			SubObj.clear();
			SubObj.end();
			Worst_p.clear();
			Worst_p.end();

			env.end();

		}//end instance rho
		
		delete[] rho;


	}

	catch (IloException &e) {
		env2.out() << "ERROR: " << e << endl;
	}
	catch (...) {
		env2.out() << "Unknown exception" << endl;
	}

	env2.end();

	getchar();
	return 0;
}
