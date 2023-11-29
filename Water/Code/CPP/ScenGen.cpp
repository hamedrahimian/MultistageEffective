#include <ilcplex/ilocplex.h>
#include <stdlib.h>   
#include <vector>
#include <time.h>
//#include <random>
#include <chrono>
#include <cmath>
using namespace std;

ILOSTLBEGIN


int main(int argc, char **argv) {
srand ((unsigned)time(NULL));
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
 IloEnv env;

 try {
	 IloInt s, t, w;
	 IloInt numYearMax = 15 ; 
	 IloInt numStage = 4;
	 IloInt numScen =25; 
	 IloNum prob = (double)1/numScen;
	 IloInt numScenElm = numYearMax+2;
	 IloNumArray SupplyThreshold(env, numStage, 0.0, 0.1, 0.25, 0.35);

	 IloNumArray numScenNode_stage_sum(env, numStage);
	 numScenNode_stage_sum[0] = 1;
	 for (s = 1; s < numStage; s++) {
		 numScenNode_stage_sum[s]=numScenNode_stage_sum[s - 1] + pow(numScen, s);
	 }

	 IloIntArray stageMap(env, numScenNode_stage_sum[numStage - 1]);
	 w = 0;
	 for (s = 0; s < numStage; s++) {
		 while (w < numScenNode_stage_sum[s]) {
			 stageMap[w] = s;
			 w++;
		 }
	 }
	 IloNumArray2  scenFac(env, numScenNode_stage_sum[numStage-1]);
	
	 for(w=0; w<numScenNode_stage_sum[numStage - 1]; w++) {
		 scenFac[w] = IloNumArray(env, numScenElm);
		 for(t = 0;t<numScenElm;t++){
			 scenFac[w][t] = 1;
		 }
	 }


	 w = 1;
	 for (s = 1; s < numStage; s++) {
		 while (w < numScenNode_stage_sum[s]) {
			 double randomNo = (double)rand() / (double)RAND_MAX;
			 if (randomNo <= SupplyThreshold[s]) {
				 scenFac[w][0] = 0.9;
			 }
			 else {
				 scenFac[w][0] = 1;
			 }

			 for (t = 1; t<numScenElm; t++) {
				 double randomNo2 = 0.9 + 0.2*(double)rand() / (double)RAND_MAX;
				 scenFac[w][t] = randomNo2;
			 }

			 w++;
		 }
		
	 }
   	
	 char resName[100];
	 const char* Scenfilename = resName;
	 sprintf_s(resName, "scenFac_%d.txt", numScen);
	 Scenfilename = resName;
	 ofstream scenTree(Scenfilename);
	 scenTree<<'[';
	 for(w=0; w<numScenNode_stage_sum[numStage-1]; ++w) {
		 scenTree<<'[';
		 for(t = 0;t<numScenElm;++t){
			 if(t+1< numScenElm)
				 scenTree<<scenFac[w][t] <<',' ;
			 else
				 scenTree<<scenFac[w][t]; 
		 }
		 if (w+1 < numScenNode_stage_sum[numStage - 1]){
			 scenTree<<']'<< ','<< endl;
		 }
		 else{
			 scenTree<<']' << ']' << endl;
		 }
	 }
	 scenTree.close();
	 cout<<"DONE";

	}

    catch(IloException &e) {
		env.out() << "ERROR: " << e << endl;
	}
	catch(...){
		env.out() << "Unknown exception" << endl;
	}
	env.end();
	
	getchar();
	return 0;

}