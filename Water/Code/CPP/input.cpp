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

#include "input.h"
#include "vars.h"

using namespace std;

void ReadCoeff(const char ambiguitytype[], int& numInstance, vector < pair<double, double> >& config) {
	//read differenet values of rho, lambda, alpha

	vector<double> vect;
	string line;
	char configName[100];
	const char* configfilename = configName;
	sprintf_s(configName, "Data\\config_%s.txt", ambiguitytype);
	configfilename = configName;
	ifstream riskpar(configfilename);
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


	if (variant.ambiguity_type == TV) {
		pair<double, double> temp;
		for (int i = 0; i < vect.size(); i++) {
			temp.first = vect[i]; //rho
			temp.second = vect[i];
			config.push_back(temp);
		}
	}
	else if (variant.ambiguity_type == EC) {
		kk = 0;
		pair<double, double> temp;
		for (int i = 0; i < vect.size(); i++) {
			if ((i + 1) % 2 != 0) {
				//read lambda
				temp.first = vect[i];
			}
			else {
				//read alpha
				temp.second = vect[i];
				config.push_back(temp);
			}
		}
	}

}//ReadCoeff

void ReadData(int& numChild, const char progname[])
{
	IloInt i, s, t;

	IloEnv env;
	IloIntArray type(env);
	const IloInt numType = 10, numNode = 62;

	IloIntArray startNodeTemp(env), endNodeTemp(env), numYearTemp(env);
	IloNumArray   lossTemp(env), costTemp(env), capacityTemp(env), Storage0Temp(env), storageUBTemp(env);


	const char* fileData0 = "Data\\input1.dat";
	ifstream Data0(fileData0);
	if (!Data0) {
		cerr << "ERROR: Cannot Open File '" << fileData0 << "' for Reading" << endl;
		cerr << "usage:   " << progname << " <file>" << endl;
		throw(-1);
	}
	Data0 >> numYearTemp >> type >> startNodeTemp >> endNodeTemp >> costTemp >> capacityTemp >> returnRate >> lossTemp >> Storage0Temp >> storageUBTemp;

	numYear = numYearTemp; startNode = startNodeTemp; endNode = endNodeTemp; cost = costTemp;
	capacity = capacityTemp; loss = lossTemp; Storage0 = Storage0Temp; storageUB = storageUBTemp;

	IloNumArray2 scenFacTemp(env);
	char fileName[100];
	const char* fileData1 = fileName;
	sprintf_s(fileName, "Data\\scenFac_%d.txt", numChild);
	ifstream Data1(fileData1);
	Data1 >> numStage >> numScen >> numScenNode >> scenFacTemp;
	Data1.close();
	scenFac = scenFacTemp;

	numYear_stage_sum.push_back(0);
	for (s = 1; s < numStage + 1; s++) {
		numYear_stage_sum.push_back(numYear[s - 1] + numYear_stage_sum[s - 1]);
	}

	for (s = 0; s < numStage; s++) {
		total_pi += pow(numScen, s)*numYear[s];
	}

	total_pi_stage_sum.push_back(0);
	for (s = 1; s < numStage + 1; s++) {
		total_pi_stage_sum.push_back(pow(numScen, s - 1)*numYear[s - 1] + total_pi_stage_sum[s - 1]);
	}


	numLink = startNode.getSize();

	IloBool consistentData = (numYear.getSize() == numStage && endNode.getSize() == numLink
		&& cost.getSize() == numLink && type.getSize() == numNode);

	if (!consistentData) {
		cerr << "ERROR: Inconsistent Data1!" << endl;
		throw(-1);
	}

	cout << "Read Data Parts 0 and 1 Successfully" << endl;

	//identify role of each node in the water network, and add node number to their corresponding role vector
	//it acts to know which constraints are needed for each node

	for (i = 0; i<numNode; i++) {
		if (type[i] == 4) {
			rechargeID.push_back(i);
		}
		else if (type[i] == 0 || type[i] == 2 || type[i] == 5 || type[i] == 6 || type[i] == 7) {
			balanceID.push_back(i);
		}
		else if (type[i] == 1) {
			capacityID.push_back(i);
		}
		else if (type[i] == 8) {
			userID.push_back(i);
			potUserID.push_back(i);
		}
	}
	for (i = 0; i<numNode; i++) {
		if (type[i] == 9) {
			userID.push_back(i);
		}
		else if (type[i] == 4) {
			capacityID.push_back(i);
		}
	}
	numUser = int(userID.size());
	numPotUser = int(potUserID.size());
	numBalance = int(balanceID.size());
	numCapacity = int(capacityID.size());
	numRecharge = int(rechargeID.size());

	IloNumArray2   population(env, numStage), TucsonPop(env, numStage);
	IloNumArray    demUnit(env, numUser);
	for (s = 0; s< numStage; s++) {
		population[s] = IloNumArray(env, numYear[s]);
		TucsonPop[s] = IloNumArray(env, numYear[s]);
	}
	const char* fileData2 = "Data\\input2.dat";
	ifstream Data2(fileData2);
	if (!Data2) {
		cerr << "ERROR: Cannot Open File '" << fileData2 << "' for Reading" << endl;
		cerr << "usage:   " << progname << " <file>" << endl;
		throw(-1);
	}

	Data2 >> population[0] >> population[1] >> population[2] >> population[3] >> demUnit >> TucsonPop[0] >> TucsonPop[1] >> TucsonPop[2] >> TucsonPop[3];
	Data2.close();
	cout << "Read Data Part 2 Successfully" << endl;

	const char* fileData3 = "Data\\DemPortion.txt";
	ifstream Data3(fileData3);
	if (!Data3) {
		cerr << "ERROR: Cannot Open File '" << fileData3 << "' for Reading" << endl;
		cerr << "usage:   " << progname << " <file>" << endl;
		throw(-1);
	}
	IloNumArray2 DembyNode(env, numYear_stage_sum[numStage]);
	for (t = 0; t<numYear_stage_sum[numStage]; t++) {
		DembyNode[t] = IloNumArray(env, numPotUser);
	}
	for (t = 0; t<numYear_stage_sum[numStage]; t++) {
		Data3 >> DembyNode[t];
	}
	Data3.close();
	cout << "Read Data Part 3 Successfully" << endl;

	IloNumArray demandTemp(env, numUser*numYear_stage_sum[numStage]);
	for (s = 0; s < numStage; s++) {
		for (t = 0; t < numYear[s]; t++) {
			for (i = 0; i < numPotUser; i++) {
				demandTemp[numYear_stage_sum[s] * numUser + t*numUser + i] = double(135 * 0.00112*population[s][t] * DembyNode[t + numYear_stage_sum[s]][i] * 0.8);
			}
		}
	}

	for (s = 0; s<numStage; s++) {
		for (t = 0; t<numYear[s]; t++) {
			for (i = numPotUser; i<numUser; i++) {
				demandTemp[numYear_stage_sum[s] * numUser + t*numUser + i] = demandTemp[numYear_stage_sum[s] * numUser + t*numUser + (i - numPotUser)] * 0.25;
			}
		}
	}
	demand = demandTemp;

	IloNumArray CAPamtTemp(env, numYear_stage_sum[numStage]);
	for (s = 0; s< numStage; s++) {
		for (t = 0; t<numYear[s]; t++) {
			CAPamtTemp[numYear_stage_sum[s] + t] = 144000 * population[s][t] / TucsonPop[s][t];
		}
	}
	CAPamt = CAPamtTemp;
}

void setTypes(char **argv) {
	if (strncmp(argv[1], "TV", 2) == 0)
		variant.ambiguity_type = TV;
	else
		variant.ambiguity_type = EC;

	if (strncmp(argv[2], "PRIMAL", 5) == 0)
		variant.problem_type = PRIMAL;
	else if (strncmp(argv[2], "DUAL", 4) == 0)
		variant.problem_type = DUAL;
	else
		variant.problem_type = RISK_NEUTRAL;

	if (strncmp(argv[3], "MULTI", 5) == 0)
		variant.cut1_type = MULTI;
	else
		variant.cut1_type = SINGLE;

	if (strncmp(argv[4], "MULTI", 5) == 0)
		variant.cut2_type = MULTI;
	else
		variant.cut2_type = SINGLE;

	if (strncmp(argv[5], "COMBINED", 8) == 0)
		variant.separation_type = COMBINED;
	else
		variant.separation_type = SEPARATED;

}//setTypes

VariantType getVariantType(void) {
	if (variant.ambiguity_type == EC && variant.problem_type == PRIMAL && variant.cut1_type == MULTI && variant.separation_type == COMBINED)
		types = VariantType::EC_P_M_C;
	else if (variant.ambiguity_type == EC && variant.problem_type == DUAL && variant.cut1_type == SINGLE && variant.separation_type == COMBINED)
		types = VariantType::EC_D_S_C;
	else if (variant.ambiguity_type == EC && variant.problem_type == DUAL && variant.cut1_type == MULTI && variant.separation_type == COMBINED)
		types = VariantType::EC_D_M_C;
	else if (variant.ambiguity_type == EC && variant.problem_type == DUAL && variant.cut1_type == SINGLE && variant.cut2_type == SINGLE && variant.separation_type == SEPARATED)
		types = VariantType::EC_D_S_S_S;
	else if (variant.ambiguity_type == EC && variant.problem_type == DUAL && variant.cut1_type == SINGLE && variant.cut2_type == MULTI && variant.separation_type == SEPARATED)
		types = VariantType::EC_D_S_M_S;
	else if (variant.ambiguity_type == EC && variant.problem_type == DUAL && variant.cut1_type == MULTI && variant.cut2_type == SINGLE && variant.separation_type == SEPARATED)
		types = VariantType::EC_D_M_S_S;
	else if (variant.ambiguity_type == EC && variant.problem_type == DUAL && variant.cut1_type == MULTI && variant.cut2_type == MULTI && variant.separation_type == SEPARATED)
		types = VariantType::EC_D_M_M_S;
	else if (variant.ambiguity_type == TV && variant.problem_type == PRIMAL && variant.cut1_type == MULTI && variant.separation_type == COMBINED)
		types = VariantType::TV_P_M_C;
	else if (variant.ambiguity_type == TV && variant.problem_type == DUAL && variant.cut1_type == SINGLE && variant.separation_type == COMBINED)
		types = VariantType::TV_D_S_C;
	else if (variant.ambiguity_type == TV && variant.problem_type == DUAL && variant.cut1_type == MULTI && variant.separation_type == COMBINED)
		types = VariantType::TV_D_M_C;
	else if (variant.ambiguity_type == TV && variant.problem_type == DUAL && variant.cut1_type == SINGLE && variant.cut2_type == SINGLE && variant.separation_type == SEPARATED)
		types = VariantType::TV_D_S_S_S;
	else if (variant.ambiguity_type == TV && variant.problem_type == DUAL && variant.cut1_type == SINGLE && variant.cut2_type == MULTI && variant.separation_type == SEPARATED)
		types = VariantType::TV_D_S_M_S;
	else if (variant.ambiguity_type == TV && variant.problem_type == DUAL && variant.cut1_type == MULTI && variant.cut2_type == SINGLE && variant.separation_type == SEPARATED)
		types = VariantType::TV_D_M_S_S;
	else if (variant.ambiguity_type == TV && variant.problem_type == DUAL && variant.cut1_type == MULTI && variant.cut2_type == MULTI && variant.separation_type == SEPARATED)
		types = VariantType::TV_D_M_M_S;

	return types;
}//getVariantType


char* setFilename(void) {
	char* configfilename = NULL;

	switch (types) {
	case  VariantType::EC_P_M_C:
		configfilename = "EC_P_M_C"; break;
	case  VariantType::EC_D_S_C:
		configfilename = "EC_D_S_C"; break;
	case  VariantType::EC_D_M_C:
		configfilename = "EC_D_M_C"; break;
	case  VariantType::EC_D_S_S_S:
		configfilename = "EC_D_S_S_S"; break;
	case  VariantType::EC_D_S_M_S:
		configfilename = "EC_D_S_M_S"; break;
	case  VariantType::EC_D_M_S_S:
		configfilename = "EC_D_M_S_S"; break;
	case  VariantType::EC_D_M_M_S:
		configfilename = "EC_D_M_M_S"; break;
	case  VariantType::TV_P_M_C:
		configfilename = "TV_P_M_C"; break;
	case  VariantType::TV_D_S_C:
		configfilename = "TV_D_S_C"; break;
	case  VariantType::TV_D_M_C:
		configfilename = "TV_D_M_C"; break;
	case  VariantType::TV_D_S_S_S:
		configfilename = "TV_D_S_S_S"; break;
	case  VariantType::TV_D_S_M_S:
		configfilename = "TV_D_S_M_S"; break;
	case  VariantType::TV_D_M_S_S:
		configfilename = "TV_D_M_S_S"; break;
	case  VariantType::TV_D_M_M_S:
		configfilename = "TV_D_M_M_S"; break;
	}
	return configfilename;
}//getnumTheta