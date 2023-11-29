Decomposition algorithm. This code implements the nested primal decomposition algorithm to solve sup-CVaR model for any general problem, given by smps files.
It also implmenets easy-to-check conditions to identify the conditional effectiveness of realizations and effectiveness of scenario paths.
This code includes source files "decomposition.cpp",  "effectiveness.cpp", "effectivenessTV.cpp", "sparseMatrix.cpp", "scenariotree.cpp", and "main.cpp". 
The header files are "decomposition.h",  "effectiveness.h", "effectivenessTV.h", "sparseMatrix.h", and "scenariotree.h". 
	
	1.a. The inputs of the code are smsp files .cor, .tim, and .sto. 
		and a series of command line arguments (see below).

	1.b. The outputs are 1) the log of the iterations, the CPU time, optimal first stage decision, optimal worst-case probabilities, and optimal costs, and 2) the information about the conditional effectiveness and effectiveness of scenario paths. 


Command Line Arguments: 
-------------------------
The user should enter the following arguments (in order) to run the code:

1. .cor file. 
2. .tim file.
3. .sto file. 
4. The type of the problem: 1) DRSO, 2) RISK_NEUTRAL.
5. The type of the cut: 1) SINGLE, 2) MULTI.


A list of arguments Data/sgpf3y3.cor Data/sgpf3y3.tim Data/sgpf3y3.sto RISK_NEUTRAL MULTI, reads 
sgpf3y3 problem smps files from folder "Data", and solves the Risk Neutral problem with multi cuts. 

Arguments 4 and 5 are optional. Default is DRSO MULTI. There is no DRSO SINGLE option. 