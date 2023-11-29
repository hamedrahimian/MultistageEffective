Decomposition algorithms and scenario generation codes in C++
---------------------------------------------------------------

1. Decomposition algorithm. This code implements the dual algorithm D5 (all variants) and the nested primal decomposition algorithm to solve mean-CVaR or sup-CVaR problem. 
It also implmenets easy-to-check conditions to identify the conditional effectiveness of realizations and effectiveness of scenario paths for sup-CVaR model ("TV").
This code includes source files "decomposition.cpp",  "effectiveness.cpp", "input.cpp", "tree.cpp", "vars.cpp", and "main.cpp". 
The header files are "decomposition.h",  "effectiveness.h", "input.h", "tree.h", "vars.h", and "types.h". 
	
	1.a. The inputs of the code are files "input1.txt", "input2.txt", "DemPortion.txt", "config_%s.txt",  where s is either "TV" (sup-CVaR) or "EC" (mean-CVaR), 
		a scenario file "scenFac_%d.txt", where d is the number of children per node (variable numScen should also be set to the number of children per node and must match %d in "scenFac_%d.txt")
		and a series of command line arguments (see below).

	1.b. The outputs are 1) the log of the iterations, the CPU time, optimal first stage decision, optimal worst-case probabilities, and optimal costs, and 2) the information about the conditional effectiveness and effectiveness of scenario paths. 



2. "ScenGen.cpp" to generate scenario files 
	1.a. The inputs are variables numStage and numScen. By default they are 4 and 25, respectively.  
	1.b. The output is "scenFac_%d.txt", where d is the number of children per node.
	



Command Line Arguments: 
-------------------------
The user should enter the following arguments (in order) to run the code:

1. The type of the ambiguity set: 1) TV, 2) EC.
2. The type of the problem: 1) PRIMAL (this means DRSO problem), 2) DUAL (this means risk-averse problem).
3. The type of the cut on sup or mean: 1) SINGLE, 2) MULTI.
4. The type of the cut on CVaR: 1) SINGLE, 2) MULTI.
5. The type of approximations separation: 1) COMBINED (meaning that mean-CVaR or sup-CVaR is approximated combined), 2)SEPARTAED (meaning that mean and CVaR, or sup and CVaR are approximated separately).

For example, a list of arguments TV PRIMAL MULTI MULTI COMBINED implements the nested primal decomposition algorithm proposed in the report. Note that, when argument 5 is COMBINED, argument 4 is meaningless, but it must be entered for consistency of the code. 
As another example, a list of arguments EC DUAL MULTI SINGLE SEPARATED implements algorithm MS_D5S. 