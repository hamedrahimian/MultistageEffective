/*
*     Nested L-shaped method for DRSO and Risk-Averse Optimization-- Built on top of SUTIL
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


#ifndef EFFECTIVENESSTV_H
#define EFFECTIVENESSTV_H

#include "effectiveness.h"


class EffectivenessTV : public Effectiveness {
	public:

		/** Forms primal categories*/
		//cIx is the overall index
		//fill categories by overall indices
		void FormPrimalCategories(const int numChild, const double * cost, const int VaR_index, const double VaR_value,
			vector<int>** categories, cost_category  *primal_category, IloIntArray* num_cost_category);

		/** Easy-to-sheck condition for the effectiveness of scenario #test_scenario with q>0 in VaR category, where |VaR|>1*/
		bool ETCCondEffVaR(const int numChild, const int test_scenario, const double * sorted_cost, const double * sorted_nominal_prob,
			const int old_VaR_index, const double old_VaR_value, const double old_beta);

		/** Checks conditional effectiveness of ALL scenario nodes with easy-to-check conditions*/
		void ETCAllCondEff(ostream& fe, const double rho, const IloNumArray* SubObj, const IloNumArray* Worst_p, OutputType** Output, vector<int>** IneffCondScen, vector<int>** EffCondScen, vector<int>** UnknownCondScen,
			IloIntArray2* num_cond_effect_category);
		
		/** Checks effectiveness of ALL scenarios paths with easy-to-check conditions, before reslving scenario nodes to identify conditional effectiveness of unknown scenario nodes*/
		void ETCAllPathEff(ostream& fs, ostream& fe, const double rho, OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen,
			IloIntArray* num_path_effect_category);

		/** Checks effectiveness of UNKNOWN scenarios paths with easy-to-check conditions, after reslving scenario nodes to identify conditional effectiveness of unknown scenario nodes*/
		void ETCUnknownPathEff(ostream& fs, ostream& fe, const double rho, OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen,
			IloIntArray* num_path_effect_category);

		/** Resolve problem to find the conditional effectiveness ONLY UNKNOWN scenario nodes*/
		void ResolveUnknownCondEff(ostream& fe, const double rho, const IloNumArray2* x_hat, const IloNumArray* optSubObj, OutputType** Output, vector<int>** IneffCondScen, vector<int>** EffCondScen,
			vector<int>** UnknownCondScen, IloIntArray2* num_cond_effect_category);

		/** Resolve problem to find the conditional effectiveness of ALL scenario nodes*/
		void ResolveAllCondEff(ostream& fe, const double rho, const IloNumArray2* x_hat, const IloNumArray* optSubObj, OutputType** Output, vector<int>** IneffCondScen, vector<int>** EffCondScen,
			vector<int>** UnknownCondScen, IloIntArray2* num_cond_effect_category);

		/** Resolve problem to find the effectiveness of ONLY UNKONWN scenario paths*/
		void ResolveUnknownPathEff(ostream& fs, ostream& fe, const double rho,
			OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen, IloIntArray* num_path_effect_category, const double objVal);

		/** Resolve problem to find the effectiveness of ALL scenario paths*/
		void ResolveAllPathEff(ostream& fs, ostream& fe, const double rho,
			OutputType** Output, vector<int>** IneffPathScen, vector<int>** EffPathScen, vector<int>** UnknownPathScen, IloIntArray* num_path_effect_category, const double objVal);

		
};


#endif