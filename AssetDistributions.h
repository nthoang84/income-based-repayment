/********************************************************************************
AssetDistributions.h

This header file simulates the loan repayment behavior of the borrowers under 
a pre-specified repayment plan.

Version 1.0.0, May 15, 2022

@author: Hoang Nguyen

********************************************************************************/
#ifndef ASSET_DISTRIBUTIONS_H
#define ASSET_DISTRIBUTIONS_H

#include <assert.h>
#include <math.h>

#include <fstream>
#include <random>
#include <vector>

#include "Calculator.h"
#include "RepaymentModel.h"
#include "StandardRepaymentModel.h"


/********************************************************************************
Purpose: This class simulates the loan repayment behavior of the borrowers.

Parameters:
	model 				(a pointer to the RepaymentModel object already created in
						 TotalUtility.h)
	beta				(discount factor for utility tradeoff across periods)
	sigma				(standard deviation of annual income shock)
	rho					(AR(1) coefficient for income shock decay)
	gamma				(curvature parameter on consumption)
	commonNumPers		(number of periods used during simulation across all 
						 repayment models)
	initPermInc 		(the initial permanent income)
	incGrowthRate		(annual income growth rate)
	numReps 			(number of replicates in the repayment model)
	calculator			(an instance of Calculator class to provide functions for 
						 common calculation)
	permIncs 			(an 1D array of permanent incomes over commonNumPers periods)
	incomes 			(a 2D array of AR(1) income draws over commonNumPers periods, 
	  					 for num_reps replicates)
	consumptions 		(a 2D array of consumptions over commonNumPers periods, for 
						 num_reps replicates)
	u 					(a 2D array of utility of consumption over commonNumPers 
						 periods, for num_reps replicates)
	totalUtility 		(total expected utility discounted to the beginning period)

********************************************************************************/
class AssetDistributions {
public:
	RepaymentModel* model;
	double beta;
	double sigma;
	double rho;
	double gamma;
	int commonNumPers;
	double initPermInc;
	double incGrowthRate;
	
	int numReps;
	Calculator calculator;
	vector<double> permIncs;
	vector<vector<double>> incomes;
	vector<vector<double>> consumptions;
	vector<vector<double>> u;
	double totalUtility;


	/********************************************************************************
	Purpose: Setters for some instance variables.

	********************************************************************************/
	void setModel(RepaymentModel* model) { this->model = model; }

	void setBeta(double beta) { this->beta = beta; }

	void setSigma(double sigma) { this->sigma = sigma; }

	void setGamma(double gamma) { this->gamma = gamma; }

	void setRho(double rho) { this->rho = rho; }

	void setCommonNumPeriods(int commonNumPers) { this->commonNumPers = commonNumPers; }

	void setInitialPermanentIncome(double initPermInc) { this->initPermInc = initPermInc; }

	void setIncomeGrowthRate(double incGrowthRate) { this->incGrowthRate = incGrowthRate; }


	/********************************************************************************
	Purpose: Calculate the permanent incomes in every period.

	********************************************************************************/
	void calcPermanentIncomes() {

		permIncs.assign(commonNumPers, 0);
		for (int i = 0; i < commonNumPers; i++) {
			permIncs[i] = (i == 0 ? initPermInc : 
						   permIncs[i - 1] * (1 + incGrowthRate));
		}
	}


	/********************************************************************************
	Purpose: Calculate the epsilons in every period, for each replicate.

	Parameters:
		epsilons 	(a matrix of epsilons for each period and replicate)

	********************************************************************************/
	void calcEpsilons(vector<vector<double>>& epsilons) {

		epsilons.assign(commonNumPers, vector<double>(numReps, 0.0));
		if (sigma == 0.0) return;

		mt19937 generator(1234567890);
		normal_distribution<double> norm(-0.5 * sigma * sigma, sigma);
		
		for (int i = 0; i < commonNumPers; i++) {
			for (int j = 0; j < numReps; j++) {
				epsilons[i][j] = norm(generator);
			}
		}

		for (int i = 0; i < commonNumPers; i++) {
			double avgEpsilon = 0.0;
			for (int j = 0; j < numReps; j++) {
				avgEpsilon += exp(epsilons[i][j]);
			}
			avgEpsilon = avgEpsilon / numReps;
			assert(avgEpsilon > 0);
			for (int j = 0; j < numReps; j++) {
				epsilons[i][j] -= log(avgEpsilon);
			}
		}
	}


	/********************************************************************************
	Purpose: Calculate the unobserved in every period, for each replicate.

	Parameters:
		unobserved 	(a matrix of unobserved for each period and replicate)
		epsilons 	(a matrix of epsilons for each period and replicate)

	********************************************************************************/
	void calcUnobserved(vector<vector<double>>& unobserved, 
						const vector<vector<double>>& epsilons) {

		unobserved.assign(commonNumPers, vector<double>(numReps, 0));
		for (int i = 0; i < commonNumPers; i++) {
			for (int j = 0; j < numReps; j++) {
				unobserved[i][j] = (i == 0 ? 0 : rho * unobserved[i - 1][j]) + 
								   epsilons[i][j];
			}
		}
	}


	/********************************************************************************
	Purpose: Calculate the real incomes in every period, for each replicate.

	Parameters:
		unobserved 	(a matrix of unobserved for each period and replicate)

	********************************************************************************/
	void calcIncomes(const vector<vector<double>>& unobserved) {

		incomes.assign(commonNumPers, vector<double>(numReps, 0));
		for (int i = 0; i < commonNumPers; i++) {
			for (int j = 0; j < numReps; j++) {
				incomes[i][j] = permIncs[i] * exp(unobserved[i][j]);
				assert(incomes[i][j] > 0.0);
			}
		}
	}


	/********************************************************************************
	Purpose: Simulate the model and calculate the consumption amount in every period,
		for each replicate.

	********************************************************************************/
	void calcConsumptions() {

		consumptions.assign(commonNumPers, vector<double>(numReps, 0));
		for (int i = 0; i < commonNumPers; i++) {
			for (int j = 0; j < numReps; j++) {
				double payment, underpaid;
				tie(payment, underpaid) = model->calcPayment(incomes[i][j]);

				consumptions[i][j] = (i < model->numPers ? 
					model->pay(i, j, incomes[i][j], payment, underpaid) : 
					incomes[i][j]);

				assert(consumptions[i][j] > 0);
			}
		}
	}


	/********************************************************************************
	Purpose: Calculate the expected utility in each period.

	********************************************************************************/
	void calcUtilities() {

		u.assign(commonNumPers, vector<double>(numReps, 0));
		for (int i = 0; i < commonNumPers; i++) {
			for (int j = 0; j < numReps; j++) {
				u[i][j] = calculator.u(consumptions[i][j]);
			}
		}
	}


	/********************************************************************************
	Purpose: Calculate the total expected utility, discounted to the present.

	********************************************************************************/
	double calcTotalUtility() {

		totalUtility = 0.0;
		for (int i = 0; i < commonNumPers; i++) {
			double avgUtil = accumulate(u[i].begin(), u[i].end(), 0.0) / u[i].size();
			totalUtility += avgUtil * pow(beta, i);
		}
		return totalUtility;
	}


	/********************************************************************************
	Purpose: Calculate the total amount forgiven (already discounted to the present).

	********************************************************************************/
	double calcAverageForgiven() {

		double totalForgiven = accumulate(model->forgiven.begin(),
										  model->forgiven.end(), 0.0);
		return totalForgiven / numReps;
	}


	/********************************************************************************
	Purpose: Initialize instance variables then do all calculations of interest.

	********************************************************************************/
	void setup() {

		calculator.setGamma(gamma);
		numReps = model->numReps;
		calcPermanentIncomes();
		
		vector<vector<double>> epsilons, unobserved;
		calcEpsilons(epsilons);
		calcUnobserved(unobserved, epsilons);
		calcIncomes(unobserved);

		calcConsumptions();
		calcUtilities();
		calcTotalUtility();
	}


	/********************************************************************************
	Purpose: Print the summary statistics of the repayment behavior of the borrowers.

	Parameters:
		dataset 	(the name of the dataset)

	********************************************************************************/
	void summarize(string dataset) {

		string modelType = model->name;
		double pctRepaid = 
			(double) accumulate(model->paidOff.begin(), model->paidOff.end(), 0) /
			numReps;
		double totalForgiven = accumulate(model->forgiven.begin(), 
										  model->forgiven.end(), 0);
		double avgForgiven = totalForgiven / numReps;
		double condAvgForgiven = (pctRepaid == 1.0 ? 0 : 
								  avgForgiven / (1 - pctRepaid));
		double govRevenue = accumulate(model->govRevenue.begin(), 
									   model->govRevenue.end(), 0);

		string strResults = dataset + "," + 
							modelType + "," +
							to_string(pctRepaid) + "," + 
							to_string(totalForgiven) + "," +
							to_string(avgForgiven) + "," + 
							to_string(condAvgForgiven) + "," +
							to_string(govRevenue) + "\n";

		string outputFile = "Results/Summaries/summary_" + dataset + "_" + modelType + ".csv";
	    ofstream ofs;
	    ofs.open(outputFile);
	    ofs << "dataset,modelType,pctRepaid,totalForgiven,avgForgiven,condAvgForgiven,govRevenue\n";
	    ofs << strResults;
	    ofs.close();
	}


	/********************************************************************************
	Purpose: Print the detailed repayment behavior of the borrowers.

	Parameters:
		dataset 	(the name of the dataset)

	********************************************************************************/
	void printSchedule(string dataset) {

		string modelType = model->name;
		string outputFile = "Results/Schedules/schedule_" + dataset + "_" + modelType + ".csv";
	    ofstream ofs;
	    ofs.open(outputFile);
	    
	    string header;
	    for (int i = 0; i < commonNumPers; i++) {
	    	header += ("Q" + to_string(i) + ",");
	    	header += ("Y" + to_string(i) + ",");
	    	header += ("q" + to_string(i) + ",");
	    	header += ("C" + to_string(i) + ",");
	    }
	    header += "paidOff,forgiven,penalty\n";
	    ofs << header;

	    for (int j = 0; j < numReps; j++) {
	    	string info;
	    	for (int i = 0; i < commonNumPers; i++) {
	    		info += (i < model->numPers ? to_string(model->balances[i][j]) : "0") + ",";
	    		info += to_string(incomes[i][j]) + ",";
	    		info += (i < model->numPers ? to_string(model->payments[i][j]) : "0") + ",";
	    		info += to_string(consumptions[i][j]) + ",";
	    	}

	    	info += (model->paidOff[j] ? "1," : "0,");
		    info += to_string(model->forgiven[j]) + ",";
		    info += to_string(model->govRevenue[j]) + "\n";
		    ofs << info;
	    }

	    ofs.close();
	}
};

#endif