/********************************************************************************
RepaymentModel.h

This header file provides a general structure for different repayment models.

Version 1.0.0, May 15, 2022

@author: Hoang Nguyen

********************************************************************************/
#ifndef REPAYMENT_MODEL_H
#define REPAYMENT_MODEL_H

#include <assert.h>

#include <vector>

using namespace std;


/********************************************************************************
Purpose: A general class of repayment model. Specific kinds of repayment model 
	would inherit variables and methods from this class.

Parameters:
	name 				(name of the repayment model, default = "gen")
	numReps 			(number of replicates)
	interestRate 		(annual interest rate of the loan)
	numPers 			(number of loan repayment periods)
	loan				(initial loan)
	graceAmount			(minimum possible amount for consumption in each period)
	balances			(a 2D array of remaining balances over num_pers periods, 
						 for numReps replicates)
	payments			(a 2D array of payments over num_pers periods, for 
						 numReps replicates)
	paidOff				(an 1D array of booleans indicating whether numReps 
						 replicates pay off the loan)
	forgiven			(an 1D array of the forgiven amounts of numReps replicates)
	penaltyTaxRate		(penalty tax rate for defaulting borrowers in the last 
						 period, default = 0.25)
	govRevenue			(an 1D array of government revenue from colleting penalty 
						 tax on defaulting borrowers)

********************************************************************************/
class RepaymentModel {
public:
	string name = "gen";
	int numReps;
	double interestRate;
	int numPers;
	double graceAmount;
	double loan;

	double penaltyTaxRate = 0.25;
	vector<vector<double>> balances;
	vector<vector<double>> payments;
	vector<bool> paidOff;
	vector<double> forgiven;
	vector<double> govRevenue;


	/********************************************************************************
	Purpose: Setters for some instance variables.

	********************************************************************************/
	void setNumReplicates(int numReps) { this->numReps = numReps; }

	void setInterestRate(double interestRate) { this->interestRate = interestRate; }

	void setNumPeriods(int numPers) { this->numPers = numPers; }

	void setGraceAmount(double graceAmount) { this->graceAmount = graceAmount; }

	void setLoan(double loan) { this->loan = loan; }


	/********************************************************************************
	Purpose: Initialize some instance variables.

	********************************************************************************/
	void setup() {

		balances.assign(numPers, vector<double>(numReps, 0));
		for (int j = 0; j < numReps; j++) {
			balances[0][j] = loan;
		}

		payments.assign(numPers, vector<double>(numReps, 0));
		paidOff.assign(numReps, true);
		forgiven.assign(numReps, 0);
		govRevenue.assign(numReps, 0);
	}


	/********************************************************************************
	Purpose: With a given level of income, determine the payment needed to make, and
		determine if income is less than required amount.

		This method in this class is just a prototype. Subclasses of specific
		repayment methods will contain detailed calculation of the required payment.

	Parameters:
		income 		(given level of income)

	Output:
					(a pair of real numbers: the required payment and the underpaid
					 amount)

	********************************************************************************/
	virtual pair<double, double> calcPayment(double income) {

		return make_pair(0, 0);
	}


	/********************************************************************************
	Purpose: Update the status of a replicate making payment in a given period. This 
		method will be inherited by all subclasses.

		Record the payment made, then update the remaining balances of the next period,
		or change the paid off status, if applicable.

		If it is the last period, determine the final paid off status and the amount 
		forgiven, if applicable. The penalty tax rate will be applied on amount 
		forgiven if it is still positive after the last period.

		Return consumption, which is the difference between income and payment made. 
		Deduct penalty tax, if applicable.

	Parameters:
		i 			(period to make payment)
		j 			(index of the replicate making payment)
		income 		(income in period i of replicate j)
		payment 	(the required payment)
		underpaid	(the difference between the required payment and standard payment)

	Returns:
					(consumption, which is the difference between income and payment 
					 made)

	********************************************************************************/
	double pay(int i, int j, double income, double payment, double underpaid) {

		if (payment < balances[i][j] * (1 + interestRate)) {
			payments[i][j] = payment;
		} else {
			underpaid = 0.0;
			payments[i][j] = balances[i][j] * (1 + interestRate);
		}

		if (underpaid > 0.0 && name == "std") {
			paidOff[j] = false;
			forgiven[j] += (underpaid * pow(1 + interestRate, -i));
		}

		double consumption = income - payments[i][j];
		assert(consumption > 0);

		if (i == numPers - 1) {
			if (fabs(payments[i][j] - balances[i][j] * (1 + interestRate)) < 1e-10) {
				if (name == "std") {
					assert(paidOff[j]);
				} else if (name == "ibr") {
					paidOff[j] = true;
				}

				return consumption;

			} else {
				if (name == "std") {
					assert(!paidOff[j]);
				} else if (name == "ibr") {
					paidOff[j] = false;
				}

				double rawForgiven;
				if (name == "ibr") {
					rawForgiven = balances[i][j] * (1 + interestRate) - 
								  payments[i][j];
					rawForgiven /= (1 + interestRate);
				} else {
					rawForgiven = underpaid;
					rawForgiven += forgiven[j] * pow(1 + interestRate, i);
				}

				double penalty = min(rawForgiven * penaltyTaxRate,
									 consumption - graceAmount);
				govRevenue[j] = penalty * pow(1 + interestRate, -i);
				forgiven[j] = (rawForgiven - penalty) * pow(1 + interestRate, -i);

				return consumption - penalty;
			}

		} else {
			balances[i + 1][j] = balances[i][j] * (1 + interestRate) - 
								 payments[i][j];
		}

		return consumption;
	}
};

#endif