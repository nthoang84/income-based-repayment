/********************************************************************************
IncomeBasedRepaymentModel.h

This header file provides the model for Standard Repayment Plan.

Version 1.0.0, May 15, 2022

@author: Hoang Nguyen

********************************************************************************/
#ifndef INCOME_BASED_REPAYMENT_MODEL_H
#define INCOME_BASED_REPAYMENT_MODEL_H

#include <math.h>

#include "RepaymentModel.h"

using namespace std;


/********************************************************************************
Purpose: A model for Income-Based Repayment Plan. This is a subclass of RepaymentModel.

Parameters:
	povertyThreshold 	(the poverty threshold)
	stdPayment			(the fixed amount that borrowers need to pay in any period,
						 calculated under the Standard Plan)

********************************************************************************/
class IncomeBasedRepaymentModel : public RepaymentModel {
public:
	double povertyThreshold;
	double stdPayment;


	/********************************************************************************
	Purpose: Initialize some instance variables.

	********************************************************************************/
	void setPovertyThreshold(double povertyThreshold) { this->povertyThreshold = povertyThreshold; }

	void setStandardPayment(double stdPayment) { this->stdPayment = stdPayment; }


	/********************************************************************************
	Purpose: Initialize some instance variables.

	********************************************************************************/
	void setup() {

		RepaymentModel::setup();
		name = "ibr";
	}


	/********************************************************************************
	Purpose: With a given level of income, determine the payment needed to make, and
		determine if income is less than required amount.
		
		In Income-Based Repayment Plan, the required payment amount in each period 
		is 10% of discretionary income or the required payment under Standard Repayment 
		Plan, whichever is less.

	Parameters:
		income 		(given level of income)

	Outputs:
					(a pair of real numbers: the required payment and the underpaid
					 amount)

	********************************************************************************/
	pair<double, double> calcPayment(double income) override {

		double discretionaryIncome = max(0.0, income - 1.5 * povertyThreshold);
		double payment = min(0.1 * discretionaryIncome, stdPayment);
		double underpaid = stdPayment - payment;
		assert(underpaid >= 0);
		return make_pair(payment, underpaid);
	}
};

#endif