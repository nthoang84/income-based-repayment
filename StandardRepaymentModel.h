/********************************************************************************
StandardRepaymentModel.h

This header file provides the model for Standard Repayment Plan.

Version 1.0.0, May 15, 2022

@author: Hoang Nguyen

********************************************************************************/
#ifndef STANDARD_REPAYMENT_MODEL_H
#define STANDARD_REPAYMENT_MODEL_H

#include <math.h>

#include "RepaymentModel.h"

using namespace std;


/********************************************************************************
Purpose: A model for Standard Repayment Plan. This is a subclass of RepaymentModel.

Parameters:
	fixed_payment		(the fixed amount that borrowers need to pay in any period)

********************************************************************************/
class StandardRepaymentModel : public RepaymentModel {
public:
	double fixedPayment;


	/********************************************************************************
	Purpose: With given interest rate, number of repayment periods, and initial loan,
		calculate the fixed payment required for each period

		This is a static method of the class StandardRepaymentModel

	Parameters:
		interestRate 		(annual interest rate of the loan)
		numPers 			(number of loan repayment periods)
		loan				(initial loan)

	Outputs:
		fixedPayment		(the fixed amount that borrowers need to pay in any period)

	********************************************************************************/
	static double calcFixedPayment(double interestRate, int numPers, double loan) {

		double fixedPayment = loan * interestRate * pow(1 + interestRate, numPers) /
					   (pow(1 + interestRate, numPers) - 1.0);

		return fixedPayment;
	}


	/********************************************************************************
	Purpose: Initialize some instance variables.

	********************************************************************************/
	void setup() {

		RepaymentModel::setup();
		name = "std";
		fixedPayment = calcFixedPayment(interestRate, numPers, loan);
	}


	/********************************************************************************
	Purpose: With a given level of income, determine the payment needed to make, and
		determine if income is less than required amount.
		
		In Standard Repayment Plan, the required payment amount each period is fixed.

	Parameters:
		income 		(given level of income)

	Outputs:
					(a pair of real numbers: the required payment and the underpaid
					 amount)

	********************************************************************************/
	pair<double, double> calcPayment(double income) override {

		double payment = min(income - graceAmount, fixedPayment);
		double underpaid = fixedPayment - payment;
		assert(underpaid >= 0);
		return make_pair(payment, underpaid);
	}
};

#endif