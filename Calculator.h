/********************************************************************************
Calculator.h

This header file provides convenient calculation methods for other programs.

Version 1.0.0, May 15, 2022

@author: Hoang Nguyen

********************************************************************************/
#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <assert.h>
#include <math.h>

using namespace std;


/********************************************************************************
Purpose: This class contains methods for common calculations.

Parameters:
	gamma				(curvature parameter on consumption)

********************************************************************************/
class Calculator {
public:
	double gamma;


	/********************************************************************************
	Purpose: Setters for some instance variables.

	********************************************************************************/
	void setGamma(double gamma) { this->gamma = gamma; }


	/********************************************************************************
	Purpose: Calculate the utility of a given consumption

	Parameters:
		consumption		(consumption level as a scalar)

	Returns:
						(the utility of the given consumption level)
	********************************************************************************/
	double u(double consumption) {
		assert(consumption > 0);
		return ((pow(consumption, 1 - gamma) - 1) / (1 - gamma));
	}
};

#endif