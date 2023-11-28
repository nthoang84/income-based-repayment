/********************************************************************************
TotalUtility.h

This header file manipulates different repayment models of interest.

Version 1.0.1, May 18, 2022

@author: Hoang Nguyen

********************************************************************************/
#ifndef TOTAL_UTILITY_H
#define TOTAL_UTILITY_H

#include <math.h>

#include <string>

#include <boost/math/tools/roots.hpp>

#include "alglib-cpp/src/stdafx.h"
#include "alglib-cpp/src/interpolation.h"

#include "AssetDistributions.h"
#include "IncomeBasedRepaymentModel.h"
#include "StandardRepaymentModel.h"

using namespace std;
using namespace alglib;
using Result = tuple<double, double, double, double, double, double>;


/********************************************************************************
Purpose: General utility class to manipulate different repayment models.

Parameters:
	name            	(name of the dataset, default = "name")
	interestRate 		(annual interest rate of the loan, default = 3.73%)
    beta            	(discount factor for utility tradeoff across periods, 
    					 default = 0.97)
    gamma 				(curvature parameter on consumption, default = 2)
    rho 				(AR(1) coefficient for income shock decay, default = 0.95)
    sigma           	(standard deviation of annual income shock, default = 0.15)
    numReps         	(number of replicates in the models, default = 1000)
	povertyThreshold 	(annual poverty threshold in thousands of dollar, 
						 default = 21.96)
	incGrowthRate 		(annual income growth rate, default = 0.01)
    stdPers    		    (number of repayment periods in Standard Plan, 
    					 default = 10)
    ibrPers        		(number of repayment periods in Income-Based Repayment 
    					 Plan, default = 20)
    graceAmount 		(minimum possible amount for consumption each period, 
    					 default = 0.001)

********************************************************************************/
class TotalUtility {
public:
	string name = "name";
	double interestRate = 0.0373;
	double beta = 0.97;
	double gamma = 2.0;
	double rho = 0.95;
	double sigma = 0.15;
	int numReps = 1000;
	double povertyThreshold = 21.96;
	double incGrowthRate = 0.01;
	int stdPers = 10;
	int ibrPers = 20;
	double graceAmount = 0.001;


	/********************************************************************************
	Purpose: Setters for some instance variables.

	********************************************************************************/
	void setName(string name) { this->name = name; }

	void setBeta(double beta) { this->beta = beta; }

	void setSigma (double sigma) { this->sigma = sigma; }

	void setNumReplicates(int numReps) { this->numReps = numReps; }
	
	void setStdPeriods(int stdPers) { this->stdPers = stdPers; }

	void setIbrPeriods(int ibrPers) { this->ibrPers = ibrPers; }


	/********************************************************************************
	Purpose: Create a Standard Repayment model with given amount of loan.

	Parameters:
	    loan            (initial loan amount)

	Output:
	    model         	(a StandardRepaymentModel object)

	********************************************************************************/
	StandardRepaymentModel createStdModel(double loan) {

		StandardRepaymentModel model;
		model.setNumReplicates(numReps);
		model.setInterestRate(interestRate);
		model.setNumPeriods(stdPers);
		model.setGraceAmount(graceAmount);
		model.setLoan(loan);

		model.setup();
		return model;
	}


	/********************************************************************************
	Purpose: Create an Income-Based Repayment model with given amount of loan and
		the fixed payment under Standard Plan.

	Parameters:
	    loan            (initial loan amount)
	    stdPayment 		(the fixed payment each period under the Standard Plan)

	Output:
	    model         	(an IncomeBasedRepaymentModel object)

	********************************************************************************/
	IncomeBasedRepaymentModel createIbrModel(double loan, double stdPayment) {

		IncomeBasedRepaymentModel model;
		model.setNumReplicates(numReps);
		model.setInterestRate(interestRate);
		model.setNumPeriods(ibrPers);
		model.setGraceAmount(graceAmount);
		model.setLoan(loan);

		model.setPovertyThreshold(povertyThreshold);
		model.setStandardPayment(stdPayment);

		model.setup();
		return model;
	}


	/********************************************************************************
	Purpose: Create the asset distributions from the given repayment model and 
		the initial permanent income amount.

	Parameters:
	    model 			(a pointer to the repayment model of interest)
	    initPermInc 	(initial permanent income)

	Output:
	    model         	(an AssetDistributions object)

	********************************************************************************/
	AssetDistributions createAssetDistributions(RepaymentModel* model, 
												double initPermInc) {

		AssetDistributions dist;
		dist.setModel(model);
		dist.setBeta(beta);
		dist.setSigma(sigma);
		dist.setRho(rho);
		dist.setGamma(gamma);
		dist.setCommonNumPeriods(max(stdPers, ibrPers));
		dist.setInitialPermanentIncome(initPermInc);
		dist.setIncomeGrowthRate(incGrowthRate);

		dist.setup();
		return dist;
	}


	/********************************************************************************
	Purpose: Calculate the total utility of a given repayment model.

	Parameters:
	    model 			(a pointer to the repayment model of interest)
	    initPermInc 	(initial permanent income)
	    isPrinted 		(a flag indicating whether outputs and schedules are printed
	    				 out, default = false)

	Output:
	    	         	(a pair of real numbers: total utility and the average
	    	         	 amount forgiven)

	********************************************************************************/
	pair<double, double> calcTotalUtility(RepaymentModel* model, 
										  double initPermInc,
										  bool isPrinted = false) {

		AssetDistributions dist = createAssetDistributions(model, initPermInc);
		if (isPrinted) {
			dist.summarize(name);
			dist.printSchedule(name);
		}
		return make_pair(dist.calcTotalUtility(), dist.calcAverageForgiven());
	}


	/********************************************************************************
	Purpose: This struct serves as a functor to feed the root-finding method later.

	********************************************************************************/
	struct CompensatingVariation {

		double ibrUtil, loan, initPermInc;
		TotalUtility* driver;

		CompensatingVariation(const double& ibrUtil, const double& loan,
							  const double& initPermInc, 
							  TotalUtility* driver) {

			this->ibrUtil = ibrUtil;
			this->loan = loan;
			this->initPermInc = initPermInc;
			this->driver = driver;
		}

		double operator()(const double& delta) {

			StandardRepaymentModel model = driver->createStdModel(loan - delta);
			double diff = ibrUtil - 
						  driver->calcTotalUtility(&model, initPermInc).first;
			return diff;
		}
	};


	/********************************************************************************
	Purpose: Compare the utility of both repayment models, then find the compensating
		variation.

	Parameters:
	    initPermInc     (initial permanent income)
    	loan            (initial loan amount)
    	driver 			(a pointer to the TotalUtility object already created in 
    					 RepaymentSimulator.cpp)

	Output:
	    	         	(a tuple of six numbers: the utility under Standard Plan,
	    	         	 the utility under Income-Based Repayment Plan, the
	    	         	 difference in utility between two plans, the compensating
	    	         	 variation, the insurance value, and the ratio CV/loan)

	********************************************************************************/
	Result compareExpectedUtilityStdAndIbr(double initPermInc, double loan, TotalUtility* driver) {

		StandardRepaymentModel stdModel = createStdModel(loan);
		double stdUtil, stdAvgForgiven;
		tie(stdUtil, stdAvgForgiven) = calcTotalUtility(&stdModel, initPermInc, true);

		IncomeBasedRepaymentModel ibrModel = createIbrModel(loan,
			StandardRepaymentModel::calcFixedPayment(interestRate, stdPers, loan));
		double ibrUtil, ibrAvgForgiven;
		tie(ibrUtil, ibrAvgForgiven) = calcTotalUtility(&ibrModel, initPermInc, true);

		double diff = stdUtil - ibrUtil;

		using namespace boost::math::tools;
		CompensatingVariation f(ibrUtil, loan, initPermInc, driver);
		boost::uintmax_t maxIter = 20;
		int digits = std::numeric_limits<double>::digits; 
		eps_tolerance<double> tol(digits - 3);

		const double tolerance = 1e-10;
		pair<double, double> bounds;
		if (fabs(diff) < tolerance) {
			bounds = make_pair(-tolerance, tolerance);
			cout << "WARNING: equal utility, not as expected!" << '\n';
		} else if (stdUtil < ibrUtil) {
			bounds = make_pair(tolerance, loan - tolerance);
			cout << "stdUtil < ibrUtil, as expected!" << '\n';
		} else {
			bounds = make_pair(-loan * 3 + tolerance, -tolerance);
			cout << "WARNING: stdUtil > ibrUtil, not as expected!" << '\n';
		}

		pair<double, double> bracket = bisect(f, bounds.first, bounds.second, tol, maxIter);
		double delta = bracket.first;
		double insurance = delta - ibrAvgForgiven;
		double ratio = delta / loan;

		/* ---- Check the utility of stdModel after accounting for CV

		StandardRepaymentModel stdModelCheck = createStdModel(loan - delta);
		double stdUtilCheck, stdAvgForgivenCheck;
		tie(stdUtilCheck, stdAvgForgivenCheck) = calcTotalUtility(&stdModelCheck, initPermInc, true);

		cout << "CHECK: " << stdUtil << ' ' << ibrUtil << ' ' << stdUtilCheck << '\n';
		*/

		return make_tuple(stdUtil, ibrUtil, diff, delta, insurance, ratio);
	}


	/********************************************************************************
	Purpose: Print a given matrix with rows and columns data.

	Parameters:
	    rows     	(the row data)
    	cols        (the column data)
    	matrix 		(the matrix)
    	filename	(name of the output file)

	Output:
	    (none)

	********************************************************************************/
	void printMatrix(const vector<double>& rows, const vector<double>& cols,
					 const vector<vector<double>>& matrix, string filename) {

		int nRows = rows.size();
		int nCols = cols.size();

		string outputFile = "Results/Matrices/" + filename + ".csv";
	    ofstream ofs;
	    ofs.open(outputFile);
	    
	    string header;
	    for (int j = 0; j < nCols; j++) {
	    	header += ",";
	    	header += to_string(cols[j]);
	    }
	    ofs << header << '\n';

	    for (int i = 0; i < nRows; i++) {
	    	string info = to_string(rows[i]);
	    	for (int j = 0; j < nCols; j++) {
	    		info += ",";
	    		info += to_string(matrix[i][j]);
	    	}
	    	ofs << info << '\n';
	    }

	    ofs.close();
	}


	/********************************************************************************
	Purpose: Create a matrix of utility of the given repayment plan across different 
		initial permanent incomes and loans.

	Parameters:
		modelType		(the type of the model, "std" or "ibr")
	    initPermIncs    (an 1D array of initial permanent incomes)
    	loans           (an 1D array of initial loans)

	Output:
	 	utilityMatrix 	(the utility matrix)

	********************************************************************************/
	vector<vector<double>> createUtilityMatrix(string modelType,
							 				   const vector<double>& initPermIncs,
											   const vector<double>& loans) {

		int nRows = initPermIncs.size();
		int nCols = loans.size();
		vector<vector<double>> utilityMatrix(nRows, vector<double>(nCols));

		for (int i = 0; i < nRows; i++) {
			for (int j = 0; j < nCols; j++) {
				double initPermInc = initPermIncs[i];
				double loan = loans[j];
				double utility;

				if (modelType == "std") {
					StandardRepaymentModel model = createStdModel(loan);
					
					utility = calcTotalUtility(&model, initPermInc).first;
				} else {
					IncomeBasedRepaymentModel model = createIbrModel(
						loan,
						StandardRepaymentModel::calcFixedPayment(interestRate, stdPers, loan));
					
					utility = calcTotalUtility(&model, initPermInc).first;
				}

				utilityMatrix[i][j] = utility;
			}
		}

		return utilityMatrix;
	}


	/********************************************************************************
	Purpose: Given row, column, and matrix data, provide all interpolants of interest.

	********************************************************************************/
	struct Interpolant {

		spline2dinterpolant interp;
		vector<double> rows;
		vector<double> cols;
		vector<vector<double>> matrix;

		Interpolant() { }

		Interpolant(const vector<double>& rows,
					const vector<double>& cols,
					const vector<vector<double>>& matrix,
					int splineDegree = 1) {

			assert(splineDegree == 1 || splineDegree == 3);

			int nRows = rows.size();
			int nCols = cols.size();

			real_1d_array x = convert(rows).c_str();
			real_1d_array y = convert(cols).c_str();
			real_1d_array f = convert(matrix).c_str();

			if (splineDegree == 1) {
				spline2dbuildbilinearv(x, nRows, y, nCols, f, 1, interp);
			} else {
				spline2dbuildbicubicv(x, nRows, y, nCols, f, 1, interp);
			}
		}

		double eval(double x, double y) {

			return spline2dcalc(interp, x, y);
		}

		double evalDiff(double x, double y, int type) {

			assert(type == 0 || type == 1);
			double f, fx, fy, fxy;
			spline2ddiff(interp, x, y, f, fx, fy, fxy);
			return (!type ? fx : fy);
		}

		string convert(const vector<double>& v) {

			string s = "[";
			for (int i = 0; i < (int) v.size(); i++) {
				s += to_string(v[i]);
				s += (i == (int) v.size() - 1 ? "]" : ", ");
			}
			return s;
		}

		string convert(const vector<vector<double>>& m) {

			assert(!m.empty());
			vector<double> v;
			for (int j = 0; j < m[0].size(); j++) {
				for (int i = 0; i < m.size(); i++) {
					v.push_back(m[i][j]);
				}
			}
			assert(v.size() == m.size() * m[0].size());
			return convert(v);
		}
	};


	/********************************************************************************
	Purpose: Print three matrices: interpolation of utility and two first order partial
		derivatives.

	Parameters:
		modelType		(the type of the model, "std" or "ibr")
	    initPermIncs    (an 1D array of initial permanent incomes)
    	loans           (an 1D array of initial loans)
    	utilityMatrix 	(a 2D pre-calculated utility matrix)

	Output:
	 	(none)

	********************************************************************************/
	void printInterpolationMatrices(string modelType,
									const vector<double>& initPermIncs,
									const vector<double>& loans,
								  	const vector<vector<double>>& utilityMatrix) {

		assert(modelType == "std" || modelType == "ibr");
		Interpolant interp(initPermIncs, loans, utilityMatrix);

		vector<vector<double>> utility(initPermIncs.size());
		vector<vector<double>> utilityWrtInitPermIncs(initPermIncs.size());
		vector<vector<double>> utilityWrtLoans(initPermIncs.size());

		for (int i = 0; i < initPermIncs.size(); i++) {
			utility[i].resize(loans.size());
			utilityWrtInitPermIncs[i].resize(loans.size());
			utilityWrtLoans[i].resize(loans.size());

			for (int j = 0; j < loans.size(); j++) {
				double initPermInc = initPermIncs[i];
				double loan = loans[j];

				utility[i][j] = interp.eval(initPermInc, loan);
				utilityWrtInitPermIncs[i][j] = interp.evalDiff(initPermInc, loan, 0);
				utilityWrtLoans[i][j] = interp.evalDiff(initPermInc, loan, 1);
			}
		}

		printMatrix(initPermIncs, loans, utility, "InterpUtility_" + modelType);
		printMatrix(initPermIncs, loans, utilityWrtInitPermIncs, "InterpUtilityWrtInitPermIncs_" + modelType);
		printMatrix(initPermIncs, loans, utilityWrtLoans, "InterpUtilityWrtLoans_" + modelType);
	}
};

#endif