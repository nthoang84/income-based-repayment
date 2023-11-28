/********************************************************************************
test.cpp

This program tests the behavior of interpolants and prints out all matrices of
interest.

Version 1.0.0, May 16, 2022

@author: Hoang Nguyen

********************************************************************************/
#include <bits/stdc++.h>

#include "TotalUtility.h"

using namespace std;

int main(int argc, char* argv[]) {

	int DEBUG_MODE = stoi(argv[1]);
	auto start = clock();

	// -------- Dataset info
	const vector<string> names = {
		"public", "privatenp", "privatefp", "dartecon", "darthist", "dartgovt", "dartpsyc", "maxloan", "unfortunate"
	};

	const vector<double> initPermIncs = {
		55.7, 55.7, 55.7, 85.2, 46.4, 59.9, 38.9, 73.25, 15
	};

	const vector<double> loans = {
		25.3, 27, 22.7, 16.83, 17.007, 14.000, 14.089, 57.5, 20
	};

	assert(names.size() == initPermIncs.size() &&
		   names.size() == loans.size());


	// -------- Matrix and interpolation parameters
	vector<double> initPermIncRange, loanRange;
	pair<int, int> initPermIncBound = (DEBUG_MODE ? make_pair(50, 60) : make_pair(15, 100));
	pair<int, int> loanBound = (DEBUG_MODE ? make_pair(20, 30) : make_pair(0, 60));
	for (int i = initPermIncBound.first; i < initPermIncBound.second; i++) {
		initPermIncRange.push_back((double) i);
	}
	for (int i = loanBound.first; i < loanBound.second; i++) {
		loanRange.push_back((double) i);
	}


	// -------- Print out utility and interpolation matrices
	TotalUtility driver;
	TotalUtility::Interpolant interpStd, interpIbr;

	for (string modelType : {"std", "ibr"}) {
		auto utilityMatrix = driver.createUtilityMatrix(modelType, initPermIncRange, loanRange);
		driver.printMatrix(initPermIncRange, loanRange, utilityMatrix, "UtilityMatrix_" + modelType);
		driver.printInterpolationMatrices(modelType, initPermIncRange, loanRange, utilityMatrix);
		auto& interp = (modelType == "std" ? interpStd : interpIbr);
		interp = TotalUtility::Interpolant(initPermIncRange, loanRange, utilityMatrix);
	}


	// -------- Testing for each dataset
	for (int i = 0; i < (DEBUG_MODE ? 1 : names.size()); i++) {
		string name = names[i];
		double initPermInc = initPermIncs[i];
		double loan = loans[i];

		ofstream ofs("Results/Consoles/console_" + name + ".txt");

		ofs << "initial permanent income = " << initPermInc << '\n';
		ofs << "loan = " << loan << '\n';
		ofs << '\n';

		ofs << "----- Standard Repayment -----" << '\n';
		auto stdModel = driver.createStdModel(loan);
		ofs << "utility (real) = " << driver.calcTotalUtility(&stdModel, initPermInc).first << '\n';
		ofs << "utility (interpolate) = " << interpStd.eval(initPermInc, loan) << '\n';
		ofs << "d_utility/d_Y0 = " << interpStd.evalDiff(initPermInc, loan, 0) << '\n';
		ofs << "d_utility/d_Q0 = " << interpStd.evalDiff(initPermInc, loan, 1) << '\n';

		ofs << '\n';

		ofs << "----- Income-Based Repayment -----" << '\n';
		auto ibrModel = driver.createIbrModel(loan, stdModel.fixedPayment);
		ofs << "utility (real) = " << driver.calcTotalUtility(&ibrModel, initPermInc).first << '\n';
		ofs << "utility (interpolate) = " << interpIbr.eval(initPermInc, loan) << '\n';
		ofs << "d_utility/d_Y0 = " << interpIbr.evalDiff(initPermInc, loan, 0) << '\n';
		ofs << "d_utility/d_Q0 = " << interpIbr.evalDiff(initPermInc, loan, 1) << '\n';

		ofs.close();

		cout << "[" << name << "] Testing finished.\n";
	}


	auto stop = clock();
    double runTime = (double) (stop - start) / CLOCKS_PER_SEC;
    cout << "Time elapsed: " << runTime << " s." << endl;
}