/********************************************************************************
RepaymentSimulator.cpp

This is the main program which calls all methods of interest.

Version 1.0.0, May 15, 2022

@author: Hoang Nguyen

********************************************************************************/
#include <time.h>

#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

#include "TotalUtility.h"

using namespace std;
using Result = tuple<double, double, double, double, double, double>;


/********************************************************************************
Purpose: With given input parameters, create a TotalUtility object to simulate
    the model.

Parameters:
    beta            (discount factor for utility tradeoff across periods)
    sigma           (standard deviation of annual income shock)
    numReps         (number of replicates in the models)
    initPermInc     (initial permanent income, in thousands of dollars)
    loan            (initial loan amount)
    stdPers         (number of repayment periods in Standard Plan)
    ibrPers         (number of repayment periods in Income-Based Repayment Plan)
    name            (name of the dataset)

Output:
    results         (a tuple of six real numbers, consisting of stdUtil, ibrUtil,
                     diff, delta, insurance, ratio, in that order)

********************************************************************************/
Result simulate(double beta, double sigma, int numReps, double initPermInc,
                double loan, int stdPers, int ibrPers, string name) {

    TotalUtility driver;
    driver.setName(name);
    driver.setBeta(beta);
    driver.setSigma(sigma);
    driver.setNumReplicates(numReps);
    driver.setStdPeriods(stdPers);
    driver.setIbrPeriods(ibrPers);

    Result results = driver.compareExpectedUtilityStdAndIbr(initPermInc, loan, &driver);
    
    return results;
}


/********************************************************************************
Purpose: Print the output csv file.

Parameters:
    name            (name of the dataset)
    initPermInc     (initial permanent income, in thousands of dollars)
    loan            (initial loan amount)
    results         (a tuple of six real numbers, consisting of stdUtil, ibrUtil,
                     diff, delta, insurance, ratio, in that order)

Output:
    (none)

********************************************************************************/
void printResult(string name, double initPermInc, double loan, Result results) {

    double stdUtil, ibrUtil, diff, delta, insurance, ratio;
    tie(stdUtil, ibrUtil, diff, delta, insurance, ratio) = results;
    string strResults = to_string(stdUtil) + "," + 
                        to_string(ibrUtil) + "," + 
                        to_string(diff) + "," + 
                        to_string(delta) + "," +
                        to_string(insurance) + "," + 
                        to_string(ratio) + "\n";

    string outputFile = "Results/Outputs/output_" + name + ".csv";
    ofstream ofs;
    ofs.open(outputFile, ios::out);
    ofs << "stdUtil,ibrUtil,diff,delta,insurance,ratio\n";
    ofs << strResults;
    ofs.close();
}


/********************************************************************************
Purpose: Parse the input parameters from the command line then execute.

Parameters:
    argc            (number of command line arguments)
    argv            (list of command line arguments)

Output:
    (none)

********************************************************************************/
int main(int argc, char* argv[]) {

    double beta = stod(argv[1]);
	double sigma = stod(argv[2]);
	int numReps = stoi(argv[3]);
	double initPermInc = stod(argv[4]);
	double loan = stod(argv[5]);
	int stdPers = stoi(argv[6]);
	int ibrPers = stoi(argv[7]);
	string name = string(argv[8]);

    auto start = clock();

	Result results = simulate(beta, sigma, numReps, initPermInc, loan, stdPers,
				              ibrPers, name);

    printResult(name, initPermInc, loan, results);

    auto stop = clock();
    double runTime = (double) (stop - start) / CLOCKS_PER_SEC;
    cout << "Time elapsed: " << runTime << " s.\n";
    cout << '[' << name << "] Repayment Simulator finished!\n" << endl;

    return 0;
}