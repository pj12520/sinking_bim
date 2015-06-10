//Header file containing functions relating to the solution of the linear system

#ifndef __SOLVE_H__
#define __SOLVE_H__
#pragma once

#include <vector>

using std::vector;

//Function to solve the linear system of equations
void Solve(int order, vector<vector<double> >* coeffs, vector<double>* known, vector<double>* unknown);


#endif /* SOLVE_H */
