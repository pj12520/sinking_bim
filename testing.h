//Header file containing functions which are used to test aspects of the model

#ifndef __TESTING_H__
#define __TESTING_H__
#pragma once

//Function to test that the dimensionless input file is read correctly
void In_test(dimless_in input);

//Function to test that the sphere and interface are produced correctly
void Config(particle sphere, surf interf);

//Function to print out the matrix of coefficients and the known vector so the elements can be compared with that calculated analytically
void Lin_sys_test(vector<vector<double> >* coeffs, vector<double>* known);

#endif /* TESTING_H */
