//Header file containing declarations for functions and structures that relate to the input and output of data

#ifndef __INOUT_H__
#define __INOUT_H__
#pragma once

#include <string>
#include <fstream>

using std::string;

//Structure to contain input parameters in case where they are dimensionless parameters
struct dimless_in 
{
  double viscos_rat; //Viscosity Ratio
  double bond; //Bond number
  double mdr; //Modified Density Ratio
  int n_sphere; //Number of elements on sphere surface
  int n_int; //Number of elements on interface
  double max_arc; //Truncation distance along interface
  double t_step; //Initial time step
  double init_height; //Initial height of sphere
};

//Function to read in the input variables for the case that they are the dimensionless parameters
void Dimless_in(string file, dimless_in input); 

#endif /* INOUT_H */
