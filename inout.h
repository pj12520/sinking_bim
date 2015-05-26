//Header file containing declarations for functions and structures that relate to the input and output of data

#ifndef __IN_OUT_H__
#define __IN_OUT_H__
#pragma once

#include <string>

using std::string;
\
struct dimless_in
{
  double viscos_rat; //Viscosity Ratio
  double bond; //Bond number
  double mdr; //Modified Density Ratio
  int n_sphere; //Number of elements on sphere surface
  int n_int; //Number of elements on interface
  double max_arc; //Truncation distance along interface
  double t_step; //Initial time step
}

void Dimless_in(string file, dimless_in input); //Function to read in the input variables for the case that they are the dimensionless parameters

#endif /* IN_OUT_H */
