#include <fstream>
#include <string>

#include "inout.h"

using std::fstream;
using std::string;

void Dimless_in(string file, dimless_in input) //Function to read in the input variables for the case that they are the dimensionless parameters
{
  fstream read; //Object for reading input parameters from file
  read.open(file.c_str()); 

  read >> input.viscos_rat; //Viscosity Ratio
  read >> input.bond; //Bond number
  read >> input.mdr; //Modified Density Ratio

  read >> input.n_sphere; //Number of elements on sphere's surface
  read >> input.n_int; //Number of elements on interface
  read >> input.max_arc; //Truncation length of interface
  read >> input.t_step; //Initial time step

  read.close();

}
