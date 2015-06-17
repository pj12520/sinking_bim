//File containing defintions for functions relating to input and output

#include <fstream>
#include <string>

#include "inout.h"

using std::fstream;
using std::string;

//Function to read in the input variables for the case that they are the dimensionless parameters
void Dimless_in(string file, dimless_in *input) 
{
  fstream read; //Object for reading input parameters from file
  read.open(file.c_str()); 

  read >> input->viscos_rat; //Viscosity Ratio
  read >> input->bond; //Bond number
  read >> input->mdr; //Modified Density Ratio

  read >> input->n_sphere; //Number of elements on sphere's surface
  read >> input->n_int; //Number of elements on interface
  read >> input->max_arc; //Truncation length of interface
  read >> input->t_step; //Initial time step
  read >> input->init_height; //Initial height of sphere
  read >> input ->max_it; //Maximum number of iterations
  read.close();

}
