//File containing defintions for functions relating to testing

#include <iostream>

#include "inout.h"

using std::cout;
using std::endl;

//Function to test that the dimensionless input file is read correctly
void In_test(dimless_in input)
{
  cout << "Viscosity Ratio = " << input.viscos_rat << endl;
  cout << "Bond number = " << input.bond << endl;
  cout << "Modified Density Ratio = " << input.mdr << endl;
  cout << "Number of intervals on sphere's surface = " << input.n_sphere << endl;
  cout << "Number of intervals on interface = " << input.n_int << endl;
  cout << "Truncation distance from axis = " << input.max_arc << endl;
  cout << "Initial time step = " << input.t_step << endl;
  cout << "Initial height of sphere = " << input.init_height << endl;
}
