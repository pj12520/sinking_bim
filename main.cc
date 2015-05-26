//Program to model the impact of a sphere onto a fluid fluid interface 

#include <string>

#include "inout.h"
#include "object.h"

using std::string;

int main()
{
  //Read in input data from input file
  dimless_in input;
  string infile = "input.dat";
  Dimless_in(infile, input);

  //Create sphere
  particle sphere;
  Create_sphere(sphere, input.init_height, input.n_sphere);

  return 0;
}

