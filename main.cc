//Program to model the impact of a sphere onto a fluid fluid interface 

#include <string>

#include "inout.h"

using std::string;

int main()
{
  //Read in input data from input file
  dimless_in input;
  string infile = "input.dat";
  Dimless_in(infile, input);

  return 0;
}

