//Program to model the impact of a sphere onto a fluid fluid interface 

#include <string>

#include <iostream> //Currently only included for debugging purposes. Can be removed when program is functional

#include "inout.h"
#include "object.h"
#include "testing.h"

using std::string;

using std::cout; //Currently only used for debugging purposes. Can be removed when program is functional
using std::endl; //Currently only used for debugging purposes. Can be removed when program is functional

int main()
{
  //Read in input data from input file
  dimless_in input;
  string infile = "input.dat";
  Dimless_in(infile, &input);

  //TESTING - Test that input data is read correctly//////////////////////
  //In_test(input);
  ///////////////////////////////////////////////////////////////////////

  //Create sphere
  particle sphere;
  sphere.intervals.resize(input.n_sphere);
  Create_sphere(&sphere, input.init_height, input.n_sphere);

  //Create interface
  surf interf;
  interf.intervals.resize(input.n_int);
  Create_interf(&interf, input.n_int, input.max_arc);

  //Testing - Test that sphere and interface are constructed correctly///////////////////////////////////

  //Produce functions r(theta), z(theta), z(r) for the sphere and r(s), z(s) and z(r) for the interface
  Config(sphere, interf);
  //Output values of the interval bounds and midpoints for both the sphere and interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  return 0;
}

