//Program to model the impact of a sphere onto a fluid fluid interface 

#include <string>

#include <iostream> //Currently only included for debugging purposes. Can be removed when program is functional
#include <iomanip> //Currently only included for debugging purposes. Can be removed when program is functional

#include "inout.h"
#include "object.h"
#include "testing.h"
#include "build.h"
#include "solve.h"
#include "vel.h"
//
using std::string;

using std::cout; //Currently only used for debugging purposes. Can be removed when program is functional
using std::endl; //Currently only used for debugging purposes. Can be removed when program is functional
using std::setw; //Currently only used for debugging purposes. Can be removed when program is functional

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
  //Config(sphere, interf);
  //Output values of the interval bounds and midpoints for both the sphere and interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  //Create the matrix of coefficients and the known vector for the linear system
  int order = 2 * (input.n_int + input.n_sphere - 1);

  vector<vector<double> > coeffs(order);
  for (int i = 0; i < order; i++)
    {
      coeffs[i].resize(order);
    }

  vector<double> known(order);

  //Start loop over iterations
  int it = 0;
  while(it < input.max_it)
    {
      //Build the linear system
      Build(&coeffs, &known, sphere, interf, input.viscos_rat, input.bond, input.mdr);

      //Testing-Test the creation of the linear system/////////////////////////////////
      Lin_sys_test(&coeffs, &known);
      //////////////////////////////////////////////////////////////////////////////////

      //Solve the linear system 
      vector<double> unknown(order);

      Solve(order, &coeffs, &known, &unknown);
      for (int i = 0; i < unknown.size(); i++)
	{
	  //cout << i << '\t' << unknown[i] << endl;
	}

      //Testing - Test the solution for the sphere velocity ///////////////////////////
      cout << setw(20) << it << setw(20) << sphere.height << setw(20) << unknown[unknown.size() - 1] << endl;
      ////////////////////////////////////////////////////////////////////////////////

      //Perform the 1st time step
      Iterate(input.n_int, &unknown, &interf.midpoints, &interf.mid_rad, &interf.mid_vert, &sphere.height, input.t_step);

      //Update the properties of the interface and sphere
      Up_interf(&interf);
      Up_sphere(&sphere);

      //Testing - Check the new configuration of the system/////////////////////////////////////////////////
      Config(sphere, interf);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      //Iterate the system
      it = it + 1;
    }
  return 0;
}

