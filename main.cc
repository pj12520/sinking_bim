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
#include "interp_1d.h"
///#include "interface_interp.h"

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

  //Build the linear system
  Build(&coeffs, &known, sphere, interf, input.viscos_rat, input.bond, input.mdr);

  //Testing-Test the creation of the linear system/////////////////////////////////
  Lin_sys_test(&coeffs, &known);
  //////////////////////////////////////////////////////////////////////////////////

  //Solve the linear system 
  vector<double> unknown(order);

  Solve(order, &coeffs, &known, &unknown);
  // for (int i = 0; i < unknown.size(); i++)
  //{
  //  cout << i << '\t' << unknown[i] << endl;
  //}

  //Testing - Test the solution for the sphere velocity ///////////////////////////
  //  cout << setw(20) << input.viscos_rat << setw(20) << input.init_height << setw(20) << unknown[unknown.size() - 1] << endl;
  ////////////////////////////////////////////////////////////////////////////////

  //Perform the 1st time step
  Iterate(input.n_int, &unknown, &interf.midpoints, &interf.mid_rad, &interf.mid_vert, &sphere.height, input.t_step);

  //Describe the interface using a cubic spline
  Spline_interp rad_spline(interf.midpoints, interf.mid_rad, 1.0, 1.0); //Structures that contain the interpolation routines
  Spline_interp vert_spline(interf.midpoints, interf.mid_vert, 0.0, 0.0);

  //Find the new maximum value of the arc-length


  //Testing - Output the new configuration of the system/////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////

  //Redistribute the points
  double max_arc = interf.midpoints[input.n_int - 1];

  vector<double> new_midpoints(input.n_int);
  vector<double> new_mid_rad(input.n_int);
  vector<double> new_mid_vert(input.n_int);

  Find_midpoints(&new_midpoints, 0.0, max_arc, input.n_int);

  for (int i = 0; i < input.n_int; i++)
    {
      if (i != 0)
	{
	  new_mid_rad[i] = rad_spline.interp(new_midpoints[i]);
	}
      else
	{
	  new_mid_rad[i] = 0.0;
	}
      new_mid_vert[i] = vert_spline.interp(new_midpoints[i]);
    }

  //Find the new intervals
  double half_width = max_arc / (2.0 * (input.n_int - 1)); //This is half the width of the intermediate intervals. 

  for (int i = 0; i < input.n_int; i++)
    {
      //Create the abscissas for the Gauss-Legendre integration in each interval
      Abscissas(&interf.intervals[i].lower, &interf.intervals[i].upper, max_arc, input.n_int, &interf.intervals[i].arc, &interf.intervals[i].width, half_width, i);

      //Set radial and vertical components of the integration points in each interval and components and divergence of normal vectors
      for (int j = 0; j < 4; j++)
	{
	  interf.intervals[i].rad[j] = rad_spline.interp(interf.intervals[i].arc[j]);
	  interf.intervals[i].vert[j] = vert_spline.interp(interf.intervals[i].arc[j]);

	  //	  (*intervals)[i].norm_rad[j] = 0.0;
	  //	  (*intervals)[i].norm_vert[j] = 1.0;
	  //	  (*intervals)[i].div_norm[j] = 0.0;
	}
    }


  //Start for loop  


  //Build linear system

  //Solve linear system

  //Iterate system in time

  return 0;
}

