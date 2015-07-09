//Program to model the impact of a sphere onto a fluid fluid interface 

#include <string>
#include <fstream>
#include <sstream>

#include <iostream> //Currently only included for debugging purposes. Can be removed when program is functional
#include <iomanip> //Currently only included for debugging purposes. Can be removed when program is functional

#include "inout.h"
#include "object.h"
#include "testing.h"
#include "build.h"
#include "solve.h"
#include "vel.h"
#include "geo.h"

using std::string;
using std::ofstream;
using std::ostringstream;

using std::cout; //Currently only used for debugging purposes. Can be removed when program is functional
using std::endl; //Currently only used for debugging purposes. Can be removed when program is functional
using std::setw; //Currently only used for debugging purposes. Can be removed when program is functional

int main(int argc, char *argv[])
{
  //Read in input data from input file
  dimless_in input;
  string infile(argv[1]);
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
  double time = 0.0;

  double t_step = min(0.01, 0.01 * input.mdr * input.bond);

  //Vectors to store the interfacial velocity
  vector<double> rad_vel(input.n_int);
  vector<double> vert_vel(input.n_int);


  //Object to output the position and velocity of the sphere
  ofstream sphere_out;
  //  string out_file = "D=" + static_cast<ostringstream*>( &(ostringstream() << input.mdr) )->str() + "/Bo=" + static_cast<ostringstream*>( &(ostringstream() << input.bond) )->str() + "/viscos_rat=" + static_cast<ostringstream*>( &(ostringstream() << input.viscos_rat) )->str() + "/sphere.dat";
  string out_file = "sphere.dat";
  sphere_out.open(out_file.c_str());
  sphere_out << setw(20) << "iteration" << setw(20) << "time" << setw(20) << "height" << setw(20) << "velocity" << endl;

  //Vectors of output data structures
  vector<out_data> output;

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

      //Update the vectors storing the interfacial velocities
      for (int i = 0; i < input.n_int; i++)
	{
	  if (i == 0)
	    {
	      rad_vel[i] = 0.0;
	      vert_vel[i] = unknown[i + input.n_int - 1];
	    }
	  else
	    {
	      rad_vel[i] = unknown[i - 1];
	      vert_vel[i] = unknown[i + input.n_int - 1];
	    }
	}

      //Output the configuration
      //      Out_sys(it, sphere, interf, input.mdr, input.bond, input.viscos_rat, &rad_vel, &vert_vel);

      //Store the configuration
      output.resize(output.size() + 1);

      output[output.size() - 1].it = it;
      output[output.size() - 1].sphere_pos = sphere.height;

      for (int i = 0; i < interf.n_int; i++)
	{
	  output[output.size() - 1].interf_rad.push_back(interf.mid_rad[i]);
	  output[output.size() - 1].interf_vert.push_back(interf.mid_vert[i]);
	}

      for (int i = 0; i < unknown.size(); i++)
	{
	  //	  	  cout << i << '\t' << unknown[i] << endl;
	}

      //Testing - Test the solution for the sphere velocity ///////////////////////////
      cout << setw(20) << it << setw(20) << time << setw(20) << sphere.height << setw(20) << unknown[unknown.size() - 1] << endl;
      ////////////////////////////////////////////////////////////////////////////////

      //Outout the sphere position and velocity
      sphere_out << setw(20) << it << setw(20) << time << setw(20) << sphere.height << setw(20) << unknown[unknown.size() - 1] << endl;
      //Perform the 1st time step
      Iterate(input.n_int, &unknown, &interf.midpoints, &interf.mid_rad, &interf.mid_vert, &sphere.height, t_step);

      //Iterate the system
      it = it + 1;
      time = time + t_step;

      //Update the properties of the interface and sphere
      Up_interf(&interf);
      Up_sphere(&sphere);

      //Testing - Check the new configuration of the system/////////////////////////////////////////////////
      Config(sphere, interf);
      ////////////////////////////////////////////////////////////////////////////////////////////////////////

      //Check the separation between the sphere and the interface isn't smaller than the distance between collocation points
      int break_criteria = Break_Crit(&interf.midpoints, &interf.mid_rad, &interf.mid_vert, sphere.height);

      if (break_criteria == 1)
	{
	  sphere_out << "# Sphere and interface collide" << endl;
	  break;
	}
      else if (break_criteria == 2)
	{
	  sphere_out << "# Tail snaps" << endl;
	  break;
	}

      //Output the system
      //    Out_sys(it, sphere, interf, input.mdr, input.bond, input.viscos_rat);
    }

  sphere_out.close();

  //Need to output the configuration at 20 iterations plus the initial and final ones
  int n_it = it - 1; //Number of iterations that occured

  if (n_it < 20)
    {
      vector<int> out_it(n_it + 1);

      for (int i = 0; i < out_it.size(); i++)
	{
	  out_it[i] = i;

	  Out_sys(output[out_it[i]], input.mdr, input.bond, input.viscos_rat);
	}

    }
  else if (n_it % 20 == 0)
    {
      vector<int> out_it(21);

      out_it[0] = 0;

      for (int i = 1; i < out_it.size(); i++)
	{
	  out_it[i] = n_it * i / 20;
	}

      for (int i = 0; i < out_it.size(); i++)
	{
	  Out_sys(output[out_it[i]], input.mdr, input.bond, input.viscos_rat);
	}
    }
  else
    {
      vector<int> out_it(22);

      out_it[0] = 0;
      out_it[21] = n_it;

      int rem = n_it % 20;

      for (int i = 1; i < out_it.size() - 1; i++)
	{
	  out_it[i] = (n_it - rem) * i / 20;
	}

      for (int i = 0; i < out_it.size(); i++)
	{
	  Out_sys(output[out_it[i]], input.mdr, input.bond, input.viscos_rat);
	}
    }

  return 0;
}

