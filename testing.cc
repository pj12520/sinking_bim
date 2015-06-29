//File containing defintions for functions relating to testing

#include <iostream>
#include <fstream>
#include <iomanip>

#include "inout.h"
#include "object.h"

using std::cout;
using std::endl;
using std:: ofstream;
using std::setw;

//Function to test that the dimensionless input file is read correctly
void In_test(dimless_in input)
{
  cout << "Viscosity Ratio = " << input.viscos_rat << endl;
  cout << "Bond number = " << input.bond << endl;
  cout << "Modified Density Ratio = " << input.mdr << endl;
  cout << "Number of intervals on sphere's surface = " << input.n_sphere << endl;
  cout << "Number of intervals on interface = " << input.n_int << endl;
  cout << "Truncation distance from axis = " << input.max_arc << endl;
  //  cout << "Initial time step = " << input.t_step << endl;
  cout << "Initial height of sphere = " << input.init_height << endl;
}


//Function to test that the sphere and interface are produced correctly
void Config(particle sphere, surf interf)
{
  ofstream sphere_out;
  sphere_out.open("testing/sphere_config.dat");

  sphere_out << setw(20) << "Theta" << setw(20) << "Radial" << setw(20) << "Vertical" << endl;

  for (int i = 0; i < sphere.n_int; i++)
    {
      for (int j = 0; j < 4; j++)
	{
	  sphere_out << setw(20) << sphere.intervals[i].theta[j] << setw(20) << sphere.intervals[i].rad[j] << setw(20) << sphere.intervals[i].vert[j] << endl;
	}
    }
  sphere_out.close();

  sphere_out.open("testing/sphere_interv.dat");

  sphere_out << setw(20) << "Interval" << setw(20) << "Lower Bound" << setw(20) << "Midpoint" << setw(20) << "Upper Bound" << endl;
  for (int i = 0; i < sphere.n_int; i++)
    {
      sphere_out << setw(20) << i << setw(20) << sphere.intervals[i].lower << setw(20) << sphere.midpoints[i] << setw(20) << sphere.intervals[i].upper << endl;
    }
  sphere_out.close();

  ofstream interf_out;
  interf_out.open("testing/interf_config.dat");

  interf_out << setw(20) << "Arc" << setw(20) << "Radial" << setw(20) << "Vertical" << setw(20) << "norm_rad" << setw(20) << "norm_vert" << setw(20) << "div_norm" << endl;

  for (int i = 0; i < interf.n_int; i++)
    {
      for (int j = 0; j < 4; j++)
	{
	  interf_out << setw(20) << interf.intervals[i].arc[j] << setw(20) << interf.intervals[i].rad[j] << setw(20) << interf.intervals[i].vert[j] << setw(20) << interf.intervals[i].norm_rad[j] << setw(20) << interf.intervals[i].norm_vert[j] << setw(20) << interf.intervals[i].div_norm[j] << endl;
	}
    }
  interf_out.close();

  interf_out.open("testing/interf_interv.dat");

  interf_out << setw(20) << "Interval" << setw(20) << "Lower Bound" << setw(20) << "Midpoint" << setw(20) << "Upper Bound" << setw(20) << "Radial" << setw(20) << "Vertical" << setw(20) << "norm_rad" << setw(20) << "norm_vert" << setw(20) << "div_norm" << endl;

  for (int i = 0; i < interf.n_int; i++)
    {
      interf_out << setw(20) << i << setw(20) << interf.intervals[i].lower << setw(20) << interf.midpoints[i] << setw(20) << interf.intervals[i].upper << setw(20) << interf.mid_rad[i] << setw(20) << interf.mid_vert[i] << setw(20) << interf.mid_norm_rad[i] << setw(20) << interf.mid_norm_vert[i] << setw(20) << interf.mid_div_norm[i] << endl;
    }
  interf_out.close();
}

//Function to print out the matrix of coefficients and the known vector so the elements can be compared with that calculated analytically
void Lin_sys_test(vector<vector<double> >* coeffs, vector<double>* known)
{
  ofstream matrix_out;
  matrix_out.open("testing/lin_sys/matrix.dat");

  for (int i = 0; i < (*coeffs).size(); i++)
    {
      for (int j = 0; j < (*coeffs)[i].size(); j++)
	{
	  matrix_out << setw(12) << (*coeffs)[i][j] << '\t';
	}
      matrix_out << endl;
    }

  matrix_out.close();

  ofstream vector_out;
  vector_out.open("testing/lin_sys/vector.dat");

  for (int i = 0; i < (*known).size(); i++)
    {
      vector_out << (*known)[i] << endl;
    }

  vector_out.close();
}
