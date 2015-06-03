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
  cout << "Initial time step = " << input.t_step << endl;
  cout << "Initial height of sphere = " << input.init_height << endl;
}


//Function to test that the sphere and interface are produced correctly
void Config(particle sphere, surf interf)
{
  ofstream sphere_out;
  sphere_out.open("testing/sphere_config.dat");

  sphere_out << "Theta" << '\t' << "Radial" << '\t' << "Vertical" << endl;

  for (int i = 0; i < sphere.n_int; i++)
    {
      for (int j = 0; j < 4; j++)
	{
	  sphere_out << sphere.intervals[i].theta[j] << '\t' << sphere.intervals[i].rad[j] << '\t' << sphere.intervals[i].vert[j] << endl;
	}
    }
  sphere_out.close();

  sphere_out.open("testing/sphere_interv.dat");

  sphere_out << "Interval" << '\t' << "Lower Bound" << '\t' << "Midpoint" << '\t' << "Upper Bound" << endl;
  for (int i = 0; i < sphere.n_int; i++)
    {
      sphere_out << i << '\t' << sphere.intervals[i].lower << '\t' << sphere.midpoints[i] << '\t' << sphere.intervals[i].upper << endl;
    }
  sphere_out.close();

  ofstream interf_out;
  interf_out.open("testing/interf_config.dat");

  interf_out << "Arc" << '\t' << "Radial" << '\t' << "Vertical" << endl;

  for (int i = 0; i < interf.n_int; i++)
    {
      for (int j = 0; j < 4; j++)
	{
	  interf_out << interf.intervals[i].arc[j] << '\t' << interf.intervals[i].rad[j] << '\t' << interf.intervals[i].vert[j] << endl;
	}
    }
  interf_out.close();

  interf_out.open("testing/interf_interv.dat");

  interf_out << "Interval" << '\t' << "Lower Bound" << '\t' << "Midpoint" << '\t' << "Upper Bound" << endl;
  for (int i = 0; i < interf.n_int; i++)
    {
      interf_out << i << '\t' << interf.intervals[i].lower << '\t' << interf.midpoints[i] << '\t' << interf.intervals[i].upper << endl;
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
