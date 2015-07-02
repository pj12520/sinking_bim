//File containing defintions for functions relating to input and output

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "inout.h"
#include "object.h"

using std::fstream;
using std::string;
using std::ofstream;
using std::ostringstream;
using std::setw;
using std::endl;

//Function to read in the input variables for the case that they are the dimensionless parameters
void Dimless_in(string file, dimless_in *input) 
{
  fstream read; //Object for reading input parameters from file
  read.open(file.c_str()); 

  read >> input->bond; //Viscosity Ratio
  read >> input->mdr; //Bond number
  read >> input->viscos_rat; //Modified Density Ratio

  read >> input->n_sphere; //Number of elements on sphere's surface
  read >> input->n_int; //Number of elements on interface
  read >> input->max_arc; //Truncation length of interface
  //  read >> input->t_step; //Initial time step
  read >> input->init_height; //Initial height of sphere
  read >> input ->max_it; //Maximum number of iterations
  read.close();

}

//Function to output the state of the system
void Out_sys(int it, particle sphere, surf interf, double mdr, double bond, double viscos_rat, vector<double>* rad_vel, vector<double>* vert_vel)
{
  ofstream write;

  //  string file = "D=" + static_cast<ostringstream*>( &(ostringstream() << mdr) )->str() + "/Bo=" + static_cast<ostringstream*>( &(ostringstream() << bond) )->str() + "/viscos_rat=" + static_cast<ostringstream*>( &(ostringstream() << viscos_rat) )->str() + "/interf_config" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";
  string file = "interf_config.dat" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";
  write.open(file.c_str());

  write << setw(20) << "Interval" << setw(20) << "Arc" << setw(20) << "Radial" << setw(20) << "Vertical" << setw(20) << "Radial_vel" << setw(20) << "Vertical_vel" << endl;

  for (int i = 0; i < interf.n_int; i++)
    {
      write << setw(20) << i << setw(20) << interf.midpoints[i] << setw(20) << interf.mid_rad[i] << setw(20) << interf.mid_vert[i] << setw(20) << (*rad_vel)[i] << setw(20) << (*vert_vel)[i] << endl;
    }

  write.close();

  //  file = "D=" + static_cast<ostringstream*>( &(ostringstream() << mdr) )->str() + "/Bo=" + static_cast<ostringstream*>( &(ostringstream() << bond) )->str() + "/viscos_rat=" + static_cast<ostringstream*>( &(ostringstream() << viscos_rat) )->str() + "/sphere_config" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";
  file = "sphere_config" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";
  write.open(file.c_str());

  write << setw(20) << "Interval" << setw(20) << "Theta" << setw(20) << "Radial" << setw(20) << "Vertical" << endl;

  for (int i = 0; i < sphere.n_int; i++)
    {
      if (i == 0)
	{
	  write << setw(20) << i << setw(20) << sphere.midpoints[i] << setw(20) << 0.0 << setw(20) << sphere.height + 1.0 << endl;
	}
      if (i == sphere.n_int - 1)
	{
	  write << setw(20) << i << setw(20) << sphere.midpoints[i] << setw(20) << 0.0 << setw(20) << sphere.height - 1.0 << endl;
	}
      else
	{
	  write << setw(20) << i << setw(20) << sphere.midpoints[i] << setw(20) << sin(sphere.midpoints[i]) << setw(20) << sphere.height + cos(sphere.midpoints[i]) << endl;
	}
    }

  write.close();
}
