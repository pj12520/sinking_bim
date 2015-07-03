//File containing defintions for functions relating to input and output

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "inout.h"
#include "object.h"
#include "const.h"

using std::fstream;
using std::string;
using std::ofstream;
using std::ostringstream;
using std::setw;
using std::endl;

using math_const::PI;

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
void Out_sys(out_data data, double mdr, double bond, double viscos_rat)
{
  ofstream write;

  //  string file = "D=" + static_cast<ostringstream*>( &(ostringstream() << mdr) )->str() + "/Bo=" + static_cast<ostringstream*>( &(ostringstream() << bond) )->str() + "/viscos_rat=" + static_cast<ostringstream*>( &(ostringstream() << viscos_rat) )->str() + "/interf_config" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";
  string file = "interf_config" + static_cast<ostringstream*>( &(ostringstream() << data.it) )->str() + ".dat";
  write.open(file.c_str());

  write << setw(20) << "Interval" << "Radial" << setw(20) << "Vertical" << setw(20) << "Radial_vel" << setw(20) << "Vertical_vel" << endl;

  for (int i = 0; i < data.interf_rad.size(); i++)
    {
      write << setw(20) << i << setw(20) << data.interf_rad[i] << setw(20) << data.interf_vert[i] << endl;
    }

  write.close();

  //  file = "D=" + static_cast<ostringstream*>( &(ostringstream() << mdr) )->str() + "/Bo=" + static_cast<ostringstream*>( &(ostringstream() << bond) )->str() + "/viscos_rat=" + static_cast<ostringstream*>( &(ostringstream() << viscos_rat) )->str() + "/sphere_config" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";
  file = "sphere_config" + static_cast<ostringstream*>( &(ostringstream() << data.it) )->str() + ".dat";
  write.open(file.c_str());

  write << setw(20) << "Interval" << setw(20) << "Radial" << setw(20) << "Vertical" << endl;

  for (int i = 0; i < 100; i++)
    {
      if (i == 0)
	{
	  write << setw(20) << i << setw(20) << 0.0 << setw(20) << data.sphere_pos + 1.0 << endl;
	}
      if (i == 99)
	{
	  write << setw(20) << i << setw(20) << 0.0 << setw(20) << data.sphere_pos - 1.0 << endl;
	}
      else
	{
	  write << setw(20) << i << setw(20) << sin(PI * i / 100.0) << setw(20) << data.sphere_pos + cos(PI * i / 100.0) << endl;
	}
    }

  write.close();
}
