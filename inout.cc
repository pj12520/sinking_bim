//File containing defintions for functions relating to input and output

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

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

  read >> input->viscos_rat; //Viscosity Ratio
  read >> input->bond; //Bond number
  read >> input->mdr; //Modified Density Ratio

  read >> input->n_sphere; //Number of elements on sphere's surface
  read >> input->n_int; //Number of elements on interface
  read >> input->max_arc; //Truncation length of interface
  read >> input->t_step; //Initial time step
  read >> input->init_height; //Initial height of sphere
  read >> input ->max_it; //Maximum number of iterations
  read.close();

}

//Function to output the state of the system
void Out_sys(int it, surf interf)
{
  ofstream write;

  string file = "results/config_" + static_cast<ostringstream*>( &(ostringstream() << it) )->str() + ".dat";

  write.open(file.c_str());

  write << setw(20) << "Interval" << setw(20) << "Arc" << setw(20) << "Radial" << setw(20) << "Vertical" << endl;

  for (int i = 0; i < interf.n_int; i++)
    {
      write << setw(20) << i << setw(20) << interf.midpoints[i] << setw(20) << interf.mid_rad[i] << setw(20) << interf.mid_vert[i] << endl;
    }

  write.close();
}
