//Program to create input files

#include <vector>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

using std::vector;
using std::ofstream;
using std::vector;
using std::endl;
using std::string;
using std::ostringstream;

int main()
{

  //Create arrays containing the different values of the dimensionless numbers over which we are sweeping

  vector<double> mod_dens_rat_data(4);
  vector<double> bond_data(7);
  vector<double> viscos_rat_data(6);

  mod_dens_rat_data[0] = 1.0;
  mod_dens_rat_data[1] = 10.0;
  mod_dens_rat_data[2] = 100.0;
  mod_dens_rat_data[3] = 1000.0;

  bond_data[0] = 0.001;
  bond_data[1] = 0.01;
  bond_data[2] = 0.1;
  bond_data[3] = 1.0;
  bond_data[4] = 10.0;
  bond_data[5] = 100.0;
  bond_data[6] = 1000.0;

  viscos_rat_data[0] = 0.001;
  viscos_rat_data[1] = 0.01;
  viscos_rat_data[2] = 0.1;
  viscos_rat_data[3] = 10;
  viscos_rat_data[4] = 100;
  viscos_rat_data[5] = 1000;

  int n_sphere = 100;
  int n_interf = 100;
  double trunc = 15.0;
  //  double t_step = 0.01;
  double height = 3.0;
  int max_it = 60000;

  //Create output directories 

  string D_dir_name;
  string Bo_dir_name;
  string viscos_dir_name;

  ostringstream D_convert;
  ostringstream Bo_convert;
  ostringstream viscos_convert;

  const char* up = "../";

  const char* data = "new_data/";

  chdir(data);

  for (int i = 0; i < mod_dens_rat_data.size(); i++)
    {
      D_convert << "D=" << mod_dens_rat_data[i];

      D_dir_name = D_convert.str();

      mkdir(D_dir_name.c_str(), S_IRWXU);

      for (int j = 0; j < bond_data.size(); j++)
	{
	  Bo_convert << "Bo=" << bond_data[j];

	  Bo_dir_name = Bo_convert.str();

	  chdir(D_dir_name.c_str());

	  mkdir(Bo_dir_name.c_str(), S_IRWXU);

	  for (int k = 0; k < viscos_rat_data.size(); k++)
	    {
	      viscos_convert << "viscos_rat=" << viscos_rat_data[k];

	      viscos_dir_name = viscos_convert.str();

	      chdir(Bo_dir_name.c_str());

	      mkdir(viscos_dir_name.c_str(), S_IRWXU);

	      chdir(viscos_dir_name.c_str());

	      ofstream fout;

	      fout.open("dimensionless_input.dat");

	      fout << bond_data[j] << endl;
	      fout << mod_dens_rat_data[i] << endl;
	      fout << viscos_rat_data[k] << endl;

	      fout << n_sphere << endl;
	      fout << n_interf << endl;
	      fout << trunc << endl;
	      //	      fout << t_step << endl;
	      fout << height << endl;
	      fout <<max_it << endl;

	      fout << "Input file containing the dimensionless numbers that characterise the system" << endl;

	      fout << "Bond Number" << endl;
	      fout << "Modified Density Ratio" << endl;
	      fout << "Viscosity Ratio" << endl;

	      fout << "Number of elements on sphere" << endl;
	      fout << "Number of elements on interface" << endl;
	      fout << "Truncation length of interface" << endl;
	      //	      fout << "Time step" << endl;
	      fout << "Initial height of sphere" << endl;
	      fout << "Maximum number of iterations" << endl;

	      fout.close();

	      chdir(up);

	      chdir(up);

	      viscos_convert.str(string());

	    }

	  chdir(up);

	  Bo_convert.str(string());

	}

      D_convert.str(string());
    }

  return 0;
}
