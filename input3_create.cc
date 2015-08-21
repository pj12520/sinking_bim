//Program to create input files for parameter sweep 3

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

  vector<double> mod_dens_rat_data(3);
  vector<double> bond_data(24);
  vector<double> viscos_rat_data(4);

  mod_dens_rat_data[0] = 10;
  mod_dens_rat_data[1] = 100;
  mod_dens_rat_data[2] = 1000;

  bond_data[0] = 0.002;
  bond_data[1] = 0.003;
  bond_data[2] = 0.004;
  bond_data[3] = 0.005;
  bond_data[4] = 0.006;
  bond_data[5] = 0.007;
  bond_data[6] = 0.008;
  bond_data[7] = 0.009;
  bond_data[8] = 0.02;
  bond_data[9] = 0.03;
  bond_data[10] = 0.04;
  bond_data[11] = 0.05;
  bond_data[12] = 0.06;
  bond_data[13] = 0.07;
  bond_data[14] = 0.08;
  bond_data[15] = 0.09;
  bond_data[16] = 0.2;
  bond_data[17] = 0.3;
  bond_data[18] = 0.4;
  bond_data[19] = 0.5;
  bond_data[20] = 0.6;
  bond_data[21] = 0.7;
  bond_data[22] = 0.8;
  bond_data[23] = 0.9;

  viscos_rat_data[0] = 0.001;
  viscos_rat_data[1] = 0.01;
  viscos_rat_data[2] = 0.1;
  viscos_rat_data[3] = 1;

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

  const char* data = "param_sweep3/";

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
