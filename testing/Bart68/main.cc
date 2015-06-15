//Program file to use the model of Bart 1968 to find the velocity of sphere moving towards a flat fluid-fluid interface as a function of position and viscosity ratio

#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>

using std::vector;
using std::ofstream;
using std::setw;
using std::endl;

int main()
{
  //Define independent variables
  vector<double> viscos_rat(7);
  viscos_rat[0] = 0.001;
  viscos_rat[1] = 0.01;
  viscos_rat[2] = 0.1;
  viscos_rat[3] = 1.0;
  viscos_rat[4] = 10.0;
  viscos_rat[5] = 100.0;
  viscos_rat[6] = 1000.0;

  vector<double> pos(9);
  pos[0] = 2.0;
  pos[1] = 3.0;
  pos[2] = 4.0;
  pos[3] = 5.0;
  pos[4] = 6.0;
  pos[5] = 7.0;
  pos[6] = 8.0;
  pos[7] = 9.0;
  pos[8] = 10.0;

  //Define dependent variables
  vector<vector<double> > sols(pos.size());
  for (int i = 0; i < sols.size(); i++)
    {
      sols[i].resize(viscos_rat.size());
    }

  //Declare temporary variables to be used (as defined in Bart 1968)
  double alpha;

  double X;
  double Y;
  double V;
  double T;

  double sum;

  double factor;
  //Define number of terms used to approximate infinite series
  int n_terms = 10;

  //Objects to output data to two data files, one sorted by viscosity ratio and one by sphere position
  ofstream viscos_sort;
  ofstream pos_sort;

  //Loop over sphere positions
  for (int i = 0; i < pos.size(); i++)
    {
      alpha = acosh(pos[i]);

      //Loop over viscosity ratios
      for (int j = 0; j < viscos_rat.size(); j++)
	{
	  sum = 0.0;

	  //Loop over param used in infinite sum
	  for (int k = 0; k < n_terms; k++)
	    {
	      X = (2.0 * k + 1.0) * (2.0 * k + 1.0) * sinh(alpha) * sinh(alpha) + 4.0 * cosh((k + 0.5) * alpha) * cosh((k + 0.5) * alpha);

	      Y = 2.0 * sinh((2.0 * k + 1.0) * alpha) + (2.0 * k + 1.0) * sinh(2.0 * alpha);

	      V = 2.0 * sinh((2.0 * k + 1.0) * alpha) - (2.0 * k + 1.0) * sinh(2.0 * alpha);

	      T = 4.0 * sinh((k + 0.5) * alpha) * sinh((k + 0.5) * alpha) - (2.0 * k + 1.0) * (2.0 * k + 1.0) * sinh(alpha) * sinh(alpha);

	      sum += k * (k + 1.0) / ((2.0 * k - 1.0) * (2.0 * k + 3.0)) * ((X + viscos_rat[j] * Y) / (V + viscos_rat[j] * T) - 1.0);
	    }

	  factor = 4.0 * sinh(alpha) * sum / 3.0;

	  sols[i][j] = 1.0 / factor;
	}
    }

  //Open data files
  viscos_sort.open("viscos_sort.dat");
  pos_sort.open("pos_sort.dat");

  //Write headings to files
  viscos_sort << setw(20) << "Viscosity Ratio" << setw(20) << "Position" << setw(20) << "Velocity" << endl;
  pos_sort << setw(20) << "Position" << setw(20) << "Viscosity Ratio" << setw(20) << "Velocity" << endl;

  //Output viscosity sorted data
  for (int i = 0; i < viscos_rat.size(); i++)
    {
      for (int j = 0; j < pos.size(); j++)
	{
	  viscos_sort << setw(20) << viscos_rat[i] << setw(20) << pos[j] << setw(20) << sols[j][i] << endl;
	}
    }

  //Output position sorted data
  for (int i = 0; i < pos.size(); i++)
    {
      for (int j = 0; j < viscos_rat.size(); j++)
	{
	  pos_sort << setw(20) << pos[i] << setw(20) << viscos_rat[j] << setw(20) << sols[i][j] << endl;
	}
    }

  //Close data files
  viscos_sort.close();
  pos_sort.close();

  return 0;
}
