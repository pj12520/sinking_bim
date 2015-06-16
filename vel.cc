//File containing function defintions relating to iteration of the system forward in time

#include <vector>
#include <math.h>

#include "geo.h"

using std::vector;

void Iterate(double n_int, vector<double>* unknown, vector<double>* arc, vector<double>* rad, vector<double>* vert, double *height, double t_step)
{
  //Extract the velocities from unknown vector

  vector<double> rad_vel(n_int - 1);
  vector<double> vert_vel(n_int);

  double sphere_vel;

  for (int i = 0; i < n_int - 1; i++)
    {
      rad_vel[i] = (*unknown)[i];
      vert_vel[i] = (*unknown)[i + n_int - 1];
    }
  vert_vel[n_int - 1] = (*unknown)[2* n_int - 2];

  sphere_vel = (*unknown)[(*unknown).size() - 1];

  //Iterate the system forward in time
  for (int i = 0; i < n_int; i++)
    {
      (*vert)[i] += vert_vel[i] * t_step;

      if (i != 0)
	{
	  (*rad)[i] += rad_vel[i - 1] * t_step;

	  //	  (*arc)[i] = sqrt(((*rad)[i] - (*rad)[i - 1]) * ((*rad)[i] - (*rad)[i - 1]) + ((*vert)[i] - (*vert)[i - 1]) * ((*vert)[i] - (*vert)[i - 1]));
	  (*arc)[i] = Pythag((*rad)[i] - (*rad)[i - 1],(*vert)[i] - (*vert)[i - 1]); 
	}
    }

  *height += sphere_vel * t_step;
}
