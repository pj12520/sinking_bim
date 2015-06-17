//File containing function defintions relating to iteration of the system forward in time

#include <vector>
#include <math.h>
#include <iostream> //Inlcuded for the purposes of debugging only

#include "geo.h"

using std::vector;

using std::cout; //Using for the purposes of debugging only
using std::endl; //Using for the purposes of debugging only

void Iterate(double n_int, vector<double>* unknown, vector<double>* arc, vector<double>* rad, vector<double>* vert, double *height, double t_step)
{
  //Extract the velocities from unknown vector

  vector<double> rad_vel(n_int);
  vector<double> vert_vel(n_int);

  double sphere_vel;

  for (int i = 0; i < n_int; i++)
    {
      if (i == 0)
	{
	  rad_vel[i] = 0.0;
	}
      else
	{
	  rad_vel[i] = (*unknown)[i - 1];
	}

     vert_vel[i] = (*unknown)[i + n_int - 1];
    }

  sphere_vel = (*unknown)[(*unknown).size() - 1];

  //Iterate the system forward in time
  for (int i = 0; i < n_int; i++)
    {
      (*rad)[i] += rad_vel[i] * t_step;
      (*vert)[i] += vert_vel[i] * t_step;

      if (i != 0)
	{
	  (*arc)[i] = (*arc)[i - 1] + Pythag((*rad)[i] - (*rad)[i - 1],(*vert)[i] - (*vert)[i - 1]); 
	}
    }

  *height += sphere_vel * t_step;
}
