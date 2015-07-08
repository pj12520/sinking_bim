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

//Function to evaluate whether or not the sphere and interface have collided yet
int Break_Crit(vector<double>* arc, vector<double>* rad, vector<double>* vert, double sphere_pos)
{
  double separation;
  int break_criteria = 0;

  double min_step = (*arc)[1] - (*arc)[0];

  for (int i = 0; i < (*arc).size(); i++)
    {
      separation = Pythag((*rad)[i], (*vert)[i] - sphere_pos);
      if (separation - 1.0 < min_step)
	{
	  break_criteria = 1;
	  cout << "Sphere and interface collide" << endl;
	  break;
	}

      if (i != 0)
	{
	  if ((*rad)[i] < (*rad)[i - 1] && (*rad)[i] < min_step)
	    {
	      break_criteria = 2;
	      cout << "Tail snaps" << endl;
	      break;
	    }
	}
    }

  return break_criteria;
}
