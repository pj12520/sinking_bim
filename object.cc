//File containing defintions for functions relating to the structures for the sphere and interface

#include <vector>
#include <math.h>
#include <new>

#include <iostream> //Currently only included for debugging purposes. Can be removed when program is functional

#include "object.h"
#include "const.h"

using std::vector;
using math_const::PI;

using std::cout; //Currently only used for debugging purposes. Can be removed when program is functional
using std::endl; //Currently only used for debugging purposes. Can be removed when program is functional

using namespace Gauss;

//Function to fill the properties of the particle structure
void Create_sphere(particle *sphere, double height, int n_int)
{
  (*sphere).height = height;
  (*sphere).n_int = n_int;

  Find_midpoints(&(*sphere).midpoints, 0.0, PI, n_int);

  Create_sphere_int(&(*sphere).intervals, n_int, height);
}

//Function to fill the properties of the surf structure
void Create_interf(surf *interf, int n_int, double max_arc)
{
  (*interf).n_int = n_int;

  Find_midpoints(&(*interf).midpoints, 0.0, max_arc, n_int);

  Create_interf_int(&(*interf).intervals, n_int, max_arc);
}

//Function to find the midpoints in a set of n_int equally spaced intervals
void Find_midpoints(vector<double>* midpoints, double start, double end, int n_int)
{
  (*midpoints).resize(n_int); 
  for (int i = 0; i < n_int; i++)
    {
      (*midpoints)[i] = start + i * (end - start) / (n_int - 1);
    }

}

//Function to create the intervals that cover the surface of the sphere
void Create_sphere_int(vector<sphere_int>* intervals, int n_int, double height)
{
  double half_width = PI / (2.0 * (n_int - 1)); //This is half the width of the intermediate intervals. 

  for (int i = 0; i < n_int; i++)
    {
      //Create the abscissas for the Gauss-Legendre integration in each interval
      Abscissas(&(*intervals)[i].lower, &(*intervals)[i].upper, PI, n_int, &(*intervals)[i].theta, &(*intervals)[i].width, half_width, i);

      //Set radial and vertical components of the integration points in each interval
      (*intervals)[i].rad.resize(4);
      (*intervals)[i].vert.resize(4);

      for (int j = 0; j < 4; j++)
	{
	  (*intervals)[i].rad[j] = sin((*intervals)[i].theta[j]);
	  (*intervals)[i].vert[j] = height + cos((*intervals)[i].theta[j]);
	}
    }
}

//Function to create the intervals that cover the interface
void Create_interf_int(vector<interf_int>* intervals, int n_int, double max_arc)
{
  double half_width = max_arc / (2.0 * (n_int - 1)); //This is half the width of the intermediate intervals. 

  for (int i = 0; i < n_int; i++)
    {
      //Create the abscissas for the Gauss-Legendre integration in each interval
      Abscissas(&(*intervals)[i].lower, &(*intervals)[i].upper, max_arc, n_int, &(*intervals)[i].arc, &(*intervals)[i].width, half_width, i);

      //Set radial and vertical components of the integration points in each interval and components of normal vectors
      (*intervals)[i].rad.resize(4);
      (*intervals)[i].vert.resize(4);

      (*intervals)[i].norm_rad.resize(4);
      (*intervals)[i].norm_vert.resize(4);

      for (int j = 0; j < 4; j++)
	{
	  (*intervals)[i].rad[j] = (*intervals)[i].arc[j];
	  (*intervals)[i].vert[j] = 0.0;

	  (*intervals)[i].norm_rad[j] = 0.0;
	  (*intervals)[i].norm_vert[j] = 1.0;
	}
    }

}

//Function to create the abscissas used for Gauss-Legendre integration in an interval
void Abscissas(double* lower, double* upper, double max, int n_int, vector<double>* points, double* width, double half_width, int interval)
{
  vector<double> Gauss_int_pts(4); //Vector to store integration points used for 4-pt Gaussian quadrature

  Gauss_int_pts[0] = Gauss_pt1;
  Gauss_int_pts[1] = Gauss_pt2;
  Gauss_int_pts[2] = Gauss_pt3;
  Gauss_int_pts[3] = Gauss_pt4;

  //Set upper and lower bounds of each interval
  if (interval == 0)
    {
      *lower = 0.0;

      *upper = half_width;
    }
  else if (interval == n_int - 1)
    {
      *lower = (2.0 * n_int - 3.0) * half_width;
      *upper = max;
    }
  else
    {
      *lower = (2.0 * interval - 1) * half_width;
      *upper = (2.0 * interval + 1) * half_width;
    }

  *width = *upper - *lower;
  //Set abscissas for each interval
  (*points).resize(4);
  for (int j = 0; j < 4; j++)
    {
      (*points)[j] = (*upper + *lower + Gauss_int_pts[j] * (*upper - *lower)) / 2.0;
    }
}
