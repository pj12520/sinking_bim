//File containing defintions for functions relating to the structures for the sphere and interface

#include <vector>

#include "object.h"
#include "const.h"

using std::vector;
using math_const::PI;

//Function to fill the properties of the particle structure
void Create_sphere(particle sphere, double height, int n_int)
{
  sphere.height = height;
  sphere.n_int = n_int;

  Find_midpoints(sphere.midpoints, 0.0, PI, n_int);
}

//Function to fill the properties of the surf structure
void Create_interf(surf interf, int n_int, double max_arc)
{
  interf.n_int = n_int;

  Find_midpoints(interf.midpoints, 0.0, max_arc, n_int);
}

//Function to find the midpoints in a set of n_int equally spaced intervals
void Find_midpoints(vector<double> midpoints, double start, double end, int n_int)
{
  midpoints.resize(n_int); 
  for (int i = 0; i < n_int; i++)
    {
      midpoints[i] = start + i * (end - start) / (n_int - 1);
    }

}
