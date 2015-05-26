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

  sphere.midpoints.resize(n_int); 
  for (int i = 0; i < n_int; i++)
    {
      sphere.midpoints[i] = i * PI / (n_int - 1);
    }
}

//Function to fill the properties of the surf structure
void Create_interf(surf interf, int n_int, double max_arc)
{
  interf.n_int = n_int;

  interf.midpoints.resize(n_int); 
  for (int i = 0; i < n_int; i++)
    {
      interf.midpoints[i] = i * max_arc / (n_int - 1);
    }
}
