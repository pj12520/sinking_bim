//File containing defintions for functions relating to the structures for the sphere and interface

#include <vector>

#include "object.h"
#include "const.h"

using std::vector;
using math_const::PI;

void Create_sphere(particle sphere, double height, double n_int)
{
  sphere.height = height;
  sphere.n_int = n_int;

  sphere.midpoints.resize(n_int); 
  for (int i = 0; i < n_int; i++)
    {
      sphere.midpoints[i] = i * PI / (n_int - 1);
    }
}
