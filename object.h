//Header file containing structures that describe the sphere and interface and related functions

#ifndef __OBJECT_H__
#define __OBJECT_H__
#pragma once

#include <vector>

using std::vector;

//Structure to contain sphere properties
struct particle
{
  double height;
  double n_int;
  vector<double> midpoints;
};

//Function to fill the properties of the particle structure
void Create_sphere(particle sphere, double height, double n_int);
#endif /* OBJECT_H */
