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
  int n_int;
  vector<double> midpoints;
};

//Structure to contain interface properties
struct surf
{
  int n_int;
  vector<double> midpoints;
};

//Function to fill the properties of the particle structure
void Create_sphere(particle sphere, double height, int n_int);

//Function to fill the properties of the surf structure
void Create_interf(surf interf, int n_int, double max_arc);

//Function to find the midpoints in a set of n_int equally spaced intervals
void Find_midpoints(vector<double> midpoints, double start, double end, int n_int);

#endif /* OBJECT_H */

