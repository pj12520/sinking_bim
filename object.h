//Header file containing structures that describe the sphere and interface and related functions

#ifndef __OBJECT_H__
#define __OBJECT_H__
#pragma once

#include <vector>

using std::vector;

//Structure to contain properties of an interval on the surface of the sphere
struct sphere_int
{
  //Theta, radial and vertical coordinates of points used in the Gauss-Legendre quadrature
  vector<double> theta; 
  vector<double> rad;
  vector<double> vert;

  double lower; //Lower bound of the interval
  double upper; //Upper bound of the interval
  double width; //Width of the interval
};

//Structure to contain properties of an interval on the interface
struct interf_int
{
  //Arclength, radial and vertical coordinates of points used in the Gauss-Legendre quadrature
  vector<double> arc;
  vector<double> rad;
  vector<double> vert;

  //Radial and vertical components of the normal vector at the points used in the Gauss-Legendre quadrature
  vector<double> norm_rad;
  vector<double> norm_vert;

  double lower; //Lower bound of the interval
  double upper; //Upper bound of the interval
  double width; //Width of the interval
};

//Structure to contain sphere properties
struct particle
{
  double height; //Height of centre of sphere above original plane of interface
  int n_int; //Number of intervals used to discretise surface of sphere
  vector<double> midpoints; //Values of theta that give the midpoints of the intervals
  vector<sphere_int> intervals; //Array of interval structures 
};

//Structure to contain interface properties
struct surf
{
  int n_int; //Number of intervals that discretise the interface
  vector<double> midpoints; //Values of arclength that give the midpoint of each interval
  vector<interf_int> intervals; //Array of interval structures
};


//Function to fill the properties of the particle structure
void Create_sphere(particle sphere, double height, int n_int);

//Function to fill the properties of the surf structure
void Create_interf(surf interf, int n_int, double max_arc);

//Function to find the midpoints in a set of n_int equally spaced intervals
void Find_midpoints(vector<double> midpoints, double start, double end, int n_int);

//Function to create the intervals that cover the surface of the sphere
void Create_sphere_int(vector<sphere_int>* intervals, int n_int, double height);

//Function to create the intervals that cover the interface
void Create_interf_int(vector<interf_int>* intervals, int n_int, double max_arc);

#endif /* OBJECT_H */

