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

  vector<double> div_norm; //Divergence of the normal vector

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
  vector<double> mid_rad; //Values of the radial coordinate at the midpoint of each interval
  vector<double> mid_vert; //Values of the vertical coordinate at the midpoint of each interval

  vector<double> mid_norm_rad; //Radial component of the normal vector at the midpoint of each interval
  vector<double> mid_norm_vert; //Vertical component of the normal vector at the midpoint of each interval
  vector<double> mid_div_norm;

  vector<interf_int> intervals; //Array of interval structures
};


//Function to fill the properties of the particle structure
void Create_sphere(particle *sphere, double height, int n_int);

//Function to fill the properties of the surf structure
void Create_interf(surf *interf, int n_int, double max_arc);

//Function to find the midpoints in a set of n_int equally spaced intervals
void Find_midpoints(vector<double>* midpoints, double start, double end, int n_int);

//Function to create the intervals that cover the surface of the sphere
void Create_sphere_int(vector<sphere_int>* intervals, int n_int, double height);

//Function to create the intervals that cover the interface
void Create_interf_int(vector<interf_int>* intervals, int n_int, double max_arc);

//Function to create the abscissas used for Gauss-Legendre integration in an interval
void Abscissas(double* lower, double* upper, double max, int n_int, vector<double>* points, double* width, double half_width, int interval, int sphere_test); 

//Function to update the properaties of the interface
void Up_interf(surf *interf);

//Function to update the properties of the sphere
void Up_sphere(particle *sphere);

#endif /* OBJECT_H */

