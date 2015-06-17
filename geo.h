//Header file containing function declarations relating to geometry

#ifndef __GEO_H__
#define __GEO_H__
#pragma once

#include <vector>

#include "interp_1d.h"

//Function to calculate hypotenuse of a right angled triangle given the length of the other sides.
double Pythag(double side1, double side2);

//Function to calculate the components of the normal vector and it's divergence at a point along the interface
void Normal(Spline_interp rad, Spline_interp height, double arc, double init_step, double *norm_rad, double *norm_vert, double *div_norm, double rad_coord);

//Function to calculate divergence of the normal of the interface at given point along it
double Div_norm(Spline_interp rad, Spline_interp height, double arc, double init_step); 

/*
//Function to calculate a tangent vector to a line at each point along it
void Tangent(vector<double>* rad, vector<double>* height, vector<double>* tangent_rad, vector<double>* tangent_height, int n_eval);
*/

//Function to calculate the radial component of the normal vector at a point along it
double Normal_rad(Spline_interp rad, Spline_interp height, double arc, double init_step);

//Function to calculate the vertical component of the normal vector at a point along it
double Normal_height(Spline_interp rad, Spline_interp height, double arc, double init_step);

/*
//Function to rotate a 2D vector in an anticlockwise direction by an angle theta
void Rotate(vector<double>* init_vector, vector<double>* final_vector, double theta);
*/



#endif /* GEO_H */
