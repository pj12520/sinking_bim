// Header file containing functions relating to iteration of the system in time

#ifndef __VELOCITY_H__
#define __VELOCITY_H__
#pragma once

#include <vector>

using std::vector;

//Function to iterate the system forward in time
void Iterate(double n_int, vector<double>* unknown, vector<double>* arc, vector<double>* rad, vector<double>* vert, double *height, double t_step);

//Function to evaluate whether or not the sphere and interface have collided yet
int Break_Crit(vector<double>* arc, vector<double>* rad, vector<double>* vert, double sphere_pos);

#endif /* VELOCITY_H */
