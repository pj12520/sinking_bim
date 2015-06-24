// Header file containing functions relating to iteration of the system in time

#ifndef __VELOCITY_H__
#define __VELOCITY_H__
#pragma once

#include <vector>

using std::vector;

//Function to iterate the system forward in time
void Iterate(double n_int, vector<double>* unknown, vector<double>* arc, vector<double>* rad, vector<double>* vert, double *height, double t_step, double *trunc_arc, double trunc_rad, double trunc_vert);

#endif /* VELOCITY_H */
