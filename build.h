//Header file containing functions relating to the construction of the linear system

#ifndef __BUILD_H__
#define __BUILD_H__
#pragma once

//Function to create the linear system
void Build(vector<vector<double> >* matrix, vector<double>* vector, particle sphere, surf interf, double viscos_rat, double bond, double mdr);

#endif /* BUILD_H */
