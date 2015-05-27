//Header file containing functions which are used to test aspects of the model

#ifndef __TESTING_H__
#define __TESTING_H__
#pragma once

//Function to test that the dimensionless input file is read correctly
void In_test(dimless_in input);

//Function to test that the sphere and interface are produced correctly
void Config(particle sphere, surf interf);

#endif /* TESTING_H */
