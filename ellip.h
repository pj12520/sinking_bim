//Header file containing functions relating to the evaluation of complete elliptic integrals of the first and second kind

#ifndef __ELLIP_H__
#define __ELLIP_H__
#pragma once

//Function to evaluate the complete elliptic integral of the 1st kind (Abramovitch and Stegun 1964 page 591)
double Ellip1(double comp_param);

//Function to evaluate the variable part of the complete elliptic integral of the second kind
double Ellip2_var(double comp_param);

//Function to evaluate the complete elliptic integral of the 2nd kind (Abramovitch and Stegun 1964 page 591)
double Ellip2(double comp_param);

//Function to evaluate the regular part of the complete elliptic integral of the first kind
double Ellip1_reg(double comp_param);

//Function to evaluate the unbounded part of the complete elliptic integral of the first kind
double Ellip1_sing(double comp_param);

//Function to evaluate the complementary parameter used as the argument for the functions which compute the elliptic integrals
double Comp_param(double beta_2, double sum);

//Function definitions

#endif /* ELLIP_H */
