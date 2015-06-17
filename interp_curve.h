//Header file containing function to find derivative at end point of spline. Modified from Press et al page 149

#ifndef __INTERP_CURVE_H__
#define __INTERP_CURVE_H__
#pragma once

//Function declarations

//Function to find derivative at end-points of spline
double fprime(double *x, double *y, int pm);


//Function definitions

//Function to find derivative at end-points of spline
inline double fprime(double *x, double *y, int pm)
{
  //Utility for estimating the derivatives at the endpoints. x and y point to the abscissa and ordinate of the endpoint. If pm is +1, points to the right will be used (left endpoint); it if is -1, points to the left will be used (rigth endpoint).

  double s1 = x[0] - x[pm * 1];
  double s2 = x[0] - x[pm * 2];
  double s3 = x[0] - x[pm * 3];

  double s12 = s1 - s2;
  double s13 = s1 - s3;
  double s23 = s2 - s3;

  double deriv = -(s1 * s2 / (s13 * s23 *s3)) * y[pm * 3] + (s1 * s3 / (s12 * s2 * s23)) * y[pm * 2] - (s2 * s3 / (s1 * s12 * s13)) * y[pm * 1] + (1.0 / s1 + 1.0 / s2 + 1.0 / s3) * y[0];

  return deriv;
}
#endif /* __INTERP_CURVE_H */
