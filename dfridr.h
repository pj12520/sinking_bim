//Header file from Numerical Recipes page 231 containing a function to evaluate the derivative of a function

//IMPORTANT: When calling any of these functions, it is imperitive that x - h and x + h lie within the interpolation range
#ifndef __DFRIDR_H__
#define __DFRIDR_H__
#pragma once

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <iostream>

#include "interp_1d.h"

using std::vector;
using std::numeric_limits;
using std::max;
using std::cout;
using std::endl;

double dfridr(double (*func)(double), const double x, const double h, double &err);

double sec_dfridr(Spline_interp spline, const double x, const double h, double &err);

//Modified version of dfridr which reads in a Spline_interp instead of a function and then interpolates the interp function of the Spline_interp.
double dfridr_interp(Spline_interp spline, const double x, const double h, double &err);

//Function to evaluate the second derivative of a spline using the array of second derivatives stored already
double sec_deriv(Spline_interp spline, const double point, vector<double> data); 

//Function to evaluate the 1st derivative of a spline using the array of second derivatives already stored
double deriv(Spline_interp spline, const double point, vector<double>* data1, vector<double>* data2);

double dfridr(double (*func)(double), const double x, const double h, double &err)
//Returns the derivativre of a function func at a point x by Riddlers' method of polynomial extrapolation. The value h is input as an estimated initial stepsize; it need not be small, but rather should be an increment in x over which func changes substantially. An estimate of the error in the derivative is returned as err.
{
  const int ntab = 10; //Set maximum size of tableau
  const double con = 1.4, con2 = con*con; //Stepsize decreased by con at each iteration
  const double big = numeric_limits<double>::max();
  const double safe = 2.0; //Return when error is safe worse than the best so far
  int i, j;
  double errt, fac, hh, ans;
  vector<vector<double> > a;

  a.resize(ntab);
  for (int k = 0; k < ntab; k++)
    {
      a[k].resize(ntab);
    }

  if (h == 0.0) 
    {
      cout << "h must be nonzero in dfridr" << endl;
      throw("h must be nonzero in dfridr.");
    }

  hh = h;
  a[0][0] = (func(x + hh) - func(x - hh)) / (2.0 * hh);
  err = big;

  for (i = 1; i <ntab; i++)
    {
      //Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
      hh /= con;
      a[0][i] = (func(x + hh) - func(x - hh)) / (2.0 * hh); //Try new, smaller stepsize.
      fac = con2;

      for (j = 1; j <= i; j++) //Compute extrapolations of various orders, requiring no new function evaluations.
	{
	  a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
	  fac = con2 * fac;
	  errt = max(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]));
	  //The error strategy is to compare each new extrapolation to one order lower, both at the present stepsize and the previous one.

	  if (errt <= err) //If the error is decreased, save the improved answer
	  {
	    err = errt;
	    ans = a[j][i];
	  }
	}

      if (abs(a[i][i] - a[i - 1][i - 1]) >= safe*err) break;
      //If higher order is worse by a significant factor safe, then quit early.
    }

  if (fabs(ans) < 1.0e-14)
    {
      ans = 0;
    } 

  return ans;
}

double sec_dfridr(Spline_interp spline, const double x, const double h, double &err)
{
  const int ntab = 10; //Set maximum size of tableau
  const double con = 1.4, con2 = con*con; //Stepsize decreased by con at each iteration
  const double big = numeric_limits<double>::max();
  const double safe = 2.0; //Return when error is safe worse than the best so far
  int i, j;
  double errt, fac, hh, ans;
  vector<vector<double> > a;

  double deriv1_plus;
  double deriv1_minus;

  double error;
  a.resize(ntab);
  for (int k = 0; k < ntab; k++)
    {
      a[k].resize(ntab);
    }

  if (h == 0.0) 
    {
      cout << "h must be nonzero in sec_dfridr" << endl;
      throw("h must be nonzero in dfridr.");
    }
  hh = h;

  deriv1_plus = dfridr_interp(spline, x + hh, h, error);

  //  if (x - hh < 0.0)
  //{
  //  deriv1_minus = dfridr_interp(spline, -(x - hh), h, error);
  //}
  //  else
  //{
  deriv1_minus = dfridr_interp(spline, x - hh, h, error);
  //  }

  a[0][0] = (deriv1_plus - deriv1_minus) / (2.0 * hh);
  err = big;

  for (i = 1; i <ntab; i++)
    {
      //Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
      hh /= con;

      deriv1_plus = dfridr_interp(spline, x + hh, h, error);

      //  if (x - hh < 0.0)
      //{
      //  deriv1_minus = dfridr_interp(spline, -(x - hh), h, error);
      //}
      //      else
      //{
      deriv1_minus = dfridr_interp(spline, x - hh, h, error);
      //	}

      a[0][i] = (deriv1_plus - deriv1_minus) / (2.0 * hh); //Try new, smaller stepsize.
      fac = con2;

      for (j = 1; j <= i; j++) //Compute extrapolations of various orders, requiring no new function evaluations.
	{
	  a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
	  fac = con2 * fac;
	  errt = max(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]));
	  //The error strategy is to compare each new extrapolation to one order lower, both at the present stepsize and the previous one.

	  if (errt <= err) //If the error is decreased, save the improved answer
	  {
	    err = errt;
	    ans = a[j][i];
	  }
	}

      if (abs(a[i][i] - a[i - 1][i - 1]) >= safe*err) break;
      //If higher order is worse by a significant factor safe, then quit early.
    }

  if (fabs(ans) < 1.0e-14)
    {
      ans = 0;
    } 

  return ans;
}

  double dfridr_interp(Spline_interp spline, const double x, const double h, double &err)
//Returns the derivativre of a function func at a point x by Riddlers' method of polynomial extrapolation. The value h is input as an estimated initial stepsize; it need not be small, but rather should be an increment in x over which func changes substantially. An estimate of the error in the derivative is returned as err.
{
  const int ntab = 10; //Set maximum size of tableau
  const double con = 1.4, con2 = con*con; //Stepsize decreased by con at each iteration
  const double big = numeric_limits<double>::max();
  const double safe = 2.0; //Return when error is safe worse than the best so far
  int i, j;
  double errt, fac, hh, ans;
  vector<vector<double> > a;

  double func_plus;
  double func_minus;

  a.resize(ntab);
  for (int k = 0; k < ntab; k++)
    {
      a[k].resize(ntab);
    }

  if (h == 0.0) 
    {
      cout << "Error in dfridr_interp" << endl;
      throw("h must be nonzero in dfridr.");
    }
  hh = h;

  //  if (x - hh < 0)
  //{
  //  cout << "Calculate function at radial coordinate less than 0" << endl;
  //}
  //  else if (x + hh > 5)
  //{
  //  cout << "Calculate function at radial coordinate greater than 5" << endl;
  //}
  func_plus = spline.interp(x + hh);

  //  if (x - hh < 0)
  //{
  //  func_minus = spline.interp(-(x - hh));
  //}
  //  else
  //{
  func_minus = spline.interp(x - hh);
  //  }

  a[0][0] = (func_plus - func_minus) / (2.0 * hh);
  err = big;

  for (i = 1; i <ntab; i++)
    {
      //Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of extrapolation.
      hh /= con;

      func_plus = spline.interp(x + hh);

      //  if (x - hh < 0)
      //{
      //  func_minus = spline.interp(-(x - hh));
      //}
      //      else
      //{
      func_minus = spline.interp(x - hh);
      //	}

      a[0][i] = (func_plus - func_minus) / (2.0 * hh); //Try new, smaller stepsize.
      fac = con2;

      for (j = 1; j <= i; j++) //Compute extrapolations of various orders, requiring no new function evaluations.
	{
	  a[j][i] = (a[j - 1][i] * fac - a[j - 1][i - 1]) / (fac - 1.0);
	  fac = con2 * fac;
	  errt = max(abs(a[j][i] - a[j - 1][i]), abs(a[j][i] - a[j - 1][i - 1]));
	  //The error strategy is to compare each new extrapolation to one order lower, both at the present stepsize and the previous one.

	  if (errt <= err) //If the error is decreased, save the improved answer
	  {
	    err = errt;
	    ans = a[j][i];
	  }
	}

      if (abs(a[i][i] - a[i - 1][i - 1]) >= safe*err) break;
      //If higher order is worse by a significant factor safe, then quit early.
    }

  if (fabs(ans) < 1.0e-14)
    {
      ans = 0;
    } 

  return ans;
}
 
//Function to evaluate the second derivative of a spline using the array of second derivatives stored already
double sec_deriv(Spline_interp spline, const double point, vector<double>* data) 
{
  double upper_deriv2;
  double lower_deriv2;

  double lower_point;
  double upper_point;

  double deriv2;

  double A;
  double B;

  for (int i = 0; i < (*data).size(); i++)
    {
      if (point == (*data)[i])
	{
	  deriv2 = spline.y2[i];
	}
      else if (point > (*data)[i])
	{
	  if (point < (*data)[i + 1])
	    {
	      lower_point = (*data)[i];
	      upper_point = (*data)[i + 1];

	      lower_deriv2 = spline.y2[i];
	      upper_deriv2 = spline.y2[i + 1];

	      A = (upper_point - point) / (upper_point - lower_point);
	      B = 1 - A;
	      
	      //	      deriv2 = lower_deriv2 + (upper_deriv2 - lower_deriv2) * (point - lower_point) / (upper_point - lower_point);
	      deriv2 = A * lower_deriv2 + B * upper_deriv2;
	    }
	}
    }

  return deriv2;
} 

//Function to evaluate the 1st derivative of a spline using the array of second derivatives already stored
double deriv(Spline_interp spline, const double point, vector<double>* data1, vector<double>* data2)
{
  double upper_deriv2;
  double lower_deriv2;

  double lower_point;
  double upper_point;

  double lower_y;
  double upper_y;

  double A;
  double B;

  double deriv1 = 0;
  for (int i = 0; i < (*data1).size(); i++)
    {
      if (point > (*data1)[i])
	{
	  if (point <= (*data1)[i + 1])
	    {
	      lower_point = (*data1)[i];
	      upper_point = (*data1)[i + 1];

	      lower_y = (*data2)[i];
	      upper_y = (*data2)[i + 1];

	      lower_deriv2 = spline.y2[i];
	      upper_deriv2 = spline.y2[i + 1];

	      A = (upper_point - point) / (upper_point - lower_point);
	      B = 1 - A;

	      deriv1 = (upper_y - lower_y) / (upper_point - lower_point) - (3.0 * A * A - 1.0) * (upper_point - lower_point) * lower_deriv2 / 6.0 + (3.0 * B * B - 1.0) * (upper_point - lower_point) * upper_deriv2 / 6.0; //From Numerical recipes equation 3.3.5 page 121
	    }
	}
    }
  return deriv1;
}
#endif /* DFRIDR_H */
