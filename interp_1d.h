//Header file containing cubic polynomial interpolation routines used. These are from Numerical Recipes 3rd edition Press et al. Pages 122-124

#ifndef __INTERP_1D_H__
#define __INTERP_1D_H__
#pragma once

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using std::vector;
using std::min;
using std::max;
using std::cout;
using std::endl;
 
struct Base_interp
//Abstract base class used by all interpolation routines in this chapter. Only the routine interp is called directly by the user.
	   {
	     int n, mm, jsav, cor, dj;
	     const double *xx, *yy;
	     Base_interp(const vector<double> &x, const double *y, int m)
	     //Constructor: Set up for interpolating on a table of x's and y's of length m. Normally called by a derived class, not by the user.
	     : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
	       dj = min(1, (int)pow((double)n,0.25));
	     }

	     double interp(double x) {
	       //Given a value x, return an interpolated value using data pointed to by xx and yy.
	       int jlo = cor ? hunt(x) : locate(x);
	       return rawinterp(jlo,x);
	     }

	     int locate(const double x); //See definitions below
	     int hunt(const double x);

	     double virtual rawinterp(int jlo, double x) = 0;
	     //Derived classes provide this as the actual interpolation method.
};

inline int Base_interp::locate(const double x)
		//Given a value x, return a value j such that x is (insofar as possible) centered in the subrange xx[i...j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
{
  int ju, jm, j1;

  if (n < 2 || mm < 2 || mm > n) throw("locate size error");

  bool ascend = (xx[n - 1] >= xx[0]); //True if ascending order of table, false otherwise.

  j1 = 0; //Initialise lower and upper limits
  ju = n - 1;

  while (ju - j1 > 1) { //If we are not yet done, 
    jm = (ju + j1) >>1; //compute a mid-point,

    if(x >= xx[jm] == ascend)
      {
      j1 = jm; //and replace either the lower limit,
      }
    else
      {
      ju = jm; //or the upper limit as appropiate.
      }
  } //Repeat until the test condition is satisfied

  cor = abs(j1 - jsav) > dj ? 0 : 1;
  jsav = j1;

  return max(0, min(n - mm, j1 - ((mm - 2) >> 1)));
}

inline int Base_interp::hunt(const double x)
		//Given a value x, return a value j such that x is (insofar as possible) centred in the subrange xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either increasing or decreasing. The returned value is not less than 0, nor greater than n - 1.
{
  int j1 = jsav, jm, ju, inc = 1;

  if (n < 2 || mm < 2 || mm > n) throw ("hunt size error");

  bool ascnd = (xx[n - 1] > xx[0]); //True is ascending order of table, false otherwise. 

  if (j1 < 0 || j1 > n - 1) { //Input guess not useful. Go immeadiately to bisection.
    j1 = 0;
    ju = n - 1;
  } else {
    if (x >= xx[j1] == ascnd) { //Hunt up:
      for (;;) {
	ju = j1 + inc;
	if (ju >= n - 1) {ju = n - 1; break;} //Off end of table.
	else if (x < xx[ju] == ascnd) break; //Found bracket.
	else { //Not done, so double the increment and try again
	  j1 = ju;
	  inc += inc;
	}
      }
    } else { //Hunt down:
      ju = j1;
      for(;;) {
	j1 = j1 - inc;
	if(j1 <= 0) {j1 = 0; break;} //Off end of table
	else if (x >= xx[j1] == ascnd) break; //Found bracket.
	else { //Not done so double the increment and try again.
	  ju = j1;
	  inc += inc;
	}
      }
    }
  }
  while (ju - j1 > 1) { //Hunt is done so begin the final bisection phase:
    jm = (ju + j1) >> 1;
    if (x >= xx[jm] == ascnd)
      j1 = jm;
    else
      ju = jm;
  }
  cor = abs(j1 - jsav) > dj ? 0 : 1;
  jsav = j1;
  return max(0, min(n - mm, j1 - ((mm - 2) >> 1)));
}

struct Spline_interp : Base_interp
//Cubic spline interpolation object. Construct with x and y vectors and (optionally) values of the 1st derivative at the end points, then call interp for interpolated values.
{
  vector<double> y2;

  Spline_interp(const vector<double> &xv, const vector<double> &yv, double yp1 = 1.e99, double ypn = 1.e99)
    : Base_interp(xv, &yv[0], 2), y2(xv.size())
    {sety2(&xv[0], &yv[0], yp1, ypn);}

  //  Spline_interp(const vector<double> &xv, const double *yv, double yp1 = 1.e99, double ypn = 1.e99)
  //: Base_interp(xv, yv, 2), y2(xv.size())
  //{sety2(&xv[0], yv, yp1, ypn);}

  void sety2(const double *xv, const double *yv, double yp1, double ypn);
  double rawinterp(int j1, double xv);
};

inline void Spline_interp::sety2(const double *xv, const double *yv, double yp1, double ypn)

		   //This routine stores an array y2[0..n-1] with second derivatives of the interpolating function at the tabulated points pointed to by xv, using function values pointed to by yv. If yp1 and/or ypn are equal to1 X 10^99 or larger, the routine is signalled to set the corresponding boundary condition for a natural spline, with zero second derivatives on that boundary; otherwise, they are the values of the first derivatives at the en  dpoints.

{
  int i, k;
  double p, qn, sig, un;

  int n = y2.size();

  vector<double> u(n - 1);

  if (yp1 > 0.99e99) //The lower boundary condition is set either to be "natural" or else to have a specified first derivative
    y2[0] = u[0] = 0.0;
  else 
    {
      y2[0] = -0.5;
      u[0] = (3.0/(xv[1] - xv[0]))*((yv[1] - yv[0])/(xv[1] - xv[0]) - yp1);
    }

  for (i = 1; i < n - 1; i++) //This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.
    {
      sig = (xv[i] - xv[i - 1])/(xv[i + 1] - xv[i - 1]);
      p = sig * y2[i - 1] + 2.0;
      y2[i] = (sig - 1.0)/p;
      u[i] = (yv[i + 1] - yv[i])/(xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1])/(xv[i] - xv[i - 1]);
      u[i] = (6.0 * u[i]/(xv[i + 1] - xv[i - 1]) - sig*u[i - 1])/p;
    }

  if (ypn > 0.99e99) //The upper boundary condition is set either to be "natural" or else to have a specified first derivative
    {
      qn = un = 0.0;
    }
  else
    {
      qn = 0.5;
      un = (3.0/(xv[n - 1] - xv[n - 2])) * (ypn - (yv[n - 1] - yv[n - 2])/(xv[n - 1] - xv[n - 2]));
    }

  y2[n - 1] = (un - qn*u[n - 2])/(qn*y2[n - 2] + 1.0);

  for (k = n - 2; k >=0; k--) //This is the backsubstitution toop of the tridiagonal algorithm
    {
      y2[k] = y2[k] * y2[k + 1] + u[k];
    }
}

inline double Spline_interp::rawinterp(int j1, double x) 
//Given a value x, and using pointers to data xx and yy, and the stored vector of second derivatives y2, this routine returns the cubic spline interpolated value y.
{
  int klo = j1, khi = j1 + 1;
  double y, h, b, a;

  h = xx[khi] - xx[klo];

  if (h == 0.0)
    {
      cout << "Error in spline_interp" << endl;
      throw("Bad input to routine splint"); //The xa's must be distinct
    }

  a = (xx[khi] - x)/h;

  b = (x - xx[klo])/h; //Cubic spline polynomial is now evaluated

  y = a * yy[klo] + b * yy[khi] + ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * (h * h)/6.0;

  return y;
}


#endif /* INTERP_1D_H */
 
