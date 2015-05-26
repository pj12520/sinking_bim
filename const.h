//Header file containing constants kept within a namespace

#ifndef __CONST_H__
#define __CONST_H__
#pragma once

namespace math_const
{
  const double PI = 3.141592653589793; //To 15dp. cite Woan 2004 The Cambridge Handbook of Physics Formulas
}

namespace Gauss
{
  //4 pt-Gaussian quadrature integration points on the interval [-1:1] (Riley, Hobson and Bence 2006 page 1008)
  const double Gauss_pt1 = -0.8611363116;
  const double Gauss_pt2 = -0.3399810436;
  const double Gauss_pt3 = 0.3399810436;
  const double Gauss_pt4 = 0.8611363116;

  //4 pt-Gaussian quadrature weights (Riley, Hobson and Bence 2006 page 1008)
  const double Gauss_wt1 = 0.3478548451;
  const double Gauss_wt2 = 0.6521451549;
  const double Gauss_wt3 = 0.6521451549;
  const double Gauss_wt4 = 0.3478548451;

}
#endif /* CONST_H */
