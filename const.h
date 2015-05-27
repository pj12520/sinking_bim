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

//Constants used in the polynomial approximations of the complete elliptic integrals (Abramovitch and Stegun 1956 pages 591 - 592)
namespace ellip_poly
{
  const double ellip1_a0 = 1.38629436112;
  const double ellip1_a1 = 0.09666344259;
  const double ellip1_a2 = 0.03590092383;
  const double ellip1_a3 = 0.03742563713;
  const double ellip1_a4 = 0.01451196212;

  const double ellip1_b0 = 0.5;
  const double ellip1_b1 = 0.12498593597;
  const double ellip1_b2 = 0.06880248576;
  const double ellip1_b3 = 0.03328355346;
  const double ellip1_b4 = 0.00441787012;

  const double ellip2_a1 = 0.44325141463;
  const double ellip2_a2 = 0.06260601220;
  const double ellip2_a3 = 0.04757383546;
  const double ellip2_a4 = 0.01736506451;

  const double ellip2_b1 = 0.24998368310;
  const double ellip2_b2 = 0.09200180037;
  const double ellip2_b3 = 0.04069697526;
  const double ellip2_b4 = 0.00526449639;
}
#endif /* CONST_H */
