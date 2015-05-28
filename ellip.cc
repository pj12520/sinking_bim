//File containing definitions for functions related to evaluation of complete elliptic integrals of the first and second kind

#include <math.h>
#include <vector>

#include "const.h"

using std::vector;

using namespace ellip_poly;

//Function to evaluate the complete elliptic integral of the 1st kind (Abramovitch and Stegun 1964 page 591)
double Ellip1(double comp_param)
{
  vector<double> a(5);
  vector<double> b(5);

  double log_factor = log(comp_param);

  a[0] = ellip1_a0;
  a[1] = ellip1_a1;
  a[2] = ellip1_a2;
  a[3] = ellip1_a3;
  a[4] = ellip1_a4;

  b[0] = ellip1_b0;
  b[1] = ellip1_b1;
  b[2] = ellip1_b2;
  b[3] = ellip1_b3;
  b[4] = ellip1_b4;

  double part1 = 0.0;
  double part2 = 0.0;

  double factor = 1.0;
  for (int i = 0; i < 5; i++)
    {
      part1 += a[i] * factor; 
      part2 += b[i] * factor;

      factor = factor * comp_param;
    }
  double ellip1 = part1 - part2 * log_factor;

  return ellip1;
}

//Function to evaluate the variable part of the complete elliptic integral of the second kind
double Ellip2_var(double comp_param)
{
  vector<double> a(5);
  vector<double> b(5);

  double log_factor = log(comp_param);

  a[1] = ellip2_a1;
  a[2] = ellip2_a2;
  a[3] = ellip2_a3;
  a[4] = ellip2_a4;

  b[1] = ellip2_b1;
  b[2] = ellip2_b2;
  b[3] = ellip2_b3;
  b[4] = ellip2_b4;

  double part1 = 0.0;
  double part2 = 0.0;
  double factor = comp_param;

  for (int i = 1; i < 5; i++)
    {
      part1 += a[i] * factor;
      part2 += b[i] * factor;

      factor = factor * comp_param;
    }
  double ellip2_var = part1 - part2 * log_factor;

  return ellip2_var;
}

//Function to evaluate the complete elliptic integral of the 2st kind (Abramovitch and Stegun 1964 page 591)
double Ellip2(double comp_param)
{
  double ellip2 = 1.0 + Ellip2_var(comp_param);

  return ellip2;
}

//Function to evaluate the regular part of the complete elliptic integral of the first kind
double Ellip1_reg(double comp_param)
{
  vector<double> a(5);
  vector<double> b(5);

  double log_factor = log(comp_param);

  a[0] = ellip1_a0;
  a[1] = ellip1_a1;
  a[2] = ellip1_a2;
  a[3] = ellip1_a3;
  a[4] = ellip1_a4;

  b[1] = ellip1_b1;
  b[2] = ellip1_b2;
  b[3] = ellip1_b3;
  b[4] = ellip1_b4;

  double part1 = a[0];
  double part2 = 0.0;

  //  double factor = 1.0;
  double factor = comp_param;
  for (int i = 1; i < 5; i++)
    {
      part1 += a[i] * factor; 
      part2 += b[i] * factor;

      factor = factor * comp_param;
    }
  double ellip1 = part1 - part2 * log_factor;

  return ellip1;
}

//Function to evaluate the unbounded part of the complete elliptic integral of the first kind
double Ellip1_sing(double comp_param)
{
  double ellip1_sing = -ellip1_b0 * log(comp_param);

  return ellip1_sing;
}


//Function to evaluate the complementary parameter used as the argument for the functions which compute the elliptic integrals
double Comp_param(double beta_2, double sum)
{
  double comp_param = 1.0 - 2.0 * beta_2 / sum;

  return comp_param;
}
