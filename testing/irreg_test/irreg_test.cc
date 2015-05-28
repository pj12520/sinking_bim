//Program file to test functions evaluated when computing irregular integrals over A11, C1 and C2

#include <math.h>
#include <iostream>

#include "../../axisym.h"
#include "../../ellip.h"

using std::cout;
using std::endl;

int main()
{
  //Set dimensionless variables to be used
  double viscos_rat = 2.0;
  double bond = 2.0;
  double mdr = 2.0;

  //Set location of source and position points
  double source_rad = 2.0;
  double source_vert = 3.0;

  double pos_rad = 4.0;
  double pos_vert = 5.0;

  //Set value of normal components at postion point
  double norm_rad = 0.6;
  double norm_vert = 0.8;
  double div_norm = 0.5;

  //Set values for quantities that depend on location of source and position points
  double source_rad_2 = source_rad * source_rad;
  double pos_rad_2 = pos_rad * pos_rad;

  double vert_diff = source_vert - pos_vert;
  double vert_diff_2 = vert_diff * vert_diff;

  double alpha_2 = source_rad_2 + pos_rad_2 + vert_diff_2;
  double beta_2 = 2.0 * source_rad * pos_rad;

  double alpha_4 = alpha_2 * alpha_2;
  double beta_4 = beta_2 * beta_2;

  double alpha_8 = alpha_4 * alpha_4;
  double beta_8 = beta_4 * beta_4;

  double sum = alpha_2 + beta_2;
  double sum_half = sqrt(sum);

  double sum_3_2 = sum * sum_half;

  double diff = alpha_2 - beta_2;
  double diff_2 = diff * diff;

  double arc_diff = sqrt((pos_rad - source_rad) * (pos_rad - source_rad) + (pos_vert - source_vert) * (pos_vert - source_vert));

  //Set values for complete elliptic integrals of the first and second kind

  //Compute the complementary parameter
  double comp_param = Comp_param(beta_2, sum);

  //Compute the complete elliptic integrals of the first and second kind
  double ellip1 = Ellip1(comp_param);

  double ellip2 = Ellip2(comp_param);

  double ellip2_var = ellip2 - 1.0;

  //Evaluate quantities a1, a2 etc.
  double a1 = A1(viscos_rat, sum_3_2, diff, beta_4, source_rad, alpha_2, alpha_4, source_rad_2, pos_rad_2, pos_rad, beta_2);

  double a2 = A2(viscos_rat, vert_diff, alpha_4, beta_4, alpha_2, vert_diff_2, sum_3_2, diff, beta_2);

  double a3 = A3(viscos_rat, sum_3_2, diff_2, beta_4, source_rad, alpha_8, alpha_4, beta_8, pos_rad_2, source_rad_2, alpha_2, pos_rad, beta_2);

  double a4 = A4(viscos_rat, vert_diff, beta_2, alpha_2, alpha_4, beta_4, vert_diff_2, sum_3_2, diff);

  //Evaluate the regular part of A11 and output for testing
  double matrix_A11_reg = Matrx_A11_reg(a1, a2, a3, a4, norm_rad, norm_vert, ellip1, ellip2, ellip2_var);
  cout << "The regular part of A11 = " << matrix_A11_reg << endl;

  return 0;
}
