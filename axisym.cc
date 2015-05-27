//File containing defintions for functions evaluating the components of A, B and C

#include "const.h"

using math_const::PI;

using namespace ellip_poly;

//Function to calculate the 11 component of B for an intermediate and regular interval
double Matrix_B11(double beta_2, double sum_half, double alpha_2, double vert_diff_2, double sum, double diff, double ellip1, double ellip2)
{
  double matrix_B11 = ((alpha_2 + vert_diff_2) * ellip1 - (sum + alpha_2 * vert_diff_2 / diff) * ellip2) / (2.0 * PI * beta_2 * sum_half);

  return matrix_B11;
}

//Function to calculate the 12 component of B for an intermediate and regular interval
double Matrix_B12(double vert_diff, double source_rad, double sum_half, double alpha_2, double diff, double ellip2, double ellip1)
{
  double matrix_B12 = vert_diff / (4.0 * PI * source_rad * sum_half) * ((2.0 * source_rad * source_rad - alpha_2) * ellip2 / diff + ellip1);

  return matrix_B12;
}

//Function to calculate the 21 component of B for an intermediate and regular interval
double Matrix_B21(double vert_diff, double pos_rad, double sum_half, double alpha_2, double diff, double ellip2, double ellip1)
{
  double matrix_B21 = vert_diff / (4.0 * PI * pos_rad * sum_half) * ((alpha_2 - 2.0 * pos_rad * pos_rad) * ellip2 / diff - ellip1);

  return matrix_B21;
}

//Function to calculate the 22 component of B for an intermediate interval
double Matrix_B22(double sum_half, double vert_diff_2, double diff, double ellip1, double ellip2)
{
  double matrix_B22 = (ellip1 + vert_diff_2 * ellip2 / diff) / (2.0 * PI * sum_half);

  return matrix_B22;
}

//Function to calculate the 21 component of B when the source point is on axis
double Matrix_B21_axisource(double vert_diff, double pos_rad, double alpha_3)
{
  double matrix_B21 = - vert_diff * pos_rad / (4.0 * alpha_3);

  return matrix_B21;
}

//Function to calculate the 22 component of B when the source point is on axis
double Matrix_B22_axisource(double vert_diff_2, double alpha, double alpha_2)
{
  double matrix_B22 = 1.0 / (4.0 * alpha) * (1.0 + vert_diff_2 / alpha_2);

  return matrix_B22;
}

//Function to calculate the regular part of the B11 function
double Matrix_B11_reg(double beta_2, double sum_half, double alpha_2, double vert_diff_2, double sum, double diff, double ellip1_reg, double ellip1_sing, double ellip2)
{
  double matrix_B11_reg = ((alpha_2 + vert_diff_2) * ellip1_reg + vert_diff_2 * ellip1_sing - (sum + alpha_2 * vert_diff_2 / diff) * ellip2) / (2.0 * PI * beta_2 * sum_half);

  return matrix_B11_reg;
}

//Function to calculate the function g1 as defined in the notes
double G1(double alpha_2, double beta_2, double sum_half)
{
  double g1 = - alpha_2 * ellip1_b0 / (2.0 * PI * beta_2 * sum_half);

  return g1;
}

//Function to calculate the regular part of the B22 function
double Matrix_B22_reg(double sum_half, double vert_diff_2, double diff, double ellip1_reg, double ellip2)
{
  double matrix_B22_reg = (ellip1_reg + vert_diff_2 * ellip2 / diff) / (2.0 * PI * sum_half);

  return matrix_B22_reg;
}

//Function to calculate the function g2 as defined in the notes
double G2(double sum_half)
{
  double g2 = -ellip1_b0 / (2.0 * PI * sum_half);

  return g2;
}
