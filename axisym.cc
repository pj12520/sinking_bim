//File containing defintions for functions evaluating the components of A, B and C

#include <iostream> //Currently only inlcuded for debugging purposes. Remove when program is operational

#include "const.h"

using math_const::PI;

using namespace ellip_poly;

using std::cout; //Currently only using for debugging purposes. Remove when program is operational
using std::endl; //Currently only using for debugging purposes. Remove when program is operational

//Function to calculate the quantity a1 as defined in the notes
double A1(double viscos_rat, double sum_3_2, double diff, double beta_4, double source_rad, double alpha_2, double alpha_4, double source_rad_2, double pos_rad_2, double pos_rad, double beta_2)
{
  double a1 = ((1.0 - viscos_rat) / (PI * sum_3_2 * diff * beta_4)) * (source_rad * alpha_2 * (4.0 * alpha_4 - 18.0 * source_rad_2 * pos_rad_2) - (2.0 * pos_rad_2 + source_rad_2) * source_rad * (2.0 * alpha_4 - 3.0 * beta_4) - (pos_rad_2 + 2.0 * source_rad_2) * pos_rad * alpha_2 * beta_2 + source_rad * pos_rad_2 * beta_4);

  //  cout << sum_3_2 << '\t' << diff << '\t' << beta_4 << '\t' << a1 << endl;
  return a1;
}

//Function to calculate the quantity a2 as defined in the notes
double A2(double viscos_rat, double vert_diff, double alpha_4, double beta_4, double alpha_2, double vert_diff_2, double sum_3_2, double diff, double beta_2)
{
  double a2 = (1.0 - viscos_rat) * vert_diff * (2.0 * alpha_4 - 2.0 * beta_4 - alpha_2 * vert_diff_2) / (PI * sum_3_2 * diff * beta_2);

  return a2;
}

//Function to calculate the quantity a3 as defined in the notes
double A3(double viscos_rat, double sum_3_2, double diff_2, double beta_4, double source_rad, double alpha_8, double alpha_4, double beta_8, double pos_rad_2, double source_rad_2, double alpha_2, double pos_rad, double beta_2)
{
  double a3 = ((1.0 - viscos_rat) / (PI * sum_3_2 * diff_2 * beta_4)) * ((source_rad / 2.0) * (-8.0 * alpha_8 + 15.0 * alpha_4 * beta_4 - 3.0 * beta_8) - 2.0 * (2.0 * pos_rad_2 + source_rad_2) * source_rad * alpha_2 * (-alpha_4 + 3.0 * beta_4) + (pos_rad_2 + 2.0 * source_rad_2) * pos_rad * beta_2 * (alpha_4 + 3.0 * beta_4) - 4.0 * source_rad * pos_rad_2 * alpha_2 * beta_4);

  return a3;
}
//Function to calculate the quantity a4 as defined in the notes
double A4(double viscos_rat, double vert_diff, double beta_2, double alpha_2, double alpha_4, double beta_4, double vert_diff_2, double sum_3_2, double diff, double diff_2)
{
  double a4 = - (1.0 - viscos_rat) * vert_diff * beta_2 * (alpha_2 * (alpha_4 - 5.0 * beta_4) + (alpha_2 - vert_diff_2)* (alpha_4 + 3.0 * beta_4)) / (PI * sum_3_2 * diff_2 * beta_4);
  //cout << a4 << '\t' << -vert_diff * beta_2 * (1.0 - viscos_rat) << endl;
  return a4;
}

//Function to calculate the quantity a6 as defined in the notes
double A6(double viscos_rat, double vert_diff_2, double source_rad_2, double alpha_2, double sum_3_2, double diff, double source_rad)
{
  double a6 = (1.0 - viscos_rat) * vert_diff_2 * (2.0 * source_rad_2 - alpha_2) / (2.0 * PI * sum_3_2 * diff * source_rad);

  return a6;
}

//Function to calculate the quantity a8 as defined in the notes
double A8(double viscos_rat, double vert_diff_2, double alpha_4, double beta_4, double source_rad_2, double alpha_2, double sum_3_2, double diff, double source_rad, double diff_2)
{
  double a8 = (1.0 - viscos_rat) * vert_diff_2 * (alpha_4 + 3.0 * beta_4 - 8.0 * source_rad_2 * alpha_2) / (2.0 * PI * sum_3_2 * diff_2 * source_rad);

  return a8;
}

//Function to calculate the quantity a9 as defined in the notes
double A9(double viscos_rat, double vert_diff, double alpha_4, double beta_4, double pos_rad_2, double alpha_2, double pos_rad_4, double sum_3_2, double diff)
{
  double a9 = (1.0 - viscos_rat) * vert_diff * (-2.0 * alpha_4 + 3.0 * beta_4 - 4.0 * pos_rad_2 * alpha_2 + 4.0 * pos_rad_4) / (4.0 * PI * sum_3_2 * diff * pos_rad_2);

  return a9;
}

//Function to calculate the quantity a10 as defined in the notes
double A10(double viscos_rat, double vert_diff_2, double alpha_2, double pos_rad_2, double sum_3_2, double diff, double pos_rad)
{
  double a10 = (1.0 - viscos_rat) * vert_diff_2 * (alpha_2 - 2.0 * pos_rad_2) / (2.0 * PI * sum_3_2 * diff * pos_rad);

  return a10;
}

//Function to calculate the quantity a11 as defined in the notes
double A11(double viscos_rat, double vert_diff, double alpha_6, double alpha_2, double beta_4, double pos_rad_2, double pos_rad_4, double sum_3_2, double diff_2, double alpha_4)
{
  double a11 = (1.0 - viscos_rat) * vert_diff * (alpha_6 - 3.0 * alpha_2 * beta_4 + 2.0 * pos_rad_2 * alpha_4 + 6.0 * pos_rad_2 * beta_4 - 8.0 * pos_rad_4 * alpha_2) / (2.0 * PI * sum_3_2 * diff_2 * pos_rad_2);

  return a11;
}

//Function to calculate the quantity a12 as defined in the notes
double A12(double viscos_rat, double vert_diff_2, double pos_rad_2, double alpha_2, double alpha_4, double beta_4, double sum_3_2, double diff_2, double pos_rad)
{
  double a12 = (1.0 - viscos_rat) * vert_diff_2 * (8.0 * pos_rad_2 * alpha_2 - alpha_4 - 3.0 * beta_4) / (2.0 * PI * sum_3_2 * diff_2 * pos_rad);

  return a12;
}

//Function to calculate the quantity a14 as defined in the notes
double A14(double viscos_rat, double vert_diff_3, double sum_3_2, double diff)
{
  double a14 = (1.0 - viscos_rat) * vert_diff_3 / (PI * sum_3_2 * diff);

  return a14;
}

//Function to calculate the quantity a16 as defined in the notes
double A16(double viscos_rat, double alpha_2, double sum_3_2, double diff_2, double vert_diff_3)
{
  double a16 = -4.0 * (1.0 - viscos_rat) * vert_diff_3 * alpha_2 / (PI * sum_3_2 * diff_2);

  return a16;
}

//Function to calculate a component of A for an intermediate and regular interval
double Matrix_A(double a1, double a2, double a3, double a4, double norm_rad, double norm_vert, double ellip1, double ellip2)
{
  double matrix_A = (a1 * norm_rad + a2 * norm_vert) * ellip1 + (a3 * norm_rad + a4 * norm_vert) * ellip2;

  return matrix_A;
}

//Function to calculate the 21 component of A when the source point is on axis
double Matrix_A21_axisource(double viscos_rat, double vert_diff, double pos_rad, double alpha_5, double norm_vert, double norm_rad)
{
  double matrix_A21 = 3.0 * (1.0 - viscos_rat) * vert_diff * pos_rad * (vert_diff * norm_vert - pos_rad * norm_rad) / (2.0 * alpha_5);

  return matrix_A21;
}

//Function to calculate the 22 component of A when the source point is on axis
double Matrix_A22_axisource(double viscos_rat, double vert_diff_2, double alpha_5, double pos_rad, double norm_rad, double vert_diff, double norm_vert)
{
  double matrix_A22 = 3.0 * (1.0 - viscos_rat) * vert_diff_2 * (pos_rad * norm_rad - vert_diff * norm_vert) / (2.0 * alpha_5);

  return matrix_A22;
}

//Function to calculate the regular part of the A11 function
double Matrix_A11_reg(double a1, double a2, double a3, double a4, double norm_rad, double norm_vert, double ellip1, double ellip2, double ellip2_var)
{
  double matrix_A11_reg = (a1 * norm_rad + a2 * norm_vert) * ellip1 + a4 * norm_vert * ellip2 + a3 * norm_rad * ellip2_var;

  return matrix_A11_reg;
}

//Function to calculate the function h as defined in the notes
double H(double viscos_rat, double sum_3_2, double vert_diff_2, double beta_4, double source_rad, double eps, double alpha_8, double alpha_4, double beta_8, double alpha_2, double pos_rad_2, double source_rad_2, double pos_rad, double norm_rad, double diff_2, double beta_2, double source_norm_vert_3, double source_norm_rad)
{
  double prefac = (1.0 - viscos_rat) / (PI * sum_3_2 * diff_2 * beta_4 * source_rad * eps);

  double term1 = source_rad * (- 8.0 * alpha_8 + 15.0 * alpha_4 * beta_4 - 3.0 * beta_8) / 2.0;

  double term2 = -2.0 * source_rad * alpha_2 * (2.0 * pos_rad_2 + source_rad_2) * (-alpha_4 + 3.0 * beta_4);

  double term3 = pos_rad * beta_2 * (pos_rad_2 + 2.0 * source_rad_2) * (alpha_4 + 3.0 * beta_4);

  double term4 = -4.0 * source_rad * pos_rad_2 * alpha_2 * beta_4;

  double term5 = - source_norm_vert_3 * source_norm_rad * sum_3_2 * diff_2 * beta_4;

  double h = prefac * (norm_rad * source_rad * eps * (term1 + term2 + term3 + term4) + term5);

  return h;
}

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
  // cout << ellip2 << " " << matrix_B11_reg <<endl;
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

//Function to calculate the prefactor in both components of C
double C_prefac(double div_norm, double bond, double pos_vert, double mdr, double beta_2, double sum_half)
{
  double prefac = 9.0 * (div_norm - bond * pos_vert) / (4.0 * PI * mdr * bond * beta_2 * sum_half);

  return prefac;
}

//Function to calculate the 1 component of C for a regular and intermediate interval
double Vector_C1(double alpha_2, double vert_diff_2, double norm_rad, double pos_rad, double vert_diff, double norm_vert, double beta_4, double source_rad, double prefac, double ellip1, double ellip2, double diff, double beta_2)
{
  double term1 = (alpha_2 + vert_diff_2) * norm_rad + pos_rad * vert_diff * norm_vert;

  double term2 = norm_rad * (beta_4 - alpha_2 * (alpha_2 + vert_diff_2));

  double term3 = norm_vert * vert_diff * (source_rad * beta_2 - pos_rad * alpha_2);

  double vector_C1 = prefac * (term1 * ellip1 + (term2 + term3) * ellip2 / diff);

  return vector_C1;
}

//Function to calculate the 2 component of C for a regular and intermediate interval
double Vector_C2(double beta_2, double norm_vert, double pos_rad, double vert_diff, double norm_rad, double source_rad, double alpha_2, double prefac, double ellip1, double ellip2, double diff)
{
  double term1 = beta_2 * norm_vert - source_rad * vert_diff * norm_rad;

  double term2 = norm_rad * (source_rad * alpha_2 - pos_rad * beta_2);

  double term3 = vert_diff * beta_2 * norm_vert;

  double vector_C2 = prefac * (term1 * ellip1 + (term2 + term3) * vert_diff * ellip2 / diff);

  return vector_C2;
}

//Function to calculate the 2 component of C when the source point is on axis
double Vector_C2_axisource(double div_norm, double bond, double pos_vert, double mdr, double alpha, double vert_diff_2, double alpha_2, double norm_vert, double pos_rad, double vert_diff, double norm_rad)
{
  double vector_C2 = 9.0 * (div_norm - bond * pos_vert) / (8.0 * mdr * bond * alpha) * (norm_vert + (vert_diff_2 * norm_vert - pos_rad * vert_diff * norm_rad) / alpha_2);

  return vector_C2;
}

//Function to calculate the regular part of the C1 function
double Vector_C1_reg(double vert_diff, double norm_rad, double pos_rad, double norm_vert, double beta_4, double alpha_2, double vert_diff_2, double source_rad, double beta_2, double prefac, double ellip1, double ellip2, double diff, double ellip1_reg)
{
  double term1 = vert_diff * (vert_diff * norm_rad + pos_rad * norm_vert);

  double term2 = norm_rad * (beta_4 - alpha_2 * (alpha_2 + vert_diff_2));

  double term3 = norm_vert * vert_diff * (source_rad * beta_2 - pos_rad * alpha_2);

  double vector_C1_reg = prefac * (term1 * ellip1 + ellip2 * (term2 + term3) / diff + alpha_2 * norm_rad * ellip1_reg);

  return vector_C1_reg;
}

//Function to calculate the function j1 as defined in the notes
double J1(double prefac, double alpha_2, double norm_rad)
{
  double j1 = prefac * alpha_2 * ellip1_b0 * norm_rad;

  return j1;
}

//Function to calculate the regular part of the C2 function
double Vector_C2_reg(double source_rad, double vert_diff, double norm_rad, double ellip1, double alpha_2, double pos_rad, double beta_2, double norm_vert, double ellip2, double diff, double ellip1_reg, double prefac)
{
  double term1 = - source_rad * vert_diff * norm_rad * ellip1;

  double term2 = vert_diff * (norm_rad * (source_rad * alpha_2 - pos_rad * beta_2) + vert_diff * beta_2 * norm_vert) * ellip2 / diff;

  double term3 = beta_2 * norm_vert * ellip1_reg;

  double vector_C2_reg = prefac * (term1 + term2 + term3);

  return vector_C2_reg;
}

//Function to calculate the function j2 as defined in the notes
double J2(double prefac, double norm_vert, double beta_2)
{
  double j2 = prefac * norm_vert * ellip1_b0 * beta_2;

  return j2;
}
