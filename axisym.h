//Header file containing function declarations to evaluate the components of the quantities A, B and C

#ifndef __AXISYM_H__
#define __AXISYM_H__
#pragma once

//Function to calculate the quantity a1 as defined in the notes
double A1(double viscos_rat, double sum_3_2, double diff, double beta_4, double source_rad, double alpha_2, double alpha_4, double source_rad_2, double pos_rad_2, double pos_rad, double beta_2);

//Function to calculate the quantity a2 as defined in the notes
double A2(double viscos_rat, double vert_diff, double alpha_4, double beta_4, double alpha_2, double vert_diff_2, double sum_3_2, double diff, double beta_2);

//Function to calculate the quantity a3 as defined in the notes
double A3(double viscos_rat, double sum_3_2, double diff_2, double beta_4, double source_rad, double alpha_8, double alpha_4, double beta_8, double pos_rad_2, double source_rad_2, double alpha_2, double pos_rad, double beta_2);

//Function to calculate the quantity a1 as defined in the notes
double A4(double viscos_rat, double vert_diff, double beta_2, double alpha_2, double alpha_4, double beta_4, double vert_diff_2, double sum_3_2, double diff, double diff_2);

//Function to calculate the quantity a6 as defined in the notes
double A6(double viscos_rat, double vert_diff_2, double source_rad_2, double alpha_2, double sum_3_2, double diff, double source_rad);

//Function to calculate the quantity a8 as defined in the notes
double A8(double viscos_rat, double vert_diff_2, double alpha_4, double beta_4, double source_rad_2, double alpha_2, double sum_3_2, double diff, double source_rad, double diff_2);

//Function to calculate the quantity a9 as defined in the notes
double A9(double viscos_rat, double vert_diff, double alpha_4, double beta_4, double pos_rad_2, double alpha_2, double pos_rad_4, double sum_3_2, double diff);

//Function to calculate the quantity a10 as defined in the notes
double A10(double viscos_rat, double vert_diff_2, double alpha_2, double pos_rad_2, double sum_3_2, double diff, double pos_rad);

//Function to calculate the quantity a11 as defined in the notes
double A11(double viscos_rat, double vert_diff, double alpha_6, double alpha_2, double beta_4, double pos_rad_2, double pos_rad_4, double sum_3_2, double diff_2, double alpha_4);

//Function to calculate the quantity a12 as defined in the notes
double A12(double viscos_rat, double vert_diff_2, double pos_rad_2, double alpha_2, double alpha_4, double beta_4, double sum_3_2, double diff_2, double pos_rad);

//Function to calculate the quantity a14 as defined in the notes
double A14(double viscos_rat, double vert_diff_3, double sum_3_2, double diff);

//Function to calculate the quantity a16 as defined in the notes
double A16(double viscos_rat, double alpha_2, double sum_3_2, double diff_2, double vert_diff_3);

//Function to calculate a 11 component of A for an intermediate and regular interval
double Matrix_A(double a1, double a2, double a3, double a4, double norm_rad, double norm_vert, double ellip1, double ellip2);

//Function to calculate the 21 component of A when the source point is on axis
double Matrix_A21_axisource(double viscos_rat, double vert_diff, double pos_rad, double alpha_5, double norm_vert, double norm_rad);

//Function to calculate the 22 component of A when the source point is on axis
double Matrix_A22_axisource(double viscos_rat, double vert_diff_2, double alpha_5, double pos_rad, double norm_rad, double vert_diff, double norm_vert);

//Function to calculate the regular part of the A11 function
double Matrix_A11_reg(double a1, double a2, double a3, double a4, double norm_rad, double norm_vert, double ellip1, double ellip2, double ellip2_var);

//Function to calculate the function h as defined in the notes
double H(double viscos_rat, double sum_3_2, double vert_diff_2, double beta_4, double source_rad, double eps, double alpha_8, double alpha_4, double beta_8, double alpha_2, double pos_rad_2, double source_rad_2, double pos_rad, double norm_rad, double diff_2, double beta_2, double source_norm_vert_3, double source_norm_rad);

//Function to calculate the 11 component of B for an intermediate and regular interval
double Matrix_B11(double beta_2, double sum_half, double alpha_2, double vert_diff_2, double sum, double diff, double ellip1, double ellip2);

//Function to calculate the 12 component of B for an intermediate interval
double Matrix_B12(double vert_diff, double source_rad, double sum_half, double alpha_2, double diff, double ellip2, double ellip1);

//Function to calculate the 21 component of B for an intermediate interval
double Matrix_B21(double vert_diff, double pos_rad, double sum_half, double alpha_2, double diff, double ellip2, double ellip1);

//Function to calculate the 22 component of B for an intermediate interval
double Matrix_B22(double sum_half, double vert_diff_2, double diff, double ellip1, double ellip2);

//Function to calculate the 21 component of B when the source point is on axis
double Matrix_B21_axisource(double vert_diff, double pos_rad, double alpha_3);

//Function to calculate the 22 component of B when the source point is on axis
double Matrix_B22_axisource(double vert_diff_2, double alpha, double alpha_2);

//Function to calculate the regular part of the B11 function
double Matrix_B11_reg(double beta_2, double sum_half, double alpha_2, double vert_diff_2, double sum, double diff, double ellip1_reg, double ellip1_sing, double ellip2);

//Function to calculate the regular part of the B22 function
double Matrix_B22_reg(double sum_half, double vert_diff_2, double diff, double ellip1_reg, double ellip2);

//Function to calculate the function g1 as defined in the notes
double G1(double alpha_2, double beta_2, double sum_half);

//Function to calculate the function g2 as defined in the notes
double G2(double sum_half);

//Function to calculate the prefactor in both components of C
double C_prefac(double div_norm, double bond, double pos_vert, double mdr, double beta_2, double sum_half);

//Function to calculate the 1 component of C for a regular and intermediate interval
double Vector_C1(double alpha_2, double vert_diff_2, double norm_rad, double pos_rad, double vert_diff, double norm_vert, double beta_4, double source_rad, double prefac, double ellip1, double ellip2, double diff, double beta_2);

//Function to calculate the 2 component of C for a regular and intermediate interval
double Vector_C2(double beta_2, double norm_vert, double pos_rad, double vert_diff, double norm_rad, double source_rad, double alpha_2, double prefac, double ellip1, double ellip2, double diff);

//Function to calculate the 2 component of C when the source point is on axis
double Vector_C2_axisource(double div_norm, double bond, double pos_vert, double mdr, double alpha, double vert_diff_2, double alpha_2, double norm_vert, double pos_rad, double vert_diff, double norm_rad);

//Function to calculate the regular part of the C1 function
double Vector_C1_reg(double vert_diff, double norm_rad, double pos_rad, double norm_vert, double beta_4, double alpha_2, double vert_diff_2, double source_rad, double beta_2, double prefac, double ellip1, double ellip2, double diff, double ellip1_reg);

//Function to calculate the function j1 as defined in the notes
double J1(double prefac, double alpha_2, double norm_rad);

//Function to calculate the regular part of the C2 function
double Vector_C2_reg(double source_rad, double vert_diff, double norm_rad, double ellip1, double alpha_2, double pos_rad, double beta_2, double norm_vert, double ellip2, double diff, double ellip1_reg, double prefac);

//Function to calculate the function j2 as defined in the notes
double J2(double prefac, double norm_vert, double beta_2);

#endif /* AXISYM_H */
