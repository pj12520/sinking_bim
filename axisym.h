//Header file containing function declarations to evaluate the components of the quantities A, B and C

#ifndef __AXISYM_H__
#define __AXISYM_H__
#pragma once

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



#endif /* AXISYM_H */
