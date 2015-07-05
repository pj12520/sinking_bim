//Header file containing declarations for supplementary functions used in building the linear system

#ifndef __BUILD_SUP_H__
#define __BUILD_SUP_H__
#pragma once

//Function to assign values to the source coordinate on the interface
void Interf_source(double *out_rad, double *out_vert, double in_rad, double in_vert, double *out_rad_2);

//Function to calculate the vertical distance between points and its square
void Vert_diff(double *vert_diff, double source_vert, double pos_vert, double *vert_diff_2);

//Function to calculate alpha_2 when the source point is on axis
double Alpha_2_axisource(double pos_rad, double vert_diff_2);

//Function to calculate alpha_2
double Alpha_2(double source_rad_2, double pos_rad_2, double vert_diff_2);

//Function to calculate beta_2
double Beta_2(double source_rad, double pos_rad);

//Function to calculate the complementary parameter and the complete elliptic integral of the 1st and 2nd kind
void Stand_ellip(double beta_2, double sum, double *comp_param, double *ellip1, double *ellip2);

//Function to calculate the component of A when the source point is on axis
void A_axisource(double source_vert, vector<double>* pos_vert, vector<double>* pos_rad, double viscos_rat, vector<double>* matrix_A21, vector<double>* matrix_A22, vector<double>* norm_rad, vector<double>* norm_vert);

//Function to calculate the components of A when the source point is off axis
void A(vector<double>* pos_rad, double source_vert, vector<double>* pos_vert, double source_rad_2, double viscos_rat, double midpoint, vector<double>* arc, vector<double>* h, vector<double>* matrix_A11, vector<double>* matrix_A12, vector<double>* matrix_A21, vector<double>* matrix_A22, double source_rad, int sing_test, vector<double>* pos_norm_rad, vector<double>* pos_norm_vert, double source_norm_vert_3, double source_norm_rad);

//Function to calculate the components of B when the source point is on axis
void B_axisource(double source_vert, vector<double>* pos_vert, vector<double>* pos_rad, vector<double>* matrix_B21, vector<double>* matrix_B22);

//Function to calculate the components of B when the source point is off axis
void B(vector<double>* pos_rad, double source_vert, vector<double>* pos_vert, vector<double>* matrix_B11, vector<double>* matrix_B12, vector<double>* matrix_B21, vector<double>* matrix_B22, double source_rad, double source_rad_2, int sing_test, vector<double>* g1, vector<double>*g2, double *sum, double *diff, vector<double>* matrix_B11_reg, vector<double>* matrix_B22_reg);

//Function to calculate the components of C when the source point is on axis
void C_axisource(double source_vert, vector<double>* pos_vert, vector<double>* pos_rad, vector<double>* pos_div_norm, double *temp2, vector<double>* Gauss_int_wts, double bond, double mdr, vector<double>* norm_vert, vector<double>* norm_rad);

//Function to calculate the components of C when the source point is off axis
void C(double source_vert, double source_rad_2, double source_rad, vector<double>* pos_vert, vector<double>* pos_rad, vector<double>* pos_div_norm, double bond, double mdr, int sing_test, vector<double>* pos_norm_rad, vector<double>* pos_norm_vert, double *temp1, double *temp2, double mid_div_norm, double mid_vert, double mid_norm_rad, double mid_norm_vert, vector<double>* pos_arc, double midpoint, vector<double>* Gauss_int_wts);

#endif /* BUILD_SUP_H */


