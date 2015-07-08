//File containing function definitions for the supplementary functions used to build the linear system

#include <vector>
#include <math.h>

#include <iostream> //Included for debugging purposes only

#include "ellip.h"
#include "axisym.h"
#include "const.h"

using std::vector;

using std:: cout; //Using for debugging purposes only
using std:: endl;

using ellip_poly::ellip1_b0;

using math_const::PI;

//Function to assign values to the source coordinate on the interface
void Interf_source(double *out_rad, double *out_vert, double in_rad, double in_vert, double *out_rad_2)
{
  *out_rad = in_rad;
  *out_vert = in_vert;
  *out_rad_2 = in_rad * in_rad;
}

//Function to calculate the vertical distance between points and its square
void Vert_diff(double *vert_diff, double source_vert, double pos_vert, double *vert_diff_2)
{
  *vert_diff = source_vert - pos_vert;
  *vert_diff_2 = *vert_diff * *vert_diff;
}

//Function to calculate alpha_2 when the source point is on axis
double Alpha_2_axisource(double pos_rad, double vert_diff_2)
{
  double alpha_2 = pos_rad * pos_rad + vert_diff_2;

  return alpha_2;
}

//Function to calculate alpha_2
double Alpha_2(double source_rad_2, double pos_rad_2, double vert_diff_2)
{
  double alpha_2 = source_rad_2 + pos_rad_2 + vert_diff_2;

  return alpha_2;
}

//Function to calculate beta_2
double Beta_2(double source_rad, double pos_rad)
{
  double beta_2 = 2.0 * source_rad * pos_rad;

  return beta_2;
}

//Function to calculate the complementary parameter and the complete elliptic integral of the 1st and 2nd kind
void Stand_ellip(double beta_2, double sum, double *comp_param, double *ellip1, double *ellip2)
{
  *comp_param = Comp_param(beta_2, sum);
  *ellip1 = Ellip1(*comp_param);
  *ellip2 = Ellip2(*comp_param);
}

//Function to calculate the components of A when the source point is on axis
void A_axisource(double source_vert, vector<double>* pos_vert, vector<double>* pos_rad, double viscos_rat, vector<double>* matrix_A21, vector<double>* matrix_A22, vector<double>* norm_rad, vector<double>* norm_vert)
{
  //Loop over the integration points in the interval and find the values of the integrands
  for (int k = 0; k < 4; k++)
    {
      double vert_diff, vert_diff_2;

      Vert_diff(&vert_diff, source_vert, (*pos_vert)[k], &vert_diff_2);
      
      double alpha_2 = Alpha_2_axisource((*pos_rad)[k], vert_diff_2);
      double alpha_5 = pow(alpha_2, 2.5);

      (*matrix_A21)[k] = Matrix_A21_axisource(viscos_rat, vert_diff, (*pos_rad)[k], alpha_5, (*norm_vert)[k], (*norm_rad)[k]);
      (*matrix_A22)[k] = Matrix_A22_axisource(viscos_rat, vert_diff_2, alpha_5, (*pos_rad)[k], (*norm_rad)[k], vert_diff, (*norm_vert)[k]);
    }
}

//Function to calculate the components of A when the source point is off axis
void A(vector<double>* pos_rad, double source_vert, vector<double>* pos_vert, double source_rad_2, double viscos_rat, double midpoint, vector<double>* arc, vector<double>* h, vector<double>* matrix_A11, vector<double>* matrix_A12, vector<double>* matrix_A21, vector<double>* matrix_A22, double source_rad, int sing_test, vector<double>* pos_norm_rad, vector<double>* pos_norm_vert, double source_norm_vert_3, double source_norm_rad)
{
  vector<double> ellip1(4);
  vector<double> ellip2(4);
  vector<double> ellip2_var(4);

  vector<double> a1(4), a2(4), a3(4), a4(4), a6(4), a8(4), a9(4), a10(4), a11(4), a12(4), a14(4), a16(4);

  //Loop over the integration points in the interval and find the values of the integrands
  for (int k = 0; k < 4; k++)
    {
      double pos_rad_2 = (*pos_rad)[k] * (*pos_rad)[k];
      double pos_rad_4 = pos_rad_2 * pos_rad_2;
      
      double vert_diff, vert_diff_2;
      Vert_diff(&vert_diff, source_vert, (*pos_vert)[k], &vert_diff_2);
      double vert_diff_3 = vert_diff_2 * vert_diff;
	      
      double alpha_2 = Alpha_2(source_rad_2, pos_rad_2, vert_diff_2);
      double beta_2 = Beta_2(source_rad, (*pos_rad)[k]);

      double alpha_4 = alpha_2 * alpha_2;
      double beta_4 = beta_2 * beta_2;

      double alpha_6 = alpha_4 * alpha_2;

      double alpha_8 = alpha_4 * alpha_4;
      double beta_8 = beta_4 * beta_4;

      double sum = alpha_2 + beta_2;
      double sum_3_2 = pow(sum, 1.5);

      double diff = alpha_2 - beta_2;
      double diff_2 = diff * diff;

      double comp_param;
      Stand_ellip(beta_2, sum, &comp_param, &ellip1[k], &ellip2[k]);
      ellip2_var[k] = Ellip2_var(comp_param);

      a1[k] = A1(viscos_rat, sum_3_2, diff, beta_4, source_rad, alpha_2, alpha_4, source_rad_2, pos_rad_2, (*pos_rad)[k], beta_2);
      a2[k] = A2(viscos_rat, vert_diff, alpha_4, beta_4, alpha_2, vert_diff_2, sum_3_2, diff, beta_2);
      a3[k] = A3(viscos_rat, sum_3_2, diff_2, beta_4, source_rad, alpha_8, alpha_4, beta_8, pos_rad_2, source_rad_2, alpha_2, (*pos_rad)[k], beta_2);
      a4[k] = A4(viscos_rat, vert_diff, beta_2, alpha_2, alpha_4, beta_4, vert_diff_2, sum_3_2, diff, diff_2);

      a6[k] = A6(viscos_rat, vert_diff_2, source_rad_2, alpha_2, sum_3_2, diff, source_rad);
      a8[k] = A8(viscos_rat, vert_diff_2, alpha_4, beta_4, source_rad_2, alpha_2, sum_3_2, diff, source_rad, diff_2);
      
      a9[k] = A9(viscos_rat, vert_diff, alpha_4, beta_4, pos_rad_2, alpha_2, pos_rad_4, sum_3_2, diff);
      a10[k] = A10(viscos_rat, vert_diff_2, alpha_2, pos_rad_2, sum_3_2, diff, (*pos_rad)[k]);
      a11[k] = A11(viscos_rat, vert_diff, alpha_6, alpha_2, beta_4, pos_rad_2, pos_rad_4, sum_3_2, diff_2, alpha_4);
      a12[k] = A12(viscos_rat, vert_diff_2, pos_rad_2, alpha_2, alpha_4, beta_4, sum_3_2, diff_2, (*pos_rad)[k]);

      a14[k] = A14(viscos_rat, vert_diff_3, sum_3_2, diff);
      a16[k] = A16(viscos_rat, alpha_2, sum_3_2, diff_2, vert_diff_3);
      //                  cout << k << '\t' << a4[k] << '\t' << source_rad << '\t' << source_vert << '\t' << (*pos_rad)[k] << '\t' << (*pos_vert)[k] << '\t' << diff << endl;
      double arc_diff;
      if (sing_test == 1) //A11 needs to be handled differently as it is singular in the range of integration
	{
	  (*matrix_A11)[k] = Matrix_A11_reg(a1[k], a2[k], a3[k], a4[k], (*pos_norm_rad)[k], (*pos_norm_vert)[k], ellip1[k], ellip2[k], ellip2_var[k]);

	  arc_diff = -(midpoint - (*arc)[k]);

	  (*h)[k] = H(viscos_rat, sum_3_2, vert_diff_2, beta_4, source_rad, arc_diff, alpha_8, alpha_4, beta_8, alpha_2, pos_rad_2, source_rad_2, (*pos_rad)[k], (*pos_norm_rad)[k], diff_2, beta_2, source_norm_vert_3, source_norm_rad);

	}
      else
	{
	  (*matrix_A11)[k] = Matrix_A(a1[k], a2[k], a3[k], a4[k], (*pos_norm_rad)[k], (*pos_norm_vert)[k], ellip1[k], ellip2[k]);

	  (*h)[k] = 0.0;
	}

      (*matrix_A12)[k] = Matrix_A(a2[k], a6[k], a4[k], a8[k], (*pos_norm_rad)[k], (*pos_norm_vert)[k], ellip1[k], ellip2[k]);
      //	  	  cout << k << '\t' << (*matrix_A12)[k] << '\t' << a2[k] << '\t' << a6[k] << '\t' << a4[k] << '\t' << a8[k] << endl;
      (*matrix_A21)[k] = Matrix_A(a9[k], a10[k], a11[k], a12[k], (*pos_norm_rad)[k], (*pos_norm_vert)[k], ellip1[k], ellip2[k]);
      (*matrix_A22)[k] = Matrix_A(a10[k], a14[k], a12[k], a16[k], (*pos_norm_rad)[k], (*pos_norm_vert)[k], ellip1[k], ellip2[k]);
    }

}

//Function to calculate the components of B when the source point is on axis
void B_axisource(double source_vert, vector<double>* pos_vert, vector<double>* pos_rad, vector<double>* matrix_B21, vector<double>* matrix_B22)
{
  //Loop over the integration points in the interval and find the values of the integrands
  for (int k = 0; k < 4; k++)
    {
      double vert_diff, vert_diff_2;
      Vert_diff(&vert_diff, source_vert, (*pos_vert)[k], &vert_diff_2);

      double alpha_2 = Alpha_2_axisource((*pos_rad)[k], vert_diff_2);
      double alpha = sqrt(alpha_2);
      double alpha_3 = alpha_2 * alpha;

      (*matrix_B21)[k] = Matrix_B21_axisource(vert_diff, (*pos_rad)[k], alpha_3);
      (*matrix_B22)[k] = Matrix_B22_axisource(vert_diff_2, alpha, alpha_2);
    }
}

//Function to calculate the components of B when the source point is off axis
void B(vector<double>* pos_rad, double source_vert, vector<double>* pos_vert, vector<double>* matrix_B11, vector<double>* matrix_B12, vector<double>* matrix_B21, vector<double>* matrix_B22, double source_rad, double source_rad_2, int sing_test, vector<double>* g1, vector<double>*g2, double *sum, double *diff, vector<double>* matrix_B11_reg, vector<double>* matrix_B22_reg)
{
  //Loop over the integration points in the interval and find the values of the integrands
  for (int k = 0; k < 4; k++)
    {
      double pos_rad_2 = (*pos_rad)[k] * (*pos_rad)[k];

      double vert_diff, vert_diff_2;
      Vert_diff(&vert_diff, source_vert, (*pos_vert)[k], &vert_diff_2);

      double alpha_2 = Alpha_2(source_rad_2, pos_rad_2, vert_diff_2);
      double beta_2 = Beta_2(source_rad, (*pos_rad)[k]);

      *sum = alpha_2 + beta_2;
      double sum_half = sqrt(*sum);

      *diff = alpha_2 - beta_2;

      double comp_param;
      vector<double> ellip1(4);
      vector<double> ellip2(4);
      Stand_ellip(beta_2, *sum, &comp_param, &ellip1[k], &ellip2[k]);

      vector<double> ellip1_reg(4);
      vector<double> ellip1_sing(4);
      if (sing_test == 1)
	{
	  ellip1_reg[k] = Ellip1_reg(comp_param);
	  ellip1_sing[k] = Ellip1_sing(comp_param);

	  (*matrix_B11_reg)[k] = Matrix_B11_reg(beta_2, sum_half, alpha_2, vert_diff_2, *sum, *diff, ellip1_reg[k], ellip1_sing[k], ellip2[k]);
	  (*matrix_B22_reg)[k] = Matrix_B22_reg(sum_half, vert_diff_2, *diff, ellip1_reg[k], ellip2[k]);

	  (*g1)[k] = G1(alpha_2, beta_2, sum_half);
	  (*g2)[k] = G2(sum_half);

	}
      else
	{
	  (*matrix_B11)[k] = Matrix_B11(beta_2, sum_half, alpha_2, vert_diff_2, *sum, *diff, ellip1[k], ellip2[k]);
	  (*matrix_B22)[k] = Matrix_B22(sum_half, vert_diff_2, *diff, ellip1[k], ellip2[k]);
	}

      (*matrix_B12)[k] = Matrix_B12(vert_diff, source_rad, sum_half, alpha_2, *diff, ellip2[k], ellip1[k]);
      (*matrix_B21)[k] = Matrix_B21(vert_diff, (*pos_rad)[k], sum_half, alpha_2, *diff, ellip2[k], ellip1[k]);

    }

}

//Function to calculate the components of C when the source point is on axis
void C_axisource(double source_vert, vector<double>* pos_vert, vector<double>* pos_rad, vector<double>* pos_div_norm, double *temp2, vector<double>* Gauss_int_wts, double bond, double mdr, vector<double>* norm_vert, vector<double>* norm_rad)
{
  vector<double> vector_C2(4);  
  //Loop over the integration points in the interval and find the values of the integrands
  for (int k = 0; k < 4; k++)
    {
      double vert_diff, vert_diff_2;
      Vert_diff(&vert_diff, source_vert, (*pos_vert)[k], &vert_diff_2);
      
      double alpha_2 = Alpha_2_axisource((*pos_rad)[k], vert_diff_2);
      double alpha = sqrt(alpha_2);

      vector_C2[k] = Vector_C2_axisource((*pos_div_norm)[k], bond, (*pos_vert)[k], mdr, alpha, vert_diff_2, alpha_2, (*norm_vert)[k], (*pos_rad)[k], vert_diff, (*norm_rad)[k]);

      *temp2 += vector_C2[k] * (*pos_rad)[k] * (*Gauss_int_wts)[k];
    }
}

//Function to calculate the components of C when the source point is off axis
void C(double source_vert, double source_rad_2, double source_rad, vector<double>* pos_vert, vector<double>* pos_rad, vector<double>* pos_div_norm, double bond, double mdr, int sing_test, vector<double>* pos_norm_rad, vector<double>* pos_norm_vert, double *temp1, double *temp2, double mid_div_norm, double mid_vert, double mid_norm_rad, double mid_norm_vert, vector<double>* pos_arc, double midpoint, vector<double>* Gauss_int_wts)
{
  vector<double> prefac(4);

  vector<double> ellip1(4);
  vector<double> ellip2(4);
  vector<double> ellip1_reg(4);

  vector<double> vector_C1(4);
  vector<double> vector_C2(4);

  vector<double> j1(4);
  vector<double> j2(4);

  //Loop over the integration points in the interval and find the values of the integrands
  for (int k = 0; k < 4; k++)
    {
      double vert_diff, vert_diff_2;
      Vert_diff(&vert_diff, source_vert, (*pos_vert)[k], &vert_diff_2);

      double alpha_2 = Alpha_2(source_rad_2, (*pos_rad)[k] * (*pos_rad)[k], vert_diff_2);
      double beta_2 = Beta_2(source_rad, (*pos_rad)[k]);

      double beta_4 = beta_2 * beta_2;

      double sum = alpha_2 + beta_2;
      double sum_half = sqrt(sum);

      double diff = alpha_2 - beta_2;

      prefac[k] = C_prefac((*pos_div_norm)[k], bond, (*pos_vert)[k], mdr, beta_2, sum_half);

      double comp_param;
      Stand_ellip(beta_2, sum, &comp_param, &ellip1[k], &ellip2[k]);

      if (sing_test == 1)
	{
	  ellip1_reg[k] = Ellip1_reg(comp_param);
      
	  vector_C1[k] = Vector_C1_reg(vert_diff, (*pos_norm_rad)[k], (*pos_rad)[k], (*pos_norm_vert)[k], beta_4, alpha_2, vert_diff_2, source_rad, beta_2, prefac[k], ellip1[k], ellip2[k], diff, ellip1_reg[k]);

	  vector_C2[k] = Vector_C2_reg(source_rad, vert_diff, (*pos_norm_rad)[k], ellip1[k], alpha_2, (*pos_rad)[k], beta_2, (*pos_norm_vert)[k], ellip2[k], diff, ellip1_reg[k], prefac[k]);

	  j1[k] = J1(prefac[k], alpha_2, (*pos_norm_rad)[k]);
	  j2[k] = J2(prefac[k], (*pos_norm_vert)[k], beta_2);

	  double arc_diff = (*pos_arc)[k] - midpoint;

	  *temp1 += (vector_C1[k] + j1[k] * log(arc_diff * arc_diff * sum / (4.0 * source_rad_2 * diff)) + 2.0 * (9.0 * (mid_div_norm - bond * mid_vert) * ellip1_b0 * mid_norm_rad / (8.0 * PI * mdr * bond * (*pos_rad)[k]) - j1[k]) * log(fabs(arc_diff) / (2.0 * source_rad))) * (*pos_rad)[k] * (*Gauss_int_wts)[k];

	  *temp2 += (vector_C2[k] + j2[k] * log(arc_diff * arc_diff * sum / (4.0 * source_rad_2  * diff)) + 2.0 * (9.0 * (mid_div_norm - bond * mid_vert) * ellip1_b0 * mid_norm_vert / (8.0 * PI * mdr * bond * (*pos_rad)[k]) - j2[k]) * log(fabs(arc_diff) / (2.0 * source_rad))) * (*pos_rad)[k] * (*Gauss_int_wts)[k];
	}

      else
	{
	  vector_C1[k] = Vector_C1(alpha_2, vert_diff_2, (*pos_norm_rad)[k], (*pos_rad)[k], vert_diff, (*pos_norm_vert)[k], beta_4, source_rad, prefac[k], ellip1[k], ellip2[k], diff, beta_2);

	  vector_C2[k] = Vector_C2(beta_2, (*pos_norm_vert)[k], (*pos_rad)[k], vert_diff, (*pos_norm_rad)[k], source_rad, alpha_2, prefac[k], ellip1[k], ellip2[k], diff);

	  *temp1 += vector_C1[k] * (*pos_rad)[k] * (*Gauss_int_wts)[k];

	  *temp2 += vector_C2[k] * (*pos_rad)[k] * (*Gauss_int_wts)[k];
	}
    }
}
