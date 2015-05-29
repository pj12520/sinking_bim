//File containing functions relating to the construction of the linear system

#include <vector>
#include <math.h>

#include "object.h"
#include "const.h"
#include "axisym.h"
#include "ellip.h"

using std::vector;

using Gauss::Gauss_wt1;
using Gauss::Gauss_wt2;
using Gauss::Gauss_wt3;
using Gauss::Gauss_wt4;

//Function to create the linear system
void Build(vector<vector<double> >* matrix, vector<double>* vec, particle sphere, surf interf, double viscos_rat)
{
  //Temporary vector object to store the matrix of coefficients and known vector
  vector<vector<double> > coeffs((*matrix).size());
  vector<double> known((*vec).size());
  for (int i = 0; i < coeffs.size(); i++)
    {
      coeffs[i].resize((*matrix)[i].size());
    }

  vector<double> Gauss_int_wts(4); //Vector to store weights used for 4-pt Gaussian quadrature

  Gauss_int_wts[0] = Gauss_wt1;
  Gauss_int_wts[1] = Gauss_wt2;
  Gauss_int_wts[2] = Gauss_wt3;
  Gauss_int_wts[3] = Gauss_wt4;

  //Variables to contain position of source when creating matrix
  double source_rad;
  double source_vert;

  //Temporary variables that are used in evaluating the matrix components
  double source_rad_2; //Square of the radial coordinate of the source point
  double source_vert_2;//Square of the vertical coordinate of the source point

  double pos_rad_2; //Square of the radial coordinate of the integration point
  double pos_rad_4; //Radial coordinate of the integration point raised to the 4th power

  double vert_diff; //Difference between vertical position of source and integration points
  double vert_diff_2; //vert_diff squared
  double vert_diff_3; //vert_diff cubed

  double alpha; 
  double alpha_2; //Alpha squared
  double alpha_3; //Alpha cubed
  double alpha_4; //etc.
  double alpha_5;
  double alpha_6;
  double alpha_8;

  double beta_2;
  double beta_4;
  double beta_6;
  double beta_8;

  double sum; //alpha_2 + beta_2
  double sum_half; //Square root of sum
  double sum_3_2; //Sum^(3/2)

  double diff; //alpha_2 - beta_2
  double diff_2; //diff squared

  double comp_param; //Argument that is passed to functions which evaluate elliptic integrals

  double theta_diff; //Separation in theta between source and integration point when they are both in the same interval on the sphere
  double theta_diff_2; //theta_diff squared

  double arc_diff; //Separation in arc length between source and integration point when they are both in the same interval on the interface
  double arc_diff_2; //arc_diff squared

  //Vectors to contain temporary values of the elliptic integrals at the integration points
  vector<double> ellip1(4);
  vector<double> ellip2(4);

  vector<double> ellip1_reg(4);
  vector<double> ellip1_sing(4);

  vector<double> ellip2_var(4);

  //Vectors to contain temporary values of the integrands at the integration points
  vector<double> a1(4);
  vector<double> a2(4);
  vector<double> a3(4);
  vector<double> a4(4);

  vector<double> a6(4);
  vector<double> a8(4);

  vector<double> a9(4);
  vector<double> a10(4);
  vector<double> a11(4);
  vector<double> a12(4);

  vector<double> a14(4);
  vector<double> a16(4);

  vector<double> matrix_A11(4);
  vector<double> matrix_A12(4);
  vector<double> matrix_A21(4);
  vector<double> matrix_A22(4);

  vector<double> matrix_A11_reg(4);

  vector<double> h(4);

  vector<double> matrix_B11(4);
  vector<double> matrix_B12(4);
  vector<double> matrix_B21(4);
  vector<double> matrix_B22(4);

  vector<double> matrix_B11_reg(4);
  vector<double> matrix_B22_reg(4);

  vector<double> g1(4);
  vector<double> g2(4);

  vector<double> vector_C1(4);
  vector<double> vector_C2(4);
 
  vector<double> vector_C1_reg(4);
  vector<double> vector_C2_reg(4);

  vector<double> j1(4);
  vector<double> j2(4);

  //Fill up the matrix as defined in the notes

  //Start looping over the source points on the interface
  for (int i = 0; i < interf.n_int; i++)
    {
      source_rad = interf.mid_rad[i];
      source_vert = interf.mid_vert[i];

      source_rad_2 = source_rad * source_rad;

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  //For the case that the source point is on axis
	  if (i == 0)
	    {
	      //Loop over the integration points in the interval and find the values of the integrands
	      for (int k = 0; k < 4; k++)
		{
		  vert_diff = source_vert - interf.intervals[j].vert[k];
		  vert_diff_2 = vert_diff * vert_diff;
		  vert_diff_3 = vert_diff_2 * vert_diff;

		  alpha_2 = interf.intervals[j].rad[k] * interf.intervals[j].rad[k] + vert_diff_2;
		  alpha_5 = pow(alpha_2, 2.5);

		  matrix_A21[k] = Matrix_A21_axisource(viscos_rat, vert_diff_2, interf.intervals[j].rad[k], alpha_5);
		  matrix_A22[k] = Matrix_A22_axisource(viscos_rat, vert_diff_3, alpha_5);
		}

	      //Perform the Gauss-Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      coeffs[i][j] = 0.0;
	      coeffs[i][j + interf.n_int] = 0.0;

	      coeffs[i + interf.n_int][j] = 0.0;
	      coeffs[i + interf.n_int][j + interf.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + interf.n_int][j] += matrix_A21[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + interf.n_int] += matrix_A22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + interf.n_int][j] = interf.intervals[j].width * coeffs[i + interf.n_int][j] / 2.0;
	      coeffs[i + interf.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + interf.n_int][j] / 2.0;

	      if (j == i)
		{
		  coeffs[i][j] -= (1.0 + viscos_rat) / 2.0;
		  coeffs[i + interf.n_int][j + interf.n_int] -= (1.0 + viscos_rat) / 2.0;
		}
	    }

	  //For the case that the source point is not on axis
	  else 
	    {
	      //Loop over the integration points in the interval and find the values of the integrands
	      for (int k = 0; k < 4; k++)
		{
		  pos_rad_2 = interf.intervals[j].rad[k];
		  pos_rad_4 = pos_rad_2 * pos_rad_2;

		  vert_diff = source_vert - interf.intervals[j].vert[k];
		  vert_diff_2 = vert_diff * vert_diff;
	      
		  alpha_2 = source_rad_2 + pos_rad_2 + vert_diff_2;
		  beta_2 = 2.0 * source_rad * interf.intervals[j].rad[k];

		  alpha_4 = alpha_2 * alpha_2;
		  beta_4 = beta_2 * beta_2;

		  alpha_6 = alpha_4 * alpha_2;

		  alpha_8 = alpha_4 * alpha_4;
		  beta_8 = beta_4 * beta_4;

		  sum = alpha_2 + beta_2;
		  sum_3_2 = pow(sum, 1.5);

		  diff = alpha_2 - beta_2;
		  diff_2 = diff * diff;

		  comp_param = Comp_param(beta_2, sum);
		  ellip1 = Ellip1(comp_param);
		  ellip2 = Ellip2(comp_param);
		  ellip2_var = Ellip2_var(comp_param);

		  a1[k] = A1(viscos_rat, sum_3_2, diff, beta_4, source_rad, alpha_2, alpha_4, source_rad_2, pos_rad_2, interf.intervals[j].rad[k], beta_2);
		  a2[k] = A2(viscos_rat, vert_diff, alpha_4, beta_4, alpha_2, vert_diff_2, sum_3_2, diff, beta_2);
		  a3[k] = A3(viscos_rat, sum_3_2, diff_2, beta_4, source_rad, alpha_8, alpha_4, beta_8, pos_rad_2, source_rad_2, alpha_2, interf.intervals[j].rad[k], beta_2);
		  a4[k] = A4(viscos_rat, vert_diff, beta_2, alpha_2, alpha_4, beta_4, vert_diff_2, sum_3_2, diff);

		  a6[k] = A6(viscos_rat, vert_diff_2, source_rad_2, alpha_2, sum_3_2, diff, source_rad);
		  a8[k] = A8(viscos_rat, vert_diff_2, alpha_4, beta_4, source_rad_2, alpha_2, sum_3_2, diff, source_rad);

		  a9[k] = A9(viscos_rat, vert_diff, alpha_4, beta_4, pos_rad_2, alpha_2, pos_rad_4, sum_3_2, diff);
		  a10[k] = A10(viscos_rat, vert_diff_2, alpha_2, pos_rad_2, sum_3_2, diff, interf.intervals[j].rad[k]);
		  a11[k] = A11(viscos_rat, vert_diff, alpha_6, alpha_2, beta_4, pos_rad_2, pos_rad_4, sum_3_2, diff_2, alpha_4);
		  a12[k] = A12(viscos_rat, vert_diff_2, pos_rad_2, alpha_2, alpha_4, beta_4, sum_3_2, diff_2, interf.intervals[j].rad[k]);

		  if (i == j) //A11 needs to be handled differently as it is singular in the range of integration
		    {
		      matrix_A11[k] = Matrix_A11_reg(a1, a2, a3, a4, interf.intervals[j].norm_rad[k], interf.intervals[j].norm_vert[k], ellip1, ellip2, ellip2_var);

		      arc_diff = interf.midpoints[i] - interf.intervals[j].arc[k];

		    }
		  else
		    {

		    }

		}
	    }
	      
	  //For the case that the source point is in the region of integration for an intermediate interval
	  else if (i == j)
	    {
	      for (int k = 0; k < 4; k++) //Loop over integration points in interval
		{
		  //Define temporary variables
		  vert_diff = source[1] - (*intervals)[j].vert[k];
		  vert_diff_2 = vert_diff * vert_diff;

		  alpha_2 = source[0] * source[0] + (*intervals)[j].rad[k] * (*intervals)[j].rad[k] + vert_diff_2;

		  beta_2 = 2.0 * source[0] * (*intervals)[j].rad[k];

		  sum = alpha_2 + beta_2;
		  sum_half = sqrt(sum);

		  diff = alpha_2 - beta_2;

		  comp_param = 1.0 - 2.0 * beta_2 / sum;

		  ellip1_reg[k] = Ellip1_reg(comp_param);
		  ellip1_sing[k] = Ellip1_sing(comp_param);

		  ellip1[k] = Ellip1(comp_param);
		  ellip2[k] = Ellip2(comp_param);

		  matrix_B11_reg[k] = Matrix_B11_reg(beta_2, sum_half, alpha_2, vert_diff_2, sum, diff, ellip1_reg[k], ellip1_sing[k], ellip2[k]);
		  matrix_B22_reg[k] = Matrix_B22_reg(sum_half, vert_diff_2, diff, ellip1_reg[k], ellip2[k]);

		  matrix_B12[k] = Matrix_B12(vert_diff, source[0], sum_half, alpha_2, diff, ellip2[k], ellip1[k]);
		  matrix_B21[k] = Matrix_B21(vert_diff, (*intervals)[j].rad[k], sum_half, alpha_2, diff, ellip2[k], ellip1[k]);

		  g1[k] = G1(alpha_2, beta_2, sum_half);
		  g2[k] = G2(sum_half);
		}

	      //Perform the Gauss Legendre integration
	      coeffs[i - 1][j - 1] = 0.0;
	      coeffs[i - 1][j + n_int - 2] = 0.0;
	      coeffs[i + n_int - 2][j - 1] = 0.0;
	      coeffs[i + n_int - 2][j + n_int - 2] = 0.0;
	      
	      for (int k = 0; k < 4; k++)
		{
		  theta_diff = (*intervals)[j].theta[k] - midpoints[i];
		  theta_diff_2 = theta_diff * theta_diff;

		  coeffs[i - 1][j - 1] += (matrix_B11_reg[k] + g1[k] * log(4.0 * source[0] * source[0] * diff / (sum * theta_diff_2)) + 2.0 * (g1[k] + ellip1_b0 / (4.0 * PI * source[0])) * log(fabs(theta_diff) / (2.0 * source[0]))) * Gauss_int_wts[k];

		  coeffs[i - 1][j + n_int - 2] += matrix_B12[k] * Gauss_int_wts[k];
		  
		  coeffs[i + n_int - 2][j - 1] += matrix_B21[k] * Gauss_int_wts[k];

		  coeffs[i + n_int - 2][j + n_int - 2] += (matrix_B22_reg[k] + g2[k] * log(4.0 * source[0] * source[0] * diff / (sum * theta_diff_2)) + 2.0 * (g2[k] + ellip1_b0 / (4.0 * PI * source[0])) * log(fabs(theta_diff) / (2.0 * source[0]))) * Gauss_int_wts[k];
		}

	      coeffs[i - 1][j - 1] = (*intervals)[j].width * coeffs[i - 1][j - 1] / 2.0 - ellip1_b0 / (2.0 * (n_int - 1.0) * source[0]) * (log(PI / (4.0 * (n_int - 1.0) * source[0])) - 1.0);

	      coeffs[i - 1][j + n_int - 2] = (*intervals)[j].width * coeffs[i - 1][j + n_int - 2] / 2.0;

	      coeffs[i + n_int - 2][j - 1] = (*intervals)[j].width * coeffs[i + n_int - 2][j - 1] / 2.0;

	      coeffs[i + n_int - 2][j + n_int - 2] = (*intervals)[j].width * coeffs[i + n_int - 2][j + n_int - 2] / 2.0 - ellip1_b0 / (2.0 * (n_int - 1.0) * source[0]) * (log(PI / (4.0 * (n_int - 1.0) * source[0])) - 1.0);
	    }

	  else
	    {
	      for (int k = 0; k < 4; k++) //Loop over integration points in interval
		{
		  //Define temporary variables
		  vert_diff = source[1] - (*intervals)[j].vert[k];
		  vert_diff_2 = vert_diff * vert_diff;

		  alpha_2 = source[0] * source[0] + (*intervals)[j].rad[k] * (*intervals)[j].rad[k] + vert_diff_2;

		  beta_2 = 2.0 * source[0] * (*intervals)[j].rad[k];

		  sum = alpha_2 + beta_2;
		  sum_half = sqrt(sum);

		  diff = alpha_2 - beta_2;

		  comp_param = 1.0 - 2.0 * beta_2 / sum;

		  ellip1[k] = Ellip1(comp_param);
		  ellip2[k] = Ellip2(comp_param);

		  matrix_B11[k] = Matrix_B11(beta_2, sum_half, alpha_2, vert_diff_2, sum, diff, ellip1[k], ellip2[k]);
		  matrix_B12[k] = Matrix_B12(vert_diff, source[0], sum_half, alpha_2, diff, ellip2[k], ellip1[k]);
		  matrix_B21[k] = Matrix_B21(vert_diff, (*intervals)[j].rad[k], sum_half, alpha_2, diff, ellip2[k], ellip1[k]);
		  matrix_B22[k] = Matrix_B22(sum_half, vert_diff_2, diff, ellip1[k], ellip2[k]);
		}

	      //Perform the Gauss Legendre integration
	      if (j != 0 && j != n_int - 1)
		{
		  coeffs[i - 1][j - 1] = 0.0;
		  coeffs[i + n_int - 2][j - 1] = 0.0;
		}
	      coeffs[i - 1][j + n_int - 2] = 0.0;
	      coeffs[i + n_int - 2][j + n_int - 2] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  if (j != 0 && j != n_int - 1)
		    {
		      coeffs[i - 1][j - 1] += matrix_B11[k] * Gauss_int_wts[k];
		      coeffs[i + n_int - 2][j - 1] += matrix_B21[k] * Gauss_int_wts[k];
		    }
		  coeffs[i - 1][j + n_int - 2] += matrix_B12[k] * Gauss_int_wts[k];
		  coeffs[i + n_int - 2][j + n_int - 2] += matrix_B22[k] * Gauss_int_wts[k];
		}
	      if (j != 0 && j != n_int - 1)
		{
		  coeffs[i - 1][j - 1] = (*intervals)[j].width * coeffs[i - 1][j - 1] / 2.0;
		  coeffs[i + n_int - 2][j - 1] = (*intervals)[j].width * coeffs[i + n_int - 2][j - 1] / 2.0;
		}
	      coeffs[i - 1][j + n_int - 2] = (*intervals)[j].width * coeffs[i - 1][j + n_int - 2] / 2.0;
	      coeffs[i + n_int - 2][j + n_int - 2] = (*intervals)[j].width * coeffs[i + n_int - 2][j + n_int - 2] / 2.0;
	    }
	}
    }

}
