//File containing functions relating to the construction of the linear system

#include <vector>
#include <math.h>

#include "object.h"
#include "const.h"
#include "axisym.h"
#include "ellip.h"
#include "build_sup.h"

using std::vector;

using Gauss::Gauss_wt1;
using Gauss::Gauss_wt2;
using Gauss::Gauss_wt3;
using Gauss::Gauss_wt4;

using ellip_poly::ellip1_b0;

using math_const::PI;

//Function to create the linear system
void Build(vector<vector<double> >* matrix, vector<double>* vec, particle sphere, surf interf, double viscos_rat, double bond, double mdr)
{
  //Temporary vector object to store the matrix of coefficients and known vector
  vector<vector<double> > coeffs((*matrix).size());
  vector<double> known((*vec).size());
  for (int i = 0; i < coeffs.size(); i++)
    {
      coeffs[i].resize((*matrix)[i].size());
    }

  vector<double> temp1(interf.n_int); //Vectors to store temporary values used in the calculation of the known vector. Needs to be resized in the calculation
  vector<double> temp2(interf.n_int); 

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

  double source_norm_vert_3; //Cube of the vertical component of the normal vector at the source point

  int sing_test;

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

  vector<double> prefac(4);

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
      Interf_source(&source_vert, &source_rad, interf.mid_rad[i], interf.mid_vert[i], &source_rad_2);

      source_norm_vert_3 = interf.mid_norm_vert[i] * interf.mid_norm_vert[i] * interf.mid_norm_vert[i];

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  //For the case that the source point is on axis
	  if (i == 0)
	    {
	      A_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), viscos_rat, &matrix_A21, &matrix_A22);

	      //Perform the Gauss-Legendre integration
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
	      coeffs[i + interf.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + interf.n_int][j + interf.n_int] / 2.0;

	      if (j == i)
		{
		  coeffs[i][j] -= (1.0 + viscos_rat) / 2.0;
		  coeffs[i + interf.n_int][j + interf.n_int] -= (1.0 + viscos_rat) / 2.0;
		}
	    }

	  //For the case that the source point is not on axis
	  else 
	    {
	      if (i == j)
		{
		  sing_test = 1;
		}
	      else
		{
		  sing_test = 0;
		}

	      A(&(interf.intervals[j].rad), source_vert, &(interf.intervals[j].vert), source_rad_2, viscos_rat, interf.midpoints[i], &(interf.intervals[j].arc), &h, &matrix_A11, &matrix_A12, &matrix_A21, &matrix_A22, source_rad, sing_test, &(interf.intervals[j].norm_rad), &(interf.intervals[j].norm_vert), source_norm_vert_3, interf.mid_norm_rad[i]);

	      //Perform the Gauss-Legendre integration
	      coeffs[i][j] = 0.0;
	      coeffs[i][j + interf.n_int] = 0.0;

	      coeffs[i + interf.n_int][j] = 0.0;
	      coeffs[i + interf.n_int][j + interf.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i][j] += (matrix_A11[k] + h[k]) * Gauss_int_wts[k];
		  coeffs[i][j + interf.n_int] += matrix_A12[k] * Gauss_int_wts[k];

		  coeffs[i + interf.n_int][j] += matrix_A21[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + interf.n_int] += matrix_A22[k] * Gauss_int_wts[k];
		}

	      coeffs[i][j] = interf.intervals[j].width * coeffs[i][j] / 2.0;
	      coeffs[i][j + interf.n_int] = interf.intervals[j].width * coeffs[i][j + interf.n_int] / 2.0;

	      coeffs[i + interf.n_int][j] = interf.intervals[j].width * coeffs[i + interf.n_int][j] / 2.0;
	      coeffs[i + interf.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + interf.n_int][j] / 2.0;

	      if (j == i)
		{
		  coeffs[i][j] -= (1.0 + viscos_rat) / 2.0;
		  coeffs[i + interf.n_int][j + interf.n_int] -= (1.0 + viscos_rat) / 2.0;
		}

	    }
	}

      //Loop over intervals on the sphere
      for (int j = 0; j < sphere.n_int; j++)
	{
	  //Consider the case that the source point is on axis
	  if (i == 0)
	    {
	      B_axisource(source_vert, &(sphere.intervals[j].vert), &(sphere.intervals[j].rad), &matrix_B21, &matrix_B22);

	      //Perform the Gauss-Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      coeffs[i][j+ 2 * interf.n_int] = 0.0;
	      coeffs[i][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      coeffs[i + interf.n_int][j + 2 * interf.n_int] = 0.0;
	      coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + interf.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + interf.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
	    }

	  else //For the case that the source point is not on axis
	    {
	      sing_test = 0;
	      B(&(sphere.intervals[j].rad), source_vert, &(sphere.intervals[j].vert), &matrix_B11, &matrix_B12, &matrix_B21, &matrix_B22, source_rad, source_rad_2, sing_test, &g1, &g2);

	      //Perform the Gauss-Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      coeffs[i][j+ 2 * interf.n_int] = 0.0;
	      coeffs[i][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      coeffs[i + interf.n_int][j + 2 * interf.n_int] = 0.0;
	      coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i][j+ 2 * interf.n_int] += matrix_B11[k] * Gauss_int_wts[k];
		  coeffs[i][j + 2 * interf.n_int + sphere.n_int] += matrix_B12[k] * Gauss_int_wts[k];

		  coeffs[i + interf.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i][j+ 2 * interf.n_int] = interf.intervals[j].width * coeffs[i][j+ 2 * interf.n_int] / 2.0;
	      coeffs[i][j + 2 * interf.n_int + sphere.n_int] = interf.intervals[j].width * coeffs[i][j + 2 * interf.n_int + sphere.n_int] / 2.0;

	      coeffs[i + interf.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
	    }  
	}
      //Complete the last column of the matrix
      coeffs[i][2 * (interf.n_int + sphere.n_int)] = 0.0;
      coeffs[i + interf.n_int][2 * (interf.n_int + sphere.n_int)] = 0.0;
    }


  //Loop over the source points on the sphere	      
  for (int i = 0; i < interf.n_int; i++)
    {
      if (i == 0)
	{
	  source_rad = 0.0;
	  source_vert = sphere.height + 1.0;

	  source_rad_2 = 0.0;
	}
      else if (i == sphere.n_int - 1)
	{
	  source_rad = 0.0;
	  source_vert = sphere.height - 1.0;

	  source_rad_2 = 0.0;
	}
      else
	{
	  source_rad = sin(sphere.midpoints[i]);
	  source_vert = sphere.height + cos(sphere.midpoints[i]);
	}

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  //For the case that the source point is on axis
	  if (i == 0 || i == sphere.n_int - 1)
	    {
	      A_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), viscos_rat, &matrix_A21, &matrix_A22);

	      //Perform the Gauss-Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      coeffs[i + 2 * interf.n_int][j] = 0.0;
	      coeffs[i + 2 * interf.n_int][j + interf.n_int] = 0.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j] = 0.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j] += matrix_A21[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] += matrix_A22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + 2* interf.n_int + sphere.n_int][j + interf.n_int] / 2.0;
	    }

	  //For the case that the source point is not on axis
	  else 
	    {
	      sing_test = 0;

	      A(&(interf.intervals[j].rad), source_vert, &(interf.intervals[j].vert), source_rad_2, viscos_rat, interf.midpoints[i], &(interf.intervals[j].arc), &h, &matrix_A11, &matrix_A12, &matrix_A21, &matrix_A22, source_rad, sing_test, &(interf.intervals[j].norm_rad), &(interf.intervals[j].norm_vert), source_norm_vert_3, interf.mid_norm_rad[i]);

	      //Perform the Gauss-Legendre integration
	      coeffs[i + 2 * interf.n_int][j] = 0.0;
	      coeffs[i + 2 * interf.n_int][j + interf.n_int] = 0.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j] = 0.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int][j] += matrix_A11[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int][j + interf.n_int] += matrix_A12[k] * Gauss_int_wts[k];

		  coeffs[i + 2 * interf.n_int + sphere.n_int][j] += matrix_A21[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] += matrix_A22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int][j] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int][j] / 2.0;
	      coeffs[i + 2 * interf.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int][j + interf.n_int] / 2.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j] / 2.0;

	    }
	
	}

      //Loop over intervals on the sphere
      for (int j = 0; j < sphere.n_int; j++) 
	{
	  //For the case that the source point is at theta = 0 or PI
	  if (i == 0 || i == sphere.n_int - 1)
	    {
	      B_axisource(source_vert, &(sphere.intervals[j].vert), &(sphere.intervals[j].rad), &matrix_B21, &matrix_B22);

	      //Perform the Gauss Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] = 0.0;
	      coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = 0.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
	    }

	  //For the case that the source point is in the region of integration for an intermediate interval
	  if (i == j)
	    {
	      sing_test = 1;
	    }
	  else
	    {
	      sing_test = 0;
	    }

	  B(&(sphere.intervals[j].rad), source_vert, &(sphere.intervals[j].vert), &matrix_B11, &matrix_B12, &matrix_B21, &matrix_B22, source_rad, source_rad_2, sing_test, &g1, &g2);

	  if (i == j)
	    {
	      //Perform the Gauss Legendre integration
	      coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] = 0.0;
	      coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = 0.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;
	      
	      for (int k = 0; k < 4; k++)
		{
		  theta_diff = sphere.intervals[j].theta[k] - sphere.midpoints[i];
		  theta_diff_2 = theta_diff * theta_diff;

		  coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] += (matrix_B11_reg[k] + g1[k] * log(4.0 * source_rad_2 * diff / (sum * theta_diff_2)) + 2.0 * (g1[k] + ellip1_b0 / (4.0 * PI * source_rad)) * log(fabs(theta_diff) / (2.0 * source_rad))) * Gauss_int_wts[k];

		  coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B12[k] * Gauss_int_wts[k];
		  
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];

		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] += (matrix_B22_reg[k] + g2[k] * log(4.0 * source_rad_2 * diff / (sum * theta_diff_2)) + 2.0 * (g2[k] + ellip1_b0 / (4.0 * PI * source_rad)) * log(fabs(theta_diff) / (2.0 * source_rad))) * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] / 2.0 - ellip1_b0 / (2.0 * (sphere.n_int - 1.0) * source_rad) * (log(PI / (4.0 * (sphere.n_int - 1.0) * source_rad)) - 1.0);

	      coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] / 2.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0 - ellip1_b0 / (2.0 * (sphere.n_int - 1.0) * source_rad) * (log(PI / (4.0 * (sphere.n_int - 1.0) * source_rad)) - 1.0);
	    }

	  else
	    {
	      //Perform the Gauss Legendre integration
	      coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] = 0.0;
	      coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = 0.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = 0.0;

	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] += matrix_B11[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B21[k] * Gauss_int_wts[k];

		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] += matrix_B12[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] / 2.0;
	      coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
	    }
	}

      //Complete the last column of the matrix
      coeffs[i + 2 * interf.n_int][2 * (interf.n_int + sphere.n_int)] = 0.0;
      coeffs[i + 2 * interf.n_int + sphere.n_int][2 * (interf.n_int + sphere.n_int)] = 1.0;
    }


  //Complete the last row of the matrix
  for (int j = 0; j < 2 * interf.n_int + sphere.n_int; j++)
    {
      coeffs[2 * (interf.n_int + sphere.n_int)][j] = 0.0;
    }

  for (int j = 2 * interf.n_int + sphere.n_int; j < 2 * (interf.n_int + sphere.n_int); j++)
    {
      coeffs[2 * (interf.n_int + sphere.n_int)][j] = 1.0;
    }

  coeffs[2 * (interf.n_int + sphere.n_int)][2 * (interf.n_int + sphere.n_int)] = 0.0;

  //Fill up the known vector as defined in the notes

  //Start looping over the source points on the interface
  for (int i = 0; i < interf.n_int; i++)
    {
      Interf_source(&source_rad, &source_vert, interf.mid_rad[i], interf.mid_vert[i], &source_rad_2);

      known[i] = 0.0;
      known[i + interf.n_int] = 0.0;

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  temp1[j] = 0.0;
	  temp2[j] = 0.0;

	  //For the case that the source point is on axis
	  if (i == 0)
	    {
	      C_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), &(interf.intervals[j].div_norm), &(temp2[j]), &Gauss_int_wts, bond, mdr);

	      known[i + interf.n_int] += interf.intervals[j].width * temp2[j] / 2.0;
	    }

	  //For the case that the integration point is in the interval
	  else if (i == j)
	    {
	      //Loop over the integration points in the interval and find the values of the integrands
	      for (int k = 0; k < 4; k++)
		{
		  Vert_diff(&vert_diff, source_vert, interf.intervals[j].vert[k], &vert_diff_2);

		  alpha_2 = Alpha_2(source_rad_2, interf.intervals[j].rad[k] * interf.intervals[j].rad[k], vert_diff_2);
		  beta_2 = Beta_2(source_rad, interf.intervals[j].rad[k]);

		  beta_4 = beta_2 * beta_2;

		  sum = alpha_2 + beta_2;
		  sum_half = sqrt(sum);

		  diff = alpha_2 - beta_2;

		  prefac[k] = C_prefac(interf.intervals[j].div_norm[k], bond, interf.intervals[j].vert[k], mdr, beta_2, sum_half);

		  Stand_ellip(beta_2, sum, &comp_param, &ellip1[k], &ellip2[k]);
		  ellip1_reg[k] = Ellip1_reg(comp_param);


		  vector_C1_reg[k] = Vector_C1_reg(vert_diff, interf.intervals[j].norm_rad[k], interf.intervals[j].rad[k], interf.intervals[j].norm_vert[k], beta_4, alpha_2, vert_diff_2, source_rad, beta_2, prefac[k], ellip1[k], ellip2[k], diff, ellip1_reg[k]);

		  vector_C2_reg[k] = Vector_C2_reg(source_rad, vert_diff, interf.intervals[j].norm_rad[k], ellip1[k], alpha_2, interf.intervals[j].rad[k], beta_2, interf.intervals[j].norm_vert[k], ellip2[k], diff, ellip1_reg[k], prefac[k]);

		  j1[k] = (prefac[k], alpha_2, interf.intervals[j].norm_rad[k]);
		  j2[k] = (prefac[k], interf.intervals[j].norm_vert[k]);

		  temp1[j] += (vector_C1_reg[k] + j1[k] * log(fabs(arc_diff) * sum / (2.0 * source_rad * diff)) + (9.0 * (interf.mid_div_norm[j] - bond* interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_rad[j] / (8.0 * PI * mdr * bond * source_rad) - j1[k]) * log(fabs(arc_diff) / 2.0 * source_rad)) * Gauss_int_wts[k];

		  temp2[j] += (vector_C2_reg[k] + j2[k] * log(fabs(arc_diff) * sum / (2.0 * source_rad * diff)) + (9.0 * (interf.mid_div_norm[j] - bond* interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_vert[j] / (8.0 * PI * mdr * bond * source_rad) - j2[k]) * log(fabs(arc_diff) / 2.0 * source_rad)) * Gauss_int_wts[k];
		}

	      known[i] += interf.intervals[j].width * temp1[j] / 2.0 - 9.0 * (interf.mid_div_norm[j] - bond* interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_rad[j] * interf.midpoints[interf.n_int - 1] / (8.0 * PI * mdr * bond * source_rad * (interf.n_int - 1.0)) * (log(interf.midpoints[interf.n_int - 1] / (4.0 * source_rad * (interf.n_int - 1.0))) - 1.0);

	      known[i + interf.n_int] += interf.intervals[j].width * temp2[j] / 2.0 - 9.0 * (interf.mid_div_norm[j] - bond* interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_rad[j] * interf.midpoints[interf.n_int - 1] / (8.0 * PI * mdr * bond * source_rad * (interf.n_int - 1.0)) * (log(interf.midpoints[interf.n_int - 1] / (4.0 * source_rad * (interf.n_int - 1.0))) - 1.0);
	    }	      

	  //For the case that the source point is off axis and the integral is regular
	  else
	    {
	      //Loop over the integration points in the interval and find the values of the integrands
	      for (int k = 0; k < 4; k++)
		{
		  //Define temporary variables
		  Vert_diff(&vert_diff, source_vert, interf.intervals[j].vert[k], &vert_diff_2);

		  alpha_2 = Alpha_2(source_rad_2, interf.intervals[j].rad[k] * interf.intervals[j].rad[k], vert_diff_2);
		  beta_2 = Beta_2(source_rad, interf.intervals[j].rad[k]);

		  beta_4 = beta_2 * beta_2;

		  sum = alpha_2 + beta_2;
		  sum_half = sqrt(sum);

		  diff = alpha_2 - beta_2;

		  prefac[k] = C_prefac(interf.intervals[j].div_norm[k], bond, interf.intervals[j].vert[k], mdr, beta_2, sum_half);

		  Stand_ellip(beta_2, sum, &comp_param, &ellip1[k], &ellip2[k]);

		  vector_C1[k] = Vector_C1(alpha_2, vert_diff_2, interf.intervals[j].norm_rad[k], interf.intervals[j].rad[k], vert_diff, interf.intervals[j].norm_vert[k], beta_4, source_rad, prefac[k], ellip1[k], ellip2[k], diff, beta_2);

		  vector_C2[k] = Vector_C2(beta_2, interf.intervals[j].norm_vert[k], interf.intervals[j].rad[k], vert_diff, interf.intervals[j].norm_rad[k], source_rad, alpha_2, prefac[k], ellip1[k], ellip2[k], diff);

		  temp1[j] += vector_C1[k] * Gauss_int_wts[k];

		  temp2[j] += vector_C2[k] * Gauss_int_wts[k];
		}

	      known[i] += interf.intervals[j].width * temp1[j] / 2.0; 
	      known[i + interf.n_int] += interf.intervals[j].width * temp1[j] / 2.0; 
	    }
	}
    }

  //Resize the temp vectors
  temp1.resize(sphere.n_int);
  temp2.resize(sphere.n_int);
  //Loop over the source points on the sphere
  for (int i = 0; i < sphere.n_int; i++)
    {
      if (i == 0)
	{
	  source_rad = 0.0;
	  source_vert = sphere.height + 1.0;

	  source_rad_2 = 0.0;
	}
      else if (i == sphere.n_int - 1)
	{
	  source_rad = 0.0;
	  source_vert = sphere.height - 1.0;

	  source_rad_2 = 0.0;
	}
      else
	{
	  source_rad = sin(sphere.midpoints[i]);
	  source_vert = sphere.height + cos(sphere.midpoints[i]);
	}

      source_rad_2 = source_rad * source_rad;

      known[i + 2 * interf.n_int] = 0.0;
      known[i + 2 * interf.n_int + sphere.n_int] = 0.0;

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  temp1[j] = 0.0;
	  temp2[j] = 0.0;

	  //For the case that the source point is on axis
	  if (i == 0 || i == sphere.n_int - 1)
	    {
	      C_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), &(interf.intervals[j].div_norm), &(temp2[j]), &Gauss_int_wts, bond, mdr);

	      known[i + 2 * interf.n_int + sphere.n_int] += interf.intervals[j].width * temp2[j] / 2.0;
	    }

	  //For the case that the source point is off axis and the integral is regular
	  else
	    {
	      //Loop over the integration points in the interval and find the values of the integrands
	      for (int k = 0; k < 4; k++)
		{
		  //Define temporary variables
		  Vert_diff(&vert_diff, source_vert, interf.intervals[j].vert[k], &vert_diff_2);

		  alpha_2 = source_rad * source_rad + interf.intervals[j].rad[k] * interf.intervals[j].rad[k] + vert_diff_2;
		  beta_2 = Beta_2(source_rad, interf.intervals[j].rad[k]);

		  beta_4 = beta_2 * beta_2;

		  sum = alpha_2 + beta_2;
		  sum_half = sqrt(sum);

		  diff = alpha_2 - beta_2;

		  prefac[k] = C_prefac(interf.intervals[j].div_norm[k], bond, interf.intervals[j].vert[k], mdr, beta_2, sum_half);

		  Stand_ellip(beta_2, sum, &comp_param, &ellip1[k], &ellip2[k]);

		  vector_C1[k] = Vector_C1(alpha_2, vert_diff_2, interf.intervals[j].norm_rad[k], interf.intervals[j].rad[k], vert_diff, interf.intervals[j].norm_vert[k], beta_4, source_rad, prefac[k], ellip1[k], ellip2[k], diff, beta_2);

		  vector_C2[k] = Vector_C2(beta_2, interf.intervals[j].norm_vert[k], interf.intervals[j].rad[k], vert_diff, interf.intervals[j].norm_rad[k], source_rad, alpha_2, prefac[k], ellip1[k], ellip2[k], diff);

		  temp1[j] += vector_C1[k] * Gauss_int_wts[k];

		  temp2[j] += vector_C2[k] * Gauss_int_wts[k];
		}

	      known[i + 2 * interf.n_int] += interf.intervals[j].width * temp1[j] / 2.0; 
	      known[i + 2 * interf.n_int + sphere.n_int] += interf.intervals[j].width * temp1[j] / 2.0; 
	    }
	}

    }

  //Move calculated matrix and vector elements into objects that were passed into function

  for (int i = 0; i < coeffs.size(); i++)
    {
      for (int j = 0; j < coeffs[i].size(); j++)
	{
	  (*matrix)[i][j] = coeffs[i][j];
	}
    }

  for (int i = 0; i < known.size(); i++)
    {
      (*vec)[i] = known[i];
    }
}
