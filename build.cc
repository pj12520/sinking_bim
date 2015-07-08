//File containing functions relating to the construction of the linear system

#include <vector>
#include <math.h>

#include <iostream> //Used for debugging purposes only. Can be removed when program is operational
#include <typeinfo> //Used for debugging purposes only. Can be removed when program is operational
#include "object.h"
#include "const.h"
#include "axisym.h"
#include "ellip.h"
#include "build_sup.h"

using std::vector;
using std::iterator;

using std::cout; //Using for debugging purposes only. Can be removed when program is operational
using std::endl;

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
  vector<vector<double> > coeffs((*matrix).size() + 3);
  vector<double> known((*vec).size() + 3);
  for (int i = 0; i < coeffs.size(); i++)
    {
      coeffs[i].resize((*matrix).size() + 3);
    }

  //Initialise all elements with zeros
  for (int i = 0; i < coeffs.size(); i++)
    {
      for (int j = 0; j < coeffs[i].size(); j++)
	{
	  coeffs[i][j] = 0.0;
	}
      known[i] = 0.0;
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

  double sum; //alpha_2 + beta_2
  double diff; //alpha_2 - beta_2

  double theta_diff; //Separation in theta between source and integration point when they are both in the same interval on the sphere
  double theta_diff_2; //theta_diff squared

  double source_norm_vert_3; //Cube of the vertical component of the normal vector at the source point

  //Variable passed to supplementary functions to indicate if integrand is regular or not
  int sing_test;

  //Vectors to contain temporary values of the integrands at the integration points
  vector<double> matrix_A11(4);
  vector<double> matrix_A12(4);
  vector<double> matrix_A21(4);
  vector<double> matrix_A22(4);

  vector<double> h(4);

  vector<double> matrix_B11(4);
  vector<double> matrix_B12(4);
  vector<double> matrix_B21(4);
  vector<double> matrix_B22(4);

  vector<double> matrix_B11_reg(4);
  vector<double> matrix_B22_reg(4);

  vector<double> g1(4);
  vector<double> g2(4);

  vector<double> temp1(interf.n_int); //Vectors to store temporary values used in the calculation of the known vector. Needs to be resized in the calculation
  vector<double> temp2(interf.n_int); 

  //Fill up the matrix as defined in the notes

  //Start looping over the source points on the interface
  for (int i = 0; i < interf.n_int; i++)
    {
      Interf_source(&source_rad, &source_vert, interf.mid_rad[i], interf.mid_vert[i], &source_rad_2);
      source_norm_vert_3 = interf.mid_norm_vert[i] * interf.mid_norm_vert[i] * interf.mid_norm_vert[i];

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  //For the case that the source point is on axis
	  if (i == 0)
	    {
	      A_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), viscos_rat, &matrix_A21, &matrix_A22, &interf.intervals[j].norm_rad, &interf.intervals[j].norm_vert);
	      //Perform the Gauss-Legendre integration
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + interf.n_int][j] -= matrix_A21[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + interf.n_int] -= matrix_A22[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
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
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i][j] -= (matrix_A11[k] + h[k]) * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  coeffs[i][j + interf.n_int] -= matrix_A12[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  
		  coeffs[i + interf.n_int][j] -= matrix_A21[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + interf.n_int] -= matrix_A22[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		}

	      coeffs[i][j] = interf.intervals[j].width * coeffs[i][j] / 2.0;
	      coeffs[i][j + interf.n_int] = interf.intervals[j].width * coeffs[i][j + interf.n_int] / 2.0;

	      coeffs[i + interf.n_int][j] = interf.intervals[j].width * coeffs[i + interf.n_int][j] / 2.0;
	      coeffs[i + interf.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + interf.n_int][j + interf.n_int] / 2.0;

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
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + interf.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + interf.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
	      //	      cout << i << " " << j << " " << coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] << " " << endl;
	    }

	  else //For the case that the source point is not on axis
	    {
	      sing_test = 0;
	      B(&(sphere.intervals[j].rad), source_vert, &(sphere.intervals[j].vert), &matrix_B11, &matrix_B12, &matrix_B21, &matrix_B22, source_rad, source_rad_2, sing_test, &g1, &g2, &sum, &diff, &matrix_B11_reg, &matrix_B22_reg);

	      //Perform the Gauss-Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i][j + 2 * interf.n_int] += matrix_B11[k] * Gauss_int_wts[k];
		  coeffs[i][j + 2 * interf.n_int + sphere.n_int] += matrix_B12[k] * Gauss_int_wts[k];

		  coeffs[i + interf.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		  coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i][j+ 2 * interf.n_int] / 2.0;
	      coeffs[i][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i][j + 2 * interf.n_int + sphere.n_int] / 2.0;

	      coeffs[i + interf.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;

	    }  
	}

    }

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

	  source_rad_2 = source_rad * source_rad;
	}

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  //For the case that the source point is on axis
	  if (i == 0 || i == sphere.n_int - 1)
	    {
	      A_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), viscos_rat, &matrix_A21, &matrix_A22, &(interf.intervals[j].norm_rad), &(interf.intervals[j].norm_vert));

	      //Perform the Gauss-Legendre integration (Riley Hobson and Bence 2006 page 1006)
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j] -= matrix_A21[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] -= matrix_A22[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] / 2.0;
	    }

	  //For the case that the source point is not on axis
	  else 
	    {
	      sing_test = 0;

	      A(&(interf.intervals[j].rad), source_vert, &(interf.intervals[j].vert), source_rad_2, viscos_rat, interf.midpoints[i], &(interf.intervals[j].arc), &h, &matrix_A11, &matrix_A12, &matrix_A21, &matrix_A22, source_rad, sing_test, &(interf.intervals[j].norm_rad), &(interf.intervals[j].norm_vert), source_norm_vert_3, interf.mid_norm_rad[i]);

	      //Perform the Gauss-Legendre integration
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int][j] -= matrix_A11[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int][j + interf.n_int] -= matrix_A12[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  //		  cout << i << '\t' << j << '\t' << k << '\t' << matrix_A11[k] << '\t' << matrix_A12[k] << endl; 
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j] -= matrix_A21[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] -= matrix_A22[k] * interf.intervals[j].rad[k] * Gauss_int_wts[k];
		}
	      //      cout << i << '\t' << j << '\t' << coeffs[i + 2 * interf.n_int][j] << '\t' << coeffs[i + 2 * interf.n_int][j + interf.n_int] << endl;

	      coeffs[i + 2 * interf.n_int][j] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int][j] / 2.0;
	      coeffs[i + 2 * interf.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int][j + interf.n_int] / 2.0;

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] = interf.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + interf.n_int] / 2.0;
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
	      for (int k = 0; k < 4; k++)
		{
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		}

	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] / 2.0;
	      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
	    }
	  else
	    {
	      //For the case that the source point is in the region of integration for an intermediate interval
	      if (i == j)
		{
		  sing_test = 1;
		}
	      else
		{
		  sing_test = 0;
		}

	      B(&(sphere.intervals[j].rad), source_vert, &(sphere.intervals[j].vert), &matrix_B11, &matrix_B12, &matrix_B21, &matrix_B22, source_rad, source_rad_2, sing_test, &g1, &g2, &sum, &diff, &matrix_B11_reg, &matrix_B22_reg);

	      if (i == j)
		{
		  //Perform the Gauss Legendre integration
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
		  for (int k = 0; k < 4; k++)
		    {
		      coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] += matrix_B11[k] * Gauss_int_wts[k];
		      coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B12[k] * Gauss_int_wts[k];

		      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] += matrix_B21[k] * Gauss_int_wts[k];
		      coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] += matrix_B22[k] * Gauss_int_wts[k];
		    }

		  coeffs[i + 2 * interf.n_int][j+ 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int] / 2.0;
		  coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;

		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int] / 2.0;
		  coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] = sphere.intervals[j].width * coeffs[i + 2 * interf.n_int + sphere.n_int][j + 2 * interf.n_int + sphere.n_int] / 2.0;
		}
	    }
	}

      //Complete the last column of the matrix
      coeffs[i + 2 * interf.n_int + sphere.n_int][2 * (interf.n_int + sphere.n_int)] = -1.0;
    }

  //Complete the last row of the matrix
  for (int j = 2 * interf.n_int + sphere.n_int; j < 2 * (interf.n_int + sphere.n_int); j++)
    {
      coeffs[2 * (interf.n_int + sphere.n_int)][j] = sphere.intervals[j - (2 * interf.n_int + sphere.n_int)].width;
    }

  //Fill up the known vector as defined in the notes

  //Start looping over the source points on the interface
  for (int i = 0; i < interf.n_int; i++)
    {
      Interf_source(&source_rad, &source_vert, interf.mid_rad[i], interf.mid_vert[i], &source_rad_2);

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  temp1[j] = 0.0;
	  temp2[j] = 0.0;

	  //For the case that the source point is on axis
	  if (i == 0)
	    {
	      C_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), &(interf.intervals[j].div_norm), &(temp2[j]), &Gauss_int_wts, bond, mdr, &(interf.intervals[j].norm_vert), &(interf.intervals[j].norm_rad));

	      known[i + interf.n_int] -= interf.intervals[j].width * temp2[j] / 2.0;
	    }
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

	      C(source_vert, source_rad_2, source_rad, &(interf.intervals[j].vert), &(interf.intervals[j].rad), &(interf.intervals[j].div_norm), bond, mdr, sing_test, &(interf.intervals[j].norm_rad), &(interf.intervals[j].norm_vert), &(temp1[j]), &(temp2[j]), interf.mid_div_norm[i], interf.mid_vert[i], interf.mid_norm_rad[i], interf.mid_norm_vert[i], &(interf.intervals[j].arc), interf.midpoints[i], &Gauss_int_wts);

	      if (i == j)
		{
		  //		  known[i] -= interf.intervals[j].width * temp1[j] / 2.0 - 9.0 * (interf.mid_div_norm[j] - bond * interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_rad[j] * interf.midpoints[interf.n_int - 1] / (8.0 * PI * mdr * bond * (interf.n_int - 1.0)) * (log(interf.midpoints[interf.n_int - 1] / (4.0 * source_rad * (interf.n_int - 1.0))) - 1.0);

		  known[i] -= interf.intervals[j].width * temp1[j] / 2.0 - 9.0 * (interf.mid_div_norm[j] - bond * interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_rad[j] * interf.intervals[j].width / (4.0 * PI * mdr * bond) * (log(interf.intervals[j].width / (4.0 * source_rad)) - 1.0);

		  //		  known[i + interf.n_int] -= interf.intervals[j].width * temp2[j] / 2.0 - 9.0 * (interf.mid_div_norm[j] - bond * interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_vert[j] * interf.midpoints[interf.n_int - 1] / (8.0 * PI * mdr * bond * (interf.n_int - 1.0)) * (log(interf.midpoints[interf.n_int - 1] / (4.0 * source_rad * (interf.n_int - 1.0))) - 1.0);

		  known[i + interf.n_int] -= interf.intervals[j].width * temp2[j] / 2.0 - 9.0 * (interf.mid_div_norm[j] - bond * interf.mid_vert[j]) * ellip1_b0 * interf.mid_norm_vert[j] * interf.intervals[j].width / (4.0 * PI * mdr * bond) * (log(interf.intervals[j].width / (4.0 * source_rad)) - 1.0);
		}	      

	      else
		{
		  known[i] -= interf.intervals[j].width * temp1[j] / 2.0; 
		  known[i + interf.n_int] -= interf.intervals[j].width * temp2[j] / 2.0; 
		}
	    }
	}
    }

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

      //Loop over the intervals on the interface
      for (int j = 0; j < interf.n_int; j++)
	{
	  temp1[j] = 0.0;
	  temp2[j] = 0.0;

	  //For the case that the source point is on axis
	  if (i == 0 || i == sphere.n_int - 1)
	    {
	      C_axisource(source_vert, &(interf.intervals[j].vert), &(interf.intervals[j].rad), &(interf.intervals[j].div_norm), &(temp2[j]), &Gauss_int_wts, bond, mdr, &(interf.intervals[j].norm_vert), &(interf.intervals[j].norm_rad));

	      known[i + 2 * interf.n_int + sphere.n_int] -= interf.intervals[j].width * temp2[j] / 2.0;
	    }

	  //For the case that the source point is off axis and the integral is regular
	  else
	    {
	      sing_test = 0;

	      C(source_vert, source_rad_2, source_rad, &(interf.intervals[j].vert), &(interf.intervals[j].rad), (&interf.intervals[j].div_norm), bond, mdr, sing_test, (&interf.intervals[j].norm_rad), (&interf.intervals[j].norm_vert), &(temp1[j]), &(temp2[j]), interf.mid_div_norm[j], interf.mid_vert[j], interf.mid_norm_rad[j], interf.mid_norm_vert[j], &(interf.intervals[j].arc), interf.midpoints[j], &Gauss_int_wts);

	      known[i + 2 * interf.n_int] -= interf.intervals[j].width * temp1[j] / 2.0; 
	      known[i + 2 * interf.n_int + sphere.n_int] -= interf.intervals[j].width * temp2[j] / 2.0; 
	    }
	}

    }
  //Final row of known vector
  known[2 * (interf.n_int + sphere.n_int)] = -3.0;

  //Now remove degenerate columns from matrix elements from vector
  vector<unsigned> to_delete(3);
  to_delete[0] = 0;
  to_delete[1] = 2 * interf.n_int - 1;
  to_delete[2] = 2 * interf.n_int + sphere.n_int - 3;

  for (int i = 0; i < 3; i++)
    {
      if (coeffs.size() > to_delete[i])
	{
	  coeffs.erase(coeffs.begin() + to_delete[i]);
	}

      if (known.size() > to_delete[i])
	{
	  known.erase(known.begin() + to_delete[i]);
	}

      for (int j = 0; j < coeffs.size(); j++)
	{
	  if (coeffs[j].size() > to_delete[i])
	    {
	      coeffs[j].erase(coeffs[j].begin() + to_delete[i]);
	    }
	}
    }

  //Remove degenerate rows
  /*
  vector<vector<double > > hold1(coeffs.size() - 1);
  vector<vector<double > > hold2(coeffs.size() - 2);
  vector<vector<double > > hold3(coeffs.size() - 3);

  for (int i = 0; i < hold1.size(); i++)
    {
      hold1[i].resize(coeffs[i + 1].size());
      
      hold1[i] = coeffs[i + 1]; 

    }

  for (int i = 0; i < 2 * interf.n_int - 1; i++)
    {
      hold2[i].resize(hold1[i].size());
      
      hold2[i] = hold1[i]; 
    }

  for (int i = 2 * interf.n_int - 1; i < hold2.size(); i++)
    {
      hold2[i].resize(hold1[i + 1].size());
      
      hold2[i] = hold1[i + 1]; 
    }

  for (int i = 0; i < 2 * interf.n_int + sphere.n_int - 3; i++)
    {
      hold3[i].resize(hold2[i].size());
      
      hold3[i] = hold2[i]; 
    }

  for (int i = 2 * interf.n_int + sphere.n_int - 3; i < hold3.size(); i++)
    {
      hold3[i].resize(hold2[i + 1].size());
      
      hold3[i] = hold2[i + 1]; 
    }
  */
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
