//Program file containing functions definitions relating to geometry

#include <math.h>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <iomanip> //Inlcuded for debugging purposes only. Remove when program is functional
#include <fstream> //Included for debugging purposes only

#include "const.h"
#include "interp_1d.h"
#include "dfridr.h"

using std::vector;
using std::cout;
using std::endl;

using std::setw; //Using for debugging purposes only. Remove when program is functional
using std::ofstream; //Using for debugging purposes only

//Function to calculate hypotenuse of a right angled triangle given the length of the other sides.
double Pythag(double side1, double side2)
{
  double hyp;

  hyp = pow(pow(side1,2.0) + pow(side2,2.0),0.5);

  return hyp;
}


//Function to calculate the components of the normal vector and it's divergence at a point along the interface
void Normal(Spline_interp rad, Spline_interp height, double arc, double init_step, double *norm_rad, double *norm_vert, double *div_norm, double rad_coord, vector<double>* midpoints, vector<double>* pos_rad, vector<double>* pos_vert, double fit_const0, double fit_const1, double fit_const2, double fit_const3, double arc_max, ofstream& out)
{
  double rad_deriv_error;
  double height_deriv_error;
  double rad_deriv2_error;
  double height_deriv2_error;

  //    double rad_deriv = dfridr_interp(rad, arc, init_step, rad_deriv_error);
  //    double height_deriv = dfridr_interp(height, arc, init_step, height_deriv_error);
  //    double rad_deriv2 = sec_dfridr(rad, arc, init_step, rad_deriv2_error);
  //    double height_deriv2 = sec_dfridr(height, arc, init_step, height_deriv2_error);
  
  double rad_deriv = My_dfridr(&Rad, arc, init_step, rad_deriv_error, rad, fit_const0, fit_const1, arc_max, 0);
  double height_deriv = My_dfridr(&Vert, arc, init_step, height_deriv_error, height, fit_const2, fit_const3, arc_max, 1);
  double rad_deriv2 = My_sec_dfridr(rad, arc, init_step, rad_deriv2_error, fit_const0, fit_const1, arc_max, &Rad, 0);
  double height_deriv2 = My_sec_dfridr(height, arc, init_step, height_deriv2_error, fit_const2, fit_const3, arc_max, &Vert, 1);

  //  out << setw(20) << arc << setw(20) << rad_deriv << setw(20) << height_deriv << setw(20) << rad_deriv2 << setw(20) << height_deriv2 << endl;


    //double rad_deriv = deriv(rad, arc, &(*midpoints), &(*pos_rad));
    //  double height_deriv = deriv(height, arc, &(*midpoints), &(*pos_vert));
    //  double rad_deriv2 = sec_deriv(rad, arc, &(*midpoints)); 
    //  double height_deriv2 = sec_deriv(height, arc, &(*midpoints));

  //    cout << setw(15) << arc << " " << setw(15) << rad_deriv << " " << setw(15) << height_deriv << " " << setw(15) << rad_deriv2 << " " << setw(15) << height_deriv2 << endl;
  *norm_rad = - height_deriv / Pythag(rad_deriv, height_deriv);

  *norm_vert =  rad_deriv / Pythag(rad_deriv, height_deriv);

  *div_norm = 1.0 * (rad_deriv2 * height_deriv - rad_deriv * height_deriv2) / pow(rad_deriv * rad_deriv + height_deriv * height_deriv, 1.5) - height_deriv / (rad_coord * pow(rad_deriv * rad_deriv + height_deriv * height_deriv, 0.5));

  *div_norm = -*div_norm;
  //  cout << arc << " " << rad_coord << " " << *norm_rad << " " << *norm_vert << " " << *div_norm << endl;
}

//Function to calculate the normal components and the divergence of the normal to the interface when it is described by the extrapolated functions
void Normal_fit(double fit_const1, double fit_const2, double fit_const3, double arc, double arc4, double *norm_rad, double *norm_vert, double *div_norm, double rad_coord)
{
  double arc5 = arc4 * arc;
  double arc6 = arc5 * arc;

  double rad_deriv = fit_const1;

  double height_deriv = -3.0 * fit_const2 / arc4 - 4.0 * fit_const3 / arc5;
  double height_deriv2 = 12.0 * fit_const2 / arc5 + 20 * fit_const3 / arc6;

  *norm_rad = - height_deriv / Pythag(rad_deriv, height_deriv);

  *norm_vert =  rad_deriv / Pythag(rad_deriv, height_deriv);

  *div_norm =  -rad_deriv * height_deriv2 / pow(rad_deriv * rad_deriv + height_deriv * height_deriv, 1.5) - height_deriv / (rad_coord * pow(rad_deriv * rad_deriv + height_deriv * height_deriv, 0.5));
}

//Function to calculate divergence of the normal of the interface at each point along it
double Div_norm(Spline_interp rad, Spline_interp height, double arc, double init_step) //When this is called arc must not be within init_step of the boundaries of the interpolation range 
{
  double curve;

  if (arc == 0)
    {
      curve = 0;
    }
  else
    {
      if (arc < init_step) init_step = arc;
      double rad_deriv_error;
      double height_deriv_error;
      double rad_deriv2_error;
      double height_deriv2_error;

      double rad_coord = rad.interp(arc);
      
      double rad_deriv = dfridr_interp(rad, arc, init_step, rad_deriv_error);
      double height_deriv = dfridr_interp(height, arc, init_step, height_deriv_error);
      double rad_deriv2 = sec_dfridr(rad, arc, init_step, rad_deriv2_error);
      double height_deriv2 = sec_dfridr(height, arc, init_step, height_deriv2_error);

      //      curve = height_deriv / (rad_coord * sqrt(pow(rad_deriv,2.0) + pow(height_deriv,2.0))) + (rad_deriv * height_deriv2 - rad_deriv2 * height_deriv) / pow(pow(rad_deriv,2.0) + pow(height_deriv,2.0), 1.5); //From Manga 94 equation A.38

      //      curve = 2.0 * (height_deriv * rad_deriv2 - rad_deriv * height_deriv2)/ pow(pow(rad_deriv, 2.0) + pow(height_deriv, 2.0), 1.5);

      //      curve = 2.0* (height_deriv * rad_deriv * (rad_deriv2 - height_deriv2) * rad_deriv2 * pow(height_deriv, 2.0) - height_deriv2 * pow(rad_deriv, 2.0)) / pow((pow(height_deriv, 2.0) + pow(rad_deriv, 2.0)), 2.5);

      //            curve = 2.0 * (height_deriv * rad_deriv2 - rad_deriv * height_deriv2)/ pow(pow(rad_deriv, 2.0) + pow(height_deriv, 2.0), 1.5);

      curve = 1.0 * (rad_deriv2 * height_deriv - rad_deriv * height_deriv2) / pow(rad_deriv * rad_deriv + height_deriv * height_deriv, 1.5) - height_deriv / (rad_coord * pow(rad_deriv * rad_deriv + height_deriv * height_deriv, 0.5));
    }

  return curve;
}

/*
//Function to calculate a tangent vector to a line at each point along it
void Tangent(vector<double>* rad, vector<double>* height, vector<double>* tangent_rad, vector<double>* tangent_height, int n_eval)
{
  //Impose that on symmetry axis and at infinity tangent vector is horizontal.

  (*tangent_rad)[0] = 1.0;
  (*tangent_height)[0] = 0.0;

  for (int i = 1; i < n_eval - 1 ; i++)
    {
      (*tangent_rad)[i] = ((*rad)[i + 1] - (*rad)[i - 1]) / Pythag((*rad)[i + 1] - (*rad)[i - 1], (*height)[i + 1] - (*height)[i - 1]);
      (*tangent_height)[i] = ((*height)[i + 1] - (*height)[i - 1]) / Pythag((*rad)[i + 1] - (*rad)[i - 1], (*height)[i + 1] - (*height)[i - 1]);
    }

  (*tangent_rad)[n_eval - 1] = 1;
  (*tangent_height)[n_eval - 1] = 0;
}
*/

//Function to calculate the radial component of the normal vector at a point along it
double Normal_rad(Spline_interp rad, Spline_interp height, double arc, double init_step)
{
  double rad_deriv_error;
  double height_deriv_error;

  double rad_deriv = dfridr_interp(rad, arc, init_step, rad_deriv_error);

  if (rad_deriv_error / rad_deriv > 0.001)
    {
      cout << "Radial derivative with respect to arc large " << arc << " " << rad.interp(arc) << " " << init_step << endl;
    }

  double height_deriv = dfridr_interp(height, arc, init_step, height_deriv_error);

  if (height_deriv_error / height_deriv > 0.001)
    {
      cout << "Vertical derivative with respect to arc large " << arc << " " << height.interp(arc) << " " << init_step << endl;
    }

  double norm_rad = - height_deriv / Pythag(rad_deriv, height_deriv);

  return norm_rad;
}

//Function to calculate the vertical component of the normal vector at a point along it
double Normal_height(Spline_interp rad, Spline_interp height, double arc, double init_step)
{
  double rad_deriv_error;
  double height_deriv_error;

  double rad_deriv = dfridr_interp(rad, arc, init_step, rad_deriv_error);

  if (rad_deriv_error / rad_deriv > 0.001)
    {
      cout << "Radial derivative with respect to arc large " << arc << " " << rad.interp(arc) << " " << init_step << endl;
    }

  double height_deriv = dfridr_interp(height, arc, init_step, height_deriv_error);

  if (height_deriv_error / height_deriv > 0.001)
    {
      cout << "Vertical derivative with respect to arc large " << arc << " " << height.interp(arc) << " " << init_step << endl;
    }

  double norm_height =  rad_deriv / Pythag(rad_deriv, height_deriv);

  return norm_height;
}

/*
//Function to rotate a 2D vector in an anticlockwise direction by an angle theta
void Rotate(vector<double>* init_vector, vector<double>* final_vector, double theta)
{
  vector<vector<double> > rotate;

  rotate.resize(2);

  for (int i = 0; i < 2; i++)
    {
      rotate[i].resize(2);
    }

  rotate[0][0] = cos(theta);
  rotate[0][1] = -sin(theta);
  rotate[1][0] = sin(theta);
  rotate[1][1] = cos(theta);

  for (int i = 0; i <2; i++)
    {
      (*final_vector)[i] = 0;

      for (int j = 0; j < 2; j++)
	{
	  (*final_vector)[i] = (*final_vector)[i] + rotate[i][j] * (*init_vector)[j];
	}
    }
}

*/

