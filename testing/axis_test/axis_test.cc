//Program file to test the evaluation of the components of A and C when the source point is on axi

#include <math.h>
#include <iostream>

#include "../../axisym.h"

using std::cout;
using std::endl;

int main()
{
  //Set dimensionless variables to be used
  double viscos_rat = 2.0;
  double bond = 2.0;
  double mdr = 2.0;

  //Set location of source and position points
  double source_vert = 3.0;

  double pos_rad = 4.0;
  double pos_vert = 5.0;

  //Calculate quantities that depend on these locations
  double vert_diff = source_vert - pos_vert;
  double vert_diff_2 = vert_diff * vert_diff;
  double vert_diff_3 = vert_diff_2 * vert_diff;

  double pos_rad_2 = pos_rad * pos_rad;

  double alpha_2 = pos_rad_2 + vert_diff_2;

  double alpha = sqrt(alpha_2);

  double alpha_5 = alpha_2 * alpha_2 * alpha;

  //Set divergence of normal vector
  double div_norm = 0.5;

  //Calculate components of A and output results for testing
  cout << "Quantities A21, A22 and C2 when the source point is on axis" << endl;

  double matrix_A21 = Matrix_A21_axisource(viscos_rat, vert_diff_2, pos_rad, alpha_5);
  cout << "A21 = " << matrix_A21 << endl;

  double matrix_A22 = Matrix_A22_axisource(viscos_rat, vert_diff_3, alpha_5);
  cout << "A22 = " << matrix_A22 << endl;

  double vector_C2 = Vector_C2_axisource(div_norm, bond, pos_vert, mdr, alpha, vert_diff_2, alpha_2);
  cout << "C2 = " << vector_C2 << endl;

  return 0;
}
