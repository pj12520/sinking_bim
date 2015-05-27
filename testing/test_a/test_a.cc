//Program file to test the evaluation of a1, a2, etc. 

#include <math.h>
#include <iostream>

#include "../../axisym.h"

using std::cout;
using std::endl;

int main()
{
  double viscos_rat = 2.0;

  double source_rad = 2.0;
  double source_vert = 3.0;

  double pos_rad = 4.0;
  double pos_vert = 5.0;

  double source_rad_2 = source_rad * source_rad;
  double pos_rad_2 = pos_rad * pos_rad;

  double pos_rad_4 = pos_rad_2 * pos_rad_2;

  double vert_diff = source_vert - pos_vert;
  double vert_diff_2 = vert_diff * vert_diff;
  double vert_diff_3 = vert_diff_2 * vert_diff;

  double alpha_2 = source_rad_2 + pos_rad_2 + vert_diff_2;
  double beta_2 = 2.0 * source_rad * pos_rad;

  double alpha_4 = alpha_2 * alpha_2;
  double beta_4 = beta_2 * beta_2;

  double alpha_6 = alpha_4 * alpha_2;
  double beta_6 = beta_4 * beta_2;

  double alpha_8 = alpha_6 * alpha_2;
  double beta_8 = beta_6 * beta_2;

  double sum = alpha_2 + beta_2;
  double sum_3_2 = pow(sum, 1.5);

  double diff = alpha_2 - beta_2;
  double diff_2 = diff * diff;

  double a1 = A1(viscos_rat, sum_3_2, diff, beta_4, source_rad, alpha_2, alpha_4, source_rad_2, pos_rad_2, pos_rad, beta_2);
  cout << "a1 = " << a1 << endl;

  double a2 = A2(viscos_rat, vert_diff, alpha_4, beta_4, alpha_2, vert_diff_2, sum_3_2, diff, beta_2);
  cout << "a2 = " << a2 << endl;

  double a3 = A3(viscos_rat, sum_3_2, diff_2, beta_4, source_rad, alpha_8, alpha_4, beta_8, pos_rad_2, source_rad_2, alpha_2, pos_rad, beta_2);
  cout << "a3 = " << a3 << endl;

  double a4 = A4(viscos_rat, vert_diff, beta_2, alpha_2, alpha_4, beta_4, vert_diff_2, sum_3_2, diff);
  cout << "a4 = " << a4 << endl;

  double a6 = A6(viscos_rat, vert_diff_2, source_rad_2, alpha_2, sum_3_2, diff, source_rad); 
  cout << "a6 = " << a6 << endl;

  double a8 = A8(viscos_rat, vert_diff_2, alpha_4, beta_4, source_rad_2, alpha_2, sum_3_2, diff, source_rad);
  cout << "a8 = " << a8 << endl;

  double a9 = A9(viscos_rat, vert_diff, alpha_4, beta_4, pos_rad_2, alpha_2, pos_rad_4, sum_3_2, diff);
  cout << "a9 = " << a9 << endl;

  double a10 = A10(viscos_rat, vert_diff_2, alpha_2, pos_rad_2, sum_3_2, diff, pos_rad);
  cout << "a10 = " << a10 << endl;

  double a11 = A11(viscos_rat, vert_diff, alpha_6, alpha_2, beta_4, pos_rad_2, pos_rad_4, sum_3_2, diff_2, alpha_4);
  cout << "a11 = " << a11 << endl;

  double a12 = A12(viscos_rat, vert_diff_2, pos_rad_2, alpha_2, alpha_4, beta_4, sum_3_2, diff_2, pos_rad);
  cout << "a12 = " << a12 << endl;

  double a14 = A14(viscos_rat, vert_diff_3, sum_3_2, diff);
  cout << "a14 = " << a14 << endl;

  double a16 = A16(viscos_rat, alpha_2, sum_3_2, diff_2, vert_diff_3);
  cout << "a16 = " << a16 << endl;

  return 0;
}
