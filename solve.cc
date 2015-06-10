#include <vector>
#include <gsl/gsl_linalg.h>

using std::vector;

void Solve(int order, vector<vector<double> >* coeffs, vector<double>* known, vector<double>* unknown)
{
  vector<double> coeffs_vector(0);
  for (int i = 0; i < order; i++)
    {
      for (int j = 0; j < order; j++)
	{
	  coeffs_vector.push_back((*coeffs)[i][j]);
	}
    }

  double coeffs_array[coeffs_vector.size()];
  for (int i = 0; i < coeffs_vector.size(); i++)
    {
      coeffs_array[i] = coeffs_vector[i];
    }

  double known_array[order];
  for (int i = 0; i < order; i++)
    {
      known_array[i] = (*known)[i];
    }

  gsl_matrix_view m = gsl_matrix_view_array(coeffs_array, order, order);

  gsl_vector_view b = gsl_vector_view_array(known_array, order);

  gsl_vector *x = gsl_vector_alloc(order);

  int s;
  
  gsl_permutation *p = gsl_permutation_alloc(order);

  //Perform the decompostion
  gsl_linalg_LU_decomp(&m.matrix, p, &s);

  //Solve
  gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);

  for (int i = 0; i < order; i++)
    {
      (*unknown)[i] = gsl_vector_get(x, i);
    }

  //Free up space
  gsl_permutation_free(p);
  gsl_vector_free(x);

}
