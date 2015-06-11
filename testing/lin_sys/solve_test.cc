//Program file to test the solution of linear systems

#include <vector>
#include <iostream>

#include "../../solve.h"

using std::vector;
using std::cout;
using std::endl;

int main()
{
  vector<double> known(5);
  vector<vector< double> > matrix(5);
  for (int i = 0; i < matrix.size(); i++)
    {
      matrix[i].resize(matrix.size());
    }
  vector<double> unknown(5);

  known[0] = 0.0;
  known[1] = 0.0;
  known[2] = 0.0;
  known[3] = 0.0;
  known[4] = 3.0;

  matrix[0][0] = -0.011073;
  matrix[0][1] = -0.0648877;
  matrix[0][2] = 3.2e-17;
  matrix[0][3] = 0.0648877;
  matrix[0][4] = 0.0;

  matrix[1][0] = -0.135299;
  matrix[1][1] = 1.06394;
  matrix[1][2] = 0.438211;
  matrix[1][3] = 0.196471;
  matrix[1][4] = -1.0;

  matrix[2][0] = -1.4e-17;
  matrix[2][1] = 0.213235;
  matrix[2][2] = 0.474671;
  matrix[2][3] = 0.213235;
  matrix[2][4] = -1.0;

  matrix[3][0] = 0.135299;
  matrix[3][1] = 0.196471;
  matrix[3][2] = 0.438211;
  matrix[3][3] = 1.06394;
  matrix[3][4] = -1.0;

  matrix[4][0] = 0.0;
  matrix[4][1] = 1.0;
  matrix[4][2] = 1.0;
  matrix[4][3] = 1.0;
  matrix[4][4] = 0.0;

  Solve(5, &matrix, &known, &unknown);

  for (int i = 0; i < unknown.size(); i++)
    {
      cout << i << " " << unknown[i] << endl;
    }

  matrix.erase(matrix.begin() + 4);
  for (int i = 0; i < matrix[i].size(); i++)
    {

    }
  return 0;
}
