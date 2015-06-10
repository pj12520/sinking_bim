#include <vector>
#include <iostream>

using namespace std;

typedef std::vector<int> V1D;
typedef std::vector<V1D> V2D;

void PrintVector2D(const V2D & v)
{
    for (int i=0; i<v.size(); ++i)
    {
        for (int j=0; j<v[i].size(); ++j)
        {
            cout << v[i][j] << " ";
        }
        cout << endl;
    }
}

void PrintVectortest(vector<vector<double> >* v)
{
  for (int i=0; i<(*v).size(); ++i)
    {
      for (int j=0; j<(*v)[i].size(); ++j)
        {
	  cout << (*v)[i][j] << " ";
        }
        cout << endl;
    }
}


int main()
{
    V2D v(3,V1D(3));  // 3 x 3
    vector<vector<double> > test;
    test.resize(v.size());
    for (int i = 0; i < test.size(); i++)
      {
	test[i].resize(v[i].size());
      }
    // give some values

    int ic = 0;

    for (int i=0; i<v.size(); ++i)
    {
        for (int j=0; j<v[i].size(); ++j)
        {
            v[i][j] = ++ic;
	    test[i][j] = 0.5 * i * j;
        }
    }

    // print out original matrix

    cout << "\n\nThe original matrix is \n\n";
    PrintVector2D(v);
    PrintVectortest(&test);
    // erase a row
    v.erase(v.begin()+1); // erase 2nd row
    test.erase(test.begin()+1);
    //print out

    cout << "\n\nThe matrix after erasing second row is \n\n";
    PrintVector2D(v);
    PrintVectortest(&test);

    return 0;
}
