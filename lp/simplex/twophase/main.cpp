#include <iostream>
using namespace std;

#include "StdLP.hpp"

int main(int argc, char *argv[])
{
  simplex::StdLP<double, int> sim;

  double c[5]   = {1,1,1,0,0};
  double A[5*2] = {1,1,0,1,0,  0,-1,1,0,1};
  double b[2]   = {1,0};

  sim.set(5,2,c,A,b);

  double x[5];
  double obj;
  sim.run(x,obj);

  cout << "--------------------------" << endl;
  cout << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << '\t' << x[4] << endl;
  cout << obj << endl;
  
  return 0;
}
