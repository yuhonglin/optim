#include <iostream>
using namespace std;

#include "mat.h"

#include "BOUNDLP.hpp"

int main(int argc, char *argv[])
{
  simplex::BoundLP <double, int> sim;

  MATFile *pmat = matOpen("./data.mat", "r");
  if (pmat == NULL) {
    cout << "Error opening file data.mat" << endl;
    return(1);
  }
  
  mxArray *c   = matGetVariable(pmat, "c");
  mxArray *Aeq = matGetVariable(pmat, "Aeq");
  mxArray *beq = matGetVariable(pmat, "beq");
  mxArray *lb  = matGetVariable(pmat, "lb");
  mxArray *ub  = matGetVariable(pmat, "ub");  

  int size = mxGetM(c);
  int m    = mxGetM(Aeq);
  
  sim.set(size, m, mxGetPr(c), mxGetPr(Aeq), mxGetPr(beq), mxGetPr(lb), mxGetPr(ub));

  double *x = new double[size];
  double obj;
  int inform;
  sim.run(x,obj, inform);

  cout << "--------------------------" << endl;
  for (int i=0; i<size; i++)
    cout << x[i] << '\t';
  cout << endl;
  cout << "obj    = " << obj << endl;
  cout << "inform = " << inform << endl;

  delete[] x;
  return 0;
}
