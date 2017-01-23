#include <cmath>

#include "mex.h"
#include "BOUNDLP.hpp"

using namespace simplex;


// prhs = [c, Aeq, beq, lb, ub]
// plhs = [x, inform]
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nlhs<1)
    return;
  
  if (nrhs<5)
    mexErrMsgTxt("Too few inputs");

  int size = max(mxGetM(prhs[0]), mxGetN(prhs[0]));
  if (mxGetN(prhs[1])!=size) {
    mexErrMsgTxt("1st and 2nd inputs do not match");
  }
  int m = max(mxGetM(prhs[2]), mxGetN(prhs[2]));
  if (m!=mxGetM(prhs[1])) {
    mexErrMsgTxt("2nd and 3rd inputs do not match");
  }
  if (max(mxGetM(prhs[3]), mxGetN(prhs[3]))!=size) {
    mexErrMsgTxt("1st and 4th inputs do not match");
  }
  if (max(mxGetM(prhs[4]), mxGetN(prhs[4]))!=size) {
    mexErrMsgTxt("1st and 5th inputs do not match");
  }

  BoundLP<double, int> sim;
  sim.set(size, m, mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]),
	  mxGetPr(prhs[3]), mxGetPr(prhs[4]));

  plhs[0] = mxCreateDoubleMatrix(size,1,mxREAL);
  double obj;
  int inform;

  sim.run(mxGetPr(plhs[0]), obj, inform);

  // other output
  if (nlhs>1) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[1])) = obj;
  }
  if (nlhs>2) {
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(plhs[2])) = static_cast<double>(inform);
  }
}
