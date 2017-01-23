#ifndef STDLP_H
#define STDLP_H

#include <iostream>
#include <limits>
#include <string>
#include <cmath>

using namespace std;

namespace simplex {

  template<typename Real, typename Integer, bool RowMajor=true>
  class StdLP {
  public:
    StdLP() : tol(1e-10) {
      holdmem  = true;
      iwork    = NULL;
      dwork    = NULL;
      phaseOne = true;
      degenerate = false;
    }

    ~StdLP() {
      if (holdmem) {
	if (iwork!=NULL) delete iwork;
	if (dwork!=NULL) delete dwork;
      }
    }

    void set(const Integer& sz, const Integer& ncon, const Real* cc,
	     const Real* aa, const Real* bb) {
      phaseOne   = true;
      degenerate = false;
      // allocate memory if necessary
      if (size!=sz || m!=ncon) {
	size = sz;
	n    = sz + ncon;
	m    = ncon;
      
	if (iwork!=NULL) delete iwork;
	if (dwork!=NULL) delete dwork;
	iwork = new Integer[n+n+m];
	dwork = new Real[n*(m+2)+m];

	pri    = iwork;
	pci    = pri+n;
	pstate = pci+m;
	pA     = dwork;
	pb     = pA+n*m;
	pc     = pb+m;
	pw     = pc+n;
      }
    
      // assign A
      int pAi = 0;
      int aai = 0;
      if (RowMajor) {
	// copy aa
	for (int i=0; i<m; i++) {
	  for (int j=0;j<sz;j++) {
	    if (bb[i]<0)
	      pA[pAi++] = -aa[aai++];
	    else
	      pA[pAi++] =  aa[aai++];
	  }
	  pAi += m;
	}
	// assign the art variables
	pAi = sz;
	for (int i=0; i<m; i++) {
	  for (int j=0; j<m; j++) {
	    if (i==j)
	      pA[pAi++] = 1;
	    else
	      pA[pAi++] = 0;
	  }
	  pAi += sz;
	}
      } else {
	// copy aa
	int idx = 0;
	for (int i=0; i<sz; i++) {
	  for (int j=0; j<m; j++) {
	    if (bb[m]<0)
	      pA[idx] = -aa[idx];
	    else
	      pA[idx] =  aa[idx];
	    idx++;
	  }
	}
	// assign the art variables
	pAi = sz*m;
	for (int j=0; j<m; j++) {
	  for (int i=0; i<m; i++) {
	    if (i==j)
	      pA[pAi++] = 1;
	    else
	      pA[pAi++] = 0;
	  }
	}
      }
      // assign b
      for (int i=0; i<n; i++) {
	pb[i] = abs(bb[i]);
      }
      // assign c
      for (int i=0; i<sz; i++) {
	pc[i] = cc[i];
      }
      for (int i=sz; i<n; i++) {
	pc[i] = 0;
      }
      // assign w
      for (int i=0; i<sz; i++) {
	w(i) = 0;
	for (int j=0; j<m; j++) {
	  w(i) -= A(j,i);
	}
      }  
      for (int i=sz; i<n; i++) {
	pw[i] = 0;
      }
      // assign wobj
      wobj = 0.0;
      for (int i=0; i<m; i++) {
	wobj += b(i);
      }
      // assign pstate
      for (int i=0; i<sz; i++) {
	pstate[i] = 0;
      }
      for (int i=sz; i<n; i++) {
	pstate[i] = 1;
      }
      // assign pci
      for (int i=0; i<ncon; i++) {
	pci[i] = sz+i;
      }
      // assign pri
      for (int i=sz; i<n; i++) {
	pri[i] = i-sz;
      }

    }
  
  private:
    // working arrays
    Integer  *iwork; // = pri(n) + pci(m) + pstate(n)
    Real     *dwork; // = pA(n*m)  + pb(m)  + pc(n) + pw(n)

    bool     holdmem;
  
    const Real  tol;
  
    // problem size
    Integer  size;

    // num of cols of A
    Integer  n;
    // num of rows of A
    Integer  m;

    // ob pointer to A
    Real     *pA;
    // ob pointer to b
    Real     *pb;
    // ob pointer to c
    Real     *pc;
    // ob pointer to w
    Real     *pw;

    // objs
    Real     wobj, cobj;

    // solution state
    Integer  solstate;

    // ob pointer to index of the nonzero rows
    // of each column (only meaningful when it
    // is active). len=n.
    Integer  *pri;
    // ob pointer to index of the basic cols
    // of each row
    Integer  *pci;

    // States of all the variables, including both admissible
    // and artificial variables. Artificial variables are
    // at the tail (len=n)
    Integer  *pstate;

    // phase
    bool     phaseOne;

    // num of pivot
    Integer  numPivot;

    // degeneracy
    bool degenerate;
  
    inline Real&  A(const Integer& ri, const Integer& ci) const {
      if (RowMajor)
	return pA[ri*n+ci];
      else
	return pA[ci*m+ri];
    };

    inline Real&  A(const Integer& li) const {
      return pA[li];
    };

  
    inline Real&  c(const Integer& i) const {
      return pc[i];
    };

    inline Real&  b(const Integer& i) const {
      return pb[i];
    };

    inline Real&  w(const Integer& i) const {
      return pw[i];
    };

    void dropCol(const Integer& ci) {
      pstate[ci] = 2;
    }

    // whether v < 0
    inline bool lt0(const Real& v) const {
      if (v < -tol) return true;
      return false;
    }

    // whether v == 0
    inline bool eq0(const Real& v) const {
      if (v > -tol and v < tol) return true;
      return false;
    }

    // whether v > 0
    inline bool gt0(const Real& v) const {
      if (v > tol) return true;
      return false;
    }
  
    void pivot(const Integer& ri, const Integer& ci) {
      Real ratio;
      // normalise A(ri,ci) and b(ri)
      ratio = 1.0/A(ri,ci);
      for (int i=0; i<n; i++) {
	A(ri,i) *= ratio;
      }
      b(ri) *= ratio;
      // process c, cobj
      ratio = c(ci);
      for (int i=0; i<n; i++) {
	c(i) -= A(ri,i)*ratio;
      }
      cobj += b(ri)*ratio;
      // process A, b
      for (int i=0; i<ri; i++) {
	ratio = A(i,ci);
	for (int j=0; j<n; j++)
	  A(i,j) -= A(ri,j)*ratio;
	// b
	b(i) -= b(ri)*ratio;
      }
      for (int i=ri+1; i<m; i++) {
	ratio = A(i,ci);
	for (int j=0; j<n; j++)
	  A(i,j) -= A(ri,j)*ratio;
	// b
	b(i) -= b(ri)*ratio;
      }
      // process w and wobj
      if (phaseOne) {
	ratio = w(ci);
	for (int i=0; i<n; i++) {
	  w(i) -= A(ri,i)*ratio;
	}
	wobj += b(ri)*ratio;
      }
    
      // update state
      pstate[ci] = 1;
      pstate[pci[ri]] = 0;
      pci[ri] = ci;
      pri[ci] = ri;
      numPivot++;
    }
  
    inline Integer getPivotCol() const {
      Real     mincoef = std::numeric_limits<Real>::max();
      Integer  minidx;

      if (phaseOne) {
	// phase I, optimize w
	for (int i=0; i<n; i++) {
	  if (pw[i]<mincoef) {
	    mincoef = pw[i];
	    minidx = i;
	  }
	}
	if (lt0(mincoef)) {
	  return minidx;
	} else {
	  return n; // means "the obj is optimal"
	}
      } else {
	// phase II, optimize c
	for (int i=0; i<n; i++) {
	  if (pstate[i]==2) continue;
	
	  if (c(i)<mincoef) {
	    mincoef = c(i);
	    minidx = i;
	  }
	}
	if (lt0(mincoef)) {
	  return minidx;
	} else {
	  return n; // means "the obj is optimal"
	}
      }
    }

    inline Integer getPivotRow(const Integer& ci) const {
      Real     minratio = std::numeric_limits<Real>::max();
      Integer  minidx   = -1;
      Real     ratio;

      for (int i=0; i<m; i++) {
	if (lt0(A(i,ci))) continue;

	ratio = b(i)/A(i,ci);
	if (ratio<minratio) {
	  minratio = ratio;
	  minidx   = i;
	}
      }
      if (minratio>=0)
	return minidx;
      else
	return m; // unbounded
    }
  
  public:
    void run(Real* x, Real& obj) {
      int ri, ci;
      phaseOne = true;

      while (true) {
	ci = getPivotCol();
	if (ci >= n) {
	  // optimal
	  if (phaseOne) {
	    // First, check feasiblity
	    if (gt0(wobj)) {
	      // infeasible
	      solstate = 1;
	      break;
	    }
	    // Second, detect degeneracy
	    int i = size;
	    for (; i<n; i++) {
	      if (pstate[i]==1) {
		degenerate = true;
		break;
	      }
	    }
	    if (degenerate) {
	      // degenerate
	      // 1. drop nonbasic artificial variables
	      for (int j=size; j<i; j++) {
		dropCol(j);
	      }
	      for (int j=i+1; j<n; j++) {
		if (pstate[i]==0)
		  dropCol(j);
	      }
	      // 2. drop non-art variables with d>0
	      for (int j=0; j<size; j++) {
		if (gt0(pw[j]))
		  dropCol(j);
	      }
	    } else {
	      // non degenerate, just remove art vars
	      for (int j=size; j<n; j++) {
		pstate[j] = 2;
	      }
	    }
	    // switch to phase II
	    phaseOne = false;
	  } else {
	    // phase II
	    ri = getPivotRow(ci);
	    if (ri>=m) {
	      // unbounded
	      solstate = 2;
	    } else {
	      // got solution
	      solstate = 0;
	      // return x
	      for (int i=0; i<size; i++) {
		if (pstate[i]==0) {
		  x[i] = 0;
		} else {
		  x[i] = b(pri[i]);
		}
	      }
	      // return obj
	      obj = cobj;
	    }
	    break;
	  }
	} else {
	  // not optimal, do the pivot
	  ri = getPivotRow(ci);
	  if (ri>=m) {
	    // unbounded
	    // (this cannot happen in phase I)
	    solstate = 2;
	    break;
	  }
	  pivot(ri, ci);
	}
      }
    }

    void dump(const string & header="") {
      if (header.length()==0)
	cout << "------------------------------------------" << endl;
      else
	cout << header << endl;
      // output w
      for (int i=0; i<n; i++) {
	cout << setw(5) << w(i) << '\t';
      }
      cout << setw(5) << wobj << endl;

      // output c
      for (int i=0; i<n; i++) {
	cout << setw(5) << c(i) << '\t';
      }
      cout << setw(5) << cobj << endl;

      // output A and b
      for (int i=0; i<m; i++) {
	for (int j=0; j<n; j++) {
	  cout << setw(5) << A(i,j) << '\t';
	}
	cout << setw(5) << b(i) << endl;
      }
    
      // output pstate
      for (int i=0; i<n; i++) {
	cout << pstate[i] << '\t';
      }    
      cout << endl;
    }
  };
}
#endif /* STDLP_H */
