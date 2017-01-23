#ifndef BOUNDLP_H
#define BOUNDLP_H

#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <cmath>

using namespace std;

namespace simplex {
  
  template<typename Real, typename Integer>
  class BoundLP {
  public:
    // use c++11 syntax
    enum VariableState : Integer {
      VS_FREE, VS_LOWER, VS_UPPER,
      VS_DROP_FROM_LOWER, VS_DROP_FROM_UPPER
    };
    
    
    BoundLP() : tol(1e-10) {
      holdmem  = true;
      iwork    = NULL;
      dwork    = NULL;
      phaseOne = true;
      degenerate = false;
    }

    ~BoundLP() {
      if (holdmem) {
	if (iwork!=NULL) delete iwork;
	if (dwork!=NULL) delete dwork;
      }
    }

    void set(const Integer& sz, const Integer& ncon, const Real* cc,
	     const Real* aa, const Real* bb, Real*llbb, Real*uubb) {
      phaseOne   = true;
      degenerate = false;
      // allocate memory if necessary
      if (size!=sz || m!=ncon) {
	size = sz;
	n    = sz + ncon;
	m    = ncon;
      
	if (iwork!=NULL) delete iwork;
	if (dwork!=NULL) delete dwork;
	iwork = new Integer[n*3+m];
	dwork = new Real[n*(m+4)+m];

	pri    = iwork;
	pci    = pri+n;
	pstate = pci+m;
	pkx    = pci+n;

	pA     = dwork;
	pb     = pA+n*m;
	pc     = pb+m;
	pw     = pc+n;
	plb    = pw+n;
	pub    = plb+n;
      }

      int pAi = 0;
      int aai = 0;

      // initialise kx
      for (int i=0; i<n; i++) {
	kx(i) = i;
      }
      
      // assign lb
      if (llbb==NULL)
	for (int i=0; i<sz; i++)
	  plb[i] = -numeric_limits<double>::max();
      else
	for (int i=0; i<sz; i++)
	  plb[i] = llbb[i];

      for (int i=sz; i<n; i++)
	plb[i] = 0;
      
      // assign ub
      if (uubb==NULL)
	for (int i=0; i<sz; i++)
	  pub[i] = numeric_limits<double>::max();
      else
	for (int i=0; i<sz; i++)
	  pub[i] = uubb[i];

      for (int i=sz; i<n; i++)
	pub[i] = numeric_limits<double>::max();

      // assign A and b
      pAi = 0;
      aai = 0;
      // copy aa
      for (int i=0; i<m*sz; i++) {
	    pA[i] = aa[i];
      }
      // assign the artificial variables
      pAi = sz*m;
      for (int i=0; i<m; i++) {
	for (int j=0; j<m; j++) {
	  if (i==j)
	    pA[pAi++] = 1;
	  else
	    pA[pAi++] = 0;
	}
      }

      // assign c
      for (int i=0; i<sz; i++) {
	pc[i] = cc[i];
      }
      for (int i=sz; i<n; i++) {
	pc[i] = 0;
      }

      for (int i=sz; i<n; i++) {
	pw[i] = 0;
      }

      // assign pstate
      for (int i=0; i<sz; i++) {
	  pstate[i] = VS_LOWER;
      }
      for (int i=sz; i<n; i++) {
	pstate[i] = VS_FREE;
      }

      // test if we need to negate the sign
      // and assign b
      for (int i=0; i<m; i++) {
	pb[i] = bb[i];
      }
      pAi = 0;

      for (int i=0; i<sz; i++)
	for (int j=0; j<m; j++) {
	  b(j) -= lb(i)*A(pAi);
	  pAi++;
	}

      for (int i=0; i<m; i++) {
	if (b(i)<0) {
	  b(i) = -bb[i];
	  pAi = i;
	  for (int i=0; i<sz;i++) {
	    A(pAi) = -A(pAi);
	    pAi += m;
	  }
	}
	else {
	  b(i) = bb[i];
	}
      }
      
      // assign w
      for (int i=0; i<sz; i++) {
	w(i) = 0;
	for (int j=0; j<m; j++) {
	  w(i) -= A(j,i);
	}
      }

      // assign bw
      bw = 0.0;
      for (int i=0; i<m; i++) {
	bw -= b(i);
      }

      // assign bc
      bc = 0.0;
      
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
    Integer  *iwork; // = pri(n) + pci(m) + pstate(n) + pkx(n)
    Real     *dwork; // = pA(n*m)  + pb(m)  + pc(n) + pw(n) +  pub(n) + plb(n)

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
    // ob pointer to lb
    Real     *plb;
    // ob pointer to ub
    Real     *pub;
    

    // objs
    Real     bw, bc;

    // solution state
    Integer  solstate;

    // index of columns
    Integer  *pkx;

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
  
    void pivot(const Integer& ri, const Integer& ci) {
      //      cout << ri << "," << ci << endl;
      Real ratio;
      // store A(ri,ci) before change the value
      // Arc is used in updating the state
      Real Arc = A(ri,ci);
      // normalise A(ri,ci) and b(ri)
      ratio = 1.0/A(ri,ci);
      for (int i=0; i<n; i++) {
	A(ri,i) *= ratio;
      }
      b(ri) *= ratio;
      // process c, bc
      ratio = c(ci);
      for (int i=0; i<n; i++) {
	c(i) -= A(ri,i)*ratio;
      }
      bc -= b(ri)*ratio;
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
      // process w and bw
      if (phaseOne) {
	ratio = w(ci);
	for (int i=0; i<n; i++) {
	  w(i) -= A(ri,i)*ratio;
	}
	bw -= b(ri)*ratio;
      }

      // update state
      if (pstate[ci]==VS_LOWER) {
	if (gt0(Arc)) {
	  pstate[pci[ri]] = VS_LOWER;
	} else {
	  pstate[pci[ri]] = VS_UPPER;
	}
      } else {
	if (gt0(Arc)) {
	  pstate[pci[ri]] = VS_UPPER;
	} else {
	  pstate[pci[ri]] = VS_LOWER;
	}
      }

      pstate[ci] = VS_FREE;
      
      pci[ri] = ci;
      pri[ci] = ri;
      numPivot++;

      //                  dump(string("Pivot (") + to_string(ri) + string(",") \
      //            	   + to_string(ci) + string(") -------------------------"));
    }
  
    inline Integer getPivotCol() const {
      Real     maxgrad = std::numeric_limits<Real>::min();
      Integer  maxidx;
      Real     *workvec;  // the object vector work work with

      // decide which phase
      if (phaseOne) {
	workvec = pw;
      } else {
	workvec = pc;
      }

      // do the work
      for (int i=0; i<n; i++) {
	//	if (pstate[i]==VS_DROP_FROM_LOWER or	\
	//	    pstate[i]==VS_DROP_FROM_UPPER)
	//	  continue;
	  
	if (pstate[i]==VS_LOWER and lt0(workvec[i])) {
	  // can increase x(i) to decrease obj
	  //	  if (-workvec[i] > maxgrad) {
	  //	    maxgrad = -workvec[i];
	  //	    maxidx  = i;
	  //      }
	  return i;
	} else if (pstate[i]==VS_UPPER and gt0(workvec[i])) {
	  // can decrease x(i) to decrease obj
	  //	  if (workvec[i] > maxgrad) {
	  //	    maxgrad = workvec[i];
	  //	    maxidx  = i;
	  //      }
	  return i;
	}
      }
      //      if (gt0(maxgrad)) {
      //	return maxidx;
      //      } else {
      //	return n; // means "the obj is optimal"
      //      }
      return n; // means "the obj is optimal"
  }

    inline Integer getPivotRow(const Integer& ci) const {
      Real     minabschg = std::numeric_limits<Real>::max();
      Integer  minidx    = std::numeric_limits<Integer>::max();
      Real     abschg    =  0; // absolute change
      Real     value     =  0; // value of basic variable
      
      Integer  pAi = m*ci;

      for (int i=0; i<m; i++) {
	if (eq0(A(pAi))) {
	  // abschg is infty
	  pAi++;
	  continue;
	} else {
	  // compute the value of current basic variable
	  value = b(i);
	  for (int j=0; j<size; j++) {
	    if (pstate[j]==VS_LOWER)
	      value -= A(i,j)*lb(j);
	    if (pstate[j]==VS_UPPER)
	      value -= A(i,j)*ub(j);
	  }

	  //	  cout << "value : " << lb(pci[i]) << " < " << value << " < " << ub(pci[i]) << endl;
	  
	  if (pstate[ci]==VS_LOWER) {
	    // x[ci] is going to increase
	    if (gt0(A(pAi))) {
	      // A(i,ci)*x[ci] is going to increase
	      // so the corresponding basic solution
	      // is going to decrease
	      abschg = (value-lb(pci[i]))/A(pAi);
	    } else {
	      // A(i,ci)*x[ci] is going to decrease
	      // so the corresponding basic solution
	      // is going to increase
	      abschg = (value-ub(pci[i]))/A(pAi);
	    }
	  } else {
	    // x[ci] is going to decrease
	    if (gt0(A(pAi))) {
	      // A(i,ci)*x[ci] is going to decrease
	      // so the corresponding basic solution
	      // is going to increase
	      abschg = (ub(pci[i])-value)/A(pAi);
	      //	      cout << abschg << "aaaaaaaaaaaaa" <<ub(pci[i])<< "aaaaaa" << value << endl;
	    } else {
	      // A(i,ci)*x[ci] is going to increase
	      // so the corresponding basic solution
	      // is going to decrease
	      abschg = (lb(pci[i])-value)/A(pAi);
	    }
	  }
	}
	if (abschg<minabschg) {
	  minabschg = abschg;
	  minidx    = i;
	}
	
	pAi++;
      }

      if (minabschg>=(ub(ci)-lb(ci))) {
	// donot do pivot, just change the current state
	//	cout << ub(ci)-lb(ci) << "," << w(ci) << "," << getObj() << endl;
	return -1;
      }
      //      cout << ci << "," << minabschg << "," << minidx << "," << w(ci) << "," << getObj() << endl;
      if (ge0(minabschg))
	return minidx;
      else
	return m; // unbounded
    }

    inline void swapCol(const Integer &ci1, const Integer &ci2) {
      int idx1 = ci1*m;
      int idx2 = ci2*m;
      Real rtmp;
      for (int i=0; i<m; i++) {
	rtmp   = A(idx1);
	A(idx1) = A(idx2);
	A(idx2) = rtmp;
	idx1++;
	idx2++;
      }
      Integer itmp;
      itmp    = kx(ci1);
      kx(ci1) = kx(ci2);
      kx(ci2) = itmp;
    }
    
  public:
    void run(Real* x, Real& obj, Integer& inform) {

      //      dump("[Init]-------------------------");
      
      int ri, ci;
      phaseOne = true;

      while (true) {
	ci = getPivotCol();
	if (ci >= n) {
	  // optimal
	  if (phaseOne) {
	    // First, check feasiblity
	    Real truebw = -bw;
	    for (int i=0;i<size;i++) {
	      if (pstate[i]==VS_LOWER)
		truebw += lb(i)*w(i);
	      if (pstate[i]==VS_UPPER)
		truebw += ub(i)*w(i);
	    }
	    if (gt0(truebw)) {
	      // infeasible
	      solstate = 1;
	      break;
	    } else {
	      //	      	      	      dump(string("Feasible : ") + to_string(truebw));
	    }
	    // Second, detect degeneracy
	    int i = size;
	    degenerate = false;
	    for (; i<n; i++) {
	      if (pstate[i]==VS_FREE) {
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
		if (pstate[i]==VS_LOWER)
		  dropCol(j);
	      }
	      // 2. drop non-art variables with
	      //    2.1 d>0 and state==VS_LOWER
	      //    2.2 d<0 and state==VS_UPPER
	      for (int j=0; j<size; j++) {
		if ( (gt0(pw[j]) and pstate[j]==VS_LOWER) or \
		     (lt0(pw[j]) and pstate[j]==VS_UPPER) )
		  dropCol(j);
	      }
	    } else {
	      // non degenerate, just remove art vars
	      for (int j=size; j<n; j++) {
		pstate[j] = VS_DROP_FROM_LOWER; // upper or lower not matter
	      }
	    }
	    // switch to phase II
	    phaseOne = false;
	  } else {
	    // got solution
	    solstate = 0;
	    // return x
	    for (int i=0; i<size; i++) {
	      if (pstate[i]==VS_LOWER or \
		  pstate[i]==VS_DROP_FROM_LOWER) {
		x[i] = lb(i);
	      } else if (pstate[i]==VS_UPPER or \
			 pstate[i]==VS_DROP_FROM_UPPER) {
		x[i] = ub(i);
	      } else {
		// basic variables
		x[i] = b(pri[i]);
		int pAi = pri[i];
		for (int j = 0; j<size; j++) {
		  if (pstate[j]==VS_LOWER or \
		      pstate[j]==VS_DROP_FROM_LOWER) 
		    x[i] -= A(pAi)*lb(j);
		  else if (pstate[j]==VS_UPPER or \
			   pstate[j]==VS_DROP_FROM_UPPER)
		    x[i] -= A(pAi)*ub(j);
		  pAi += m;
		}
	      }
	    }
	    // return obj
	    obj = -bc;
	    for (int i=0; i<size; i++) {
	      obj += x[i]*c(i);
	    }
	    break;
	  }
	} else {
	  // not optimal, do the pivot
	  ri = getPivotRow(ci);
	  if (ri<0) {
	      // no need to do pivot
	      if (pstate[ci]==VS_UPPER) {
		pstate[ci] = VS_LOWER;
	      } else {
		pstate[ci] = VS_UPPER;
	      }
	      //	      dump(string("Pivot (-,") +
	      //		   to_string(ci) + string(") -------------------------"));
	      continue;
	    }

	  if (ri>=m) {
	    // unbounded
	    // (this cannot happen in phase I)
	    solstate = 2;
	    break;
	  }
	  pivot(ri, ci);
	}
      }
      inform = solstate;
    }

    void dump(const string & inhead="") const {
      string header = inhead;
      if (header.length()==0)
	header = "------------------------------------------";
	
      if (phaseOne)
	header = header + string(" Phase I");
      else
	header = header + string(" Phase II");

      cout << header << endl;
      // output w
      cout << "w:\t";
      for (int i=0; i<n; i++) {
	cout << fixed << setw(5) << setprecision(3) << w(i) << '\t';
      }
      cout << fixed << setw(5) << setprecision(3) << bw << endl;

      // output c
      cout << "c:\t";
      for (int i=0; i<n; i++) {
	cout << fixed << setw(5) << setprecision(3) << c(i) << '\t';
      }
      cout << fixed << setw(5) << setprecision(3) << bc << endl;

      // output A and b
      for (int i=0; i<m; i++) {
	cout << setprecision(5) << "A:\t";	
	for (int j=0; j<n; j++) {
	  cout << fixed << setw(5) << setprecision(3) << A(i,j) << '\t';
	}
	cout << fixed << setw(5) << setprecision(3) << b(i) << endl;
      }
    
      // output pstate
      cout << "s:\t";
      for (int i=0; i<n; i++) {
	switch (pstate[i]) {
	case VS_FREE:
	  cout << fixed << setw(5) << "   F\t";
	  break;
	case VS_LOWER:
	  cout << fixed << setw(5) << "   L\t";
	  break;
	case VS_UPPER:
	  cout << fixed << setw(5) << "   U\t";
	  break;
	case VS_DROP_FROM_LOWER:
	  cout << fixed << setw(5) << "  Dl\t";
	  break;
	case VS_DROP_FROM_UPPER:
	  cout << fixed << setw(5) << "  Du\t";
	  break;
	};
	
      }    
      cout << endl;
    }

    // Helper functions
    inline Real&  A(const Integer& ri, const Integer& ci) const {
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

    inline Real&  lb(const Integer& i) const {
      return plb[i];
    };

    inline Real&  ub(const Integer& i) const {
      return pub[i];
    };

    inline Integer&  kx(const Integer& i) const {
      return pkx[i];
    };


    
    void dropCol(const Integer& ci) {
      if (pstate[ci] == VS_LOWER)
	pstate[ci] = VS_DROP_FROM_LOWER;
      else
	pstate[ci] = VS_DROP_FROM_UPPER;
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

    // whether v > 0
    inline bool ge0(const Real& v) const {
      if (v > -tol) return true;
      return false;
    }

    // compute obj
    Real getObj() const {
      Real *workobj;
      Real ret;
      if (phaseOne) {
	workobj =  pw;
	ret     = -bw;
      }
      else {
	workobj =  pc;
	ret     = -bc;
      }
      
      for (int i=0; i<size; i++) {
	if (pstate[i]==VS_LOWER or pstate[i]==VS_DROP_FROM_LOWER)
	  ret += workobj[i]*lb(i);
	else if (pstate[i]==VS_UPPER or pstate[i]==VS_DROP_FROM_UPPER)
	  ret += workobj[i]*ub(i);
      }

      return ret;
    }
    
  };
}
#endif /* BOUNDLP_H */
