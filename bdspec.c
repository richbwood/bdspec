// Copyright 2010 Richard Wood. All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
//    1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY Richard Wood ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL <COPYRIGHT HOLDER> OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// The views and conclusions contained in the software and documentation are
// those of the authors and should not be interpreted as representing official
// policies, either expressed or implied, of Richard Wood.

#include <stdio.h>
#include <stdlib.h>

// Required by bdspec:
#include <math.h>
#include <float.h>
#define SWAP(a, b) {double dum = (a); (a) = (b); (b) = dum;}
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Required for testing only:
#include <time.h>
#include <stdbool.h>
#include <meschach/matrix.h>
#include <meschach/matrix2.h>
#ifndef NAN
#define NAN 0.0/0.0
#endif


// Required for testing only:
// Dump matrix 'a' to stdout
void m_out(double *a, int n, int m) {
  const int d = 3; //digits of accuracy
  double(*aa)[m] = (double(*)[m])a;
  for(int i = 0; i < n; i++) {
    printf("[%2d]", i);
    for(int j = 0; j < m; j++) {
      if(FP_ZERO == fpclassify(aa[i][j]))
        printf(" %.*s", d+3, "                 ");
      else
        printf(" %*.*F", d+3, d, aa[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

// Required for testing only:
// Dump vector 'a' to stdout
void v_out(double a[], int len) {
  const int d = 3; //digits of accuracy
  for(int i = 0; i < len; i++) {
    printf("[%2d] %*.*f\n", i, d+3, d, a[i]);
  }
  printf("\n");
}

// Required for testing only:
// Dump permutation 'p' to stdout
void p_out(int p[], int len) {
  for(int i = 0; i < len; i++) {
    printf("[%2d] %2d\n", i, p[i]);
  }
  printf("\n");
}

// Required for testing only:
// Fill 'a' with zeros
void zero(double a[], int len) {
  for(int i = 0; i < len; i++) {
    a[i] = 0.0;
  }
}

// Required for testing only:
// Fill 'a' with random numbers uniformly distributed between -1 and 1
void  randlist(double a[], int len) {
  for(int i = 0; i < len; i++) {
    a[i] = 2.0*drand48()-1.0;
  }
}

// Required for testing only:
// Calculate the L_2 norm of a-b
double v_diff(double a[], double b[], int len) {
  double acc = 0.0;
  for(int i = 0; i < len; i++) {
    double diff = a[i]-b[i];
    acc += diff*diff;
  }
  return sqrt(acc);
}

// Calculates bA*x for (n*n) special-band-matrix bA and (n) vector x.
// -- the (n*n) matrix A is stored in band storage array bA of dimension
// (n*[2*lb+1+ub]). Elements are stored in columns lb to 2*lb+ub.  The jth row
// of bA is stored in the jth row of band A (bA) as follows: bA[i][2*lb+j-i] =
// A[i][j] for max(0, i-lb) <= j <= min(A->n-1, i+ub); Note that the first lb
// columns of bA are unused!
void bdspecLUmlt(double *bA, int n, int lb, int ub, double x[], double out[]) {
  int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
  double acc = 0.0;
  for(int i = n-1; i >= 0; i--) {
    {
      int k = i-2*lb;
      out[i] = 0.0;
      int j0 = MAX(lb, -k);
      int j1 = MIN(mm-1, n-k);
      for(int j = j0, l = j0+k; j < j1; j++, l++)
        out[i] += a[i][j]*x[l];
    }
    {
      int k = ub+i;
      if(k < n) {
        acc += x[k];
        out[i] += a[i][mm-1]*acc;
      }
    }
  }
}

// bdspecLUfactor -- Gaussian elimination with partial pivoting
// -- on entry, the (n*n) special-band-matrix A is stored in band storage array
// bA of dimension (n*[2*lb+1+ub]). Elements are stored in columns lb to
// 2*lb+ub.  The jth row of bA is stored in the jth row of A  as follows:
// bA[i][2*lb+j-i] = A[i][j] for max(0, i-lb) <= j <= min(A->n-1, i+ub); Note
// that the first lb columns of bA are unused!
// -- on exit: Upper triangular matrix U is stored with diagonal and lb+ub
// super-diagonals in columns lb to 2*lb+ub. The diagonals of lower triangular
// matrix L are ones, other nonzero elements of lower triangular matrix L is
// stored in columns 0 to lb-1.
// -- Array indx stores the row re-ordering that takes place
void bdspecLUfactor(double *bA, int n, int lb, int ub, int indx[]) {
  int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
  // Ensure repeated columns are repeated
  for(int k = 0; k < n; k++)
    for(int j = mm; j < mm; j++)
      a[k][j] = a[k][mm-1];
  // Rearrange the storage:
  for(int i = 0, l = lb; i < lb; i++, l--) {
    for(int j = lb+lb-i; j < mm; j++)
      a[i][j-l] = a[i][j];
    for(int j = mm-l; j < mm; j++)
      a[i][j] = a[i][mm-l-1];
  }
  // For each row
  for(int k = 0, l = lb; k < n; k++) {
    if(l < n) l++;
    // Find the pivot element:
    double maxa0 = fabs(a[k][lb]);
    int i = k;
    for(int j = k+1; j < l; j++) {
      double fabsa = fabs(a[j][lb]);
      if(fabsa > maxa0) {
        maxa0 = fabsa;
        i = j;
      }
    }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(maxa0 == 0.0) // If singular, give up and go to the next row
      continue;
    // Interchange rows
    if(i != k) {
      for(int j = lb; j < mm; j++)
        SWAP(a[k][j], a[i][j]);
    }
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      double fac = a[i][lb]/a[k][lb];
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++)
        a[i][j-1] = a[i][j]-fac*a[k][j];
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// Same as bdspecLUfactor, but pivots on rows that are normalized such that the
// greatest element in each row is +1 or -1. Note that this only affects which
// elements are pivoted: the resulting LU factorization is not normalized.
void bdspecLUfactorscale(double *bA, int n, int lb, int ub, int indx[]) {
  int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
  double scale[n];
  // Ensure repeated columns are repeated
  for(int k = 0; k < n; k++)
    for(int j = mm; j < mm; j++)
      a[k][j] = a[k][mm-1];
  // Calculate the scale of each row
  for(int i = 0; i < n; i++) {
    double max = 0.0;
    int k = i-2*lb;
    int j0 = MAX(lb, -k);
    int j1 = MIN(mm, n-k);
    for(int j = j0; j < j1; j++)
      max = fmax(fabs(a[i][j]), max);
    scale[i] = max;
  }
  // Rearrange the storage:
  for(int i = 0, l = lb; i < lb; i++, l--) {
    for(int j = lb+lb-i; j < mm; j++)
      a[i][j-l] = a[i][j];
    for(int j = mm-l; j < mm; j++)
      a[i][j] = a[i][mm-l-1];
  }
  // For each row
  for(int k = 0, l = lb; k < n; k++) {
    if(l < n) l++;
    // Find the pivot element:
    double maxa0 = 0.0;
    int i = -1;
    for(int j = k; j < l; j++) {
      if(scale[j] >= DBL_EPSILON*fabs(a[j][lb])) {
        double temp = fabs(a[j][lb])/scale[j];
        if(temp > maxa0) {
          maxa0 = temp;
          i = j;
        }
      }
    }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(i == -1) { // If singular continue to next row
      a[k][lb] = 0.0;
      continue;
    }
    // Interchange rows
    if(i != k) {
      SWAP(scale[k], scale[i]);
      for(int j = lb; j < mm; j++)
        SWAP(a[k][j], a[i][j]);
    }
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      double fac = a[i][lb]/a[k][lb];
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++)
        a[i][j-1] = a[i][j]-fac*a[k][j];
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// Same as bdspecLUfactorscale, but includes a big found in meschach 1.2.
// Having identical bugs allows testing of the LU factorization
void bdspecLUfactormeschscale(double *bA, int n, int lb, int ub, int indx[]) {
  int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
  double scale[n];
  // Ensure repeated columns are repeated
  for(int k = 0; k < n; k++)
    for(int j = mm; j < mm; j++)
      a[k][j] = a[k][mm-1];
  // Calculate the scale of each row
  for(int i = 0; i < n; i++) {
    double max = 0.0;
    int k = i-2*lb;
    int j0 = MAX(lb, -k);
    int j1 = MIN(mm, n-k);
    for(int j = j0; j < j1; j++)
      max = fmax(fabs(a[i][j]), max);
    scale[i] = max;
  }
  // Rearrange the storage:
  for(int i = 0, l = lb; i < lb; i++, l--) {
    for(int j = lb+lb-i; j < mm; j++)
      a[i][j-l] = a[i][j];
    for(int j = mm-l; j < mm; j++)
      a[i][j] = a[i][mm-l-1];
  }
  // For each row
  for(int k = 0, l = lb; k < n; k++) {
    if(l < n) l++;
    // Find the pivot element:
    double maxa0 = 0.0;
    int i = -1;
    for(int j = k; j < l; j++) {
      if(scale[j] >= DBL_EPSILON*fabs(a[j][lb])) {
        double temp = fabs(a[j][lb])/scale[j];
        if(temp > maxa0) {
          maxa0 = temp;
          i = j;
        }
      }
    }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(i == -1) { // If singular continue to next row
      a[k][lb] = 0.0;
      continue;
    }
    // Interchange rows
    if(i != k) {
      // Meschach doesn't have this:
      //SWAP(scale[k], scale[i]);
      for(int j = lb; j < mm; j++)
        SWAP(a[k][j], a[i][j]);
    }
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      double fac = a[i][lb]/a[k][lb];
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++)
        a[i][j-1] = a[i][j]-fac*a[k][j];
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// bdLUsolve -- given an LU factorisation in bA, solve bA*x=b for x
// -- Solution is stored in-place in array b.
// -- L and U are stored together in array bLU.  Upper triangular matrix U is
// stored with diagonal and lb+ub super-diagonals in columns lb to 2*lb+ub. The
// diagonals of lower triangular matrix L are ones, other nonzero elements of
// lower triangular matrix L is stored in columns 0 to lb-1.
// -- Array indx stores the row reordering required
void bdspecLUsolve(double *bLU, int n, int lb, int ub, int indx[], double b[]) {
  int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bLU;
  // Forward substitution:
  for(int k = 0, l = lb+1; k < n; k++, l++) {
    l = MIN(l, n);
    // b must be reordered during forward substitution, not before
    // quicker to swap without testing
    SWAP(b[k], b[indx[k]]);
    for(int i = k+1, j = 0; i < l; i++, j++)
      b[i] -= a[k][j]*b[k];
  }
  // Backward substitution:
  double acc2 = 0.0;
  for(int i = n-1, l = 1+lb; i >= 0; i--, l++) {
    l = MIN(l, mm-1);
    double acc = b[i];
    for(int k = lb+1, j = i+1; k < l; k++, j++) {
      acc -= a[i][k]*b[j];
    }
    {
      int k = MAX(ub+lb, 1)+i;
      if(k < n) {
        acc2 += b[k];
        acc -= a[i][mm-1]*acc2;
      }
    }
    b[i] = acc/a[i][lb];
  }
}

void test(int n, int lb, int ub, bool dumpfull) {
  // Generate a special-band-matrix mfull With lb lower diagonals, ub upper
  // diagonals and the elements of the highest upper diagonal extended across
  // each row
  MAT *mfull = m_get(n, n);
  //m_rand(mfull);
  randlist(mfull->base, (mfull->n)*(mfull->n));
  double **me = mfull->me;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < i - lb; j++)
      me[i][j] = 0.0;
    for(int j = i+ub+1; j < n; j++)
      me[i][j] = me[i][i+ub];
  }

  // Copy matrix mfull to a compactly stored version
  // mcmpct
  // First lb columns padding for later use
  // Next lb columns for lower diagonals
  // Next column for diagonal
  // Next ub columns for upper diagonals
  // as the highest upper diagonal of the same row
  int mm = 2*lb+ub+1;
  double mcmpct[n][mm];
  zero(&mcmpct[0][0], n*mm);
  for(int i = 0; i < n; i++) {
    for(int j = MAX(i-lb, 0); j < MIN(i+ub+1, n); j++) {
      mcmpct[i][j-i+lb+lb] = me[i][j];
    }
  }
  // Replace unused values with NAN to be sure they aren't used
  for(int k = 0; k < n; k++)
    for(int i = 0; i < lb; i++)
      mcmpct[k][i] = NAN;
  for(int k = 0; k < lb; k++) {
    for(int i = 0; i < lb-k; i++) {
      mcmpct[k][i+lb] = NAN;
    }
  }
  for(int k=n-1; k >= n-ub; k--) {
    for(int i = n-1-k+1+lb; i < mm; i++)
      mcmpct[k][i+lb] = NAN;
  }

  // Generate start vector x1 for test
  VEC *x1 = v_get(n);
  randlist(x1->ve, n);

  // Calculate mfull*x1 = dfull
  VEC *dfull = v_get(n);
  mv_mlt(mfull, x1, dfull);

  // Calculate mcmpct*x1 = dcmpct
  double dcmpct[n];
  bdspecLUmlt(&mcmpct[0][0], n, lb, ub, x1->ve, dcmpct);

  if(dumpfull) {
    printf("Vector x (random values)\n");
    printf("=======================\n");
    v_out(x1->ve, n);

    printf("Matrix A (random values)\n");
    printf("========================\n");
    printf("Full NxN Meschach Matrix\n");
    m_out(mfull->base, n, n);
    printf("Compact bdspec Array\n");
    m_out(&mcmpct[0][0], n, mm);

    printf("Vector d = A*x\n");
    printf("==============\n");
    printf("Calculated from Full Meschach Matrix:\n");
    v_out(dfull->ve, n);
    printf("Calculated from Compact bdspec Array:\n");
    v_out(dcmpct, n);
    printf("L2 norm of difference between Meschach and bdspec calculations of d\n");
  }
  printf("d diff:%12E\t", v_diff(dfull->ve, dcmpct, n));

  PERM  *p = px_get(n);
  LUfactor(mfull, p);

  int indx[n];
  bdspecLUfactormeschscale(&mcmpct[0][0], n, lb, ub, indx);

  VEC *yfull = v_get(n);
  catchall(LUsolve(mfull, p, dfull, yfull), printf("--matrix singular--:\n"));
  double ycmpct[n];
  for(int i = 0; i < n; i++)
    ycmpct[i] = dcmpct[i];
  bdspecLUsolve(&mcmpct[0][0], n, lb, ub, indx, ycmpct);

  if(dumpfull) {
    printf("\n\n");
    printf("LU Factorization\n");
    printf("================\n");
    printf("Meschach LU Array\n");
    m_out(mfull->base, n, n);
    printf("Compact bdspec LU Array\n");
    m_out(&mcmpct[0][0], n, mm);
    printf("Permutation\n");
    printf("===========\n");
    printf("Meschach permutation vector\n");
    p_out((int *) p->pe, n);
    printf("bdspec indx vector\n");
    p_out(indx, n);

    printf("A*y = d Solved for y\n");
    printf("====================\n");
    printf("Meschach result\n");
    v_out(yfull->ve, n);
    printf("bdspec result\n");
    v_out(ycmpct, n);
    printf("L2 norm of difference between Meschach and bdspec calculations of y:\n");
  }

  printf("y diff:%12E\t", v_diff(yfull->ve, ycmpct, n));
  if(dumpfull) {
    printf("\n\n");
    printf("L2 norm of error = y-x\n");
    printf("======================\n");
  }
  printf("mesch err:%12E\t", v_diff(yfull->ve, x1->ve, n));
  printf("bdspec y err:%12E", v_diff(ycmpct, x1->ve, n));
  if(dumpfull) {
    printf("\n\n");
  }
  fflush(stdout);
}


int main() {
  srand48(time(NULL));

  test(8,1,2,true);

  printf("Testing with a whole lot of numbers...\n");
  printf("(Note that bdspec error is only large when Meschach error is large)\n\n");
  for(int i = 0; i < 4; i++)
    for(int j = 0; j <= i; j++)
      for(int n = 8; n < 2048; n *= 2) {
        printf("n=%4d, lb=%d, ub=%d  ", n, i-j, j);
        test(n, i-j, j, false);
        printf("\n");
      }

  return 0;
}

// vim: ft=cpp:ts=2:sw=2:et:sta:ai:si
