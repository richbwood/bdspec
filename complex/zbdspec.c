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

#include "zbdspec.h"

#define SWAP(a, b) do {const double dum = (a); (a) = (b); (b) = dum;} while(0)
#define ZSWAP(a, b) do {const bdspec_complex dum = (a); (a) = (b); (b) = dum;} while(0)
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define cabssqr(a) ((a).r*(a).r+(a).i*(a).i)

// Calculates bA*x for (n*n) special-band-matrix bA and (n) vector x.
// -- the (n*n) matrix A is stored in band storage array bA of dimension
// (n*[2*lb+1+ub]). Elements are stored in columns lb to 2*lb+ub.  The jth row
// of bA is stored in the jth row of band A (bA) as follows: bA[i][2*lb+j-i] =
// A[i][j] for max(0, i-lb) <= j <= min(A->n-1, i+ub); Note that the first lb
// columns of bA are unused!
void
zbdspecLUmlt(const bdspec_complex *bA, const int n, const int lb, const int ub,
            const bdspec_complex x[], bdspec_complex out[]) {
  const int mm = 2*lb+ub+1;
  const bdspec_complex(*a)[mm] = (const bdspec_complex(*)[mm])bA;
  bdspec_complex acc = {0};
  for(int i = n-1; i >= 0; i--) {
    {
      const int k = i-2*lb;
      const int j0 = MAX(lb, -k);
      const int j1 = MIN(mm-1, n-k);
      out[i] = (bdspec_complex){0};
      for(int j = j0, l = j0+k; j < j1; j++, l++) {
        out[i].r += a[i][j].r*x[l].r-a[i][j].i*x[l].i;
        out[i].i += a[i][j].r*x[l].i+a[i][j].i*x[l].r;
      }
    }
    {
      const int k = ub+i;
      if(k < n) {
        acc.r += x[k].r;
        acc.i += x[k].i;
        out[i].r += a[i][mm-1].r*acc.r-a[i][mm-1].i*acc.i;
        out[i].i += a[i][mm-1].r*acc.i+a[i][mm-1].i*acc.r;
      }
    }
  }
}

// Gaussian elimination with partial pivoting
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
void
zbdspecLUfactor(bdspec_complex *bA, const int n, const int lb, const int ub, int indx[]) {
  const int mm = 2*lb+ub+1;
  bdspec_complex(*a)[mm] = (bdspec_complex(*)[mm])bA;
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
    double maxa0 = cabssqr(a[k][lb]);
    int i = k;
    for(int j = k+1; j < l; j++) {
      const double fabsa = cabssqr(a[j][lb]);
      if(fabsa > maxa0) {
        maxa0 = fabsa;
        i = j;
      }
    }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(0.0 == maxa0) // If singular, give up and go to the next row
      continue;
    // Interchange rows
    if(i != k)
      for(int j = lb; j < mm; j++) {
        ZSWAP(a[k][j], a[i][j]);
      }
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      bdspec_complex fac;
      fac.r = (a[i][lb].r*a[k][lb].r+a[i][lb].i*a[k][lb].i)/maxa0;
      fac.i = (a[i][lb].i*a[k][lb].r-a[i][lb].r*a[k][lb].i)/maxa0;
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++) {
        a[i][j-1].r = a[i][j].r-fac.r*a[k][j].r+fac.i*a[k][j].i;
        a[i][j-1].i = a[i][j].i-fac.r*a[k][j].i-fac.i*a[k][j].r;
      }
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// Same as bdspecLUfactor, but pivots on rows that are normalized such that the
// greatest element in each row is +1 or -1. Note that this only affects which
// elements are pivoted: the resulting LU factorization is not normalized.
void
zbdspecLUfactorscale(bdspec_complex *bA, const int n, const int lb, const int ub,
                    int indx[]) {
  const int mm = 2*lb+ub+1;
  bdspec_complex(*a)[mm] = (bdspec_complex(*)[mm])bA;
  double scale[n];
  // Calculate the scale of each row
  for(int i = 0; i < n; i++) {
    double max = 0.0;
    const int k = i-2*lb;
    const int j0 = MAX(lb, -k);
    const int j1 = MIN(mm, n-k);
    for(int j = j0; j < j1; j++)
      max = fmax(cabssqr(a[i][j]), max);
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
      const double fabsa = cabssqr(a[j][lb]);
      if(scale[j] >= DBL_EPSILON*DBL_EPSILON*fabsa) {
        const double temp = fabsa/scale[j];
        if(temp > maxa0) {
          maxa0 = temp;
          i = j;
        }
      }
    }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(-1 == i) { // If singular continue to next row
      indx[k] = k;
      a[k][lb] = (bdspec_complex){0.0};
      continue;
    }
    // Interchange rows
    if(i != k) {
      SWAP(scale[k], scale[i]);
      for(int j = lb; j < mm; j++) {
        ZSWAP(a[k][j], a[i][j]);
      }
    }
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      const double denom = maxa0*scale[k];
      bdspec_complex fac;
      fac.r = (a[i][lb].r*a[k][lb].r+a[i][lb].i*a[k][lb].i)/denom;
      fac.i = (a[i][lb].i*a[k][lb].r-a[i][lb].r*a[k][lb].i)/denom;
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++) {
        a[i][j-1].r = a[i][j].r-fac.r*a[k][j].r+fac.i*a[k][j].i;
        a[i][j-1].i = a[i][j].i-fac.r*a[k][j].i-fac.i*a[k][j].r;
      }
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// Same as bdspecLUfactorscale, but includes a bug found in meschach 1.2.
// Having identical bugs allows testing of the LU factorization
void
zbdspecLUfactormeschscale(bdspec_complex *bA, const int n, const int lb, const int ub,
                         int indx[]) {
  const int mm = 2*lb+ub+1;
  bdspec_complex(*a)[mm] = (bdspec_complex(*)[mm])bA;
  double scale[n];
  // Calculate the scale of each row
  for(int i = 0; i < n; i++) {
    double max = 0.0;
    const int k = i-2*lb;
    const int j0 = MAX(lb, -k);
    const int j1 = MIN(mm, n-k);
    for(int j = j0; j < j1; j++)
      max = fmax(cabssqr(a[i][j]), max);
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
      const double fabsa = cabssqr(a[j][lb]);
      if(scale[j] >= DBL_EPSILON*DBL_EPSILON*fabsa) {
        const double temp = fabsa/scale[j];
        if(temp > maxa0) {
          maxa0 = temp;
          i = j;
        }
      }
    }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(-1 == i) { // If singular continue to next row
      indx[k] = k;
      a[k][lb] = (bdspec_complex){0.0};
      continue;
    }
    // Interchange rows
    if(i != k) {
      // Meschach doesn't have this:
      //SWAP(scale[k], scale[i]);
      for(int j = lb; j < mm; j++) {
        ZSWAP(a[k][j], a[i][j]);
      }
    }
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      // Meschach above requires modified code here:
      //const double denom = maxa0*scale[k];
      const double denom = cabssqr(a[k][lb]);
      bdspec_complex fac;
      fac.r = (a[i][lb].r*a[k][lb].r+a[i][lb].i*a[k][lb].i)/denom;
      fac.i = (a[i][lb].i*a[k][lb].r-a[i][lb].r*a[k][lb].i)/denom;
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++) {
        a[i][j-1].r = a[i][j].r-fac.r*a[k][j].r+fac.i*a[k][j].i;
        a[i][j-1].i = a[i][j].i-fac.r*a[k][j].i-fac.i*a[k][j].r;
      }
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// given an LU factorisation in bA, solve bA*x=b for x
// -- Solution is stored in-place in array b.
// -- L and U are stored together in array bLU.  Upper triangular matrix U is
// stored with diagonal and lb+ub super-diagonals in columns lb to 2*lb+ub. The
// diagonals of lower triangular matrix L are ones, other nonzero elements of
// lower triangular matrix L is stored in columns 0 to lb-1.
// -- Array indx stores the row reordering required
void
zbdspecLUsolve(const bdspec_complex *bLU, const int n, const int lb, const int ub,
              const int indx[], bdspec_complex b[]) {
  const int mm = 2*lb+ub+1;
  const bdspec_complex(*a)[mm] = (const bdspec_complex(*)[mm])bLU;
  // Forward substitution:
  for(int k = 0, l = lb+1; k < n; k++, l++) {
    l = MIN(l, n);
    // b must be reordered during forward substitution, not before
    // quicker to swap without testing
    ZSWAP(b[k], b[indx[k]]);
    for(int i = k+1, j = 0; i < l; i++, j++) {
      b[i].r -= a[k][j].r*b[k].r-a[k][j].i*b[k].i;;
      b[i].i -= a[k][j].r*b[k].i+a[k][j].i*b[k].r;;
    }
  }
  // Backward substitution:
  bdspec_complex acc2 = {0};
  for(int i = n-1, l = 1+lb; i >= 0; i--, l++) {
    l = MIN(l, mm-1);
    bdspec_complex acc;
    acc = b[i];
    for(int k = lb+1, j = i+1; k < l; k++, j++) {
      acc.r -= a[i][k].r*b[j].r-a[i][k].i*b[j].i;
      acc.i -= a[i][k].r*b[j].i+a[i][k].i*b[j].r;
    }
    {
      const int k = MAX(ub+lb, 1)+i;
      if(k < n) {
        acc2.r += b[k].r;
        acc2.i += b[k].i;
        acc.r -= a[i][mm-1].r*acc2.r-a[i][mm-1].i*acc2.i;
        acc.i -= a[i][mm-1].r*acc2.i+a[i][mm-1].i*acc2.r;
      }
    }
    const double denom = a[i][lb].r*a[i][lb].r+a[i][lb].i*a[i][lb].i;
    b[i].r = (acc.r*a[i][lb].r+acc.i*a[i][lb].i)/denom;
    b[i].i = (acc.i*a[i][lb].r-acc.r*a[i][lb].i)/denom;
  }
}

// vim: ft=cpp:ts=2:sw=2:et:sta:ai:si
