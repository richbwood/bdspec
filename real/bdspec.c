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
// This software is provided by Richard Wood ``as is'' and any express or
// implied warranties, including, but not limited to, the implied warranties of
// merchantability and fitness for a particular purpose are disclaimed. In no
// event shall Richard Wood or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused and
// on any theory of liability, whether in contract, strict liability, or tort
// (including negligence or otherwise) arising in any way out of the use of
// this software, even if advised of the possibility of such damage.
// 
// The views and conclusions contained in the software and documentation are
// those of the authors and should not be interpreted as representing official
// policies, either expressed or implied, of Richard Wood.

#include "bdspec.h"

// Calculates bA*x for (n*n) special-band-matrix bA and (n) vector x.
// -- the (n*n) matrix A is stored in band storage array bA of dimension
// (n*[2*lb+1+ub]). Elements are stored in columns lb to 2*lb+ub.  The jth row
// of bA is stored in the jth row of band A (bA) as follows: bA[i][2*lb+j-i] =
// A[i][j] for max(0, i-lb) <= j <= min(A->n-1, i+ub); Note that the first lb
// columns of bA are unused!
void
bdspecLUmlt(const double *bA, const int n, const int lb, const int ub,
            const double x[], double out[]) {
  const int mm = 2*lb+ub+1;
  const double(*a)[mm] = (const double(*)[mm])bA;
  double acc = 0.0;
  for(int i = n-1; i >= 0; i--) {
    {
      const int k = i-2*lb;
      const int j0 = MAX(lb, -k);
      const int j1 = MIN(mm-1, n-k);
      out[i] = 0.0;
      for(int j = j0, l = j0+k; j < j1; j++, l++)
        out[i] += a[i][j]*x[l];
    }
    {
      const int k = ub+i;
      if(k < n) {
        acc += x[k];
        out[i] += a[i][mm-1]*acc;
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
bdspecLUfactor(double *bA, const int n, const int lb, const int ub, int indx[]) {
  const int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
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
      const double fabsa = fabs(a[j][lb]);
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
      for(int j = lb; j < mm; j++)
        SWAP(a[k][j], a[i][j]);
    // Perform LU factorisation
    for(i = k+1; i < l; i++) {
      const double fac = a[i][lb]/a[k][lb];
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
void
bdspecLUfactorscale(double *bA, const int n, const int lb, const int ub,
                    int indx[]) {
  const int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
  double scale[n];
  // Ensure repeated columns are repeated
  for(int k = 0; k < n; k++)
    for(int j = mm; j < mm; j++)
      a[k][j] = a[k][mm-1];
  // Calculate the scale of each row
  for(int i = 0; i < n; i++) {
    double max = 0.0;
    const int k = i-2*lb;
    const int j0 = MAX(lb, -k);
    const int j1 = MIN(mm, n-k);
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
    for(int j = k; j < l; j++)
      if(scale[j] >= DBL_EPSILON*fabs(a[j][lb])) {
        const double temp = fabs(a[j][lb])/scale[j];
        if(temp > maxa0) {
          maxa0 = temp;
          i = j;
        }
      }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(-1 == i) { // If singular continue to next row
      indx[k] = k;
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
      const double fac = a[i][lb]/a[k][lb];
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++)
        a[i][j-1] = a[i][j]-fac*a[k][j];
      a[i][mm-1] = a[i][mm-2];
    }
  }
}

// Same as bdspecLUfactorscale, but includes a big found in meschach 1.2.
// Having identical bugs allows testing of the LU factorization
void
bdspecLUfactormeschscale(double *bA, const int n, const int lb, const int ub,
                         int indx[]) {
  const int mm = 2*lb+ub+1;
  double(*a)[mm] = (double(*)[mm])bA;
  double scale[n];
  // Calculate the scale of each row
  for(int i = 0; i < n; i++) {
    double max = 0.0;
    const int k = i-2*lb;
    const int j0 = MAX(lb, -k);
    const int j1 = MIN(mm, n-k);
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
    for(int j = k; j < l; j++)
      if(scale[j] >= DBL_EPSILON*fabs(a[j][lb])) {
        const double temp = fabs(a[j][lb])/scale[j];
        if(temp > maxa0) {
          maxa0 = temp;
          i = j;
        }
      }
    // Store pivot row index:
    indx[k] = i;
    // Check if singular
    if(-1 == i) { // If singular continue to next row
      indx[k] = k;
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
      const double fac = a[i][lb]/a[k][lb];
      a[k][i-k-1] = fac;
      for(int j = 1+lb; j < mm; j++)
        a[i][j-1] = a[i][j]-fac*a[k][j];
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
bdspecLUsolve(const double *bLU, const int n, const int lb, const int ub,
              const int indx[], double b[]) {
  const int mm = 2*lb+ub+1;
  const double(*a)[mm] = (const double(*)[mm])bLU;
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
    //for(int k = lb+1, j = i+1; j < n; k++, j++) {
      acc -= a[i][k]*b[j];
    }
    {
      const int j = MAX(ub+lb, 1)+i;
      if(j < n) {
        acc2 += b[j];
        acc -= a[i][mm-1]*acc2;
      }
    }
    b[i] = acc/a[i][lb];
  }
}

// vim: ft=cpp:ts=2:sw=2:et:sta:ai:si
