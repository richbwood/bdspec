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

#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <stdbool.h>
#include <meschach/matrix.h>
#include <meschach/matrix2.h>
#ifndef NAN
#define NAN 0.0/0.0
#endif

#include "bdspec.h"

// Required for testing only:
// Dump matrix 'a' to stdout
void
m_out(const double *a, const int n, const int m) {
  const int d = 3; //digits of accuracy
  const double(*aa)[m] = (const double(*)[m])a;
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
void
v_out(const double a[], const int len) {
  const int d = 3; //digits of accuracy
  for(int i = 0; i < len; i++)
    printf("[%2d] %*.*f\n", i, d+3, d, a[i]);
  printf("\n");
}

// Required for testing only:
// Dump permutation 'p' to stdout
void
p_out(const int p[], const int len) {
  for(int i = 0; i < len; i++)
    printf("[%2d] %2d\n", i, p[i]);
  printf("\n");
}

// Required for testing only:
// Fill 'a' with zeros
void
zero(double a[], const int len) {
  for(int i = 0; i < len; i++)
    a[i] = 0.0;
}

// Required for testing only:
// Fill 'a' with random numbers uniformly distributed between -1 and 1
void
randlist(double a[], const int len) {
  for(int i = 0; i < len; i++)
    a[i] = 2.0*drand48()-1.0;
}

// Required for testing only:
// Calculate the L_2 norm of a-b
double v_diff(const double a[], const double b[], const int len) {
  double acc = 0.0;
  for(int i = 0; i < len; i++) {
    const double diff = a[i]-b[i];
    acc += diff*diff;
  }
  return sqrt(acc);
}

// Required for testing only:
// Calculate the L_2 norm squared of a
double v_normsq(const double a[], const int len) {
  double acc = 0.0;
  for(int i = 0; i < len; i++) {
    acc += a[i]*a[i];
  }
  return acc;
}

// Required for testing only:
// Calculate the L_2 norm of a
double v_norm(const double a[], const int len) {
  return sqrt(v_normsq(a, len));
}

void
test(const int n, const int lb, const int ub, bool dumpfull) {
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
  const int mm = 2*lb+ub+1;
  double mcmpct[n][mm];
  zero(&mcmpct[0][0], n*mm);
  for(int i = 0; i < n; i++)
    for(int j = MAX(i-lb, 0); j < MIN(i+ub+1, n); j++)
      mcmpct[i][j-i+lb+lb] = me[i][j];
  // Replace unused values with NAN to be sure they aren't used
  for(int k = 0; k < n; k++)
    for(int i = 0; i < lb; i++)
      mcmpct[k][i] = NAN;
  for(int k = 0; k < lb; k++)
    for(int i = 0; i < lb-k; i++)
      mcmpct[k][i+lb] = NAN;
  for(int k=n-1; k >= n-ub; k--)
    for(int i = n-1-k+1+lb; i < mm; i++)
      mcmpct[k][i+lb] = NAN;

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
    printf("========================\n");
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
  double ddiff = v_diff(dfull->ve, dcmpct, n);
  printf("d diff=%6.0E ", ddiff);
  if(ddiff*ddiff > DBL_EPSILON*v_normsq(dfull->ve, n))
    printf("FAIL,");
  else
    printf("PASS,");

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

  double ydiff = v_diff(yfull->ve, ycmpct, n);
  printf("y diff=%6.0E ", ydiff);
  if(ydiff*ydiff > DBL_EPSILON*v_normsq(yfull->ve, n))
    printf("FAIL,");
  else
    printf("PASS,");
  if(dumpfull) {
    printf("\n\n");
    printf("L2 norm of error = y-x\n");
    printf("======================\n");
  }
  double x1normsq = v_normsq(x1->ve, n);
  double mescherr = v_diff(yfull->ve, x1->ve, n);
  printf("mesch err=%6.0E ", mescherr);
  if(mescherr*mescherr > DBL_EPSILON*x1normsq)
    printf("FAIL,");
  else
    printf("PASS,");
  double bdspecerr = v_diff(ycmpct, x1->ve, n);
  printf("bdspec err=%6.0E ", bdspecerr);
  if(bdspecerr*bdspecerr > DBL_EPSILON*x1normsq)
    printf("FAIL ");
  else
    printf("PASS ");

  if(dumpfull) {
    printf("\n\n");
  }
  fflush(stdout);
}


int
main(void) {
  srand48(time(NULL));

  test(8, 1, 2, true);

  printf("Testing with a whole lot of numbers...\n");
  printf("(bdspec should only fail when Meschac does)\n\n");
  for(int i = 0; i < 4; i++)
    for(int j = 0; j <= i; j++)
      for(int n = 8; n < 2048; n *= 2) {
        printf("n=%4d,lb=%d,ub=%d,", n, i-j, j);
        test(n, i-j, j, false);
        printf("\n");
      }

  return 0;
}

// vim: ft=cpp:ts=2:sw=2:et:sta:ai:si
