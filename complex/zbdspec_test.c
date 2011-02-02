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

#include <time.h>
#include <stdbool.h>
#include <meschach/zmatrix.h>
#include <meschach/zmatrix2.h>
#ifndef NAN
#define NAN 0.0/0.0
#endif
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#include "zbdspec.h"

// Required for testing only:
// Dump matrix 'a' to stdout
void
zm_out(const double *a, const int n, const int m) {
  const int d = 3; //digits of accuracy
  const double(*aa)[m][2] = (const double(*)[m][2])a;
  for(int i = 0; i < n; i++) {
    printf("[%2d]", i);
    for(int j = 0; j < m; j++) {
      if(FP_ZERO == fpclassify(aa[i][j][0]*aa[i][j][0]+aa[i][j][1]*aa[i][j][1]))
        printf(" %.*s%.*s ", d+3, "                 ", d+3, "                 ");
      else
        printf(" %*.*F%+*.*FI", d+3, d, aa[i][j][0], d+3, d, aa[i][j][1]);
    }
    printf("\n");
  }
  printf("\n");
}

// Required for testing only:
// Dump vector 'a' to stdout
void
zv_out(const double a[], const int len) {
  const int d = 3; //digits of accuracy
  const double(*aa)[2] = (const double(*)[2])a;
  for(int i = 0; i < len; i++)
    printf("[%2d] %*.*F%+*.*FI\n", i, d+3, d, aa[i][0], d+3, d, aa[i][1]);
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
  ZMAT *mfull = zm_get(n, n);
  randlist((double *)mfull->base, (mfull->n)*(mfull->n)*2);
  complex (**me) = mfull->me;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < i - lb; j++) {
      me[i][j].re = 0.0;
      me[i][j].im = 0.0;
    }
    for(int j = i+ub+1; j < n; j++) {
      me[i][j].re = me[i][i+ub].re;
      me[i][j].im = me[i][i+ub].im;
    }
  }

  // Copy matrix mfull to a compactly stored version
  // mcmpct
  // First lb columns padding for later use
  // Next lb columns for lower diagonals
  // Next column for diagonal
  // Next ub columns for upper diagonals
  // as the highest upper diagonal of the same row
  const int mm = 2*lb+ub+1;
  double mcmpct[n][mm][2];
  zero((double *)&mcmpct[0][0], n*mm*2);
  for(int i = 0; i < n; i++)
    for(int j = MAX(i-lb, 0); j < MIN(i+ub+1, n); j++) {
      mcmpct[i][j-i+lb+lb][0] = me[i][j].re;
      mcmpct[i][j-i+lb+lb][1] = me[i][j].im;
    }
  // Replace unused values with NAN to be sure they aren't used
  for(int k = 0; k < n; k++)
    for(int i = 0; i < lb; i++) {
      mcmpct[k][i][0] = NAN;
      mcmpct[k][i][1] = NAN;
    }
  for(int k = 0; k < lb; k++)
    for(int i = 0; i < lb-k; i++) {
      mcmpct[k][i+lb][0] = NAN;
      mcmpct[k][i+lb][1] = NAN;
    }
  for(int k=n-1; k >= n-ub; k--)
    for(int i = n-1-k+1+lb; i < mm; i++) {
      mcmpct[k][i+lb][0] = NAN;
      mcmpct[k][i+lb][1] = NAN;
    }

  // Generate start vector x1 for test
  ZVEC *x1 = zv_get(n);
  randlist((double *)x1->ve, n*2);

  // Calculate mfull*x1 = dfull
  ZVEC *dfull = zv_get(n);
  zmv_mlt(mfull, x1, dfull);

  // Calculate mcmpct*x1 = dcmpct
  double dcmpct[n][2];
  zbdspecLUmlt(&mcmpct[0][0][0], n, lb, ub, (double *)x1->ve, &dcmpct[0][0]);

  if(dumpfull) {
    printf("Vector x (random values)\n");
    printf("========================\n");
    zv_out((double *)x1->ve, n);

    printf("Matrix A (random values)\n");
    printf("========================\n");
    printf("Full NxN Meschach Matrix\n");
    zm_out((double *)mfull->base, n, n);
    printf("Compact bdspec Array\n");
    zm_out(&mcmpct[0][0][0], n, mm);

    printf("Vector d = A*x\n");
    printf("==============\n");
    printf("Calculated from Full Meschach Matrix:\n");
    zv_out((double *)dfull->ve, n);
    printf("Calculated from Compact bdspec Array:\n");
    zv_out((double *)dcmpct, n);
    printf("L2 norm of difference between Meschach and bdspec calculations of d\n");
  }
  double ddiff = v_diff((double *)dfull->ve, &dcmpct[0][0], 2*n);
  printf("d diff=%6.0E ", ddiff);
  if(ddiff*ddiff > DBL_EPSILON*v_normsq((double *)dfull->ve, 2*n))
    printf("FAIL,");
  else
    printf("PASS,");

  PERM  *p = px_get(n);
  zLUfactor(mfull, p);

  int indx[n];
  zbdspecLUfactormeschscale(&mcmpct[0][0][0], n, lb, ub, indx);
  //zbdspecLUfactorscale(&mcmpct[0][0][0], n, lb, ub, indx);
  //zbdspecLUfactor(&mcmpct[0][0][0], n, lb, ub, indx);

  ZVEC *yfull = zv_get(n);
  catchall(zLUsolve(mfull, p, dfull, yfull), printf("--matrix singular--:\n"));
  double ycmpct[n][2];
  for(int i = 0; i < n; i++) {
    ycmpct[i][0] = dcmpct[i][0];
    ycmpct[i][1] = dcmpct[i][1];
  }
  zbdspecLUsolve(&mcmpct[0][0][0], n, lb, ub, indx, &ycmpct[0][0]);

  if(dumpfull) {
    printf("\n\n");
    printf("LU Factorization\n");
    printf("================\n");
    printf("Meschach LU Array\n");
    zm_out((double *)mfull->base, n, n);
    printf("Compact bdspec LU Array\n");
    zm_out(&mcmpct[0][0][0], n, mm);
    printf("Permutation\n");
    printf("===========\n");
    printf("Meschach permutation vector\n");
    p_out((int *) p->pe, n);
    printf("bdspec indx vector\n");
    p_out(indx, n);

    printf("A*y = d Solved for y\n");
    printf("====================\n");
    printf("Meschach result\n");
    zv_out((double *)yfull->ve, n);
    printf("bdspec result\n");
    zv_out(&ycmpct[0][0], n);
    printf("L2 norm of difference between Meschach and bdspec calculations of y:\n");
  }

  double ydiff = v_diff((double *)yfull->ve, &ycmpct[0][0], 2*n);
  printf("y diff=%6.0E ", ydiff);
  if(ydiff*ydiff > DBL_EPSILON*v_normsq((double *)yfull->ve, 2*n))
    printf("FAIL,");
  else
    printf("PASS,");
  if(dumpfull) {
    printf("\n\n");
    printf("L2 norm of error = y-x\n");
    printf("======================\n");
  }
  double x1normsq = v_normsq((double *)x1->ve, 2*n);
  double mescherr = v_diff((double *)yfull->ve, (double *)x1->ve, 2*n);
  printf("mesch err=%6.0E ", mescherr);
  if(mescherr*mescherr > DBL_EPSILON*x1normsq)
    printf("FAIL,");
  else
    printf("PASS,");
  double bdspecerr = v_diff(&ycmpct[0][0], (double *)x1->ve, 2*n);
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
