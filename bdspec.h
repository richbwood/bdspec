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
#define SWAP(a, b) do {const double dum = (a); (a) = (b); (b) = dum;} while(0)
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void bdspecLUmlt(const double *bA, const int n, const int lb, const int ub,
            const double x[], double out[]);
void bdspecLUfactor(double *bA, const int n, const int lb, const int ub, int indx[]);
void bdspecLUfactorscale(double *bA, const int n, const int lb, const int ub,
                    int indx[]);
void bdspecLUfactormeschscale(double *bA, const int n, const int lb, const int ub,
                         int indx[]);
void bdspecLUsolve(const double *bLU, const int n, const int lb, const int ub,
              const int indx[], double b[]);

// vim: ft=cpp:ts=2:sw=2:et:sta:ai:si
