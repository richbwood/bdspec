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

// Required by bdspec:
#include <math.h>
#include <float.h>

typedef struct  {
  double r, i;
} bdspec_complex;

void zbdspecLUmlt(const bdspec_complex *bA, const int n, const int lb, const int ub,
            const bdspec_complex x[], bdspec_complex out[]);
void zbdspecLUfactor(bdspec_complex *bA, const int n, const int lb, const int ub, int indx[]);
void zbdspecLUfactorscale(bdspec_complex *bA, const int n, const int lb, const int ub,
                    int indx[]);
void zbdspecLUfactormeschscale(bdspec_complex *bA, const int n, const int lb, const int ub,
                         int indx[]);
void zbdspecLUsolve(const bdspec_complex *bLU, const int n, const int lb, const int ub,
              const int indx[], bdspec_complex b[]);

// vim: ft=cpp:ts=2:sw=2:et:sta:ai:si
