//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fullHybridClassifier2.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 16-Mar-2017 10:57:23
//
#ifndef FULLHYBRIDCLASSIFIER2_H
#define FULLHYBRIDCLASSIFIER2_H

// Include Files
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "fullHybridClassifier2_types.h"

// Variable Declarations
extern omp_nest_lock_t emlrtNestLockGlobal;

// Function Declarations
extern void fullHybridClassifier2(const double ch1_data[], const int ch1_size[1],
  const double ch2_data[], const int ch2_size[1], const double ch3_data[], const
  int ch3_size[1], const double ch4_data[], const int ch4_size[1], double Fs,
  double Y[5]);
extern void fullHybridClassifier2_initialize();
extern void fullHybridClassifier2_terminate();

#endif

//
// File trailer for fullHybridClassifier2.h
//
// [EOF]
//
