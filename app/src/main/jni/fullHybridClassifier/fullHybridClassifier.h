//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fullHybridClassifier.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Mar-2017 00:59:20
//
#ifndef FULLHYBRIDCLASSIFIER_H
#define FULLHYBRIDCLASSIFIER_H

// Include Files
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "fullHybridClassifier_types.h"

// Variable Declarations
extern omp_nest_lock_t emlrtNestLockGlobal;

// Function Declarations
extern emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);
extern emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
extern emxArray_real_T *emxCreate_real_T(int rows, int cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
extern void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern double fullHybridClassifier(const double ch1_data[], const int ch1_size[1],
  const double ch2_data[], const int ch2_size[1], const double ch3_data[], const
  int ch3_size[1], const double ch4_data[], const int ch4_size[1], const
  emxArray_real_T *tXEOG, const emxArray_real_T *tYEOG, const emxArray_real_T
  *tXSSVEP, const emxArray_real_T *tYSSVEP, double Fs);
extern void fullHybridClassifier_initialize();
extern void fullHybridClassifier_terminate();

#endif

//
// File trailer for fullHybridClassifier.h
//
// [EOF]
//
