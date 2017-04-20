//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eog_knn.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 20-Apr-2017 10:40:49
//
#ifndef EOG_KNN_H
#define EOG_KNN_H

// Include Files
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "eog_knn_types.h"

// Variable Declarations
extern omp_nest_lock_t emlrtNestLockGlobal;

// Function Declarations
extern double eog_knn(const double ch1[250], const double ch2[250], const double
                      ch3[250]);
extern void eog_knn_initialize();
extern void eog_knn_terminate();

#endif

//
// File trailer for eog_knn.h
//
// [EOF]
//
