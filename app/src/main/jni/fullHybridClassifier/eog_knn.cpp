//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eog_knn.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 20-Apr-2017 10:40:49
//

// Include Files
#include "rt_nonfinite.h"
#include "eog_knn.h"
#include <cstdlib>
#include <math.h>
#include <android/log.h>

#define  LOG_TAG "eog_knn-cpp"
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)
// Type Definitions
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common {
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray__common

#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

struct emxArray_boolean_T
{
  boolean_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_boolean_T

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_int32_T

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T

// Variable Definitions
//omp_nest_lock_t emlrtNestLockGlobal;

// Function Declarations
static void b_diff(const double x[250], double y[249]);
static void b_findLocalMaxima(const double yTemp[250], emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void b_findpeaks(const double Yin[250], emxArray_real_T *Ypk,
  emxArray_real_T *Xpk);
static void b_floor(emxArray_real_T *x);
static void b_merge(int idx[123], double x[123], int offset, int np, int nq, int
                    iwork[123], double xwork[123]);
static void b_merge_block(int idx[123], double x[123], int offset, int n, int
  preSortLevel, int iwork[123], double xwork[123]);
static void b_sign(emxArray_real_T *x);
static void b_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx);
static double b_xnrm2(int n, const emxArray_real_T *x, int ix0);
static void c_diff(const emxArray_real_T *x, emxArray_real_T *y);
static void c_sort(double x[123], int idx[123]);
static void combinePeaks(const emxArray_real_T *iPk, const emxArray_real_T *iInf,
  emxArray_real_T *iPkOut);
static void customFilt(const double X[250], double Y[250]);
static void cwt_haar(const double SIG[250], double coefs[5000]);
static void diff(const double x[250], double y[249]);
static void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static double eps(double x);
static void featureExtractionEOG(const double samplesX[250], emxArray_real_T *F);
static void filter(double b[5], double a[5], const double x[274], const double
                   zi[4], double y[274]);
static void findLocalMaxima(const double yTemp[249], emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void findpeaks(const double Yin[249], emxArray_real_T *Ypk,
                      emxArray_real_T *Xpk);
static int ixamax(int n, double x);
static double knn(const emxArray_real_T *tsX);
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void mrdivide(emxArray_real_T *A, const emxArray_real_T *B);
static void power(const double a[249], double y[249]);
static double rt_hypotd_snf(double u0, double u1);
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);
static void sort(emxArray_real_T *x, emxArray_int32_T *idx);
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
static double trapz(const double x[250]);
static void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, int *jpvt);
static double xnrm2(int n, const emxArray_real_T *x, int ix0);
static void xscal(int n, double a, emxArray_real_T *x, int ix0);

// Function Definitions

//
// Arguments    : const double x[250]
//                double y[249]
// Return Type  : void
//
static void b_diff(const double x[250], double y[249])
{
  int ixLead;
  int iyLead;
  double work;
  int m;
  double tmp2;
  ixLead = 1;
  iyLead = 0;
  work = x[0];
  for (m = 0; m < 249; m++) {
    tmp2 = work;
    work = x[ixLead];
    tmp2 = x[ixLead] - tmp2;
    ixLead++;
    y[iyLead] = tmp2;
    iyLead++;
  }
}

//
// Arguments    : const double yTemp[250]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void b_findLocalMaxima(const double yTemp[250], emxArray_real_T *iPk,
  emxArray_real_T *iInflect)
{
  double b_yTemp[252];
  boolean_T yFinite[252];
  int ii;
  boolean_T x[251];
  emxArray_int32_T *b_ii;
  int idx;
  int i1;
  boolean_T exitg3;
  emxArray_int32_T *r2;
  boolean_T guard3 = false;
  emxArray_real_T *iTemp;
  emxArray_real_T *c_yTemp;
  emxArray_real_T *s;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r3;
  int nx;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *c_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  b_yTemp[0] = rtNaN;
  memcpy(&b_yTemp[1], &yTemp[0], 250U * sizeof(double));
  b_yTemp[251] = rtNaN;
  for (ii = 0; ii < 252; ii++) {
    yFinite[ii] = !rtIsNaN(b_yTemp[ii]);
  }

  for (ii = 0; ii < 251; ii++) {
    x[ii] = ((b_yTemp[ii] != b_yTemp[ii + 1]) && (yFinite[ii] || yFinite[ii + 1]));
  }

  emxInit_int32_T(&b_ii, 1);
  idx = 0;
  i1 = b_ii->size[0];
  b_ii->size[0] = 251;
  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 252)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= 251) {
        exitg3 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3) {
      ii++;
    }
  }

  emxInit_int32_T(&r2, 1);
  i1 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  i1 = r2->size[0];
  r2->size[0] = 1 + b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r2, i1, (int)sizeof(int));
  r2->data[0] = 1;
  ii = b_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    r2->data[i1 + 1] = b_ii->data[i1] + 1;
  }

  emxInit_real_T1(&iTemp, 1);
  i1 = iTemp->size[0];
  iTemp->size[0] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i1, (int)sizeof(double));
  ii = r2->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iTemp->data[i1] = 1.0 + (double)(r2->data[i1] - 1);
  }

  emxFree_int32_T(&r2);
  emxInit_real_T1(&c_yTemp, 1);
  i1 = c_yTemp->size[0];
  c_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)c_yTemp, i1, (int)sizeof(double));
  ii = iTemp->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    c_yTemp->data[i1] = b_yTemp[(int)iTemp->data[i1] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_boolean_T(&b_x, 1);
  emxInit_real_T1(&r3, 1);
  c_diff(c_yTemp, s);
  b_sign(s);
  c_diff(s, r3);
  i1 = b_x->size[0];
  b_x->size[0] = r3->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i1, (int)sizeof(boolean_T));
  ii = r3->size[0];
  emxFree_real_T(&c_yTemp);
  for (i1 = 0; i1 < ii; i1++) {
    b_x->data[i1] = (r3->data[i1] < 0.0);
  }

  emxFree_real_T(&r3);
  nx = b_x->size[0];
  idx = 0;
  i1 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  ii = 1;
  exitg2 = false;
  while ((!exitg2) && (ii <= nx)) {
    guard2 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg2 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i1 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
    }
  } else {
    i1 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  }

  if (1.0 > (double)s->size[0] - 1.0) {
    ii = 0;
  } else {
    ii = (int)((double)s->size[0] - 1.0);
  }

  if (2 > s->size[0]) {
    i1 = 0;
  } else {
    i1 = 1;
  }

  idx = b_x->size[0];
  b_x->size[0] = ii;
  emxEnsureCapacity((emxArray__common *)b_x, idx, (int)sizeof(boolean_T));
  for (idx = 0; idx < ii; idx++) {
    b_x->data[idx] = (s->data[idx] != s->data[i1 + idx]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i1 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i1, (int)sizeof(int));
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx)) {
    guard1 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      c_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i1 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i1, (int)sizeof(int));
    }
  } else {
    i1 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i1, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  i1 = iInflect->size[0];
  iInflect->size[0] = c_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i1, (int)sizeof(double));
  ii = c_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iInflect->data[i1] = iTemp->data[c_ii->data[i1]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  i1 = iPk->size[0];
  iPk->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i1, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iPk->data[i1] = iTemp->data[b_ii->data[i1]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : const double Yin[250]
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void b_findpeaks(const double Yin[250], emxArray_real_T *Ypk,
  emxArray_real_T *Xpk)
{
  boolean_T x[250];
  int k;
  emxArray_int32_T *ii;
  int idx;
  int cdiff;
  boolean_T exitg1;
  emxArray_real_T *iInfite;
  boolean_T guard1 = false;
  double yTemp[250];
  emxArray_real_T *iPk;
  emxArray_real_T *b_idx;
  int ndbl;
  emxArray_real_T *base;
  double extremum;
  emxArray_real_T *varargin_2;
  int apnd;
  emxArray_real_T *y;
  emxArray_real_T *c_idx;
  for (k = 0; k < 250; k++) {
    x[k] = (rtIsInf(Yin[k]) && (Yin[k] > 0.0));
  }

  emxInit_int32_T(&ii, 1);
  idx = 0;
  k = ii->size[0];
  ii->size[0] = 250;
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = 1;
  exitg1 = false;
  while ((!exitg1) && (cdiff < 251)) {
    guard1 = false;
    if (x[cdiff - 1]) {
      idx++;
      ii->data[idx - 1] = cdiff;
      if (idx >= 250) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      cdiff++;
    }
  }

  emxInit_real_T1(&iInfite, 1);
  k = ii->size[0];
  if (1 > idx) {
    ii->size[0] = 0;
  } else {
    ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  k = iInfite->size[0];
  iInfite->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInfite, k, (int)sizeof(double));
  cdiff = ii->size[0];
  for (k = 0; k < cdiff; k++) {
    iInfite->data[k] = ii->data[k];
  }

  memcpy(&yTemp[0], &Yin[0], 250U * sizeof(double));
  k = ii->size[0];
  ii->size[0] = iInfite->size[0];
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = iInfite->size[0];
  for (k = 0; k < cdiff; k++) {
    ii->data[k] = (int)iInfite->data[k];
  }

  cdiff = ii->size[0];
  for (k = 0; k < cdiff; k++) {
    yTemp[ii->data[k] - 1] = rtNaN;
  }

  emxFree_int32_T(&ii);
  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&b_idx, 1);
  b_findLocalMaxima(yTemp, iPk, b_idx);
  if (!(iPk->size[0] == 0)) {
    cdiff = iPk->size[0] - 1;
    ndbl = 0;
    for (idx = 0; idx <= cdiff; idx++) {
      if (Yin[(int)iPk->data[idx] - 1] > 0.001) {
        ndbl++;
      }
    }

    k = 0;
    for (idx = 0; idx <= cdiff; idx++) {
      if (Yin[(int)iPk->data[idx] - 1] > 0.001) {
        iPk->data[k] = iPk->data[idx];
        k++;
      }
    }

    k = iPk->size[0];
    iPk->size[0] = ndbl;
    emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  }

  cdiff = iPk->size[0];
  emxInit_real_T1(&base, 1);
  k = base->size[0];
  base->size[0] = cdiff;
  emxEnsureCapacity((emxArray__common *)base, k, (int)sizeof(double));
  for (k = 0; k + 1 <= cdiff; k++) {
    if ((Yin[(int)(iPk->data[k] - 1.0) - 1] >= Yin[(int)(iPk->data[k] + 1.0) - 1])
        || rtIsNaN(Yin[(int)(iPk->data[k] + 1.0) - 1])) {
      extremum = Yin[(int)(iPk->data[k] - 1.0) - 1];
    } else {
      extremum = Yin[(int)(iPk->data[k] + 1.0) - 1];
    }

    base->data[k] = extremum;
  }

  cdiff = iPk->size[0] - 1;
  ndbl = 0;
  for (idx = 0; idx <= cdiff; idx++) {
    if (Yin[(int)iPk->data[idx] - 1] - base->data[idx] >= 0.0) {
      ndbl++;
    }
  }

  k = 0;
  for (idx = 0; idx <= cdiff; idx++) {
    if (Yin[(int)iPk->data[idx] - 1] - base->data[idx] >= 0.0) {
      iPk->data[k] = iPk->data[idx];
      k++;
    }
  }

  emxFree_real_T(&base);
  emxInit_real_T1(&varargin_2, 1);
  k = iPk->size[0];
  iPk->size[0] = ndbl;
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  combinePeaks(iPk, iInfite, varargin_2);
  emxFree_real_T(&iInfite);
  if (varargin_2->size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)varargin_2->size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - varargin_2->size[0]) + 1;
    idx = varargin_2->size[0];
    if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)idx) {
      ndbl++;
      apnd = varargin_2->size[0];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }
  }

  emxInit_real_T(&y, 2);
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  if (ndbl > 0) {
    y->data[0] = 1.0;
    if (ndbl > 1) {
      y->data[ndbl - 1] = apnd;
      idx = (ndbl - 1) / 2;
      for (k = 1; k < idx; k++) {
        y->data[k] = 1.0 + (double)k;
        y->data[(ndbl - k) - 1] = apnd - k;
      }

      if (idx << 1 == ndbl - 1) {
        y->data[idx] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[idx] = 1.0 + (double)idx;
        y->data[idx + 1] = apnd - idx;
      }
    }
  }

  k = b_idx->size[0];
  b_idx->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_idx, k, (int)sizeof(double));
  cdiff = y->size[1];
  for (k = 0; k < cdiff; k++) {
    b_idx->data[k] = y->data[y->size[0] * k];
  }

  emxFree_real_T(&y);
  if (b_idx->size[0] > 250) {
    emxInit_real_T1(&c_idx, 1);
    k = c_idx->size[0];
    c_idx->size[0] = 250;
    emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
    for (k = 0; k < 250; k++) {
      c_idx->data[k] = b_idx->data[k];
    }

    k = b_idx->size[0];
    b_idx->size[0] = c_idx->size[0];
    emxEnsureCapacity((emxArray__common *)b_idx, k, (int)sizeof(double));
    cdiff = c_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      b_idx->data[k] = c_idx->data[k];
    }

    emxFree_real_T(&c_idx);
  }

  k = iPk->size[0];
  iPk->size[0] = b_idx->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  cdiff = b_idx->size[0];
  for (k = 0; k < cdiff; k++) {
    iPk->data[k] = varargin_2->data[(int)b_idx->data[k] - 1];
  }

  emxFree_real_T(&varargin_2);
  emxFree_real_T(&b_idx);
  k = Ypk->size[0] * Ypk->size[1];
  Ypk->size[0] = 1;
  Ypk->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Ypk, k, (int)sizeof(double));
  cdiff = iPk->size[0];
  for (k = 0; k < cdiff; k++) {
    Ypk->data[Ypk->size[0] * k] = Yin[(int)iPk->data[k] - 1];
  }

  k = Xpk->size[0] * Xpk->size[1];
  Xpk->size[0] = 1;
  Xpk->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Xpk, k, (int)sizeof(double));
  cdiff = iPk->size[0];
  for (k = 0; k < cdiff; k++) {
    Xpk->data[Xpk->size[0] * k] = 1.0 + (double)((int)iPk->data[k] - 1);
  }

  emxFree_real_T(&iPk);
}

//
// Arguments    : emxArray_real_T *x
// Return Type  : void
//
static void b_floor(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = std::floor(x->data[k]);
  }
}

//
// Arguments    : int idx[123]
//                double x[123]
//                int offset
//                int np
//                int nq
//                int iwork[123]
//                double xwork[123]
// Return Type  : void
//
static void b_merge(int idx[123], double x[123], int offset, int np, int nq, int
                    iwork[123], double xwork[123])
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if ((np == 0) || (nq == 0)) {
  } else {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork[qend] = idx[offset + qend];
      xwork[qend] = x[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[n]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[n];
        x[iout] = xwork[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = iout - p;
          while (p + 1 <= np) {
            idx[(n + p) + 1] = iwork[p];
            x[(n + p) + 1] = xwork[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : int idx[123]
//                double x[123]
//                int offset
//                int n
//                int preSortLevel
//                int iwork[123]
//                double xwork[123]
// Return Type  : void
//
static void b_merge_block(int idx[123], double x[123], int offset, int n, int
  preSortLevel, int iwork[123], double xwork[123])
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        b_merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      b_merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork,
              xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    b_merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : emxArray_real_T *x
// Return Type  : void
//
static void b_sign(emxArray_real_T *x)
{
  int nx;
  int k;
  double b_x;
  nx = x->size[0];
  for (k = 0; k + 1 <= nx; k++) {
    if (x->data[k] < 0.0) {
      b_x = -1.0;
    } else if (x->data[k] > 0.0) {
      b_x = 1.0;
    } else if (x->data[k] == 0.0) {
      b_x = 0.0;
    } else {
      b_x = x->data[k];
    }

    x->data[k] = b_x;
  }
}

//
// Arguments    : emxArray_real_T *x
//                int dim
//                emxArray_int32_T *idx
// Return Type  : void
//
static void b_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
  int i3;
  emxArray_real_T *vwork;
  int k;
  unsigned int unnamed_idx_0;
  int vstride;
  emxArray_int32_T *iidx;
  int j;
  if (dim <= 1) {
    i3 = x->size[0];
  } else {
    i3 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  k = vwork->size[0];
  vwork->size[0] = i3;
  emxEnsureCapacity((emxArray__common *)vwork, k, (int)sizeof(double));
  unnamed_idx_0 = (unsigned int)x->size[0];
  k = idx->size[0];
  idx->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, k, (int)sizeof(int));
  if (dim > 2) {
    vstride = x->size[0];
  } else {
    vstride = 1;
    k = 1;
    while (k <= dim - 1) {
      k = x->size[0];
      vstride *= k;
      k = 2;
    }
  }

  emxInit_int32_T(&iidx, 1);
  for (j = 0; j + 1 <= vstride; j++) {
    for (k = 0; k + 1 <= i3; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, iidx);
    for (k = 0; k + 1 <= i3; k++) {
      x->data[j + k * vstride] = vwork->data[k];
      idx->data[j + k * vstride] = iidx->data[k];
    }
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

//
// Arguments    : int n
//                const emxArray_real_T *x
//                int ix0
// Return Type  : double
//
static double b_xnrm2(int n, const emxArray_real_T *x, int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x->data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x->data[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * std::sqrt(y);
  }

  return y;
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_real_T *y
// Return Type  : void
//
static void c_diff(const emxArray_real_T *x, emxArray_real_T *y)
{
  int iyLead;
  int orderForDim;
  emxArray_real_T *work;
  int ySize_idx_0;
  int m;
  double tmp1;
  int k;
  double tmp2;
  if (x->size[0] == 0) {
    iyLead = y->size[0];
    y->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
  } else {
    if (x->size[0] - 1 <= 1) {
      orderForDim = x->size[0] - 1;
    } else {
      orderForDim = 1;
    }

    if (orderForDim < 1) {
      iyLead = y->size[0];
      y->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
    } else {
      emxInit_real_T1(&work, 1);
      ySize_idx_0 = x->size[0] - orderForDim;
      iyLead = work->size[0];
      work->size[0] = orderForDim;
      emxEnsureCapacity((emxArray__common *)work, iyLead, (int)sizeof(double));
      iyLead = y->size[0];
      y->size[0] = ySize_idx_0;
      emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
      if (!(y->size[0] == 0)) {
        ySize_idx_0 = 1;
        iyLead = 0;
        work->data[0] = x->data[0];
        if (orderForDim >= 2) {
          for (m = 1; m < orderForDim; m++) {
            tmp1 = x->data[ySize_idx_0];
            for (k = 0; k + 1 <= m; k++) {
              tmp2 = work->data[k];
              work->data[k] = tmp1;
              tmp1 -= tmp2;
            }

            work->data[m] = tmp1;
            ySize_idx_0++;
          }
        }

        for (m = orderForDim + 1; m <= x->size[0]; m++) {
          tmp1 = x->data[ySize_idx_0];
          for (k = 0; k + 1 <= orderForDim; k++) {
            tmp2 = work->data[k];
            work->data[k] = tmp1;
            tmp1 -= tmp2;
          }

          ySize_idx_0++;
          y->data[iyLead] = tmp1;
          iyLead++;
        }
      }

      emxFree_real_T(&work);
    }
  }
}

//
// Arguments    : double x[123]
//                int idx[123]
// Return Type  : void
//
static void c_sort(double x[123], int idx[123])
{
  double x4[4];
  signed char idx4[4];
  int m;
  double xwork[123];
  int nNaNs;
  int ib;
  int k;
  signed char perm[4];
  int i2;
  int iwork[123];
  int i3;
  int i4;
  memset(&idx[0], 0, 123U * sizeof(int));
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  memset(&xwork[0], 0, 123U * sizeof(double));
  nNaNs = -122;
  ib = 0;
  for (k = 0; k < 123; k++) {
    if (rtIsNaN(x[k])) {
      idx[-nNaNs] = k + 1;
      xwork[-nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = (signed char)(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        ib = (k - nNaNs) - 125;
        if (x4[0] <= x4[1]) {
          m = 1;
          i2 = 2;
        } else {
          m = 2;
          i2 = 1;
        }

        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }

        if (x4[m - 1] <= x4[i3 - 1]) {
          if (x4[i2 - 1] <= x4[i3 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i2;
            perm[2] = (signed char)i3;
            perm[3] = (signed char)i4;
          } else if (x4[i2 - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i2;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)i2;
          }
        } else if (x4[m - 1] <= x4[i4 - 1]) {
          if (x4[i2 - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i2;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)i2;
          }
        } else {
          perm[0] = (signed char)i3;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)m;
          perm[3] = (signed char)i2;
        }

        idx[ib] = idx4[perm[0] - 1];
        idx[ib + 1] = idx4[perm[1] - 1];
        idx[ib + 2] = idx4[perm[2] - 1];
        idx[ib + 3] = idx4[perm[3] - 1];
        x[ib] = x4[perm[0] - 1];
        x[ib + 1] = x4[perm[1] - 1];
        x[ib + 2] = x4[perm[2] - 1];
        x[ib + 3] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  if (ib > 0) {
    for (m = 0; m < 4; m++) {
      perm[m] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx[(k - nNaNs) - ib] = idx4[perm[k - 1] - 1];
      x[(k - nNaNs) - ib] = x4[perm[k - 1] - 1];
    }
  }

  m = (nNaNs + 122) >> 1;
  for (k = 1; k <= m; k++) {
    ib = idx[k - nNaNs];
    idx[k - nNaNs] = idx[123 - k];
    idx[123 - k] = ib;
    x[k - nNaNs] = xwork[123 - k];
    x[123 - k] = xwork[k - nNaNs];
  }

  if (((nNaNs + 122) & 1) != 0) {
    x[(m - nNaNs) + 1] = xwork[(m - nNaNs) + 1];
  }

  if (1 - nNaNs > 1) {
    memset(&iwork[0], 0, 123U * sizeof(int));
    b_merge_block(idx, x, 0, 1 - nNaNs, 2, iwork, xwork);
  }
}

//
// Arguments    : const emxArray_real_T *iPk
//                const emxArray_real_T *iInf
//                emxArray_real_T *iPkOut
// Return Type  : void
//
static void combinePeaks(const emxArray_real_T *iPk, const emxArray_real_T *iInf,
  emxArray_real_T *iPkOut)
{
  emxArray_int32_T *ia;
  emxArray_int32_T *ib;
  emxInit_int32_T(&ia, 1);
  emxInit_int32_T(&ib, 1);
  do_vectors(iPk, iInf, iPkOut, ia, ib);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&ia);
}

//
// FILT_CUSTOM Allows for quick customization of bandpass filter parameters
// Arguments    : const double X[250]
//                double Y[250]
// Return Type  : void
//
static void customFilt(const double X[250], double Y[250])
{
  double xtmp;
  double d0;
  int i;
  double y[274];
  double dv0[5];
  double dv1[5];
  double a[4];
  static const double dv2[5] = { 0.013356737905790865, 0.0,
    -0.026713475811581729, 0.0, 0.013356737905790865 };

  static const double dv3[5] = { 1.0, -3.6474831484713031, 4.9958898770882474,
    -3.0493284095486253, 0.70092168096543528 };

  double b_y[274];
  static const double b_a[4] = { -0.013356620059890656, -0.013357049900825781,
    0.013357014655895857, 0.013356655305044398 };

  double c_y[274];

  // vectorize
  //  cut off based on Fs
  xtmp = 2.0 * X[0];
  d0 = 2.0 * X[249];
  for (i = 0; i < 12; i++) {
    y[i] = xtmp - X[12 - i];
  }

  memcpy(&y[12], &X[0], 250U * sizeof(double));
  for (i = 0; i < 12; i++) {
    y[i + 262] = d0 - X[248 - i];
  }

  for (i = 0; i < 5; i++) {
    dv0[i] = dv2[i];
    dv1[i] = dv3[i];
  }

  for (i = 0; i < 4; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 274U * sizeof(double));
  filter(dv0, dv1, b_y, a, y);
  for (i = 0; i < 137; i++) {
    xtmp = y[i];
    y[i] = y[273 - i];
    y[273 - i] = xtmp;
  }

  for (i = 0; i < 5; i++) {
    dv0[i] = dv2[i];
    dv1[i] = dv3[i];
  }

  for (i = 0; i < 4; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&c_y[0], &y[0], 274U * sizeof(double));
  filter(dv0, dv1, c_y, a, y);
  for (i = 0; i < 137; i++) {
    xtmp = y[i];
    y[i] = y[273 - i];
    y[273 - i] = xtmp;
  }

  memcpy(&Y[0], &y[12], 250U * sizeof(double));
}

//
// CWT Continuous 1-D wavelet transform.
// Arguments    : const double SIG[250]
//                double coefs[5000]
// Return Type  : void
//
static void cwt_haar(const double SIG[250], double coefs[5000])
{
  int ind;
  emxArray_real_T *j;
  emxArray_real_T *x;
  emxArray_real_T *y;
  emxArray_real_T *work;
  emxArray_real_T *a;
  int k;
  double tmp2;
  int ndbl;
  int apnd;
  double cdiff;
  int iyLead;
  int b_k;
  int orderForDim;
  static const double dv7[1026] = { 0.0, 0.000976562500000001, 0.001953125,
    0.0029296875, 0.00390625, 0.0048828125, 0.00585937500000001,
    0.00683593750000001, 0.00781250000000001, 0.00878906250000001,
    0.00976562500000001, 0.0107421875, 0.01171875, 0.0126953125, 0.013671875,
    0.0146484375, 0.015625, 0.0166015625, 0.017578125, 0.0185546875, 0.01953125,
    0.0205078125, 0.021484375, 0.0224609375, 0.0234375, 0.0244140625,
    0.025390625, 0.0263671875, 0.02734375, 0.0283203125, 0.029296875,
    0.0302734375, 0.03125, 0.0322265625, 0.033203125, 0.0341796875, 0.03515625,
    0.0361328125, 0.037109375, 0.0380859375, 0.0390625, 0.0400390625,
    0.041015625, 0.0419921875, 0.04296875, 0.0439453125, 0.044921875,
    0.0458984375, 0.046875, 0.0478515625, 0.048828125, 0.0498046875, 0.05078125,
    0.0517578125, 0.052734375, 0.0537109375, 0.0546875, 0.0556640625,
    0.056640625, 0.0576171875, 0.05859375, 0.0595703125, 0.060546875,
    0.0615234375, 0.0625, 0.0634765625, 0.064453125, 0.0654296875, 0.06640625,
    0.0673828125, 0.068359375, 0.0693359375, 0.0703125, 0.0712890625,
    0.072265625, 0.0732421875, 0.07421875, 0.0751953125, 0.076171875,
    0.0771484375, 0.078125, 0.0791015625, 0.080078125, 0.0810546875, 0.08203125,
    0.0830078125, 0.083984375, 0.0849609375, 0.0859375, 0.0869140625,
    0.087890625, 0.0888671875, 0.08984375, 0.0908203125, 0.091796875,
    0.0927734375, 0.09375, 0.0947265625, 0.095703125, 0.0966796875, 0.09765625,
    0.0986328125, 0.099609375, 0.1005859375, 0.1015625, 0.1025390625,
    0.103515625, 0.1044921875, 0.10546875, 0.1064453125, 0.107421875,
    0.1083984375, 0.109375, 0.1103515625, 0.111328125, 0.1123046875, 0.11328125,
    0.1142578125, 0.115234375, 0.1162109375, 0.1171875, 0.1181640625,
    0.119140625, 0.1201171875, 0.12109375, 0.1220703125, 0.123046875,
    0.1240234375, 0.125, 0.1259765625, 0.126953125, 0.1279296875, 0.12890625,
    0.1298828125, 0.130859375, 0.1318359375, 0.1328125, 0.1337890625,
    0.134765625, 0.1357421875, 0.13671875, 0.1376953125, 0.138671875,
    0.1396484375, 0.140625, 0.1416015625, 0.142578125, 0.1435546875, 0.14453125,
    0.1455078125, 0.146484375, 0.1474609375, 0.1484375, 0.1494140625,
    0.150390625, 0.1513671875, 0.15234375, 0.1533203125, 0.154296875,
    0.1552734375, 0.15625, 0.1572265625, 0.158203125, 0.1591796875, 0.16015625,
    0.1611328125, 0.162109375, 0.1630859375, 0.1640625, 0.1650390625,
    0.166015625, 0.1669921875, 0.16796875, 0.1689453125, 0.169921875,
    0.1708984375, 0.171875, 0.1728515625, 0.173828125, 0.1748046875, 0.17578125,
    0.1767578125, 0.177734375, 0.1787109375, 0.1796875, 0.1806640625,
    0.181640625, 0.1826171875, 0.18359375, 0.1845703125, 0.185546875,
    0.1865234375, 0.1875, 0.1884765625, 0.189453125, 0.1904296875, 0.19140625,
    0.1923828125, 0.193359375, 0.1943359375, 0.1953125, 0.1962890625,
    0.197265625, 0.1982421875, 0.19921875, 0.2001953125, 0.201171875,
    0.2021484375, 0.203125, 0.2041015625, 0.205078125, 0.2060546875, 0.20703125,
    0.2080078125, 0.208984375, 0.2099609375, 0.2109375, 0.2119140625,
    0.212890625, 0.2138671875, 0.21484375, 0.2158203125, 0.216796875,
    0.2177734375, 0.21875, 0.2197265625, 0.220703125, 0.2216796875, 0.22265625,
    0.2236328125, 0.224609375, 0.2255859375, 0.2265625, 0.2275390625,
    0.228515625, 0.2294921875, 0.23046875, 0.2314453125, 0.232421875,
    0.2333984375, 0.234375, 0.2353515625, 0.236328125, 0.2373046875, 0.23828125,
    0.2392578125, 0.240234375, 0.2412109375, 0.2421875, 0.2431640625,
    0.244140625, 0.2451171875, 0.24609375, 0.2470703125, 0.248046875,
    0.2490234375, 0.25, 0.2509765625, 0.251953125, 0.2529296875, 0.25390625,
    0.2548828125, 0.255859375, 0.2568359375, 0.2578125, 0.2587890625,
    0.259765625, 0.2607421875, 0.26171875, 0.2626953125, 0.263671875,
    0.2646484375, 0.265625, 0.2666015625, 0.267578125, 0.2685546875, 0.26953125,
    0.2705078125, 0.271484375, 0.2724609375, 0.2734375, 0.2744140625,
    0.275390625, 0.2763671875, 0.27734375, 0.2783203125, 0.279296875,
    0.2802734375, 0.28125, 0.2822265625, 0.283203125, 0.2841796875, 0.28515625,
    0.2861328125, 0.287109375, 0.2880859375, 0.2890625, 0.2900390625,
    0.291015625, 0.2919921875, 0.29296875, 0.2939453125, 0.294921875,
    0.2958984375, 0.296875, 0.2978515625, 0.298828125, 0.2998046875, 0.30078125,
    0.3017578125, 0.302734375, 0.3037109375, 0.3046875, 0.3056640625,
    0.306640625, 0.3076171875, 0.30859375, 0.3095703125, 0.310546875,
    0.3115234375, 0.3125, 0.3134765625, 0.314453125, 0.3154296875, 0.31640625,
    0.3173828125, 0.318359375, 0.3193359375, 0.3203125, 0.3212890625,
    0.322265625, 0.3232421875, 0.32421875, 0.3251953125, 0.326171875,
    0.3271484375, 0.328125, 0.3291015625, 0.330078125, 0.3310546875, 0.33203125,
    0.3330078125, 0.333984375, 0.3349609375, 0.3359375, 0.3369140625,
    0.337890625, 0.3388671875, 0.33984375, 0.3408203125, 0.341796875,
    0.3427734375, 0.34375, 0.3447265625, 0.345703125, 0.3466796875, 0.34765625,
    0.3486328125, 0.349609375, 0.3505859375, 0.3515625, 0.3525390625,
    0.353515625, 0.3544921875, 0.35546875, 0.3564453125, 0.357421875,
    0.3583984375, 0.359375, 0.3603515625, 0.361328125, 0.3623046875, 0.36328125,
    0.3642578125, 0.365234375, 0.3662109375, 0.3671875, 0.3681640625,
    0.369140625, 0.3701171875, 0.37109375, 0.3720703125, 0.373046875,
    0.3740234375, 0.375, 0.3759765625, 0.376953125, 0.3779296875, 0.37890625,
    0.3798828125, 0.380859375, 0.3818359375, 0.3828125, 0.3837890625,
    0.384765625, 0.3857421875, 0.38671875, 0.3876953125, 0.388671875,
    0.3896484375, 0.390625, 0.3916015625, 0.392578125, 0.3935546875, 0.39453125,
    0.3955078125, 0.396484375, 0.3974609375, 0.3984375, 0.3994140625,
    0.400390625, 0.4013671875, 0.40234375, 0.4033203125, 0.404296875,
    0.4052734375, 0.40625, 0.4072265625, 0.408203125, 0.4091796875, 0.41015625,
    0.4111328125, 0.412109375, 0.4130859375, 0.4140625, 0.4150390625,
    0.416015625, 0.4169921875, 0.41796875, 0.4189453125, 0.419921875,
    0.4208984375, 0.421875, 0.4228515625, 0.423828125, 0.4248046875, 0.42578125,
    0.4267578125, 0.427734375, 0.4287109375, 0.4296875, 0.4306640625,
    0.431640625, 0.4326171875, 0.43359375, 0.4345703125, 0.435546875,
    0.4365234375, 0.4375, 0.4384765625, 0.439453125, 0.4404296875, 0.44140625,
    0.4423828125, 0.443359375, 0.4443359375, 0.4453125, 0.4462890625,
    0.447265625, 0.4482421875, 0.44921875, 0.4501953125, 0.451171875,
    0.4521484375, 0.453125, 0.4541015625, 0.455078125, 0.4560546875, 0.45703125,
    0.4580078125, 0.458984375, 0.4599609375, 0.4609375, 0.4619140625,
    0.462890625, 0.4638671875, 0.46484375, 0.4658203125, 0.466796875,
    0.4677734375, 0.46875, 0.4697265625, 0.470703125, 0.4716796875, 0.47265625,
    0.4736328125, 0.474609375, 0.4755859375, 0.4765625, 0.4775390625,
    0.478515625, 0.4794921875, 0.48046875, 0.4814453125, 0.482421875,
    0.4833984375, 0.484375, 0.4853515625, 0.486328125, 0.4873046875, 0.48828125,
    0.4892578125, 0.490234375, 0.4912109375, 0.4921875, 0.4931640625,
    0.494140625, 0.4951171875, 0.49609375, 0.4970703125, 0.498046875,
    0.4990234375, 0.5, 0.4990234375, 0.498046875, 0.4970703125, 0.49609375,
    0.4951171875, 0.494140625, 0.4931640625, 0.4921875, 0.4912109375,
    0.490234375, 0.4892578125, 0.48828125, 0.4873046875, 0.486328125,
    0.4853515625, 0.484375, 0.4833984375, 0.482421875, 0.4814453125, 0.48046875,
    0.4794921875, 0.478515625, 0.4775390625, 0.4765625, 0.4755859375,
    0.474609375, 0.4736328125, 0.47265625, 0.4716796875, 0.470703125,
    0.4697265625, 0.46875, 0.4677734375, 0.466796875, 0.4658203125, 0.46484375,
    0.4638671875, 0.462890625, 0.4619140625, 0.4609375, 0.4599609375,
    0.458984375, 0.4580078125, 0.45703125, 0.4560546875, 0.455078125,
    0.4541015625, 0.453125, 0.4521484375, 0.451171875, 0.4501953125, 0.44921875,
    0.4482421875, 0.447265625, 0.4462890625, 0.4453125, 0.4443359375,
    0.443359375, 0.4423828125, 0.44140625, 0.4404296875, 0.439453125,
    0.4384765625, 0.4375, 0.4365234375, 0.435546875, 0.4345703125, 0.43359375,
    0.4326171875, 0.431640625, 0.4306640625, 0.4296875, 0.4287109375,
    0.427734375, 0.4267578125, 0.42578125, 0.4248046875, 0.423828125,
    0.4228515625, 0.421875, 0.4208984375, 0.419921875, 0.4189453125, 0.41796875,
    0.4169921875, 0.416015625, 0.4150390625, 0.4140625, 0.4130859375,
    0.412109375, 0.4111328125, 0.41015625, 0.4091796875, 0.408203125,
    0.4072265625, 0.40625, 0.4052734375, 0.404296875, 0.4033203125, 0.40234375,
    0.4013671875, 0.400390625, 0.3994140625, 0.3984375, 0.3974609375,
    0.396484375, 0.3955078125, 0.39453125, 0.3935546875, 0.392578125,
    0.3916015625, 0.390625, 0.3896484375, 0.388671875, 0.3876953125, 0.38671875,
    0.3857421875, 0.384765625, 0.3837890625, 0.3828125, 0.3818359375,
    0.380859375, 0.3798828125, 0.37890625, 0.3779296875, 0.376953125,
    0.3759765625, 0.375, 0.3740234375, 0.373046875, 0.3720703125, 0.37109375,
    0.3701171875, 0.369140625, 0.3681640625, 0.3671875, 0.3662109375,
    0.365234375, 0.3642578125, 0.36328125, 0.3623046875, 0.361328125,
    0.3603515625, 0.359375, 0.3583984375, 0.357421875, 0.3564453125, 0.35546875,
    0.3544921875, 0.353515625, 0.3525390625, 0.3515625, 0.3505859375,
    0.349609375, 0.3486328125, 0.34765625, 0.3466796875, 0.345703125,
    0.3447265625, 0.34375, 0.3427734375, 0.341796875, 0.3408203125, 0.33984375,
    0.3388671875, 0.337890625, 0.3369140625, 0.3359375, 0.3349609375,
    0.333984375, 0.3330078125, 0.33203125, 0.3310546875, 0.330078125,
    0.3291015625, 0.328125, 0.3271484375, 0.326171875, 0.3251953125, 0.32421875,
    0.3232421875, 0.322265625, 0.3212890625, 0.3203125, 0.3193359375,
    0.318359375, 0.3173828125, 0.31640625, 0.3154296875, 0.314453125,
    0.3134765625, 0.3125, 0.3115234375, 0.310546875, 0.3095703125, 0.30859375,
    0.3076171875, 0.306640625, 0.3056640625, 0.3046875, 0.3037109375,
    0.302734375, 0.3017578125, 0.30078125, 0.2998046875, 0.298828125,
    0.2978515625, 0.296875, 0.2958984375, 0.294921875, 0.2939453125, 0.29296875,
    0.2919921875, 0.291015625, 0.2900390625, 0.2890625, 0.2880859375,
    0.287109375, 0.2861328125, 0.28515625, 0.2841796875, 0.283203125,
    0.2822265625, 0.28125, 0.2802734375, 0.279296875, 0.2783203125, 0.27734375,
    0.2763671875, 0.275390625, 0.2744140625, 0.2734375, 0.2724609375,
    0.271484375, 0.2705078125, 0.26953125, 0.2685546875, 0.267578125,
    0.2666015625, 0.265625, 0.2646484375, 0.263671875, 0.2626953125, 0.26171875,
    0.2607421875, 0.259765625, 0.2587890625, 0.2578125, 0.2568359375,
    0.255859375, 0.2548828125, 0.25390625, 0.2529296875, 0.251953125,
    0.2509765625, 0.25, 0.2490234375, 0.248046875, 0.2470703125, 0.24609375,
    0.2451171875, 0.244140625, 0.2431640625, 0.2421875, 0.2412109375,
    0.240234375, 0.2392578125, 0.23828125, 0.2373046875, 0.236328125,
    0.2353515625, 0.234375, 0.2333984375, 0.232421875, 0.2314453125, 0.23046875,
    0.2294921875, 0.228515625, 0.2275390625, 0.2265625, 0.2255859375,
    0.224609375, 0.2236328125, 0.22265625, 0.2216796875, 0.220703125,
    0.2197265625, 0.21875, 0.2177734375, 0.216796875, 0.2158203125, 0.21484375,
    0.2138671875, 0.212890625, 0.2119140625, 0.2109375, 0.2099609375,
    0.208984375, 0.2080078125, 0.20703125, 0.2060546875, 0.205078125,
    0.2041015625, 0.203125, 0.2021484375, 0.201171875, 0.2001953125, 0.19921875,
    0.1982421875, 0.197265625, 0.1962890625, 0.1953125, 0.1943359375,
    0.193359375, 0.1923828125, 0.19140625, 0.1904296875, 0.189453125,
    0.1884765625, 0.1875, 0.1865234375, 0.185546875, 0.1845703125, 0.18359375,
    0.1826171875, 0.181640625, 0.1806640625, 0.1796875, 0.1787109375,
    0.177734375, 0.1767578125, 0.17578125, 0.1748046875, 0.173828125,
    0.1728515625, 0.171875, 0.1708984375, 0.169921875, 0.1689453125, 0.16796875,
    0.1669921875, 0.166015625, 0.1650390625, 0.1640625, 0.1630859375,
    0.162109375, 0.1611328125, 0.16015625, 0.1591796875, 0.158203125,
    0.1572265625, 0.15625, 0.1552734375, 0.154296875, 0.1533203125, 0.15234375,
    0.1513671875, 0.150390625, 0.1494140625, 0.1484375, 0.1474609375,
    0.146484375, 0.1455078125, 0.14453125, 0.1435546875, 0.142578125,
    0.1416015625, 0.140625, 0.1396484375, 0.138671875, 0.1376953125, 0.13671875,
    0.1357421875, 0.134765625, 0.1337890625, 0.1328125, 0.1318359375,
    0.130859375, 0.1298828125, 0.12890625, 0.1279296875, 0.126953125,
    0.1259765625, 0.125, 0.1240234375, 0.123046875, 0.1220703125, 0.12109375,
    0.1201171875, 0.119140625, 0.1181640625, 0.1171875, 0.1162109375,
    0.115234375, 0.1142578125, 0.11328125, 0.1123046875, 0.111328125,
    0.1103515625, 0.109375, 0.1083984375, 0.107421875, 0.1064453125, 0.10546875,
    0.1044921875, 0.103515625, 0.1025390625, 0.1015625, 0.1005859375,
    0.0996093750000001, 0.0986328125000001, 0.0976562500000001,
    0.0966796875000001, 0.0957031250000001, 0.0947265625000001,
    0.0937500000000001, 0.0927734375000001, 0.0917968750000001,
    0.0908203125000001, 0.0898437500000001, 0.0888671875000001,
    0.0878906250000001, 0.0869140625000001, 0.0859375000000001,
    0.0849609375000001, 0.0839843750000001, 0.0830078125000001,
    0.0820312500000001, 0.0810546875000001, 0.0800781250000001,
    0.0791015625000001, 0.0781250000000001, 0.0771484375000001,
    0.0761718750000001, 0.0751953125000001, 0.0742187500000001,
    0.0732421875000001, 0.0722656250000001, 0.0712890625000001,
    0.0703125000000001, 0.0693359375000001, 0.0683593750000001,
    0.0673828125000001, 0.0664062500000001, 0.0654296875000001,
    0.0644531250000001, 0.0634765625000001, 0.0625000000000001,
    0.0615234375000001, 0.0605468750000001, 0.0595703125000001,
    0.0585937500000001, 0.0576171875000001, 0.0566406250000001,
    0.0556640625000001, 0.0546875000000001, 0.0537109375000001,
    0.0527343750000001, 0.0517578125000001, 0.0507812500000001,
    0.0498046875000001, 0.0488281250000001, 0.0478515625000001,
    0.0468750000000001, 0.0458984375000001, 0.0449218750000001,
    0.0439453125000001, 0.0429687500000001, 0.0419921875000001,
    0.0410156250000001, 0.0400390625000001, 0.0390625000000001,
    0.0380859375000001, 0.0371093750000001, 0.0361328125000001,
    0.0351562500000001, 0.0341796875000001, 0.0332031250000001,
    0.0322265625000001, 0.0312500000000001, 0.0302734375000001,
    0.0292968750000001, 0.0283203125000001, 0.0273437500000001,
    0.0263671875000001, 0.0253906250000001, 0.0244140625000001,
    0.0234375000000001, 0.0224609375000001, 0.0214843750000001,
    0.0205078125000001, 0.0195312500000001, 0.0185546875000001,
    0.0175781250000001, 0.0166015625000001, 0.0156250000000001,
    0.0146484375000001, 0.0136718750000001, 0.0126953125000001,
    0.0117187500000001, 0.0107421875000001, 0.00976562500000011,
    0.00878906250000011, 0.00781250000000011, 0.00683593750000011,
    0.00585937500000011, 0.00488281250000011, 0.00390625000000011,
    0.00292968750000011, 0.00195312500000011, 0.000976562500000105,
    1.04083408558608E-16, 1.04083408558608E-16 };

  int nc;
  boolean_T b2;
  double d;
  double last;

  //  xSIG    = (1:lenSIG);
  memset(&coefs[0], 0, 5000U * sizeof(double));
  ind = 0;

  //  [val_WAV,xWAV] = intwave(WAV,precis);
  //  wtype = 1;
  emxInit_real_T(&j, 2);
  emxInit_real_T(&x, 2);
  emxInit_real_T(&y, 2);
  emxInit_real_T1(&work, 1);
  emxInit_real_T(&a, 2);
  for (k = 0; k < 20; k++) {
    tmp2 = (1.0 + (double)k) * 1.0009765625;
    ndbl = (int)std::floor(tmp2 + 0.5);
    apnd = ndbl;
    cdiff = (double)ndbl - tmp2;
    if (fabs(cdiff) < 4.4408920985006262E-16 * tmp2) {
      ndbl++;
    } else if (cdiff > 0.0) {
      tmp2 = (double)ndbl - 1.0;
    } else {
      ndbl++;
      tmp2 = apnd;
    }

    iyLead = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = ndbl;
    emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
    j->data[0] = 0.0;
    if (ndbl > 1) {
      j->data[ndbl - 1] = tmp2;
      apnd = (ndbl - 1) / 2;
      for (b_k = 1; b_k < apnd; b_k++) {
        j->data[b_k] = b_k;
        j->data[(ndbl - b_k) - 1] = tmp2 - (double)b_k;
      }

      if (apnd << 1 == ndbl - 1) {
        j->data[apnd] = tmp2 / 2.0;
      } else {
        j->data[apnd] = apnd;
        j->data[apnd + 1] = tmp2 - (double)apnd;
      }
    }

    cdiff = (1.0 + (double)k) * 0.0009765625;
    iyLead = j->size[0] * j->size[1];
    j->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
    orderForDim = j->size[0];
    apnd = j->size[1];
    apnd *= orderForDim;
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      j->data[iyLead] /= cdiff;
    }

    b_floor(j);
    iyLead = j->size[0] * j->size[1];
    j->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
    orderForDim = j->size[0];
    apnd = j->size[1];
    apnd *= orderForDim;
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      j->data[iyLead]++;
    }

    if (j->size[1] == 1) {
      iyLead = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = 2;
      emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
      for (iyLead = 0; iyLead < 2; iyLead++) {
        j->data[iyLead] = 1.0;
      }
    }

    iyLead = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = j->size[1];
    emxEnsureCapacity((emxArray__common *)x, iyLead, (int)sizeof(double));
    apnd = j->size[0] * j->size[1];
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      x->data[iyLead] = dv7[(int)j->data[iyLead] - 1];
    }

    iyLead = j->size[0] * j->size[1];
    j->size[0] = 1;
    j->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
    apnd = x->size[0] * x->size[1];
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      j->data[iyLead] = x->data[iyLead];
    }

    apnd = x->size[1] >> 1;
    for (ndbl = 1; ndbl <= apnd; ndbl++) {
      iyLead = x->size[1] - ndbl;
      cdiff = j->data[ndbl - 1];
      j->data[ndbl - 1] = j->data[iyLead];
      j->data[iyLead] = cdiff;
    }

    //      coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
    // UNTITLED5 mod of wconv1
    apnd = j->size[1];
    iyLead = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = apnd;
    emxEnsureCapacity((emxArray__common *)x, iyLead, (int)sizeof(double));
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      x->data[x->size[0] * iyLead] = j->data[iyLead];
    }

    if (x->size[1] == 0) {
      nc = 250;
    } else {
      nc = x->size[1] + 249;
    }

    iyLead = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nc;
    emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
    for (iyLead = 0; iyLead < nc; iyLead++) {
      y->data[iyLead] = 0.0;
    }

    if ((x->size[1] == 0) || (nc == 0)) {
      b2 = true;
    } else {
      b2 = false;
    }

    if (!b2) {
      if (x->size[1] <= nc - 1) {
        apnd = x->size[1];
      } else {
        apnd = nc;
      }

      for (orderForDim = 0; orderForDim < apnd; orderForDim++) {
        if (orderForDim + 249 < nc - 1) {
          ndbl = 249;
        } else {
          ndbl = (nc - orderForDim) - 1;
        }

        for (b_k = 0; b_k <= ndbl; b_k++) {
          iyLead = orderForDim + b_k;
          if (iyLead >= 0) {
          } else {
            iyLead = 0;
          }

          y->data[iyLead] += x->data[orderForDim] * SIG[b_k];
        }
      }
    }

    if (y->size[1] == 0) {
      iyLead = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
    } else {
      if (y->size[1] - 1 <= 1) {
        orderForDim = y->size[1] - 1;
      } else {
        orderForDim = 1;
      }

      if (orderForDim < 1) {
        iyLead = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
      } else {
        apnd = y->size[1] - orderForDim;
        iyLead = work->size[0];
        work->size[0] = orderForDim;
        emxEnsureCapacity((emxArray__common *)work, iyLead, (int)sizeof(double));
        iyLead = j->size[0] * j->size[1];
        j->size[0] = 1;
        j->size[1] = apnd;
        emxEnsureCapacity((emxArray__common *)j, iyLead, (int)sizeof(double));
        if (!(j->size[1] == 0)) {
          apnd = 1;
          iyLead = 0;
          work->data[0] = y->data[0];
          if (orderForDim >= 2) {
            for (ndbl = 1; ndbl < orderForDim; ndbl++) {
              cdiff = y->data[apnd];
              for (b_k = 0; b_k + 1 <= ndbl; b_k++) {
                tmp2 = work->data[b_k];
                work->data[b_k] = cdiff;
                cdiff -= tmp2;
              }

              work->data[ndbl] = cdiff;
              apnd++;
            }
          }

          for (ndbl = orderForDim + 1; ndbl <= y->size[1]; ndbl++) {
            cdiff = y->data[apnd];
            for (b_k = 0; b_k + 1 <= orderForDim; b_k++) {
              tmp2 = work->data[b_k];
              work->data[b_k] = cdiff;
              cdiff -= tmp2;
            }

            apnd++;
            j->data[iyLead] = cdiff;
            iyLead++;
          }
        }
      }
    }

    // UNTITLED4 mod of wkeep1
    iyLead = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = j->size[1];
    emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
    apnd = j->size[0] * j->size[1];
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      y->data[iyLead] = j->data[iyLead];
    }

    b2 = (250 < j->size[1]);
    if (!b2) {
    } else {
      d = ((double)j->size[1] - 250.0) / 2.0;
      apnd = (int)std::floor(d);
      last = (double)j->size[1] - std::ceil(d);
      if (1.0 + (double)apnd > last) {
        iyLead = 1;
        ndbl = 0;
      } else {
        iyLead = apnd + 1;
        ndbl = (int)last;
      }

      apnd = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (ndbl - iyLead) + 1;
      emxEnsureCapacity((emxArray__common *)y, apnd, (int)sizeof(double));
      apnd = (ndbl - iyLead) + 1;
      for (ndbl = 0; ndbl < apnd; ndbl++) {
        y->data[y->size[0] * ndbl] = j->data[(iyLead + ndbl) - 1];
      }
    }

    cdiff = -std::sqrt(1.0 + (double)k);
    iyLead = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)a, iyLead, (int)sizeof(double));
    apnd = y->size[1];
    for (iyLead = 0; iyLead < apnd; iyLead++) {
      a->data[a->size[0] * iyLead] = cdiff * y->data[y->size[0] * iyLead];
    }

    for (iyLead = 0; iyLead < 250; iyLead++) {
      coefs[ind + 20 * iyLead] = a->data[iyLead];
    }

    ind++;
  }

  emxFree_real_T(&a);
  emxFree_real_T(&work);
  emxFree_real_T(&y);
  emxFree_real_T(&x);
  emxFree_real_T(&j);
}

//
// Arguments    : const double x[250]
//                double y[249]
// Return Type  : void
//
static void diff(const double x[250], double y[249])
{
  int ixLead;
  int iyLead;
  double work;
  int m;
  double tmp2;
  ixLead = 1;
  iyLead = 0;
  work = x[0];
  for (m = 0; m < 249; m++) {
    tmp2 = work;
    work = x[ixLead];
    tmp2 = x[ixLead] - tmp2;
    ixLead++;
    y[iyLead] = tmp2;
    iyLead++;
  }
}

//
// Arguments    : const emxArray_real_T *a
//                const emxArray_real_T *b
//                emxArray_real_T *c
//                emxArray_int32_T *ia
//                emxArray_int32_T *ib
// Return Type  : void
//
static void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib)
{
  int na;
  int nb;
  int ncmax;
  int ibfirst;
  int nc;
  int nia;
  int nib;
  int iafirst;
  int ialast;
  int iblast;
  int b_ialast;
  double ak;
  int b_iblast;
  double bk;
  boolean_T p;
  emxArray_int32_T *b_ia;
  emxArray_int32_T *b_ib;
  emxArray_real_T *b_c;
  na = a->size[0];
  nb = b->size[0];
  ncmax = a->size[0] + b->size[0];
  ibfirst = c->size[0];
  c->size[0] = ncmax;
  emxEnsureCapacity((emxArray__common *)c, ibfirst, (int)sizeof(double));
  ibfirst = ia->size[0];
  ia->size[0] = a->size[0];
  emxEnsureCapacity((emxArray__common *)ia, ibfirst, (int)sizeof(int));
  ibfirst = ib->size[0];
  ib->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)ib, ibfirst, (int)sizeof(int));
  nc = -1;
  nia = -1;
  nib = 0;
  iafirst = 1;
  ialast = 1;
  ibfirst = 0;
  iblast = 1;
  while ((ialast <= na) && (iblast <= nb)) {
    b_ialast = ialast;
    ak = skip_to_last_equal_value(&b_ialast, a);
    ialast = b_ialast;
    b_iblast = iblast;
    bk = skip_to_last_equal_value(&b_iblast, b);
    iblast = b_iblast;
    if ((fabs(bk - ak) < eps(bk / 2.0)) || (rtIsInf(ak) && rtIsInf(bk) &&
         ((ak > 0.0) == (bk > 0.0)))) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      nc++;
      c->data[nc] = ak;
      nia++;
      ia->data[nia] = iafirst;
      ialast = b_ialast + 1;
      iafirst = b_ialast + 1;
      iblast = b_iblast + 1;
      ibfirst = b_iblast;
    } else {
      if ((ak < bk) || rtIsNaN(bk)) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        nc++;
        nia++;
        c->data[nc] = ak;
        ia->data[nia] = iafirst;
        ialast = b_ialast + 1;
        iafirst = b_ialast + 1;
      } else {
        nc++;
        nib++;
        c->data[nc] = bk;
        ib->data[nib - 1] = ibfirst + 1;
        iblast = b_iblast + 1;
        ibfirst = b_iblast;
      }
    }
  }

  while (ialast <= na) {
    iafirst = ialast;
    ak = skip_to_last_equal_value(&iafirst, a);
    nc++;
    nia++;
    c->data[nc] = ak;
    ia->data[nia] = ialast;
    ialast = iafirst + 1;
  }

  while (iblast <= nb) {
    iafirst = iblast;
    bk = skip_to_last_equal_value(&iafirst, b);
    nc++;
    nib++;
    c->data[nc] = bk;
    ib->data[nib - 1] = iblast;
    iblast = iafirst + 1;
  }

  if (a->size[0] > 0) {
    if (1 > nia + 1) {
      iafirst = -1;
    } else {
      iafirst = nia;
    }

    emxInit_int32_T(&b_ia, 1);
    ibfirst = b_ia->size[0];
    b_ia->size[0] = iafirst + 1;
    emxEnsureCapacity((emxArray__common *)b_ia, ibfirst, (int)sizeof(int));
    for (ibfirst = 0; ibfirst <= iafirst; ibfirst++) {
      b_ia->data[ibfirst] = ia->data[ibfirst];
    }

    ibfirst = ia->size[0];
    ia->size[0] = b_ia->size[0];
    emxEnsureCapacity((emxArray__common *)ia, ibfirst, (int)sizeof(int));
    iafirst = b_ia->size[0];
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      ia->data[ibfirst] = b_ia->data[ibfirst];
    }

    emxFree_int32_T(&b_ia);
  }

  if (b->size[0] > 0) {
    if (1 > nib) {
      iafirst = 0;
    } else {
      iafirst = nib;
    }

    emxInit_int32_T(&b_ib, 1);
    ibfirst = b_ib->size[0];
    b_ib->size[0] = iafirst;
    emxEnsureCapacity((emxArray__common *)b_ib, ibfirst, (int)sizeof(int));
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      b_ib->data[ibfirst] = ib->data[ibfirst];
    }

    ibfirst = ib->size[0];
    ib->size[0] = b_ib->size[0];
    emxEnsureCapacity((emxArray__common *)ib, ibfirst, (int)sizeof(int));
    iafirst = b_ib->size[0];
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      ib->data[ibfirst] = b_ib->data[ibfirst];
    }

    emxFree_int32_T(&b_ib);
  }

  if (ncmax > 0) {
    if (1 > nc + 1) {
      iafirst = -1;
    } else {
      iafirst = nc;
    }

    emxInit_real_T1(&b_c, 1);
    ibfirst = b_c->size[0];
    b_c->size[0] = iafirst + 1;
    emxEnsureCapacity((emxArray__common *)b_c, ibfirst, (int)sizeof(double));
    for (ibfirst = 0; ibfirst <= iafirst; ibfirst++) {
      b_c->data[ibfirst] = c->data[ibfirst];
    }

    ibfirst = c->size[0];
    c->size[0] = b_c->size[0];
    emxEnsureCapacity((emxArray__common *)c, ibfirst, (int)sizeof(double));
    iafirst = b_c->size[0];
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      c->data[ibfirst] = b_c->data[ibfirst];
    }

    emxFree_real_T(&b_c);
  }
}

//
// Arguments    : emxArray__common *emxArray
//                int oldNumel
//                int elementSize
// Return Type  : void
//
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize)
{
  int newNumel;
  int i;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, (unsigned int)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_boolean_T **pEmxArray
// Return Type  : void
//
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
// Return Type  : void
//
static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
// Return Type  : void
//
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

//
// Arguments    : emxArray_boolean_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
  emxArray_boolean_T *emxArray;
  int i;
  *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : double x
// Return Type  : double
//
static double eps(double x)
{
  double r;
  double absxk;
  int exponent;
  absxk = fabs(x);
  if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
    if (absxk <= 2.2250738585072014E-308) {
      r = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      r = std::ldexp(1.0, exponent - 53);
    }
  } else {
    r = rtNaN;
  }

  return r;
}

//
// featureExtraction Summary of this function goes here if I ever feel like
// writing one up.
// Arguments    : const double samplesX[250]
//                emxArray_real_T *F
// Return Type  : void
//
static void featureExtractionEOG(const double samplesX[250], emxArray_real_T *F)
{
  double waveletCoef[5000];
  double Scalogram[5000];
  int ixstart;
  double E;
  emxArray_real_T *Max;
  emxArray_real_T *Imax;
  double dv4[249];
  double dv5[249];
  double dv6[249];
  int ix;
  emxArray_real_T *Min;
  emxArray_real_T *Imin;
  unsigned int uv0[2];
  emxArray_real_T *Amplitude;
  emxArray_real_T *Velocity;
  emxArray_real_T *b_Imax;
  double y;
  double xbar;
  double r;
  double b_y;
  boolean_T exitg2;
  boolean_T exitg1;
  emxArray_real_T *T_findpeaks_distX;
  int T_count_findpeaks;

  //  figure(20);hold on; plot(samplesX);
  //  waveletCoef = cwt(samplesX,1:20,'haar');
  cwt_haar(samplesX, waveletCoef);
  for (ixstart = 0; ixstart < 5000; ixstart++) {
    Scalogram[ixstart] = fabs(waveletCoef[ixstart] * waveletCoef[ixstart]);
  }

  E = Scalogram[0];
  for (ixstart = 0; ixstart < 4999; ixstart++) {
    E += Scalogram[ixstart + 1];
  }

  emxInit_real_T(&Max, 2);
  emxInit_real_T(&Imax, 2);

  // Energy Level.
  b_diff(samplesX, dv4);
  findpeaks(dv4, Max, Imax);
  b_diff(samplesX, dv5);
  for (ix = 0; ix < 249; ix++) {
    dv6[ix] = -dv5[ix];
  }

  emxInit_real_T(&Min, 2);
  emxInit_real_T(&Imin, 2);
  findpeaks(dv6, Min, Imin);
  ix = Max->size[0] * Max->size[1];
  Max->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)Max, ix, (int)sizeof(double));
  ixstart = Max->size[0];
  ix = Max->size[1];
  ixstart *= ix;
  for (ix = 0; ix < ixstart; ix++) {
    Max->data[ix] -= Min->data[ix];
  }

  for (ix = 0; ix < 2; ix++) {
    uv0[ix] = (unsigned int)Max->size[ix];
  }

  emxInit_real_T(&Amplitude, 2);
  ix = Amplitude->size[0] * Amplitude->size[1];
  Amplitude->size[0] = 1;
  Amplitude->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)Amplitude, ix, (int)sizeof(double));
  for (ixstart = 0; ixstart + 1 <= Max->size[1]; ixstart++) {
    Amplitude->data[ixstart] = fabs(Max->data[ixstart]);
  }

  emxInit_real_T(&Velocity, 2);
  ix = Velocity->size[0] * Velocity->size[1];
  Velocity->size[0] = 1;
  Velocity->size[1] = Amplitude->size[1];
  emxEnsureCapacity((emxArray__common *)Velocity, ix, (int)sizeof(double));
  ixstart = Amplitude->size[0] * Amplitude->size[1];
  for (ix = 0; ix < ixstart; ix++) {
    Velocity->data[ix] = Amplitude->data[ix];
  }

  emxInit_real_T(&b_Imax, 2);
  ix = b_Imax->size[0] * b_Imax->size[1];
  b_Imax->size[0] = 1;
  b_Imax->size[1] = Imax->size[1];
  emxEnsureCapacity((emxArray__common *)b_Imax, ix, (int)sizeof(double));
  ixstart = Imax->size[0] * Imax->size[1];
  for (ix = 0; ix < ixstart; ix++) {
    b_Imax->data[ix] = Imax->data[ix] - Imin->data[ix];
  }

  emxFree_real_T(&Imin);
  emxFree_real_T(&Imax);
  mrdivide(Velocity, b_Imax);
  y = samplesX[0];
  emxFree_real_T(&b_Imax);
  ix = 0;
  xbar = samplesX[0];
  for (ixstart = 0; ixstart < 249; ixstart++) {
    y += samplesX[ixstart + 1];
    ix++;
    xbar += samplesX[ix];
  }

  xbar /= 250.0;
  ix = 0;
  r = samplesX[0] - xbar;
  b_y = r * r;
  for (ixstart = 0; ixstart < 249; ixstart++) {
    ix++;
    r = samplesX[ix] - xbar;
    b_y += r * r;
  }

  b_y /= 249.0;
  ixstart = 1;
  xbar = samplesX[0];
  if (rtIsNaN(samplesX[0])) {
    ix = 2;
    exitg2 = false;
    while ((!exitg2) && (ix < 251)) {
      ixstart = ix;
      if (!rtIsNaN(samplesX[ix - 1])) {
        xbar = samplesX[ix - 1];
        exitg2 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 250) {
    while (ixstart + 1 < 251) {
      if (samplesX[ixstart] > xbar) {
        xbar = samplesX[ixstart];
      }

      ixstart++;
    }
  }

  ixstart = 1;
  r = samplesX[0];
  if (rtIsNaN(samplesX[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 251)) {
      ixstart = ix;
      if (!rtIsNaN(samplesX[ix - 1])) {
        r = samplesX[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 250) {
    while (ixstart + 1 < 251) {
      if (samplesX[ixstart] < r) {
        r = samplesX[ixstart];
      }

      ixstart++;
    }
  }

  b_findpeaks(samplesX, Min, Max);
  emxInit_real_T(&T_findpeaks_distX, 2);
  if (Min->size[1] == 0) {
    T_count_findpeaks = 0;
    ix = T_findpeaks_distX->size[0] * T_findpeaks_distX->size[1];
    T_findpeaks_distX->size[0] = 1;
    T_findpeaks_distX->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)T_findpeaks_distX, ix, (int)sizeof
                      (double));
    T_findpeaks_distX->data[0] = 0.0;
  } else {
    ix = Min->size[1];
    if (1 >= ix) {
      ix = 1;
    }

    if (Min->size[1] == 0) {
      ix = 0;
    }

    T_count_findpeaks = ix;
    ix = Min->size[1];
    if (1 >= ix) {
      ix = 1;
    }

    if (Min->size[1] == 0) {
      ix = 0;
    }

    if (ix > 1) {
      ix = T_findpeaks_distX->size[0] * T_findpeaks_distX->size[1];
      T_findpeaks_distX->size[0] = 1;
      T_findpeaks_distX->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)T_findpeaks_distX, ix, (int)sizeof
                        (double));
      T_findpeaks_distX->data[0] = Max->data[Max->size[1] - 1] - Max->data[0];

      // TODO: TAKE AVG, NOT MAX-MIN
    } else {
      ix = T_findpeaks_distX->size[0] * T_findpeaks_distX->size[1];
      T_findpeaks_distX->size[0] = 1;
      T_findpeaks_distX->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)T_findpeaks_distX, ix, (int)sizeof
                        (double));
      T_findpeaks_distX->data[0] = 0.0;
    }
  }

  emxFree_real_T(&Min);
  emxFree_real_T(&Max);

  // %? Threshold Stuff:?
  //  T_countmin_1 = sum(samplesX<-1*10^-5 & samplesX>(-1*10^-4),2);
  //  T_countmin_2 = WCountMin(samplesX,-1*10^-4);
  //  T_countmax = WCountMax(samplesX,8.5*10^-5);
  ix = F->size[0] * F->size[1];
  F->size[0] = 1;
  F->size[1] = (Velocity->size[1] + Amplitude->size[1]) + 9;
  emxEnsureCapacity((emxArray__common *)F, ix, (int)sizeof(double));
  F->data[0] = E;
  F->data[F->size[0]] = trapz(samplesX);
  ixstart = Velocity->size[1];
  for (ix = 0; ix < ixstart; ix++) {
    F->data[F->size[0] * (ix + 2)] = Velocity->data[Velocity->size[0] * ix];
  }

  ixstart = Amplitude->size[1];
  for (ix = 0; ix < ixstart; ix++) {
    F->data[F->size[0] * ((ix + Velocity->size[1]) + 2)] = Amplitude->
      data[Amplitude->size[0] * ix];
  }

  F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 2)] = 0.004 *
    y;
  F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 3)] = std::
    sqrt(b_y);
  F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 4)] = xbar;
  F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 5)] = r;
  F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 6)] = trapz
    (samplesX);
  F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 7)] =
    T_count_findpeaks;
  for (ix = 0; ix < 1; ix++) {
    F->data[F->size[0] * ((Velocity->size[1] + Amplitude->size[1]) + 8)] =
      T_findpeaks_distX->data[0];
  }

  emxFree_real_T(&T_findpeaks_distX);
  emxFree_real_T(&Velocity);
  emxFree_real_T(&Amplitude);

  //  hold off;
}

//
// Arguments    : double b[5]
//                double a[5]
//                const double x[274]
//                const double zi[4]
//                double y[274]
// Return Type  : void
//
static void filter(double b[5], double a[5], const double x[274], const double
                   zi[4], double y[274])
{
  double a1;
  int k;
  double dbuffer[5];
  int j;
  a1 = a[0];
  if ((!((!rtIsInf(a[0])) && (!rtIsNaN(a[0])))) || (a[0] == 0.0) || (!(a[0] !=
        1.0))) {
  } else {
    for (k = 0; k < 5; k++) {
      b[k] /= a1;
    }

    for (k = 0; k < 4; k++) {
      a[k + 1] /= a1;
    }

    a[0] = 1.0;
  }

  for (k = 0; k < 4; k++) {
    dbuffer[k + 1] = zi[k];
  }

  for (j = 0; j < 274; j++) {
    for (k = 0; k < 4; k++) {
      dbuffer[k] = dbuffer[k + 1];
    }

    dbuffer[4] = 0.0;
    for (k = 0; k < 5; k++) {
      dbuffer[k] += x[j] * b[k];
    }

    for (k = 0; k < 4; k++) {
      dbuffer[k + 1] -= dbuffer[0] * a[k + 1];
    }

    y[j] = dbuffer[0];
  }
}

//
// Arguments    : const double yTemp[249]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void findLocalMaxima(const double yTemp[249], emxArray_real_T *iPk,
  emxArray_real_T *iInflect)
{
  double b_yTemp[251];
  boolean_T yFinite[251];
  int ii;
  boolean_T x[250];
  emxArray_int32_T *b_ii;
  int idx;
  int i0;
  boolean_T exitg3;
  emxArray_int32_T *r0;
  boolean_T guard3 = false;
  emxArray_real_T *iTemp;
  emxArray_real_T *c_yTemp;
  emxArray_real_T *s;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r1;
  int nx;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *c_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  b_yTemp[0] = rtNaN;
  memcpy(&b_yTemp[1], &yTemp[0], 249U * sizeof(double));
  b_yTemp[250] = rtNaN;
  for (ii = 0; ii < 251; ii++) {
    yFinite[ii] = !rtIsNaN(b_yTemp[ii]);
  }

  for (ii = 0; ii < 250; ii++) {
    x[ii] = ((b_yTemp[ii] != b_yTemp[ii + 1]) && (yFinite[ii] || yFinite[ii + 1]));
  }

  emxInit_int32_T(&b_ii, 1);
  idx = 0;
  i0 = b_ii->size[0];
  b_ii->size[0] = 250;
  emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 251)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= 250) {
        exitg3 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3) {
      ii++;
    }
  }

  emxInit_int32_T(&r0, 1);
  i0 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
  i0 = r0->size[0];
  r0->size[0] = 1 + b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r0, i0, (int)sizeof(int));
  r0->data[0] = 1;
  ii = b_ii->size[0];
  for (i0 = 0; i0 < ii; i0++) {
    r0->data[i0 + 1] = b_ii->data[i0] + 1;
  }

  emxInit_real_T1(&iTemp, 1);
  i0 = iTemp->size[0];
  iTemp->size[0] = r0->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i0, (int)sizeof(double));
  ii = r0->size[0];
  for (i0 = 0; i0 < ii; i0++) {
    iTemp->data[i0] = 1.0 + (double)(r0->data[i0] - 1);
  }

  emxFree_int32_T(&r0);
  emxInit_real_T1(&c_yTemp, 1);
  i0 = c_yTemp->size[0];
  c_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)c_yTemp, i0, (int)sizeof(double));
  ii = iTemp->size[0];
  for (i0 = 0; i0 < ii; i0++) {
    c_yTemp->data[i0] = b_yTemp[(int)iTemp->data[i0] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_boolean_T(&b_x, 1);
  emxInit_real_T1(&r1, 1);
  c_diff(c_yTemp, s);
  b_sign(s);
  c_diff(s, r1);
  i0 = b_x->size[0];
  b_x->size[0] = r1->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i0, (int)sizeof(boolean_T));
  ii = r1->size[0];
  emxFree_real_T(&c_yTemp);
  for (i0 = 0; i0 < ii; i0++) {
    b_x->data[i0] = (r1->data[i0] < 0.0);
  }

  emxFree_real_T(&r1);
  nx = b_x->size[0];
  idx = 0;
  i0 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
  ii = 1;
  exitg2 = false;
  while ((!exitg2) && (ii <= nx)) {
    guard2 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg2 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i0 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
    }
  } else {
    i0 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i0, (int)sizeof(int));
  }

  if (1.0 > (double)s->size[0] - 1.0) {
    ii = 0;
  } else {
    ii = (int)((double)s->size[0] - 1.0);
  }

  if (2 > s->size[0]) {
    i0 = 0;
  } else {
    i0 = 1;
  }

  idx = b_x->size[0];
  b_x->size[0] = ii;
  emxEnsureCapacity((emxArray__common *)b_x, idx, (int)sizeof(boolean_T));
  for (idx = 0; idx < ii; idx++) {
    b_x->data[idx] = (s->data[idx] != s->data[i0 + idx]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i0 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i0, (int)sizeof(int));
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx)) {
    guard1 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      c_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i0 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i0, (int)sizeof(int));
    }
  } else {
    i0 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i0, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  i0 = iInflect->size[0];
  iInflect->size[0] = c_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i0, (int)sizeof(double));
  ii = c_ii->size[0];
  for (i0 = 0; i0 < ii; i0++) {
    iInflect->data[i0] = iTemp->data[c_ii->data[i0]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  i0 = iPk->size[0];
  iPk->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i0, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i0 = 0; i0 < ii; i0++) {
    iPk->data[i0] = iTemp->data[b_ii->data[i0]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : const double Yin[249]
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void findpeaks(const double Yin[249], emxArray_real_T *Ypk,
                      emxArray_real_T *Xpk)
{
  boolean_T x[249];
  int k;
  emxArray_int32_T *ii;
  int idx;
  int cdiff;
  boolean_T exitg1;
  emxArray_real_T *iInfite;
  boolean_T guard1 = false;
  double yTemp[249];
  emxArray_real_T *iPk;
  emxArray_real_T *iInflect;
  int ndbl;
  emxArray_real_T *base;
  double b_idx;
  int apnd;
  emxArray_real_T *y;
  emxArray_real_T *c_idx;
  emxArray_real_T *d_idx;
  for (k = 0; k < 249; k++) {
    x[k] = (rtIsInf(Yin[k]) && (Yin[k] > 0.0));
  }

  emxInit_int32_T(&ii, 1);
  idx = 0;
  k = ii->size[0];
  ii->size[0] = 249;
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = 1;
  exitg1 = false;
  while ((!exitg1) && (cdiff < 250)) {
    guard1 = false;
    if (x[cdiff - 1]) {
      idx++;
      ii->data[idx - 1] = cdiff;
      if (idx >= 249) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      cdiff++;
    }
  }

  emxInit_real_T1(&iInfite, 1);
  k = ii->size[0];
  if (1 > idx) {
    ii->size[0] = 0;
  } else {
    ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  k = iInfite->size[0];
  iInfite->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInfite, k, (int)sizeof(double));
  cdiff = ii->size[0];
  for (k = 0; k < cdiff; k++) {
    iInfite->data[k] = ii->data[k];
  }

  memcpy(&yTemp[0], &Yin[0], 249U * sizeof(double));
  k = ii->size[0];
  ii->size[0] = iInfite->size[0];
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = iInfite->size[0];
  for (k = 0; k < cdiff; k++) {
    ii->data[k] = (int)iInfite->data[k];
  }

  cdiff = ii->size[0];
  for (k = 0; k < cdiff; k++) {
    yTemp[ii->data[k] - 1] = rtNaN;
  }

  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&iInflect, 1);
  findLocalMaxima(yTemp, iPk, iInflect);
  if (!(iPk->size[0] == 0)) {
    cdiff = iPk->size[0] - 1;
    ndbl = 0;
    for (idx = 0; idx <= cdiff; idx++) {
      if (Yin[(int)iPk->data[idx] - 1] > rtMinusInf) {
        ndbl++;
      }
    }

    k = 0;
    for (idx = 0; idx <= cdiff; idx++) {
      if (Yin[(int)iPk->data[idx] - 1] > rtMinusInf) {
        iPk->data[k] = iPk->data[idx];
        k++;
      }
    }

    k = iPk->size[0];
    iPk->size[0] = ndbl;
    emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  }

  cdiff = iPk->size[0];
  emxInit_real_T1(&base, 1);
  k = base->size[0];
  base->size[0] = cdiff;
  emxEnsureCapacity((emxArray__common *)base, k, (int)sizeof(double));
  for (k = 0; k + 1 <= cdiff; k++) {
    if ((Yin[(int)(iPk->data[k] - 1.0) - 1] >= Yin[(int)(iPk->data[k] + 1.0) - 1])
        || rtIsNaN(Yin[(int)(iPk->data[k] + 1.0) - 1])) {
      b_idx = Yin[(int)(iPk->data[k] - 1.0) - 1];
    } else {
      b_idx = Yin[(int)(iPk->data[k] + 1.0) - 1];
    }

    base->data[k] = b_idx;
  }

  cdiff = iPk->size[0] - 1;
  ndbl = 0;
  for (idx = 0; idx <= cdiff; idx++) {
    if (Yin[(int)iPk->data[idx] - 1] - base->data[idx] >= 0.0) {
      ndbl++;
    }
  }

  k = 0;
  for (idx = 0; idx <= cdiff; idx++) {
    if (Yin[(int)iPk->data[idx] - 1] - base->data[idx] >= 0.0) {
      iPk->data[k] = iPk->data[idx];
      k++;
    }
  }

  k = iPk->size[0];
  iPk->size[0] = ndbl;
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  combinePeaks(iPk, iInfite, base);
  emxFree_real_T(&iInfite);
  if (base->size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)base->size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - base->size[0]) + 1;
    idx = base->size[0];
    if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)idx) {
      ndbl++;
      apnd = base->size[0];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }
  }

  emxInit_real_T(&y, 2);
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  if (ndbl > 0) {
    y->data[0] = 1.0;
    if (ndbl > 1) {
      y->data[ndbl - 1] = apnd;
      idx = (ndbl - 1) / 2;
      for (k = 1; k < idx; k++) {
        y->data[k] = 1.0 + (double)k;
        y->data[(ndbl - k) - 1] = apnd - k;
      }

      if (idx << 1 == ndbl - 1) {
        y->data[idx] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[idx] = 1.0 + (double)idx;
        y->data[idx + 1] = apnd - idx;
      }
    }
  }

  emxInit_real_T1(&c_idx, 1);
  k = c_idx->size[0];
  c_idx->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
  cdiff = y->size[1];
  for (k = 0; k < cdiff; k++) {
    c_idx->data[k] = y->data[y->size[0] * k];
  }

  emxFree_real_T(&y);
  if (c_idx->size[0] == 0) {
  } else {
    k = iInflect->size[0];
    iInflect->size[0] = c_idx->size[0];
    emxEnsureCapacity((emxArray__common *)iInflect, k, (int)sizeof(double));
    cdiff = c_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      iInflect->data[k] = Yin[(int)base->data[(int)c_idx->data[k] - 1] - 1];
    }

    emxInit_real_T1(&d_idx, 1);
    sort(iInflect, ii);
    k = d_idx->size[0];
    d_idx->size[0] = ii->size[0];
    emxEnsureCapacity((emxArray__common *)d_idx, k, (int)sizeof(double));
    cdiff = ii->size[0];
    for (k = 0; k < cdiff; k++) {
      d_idx->data[k] = c_idx->data[ii->data[k] - 1];
    }

    k = c_idx->size[0];
    c_idx->size[0] = d_idx->size[0];
    emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
    cdiff = d_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      c_idx->data[k] = d_idx->data[k];
    }

    emxFree_real_T(&d_idx);
  }

  emxFree_int32_T(&ii);
  emxFree_real_T(&iInflect);
  if (c_idx->size[0] > 1) {
    b_idx = c_idx->data[0];
    k = c_idx->size[0];
    c_idx->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
    c_idx->data[0] = b_idx;
  }

  k = iPk->size[0];
  iPk->size[0] = c_idx->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  cdiff = c_idx->size[0];
  for (k = 0; k < cdiff; k++) {
    iPk->data[k] = base->data[(int)c_idx->data[k] - 1];
  }

  emxFree_real_T(&base);
  emxFree_real_T(&c_idx);
  k = Ypk->size[0] * Ypk->size[1];
  Ypk->size[0] = 1;
  Ypk->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Ypk, k, (int)sizeof(double));
  cdiff = iPk->size[0];
  for (k = 0; k < cdiff; k++) {
    Ypk->data[Ypk->size[0] * k] = Yin[(int)iPk->data[k] - 1];
  }

  k = Xpk->size[0] * Xpk->size[1];
  Xpk->size[0] = 1;
  Xpk->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Xpk, k, (int)sizeof(double));
  cdiff = iPk->size[0];
  for (k = 0; k < cdiff; k++) {
    Xpk->data[Xpk->size[0] * k] = 1.0 + (double)((int)iPk->data[k] - 1);
  }

  emxFree_real_T(&iPk);
}

//
// Arguments    : int n
//                double x
// Return Type  : int
//
static int ixamax(int n, double x)
{
  int idxmax;
  double smax;
  int k;
  double s;
  if (n < 1) {
    idxmax = 0;
  } else {
    idxmax = 1;
    if (n > 1) {
      smax = fabs(x);
      for (k = 2; k <= n; k++) {
        s = fabs(x);
        if (s > smax) {
          idxmax = k;
          smax = s;
        }
      }
    }
  }

  return idxmax;
}

//
// function yfit = knnclassification(testsamplesX,samplesX, samplesY, Knn, type)
//  Classify using the Nearest neighbor algorithm
//  Inputs:
//   tX    - Train samples
//  tY    - Train labels
//    tsX (testsamplesX) - Test  samples to classify
//  Knn         - Number of nearest neighbors
//
//  Outputs
//  result - Predicted targets
// if nargin < 5
//     type = '2norm';
// end
// Arguments    : const emxArray_real_T *tsX
// Return Type  : double
//
static double knn(const emxArray_real_T *tsX)
{
  double yfit;
  emxArray_real_T *Uc;
  int high_i;
  int nb;
  static const signed char a[123] = { 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 5, 5, 5, 4,
    4, 4, 5, 5, 5, 3, 3, 3, 4, 4, 4, 5, 5, 5, 4, 5, 5, 4, 4, 4, 5, 5, 5, 5, 6, 6,
    6, 4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 5, 5, 3, 4, 4, 4, 5, 5, 5, 3, 3, 3, 4, 4, 4,
    4, 5, 5, 5, 6, 6, 6, 4, 4, 4, 5, 5, 5, 3, 3, 3, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4,
    4, 5, 5, 5, 3, 3, 3, 6, 6, 6, 4, 4, 4, 4, 5, 5, 5, 3, 3, 3, 3, 5, 5, 5, 5, 4,
    4, 4, 3, 3, 3 };

  static const signed char idx[123] = { 20, 21, 22, 54, 61, 62, 63, 80, 81, 82,
    97, 98, 99, 110, 111, 112, 113, 121, 122, 123, 1, 2, 3, 14, 15, 16, 23, 24,
    25, 29, 32, 33, 34, 42, 43, 44, 45, 49, 50, 51, 55, 56, 57, 64, 65, 66, 67,
    74, 75, 76, 87, 88, 89, 90, 91, 92, 93, 103, 104, 105, 106, 118, 119, 120, 4,
    5, 6, 11, 12, 13, 17, 18, 19, 26, 27, 28, 30, 31, 35, 36, 37, 38, 46, 47, 48,
    52, 53, 58, 59, 60, 68, 69, 70, 77, 78, 79, 83, 84, 85, 86, 94, 95, 96, 107,
    108, 109, 114, 115, 116, 117, 7, 8, 9, 10, 39, 40, 41, 71, 72, 73, 100, 101,
    102 };

  emxArray_real_T *b_tsX;
  int nbins;
  int exitg5;
  int i2;
  boolean_T eok;
  emxArray_real_T *r4;
  double x[4059];
  double y[4059];
  int k;
  static const double tX[4059] = { 1.95455771918173E-5, 2.2138063922655E-5,
    2.77379073614066E-5, 6.94930556263154E-5, 7.28306583407363E-5,
    1.39354843778427E-5, 1.51532915669433E-5, 1.65161641475842E-5,
    2.63196101502014E-5, 2.17182027946734E-5, 5.99334338195675E-5,
    7.7474407548572E-5, 7.45820909595681E-5, 1.35735660191415E-5,
    1.16302466552936E-5, 1.34161271548026E-5, 3.94954875295947E-5,
    4.08132336848847E-5, 4.58555087583478E-5, 1.04664117071671E-5,
    9.21609393821206E-6, 1.29000814507056E-5, 1.12231596773552E-5,
    1.81629517599112E-5, 1.50883422830509E-5, 3.67858185768607E-5,
    3.66010470884004E-5, 4.38094724181388E-5, 6.5994017043041E-5,
    7.14981052465951E-5, 6.85894956269728E-5, 1.22756371041366E-5,
    9.33420560980704E-6, 4.42077726896543E-5, 7.44143130932922E-5,
    9.33216911142712E-5, 3.76851233741965E-5, 6.87530101177864E-5,
    1.09091231669032E-5, 1.10354691340417E-5, 1.14495477039844E-5,
    9.64176542131541E-6, 6.85188930862288E-6, 6.37131511223554E-6,
    6.87363625323037E-6, 5.61201326003545E-5, 4.78779113183162E-5,
    6.35552473861503E-5, 9.83760983510786E-6, 1.15908475679972E-5,
    1.18271189333382E-5, 4.17027829542948E-5, 4.86869168308431E-5,
    1.75196241039786E-5, 1.28690874721099E-5, 1.29302549077642E-5,
    1.67813368616538E-5, 6.25120334869243E-5, 4.63326149259269E-5,
    5.07655145550723E-5, 2.18564357241656E-5, 1.78138010470851E-5,
    2.55291438987319E-5, 1.75871577807866E-5, 1.30070330649033E-5,
    1.75952409145428E-5, 1.31288534502293E-5, 5.82992684648713E-5,
    6.08324709388825E-5, 4.43651515272214E-5, 1.81056981376546E-5,
    2.03579845605977E-5, 1.98703897072346E-5, 1.46043849387616E-5,
    1.12877399894156E-5, 8.98349754920065E-6, 7.67547447505327E-5,
    7.33454572750661E-5, 6.09995794089302E-5, 1.34496229850892E-5,
    1.08884585426564E-5, 1.0149406512475E-5, 9.49964865539569E-5,
    0.000110350915515455, 8.00969055100699E-5, 9.85491395931919E-5,
    1.02277187762546E-5, 1.37977587867506E-5, 1.05347660219272E-5,
    2.5884148373481E-5, 1.24374522898032E-5, 1.46962208975451E-5,
    2.24221738528087E-5, 7.50505906857337E-5, 5.94216315668637E-5,
    6.94286413649564E-5, 9.69603258674445E-6, 1.45727375259705E-5,
    1.89778394652199E-5, 1.91137775376173E-5, 1.54235087679355E-5,
    1.70705719998979E-5, 1.1250545548608E-5, 4.70712302354864E-5,
    4.60248609762129E-5, 3.03080206716878E-5, 0.000121073017759477,
    0.000148506231213222, 0.000156628339994503, 1.35336776065926E-5,
    2.57454517213571E-5, 1.21274563261983E-5, 1.50743187145708E-5,
    9.7722428264649E-5, 0.000115975366209266, 0.000108400736408943,
    0.000131344970724635, 8.23139251686005E-6, 7.63040364195976E-6,
    1.31475804583592E-5, 2.32156043753123E-5, 1.81112605887229E-5,
    2.76143796573164E-5, -0.0352035228685751, -0.00930365463604461,
    -0.0292385085457461, 0.0660307903173565, 0.0375617667168115,
    -0.00765906567050556, 0.00800109960989326, 0.0355260044136835,
    0.0318430454740929, -0.00365955425768218, 0.0638847045615644,
    0.0544767898877195, 0.0083954637393396, -0.0173315899445264,
    0.00541975374741047, -0.0155303713068456, 0.0674922699305101,
    0.0370445681941044, 0.00651972213820795, -0.00534928056508519,
    -0.0130385936411617, -0.0198862352748375, 0.00966630534126422,
    -0.0321838611668394, 0.000515677551483767, 0.029834395182133,
    -0.00726554549647225, -0.00344988538368944, -0.108076102862316,
    0.0286493593922848, -0.0154250371835571, -0.0227951990495208,
    -0.00544498533953719, -0.0755900327136059, 0.109865369605685,
    0.105213642414016, 0.00599372596159434, -0.116446852161714,
    0.000986057623137879, 0.00770309361680654, -0.00765706136974427,
    -0.0216434128428951, -0.0103245752290101, -0.0133565940294456,
    0.0259366787786986, 0.0606434419607969, 0.00591634099133887,
    0.00895724863119264, 0.00353807221238961, 0.0195258901820972,
    -0.00224165432969467, -0.0238179872022754, -0.00999645840744452,
    -0.00594736841588412, 0.00962388826427952, 0.00717208371744264,
    -0.012399849816307, 0.105721043669701, 0.0417430839622213,
    0.0192455127261545, -0.0237441300898348, 0.0346200850435575,
    -0.0279461876978678, -0.0282654891234318, -0.0274183574350524,
    -0.0177820915680353, -0.0204776139496475, 0.0836870885242251,
    0.0584289266050684, -0.0024590446368482, -0.00696370452593977,
    -0.0231716664740854, -0.0306926671317371, -0.040047379332178,
    0.00711136489240286, 0.0166609335332634, 0.0998536930289322,
    0.0594300360868701, 0.00946809070215567, -0.0263154819199412,
    0.00464351607440533, -0.00991079627830764, 0.0756832932440384,
    0.0713122051202202, 0.0428128065746302, 0.0107045240766212,
    -0.00946471745502574, 0.036792639746054, -0.000301862238573381,
    -0.0689477460319714, -0.0148180503527602, 0.00424537779496658,
    0.0270386138413965, 0.0779657951746449, -0.0165232817140292,
    -0.00271602193298637, 0.030560349845006, -0.0270833237588473,
    0.0262613134935012, -0.0322499365277371, -0.0100147712870111,
    -0.0131226792881083, 0.00838703220818822, -0.035979428422847,
    -0.0117875039374557, -0.00592852774646735, 0.0891147462832506,
    0.080131080898961, 0.0426457872950767, 0.00650084083172434,
    0.0591746638648513, -0.0149154207375208, 0.0192017265718292,
    0.0893268138407176, 0.0795959383192876, 0.0600327371576447,
    0.0276276274255963, 0.00312932232937147, 0.0127738079191165,
    -0.00498581144915601, 0.0392825872865747, -0.0192095025913887,
    0.0441975384644062, -1.00999624352613E-7, -1.01125881042955E-7,
    -1.01763825185602E-7, 1.55210512352928E-6, 1.55151264674353E-6,
    -1.6005373858571E-8, 4.49386262636041E-7, 4.48992783767647E-7,
    -2.30650735780585E-8, -7.68082373417693E-8, -8.95947280760499E-7,
    -8.96024325095086E-7, 2.45610354186612E-7, -2.33992012331545E-8,
    -2.24701421981153E-8, -4.01798300933377E-8, -2.71277887319069E-7,
    -2.712315441375E-7, 9.93710163737843E-8, 2.24999701295114E-8,
    8.99884905846022E-8, 9.06241582416367E-8, 5.75806619233199E-7,
    -1.55149226949429E-8, -1.63792378567277E-8, -2.44656961433188E-7,
    1.68912464314577E-7, 1.72291080638811E-7, 3.92523359769513E-8,
    5.52163805256781E-7, 5.52283844630391E-7, 4.63046635317125E-8,
    4.6137199545859E-8, -1.33091001061121E-8, 6.20670668101414E-7,
    1.6038588067161E-7, 1.17025143182494E-7, -1.05270952729294E-9,
    -3.22276413799676E-7, -3.22344875096267E-7, 9.6310847761603E-8,
    -2.41176531742534E-8, -9.74554821712019E-8, -9.73414464938316E-8,
    -2.0093058631929E-8, -7.63071961378184E-7, -7.63586472063029E-7,
    1.6406678279032E-7, 5.03992047579605E-7, -2.76510658321468E-8,
    -2.33737246471059E-8, 1.10918679190009E-7, -7.09800057059272E-7,
    1.21098877247847E-7, -6.88940967967523E-7, -9.6140038512677E-8,
    -9.6137760986346E-8, 2.54835985743481E-7, 2.55259170459348E-7,
    2.55147432664741E-7, -2.01390791544531E-7, -2.00807369074802E-7,
    1.31884665165438E-8, 1.64285749121774E-7, 1.64481767588819E-7,
    -3.66354194783804E-8, -5.371901567773E-8, -1.92583559608072E-6,
    -1.9257869561872E-6, 1.22667204052934E-7, 5.19386254686778E-7,
    5.19953368337663E-7, 5.09074669671671E-7, 1.59106870700548E-7,
    -2.88995184153054E-8, -2.18293360210844E-8, -1.82834187780998E-7,
    1.40676732404824E-6, 1.40692818830142E-6, 7.51442637374967E-9,
    7.81045832039564E-9, 7.41800621048218E-8, -2.78602118536801E-7,
    -2.78572639475596E-7, -2.86529561573017E-7, 1.25633225554234E-6,
    -7.866061601903E-7, -3.01851344259221E-8, -4.16114052364331E-8,
    8.12360709071894E-8, -4.42625296751771E-8, -4.47412097121608E-8,
    -3.62265204027752E-9, 4.23319556325471E-7, 4.23732443649291E-7,
    4.23705267647141E-7, 1.79677530404563E-7, 5.1467983437202E-8,
    2.67271833537165E-8, 6.77755047042431E-7, 6.7719603131997E-7,
    -6.29119659220251E-8, -3.57707152112893E-8, 1.57413141953656E-7,
    -1.08633496169796E-7, -1.08432891508263E-7, 4.05562671666737E-7,
    4.12527071143414E-7, 4.15556586825205E-7, 6.86913306034149E-8,
    -5.13394398563146E-9, 1.97470103785523E-8, 1.97620361829363E-8,
    -4.92567652776819E-7, -4.9266929326384E-7, -4.93016994603801E-7,
    2.19178061757997E-6, -1.42723199659035E-7, -6.69942201217349E-9,
    -8.40421566934555E-9, 5.08260332974378E-8, 5.57600010964711E-8,
    5.49401961968767E-8, 7.97897032385645E-6, 7.98894460239342E-6,
    8.03934218966255E-6, 2.94899973470563E-5, 2.94787402881271E-5,
    2.72091355595707E-6, 1.30322016164452E-5, 1.30207907292618E-5,
    9.45668016700399E-7, 2.45786359493662E-6, 2.41905765805335E-5,
    2.41926567775673E-5, 2.48066457728478E-5, 2.26972251961599E-6,
    2.17960379321718E-6, 2.1295309949469E-6, 2.17022309855255E-5,
    2.1698523531E-5, 1.77874119309074E-5, 3.89249483240547E-6,
    1.43981584935364E-6, 1.44998653186619E-6, 6.90967943079839E-6,
    2.26517871346167E-6, 2.39136872708224E-6, 2.17744695675538E-5,
    2.06073206463784E-5, 2.1019511837935E-5, 9.42056063446832E-7,
    2.98168454838662E-5, 2.98233276100411E-5, 2.03740519539535E-6,
    2.0300367800178E-6, 2.36901981888796E-6, 1.86201200430424E-5,
    1.78028327545487E-5, 1.65005451887317E-5, 8.84276002926067E-8,
    1.22465037243877E-5, 1.22491052536582E-5, 1.31945861433396E-5,
    4.10000103962308E-7, 1.75419867908163E-6, 1.75214603688897E-6,
    2.02939892182483E-6, 3.28120943392619E-5, 3.28342182987103E-5,
    3.2321156209693E-5, 6.04790457095527E-6, 3.45638322901835E-6,
    2.92171558088824E-6, 1.65268831993113E-5, 1.77450014264818E-5,
    7.99252589835793E-6, 8.9562325835778E-6, 8.55646342762825E-6,
    8.55626072778479E-6, 1.50353231588654E-5, 1.50602910571015E-5,
    1.50536985272197E-5, 3.2222526647125E-6, 3.21291790519683E-6,
    2.84870876757345E-6, 9.52857344906287E-6, 9.53994252015149E-6,
    5.45867750227868E-6, 1.28925637626552E-6, 2.88875339412109E-5,
    2.8886804342808E-5, 2.61281144632749E-5, 1.19458838577959E-5,
    1.19589274717663E-5, 1.17087174024484E-5, 4.6140992503159E-6,
    3.7569373939897E-6, 2.81598434671988E-6, 2.66937914160258E-5,
    2.67285791569167E-5, 2.67316355777271E-5, 5.18495419788727E-7,
    5.38921624107299E-7, 1.33524111788679E-6, 1.72733313492816E-5,
    1.7271503647487E-5, 1.77648328175271E-5, 3.26646386441008E-5,
    8.6526677620933E-6, 4.55795529831424E-6, 6.20009938022853E-6,
    8.36731530344051E-6, 4.51477802686806E-6, 5.45842758488362E-6,
    8.25964665183275E-7, 2.11659778162735E-5, 2.11866221824646E-5,
    2.1185263382357E-5, 5.57000344254145E-6, 5.35267027746901E-6,
    4.16944060317978E-6, 1.21995908467638E-5, 1.21895285637595E-5,
    3.64889402347746E-6, 3.29090579943861E-6, 1.29078776401998E-5,
    3.25900488509388E-6, 3.2529867452479E-6, 3.81228911366733E-5,
    3.87775446874809E-5, 3.90623191615693E-5, 1.03036995905122E-6,
    3.90179742907991E-7, 2.07343608974799E-6, 2.07501379920831E-6,
    3.5464870999931E-5, 3.54721891149965E-5, 3.54972236114736E-5,
    3.94520511164395E-5, 6.99343678329271E-6, 1.07190752194776E-6,
    1.33627029142594E-6, 4.21856076368734E-6, 4.6280800910071E-6,
    4.56003628434076E-6, -0.000141225496628324, -3.74098953540018E-5,
    -0.00011747156197138, 0.000265064285802323, 0.000151425952274796,
    -3.02372996631618E-5, 3.21345985600373E-5, 0.000142566688492709,
    0.00012749124592818, -1.47964186184976E-5, 0.000256432589272023,
    0.000219095022062852, 3.47493038634794E-5, -6.96372064998863E-5,
    2.1630464548038E-5, -6.2617055120801E-5, 0.00027074611617798,
    0.000148826461981352, 2.67080338903016E-5, -2.10039730742304E-5,
    -5.24363129980531E-5, -8.00392984559393E-5, 3.86409605545423E-5,
    -0.000129327298605874, 1.77438878842144E-6, 0.000119793440472769,
    -2.86972244247709E-5, -1.3219804307503E-5, -0.000433669797943145,
    0.000115307137164189, -6.13177933629013E-5, -9.15681655082016E-5,
    -2.21528309664511E-5, -0.000303263529454968, 0.000440787169439438,
    0.000422459075458157, 2.46103494364384E-5, -0.000467085514999692,
    4.10239809078529E-6, 3.09851611570627E-5, -3.06608041946383E-5,
    -8.68571905604252E-5, -4.15076206206351E-5, -5.36507572547649E-5,
    0.000103909759549947, 0.000243552035161577, 2.41819706381677E-5,
    3.67351067478978E-5, 1.4108541808113E-5, 7.82815947902428E-5,
    -9.04050848226879E-6, -9.52149120101837E-5, -3.97818864305285E-5,
    -2.38022971417799E-5, 3.85780929104775E-5, 2.86082197113229E-5,
    -5.00002132145179E-5, 0.000424023751055945, 0.000167799597468872,
    7.79386766374599E-5, -9.52642557420943E-5, 0.000138774418760872,
    -0.00011223169449627, -0.000113565880481503, -0.000109956074297105,
    -7.14021635141645E-5, -8.2353599820528E-5, 0.000335801866907511,
    0.000234699022019814, -9.25644799152158E-6, -2.79180880181872E-5,
    -9.30166913618746E-5, -0.000122811633013632, -0.000160809710358544,
    2.83097345797412E-5, 6.64772657466364E-5, 0.000400646273923184,
    0.000238953866223159, 3.88464024810533E-5, -0.000105584331601715,
    1.86458993529163E-5, -4.00245289179074E-5, 0.000303337188197778,
    0.000286204515366841, 0.000171946296464036, 4.40063510568519E-5,
    -3.79508488752056E-5, 0.000147531748662561, -1.39098805149254E-6,
    -0.000276641959066147, -5.95547258071779E-5, 1.68895885986709E-5,
    0.000108054164956654, 0.000312955152938512, -6.58014046277959E-5,
    -1.00993269830171E-5, 0.000122313963477163, -0.000108685471383953,
    0.000105536690661422, -0.000129432494983051, -4.01827508406443E-5,
    -5.23594001689003E-5, 3.37622011365163E-5, -0.00014431316571608,
    -4.76213854237101E-5, -2.39289700982686E-5, 0.00035770012810607,
    0.000322073027357932, 0.00017226542562978, 2.60256603088119E-5,
    0.000237626709809067, -5.95434296560964E-5, 7.70482438737563E-5,
    0.000357995585326268, 0.000319321896589479, 0.000241263272370076,
    0.000111847955101988, 1.24245874336736E-5, 5.10562860758066E-5,
    -2.01311286871191E-5, 0.000157987785893848, -7.70586779177347E-5,
    0.00017750385376542, 7.60144610956892E-5, 6.80704316927643E-5,
    4.98640774496691E-5, 0.000235971581626938, 0.000191872061726604,
    4.52603144616317E-5, 5.96320798901886E-5, 5.0148123743089E-5,
    8.72617981042286E-5, 8.85933432497177E-5, 0.00024566387355891,
    0.000254702790692719, 0.000170586874493444, 5.6719026048213E-5,
    4.92333834886759E-5, 4.12190552766119E-5, 0.000178540930312285,
    0.000188667925344535, 0.000132884379807641, 6.27115993091703E-5,
    4.97560349144256E-5, 3.78767768813491E-5, 0.000108502914826101,
    0.000109571913175819, 8.1890599465529E-5, 0.000184731858590446,
    0.00018091290694942, 0.000103832068159398, 4.63849119500226E-5,
    0.000232803952560384, 0.000178481362299659, 7.73088960654462E-5,
    6.11053521318118E-5, 5.2501838518092E-5, 0.000178887779277615,
    0.000198396491768265, 0.000154589610656714, 4.93195412166872E-5,
    4.95584399861846E-5, 5.01863462064172E-5, 5.03871955938129E-5,
    4.45270277518637E-5, 4.88890112672372E-5, 3.21545349514287E-5,
    3.00630429830427E-5, 0.000189089065724447, 0.000182628264778521,
    0.000107722446471796, 6.14122604484729E-5, 5.45320616628577E-5,
    4.63050627496899E-5, 0.000169437921579654, 0.000134236516513508,
    5.57261735572326E-5, 7.44806359046622E-5, 7.41127260574926E-5,
    6.33779588820073E-5, 0.000202140453959063, 0.00022473170530037,
    0.000172112813030656, 7.06383110266233E-5, 7.74931862295742E-5,
    8.58430020997106E-5, 7.43897895970622E-5, 6.90750788616195E-5,
    8.0240751189685E-5, 4.14410699098856E-5, 0.000168977479748565,
    0.000185475442613708, 0.000143742622947479, 7.3073628391009E-5,
    7.33947954873543E-5, 6.17188549332339E-5, 7.93222292352668E-5,
    8.25937000982814E-5, 5.14408473513303E-5, 0.000242209485088572,
    0.000253630908256094, 0.000185967921943903, 4.95923851200823E-5,
    5.90815425300268E-5, 3.95488806748276E-5, 0.000194927724307709,
    0.00025260815880175, 0.000242279410589606, 0.000156444924807736,
    7.5086579014908E-5, 7.13464854605577E-5, 4.96712306642355E-5,
    8.24900120272171E-5, 9.20275400320398E-5, 8.59883359955877E-5,
    6.06318371479049E-5, 0.000241557217397187, 0.000254955258716394,
    0.000180287113509962, 6.00434432966795E-5, 6.83258824160898E-5,
    5.74147354376131E-5, 5.7203822637635E-5, 5.63184667202832E-5,
    4.93820325234296E-5, 8.81654110425664E-5, 0.000156717974058166,
    0.000138744810214401, 0.000104772200837977, 0.000307189705977877,
    0.000284058904211272, 0.000150347635305152, 5.30926576086303E-5,
    5.95393873802421E-5, 5.8893400367096E-5, 4.37972193443924E-5,
    0.000172361070014522, 0.000280219129683511, 0.000280798772686014,
    0.000220485981424967, 8.70383093094842E-5, 7.44937550777144E-5,
    4.76321879042319E-5, 4.39912645161917E-5, 5.4643989678002E-5,
    5.14196423691602E-5, -2.39489454204274E-6, 7.36583568308204E-5,
    -3.56507215568974E-5, 0.000554240051338923, 0.000551162005420885,
    0.000226181624503806, 0.000129982844353953, 0.000240143260595644,
    0.000366905094209753, 0.000234676363425271, 0.000520263884994931,
    0.000610426434346954, 0.000616727372093112, 5.45601795024205E-5,
    0.000124687182141101, 2.21041365422768E-5, 0.000475610376622003,
    0.000432349516323667, 0.000421384134198387, 0.000167619189765204,
    4.98406144517206E-5, 3.71655760851319E-6, 0.000252651553592802,
    2.59652424056452E-5, 0.00010651844372255, 0.000318017984666425,
    0.000258587962925421, 0.000351149105189178, -0.000223214159852041,
    0.000456096083801337, 0.000381320625949701, 5.57877481278592E-5,
    8.77430443806099E-5, -0.000144600408505331, 0.000596208578554291,
    0.000679621036299319, 0.000361071211499067, -0.000220484607638359,
    8.70590059929267E-5, 0.000126825494756761, 7.67208304263206E-5,
    -9.88631421174991E-6, 2.68867690277569E-5, -8.40701535206414E-6,
    0.000168473386353756, 0.000486183915624012, 0.000350730755232484,
    0.000463502171330657, 0.000144672028684931, 0.00018721879649942,
    8.39119960634705E-5, 0.00016738934712277, 0.000299835402744128,
    8.33605554781516E-5, 0.000198458662474087, 0.000187734784106259,
    8.27002305179906E-5, 0.000617296350307554, 0.000470070745753443,
    0.000480548676782352, 0.000139533751912849, 0.000391182642337027,
    0.000154217355848611, 3.81602321147757E-5, 6.95984594286615E-7,
    8.43199806119701E-5, 3.49270107128125E-5, 0.000520036544437636,
    0.000514406259058611, 0.000379434358061064, 0.000104786380244895,
    4.26929394771736E-5, 3.93030112729765E-5, -4.18308098080258E-5,
    0.000160330255666997, 0.00017191884137523, 0.0006275563161439,
    0.000604836708193692, 0.000513656597760176, -1.5480938000214E-5,
    0.000106979583493518, 2.34333008008773E-5, 0.000520379643936414,
    0.000603369517884235, 0.000612875455228354, 0.000607137578917906,
    0.000113538940097419, 0.000278217136900878, 0.000101814577468192,
    -3.55067375606934E-5, 0.000108841334407284, 0.000144673612956474,
    0.000218434742552691, 0.000557962340115319, 0.000296111193404765,
    0.000479282698835517, 0.000225390591668283, 1.53897789885223E-5,
    0.000260612003138694, -4.1994673586326E-5, 4.84771033843959E-5,
    7.89392242525019E-5, 0.000204328057790263, 0.000292445028579414,
    0.000327530407752818, 0.000311001094163421, 0.000714141925997585,
    0.000796440718109164, 0.000782206049337098, 0.000142489000984657,
    0.000390094160953437, 8.94594701814582E-5, 0.000141850781052946,
    0.000572297551073676, 0.000677032008174562, 0.000739444051107545,
    0.000748686064891574, 0.000146873003844044, 0.000152881082029349,
    5.2937280043119E-5, 0.000341144308790586, 2.13803168552318E-5,
    0.000292497514831249, -0.000250888863489043, -0.00017662965170405,
    -0.000251234571141188, -2.24628657474001E-5, -1.91250456527387E-5,
    -0.0001091348674126, -0.000119098727958809, 2.75089101204168E-5,
    -3.06535180872256E-5, -0.000111353042269965, -0.000101612796921451,
    -1.47176050962403E-5, -9.14812438787409E-5, -0.000177659530745659,
    -6.28796695071611E-5, -0.000207737486357057, 1.34878226815748E-5,
    -4.37556803385163E-5, -9.24578499467293E-5, -0.000116048781274673,
    -0.000138003273036173, -0.000195533208068051, -0.000110499619793533,
    -0.00028659310448371, -0.000216813582484836, -0.000134638218290803,
    -0.000219726517967765, -0.000117301435543516, -0.000505873268805151,
    -0.00013641705409746, -0.000211184304575878, -0.000217458460819918,
    -0.000182213843831372, -0.000435902802927929, 8.40656949079584E-5,
    0.000122631864747293, -0.000125873442480518, -0.000552407665658911,
    -0.000132726007300897, -9.29879287844369E-5, -0.000143158930769828,
    -0.000172867176726003, -0.00014332824239424, -0.000128433383882077,
    3.50528168121648E-6, -2.63762152137727E-5, -0.000161987819545819,
    -5.34932485684835E-5, -7.41542808805609E-5, -1.8888412683523E-5,
    -0.000122162122820466, -0.000356511906129131, -0.000224650898745098,
    -0.000134590125385079, -6.39832724982459E-5, -7.45377140397396E-5,
    -0.000193939412234306, 2.88255720799177E-5, -9.24796511058289E-5,
    -6.91414409761603E-5, -0.0002045484609348, 3.20406440249852E-5,
    -0.000245227333207782, -0.000229013444926617, -0.000218589993102209,
    -0.000223404482067867, -0.000155710189494924, 6.13254701104065E-5,
    1.47997612238928E-5, -0.000134903122855828, -0.000237801036517877,
    -0.000300117389848153, -0.000311903509925639, -0.000294379759656583,
    -0.000142513206618297, -7.56894498465166E-5, -3.26844025395545E-5,
    -5.41300905870505E-5, -0.000130790905304362, -0.000211077175842497,
    -0.000103388992583931, -0.000160045850104701, -0.000176475642799815,
    -9.30063710683167E-5, -8.35565491127343E-5, -9.6273005595465E-5,
    -0.000160897065084074, 3.21259108958536E-6, -0.000118283625592834,
    -0.000389980731569933, -0.00017534583968609, -0.000126315117511708,
    -8.11099900691066E-5, -0.000102301240879159, -0.00036119058335127,
    -0.000174928690177443, -4.698363091868E-6, -0.00021434969742761,
    2.39037580981916E-5, -0.00030373532691426, -0.000213090785688218,
    -0.000226831633507268, -9.03795677310664E-5, -0.00034510348493,
    -0.000301728396355693, -0.00013879240236041, -0.000107228587949264,
    -2.48777705540085E-5, -2.99819749078424E-5, -0.000102636120877175,
    7.39330138774942E-5, -0.000194725209662085, -2.60981277281477E-5,
    -0.000138327547407774, -0.000127814143533401, -6.46616208569612E-5,
    -5.35604690554005E-5, -0.0001091913008881, -8.25455386921619E-5,
    -0.000127714915246552, 5.95480423327915E-5, -0.000188021832047069,
    8.7788452231908E-5, -0.0352035228685751, -0.00930365463604461,
    -0.0292385085457461, 0.0660307903173565, 0.0375617667168115,
    -0.00765906567050556, 0.00800109960989326, 0.0355260044136835,
    0.0318430454740929, -0.00365955425768218, 0.0638847045615644,
    0.0544767898877195, 0.0083954637393396, -0.0173315899445264,
    0.00541975374741047, -0.0155303713068456, 0.0674922699305101,
    0.0370445681941044, 0.00651972213820795, -0.00534928056508519,
    -0.0130385936411617, -0.0198862352748375, 0.00966630534126422,
    -0.0321838611668394, 0.000515677551483767, 0.029834395182133,
    -0.00726554549647225, -0.00344988538368944, -0.108076102862316,
    0.0286493593922848, -0.0154250371835571, -0.0227951990495208,
    -0.00544498533953719, -0.0755900327136059, 0.109865369605685,
    0.105213642414016, 0.00599372596159434, -0.116446852161714,
    0.000986057623137879, 0.00770309361680654, -0.00765706136974427,
    -0.0216434128428951, -0.0103245752290101, -0.0133565940294456,
    0.0259366787786986, 0.0606434419607969, 0.00591634099133887,
    0.00895724863119264, 0.00353807221238961, 0.0195258901820972,
    -0.00224165432969467, -0.0238179872022754, -0.00999645840744452,
    -0.00594736841588412, 0.00962388826427952, 0.00717208371744264,
    -0.012399849816307, 0.105721043669701, 0.0417430839622213,
    0.0192455127261545, -0.0237441300898348, 0.0346200850435575,
    -0.0279461876978678, -0.0282654891234318, -0.0274183574350524,
    -0.0177820915680353, -0.0204776139496475, 0.0836870885242251,
    0.0584289266050684, -0.0024590446368482, -0.00696370452593977,
    -0.0231716664740854, -0.0306926671317371, -0.040047379332178,
    0.00711136489240286, 0.0166609335332634, 0.0998536930289322,
    0.0594300360868701, 0.00946809070215567, -0.0263154819199412,
    0.00464351607440533, -0.00991079627830764, 0.0756832932440384,
    0.0713122051202202, 0.0428128065746302, 0.0107045240766212,
    -0.00946471745502574, 0.036792639746054, -0.000301862238573381,
    -0.0689477460319714, -0.0148180503527602, 0.00424537779496658,
    0.0270386138413965, 0.0779657951746449, -0.0165232817140292,
    -0.00271602193298637, 0.030560349845006, -0.0270833237588473,
    0.0262613134935012, -0.0322499365277371, -0.0100147712870111,
    -0.0131226792881083, 0.00838703220818822, -0.035979428422847,
    -0.0117875039374557, -0.00592852774646735, 0.0891147462832506,
    0.080131080898961, 0.0426457872950767, 0.00650084083172434,
    0.0591746638648513, -0.0149154207375208, 0.0192017265718292,
    0.0893268138407176, 0.0795959383192876, 0.0600327371576447,
    0.0276276274255963, 0.00312932232937147, 0.0127738079191165,
    -0.00498581144915601, 0.0392825872865747, -0.0192095025913887,
    0.0441975384644062, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.65865165661436E-5,
    2.40153772425395E-5, 3.99060749628191E-5, 0.000107736223023765,
    0.000116148076932817, 1.89620256697964E-5, 9.85445743055183E-5,
    8.37249477808528E-5, 0.000282191559232815, 0.000267317427022238,
    0.000102619148854497, 0.000123439964980764, 0.000141387825664208,
    9.63997249137935E-5, 9.77863378100776E-5, 7.10756013802924E-5,
    6.3622534844662E-5, 7.26578060918617E-5, 8.29371649654145E-5,
    9.31665003428367E-5, 0.000114826474915541, 0.000115023803961822,
    5.41580307219418E-6, 9.68670429443028E-6, 1.30513495045229E-5,
    7.51430073843942E-5, 6.11420956374875E-5, 6.36858443227361E-5,
    6.48638284643969E-5, 9.37435428159264E-5, 7.65643916069822E-5,
    1.00277143003455E-5, 1.20257697910472E-5, 4.04168879822613E-5,
    0.000149650803131344, 0.000137460890298612, 5.4271134856283E-5,
    5.38153145639022E-5, 8.28935478269036E-5, 7.71311194596297E-5,
    7.30191975380606E-5, 1.50118099723519E-5, 7.93889283372096E-6,
    8.08147654707132E-6, 1.27291308134187E-5, 7.56932340358719E-5,
    6.83742698695443E-5, 8.21684027414851E-5, 1.78612846331638E-5,
    2.26206227809909E-5, 2.00425861358549E-5, 5.49135986442197E-5,
    5.56235165242545E-5, 2.65572555322072E-5, 2.46698723483739E-5,
    2.13703195332401E-5, 2.28334905286641E-5, 0.000101071555270702,
    6.29957413716281E-5, 6.25070530905667E-5, 0.000393811946326805,
    0.000372507967479346, 0.000370481759183214, 2.05195642747289E-5,
    2.00272225706231E-5, 1.90924302165525E-5, 2.33249312655726E-5,
    5.92176730110708E-5, 5.55581722771687E-5, 2.96192206577505E-5,
    7.46805183972186E-5, 6.61380655314053E-5, 6.44084583550461E-5,
    9.18605483895596E-6, 7.34126533999439E-6, 1.13725875648098E-5,
    5.94136642735199E-5, 5.96401309434914E-5, 5.88724448578801E-5,
    0.000115234022151584, 0.000118755657370993, 0.000129789497771458,
    5.95253193472135E-5, 7.6568979907338E-5, 7.03139630632838E-5,
    7.02848438639128E-5, 2.75759197171858E-5, 2.34459989657249E-5,
    1.50082054974626E-5, 2.5690252886799E-5, 3.90619219096541E-5,
    3.49493945804082E-5, 4.78804414655889E-5, 0.000104585116761163,
    7.89937336534583E-5, 0.000113633500276127, 5.86257379331023E-5,
    9.3721849495361E-5, 5.7737614091802E-5, 7.89460357466762E-5,
    8.29653830053082E-5, 2.7908642134631E-5, 3.41099009394615E-5,
    0.000466918142981949, 0.000464668030978896, 0.000451611237642554,
    0.000144045258623626, 0.000156755838314871, 0.000173482892529555,
    0.000109464939958059, 8.89774775710966E-5, 9.80697806617638E-5,
    0.000111888863954397, 0.000106086081966578, 0.00013273460233453,
    0.000136362588342255, 0.000134492057558653, 1.163815506755E-5,
    1.88393486456042E-5, 1.31981703279784E-5, 8.97566285039322E-5,
    0.000115444738298557, 9.02213904925726E-5, -0.0579477137801724,
    -0.0282006763942179, -0.064470075933931, 0.0904228023440486,
    0.0510905893475633, 0.0105570470358791, 0.0648539350816602,
    0.0387194595669229, -0.0237538398198801, 0.0195856548259627,
    0.109854172788872, 0.067007034652225, 0.0250644303994382, 0.0100059821849891,
    0.008084475637745, -0.0234129812614853, 0.0961882995634636,
    0.0671052381666175, 0.0384512622518695, -0.058349105318089,
    -0.0361178786477529, -0.0371492096066879, 0.00497427887533475,
    -0.0189991085337769, -0.0155089069584194, 0.0646202870637739,
    0.000923865687138569, -0.0159288057045114, -0.112783318504642,
    0.0620213645518598, 0.00919319129099153, 0.0383239487740188,
    0.0301181137780944, -0.0768369164209851, 0.175540672351677, 0.13475285793464,
    0.000189615531421491, -0.102980894095479, 0.0742957369995521,
    0.0557548568058752, 0.0178339166108144, -0.0305463793512674,
    0.0039686270311005, -0.0250593646815294, 0.0396384933221974,
    0.0853946219492303, 0.0245897898424916, 0.0158811507676209,
    0.0543299394273095, 0.0545804463510483, 0.0238340652959139,
    -0.0262862084543839, -0.0254733508597779, -0.0203283489551911,
    0.0579831557618695, 0.0297773906277193, -0.032771571037686,
    0.147560259093868, 0.0513767102085904, 0.0151468979408811,
    -0.0602677119935489, 0.0220590937381448, -0.0193773726407221,
    0.0398574757825142, 0.00480917239482275, -0.00112674137095888,
    -0.0113284291055293, 0.114033607127821, 0.083239326470875,
    -5.13060664654533E-5, 0.049525898133738, 0.00458399840436961,
    -0.0193766783835538, -0.013836775732207, 0.0115277313083099,
    0.0321819618904084, 0.094137245090524, 0.0551378347830873,
    0.0246931811889634, -0.0337801902289907, 0.00626495602354528,
    0.00130734674683958, 0.0858811046955239, 0.085038926205366,
    0.0631713056752069, 0.0128828024987608, 0.0475441428069118,
    0.0521399606889333, 0.0178639832184625, -0.0514352807682576,
    0.0651257966829733, 0.00512732025771713, 0.0321039960202167,
    0.100841807879346, -0.00830187994385378, 0.0358189393323283,
    -0.0211075913641565, -0.0673489889624698, 0.00602810875895397,
    0.0305453666616202, 0.0208885133162025, -0.0216505694774143,
    0.0593513580067027, -0.100988413246856, -0.00193360130178145,
    0.0167091075585676, 0.119695547666538, 0.0776299178082854,
    0.0444506535344717, -0.0677997055479778, -0.00188284389994988,
    -0.0370667179456813, -0.0286063240725531, 0.107096170300044,
    0.0977714522800726, 0.0969816982774425, 0.0347045126307449,
    0.0165858464934175, 0.0337496744510967, -0.00541739107378453,
    -0.00218126461476586, -0.049119257500605, 0.0605234330574666,
    -4.13957267697098E-7, -4.15596398209966E-7, -4.20722656411001E-7,
    -3.34806720207738E-7, 2.94649106218596E-7, -5.83795592728892E-8,
    1.38820301643151E-6, 1.38815265038675E-6, -9.0554712764586E-7,
    -9.13774272580367E-7, -8.21706909587072E-7, -8.21843872870863E-7,
    4.03670021920701E-7, -7.55227433291199E-9, -7.30560105514122E-9,
    4.98870303263095E-7, -1.70416626905042E-7, 2.4687438963775E-7,
    1.0356847580076E-7, 1.95303803193171E-7, -1.81522508839298E-6,
    -1.81566739151343E-6, -1.20509311747753E-7, -3.16697435488868E-8,
    -3.2694423198195E-8, -2.29271653717511E-7, 1.83859138442437E-7,
    1.31522923746956E-7, 4.36806227483454E-8, 4.92861050328755E-7,
    1.26500167875119E-7, -5.49481963283162E-7, -5.49448367011527E-7,
    -5.48592323893166E-7, 5.50438003075434E-7, 5.49734833577458E-7,
    1.09444196023703E-7, -2.32611201306385E-8, 2.39405873866014E-6,
    2.39407809264449E-6, 2.39510182370333E-6, 1.96381515989365E-8,
    9.31886026763582E-8, 9.29474996645121E-8, -5.3847822285792E-8,
    -2.90591725894723E-7, -6.15184659351018E-7, 1.36537107879493E-7,
    -1.77280436063455E-7, -1.78146092165772E-7, -1.7849940927284E-7,
    5.68623303402118E-8, 1.05466922928125E-7, 4.26906108699402E-8,
    -3.3218670167598E-7, -3.32105988311833E-7, -3.31245387978048E-7,
    2.17463240370438E-7, 1.3986703437631E-7, 1.39851811193439E-7,
    -2.35120597564933E-6, -2.35206263853503E-6, -2.30931207959508E-6,
    -8.38786572810205E-8, -6.90294330566946E-8, -6.87009330047466E-8,
    8.56636243409768E-8, -6.83197946028809E-7, -6.83666592412821E-7,
    3.87361137644219E-8, -8.63752346524672E-7, -8.64216511284355E-7,
    1.36888796134702E-6, -3.4778856101017E-9, -9.22086145099873E-8,
    1.69369602365574E-8, -1.00924648425595E-7, -2.68363745108524E-7,
    1.20403611576593E-7, -1.85801934606501E-6, -1.85821867761994E-6,
    -1.85925686982076E-6, -7.01678398523397E-8, -6.99237559090788E-8,
    -7.46899569951382E-8, 8.47985704046574E-7, 8.46503815046294E-8,
    -1.09257723352849E-6, -1.09210192438613E-6, -1.29868610698571E-6,
    -2.35951140620613E-7, -2.35457996490795E-7, -5.16904594983093E-8,
    4.33029113615919E-7, 4.33306499351096E-7, 4.33017888657383E-7,
    -1.4856162144771E-6, -1.48517097443598E-6, -1.4882382975745E-6,
    2.98714460774002E-7, 2.69792686671816E-6, -4.93521230487975E-8,
    -2.84318514321072E-7, 7.20815399279208E-7, -1.45156606277851E-6,
    -1.45130913985734E-6, -2.4130673240206E-7, 2.90797272727267E-7,
    2.90467974294587E-7, -1.61925554856831E-6, -1.62058998706286E-6,
    -1.61853773997857E-6, -1.60581914715195E-6, -1.90772864219908E-6,
    -1.90918854152223E-6, -1.90832684474747E-6, -1.88328132509896E-6,
    -1.72875427553223E-8, -1.69186787604419E-8, -1.62749634028254E-8,
    -2.35920056135846E-6, -2.35754586645657E-6, -2.34053835650411E-6,
    7.03727355085067E-6, 7.06513876956942E-6, 7.15228515898702E-6,
    2.9462991378281E-5, 3.12328052591712E-5, 4.96226253819558E-6,
    3.60932784272193E-5, 3.60919689100554E-5, 2.62608667017299E-5,
    2.64994539048306E-5, 2.38295003780251E-5, 2.3833472313255E-5,
    2.58348814029249E-5, 3.92718265311424E-7, 3.79891254867343E-7,
    2.84356072859964E-5, 2.01091619747949E-5, 2.09843231192087E-5,
    1.84351886925353E-5, 2.79284438566234E-5, 2.90436014142877E-5,
    2.90506782642148E-5, 6.14597489913542E-6, 4.62378255813747E-6,
    4.77338578693648E-6, 2.0634448834576E-5, 2.20630966130925E-5,
    2.1043667799513E-5, 5.02327161605972E-6, 2.66144967177528E-5,
    2.59325344143995E-5, 7.14326552268111E-6, 7.14282877114985E-6,
    7.13170021061116E-6, 1.54122640861122E-5, 1.53925753401688E-5,
    1.48844106592237E-5, 1.95393409097364E-6, 3.59108810799021E-5,
    3.59111713896673E-5, 3.59265273555499E-5, 3.06355164943409E-6,
    3.35478969634889E-6, 3.34610998792243E-6, 5.86941262915133E-6,
    2.47002967010514E-5, 2.58377556927428E-5, 2.70343473601397E-5,
    8.33218049498238E-6, 8.37286633179128E-6, 8.38947223582347E-6,
    8.35876256001113E-6, 1.50817699787219E-5, 4.5678953630836E-6,
    4.65061382346372E-6, 4.64948383636566E-6, 4.63743543169267E-6,
    1.21779414607445E-5, 1.20285649563627E-5, 1.20272557626358E-5,
    3.99705015860387E-5, 3.99850648550954E-5, 3.92583053531164E-5,
    4.36169017861307E-6, 3.58953051894812E-6, 3.57244851624682E-6,
    7.28140806898303E-6, 1.16143650824898E-5, 1.1622332071018E-5,
    8.13458389052859E-6, 3.10950844748882E-5, 3.11117944062368E-5,
    3.28533110723286E-5, 5.04293413464746E-7, 1.19871198862984E-6,
    1.71063298389229E-6, 1.52396219122648E-5, 1.6906915941837E-5,
    1.72177164554528E-5, 2.97283095370402E-5, 2.9731498841919E-5,
    2.97481099171322E-5, 4.49074175054974E-6, 4.47512037818104E-6,
    4.78015724768885E-6, 2.20476283052109E-5, 1.16817526476389E-5,
    1.20183495688134E-5, 1.20131211682474E-5, 2.07789777117714E-5,
    1.39211172966162E-5, 1.38920217929569E-5, 1.18371152251128E-5,
    2.07853974535641E-5, 2.07987119688526E-5, 2.07848586555544E-5,
    2.22842432171565E-5, 2.22775646165397E-5, 2.23235744636175E-5,
    3.58457352928803E-5, 4.04689030007724E-5, 2.91177525987905E-6,
    3.98045920049501E-6, 8.57770325142258E-5, 2.90313212555702E-5,
    2.90261827971468E-5, 2.89568078882471E-5, 3.28600918181812E-5,
    3.28228810952883E-5, 2.5908088777093E-5, 2.59294397930058E-5,
    2.58966038396571E-5, 2.56931063544312E-5, 3.43391155595834E-5,
    3.43653937474002E-5, 3.43498832054544E-5, 3.38990638517813E-5,
    8.81664680521439E-7, 8.62852616782535E-7, 8.30023133544097E-7,
    3.77472089817354E-5, 3.77207338633052E-5, 3.74486137040658E-5,
    -0.000232560077234368, -0.000113029895615021, -0.000258768533944763,
    0.000363046906170551, 0.000206012026936151, 4.30457039180172E-5,
    0.000260495313492996, 0.000155939772928913, -9.53239733853024E-5,
    7.81049692687825E-5, 0.000440649349104312, 0.000269438438984996,
    0.000101936819092025, 3.98957681948904E-5, 3.20170337550346E-5,
    -9.32707235578094E-5, 0.000385955822941758, 0.000269591586779118,
    0.00015504625556969, -0.000233615943719584, -0.000145642536906716,
    -0.000149851312392806, 1.9906880506949E-5, -7.62131791670484E-5,
    -6.22849840086286E-5, 0.000259428980105814, 4.41616775903821E-6,
    -6.29643576324839E-5, -0.000452436581084312, 0.000249311117466098,
    3.76455306114018E-5, 0.00015347058280971, 0.000120652578926463,
    -0.000308208092584093, 0.000704219097755593, 0.000541037575183327,
    1.56368426773048E-6, -0.000412874323852592, 0.00029828379259623,
    0.000224156583642728, 7.21467879842046E-5, -0.000122582143706381,
    1.59012426587908E-5, -0.000100377491454773, 0.000158940966271337,
    0.000342860213256119, 9.92547661761174E-5, 6.47204726664361E-5,
    0.000217492020896702, 0.000218783678382637, 9.58669572704591E-5,
    -0.000104905009875075, -0.00010153556253641, -8.18916474011239E-5,
    0.000232477712796103, 0.00011943872967618, -0.000131494548387976,
    0.000591869058690742, 0.000206527682341303, 6.16196985746546E-5,
    -0.000242169057505078, 8.7755313167605E-5, -7.87026744776149E-5,
    0.000159739308584845, 1.9511716768412E-5, -4.55475554481118E-6,
    -4.54304730862552E-5, 0.00045749774837913, 0.000334149164490331,
    3.0376796440694E-7, 0.000198991116776977, 1.88125381779997E-5,
    -7.67918448573341E-5, -5.55781017688256E-5, 4.63429257841899E-5,
    0.000129217175938141, 0.000377619340231936, 0.000221712494178805,
    0.000100038173416382, -0.000135482243315929, 2.47547646355924E-5,
    4.3734532546281E-6, 0.000344248104584098, 0.000341240015266654,
    0.000253710346182815, 5.26475350428955E-5, 0.000191024994796131,
    0.000209192744851552, 7.16867704570485E-5, -0.00020613825093951,
    0.000261087849663537, 2.0424993018914E-5, 0.000128569713374593,
    0.000404728006351656, -3.25600601578329E-5, 0.000144748974618105,
    -8.51687523103457E-5, -0.000270611565671975, 2.38130905820195E-5,
    0.000122783738237218, 8.43937783759415E-5, -8.57851769881657E-5,
    0.000238354281438978, -0.000403991642422091, -8.2235056156219E-6,
    6.64365046810079E-5, 0.000480377396712603, 0.000312107133741891,
    0.000179768550637428, -0.000271884631695234, -7.46100854939777E-6,
    -0.00014903762961559, -0.000115592499817648, 0.000429533212764769,
    0.000392403320699747, 0.000389641322247356, 0.000140397540492798,
    6.65850380488566E-5, 0.000135668287812833, -2.15221781606881E-5,
    -8.60756554370946E-6, -0.000197391362286875, 0.000242200978914661,
    6.11037277850739E-5, 6.49554451363679E-5, 5.34333527318442E-5,
    0.000329933155295321, 0.000273756893628328, 6.01017880362342E-5,
    0.000249735958214659, 0.000190637701777342, 0.000377552986843078,
    0.00037174531226661, 0.000318818221193495, 0.000348863177172354,
    0.000247399239341213, 0.000196323881437421, 0.000189358624852576,
    0.000149172690794938, 0.000240679318758153, 0.000259992284597238,
    0.00018917348903617, 0.000254393014481343, 0.000236365635344933,
    0.000161480343114932, 5.50106563493903E-5, 6.79896008505217E-5,
    5.96034954068336E-5, 0.0002808400142252, 0.000268800078226451,
    0.000164284110407211, 7.29377554450619E-5, 0.000296217927448807,
    0.000219358610155699, 6.15694326041223E-5, 6.05449436724411E-5,
    6.40154880860365E-5, 0.000253661844613221, 0.000280565830352986,
    0.000216502980087276, 5.08307340690407E-5, 0.000189290927758254,
    0.000194563008576889, 0.000154746524209869, 4.4617777233326E-5,
    4.10859817236284E-5, 3.90939548030956E-5, 4.15551238949216E-5,
    0.000260415167339421, 0.000261827264104801, 0.000151724404607328,
    7.72327205051584E-5, 0.000108286841178761, 0.000125573235779089,
    0.00025987203534711, 0.000185923539666226, 7.3209995663941E-5,
    5.36006518212027E-5, 9.55436865121046E-5, 9.67921387222835E-5,
    0.000244040827375793, 0.000278186387929574, 0.000211996822591379,
    0.000356115680217206, 0.000330192116545714, 0.000279915473839321,
    6.65557097898315E-5, 9.21556662341637E-5, 9.05904224692437E-5,
    8.88984718587565E-5, 0.00015466110472311, 0.000190676138984423,
    0.000165063368000157, 0.000234380732581008, 0.000217502916889186,
    0.000121102771767721, 5.3855065890289E-5, 4.08653031033384E-5,
    4.78314953720029E-5, 0.000236138647425961, 0.000265473969298711,
    0.000200762027638412, 0.000232254433323624, 0.000230620734176336,
    0.000143736493872591, 0.00015641580581242, 0.000234024462385258,
    0.000248619367826202, 0.000160302819688846, 6.06853335226165E-5,
    6.35944440453106E-5, 6.35674410589141E-5, 9.06791488860068E-5,
    0.00010785620776048, 0.000109732981663286, 0.000114053190445487,
    0.00032569420200185, 0.000341227244504486, 0.000245644596646796,
    0.000179416338938757, 0.000151382185946961, 0.000112868834139642,
    0.00028489933107802, 0.000249247860258521, 6.31704569909779E-5,
    7.73125437179093E-5, 0.000406753645184541, 0.000425211770783659,
    0.000413356612799559, 0.000342493419238452, 0.000325569704424527,
    0.000177877716004313, 0.000179543348349264, 0.00020090483912727,
    0.000180110057970852, 0.000109284028081566, 0.000207200141752574,
    0.000332813006647391, 0.000329112747868577, 0.00024029805841287,
    6.0732264280463E-5, 6.41664958574703E-5, 5.77278441239008E-5,
    0.000185608639957113, 0.000171335202804113, 0.000119838863622781,
    -7.64631218328996E-5, 1.16667242796874E-5, -9.52501813846491E-5,
    0.000768394512781531, 0.000787769459695844, 0.000412269304747026,
    0.000595228610206324, 0.000552675478865765, 0.00102074871941566,
    0.00115955301963624, 0.000808065195006071, 0.000809596608480886,
    0.000877588454235295, 0.000477088677378848, 0.000453459600831625,
    0.000330778052695393, 0.000623478863047252, 0.000635237824710583,
    0.000665751564751684, 0.000324869255372279, 0.000332383449587273,
    0.000221139039950929, 0.00015318519246329, 4.76419358822869E-5,
    3.25723921285128E-5, 0.000566565946409423, 0.000464567909773758,
    0.000471950788529319, -0.000260800933550311, 0.000651159052080175,
    0.000540640693124184, 0.000305730979183301, 0.000272120727051274,
    -0.000149956907177045, 0.000939514413762414, 0.000900548956442043,
    0.000431234260258909, -0.000205144274429848, 0.00052878761700521,
    0.000515864155673889, 0.000456424591483054, -3.82493090787235E-5,
    0.000107549426914855, -1.92624964627946E-5, 0.000225926902973652,
    0.000630326760242222, 0.000525196491228356, 0.000618774659144668,
    0.00039594128069423, 0.000424694993381098, 0.000347910831435995,
    0.000252112587554352, 0.00036616449490713, 0.000100113035941088,
    0.000362221776106812, 0.000302980481327333, 5.71441671936467E-5,
    0.000805984709398233, 0.000566214592783747, 0.000534514280315835,
    0.000941541173470191, 0.00121433159128565, 0.000972739026081486,
    0.000369851667967869, 0.000261483921243951, 0.000241947227936135,
    0.000219510316517322, 0.000622134291465036, 0.000612709794428746,
    0.00038628421576513, 0.000463802898573985, 0.000379522876823182,
    0.000398099505848575, 6.15559805253535E-5, 0.000134680457189252,
    0.00023424740655026, 0.000603497554593567, 0.00059538467953879,
    0.000635822720371434, 0.000465105206325608, 0.000543619175478645,
    0.000418542566744028, 0.000556146310838728, 0.000658464679564162,
    0.000702094801608685, 0.000576766980999569, 0.000377419842835331,
    0.00036081165129573, 0.000229949224447046, 1.60252428247534E-5,
    0.000490711834054651, 0.000243682103608662, 0.000330192140313204,
    0.000727734862178332, 0.000493767896265898, 0.00075695854996452,
    0.000379955117691013, 0.000143747699265518, 0.000375337253704479,
    0.000444280449190752, 0.000546677476959082, 0.00044424970424964,
    0.000448164000290906, 0.00106204917114106, 0.00137876688088607,
    0.00141816302347269, 0.000879053408469296, 0.000866227704593365,
    0.000903809681214986, 0.000262895148017319, 0.000454903952712969,
    0.000235027599221024, 0.000187513581231546, 0.000656080808292376,
    0.00072430288435871, 0.000889065544015517, 0.000813890398742032,
    0.000171935492434203, 0.000316972582967517, 0.000121088160695309,
    0.000449523946866914, 0.000223171649162657, 0.000606020809000518,
    -0.00033047755645676, -0.000242779080810879, -0.000360563613541172,
    -4.69896493834855E-5, -3.50992784421387E-5, -3.01498875418508E-5,
    -6.04013129074278E-5, -2.91868972466963E-5, -0.000358608948338942,
    -0.00019545643075132, -8.53140295066541E-5, -8.09625095379215E-5,
    -9.188823091221E-5, -0.000188402357843168, -0.000212397910810509,
    -0.000324116894518528, -5.83932187897801E-6, -8.00325747667781E-6,
    -7.6745455424875E-6, -0.000550445960812878, -0.000576930426077899,
    -0.000540264483168983, -5.05542887778946E-5, -0.000176820789539102,
    -0.000212453415879539, -0.000154986021793856, -0.000251690401627693,
    -0.000178605732561326, -0.000557089512429277, -8.21871936027327E-5,
    -0.000143864274979405, 2.5542968006676E-5, 4.76715545547847E-7,
    -0.00042219906523264, 0.000127229122861326, 0.000131991974544707,
    -0.000200078277571573, -0.000507725240896921, -7.60252586420383E-5,
    -4.30892259972089E-5, -0.000102730831337383, -0.000212864830529643,
    -6.8045366370175E-5, -0.000171919466109965, 7.58557936666707E-5,
    1.27504522250187E-6, -0.000165381306765558, -5.60323329307168E-5,
    4.17425419052848E-5, 2.72404444442008E-5, -7.04702885276227E-5,
    -0.000418014409338704, -0.000301536533182936, -0.000271872606693636,
    5.47543247290441E-5, -7.89951118420805E-5, -0.0003260862579331,
    6.97501114091916E-5, -0.000140420474464971, -0.000141858076965093,
    -0.000586348135700866, -0.000287954145649329, -0.000502230919776719,
    1.04533927302106E-5, -0.000134213476048181, -0.000154379822218276,
    -0.000185212343005173, 8.13958729618893E-5, 4.70970745597403E-5,
    -0.000142221283436232, -0.000207738905672214, -0.000281259640687181,
    -0.000276171629069309, -0.000166039421017014, -3.34051384852718E-5,
    1.04167817036224E-5, -9.1700556828765E-5, -8.86792259007489E-5,
    -7.62394665214838E-5, -0.000392506742646888, -0.000314054966700326,
    -0.000421642766159432, -5.54482572892495E-5, 6.32428618725659E-6,
    4.81516875367441E-6, -6.47497079659003E-5, 4.67919414064225E-5,
    7.62503334909591E-5, -1.27678551512641E-5, -0.0002985678905096,
    4.44858984269253E-5, -0.000144511578228341, -5.83937692071984E-5,
    -0.000132942028066668, -0.000395257250182771, -4.29453280854291E-5,
    -0.000329334225692915, -0.000485390084306526, -0.000164238453444466,
    -0.000202745252838817, -0.000100158142698315, -0.000127210417265258,
    2.62607057923992E-5, -0.000695728047240271, -0.000377121104693356,
    -0.000337582897362302, 2.62569792575038E-5, 2.06171407108536E-5,
    5.5010035569513E-5, -0.000454878338675177, -0.00026323974878008,
    -0.000482688898199322, -0.000549771924075491, -0.000138806352834896,
    -0.000151925283684641, 3.32945690017148E-5, -4.08675178940395E-5,
    -6.21632453337867E-5, 1.78224212555643E-5, -9.70807536689145E-5,
    -0.000245093622587588, -0.000471154635938218, -7.97980352120904E-5,
    -0.0579477137801724, -0.0282006763942179, -0.064470075933931,
    0.0904228023440486, 0.0510905893475633, 0.0105570470358791,
    0.0648539350816602, 0.0387194595669229, -0.0237538398198801,
    0.0195856548259627, 0.109854172788872, 0.067007034652225, 0.0250644303994382,
    0.0100059821849891, 0.008084475637745, -0.0234129812614853,
    0.0961882995634636, 0.0671052381666175, 0.0384512622518695,
    -0.058349105318089, -0.0361178786477529, -0.0371492096066879,
    0.00497427887533475, -0.0189991085337769, -0.0155089069584194,
    0.0646202870637739, 0.000923865687138569, -0.0159288057045114,
    -0.112783318504642, 0.0620213645518598, 0.00919319129099153,
    0.0383239487740188, 0.0301181137780944, -0.0768369164209851,
    0.175540672351677, 0.13475285793464, 0.000189615531421491,
    -0.102980894095479, 0.0742957369995521, 0.0557548568058752,
    0.0178339166108144, -0.0305463793512674, 0.0039686270311005,
    -0.0250593646815294, 0.0396384933221974, 0.0853946219492303,
    0.0245897898424916, 0.0158811507676209, 0.0543299394273095,
    0.0545804463510483, 0.0238340652959139, -0.0262862084543839,
    -0.0254733508597779, -0.0203283489551911, 0.0579831557618695,
    0.0297773906277193, -0.032771571037686, 0.147560259093868,
    0.0513767102085904, 0.0151468979408811, -0.0602677119935489,
    0.0220590937381448, -0.0193773726407221, 0.0398574757825142,
    0.00480917239482275, -0.00112674137095888, -0.0113284291055293,
    0.114033607127821, 0.083239326470875, -5.13060664654533E-5,
    0.049525898133738, 0.00458399840436961, -0.0193766783835538,
    -0.013836775732207, 0.0115277313083099, 0.0321819618904084,
    0.094137245090524, 0.0551378347830873, 0.0246931811889634,
    -0.0337801902289907, 0.00626495602354528, 0.00130734674683958,
    0.0858811046955239, 0.085038926205366, 0.0631713056752069,
    0.0128828024987608, 0.0475441428069118, 0.0521399606889333,
    0.0178639832184625, -0.0514352807682576, 0.0651257966829733,
    0.00512732025771713, 0.0321039960202167, 0.100841807879346,
    -0.00830187994385378, 0.0358189393323283, -0.0211075913641565,
    -0.0673489889624698, 0.00602810875895397, 0.0305453666616202,
    0.0208885133162025, -0.0216505694774143, 0.0593513580067027,
    -0.100988413246856, -0.00193360130178145, 0.0167091075585676,
    0.119695547666538, 0.0776299178082854, 0.0444506535344717,
    -0.0677997055479778, -0.00188284389994988, -0.0370667179456813,
    -0.0286063240725531, 0.107096170300044, 0.0977714522800726,
    0.0969816982774425, 0.0347045126307449, 0.0165858464934175,
    0.0337496744510967, -0.00541739107378453, -0.00218126461476586,
    -0.049119257500605, 0.0605234330574666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.000118830099578963, 9.42867497061019E-5, 0.000108431431203111,
    0.000319146600094842, 0.000315285265719292, 2.28718105475272E-5,
    3.02963783202455E-5, 2.16002107669932E-5, 5.15147731110411E-5,
    3.32053726713654E-5, 0.000252590742845638, 0.000260223267181229,
    0.000258783230661902, 0.000215666587209383, 0.000241500459005499,
    0.000245014747298, 0.000240360607688649, 0.000271013581614261,
    0.000247453949376365, 2.94234848098821E-5, 2.78502112998125E-5,
    5.88578765203214E-6, 0.00018493940600373, 0.000155068609834147,
    0.000182543275547097, 0.000251689513606677, 0.000233521938626456,
    0.000237898794080971, 4.11674635051763E-5, 0.000195024884919976,
    0.000181811607717003, 0.000136733688857354, 0.000152522681082881,
    0.000115299022539122, 0.000201379997220528, 0.000176270659156836,
    0.000188699034957903, 3.24184696660712E-5, 1.88357674499489E-5,
    2.91240773162379E-5, 1.9241869172531E-5, 0.000154781090195011,
    0.000165247601511466, 0.000158085512693072, 0.000184580118113768,
    0.000204298017621671, 0.000245208408960445, 0.000207356569180511,
    0.000225972406234192, 0.000220283772032033, 0.00022354508456646,
    0.000243146128465111, 0.000228867327336145, 1.28449675748961E-5,
    0.000159163051704246, 0.000148666764487295, 0.000111089793023978,
    0.000216500353994512, 0.000208218039172982, 0.000183081961536816,
    2.18085726041116E-5, 2.10687396985976E-5, 2.66822982825485E-5,
    0.000138130489981132, 0.000171585315603236, 0.000170360927140262,
    7.37557739843605E-5, 9.59757116689576E-5, 0.00012015845596915,
    0.000122394274849401, 1.93764044323189E-5, 1.34340517305741E-5,
    1.51944484631687E-5, 0.000114442113052133, 0.000152498102014993,
    0.000112015891899157, 0.000188487105677467, 0.000208981700207685,
    0.000199387058015208, 1.02037229507265E-5, 9.03868992654702E-6,
    1.19935935007523E-5, 0.000539672247963784, 0.000524487709216202,
    0.000445763045671255, 0.000423216497322893, 0.000172394704882535,
    0.000127740646000638, 0.000126778773420547, 0.000182536489811769,
    0.000162188683253416, 0.00015862063916518, 0.000224942234267408,
    0.000177363260752988, 0.000300969401811449, 0.000180591044829688,
    8.50174309043072E-6, 6.77361009539538E-6, 4.16701577410318E-6,
    3.42194093278867E-5, 3.25538655544543E-5, 3.4888361360084E-5,
    0.000128472868450348, 0.000199526601194682, 0.000180022737109761,
    2.58522723230294E-5, 0.000269295860722968, 0.000307932196625213,
    0.000237630204856372, 9.93487530128404E-6, 7.93832985214355E-6,
    1.17801737985297E-5, 1.96844083941151E-5, 0.000474696165094188,
    0.000468867429460488, 0.000448305852164607, 0.000590904028819915,
    0.000142977048306222, 0.000156806845152992, 0.000161045296567369,
    1.20400059845105E-5, 1.2990524264968E-5, 1.44479409805678E-5,
    0.0953122764541685, 0.0656808515506427, 0.0407700269737749,
    -0.150348441503025, -0.0789910015600634, -0.00446903382787341,
    -0.0519630312982012, -0.0203275759983725, 0.0495350735500977,
    0.00480124575033053, -0.127558367593738, -0.0851554021545155,
    -0.0194339370363537, 0.117980215524134, 0.0965523146076427,
    0.00807199183670376, -0.166715409144615, -0.0970207766918774,
    -0.0310629109762221, -0.0466629912846003, -0.0315315206313729,
    0.00754505772790597, 0.133316142607978, 0.0568203114482877,
    0.0269140846216206, -0.140989260534467, -0.0898298822081123,
    -0.038624030884947, -0.0798831619515745, -0.0718786480509995,
    -0.0166220672954587, 0.148144929288623, 0.10041154709543,
    -0.0510319592979699, -0.136907864108615, -0.0503659510191789,
    -0.0285450071827639, -0.0857050163525366, -0.0477015411196316,
    -0.0486607900967999, -0.0178601263888384, 0.144651180527571,
    0.106909142579135, 0.0326677618942765, 0.0245029019687406,
    -0.100774847958715, -0.0759183477185552, -0.0204825698283515,
    0.178997476464521, 0.1097291748174, 0.0433254528619748, -0.0847403003413587,
    -0.0131185809981979, 0.0220168279783755, 0.161341495128655,
    0.0975845831787999, 0.0142941197360796, -0.15307176308466,
    -0.0891005807132297, -0.0264489014792431, 0.0290918438376946,
    0.00717878972001297, 0.0229388613730611, 0.146990367763805,
    0.132162695216896, 0.0618492904882786, 0.00473094436048555,
    -0.102104507602664, -0.0699490898516901, -0.0233653416276185,
    -0.0494542269246884, -0.00712945502739128, -0.0061006874887681,
    0.0957691828150444, 0.096103919410573, 0.0188296991340803,
    -0.145335050952785, -0.0976127693698435, -0.0334256102197836,
    0.0102849377547114, 0.0165066385020909, 0.0274941613282052,
    -0.150555055785698, -0.101360615165296, -0.0750272635171702,
    0.000564805087702514, 0.1359911967545, 0.0718231222260422, 0.033945840660848,
    0.209217360423698, 0.126881631034665, 0.0665686841654347, 0.0598170844032163,
    -0.105501531602768, -0.150535491935129, -0.00153440905854427,
    -0.0105049632862818, -0.00801552883154009, 0.00934612804767871,
    -0.0710415582061224, -0.0313115087274082, -0.0189510689664915,
    0.115544192346413, 0.107692016481558, 0.0149161107085507, 0.0189519205302542,
    -0.118543899989235, -0.0597282783667234, 0.0225323381280878,
    0.0124589763112464, 0.00534363957366982, 0.0177599214346167,
    0.046283095566213, -0.184881243273735, -0.168850157768001,
    -0.105342800685828, -0.0374693776660446, 0.120016087037157,
    0.084378396673115, 0.0340842796081616, 0.00251973901280574,
    0.0246877114758373, 0.0139915219739796, 2.46353090638117E-6,
    2.46386207343214E-6, 2.46084580382461E-6, -5.8998079923166E-7,
    -3.96277926658516E-7, -2.89522097125457E-8, 5.94627998905793E-8,
    -5.7619105082038E-8, -1.10705154597298E-7, -1.51842348365372E-8,
    2.73806178712843E-6, 2.73848125047993E-6, -2.76166533063913E-6,
    -1.2660657778867E-6, -1.2660027140864E-6, -3.01786974695435E-6,
    4.23548922643877E-7, -5.39177969648162E-7, -2.71506524374674E-7,
    9.32942149289486E-8, -2.86239011702028E-8, 8.53118926957893E-8,
    1.65592971197694E-6, 1.65599954131157E-6, 2.78942036507126E-7,
    -2.80598300072361E-6, -2.8069134907853E-6, -2.77825897003731E-6,
    2.57724873964094E-7, -1.27723968870768E-6, -1.27704480689098E-6,
    -2.28850899092727E-6, -2.2887464952658E-6, 2.50568322961496E-7,
    -1.15709623343544E-6, -1.15745225781771E-6, -1.15710254475001E-6,
    2.60900313634796E-8, 2.9524595723361E-8, 2.2996367341385E-7,
    2.30755998485729E-7, -5.89428781522879E-7, -5.89377042823957E-7,
    2.84746499292334E-7, 3.10610728671289E-7, 7.34608172745416E-7,
    -4.56423345032678E-7, -4.55472579143643E-7, -1.08864971964363E-6,
    -1.08872744184204E-6, -1.16081370795044E-6, -8.81910206960401E-7,
    -2.9841013546468E-7, -3.44539340951431E-8, -1.714892986992E-6,
    -1.71544773670647E-6, -1.72922454460937E-6, -7.65236566228253E-7,
    -7.64718996818473E-7, -7.6468575838431E-7, -2.14573649532233E-7,
    -2.15455532332042E-7, -1.96507693849703E-7, -4.52122046753072E-7,
    -4.52194501949934E-7, 3.06989260086679E-7, 6.29875553453854E-8,
    2.8131133300291E-7, -1.65910122439566E-6, -1.65919916243205E-6,
    -8.61633546129063E-8, -4.16404781911021E-8, -3.9967037560356E-8,
    4.22549659463173E-7, 4.20274533358338E-7, 2.30010465252064E-7,
    3.42085277713833E-7, -1.15800727203663E-6, -3.4674221943658E-7,
    -5.5163484583543E-8, -5.5152020865106E-8, -4.45542419990922E-8,
    1.04359645090112E-7, 1.04244716539359E-7, 1.03176390459778E-7,
    -2.06248620851814E-6, -3.01740968086436E-7, -1.89114255606426E-6,
    2.69813929026024E-7, -4.05939632305533E-7, -4.05731051427464E-7,
    -3.02662502997775E-6, 2.22852236136514E-7, -2.34867061361834E-6,
    -2.34679817032999E-6, -2.3483226224943E-6, -4.01494468869531E-8,
    -1.7539612212502E-10, -3.99324760436126E-8, 4.05209377080831E-7,
    -6.31470774628502E-8, -6.19957555493853E-8, -2.91860903284846E-7,
    3.59501841697162E-7, 3.59766217775768E-7, 3.93398043286199E-8,
    4.16349482288203E-7, -1.0343842083293E-6, -1.02869320670165E-6,
    -1.95501990616216E-7, -1.94693298413174E-7, -1.95105640886427E-7,
    -1.95873028024245E-7, 1.35386827032806E-6, 1.35436479899316E-6,
    1.35410590458555E-6, -7.80584534677463E-7, -7.7946587694538E-7,
    -7.7835809505105E-7, 7.60503485447878E-7, 4.99056192843825E-7,
    4.96692318122845E-7, 4.96427107230663E-7, 3.94164945020987E-5,
    3.94217931749142E-5, 3.93735328611937E-5, 6.19479839193243E-5,
    6.10268007054114E-5, 2.31617677700365E-6, 5.41111479004271E-6,
    7.60572187082902E-6, 3.87468041090544E-6, 4.40342810259579E-7,
    5.74992975296971E-5, 5.75081062600785E-5, 6.07566372740609E-5,
    6.45693546722216E-5, 6.45661384184064E-5, 6.63931344329957E-5,
    4.70139304134704E-5, 5.49961529041125E-5, 5.10432265824387E-5,
    5.78424132559481E-6, 6.068267048083E-6, 1.19436649774105E-6,
    5.96134696311697E-5, 5.96159834872164E-5, 5.66252334109465E-5,
    5.61196600144722E-5, 5.6138269815706E-5, 5.55651794007463E-5,
    2.11334396650557E-5, 5.2366827237015E-5, 5.23588370825302E-5,
    4.34816708276181E-5, 4.34861834100502E-5, 4.13437732886469E-5,
    5.20693305045946E-5, 5.20853516017968E-5, 5.20696145137504E-5,
    3.65260439088715E-6, 3.86772203976029E-6, 9.19854693655399E-6,
    9.23023993942916E-6, 5.01014464294447E-5, 5.00970486400364E-5,
    4.81221583804044E-5, 4.87658844013924E-5, 6.02378701651241E-5,
    6.48121149946402E-5, 6.46771062383973E-5, 5.44324859821814E-5,
    5.44363720921022E-5, 5.57190579816211E-5, 5.64422532454657E-5,
    5.58026953318952E-5, 2.44622932075516E-6, 3.42978597398401E-5,
    3.43089547341295E-5, 3.45844908921874E-5, 4.66794305399234E-5,
    4.66478588059269E-5, 4.66458312614429E-5, 3.86232569158019E-6,
    3.87819958197675E-6, 3.53713848929465E-6, 4.47600826285541E-5,
    4.47672556930435E-5, 4.02155930713549E-5, 5.73186753643007E-6,
    3.1225557963323E-5, 4.31366318342871E-5, 4.31391782232332E-5,
    5.85910811367763E-6, 5.12177881750556E-6, 4.95591265748414E-6,
    4.31000652652436E-5, 4.28680024025504E-5, 4.23219256063798E-5,
    5.19969622125026E-5, 5.21103272416484E-5, 5.16645906960504E-5,
    8.82615753336689E-7, 8.82432333841696E-7, 7.12867871985475E-7,
    6.05285941522647E-6, 6.04619355928281E-6, 5.98423064666713E-6,
    8.66244207577619E-5, 3.37949884256808E-5, 3.97139936773494E-5,
    3.91230197087736E-5, 4.70889973474419E-5, 4.70648019655858E-5,
    4.8426000479644E-5, 4.88046397138967E-5, 4.93220828859852E-5,
    4.92827615769298E-5, 4.93147750723802E-5, 1.00373617217383E-6,
    2.6484814440878E-8, 9.98311901090316E-7, 1.01302344270208E-5,
    9.34576746450183E-6, 9.17537182130903E-6, 3.56070302007512E-5,
    4.42187265287509E-5, 4.42512447864195E-5, 2.79312610733201E-6,
    4.70474914985669E-5, 6.9303741958063E-5, 6.78937516423092E-5,
    3.32353384047568E-6, 3.30978607302395E-6, 3.31679589506926E-6,
    3.32984147641217E-6, 4.46776529208258E-5, 4.46940383667742E-5,
    4.46854948513232E-5, 9.28895596266181E-5, 4.20911573550505E-5,
    4.20313371327567E-5, 4.41092021559769E-5, 7.48584289265737E-6,
    7.45038477184267E-6, 7.44640660845995E-6, 0.000382752106456979,
    0.000264071944283501, 0.000164533060349642, -0.000603698754344476,
    -0.000318357239172927, -1.88237820104731E-5, -0.000208579637889989,
    -8.19614818965966E-5, 0.000198074407161706, 1.9076300289834E-5,
    -0.000512206000357525, -0.000342556643353312, -7.96721518422797E-5,
    0.000473335777902639, 0.00038785742971784, 3.37111835509611E-5,
    -0.000668910010471307, -0.00039034379101188, -0.00012633130400178,
    -0.000187393236247721, -0.00012695346590492, 3.01888288014222E-5,
    0.000534763660280542, 0.000228397476978336, 0.000109258598955787,
    -0.000566048667290353, -0.000361199528181256, -0.000156691727907427,
    -0.000319929310919246, -0.000289156785351903, -6.79676673065141E-5,
    0.00059401014671937, 0.00040305833649398, -0.000203569133004302,
    -0.000549267080475896, -0.000202688137812675, -0.000115710088761428,
    -0.000343792622096549, -0.000191190569701329, -0.000195438634206962,
    -7.20315047962243E-5, 0.000579973631040588, 0.000428963599189485,
    0.000131822528432629, 9.96777311839611E-5, -0.000404565208393161,
    -0.000305749900953449, -8.38342753659654E-5, 0.000717831502169974,
    0.000440782906148366, 0.000175048060340606, -0.000340970532673593,
    -5.41210619811867E-5, 8.82509811784713E-5, 0.000647347225632431,
    0.000392099944503235, 5.84312159914034E-5, -0.000614245853148331,
    -0.000358228054146354, -0.000107447355358047, 0.000116543740817612,
    2.86387093424941E-5, 9.21802061331236E-5, 0.000589382437620292,
    0.000530525575176978, 0.000249220819498683, 2.06690267555507E-5,
    -0.000409518599115443, -0.000281150391964644, -9.48675510096378E-5,
    -0.00019846752323363, -2.89461309690558E-5, -2.48216400889187E-5,
    0.000384458596954684, 0.000386317913370411, 7.67626352907389E-5,
    -0.000582950617735156, -0.000392341057381362, -0.000135571733686533,
    4.13278728218442E-5, 6.61880655423934E-5, 0.000110331986975041,
    -0.000603465806734043, -0.000406775780111345, -0.00030090746349094,
    6.04132006636865E-8, 0.000545979389128485, 0.000288830678358521,
    0.000137277209498088, 0.000838843410894206, 0.000509123016847528,
    0.000267847174668345, 0.000241461923756936, -0.000423232302716889,
    -0.000604814487524714, -7.67451887543122E-6, -4.22356168396647E-5,
    -3.24066850950066E-5, 3.74769773805221E-5, -0.00028516646362569,
    -0.000125988536958167, -7.67802535226699E-5, 0.000463261913647381,
    0.000432704628779884, 6.11687498196116E-5, 7.61973427197873E-5,
    -0.000475556729652464, -0.000240915600597644, 8.8275034388294E-5,
    5.00033965302164E-5, 2.15837676824824E-5, 7.11375326637471E-5,
    0.000185694595265366, -0.000741360309124064, -0.00067718113696201,
    -0.00042307272782562, -0.000152638203608045, 0.000481549902413192,
    0.000339234604948969, 0.000138156085999502, 9.63775774996229E-6,
    9.90245803265093E-5, 5.62712712408816E-5, 0.000301718001893731,
    0.00025594285352202, 0.000153917987406359, 0.000547256861550094,
    0.00047198718695627, 7.38943535316071E-5, 0.000170921257528357,
    0.000154873953842744, 0.000147259257125295, 0.000115227413225606,
    0.000452643540344907, 0.000424116889222981, 0.000245345573824996,
    0.000379115174766691, 0.000401686152616924, 0.000288210265382201,
    0.000478321537254203, 0.00048108597864469, 0.000295948030745857,
    0.000132508149322709, 8.29982329643298E-5, 2.55660092041618E-5,
    0.00030420918601151, 0.000314429246944214, 0.000218386698788952,
    0.000441114747799592, 0.000402912667983546, 0.000221374670922135,
    0.000198433912448256, 0.000413745425533213, 0.000278233781704751,
    0.000304394688729298, 0.00036882668289744, 0.000299658438441149,
    0.000407689581802767, 0.000415696018001128, 0.000327223591722177,
    4.45575186476226E-5, 0.00012254632295692, 0.000151223569395584,
    0.000120334520367394, 0.000279553302189059, 0.000362435152225655,
    0.000359291803527093, 0.000191712941917795, 0.000439842371428468,
    0.000411989595012479, 0.000190814683869282, 0.000379225389077446,
    0.000457826080988526, 0.000397687652917474, 0.000450456858841406,
    0.000304089746216149, 5.05184876872357E-5, 0.000341444278357966,
    0.000367845369347319, 0.000271056976211636, 0.000417546900581164,
    0.000426269250728998, 0.00030028195849009, 8.64512364673746E-5,
    9.05538470688694E-5, 9.43857953558314E-5, 0.00030132271223912,
    0.000396652151484396, 0.000361765474312722, 0.000124492846959338,
    0.000287983860160082, 0.000318070270800761, 0.000224150497366114,
    0.000134336433263746, 0.000136043674691244, 9.43483645159842E-5,
    0.000325045667057044, 0.000311056691106834, 0.00015656890916083,
    0.000408130541208982, 0.000425903350004138, 0.000308802094726734,
    4.38235590317305E-5, 5.50774123116489E-5, 5.60698520950927E-5,
    0.000425607097993683, 0.000483391773497182, 0.000496633159893373,
    0.000319610224091737, 0.000312754043487601, 0.000340872172933542,
    0.000250227862745889, 0.000251327989892321, 0.000376659895917706,
    0.000377608590144629, 0.00023489314803782, 0.000440452073838323,
    0.000429646168116861, 0.00026762818848527, 6.22872397495189E-5,
    3.80069978115044E-5, 2.80686567193962E-5, 0.000184453423525063,
    0.000193137403147084, 0.000113277751076387, 0.000344402255397426,
    0.000374524898874369, 0.00029177494918451, 7.97078564404812E-5,
    0.000548306537818824, 0.000507262753153302, 0.000217078530850612,
    3.27070648458392E-5, 3.09606616775286E-5, 5.70444736757934E-5,
    6.96318545090334E-5, 0.000353572498628395, 0.000504711394491594,
    0.000545725729156156, 0.00049181454544673, 0.000359755692998569,
    0.000363583331254663, 0.000207169422589839, 5.16509835202149E-5,
    5.52958644575318E-5, 4.84684964515462E-5, 0.000755359173225023,
    0.000658709691241291, 0.000716407629497978, 5.78877535779075E-5,
    5.90067930323719E-5, 8.53653113453657E-5, 7.89808930564934E-5,
    9.99806881807596E-5, 0.000584955505856805, 0.000375551872062656,
    0.000101427165026454, 4.11922337730361E-5, 6.84385224401164E-5,
    0.000972956224172803, 0.00106022352399068, 0.000907001531072937,
    -6.55432594258678E-6, 3.61968274752575E-5, 3.62297086313629E-5,
    -3.2817437682983E-5, -4.65019454324183E-5, 7.45906779363219E-5,
    0.000899862002273611, 0.000734639447476746, 0.000758382398999678,
    3.2848683897515E-5, 1.69997951727887E-6, -1.26425255537999E-5,
    3.79631586173949E-5, 0.000174615502745947, 0.00018710561515928,
    0.000886291489315784, 0.000880516049772086, 0.000470038139361623,
    0.000124810572069464, 0.000266337115776505, 0.000146920027508906,
    -0.000184366299803602, 6.84627218804613E-5, -9.34846733226775E-7,
    3.71135427258973E-5, 0.000787149469221217, 0.000805906114069463,
    0.000704592513043838, 0.000859205252909083, 0.000164586326458022,
    5.24341377000444E-5, 5.29504665119023E-5, 0.00106807675738134,
    0.00101713713943374, 0.000964154137588878, 9.23587719451876E-5,
    0.00014751205516552, 0.000199593825699964, 0.000941780717806157,
    0.000871772285394588, 0.000738570755096358, 3.74019908707538E-5,
    6.23895755512627E-5, 9.4059103490837E-5, 0.000402149819911724,
    0.000322171193305693, 0.000410121605262804, 0.000802352739812971,
    0.00094717128380593, 0.000894915749257183, 0.000864538444773375,
    2.35930259893133E-5, 3.21773043148623E-5, 4.48019316435363E-5,
    3.45596952988238E-5, 0.000149157693674622, 9.00752246693207E-5,
    0.000718388418215764, 0.000893822564317687, 0.000743283621069255,
    5.38844420826629E-5, 2.58812266733562E-5, 7.3567442411976E-5,
    0.00015920445902748, 0.000211148397335097, 0.000279042769545685,
    0.00014168747246897, 0.00018512077472567, 0.000152896569641086,
    0.000266190581586479, 0.00100093507277309, 0.000776186431552911,
    0.000782028633188605, 0.00106623952289007, 0.000920529790695767,
    0.000891901582721868, 0.00109431512970648, 0.000209525360639591,
    -0.000147520905739156, 0.000161869391342701, 8.18290364141335E-5,
    6.07329040237816E-5, 0.000121080049594614, 5.47650170434365E-6,
    9.03568966213395E-5, 3.78240775292962E-5, 0.000859405300564658,
    0.000959441142243204, 0.000785787128954371, 0.00028538782220271,
    0.000228035655186529, 0.000217636772957426, 0.000273343507650793,
    0.000137638542310521, 0.000112857957579468, 0.000188138251502066,
    0.000341309732210291, 6.50449812612223E-5, 1.71174612964523E-5,
    4.6126475495073E-5, 9.93665376498025E-5, 0.000854701730053457,
    0.000899370798972247, 0.000903463168769567, 0.000110161627294974,
    0.000209379772689451, 0.000194129746519134, -3.01359887669774E-5,
    1.69084673953959E-5, 1.70650018527301E-5, -0.00118630453014227,
    -0.00118424870928789, -0.0004803891028962, -0.000412833753042408,
    -0.000342993885964701, -0.000177867904068145, -7.82635753991569E-5,
    -0.000990327950911438, -0.00102648497529612, -0.00100051799922159,
    -0.000104812677404259, -1.75218486638169E-5, -0.000173053830270444,
    -0.00115079407285502, -0.00113831364312573, -0.00105014222039853,
    -0.00045308906498031, -0.000366622230879916, -2.89143056598972E-5,
    2.355109183051E-5, -0.000141698692295282, -0.000116406596832681,
    -0.00102980021998476, -0.000987570641604159, -0.00105048944085769,
    -0.000507636549491189, -0.000868997634627863, -0.000768852049861936,
    1.78860929784877E-5, -1.67041156340299E-5, -0.000440025435657483,
    -0.000863223818963389, -0.000722267890154075, -0.000841215944597356,
    -0.000425206370099131, -0.00033692791874419, -0.000404213088896264,
    -0.000325251630280664, -6.83719783504923E-5, -4.96012717571638E-5,
    -0.000219106049042222, -4.98031927927823E-5, -0.000869247149127775,
    -0.00102368528157391, -0.00092492527924105, 8.65717533548863E-7,
    -9.46416574173234E-5, -0.000147463016755821, -0.000994457049547743,
    -0.000891907532314718, -3.69381948998536E-6, 4.88418410982121E-5,
    6.88960091023317E-6, -0.000128306649467623, -0.000991463893265437,
    -0.000947650528419984, -0.000852587689505263, 1.70244724895976E-5,
    -5.17075032887164E-5, -2.57530911114936E-6, -5.89981662238746E-5,
    2.69054062600653E-5, -4.71320545410765E-5, -7.5954347343037E-5,
    -0.000719505410292904, -0.00072967130943244, -0.000709405468150263,
    -0.00037539120172542, -0.000243625107585021, -0.000255197665854148,
    -5.46376511579067E-5, 7.27235301937431E-5, -3.45694468491223E-5,
    -0.000911592820988773, -0.000936193021450031, -0.000910179262465929,
    -2.67023148948449E-5, -2.17429690044456E-5, 4.02954472178921E-5,
    -0.00130672992102559, -0.0012631320938846, -0.00129561038477857,
    -0.00116967217737745, 6.36598246983915E-6, -8.7399668798949E-5,
    -6.27548413820733E-5, 0.000135860523357692, -2.32008430130184E-5,
    -7.64350135723489E-5, 1.14054001650242E-5, -0.000847040956974396,
    -0.0012326355343471, -0.000799653320633223, -0.000160577823736911,
    -0.000155424036089776, -5.53184198453717E-6, -0.000500172806444677,
    -0.000412342123123543, -0.00045525365856923, -0.000118940539791969,
    -1.78401060540429E-5, -0.000196976098675094, -5.12482956557939E-5,
    -0.00108756220751823, -0.00107951583931925, -0.000959600523663226,
    -2.68574644381588E-5, -5.16463822400095E-5, -7.52998926080929E-5,
    3.90311430681926E-5, -0.00148040457098098, -0.0016722199415527,
    -0.00164291658164034, -0.00158025413052086, -2.26468523238464E-5,
    2.07157130381716E-5, 1.3477971074244E-5, -0.00020797866516254,
    -4.15434057325223E-6, -2.00010638991177E-5, 0.0953122764541685,
    0.0656808515506427, 0.0407700269737749, -0.150348441503025,
    -0.0789910015600634, -0.00446903382787341, -0.0519630312982012,
    -0.0203275759983725, 0.0495350735500977, 0.00480124575033053,
    -0.127558367593738, -0.0851554021545155, -0.0194339370363537,
    0.117980215524134, 0.0965523146076427, 0.00807199183670376,
    -0.166715409144615, -0.0970207766918774, -0.0310629109762221,
    -0.0466629912846003, -0.0315315206313729, 0.00754505772790597,
    0.133316142607978, 0.0568203114482877, 0.0269140846216206,
    -0.140989260534467, -0.0898298822081123, -0.038624030884947,
    -0.0798831619515745, -0.0718786480509995, -0.0166220672954587,
    0.148144929288623, 0.10041154709543, -0.0510319592979699, -0.136907864108615,
    -0.0503659510191789, -0.0285450071827639, -0.0857050163525366,
    -0.0477015411196316, -0.0486607900967999, -0.0178601263888384,
    0.144651180527571, 0.106909142579135, 0.0326677618942765, 0.0245029019687406,
    -0.100774847958715, -0.0759183477185552, -0.0204825698283515,
    0.178997476464521, 0.1097291748174, 0.0433254528619748, -0.0847403003413587,
    -0.0131185809981979, 0.0220168279783755, 0.161341495128655,
    0.0975845831787999, 0.0142941197360796, -0.15307176308466,
    -0.0891005807132297, -0.0264489014792431, 0.0290918438376946,
    0.00717878972001297, 0.0229388613730611, 0.146990367763805,
    0.132162695216896, 0.0618492904882786, 0.00473094436048555,
    -0.102104507602664, -0.0699490898516901, -0.0233653416276185,
    -0.0494542269246884, -0.00712945502739128, -0.0061006874887681,
    0.0957691828150444, 0.096103919410573, 0.0188296991340803,
    -0.145335050952785, -0.0976127693698435, -0.0334256102197836,
    0.0102849377547114, 0.0165066385020909, 0.0274941613282052,
    -0.150555055785698, -0.101360615165296, -0.0750272635171702,
    0.000564805087702514, 0.1359911967545, 0.0718231222260422, 0.033945840660848,
    0.209217360423698, 0.126881631034665, 0.0665686841654347, 0.0598170844032163,
    -0.105501531602768, -0.150535491935129, -0.00153440905854427,
    -0.0105049632862818, -0.00801552883154009, 0.00934612804767871,
    -0.0710415582061224, -0.0313115087274082, -0.0189510689664915,
    0.115544192346413, 0.107692016481558, 0.0149161107085507, 0.0189519205302542,
    -0.118543899989235, -0.0597282783667234, 0.0225323381280878,
    0.0124589763112464, 0.00534363957366982, 0.0177599214346167,
    0.046283095566213, -0.184881243273735, -0.168850157768001,
    -0.105342800685828, -0.0374693776660446, 0.120016087037157,
    0.084378396673115, 0.0340842796081616, 0.00251973901280574,
    0.0246877114758373, 0.0139915219739796, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 28.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 98.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  double dist[123];
  int b_k;
  emxArray_real_T *edges;
  double s;
  int b_idx[123];
  int exitg4;
  int exitg3;
  emxArray_real_T *nn;
  signed char outsize_idx_0;
  boolean_T guard1 = false;
  int exitg2;
  int mid_i;
  boolean_T exitg1;
  emxInit_real_T1(&Uc, 1);
  high_i = Uc->size[0];
  Uc->size[0] = 123;
  emxEnsureCapacity((emxArray__common *)Uc, high_i, (int)sizeof(double));
  for (high_i = 0; high_i < 123; high_i++) {
    Uc->data[high_i] = a[idx[high_i] - 1];
  }

  nb = -1;
  high_i = 0;
  while (high_i + 1 <= 123) {
    nbins = (int)Uc->data[high_i];
    do {
      exitg5 = 0;
      high_i++;
      if (high_i + 1 > 123) {
        exitg5 = 1;
      } else {
        eok = (fabs((double)nbins - Uc->data[high_i]) < eps((double)nbins /
                2.0));
        if (!eok) {
          exitg5 = 1;
        }
      }
    } while (exitg5 == 0);

    nb++;
    Uc->data[nb] = nbins;
  }

  emxInit_real_T(&b_tsX, 2);
  high_i = Uc->size[0];
  if (1 > nb + 1) {
    i2 = -1;
  } else {
    i2 = nb;
  }

  Uc->size[0] = i2 + 1;
  emxEnsureCapacity((emxArray__common *)Uc, high_i, (int)sizeof(double));
  nbins = tsX->size[1];
  high_i = b_tsX->size[0] * b_tsX->size[1];
  b_tsX->size[0] = 1;
  b_tsX->size[1] = nbins;
  emxEnsureCapacity((emxArray__common *)b_tsX, high_i, (int)sizeof(double));
  for (high_i = 0; high_i < nbins; high_i++) {
    b_tsX->data[b_tsX->size[0] * high_i] = tsX->data[tsX->size[0] * high_i];
  }

  emxInit_real_T(&r4, 2);
  high_i = r4->size[0] * r4->size[1];
  r4->size[0] = 123;
  r4->size[1] = b_tsX->size[1];
  emxEnsureCapacity((emxArray__common *)r4, high_i, (int)sizeof(double));
  for (high_i = 0; high_i < 123; high_i++) {
    nbins = b_tsX->size[1];
    for (nb = 0; nb < nbins; nb++) {
      r4->data[high_i + r4->size[0] * nb] = b_tsX->data[b_tsX->size[0] * nb];
    }
  }

  emxFree_real_T(&b_tsX);
  for (high_i = 0; high_i < 33; high_i++) {
    for (nb = 0; nb < 123; nb++) {
      x[nb + 123 * high_i] = tX[nb + 123 * high_i] - r4->data[nb + 123 * high_i];
    }
  }

  emxFree_real_T(&r4);

//#pragma omp parallel for \
// num_threads(omp_get_max_threads()) \
// private(b_k)

  for (k = 1; k < 4060; k++) {
    b_k = k;
    y[b_k - 1] = x[b_k - 1] * x[b_k - 1];
  }

  for (nbins = 0; nbins < 123; nbins++) {
    s = y[nbins];
    for (high_i = 0; high_i < 32; high_i++) {
      s += y[nbins + (high_i + 1) * 123];
    }

    dist[nbins] = s;
  }

  emxInit_real_T(&edges, 2);
  c_sort(dist, b_idx);
  nbins = Uc->size[0];
  high_i = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = (signed char)(nbins + 1);
  emxEnsureCapacity((emxArray__common *)edges, high_i, (int)sizeof(double));
  high_i = 0;
  do {
    exitg4 = 0;
    nbins = Uc->size[0];
    if (high_i <= nbins - 2) {
      edges->data[1 + high_i] = Uc->data[high_i] + (Uc->data[1 + high_i] -
        Uc->data[high_i]) / 2.0;
      high_i++;
    } else {
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  edges->data[0] = rtMinusInf;
  edges->data[edges->size[1] - 1] = rtInf;
  high_i = 1;
  do {
    exitg3 = 0;
    nbins = Uc->size[0];
    if (high_i - 1 <= nbins - 2) {
      edges->data[high_i] += eps(edges->data[high_i]);
      high_i++;
    } else {
      exitg3 = 1;
    }
  } while (exitg3 == 0);

  emxInit_real_T1(&nn, 1);
  outsize_idx_0 = (signed char)edges->size[1];
  high_i = nn->size[0];
  nn->size[0] = outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)nn, high_i, (int)sizeof(double));
  nbins = outsize_idx_0;
  for (high_i = 0; high_i < nbins; high_i++) {
    nn->data[high_i] = 0.0;
  }

  nbins = edges->size[1];
  guard1 = false;
  if (nbins > 1) {
    nb = 1;
    do {
      exitg2 = 0;
      if (nb + 1 <= nbins) {
        if (!(edges->data[nb] >= edges->data[nb - 1])) {
          eok = false;
          exitg2 = 1;
        } else {
          nb++;
        }
      } else {
        guard1 = true;
        exitg2 = 1;
      }
    } while (exitg2 == 0);
  } else {
    guard1 = true;
  }

  if (guard1) {
    eok = true;
  }

  if (!eok) {
    high_i = nn->size[0];
    nn->size[0] = outsize_idx_0;
    emxEnsureCapacity((emxArray__common *)nn, high_i, (int)sizeof(double));
    nbins = outsize_idx_0;
    for (high_i = 0; high_i < nbins; high_i++) {
      nn->data[high_i] = rtNaN;
    }
  } else {
    nbins = 0;
    if ((a[b_idx[0] - 1] >= edges->data[0]) && (a[b_idx[0] - 1] < edges->
         data[edges->size[1] - 1])) {
      nbins = 1;
      nb = 2;
      high_i = edges->size[1];
      while (high_i > nb) {
        mid_i = (nbins >> 1) + (high_i >> 1);
        if (((nbins & 1) == 1) && ((high_i & 1) == 1)) {
          mid_i++;
        }

        if (a[b_idx[0] - 1] >= edges->data[mid_i - 1]) {
          nbins = mid_i;
          nb = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }
    }

    if (a[b_idx[0] - 1] == edges->data[edges->size[1] - 1]) {
      nbins = edges->size[1];
    }

    if (nbins > 0) {
      nn->data[nbins - 1] = 1.0;
    }
  }

  high_i = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = nn->size[0] - 1;
  emxEnsureCapacity((emxArray__common *)edges, high_i, (int)sizeof(double));
  for (high_i = 0; high_i <= nn->size[0] - 2; high_i++) {
    edges->data[high_i] = nn->data[high_i];
  }

  if (nn->size[0] - 1 > 0) {
    edges->data[edges->size[1] - 1] += nn->data[nn->size[0] - 1];
  }

  emxFree_real_T(&nn);
  nbins = 1;
  nb = edges->size[1];
  s = edges->data[0];
  high_i = 0;
  if (edges->size[1] > 1) {
    if (rtIsNaN(edges->data[0])) {
      mid_i = 2;
      exitg1 = false;
      while ((!exitg1) && (mid_i <= nb)) {
        nbins = mid_i;
        if (!rtIsNaN(edges->data[mid_i - 1])) {
          s = edges->data[mid_i - 1];
          high_i = mid_i - 1;
          exitg1 = true;
        } else {
          mid_i++;
        }
      }
    }

    if (nbins < edges->size[1]) {
      while (nbins + 1 <= nb) {
        if (edges->data[nbins] > s) {
          s = edges->data[nbins];
          high_i = nbins;
        }

        nbins++;
      }
    }
  }

  emxFree_real_T(&edges);
  yfit = Uc->data[high_i];
  emxFree_real_T(&Uc);
  return yfit;
}

//
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
//                int np
//                int nq
//                emxArray_int32_T *iwork
//                emxArray_real_T *xwork
// Return Type  : void
//
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if ((np == 0) || (nq == 0)) {
  } else {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork->data[qend] = idx->data[offset + qend];
      xwork->data[qend] = x->data[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork->data[p] >= xwork->data[n]) {
        idx->data[iout] = iwork->data[p];
        x->data[iout] = xwork->data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx->data[iout] = iwork->data[n];
        x->data[iout] = xwork->data[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = iout - p;
          while (p + 1 <= np) {
            idx->data[(n + p) + 1] = iwork->data[p];
            x->data[(n + p) + 1] = xwork->data[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
//                int n
//                int preSortLevel
//                emxArray_int32_T *iwork
//                emxArray_real_T *xwork
// Return Type  : void
//
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : emxArray_real_T *A
//                const emxArray_real_T *B
// Return Type  : void
//
static void mrdivide(emxArray_real_T *A, const emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_real_T *b_B;
  int minmn;
  int maxmn;
  int rankR;
  double tol;
  int mn;
  double wj;
  emxInit_real_T1(&b_A, 1);
  emxInit_real_T1(&tau, 1);
  emxInit_real_T1(&b_B, 1);
  if ((A->size[1] == 0) || (B->size[1] == 0)) {
    minmn = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)A, minmn, (int)sizeof(double));
    A->data[0] = 0.0;
  } else if (1 == B->size[1]) {
    if (A->size[1] == 0) {
    } else {
      A->data[0] *= 1.0 / B->data[0];
    }
  } else {
    minmn = b_A->size[0];
    b_A->size[0] = B->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, minmn, (int)sizeof(double));
    maxmn = B->size[1];
    for (minmn = 0; minmn < maxmn; minmn++) {
      b_A->data[minmn] = B->data[B->size[0] * minmn];
    }

    xgeqp3(b_A, tau, &minmn);
    rankR = 0;
    if (b_A->size[0] < 1) {
      minmn = 0;
      maxmn = 1;
    } else {
      minmn = 1;
      maxmn = b_A->size[0];
    }

    if (minmn > 0) {
      tol = (double)maxmn * fabs(b_A->data[0]) * 2.2204460492503131E-16;
      while ((rankR < 1) && (fabs(b_A->data[0]) >= tol)) {
        rankR = 1;
      }
    }

    tol = 0.0;
    minmn = b_B->size[0];
    b_B->size[0] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_B, minmn, (int)sizeof(double));
    maxmn = A->size[1];
    for (minmn = 0; minmn < maxmn; minmn++) {
      b_B->data[minmn] = A->data[A->size[0] * minmn];
    }

    maxmn = b_A->size[0];
    mn = !(b_A->size[0] < 1);
    minmn = 1;
    while (minmn <= mn) {
      if (tau->data[0] != 0.0) {
        wj = b_B->data[0];
        for (minmn = 1; minmn + 1 <= maxmn; minmn++) {
          wj += b_A->data[minmn] * b_B->data[minmn];
        }

        wj *= tau->data[0];
        if (wj != 0.0) {
          b_B->data[0] -= wj;
          for (minmn = 1; minmn + 1 <= maxmn; minmn++) {
            b_B->data[minmn] -= b_A->data[minmn] * wj;
          }
        }
      }

      minmn = 2;
    }

    minmn = 1;
    while (minmn <= rankR) {
      tol = b_B->data[0];
      minmn = 2;
    }

    while (rankR > 0) {
      tol /= b_A->data[0];
      rankR = 0;
    }

    minmn = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)A, minmn, (int)sizeof(double));
    A->data[0] = tol;
  }

  emxFree_real_T(&b_B);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

//
// Arguments    : const double a[249]
//                double y[249]
// Return Type  : void
//
static void power(const double a[249], double y[249])
{
  int k;
  for (k = 0; k < 249; k++) {
    y[k] = a[k] * a[k];
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

//
// Arguments    : int *k
//                const emxArray_real_T *x
// Return Type  : double
//
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x)
{
  double xk;
  boolean_T exitg1;
  boolean_T p;
  xk = x->data[*k - 1];
  exitg1 = false;
  while ((!exitg1) && (*k < x->size[0])) {
    if ((fabs(xk - x->data[*k]) < eps(xk / 2.0)) || (rtIsInf(x->data[*k]) &&
         rtIsInf(xk) && ((x->data[*k] > 0.0) == (xk > 0.0)))) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }

  return xk;
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
static void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  b_sort(x, dim, idx);
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_real_T *b_x;
  unsigned int unnamed_idx_0;
  int ib;
  int m;
  int n;
  double x4[4];
  int idx4[4];
  emxArray_int32_T *iwork;
  emxArray_real_T *xwork;
  int nNaNs;
  int k;
  int wOffset;
  signed char perm[4];
  int nNonNaN;
  int p;
  int i4;
  int nBlocks;
  int b_iwork[256];
  double b_xwork[256];
  int b;
  int bLen;
  int bLen2;
  int nPairs;
  int exitg1;
  emxInit_real_T1(&b_x, 1);
  unnamed_idx_0 = (unsigned int)x->size[0];
  ib = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, ib, (int)sizeof(double));
  m = x->size[0];
  for (ib = 0; ib < m; ib++) {
    b_x->data[ib] = x->data[ib];
  }

  ib = idx->size[0];
  idx->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, ib, (int)sizeof(int));
  m = (int)unnamed_idx_0;
  for (ib = 0; ib < m; ib++) {
    idx->data[ib] = 0;
  }

  n = x->size[0];
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  emxInit_int32_T(&iwork, 1);
  ib = iwork->size[0];
  iwork->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)iwork, ib, (int)sizeof(int));
  m = iwork->size[0];
  ib = iwork->size[0];
  iwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)iwork, ib, (int)sizeof(int));
  for (ib = 0; ib < m; ib++) {
    iwork->data[ib] = 0;
  }

  emxInit_real_T1(&xwork, 1);
  unnamed_idx_0 = (unsigned int)x->size[0];
  ib = xwork->size[0];
  xwork->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)xwork, ib, (int)sizeof(double));
  m = xwork->size[0];
  ib = xwork->size[0];
  xwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)xwork, ib, (int)sizeof(double));
  for (ib = 0; ib < m; ib++) {
    xwork->data[ib] = 0.0;
  }

  nNaNs = 0;
  ib = 0;
  for (k = 0; k + 1 <= n; k++) {
    if (rtIsNaN(b_x->data[k])) {
      idx->data[(n - nNaNs) - 1] = k + 1;
      xwork->data[(n - nNaNs) - 1] = b_x->data[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = b_x->data[k];
      if (ib == 4) {
        ib = k - nNaNs;
        if (x4[0] >= x4[1]) {
          m = 1;
          wOffset = 2;
        } else {
          m = 2;
          wOffset = 1;
        }

        if (x4[2] >= x4[3]) {
          p = 3;
          i4 = 4;
        } else {
          p = 4;
          i4 = 3;
        }

        if (x4[m - 1] >= x4[p - 1]) {
          if (x4[wOffset - 1] >= x4[p - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)p;
            perm[3] = (signed char)i4;
          } else if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)p;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)m;
            perm[1] = (signed char)p;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else if (x4[m - 1] >= x4[i4 - 1]) {
          if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)p;
            perm[1] = (signed char)m;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)p;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else {
          perm[0] = (signed char)p;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)m;
          perm[3] = (signed char)wOffset;
        }

        idx->data[ib - 3] = idx4[perm[0] - 1];
        idx->data[ib - 2] = idx4[perm[1] - 1];
        idx->data[ib - 1] = idx4[perm[2] - 1];
        idx->data[ib] = idx4[perm[3] - 1];
        b_x->data[ib - 3] = x4[perm[0] - 1];
        b_x->data[ib - 2] = x4[perm[1] - 1];
        b_x->data[ib - 1] = x4[perm[2] - 1];
        b_x->data[ib] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  wOffset = (x->size[0] - nNaNs) - 1;
  if (ib > 0) {
    for (m = 0; m < 4; m++) {
      perm[m] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] >= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] >= x4[1]) {
      if (x4[1] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
      b_x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
    }
  }

  m = nNaNs >> 1;
  for (k = 1; k <= m; k++) {
    ib = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[n - k];
    idx->data[n - k] = ib;
    b_x->data[wOffset + k] = xwork->data[n - k];
    b_x->data[n - k] = xwork->data[wOffset + k];
  }

  if ((nNaNs & 1) != 0) {
    b_x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
  }

  nNonNaN = x->size[0] - nNaNs;
  m = 2;
  if (nNonNaN > 1) {
    if (x->size[0] >= 256) {
      nBlocks = nNonNaN >> 8;
      if (nBlocks > 0) {
        for (i4 = 1; i4 <= nBlocks; i4++) {
          n = (i4 - 1) << 8;
          for (b = 0; b < 6; b++) {
            bLen = 1 << (b + 2);
            bLen2 = bLen << 1;
            nPairs = 256 >> (b + 3);
            for (k = 1; k <= nPairs; k++) {
              m = n + (k - 1) * bLen2;
              for (ib = 0; ib + 1 <= bLen2; ib++) {
                b_iwork[ib] = idx->data[m + ib];
                b_xwork[ib] = b_x->data[m + ib];
              }

              p = 0;
              wOffset = bLen;
              ib = m - 1;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[p] >= b_xwork[wOffset]) {
                  idx->data[ib] = b_iwork[p];
                  b_x->data[ib] = b_xwork[p];
                  if (p + 1 < bLen) {
                    p++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx->data[ib] = b_iwork[wOffset];
                  b_x->data[ib] = b_xwork[wOffset];
                  if (wOffset + 1 < bLen2) {
                    wOffset++;
                  } else {
                    ib = (ib - p) + 1;
                    while (p + 1 <= bLen) {
                      idx->data[ib + p] = b_iwork[p];
                      b_x->data[ib + p] = b_xwork[p];
                      p++;
                    }

                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }

        m = nBlocks << 8;
        ib = nNonNaN - m;
        if (ib > 0) {
          merge_block(idx, b_x, m, ib, 2, iwork, xwork);
        }

        m = 8;
      }
    }

    merge_block(idx, b_x, 0, nNonNaN, m, iwork, xwork);
  }

  if ((nNaNs > 0) && (nNonNaN > 0)) {
    for (k = 0; k + 1 <= nNaNs; k++) {
      xwork->data[k] = b_x->data[nNonNaN + k];
      iwork->data[k] = idx->data[nNonNaN + k];
    }

    for (k = nNonNaN - 1; k + 1 > 0; k--) {
      b_x->data[nNaNs + k] = b_x->data[k];
      idx->data[nNaNs + k] = idx->data[k];
    }

    for (k = 0; k + 1 <= nNaNs; k++) {
      b_x->data[k] = xwork->data[k];
      idx->data[k] = iwork->data[k];
    }
  }

  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
  ib = x->size[0];
  x->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)x, ib, (int)sizeof(double));
  m = b_x->size[0];
  for (ib = 0; ib < m; ib++) {
    x->data[ib] = b_x->data[ib];
  }

  emxFree_real_T(&b_x);
}

//
// Arguments    : const double x[250]
// Return Type  : double
//
static double trapz(const double x[250])
{
  double z;
  int iy;
  double ylast;
  int k;
  z = 0.0;
  iy = 0;
  ylast = x[0];
  for (k = 0; k < 249; k++) {
    iy++;
    z += (ylast + x[iy]) / 2.0;
    ylast = x[iy];
  }

  return z;
}

//
// Arguments    : emxArray_real_T *A
//                emxArray_real_T *tau
//                int *jpvt
// Return Type  : void
//
static void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, int *jpvt)
{
  int m;
  int mn;
  int pvt;
  double vn1;
  double vn2;
  int i;
  int i_i;
  int mmi;
  int iy;
  double t;
  int k;
  double d1;
  double xnorm;
  double beta1;
  m = A->size[0];
  if (A->size[0] <= 1) {
    mn = A->size[0];
  } else {
    mn = 1;
  }

  pvt = tau->size[0];
  tau->size[0] = mn;
  emxEnsureCapacity((emxArray__common *)tau, pvt, (int)sizeof(double));
  if (A->size[0] == 0) {
  } else {
    vn1 = xnrm2(A->size[0], A, 1);
    vn2 = vn1;
    for (i = 1; i <= mn; i++) {
      i_i = (i + (i - 1) * m) - 1;
      mmi = m - i;
      pvt = (i + ixamax(2 - i, vn1)) - 2;
      if (pvt + 1 != i) {
        pvt *= m;
        iy = m * (i - 1);
        for (k = 1; k <= m; k++) {
          xnorm = A->data[pvt];
          A->data[pvt] = A->data[iy];
          A->data[iy] = xnorm;
          pvt++;
          iy++;
        }
      }

      if (i < m) {
        t = A->data[i_i];
        d1 = 0.0;
        if (1 + mmi <= 0) {
        } else {
          xnorm = b_xnrm2(mmi, A, i_i + 2);
          if (xnorm != 0.0) {
            beta1 = rt_hypotd_snf(A->data[i_i], xnorm);
            if (A->data[i_i] >= 0.0) {
              beta1 = -beta1;
            }

            if (fabs(beta1) < 1.0020841800044864E-292) {
              pvt = 0;
              do {
                pvt++;
                xscal(mmi, 9.9792015476736E+291, A, i_i + 2);
                beta1 *= 9.9792015476736E+291;
                t *= 9.9792015476736E+291;
              } while (!(fabs(beta1) >= 1.0020841800044864E-292));

              xnorm = b_xnrm2(mmi, A, i_i + 2);
              beta1 = rt_hypotd_snf(t, xnorm);
              if (t >= 0.0) {
                beta1 = -beta1;
              }

              d1 = (beta1 - t) / beta1;
              xscal(mmi, 1.0 / (t - beta1), A, i_i + 2);
              for (k = 1; k <= pvt; k++) {
                beta1 *= 1.0020841800044864E-292;
              }

              t = beta1;
            } else {
              d1 = (beta1 - A->data[i_i]) / beta1;
              xscal(mmi, 1.0 / (A->data[i_i] - beta1), A, i_i + 2);
              t = beta1;
            }
          }
        }

        tau->data[i - 1] = d1;
        A->data[i_i] = t;
      } else {
        tau->data[i - 1] = 0.0;
      }

      pvt = i;
      while (pvt + 1 < 2) {
        pvt = i + m * pvt;
        if (vn1 != 0.0) {
          xnorm = fabs(A->data[i - 1]) / vn1;
          xnorm = 1.0 - xnorm * xnorm;
          if (xnorm < 0.0) {
            xnorm = 0.0;
          }

          beta1 = vn1 / vn2;
          beta1 = xnorm * (beta1 * beta1);
          if (beta1 <= 1.4901161193847656E-8) {
            if (i < m) {
              vn1 = 0.0;
              if (mmi < 1) {
              } else if (mmi == 1) {
                vn1 = fabs(A->data[pvt]);
              } else {
                xnorm = 2.2250738585072014E-308;
                iy = pvt + mmi;
                while (pvt + 1 <= iy) {
                  beta1 = fabs(A->data[pvt]);
                  if (beta1 > xnorm) {
                    t = xnorm / beta1;
                    vn1 = 1.0 + vn1 * t * t;
                    xnorm = beta1;
                  } else {
                    t = beta1 / xnorm;
                    vn1 += t * t;
                  }

                  pvt++;
                }

                vn1 = xnorm * std::sqrt(vn1);
              }

              vn2 = vn1;
            } else {
              vn1 = 0.0;
              vn2 = 0.0;
            }
          } else {
            vn1 *= std::sqrt(xnorm);
          }
        }

        pvt = 1;
      }
    }
  }

  *jpvt = 1;
}

//
// Arguments    : int n
//                const emxArray_real_T *x
//                int ix0
// Return Type  : double
//
static double xnrm2(int n, const emxArray_real_T *x, int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = fabs(x->data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = fabs(x->data[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * std::sqrt(y);
  }

  return y;
}

//
// Arguments    : int n
//                double a
//                emxArray_real_T *x
//                int ix0
// Return Type  : void
//
static void xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i4;
  int k;
  i4 = (ix0 + n) - 1;
  for (k = ix0; k <= i4; k++) {
    x->data[k - 1] *= a;
  }
}

//
// Arguments    : const double ch1[250]
//                const double ch2[250]
//                const double ch3[250]
// Return Type  : double
//
double eog_knn(const double ch1[250], const double ch2[250], const double ch3
               [250])
{
  double Y;
  double ch1f[250];
  double ch2f[250];
  double ch3f[250];
  double adchf3[249];
  double y[249];
  int ixstart;
  double adchf1[249];
  double adchf2[249];
  double mtmp;
  int ix;
  boolean_T exitg6;
  boolean_T guard2 = false;
  boolean_T b0;
  boolean_T exitg5;
  boolean_T exitg3;
  boolean_T guard1 = false;
  boolean_T exitg4;
  boolean_T b1;
  boolean_T exitg2;
  emxArray_real_T *fch1f;
  emxArray_real_T *fch2f;
  emxArray_real_T *fch3f;
  emxArray_real_T *b_fch1f;
  boolean_T exitg1;
  Y = 0.0;
  customFilt(ch1, ch1f);
  customFilt(ch2, ch2f);
  customFilt(ch3, ch3f);
  diff(ch1f, adchf3);
  for (ixstart = 0; ixstart < 249; ixstart++) {
    y[ixstart] = fabs(adchf3[ixstart]);
  }

  power(y, adchf1);
  diff(ch2f, adchf3);
  for (ixstart = 0; ixstart < 249; ixstart++) {
    y[ixstart] = fabs(adchf3[ixstart]);
  }

  power(y, adchf2);
  diff(ch3f, adchf3);
  for (ixstart = 0; ixstart < 249; ixstart++) {
    y[ixstart] = fabs(adchf3[ixstart]);
  }

  power(y, adchf3);
  ixstart = 1;
  mtmp = adchf1[0];
  if (rtIsNaN(adchf1[0])) {
    ix = 2;
    exitg6 = false;
    while ((!exitg6) && (ix < 250)) {
      ixstart = ix;
      if (!rtIsNaN(adchf1[ix - 1])) {
        mtmp = adchf1[ix - 1];
        exitg6 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 249) {
    while (ixstart + 1 < 250) {
      if (adchf1[ixstart] > mtmp) {
        mtmp = adchf1[ixstart];
      }

      ixstart++;
    }
  }

  guard2 = false;
    LOGE("CH1 MTMP1 [%f]",mtmp);
  if (mtmp > 1.375E-9) {
    guard2 = true;
  } else {
    ixstart = 1;
    mtmp = adchf2[0];
    if (rtIsNaN(adchf2[0])) {
      ix = 2;
      exitg5 = false;
      while ((!exitg5) && (ix < 250)) {
        ixstart = ix;
        if (!rtIsNaN(adchf2[ix - 1])) {
          mtmp = adchf2[ix - 1];
          exitg5 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 249) {
      while (ixstart + 1 < 250) {
        if (adchf2[ixstart] > mtmp) {
          mtmp = adchf2[ixstart];
        }

        ixstart++;
      }
    }
      LOGE("CH2 MTMP1 [%f]",mtmp);
    if (mtmp > 1.375E-9) {
      guard2 = true;
    } else {
      ixstart = 1;
      mtmp = adchf3[0];
      if (rtIsNaN(adchf3[0])) {
        ix = 2;
        exitg4 = false;
        while ((!exitg4) && (ix < 250)) {
          ixstart = ix;
          if (!rtIsNaN(adchf3[ix - 1])) {
            mtmp = adchf3[ix - 1];
            exitg4 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < 249) {
        while (ixstart + 1 < 250) {
          if (adchf3[ixstart] > mtmp) {
            mtmp = adchf3[ixstart];
          }

          ixstart++;
        }
      }
        LOGE("CH3 MTMP1 [%f]",mtmp);
      if (mtmp > 1.375E-9) {
        guard2 = true;
      } else {
        b0 = false;
      }
    }
  }

  if (guard2) {
    b0 = true;
  }

  ixstart = 1;
  mtmp = adchf1[0];
  if (rtIsNaN(adchf1[0])) {
    ix = 2;
    exitg3 = false;
    while ((!exitg3) && (ix < 250)) {
      ixstart = ix;
      if (!rtIsNaN(adchf1[ix - 1])) {
        mtmp = adchf1[ix - 1];
        exitg3 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 249) {
    while (ixstart + 1 < 250) {
      if (adchf1[ixstart] > mtmp) {
        mtmp = adchf1[ixstart];
      }

      ixstart++;
    }
  }

  guard1 = false;
  if (mtmp > 1.25E-8) {
    guard1 = true;
  } else {
    ixstart = 1;
    mtmp = adchf2[0];
    if (rtIsNaN(adchf2[0])) {
      ix = 2;
      exitg2 = false;
      while ((!exitg2) && (ix < 250)) {
        ixstart = ix;
        if (!rtIsNaN(adchf2[ix - 1])) {
          mtmp = adchf2[ix - 1];
          exitg2 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 249) {
      while (ixstart + 1 < 250) {
        if (adchf2[ixstart] > mtmp) {
          mtmp = adchf2[ixstart];
        }

        ixstart++;
      }
    }

    if (mtmp > 1.25E-8) {
      guard1 = true;
    } else {
      ixstart = 1;
      mtmp = adchf3[0];
      if (rtIsNaN(adchf3[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 250)) {
          ixstart = ix;
          if (!rtIsNaN(adchf3[ix - 1])) {
            mtmp = adchf3[ix - 1];
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < 249) {
        while (ixstart + 1 < 250) {
          if (adchf3[ixstart] > mtmp) {
            mtmp = adchf3[ixstart];
          }

          ixstart++;
        }
      }

      if (mtmp > 1.25E-8) {
        guard1 = true;
      } else {
        b1 = false;
      }
    }
  }

  if (guard1) {
    b1 = true;
  }

  emxInit_real_T(&fch1f, 2);
  emxInit_real_T(&fch2f, 2);
  emxInit_real_T(&fch3f, 2);
  emxInit_real_T(&b_fch1f, 2);
//    LOGE("b0 = [%d] ; b1 = [%d]",b0, b1);
  if (b0 && (!b1)) {
    // load data:
    featureExtractionEOG(ch1f, fch1f);
    featureExtractionEOG(ch2f, fch2f);
    featureExtractionEOG(ch3f, fch3f);
    ixstart = b_fch1f->size[0] * b_fch1f->size[1];
    b_fch1f->size[0] = 1;
    b_fch1f->size[1] = (fch1f->size[1] + fch2f->size[1]) + fch3f->size[1];
    emxEnsureCapacity((emxArray__common *)b_fch1f, ixstart, (int)sizeof(double));
    ix = fch1f->size[1];
    for (ixstart = 0; ixstart < ix; ixstart++) {
      b_fch1f->data[b_fch1f->size[0] * ixstart] = fch1f->data[fch1f->size[0] *
        ixstart];
    }

    ix = fch2f->size[1];
    for (ixstart = 0; ixstart < ix; ixstart++) {
      b_fch1f->data[b_fch1f->size[0] * (ixstart + fch1f->size[1])] = fch2f->
        data[fch2f->size[0] * ixstart];
    }

    ix = fch3f->size[1];
    for (ixstart = 0; ixstart < ix; ixstart++) {
      b_fch1f->data[b_fch1f->size[0] * ((ixstart + fch1f->size[1]) + fch2f->
        size[1])] = fch3f->data[fch3f->size[0] * ixstart];
    }

    Y = knn(b_fch1f);
  }

  emxFree_real_T(&b_fch1f);
  emxFree_real_T(&fch3f);
  emxFree_real_T(&fch2f);
  emxFree_real_T(&fch1f);
  return Y;
}

//
// Arguments    : void
// Return Type  : void
//
void eog_knn_initialize()
{
  rt_InitInfAndNaN(8U);
//  omp_init_nest_lock(&emlrtNestLockGlobal);
}

//
// Arguments    : void
// Return Type  : void
//
void eog_knn_terminate()
{
//  omp_destroy_nest_lock(&emlrtNestLockGlobal);
}

//
// File trailer for eog_knn.cpp
//
// [EOF]
//
