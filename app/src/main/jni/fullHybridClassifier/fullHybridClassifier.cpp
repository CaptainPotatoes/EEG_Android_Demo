//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fullHybridClassifier.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 22-Mar-2017 11:10:49
//

// Include Files
#include "rt_nonfinite.h"
#include "fullHybridClassifier.h"
#include <cstdlib>

// Type Definitions
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
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

#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T

struct emxArray_creal_T
{
  creal_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_creal_T

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

// Variable Definitions
//omp_nest_lock_t emlrtNestLockGlobal;

// Function Declarations
static void b_abs(const creal_T x[2048], double y[2048]);
static void b_filter(const emxArray_real_T *x, const double zi[10],
                     emxArray_real_T *y);
static void b_filtfilt(const double x_in[250], double y_out[250]);
static void b_findLocalMaxima(const double yTemp[1025], emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void b_findpeaks(const double Yin[1025], emxArray_real_T *Ypk,
  emxArray_real_T *Xpk);
static boolean_T b_isequal(const emxArray_real_T *varargin_1, const
  emxArray_real_T *varargin_2);
static void b_log10(emxArray_real_T *x);
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
                    np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static double b_mod(double x, double y);
static void b_sign(emxArray_real_T *x);
static void b_sort(emxArray_real_T *x, emxArray_int32_T *idx);
static double b_sum(const double x[4]);
static void bluestein_setup(int nRows, emxArray_creal_T *wwc);
static void c_abs(const emxArray_creal_T *x, emxArray_real_T *y);
static void c_findLocalMaxima(emxArray_real_T *yTemp, emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void c_findPeaksSeparatedByMoreThanM(const emxArray_real_T *iPk,
  emxArray_real_T *idx);
static void c_findpeaks(const emxArray_real_T *Yin, emxArray_real_T *Ypk,
  emxArray_real_T *Xpk);
static void c_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx);
static void c_sum(const boolean_T x[16], double y[4]);
static void combinePeaks(const emxArray_real_T *iPk, const emxArray_real_T *iInf,
  emxArray_real_T *iPkOut);
static void d_findpeaks(const emxArray_real_T *Yin, emxArray_real_T *Ypk,
  emxArray_real_T *Xpk);
static double d_sum(const boolean_T x[4]);
static void diff(const emxArray_real_T *x, emxArray_real_T *y);
static void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib);
static void dobluesteinfft(const emxArray_real_T *x, int N2, int n1, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, const emxArray_real_T *
  sintabinv, emxArray_creal_T *y);
static void e_sum(const double x[16], double y[4]);
static void eegcfilt(emxArray_real_T *X, emxArray_real_T *Y);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
static void emxFree_creal_T(emxArray_creal_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
static void emxInit_boolean_T1(emxArray_boolean_T **pEmxArray, int numDimensions);
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
static void emxInit_creal_T1(emxArray_creal_T **pEmxArray, int numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static double eps(double x);
static void fSSVEPnew(const emxArray_real_T *fch1, const emxArray_real_T *fch2,
                      const emxArray_real_T *fch3, double Fs, emxArray_real_T *F);
static void featureExtractionEOG(const double samplesX[250], emxArray_real_T *F);
static void featureExtractionSSVEP(const emxArray_real_T *fch1, const
  emxArray_real_T *fch2, const emxArray_real_T *fch3, double Fs, double F[30]);
static void fft(const emxArray_real_T *x, emxArray_creal_T *y);
static void filter(double b[4], double a[4], const double x[268], const double
                   zi[3], double y[268]);
static void filtfilt(const double x_in[250], double y_out[250]);
static void findLocalMaxima(const double yTemp[250], emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static int findbin(double x, const emxArray_real_T *bin_edges);
static void findpeaks(const double Yin[250], emxArray_real_T *Ypk,
                      emxArray_real_T *Xpk);
static void getAllPeaks(const emxArray_real_T *y, emxArray_real_T *iPk,
  emxArray_real_T *iInf, emxArray_real_T *iInflect);
static void get_nfft_data(const emxArray_real_T *X, double Fs, double f[1025],
  double C[1025]);
static void hannWin(double x, emxArray_real_T *w);
static boolean_T isequal(const double varargin_1[4], const double varargin_2[4]);
static void keepAtMostNpPeaks(emxArray_real_T *idx, double Np);
static double knn(const double tsX[40]);
static double mean(const double x[4]);
static void merge(int idx[712], double x[712], int offset, int np, int nq, int
                  iwork[712], double xwork[712]);
static void merge_block(int idx[712], double x[712], int offset, int n, int
  preSortLevel, int iwork[712], double xwork[712]);
static void merge_pow2_block(int idx[712], double x[712], int offset);
static void orderPeaks(const emxArray_real_T *Y, const emxArray_real_T *iPk,
  emxArray_real_T *idx);
static void r2br_r2dit_trig_impl(const emxArray_creal_T *x, int xoffInit, int
  unsigned_nRows, const emxArray_real_T *costab, const emxArray_real_T *sintab,
  emxArray_creal_T *y);
static void removePeaksBelowMinPeakHeight(const emxArray_real_T *Y,
  emxArray_real_T *iPk, double Ph);
static void removePeaksBelowThreshold(const emxArray_real_T *Y, emxArray_real_T *
  iPk, double Th);
static double rt_hypotd_snf(double u0, double u1);
static double rt_roundd_snf(double u);
static void scaleAbs(emxArray_real_T *X, emxArray_real_T *Y);
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);
static void sort(double x[712], int idx[712]);
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
static void stft(const emxArray_real_T *x, double fs, emxArray_creal_T *s,
                 double f[1025], emxArray_real_T *t);
static void sum(const emxArray_real_T *x, emxArray_real_T *y);
static void treeClassifier(const double F[30], const emxArray_real_T *F2, double
  Y[7]);
static void welch_psd(const emxArray_real_T *signals, double fs, emxArray_real_T
                      *window, emxArray_real_T *CSM, emxArray_real_T
                      *frequencies);

// Function Definitions

//
// Arguments    : const creal_T x[2048]
//                double y[2048]
// Return Type  : void
//
static void b_abs(const creal_T x[2048], double y[2048])
{
  int k;
  for (k = 0; k < 2048; k++) {
    y[k] = rt_hypotd_snf(x[k].re, x[k].im);
  }
}

//
// Arguments    : const emxArray_real_T *x
//                const double zi[10]
//                emxArray_real_T *y
// Return Type  : void
//
static void b_filter(const emxArray_real_T *x, const double zi[10],
                     emxArray_real_T *y)
{
  unsigned int unnamed_idx_0;
  int j;
  double dbuffer[11];
  int k;
  double b_dbuffer;
  static const double dv14[11] = { 2.13961520749732E-5, 0.0,
    -0.000106980760374866, 0.0, 0.000213961520749732, 0.0, -0.000213961520749732,
    0.0, 0.000106980760374866, 0.0, -2.13961520749732E-5 };

  static const double dv15[11] = { 1.0, -8.77043379286888, 35.0068378010024,
    -83.7229808056309, 132.845833785487, -146.117834417428, 112.823239428442,
    -60.389449129414, 21.4471017127118, -4.56451967201817, 0.442209182399621 };

  unnamed_idx_0 = (unsigned int)x->size[0];
  j = y->size[0];
  y->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)y, j, (int)sizeof(double));
  memcpy(&dbuffer[1], &zi[0], 10U * sizeof(double));
  for (j = 0; j + 1 <= x->size[0]; j++) {
    for (k = 0; k < 10; k++) {
      dbuffer[k] = dbuffer[k + 1];
    }

    dbuffer[10] = 0.0;
    for (k = 0; k < 11; k++) {
      b_dbuffer = dbuffer[k] + x->data[j] * dv14[k];
      dbuffer[k] = b_dbuffer;
    }

    for (k = 0; k < 10; k++) {
      dbuffer[k + 1] -= dbuffer[0] * dv15[k + 1];
    }

    y->data[j] = dbuffer[0];
  }
}

//
// Arguments    : const double x_in[250]
//                double y_out[250]
// Return Type  : void
//
static void b_filtfilt(const double x_in[250], double y_out[250])
{
  double xtmp;
  double d1;
  int i;
  double y[268];
  double dv10[4];
  double dv11[4];
  double a[3];
  static const double dv12[4] = { 0.950971887923409, -2.85291566377023,
    2.85291566377023, -0.950971887923409 };

  static const double dv13[4] = { 1.0, -2.89947959461186, 2.803947977383,
    -0.904347531392409 };

  double b_y[268];
  static const double b_a[3] = { -0.95097188792826548, 1.9019437758560462,
    -0.95097188792780118 };

  double c_y[268];
  xtmp = 2.0 * x_in[0];
  d1 = 2.0 * x_in[249];
  for (i = 0; i < 9; i++) {
    y[i] = xtmp - x_in[9 - i];
  }

  memcpy(&y[9], &x_in[0], 250U * sizeof(double));
  for (i = 0; i < 9; i++) {
    y[i + 259] = d1 - x_in[248 - i];
  }

  for (i = 0; i < 4; i++) {
    dv10[i] = dv12[i];
    dv11[i] = dv13[i];
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 268U * sizeof(double));
  filter(dv10, dv11, b_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  for (i = 0; i < 4; i++) {
    dv10[i] = dv12[i];
    dv11[i] = dv13[i];
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&c_y[0], &y[0], 268U * sizeof(double));
  filter(dv10, dv11, c_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  memcpy(&y_out[0], &y[9], 250U * sizeof(double));
}

//
// Arguments    : const double yTemp[1025]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void b_findLocalMaxima(const double yTemp[1025], emxArray_real_T *iPk,
  emxArray_real_T *iInflect)
{
  double b_yTemp[1027];
  boolean_T yFinite[1027];
  int ii;
  boolean_T x[1026];
  emxArray_int32_T *b_ii;
  int idx;
  int i4;
  boolean_T exitg3;
  emxArray_int32_T *r9;
  boolean_T guard3 = false;
  emxArray_real_T *iTemp;
  emxArray_real_T *c_yTemp;
  emxArray_real_T *s;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r10;
  int nx;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *c_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  b_yTemp[0] = rtNaN;
  memcpy(&b_yTemp[1], &yTemp[0], 1025U * sizeof(double));
  b_yTemp[1026] = rtNaN;
  for (ii = 0; ii < 1027; ii++) {
    yFinite[ii] = !rtIsNaN(b_yTemp[ii]);
  }

  for (ii = 0; ii < 1026; ii++) {
    x[ii] = ((b_yTemp[ii] != b_yTemp[ii + 1]) && (yFinite[ii] || yFinite[ii + 1]));
  }

  emxInit_int32_T(&b_ii, 1);
  idx = 0;
  i4 = b_ii->size[0];
  b_ii->size[0] = 1026;
  emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 1027)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= 1026) {
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

  emxInit_int32_T(&r9, 1);
  i4 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  i4 = r9->size[0];
  r9->size[0] = 1 + b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r9, i4, (int)sizeof(int));
  r9->data[0] = 1;
  ii = b_ii->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    r9->data[i4 + 1] = b_ii->data[i4] + 1;
  }

  emxInit_real_T1(&iTemp, 1);
  i4 = iTemp->size[0];
  iTemp->size[0] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i4, (int)sizeof(double));
  ii = r9->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    iTemp->data[i4] = 1.0 + (double)(r9->data[i4] - 1);
  }

  emxFree_int32_T(&r9);
  emxInit_real_T1(&c_yTemp, 1);
  i4 = c_yTemp->size[0];
  c_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)c_yTemp, i4, (int)sizeof(double));
  ii = iTemp->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    c_yTemp->data[i4] = b_yTemp[(int)iTemp->data[i4] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_boolean_T(&b_x, 1);
  emxInit_real_T1(&r10, 1);
  diff(c_yTemp, s);
  b_sign(s);
  diff(s, r10);
  i4 = b_x->size[0];
  b_x->size[0] = r10->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i4, (int)sizeof(boolean_T));
  ii = r10->size[0];
  emxFree_real_T(&c_yTemp);
  for (i4 = 0; i4 < ii; i4++) {
    b_x->data[i4] = (r10->data[i4] < 0.0);
  }

  emxFree_real_T(&r10);
  nx = b_x->size[0];
  idx = 0;
  i4 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
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
      i4 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
    }
  } else {
    i4 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  }

  if (1.0 > (double)s->size[0] - 1.0) {
    ii = 0;
  } else {
    ii = (int)((double)s->size[0] - 1.0);
  }

  if (2 > s->size[0]) {
    i4 = 0;
  } else {
    i4 = 1;
  }

  idx = b_x->size[0];
  b_x->size[0] = ii;
  emxEnsureCapacity((emxArray__common *)b_x, idx, (int)sizeof(boolean_T));
  for (idx = 0; idx < ii; idx++) {
    b_x->data[idx] = (s->data[idx] != s->data[i4 + idx]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i4 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i4, (int)sizeof(int));
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
      i4 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i4, (int)sizeof(int));
    }
  } else {
    i4 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i4, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  i4 = iInflect->size[0];
  iInflect->size[0] = c_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i4, (int)sizeof(double));
  ii = c_ii->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    iInflect->data[i4] = iTemp->data[c_ii->data[i4]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  i4 = iPk->size[0];
  iPk->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i4, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    iPk->data[i4] = iTemp->data[b_ii->data[i4]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : const double Yin[1025]
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void b_findpeaks(const double Yin[1025], emxArray_real_T *Ypk,
  emxArray_real_T *Xpk)
{
  boolean_T x[1025];
  int k;
  emxArray_int32_T *ii;
  int idx;
  int cdiff;
  boolean_T exitg1;
  emxArray_real_T *iInfite;
  boolean_T guard1 = false;
  double yTemp[1025];
  emxArray_real_T *iPk;
  emxArray_real_T *iInflect;
  int ndbl;
  emxArray_real_T *base;
  double extremum;
  int apnd;
  emxArray_real_T *y;
  emxArray_real_T *b_idx;
  emxArray_real_T *c_idx;
  emxArray_real_T *d_idx;
  for (k = 0; k < 1025; k++) {
    x[k] = (rtIsInf(Yin[k]) && (Yin[k] > 0.0));
  }

  emxInit_int32_T(&ii, 1);
  idx = 0;
  k = ii->size[0];
  ii->size[0] = 1025;
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = 1;
  exitg1 = false;
  while ((!exitg1) && (cdiff < 1026)) {
    guard1 = false;
    if (x[cdiff - 1]) {
      idx++;
      ii->data[idx - 1] = cdiff;
      if (idx >= 1025) {
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

  memcpy(&yTemp[0], &Yin[0], 1025U * sizeof(double));
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
  b_findLocalMaxima(yTemp, iPk, iInflect);
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

  emxInit_real_T1(&b_idx, 1);
  k = b_idx->size[0];
  b_idx->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_idx, k, (int)sizeof(double));
  cdiff = y->size[1];
  for (k = 0; k < cdiff; k++) {
    b_idx->data[k] = y->data[y->size[0] * k];
  }

  emxFree_real_T(&y);
  if (b_idx->size[0] == 0) {
  } else {
    k = iInflect->size[0];
    iInflect->size[0] = b_idx->size[0];
    emxEnsureCapacity((emxArray__common *)iInflect, k, (int)sizeof(double));
    cdiff = b_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      iInflect->data[k] = Yin[(int)base->data[(int)b_idx->data[k] - 1] - 1];
    }

    emxInit_real_T1(&d_idx, 1);
    b_sort(iInflect, ii);
    k = d_idx->size[0];
    d_idx->size[0] = ii->size[0];
    emxEnsureCapacity((emxArray__common *)d_idx, k, (int)sizeof(double));
    cdiff = ii->size[0];
    for (k = 0; k < cdiff; k++) {
      d_idx->data[k] = b_idx->data[ii->data[k] - 1];
    }

    k = b_idx->size[0];
    b_idx->size[0] = d_idx->size[0];
    emxEnsureCapacity((emxArray__common *)b_idx, k, (int)sizeof(double));
    cdiff = d_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      b_idx->data[k] = d_idx->data[k];
    }

    emxFree_real_T(&d_idx);
  }

  emxFree_int32_T(&ii);
  emxFree_real_T(&iInflect);
  if (b_idx->size[0] > 1025) {
    emxInit_real_T1(&c_idx, 1);
    k = c_idx->size[0];
    c_idx->size[0] = 1025;
    emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
    for (k = 0; k < 1025; k++) {
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
    iPk->data[k] = base->data[(int)b_idx->data[k] - 1];
  }

  emxFree_real_T(&base);
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
// Arguments    : const emxArray_real_T *varargin_1
//                const emxArray_real_T *varargin_2
// Return Type  : boolean_T
//
static boolean_T b_isequal(const emxArray_real_T *varargin_1, const
  emxArray_real_T *varargin_2)
{
  boolean_T p;
  boolean_T b_p;
  int k;
  int exitg2;
  int i11;
  int i12;
  boolean_T exitg1;
  p = false;
  b_p = false;
  k = 0;
  do {
    exitg2 = 0;
    if (k < 2) {
      i11 = varargin_1->size[k];
      i12 = varargin_2->size[k];
      if (i11 != i12) {
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      b_p = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_p && (!(varargin_1->size[1] == 0)) && (!(varargin_2->size[1] == 0))) {
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= varargin_2->size[1] - 1)) {
      if (!(varargin_1->data[k] == varargin_2->data[k])) {
        b_p = false;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (!b_p) {
  } else {
    p = true;
  }

  return p;
}

//
// Arguments    : emxArray_real_T *x
// Return Type  : void
//
static void b_log10(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0] * x->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = std::log10(x->data[k]);
  }
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
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
                    np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
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
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
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
// Arguments    : double x
//                double y
// Return Type  : double
//
static double b_mod(double x, double y)
{
  double r;
  if (y == 0.0) {
    r = x;
  } else if (y == std::floor(y)) {
    r = x - std::floor(x / y) * y;
  } else {
    r = x / y;
    if (fabs(r - rt_roundd_snf(r)) <= 2.2204460492503131E-16 * fabs(r))
    {
      r = 0.0;
    } else {
      r = (r - std::floor(r)) * y;
    }
  }

  return r;
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
//                emxArray_int32_T *idx
// Return Type  : void
//
static void b_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  c_sort(x, dim, idx);
}

//
// Arguments    : const double x[4]
// Return Type  : double
//
static double b_sum(const double x[4])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 3; k++) {
    y += x[k + 1];
  }

  return y;
}

//
// Arguments    : int nRows
//                emxArray_creal_T *wwc
// Return Type  : void
//
static void bluestein_setup(int nRows, emxArray_creal_T *wwc)
{
  int nInt2m1;
  int idx;
  int rt;
  int nInt2;
  int k;
  int y;
  double nt_im;
  nInt2m1 = (nRows + nRows) - 1;
  idx = wwc->size[0];
  wwc->size[0] = nInt2m1;
  emxEnsureCapacity((emxArray__common *)wwc, idx, (int)sizeof(creal_T));
  idx = nRows;
  rt = 0;
  wwc->data[nRows - 1].re = 1.0;
  wwc->data[nRows - 1].im = 0.0;
  nInt2 = nRows << 1;
  for (k = 1; k < nRows; k++) {
    y = (k << 1) - 1;
    if (nInt2 - rt <= y) {
      rt += y - nInt2;
    } else {
      rt += y;
    }

    nt_im = -3.1415926535897931 * (double)rt / (double)nRows;
    wwc->data[idx - 2].re = std::cos(nt_im);
    wwc->data[idx - 2].im = -std::sin(nt_im);
    idx--;
  }

  idx = 0;
  for (k = nInt2m1 - 1; k >= nRows; k--) {
    wwc->data[k] = wwc->data[idx];
    idx++;
  }
}

//
// Arguments    : const emxArray_creal_T *x
//                emxArray_real_T *y
// Return Type  : void
//
static void c_abs(const emxArray_creal_T *x, emxArray_real_T *y)
{
  unsigned int uv0[2];
  int n;
  int k;
  for (n = 0; n < 2; n++) {
    uv0[n] = (unsigned int)x->size[n];
  }

  n = y->size[0] * y->size[1];
  y->size[0] = (int)uv0[0];
  y->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)y, n, (int)sizeof(double));
  n = x->size[0] * x->size[1];
  for (k = 0; k + 1 <= n; k++) {
    y->data[k] = rt_hypotd_snf(x->data[k].re, x->data[k].im);
  }
}

//
// Arguments    : emxArray_real_T *yTemp
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void c_findLocalMaxima(emxArray_real_T *yTemp, emxArray_real_T *iPk,
  emxArray_real_T *iInflect)
{
  emxArray_real_T *r12;
  int i7;
  int cdiff;
  int ndbl;
  int apnd;
  int absb;
  emxArray_real_T *y;
  emxArray_real_T *iTemp;
  emxArray_boolean_T *yFinite;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  boolean_T exitg3;
  boolean_T guard3 = false;
  emxArray_int32_T *r13;
  emxArray_real_T *b_iTemp;
  emxArray_real_T *b_yTemp;
  emxArray_real_T *s;
  emxArray_real_T *r14;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *b_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxInit_real_T1(&r12, 1);
  i7 = r12->size[0];
  r12->size[0] = 2 + yTemp->size[0];
  emxEnsureCapacity((emxArray__common *)r12, i7, (int)sizeof(double));
  r12->data[0] = rtNaN;
  cdiff = yTemp->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    r12->data[i7 + 1] = yTemp->data[i7];
  }

  r12->data[1 + yTemp->size[0]] = rtNaN;
  i7 = yTemp->size[0];
  yTemp->size[0] = r12->size[0];
  emxEnsureCapacity((emxArray__common *)yTemp, i7, (int)sizeof(double));
  cdiff = r12->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    yTemp->data[i7] = r12->data[i7];
  }

  emxFree_real_T(&r12);
  ndbl = (int)std::floor(((double)yTemp->size[0] - 1.0) + 0.5);
  apnd = ndbl + 1;
  cdiff = (ndbl - yTemp->size[0]) + 1;
  absb = yTemp->size[0];
  if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
    ndbl++;
    apnd = yTemp->size[0];
  } else if (cdiff > 0) {
    apnd = ndbl;
  } else {
    ndbl++;
  }

  emxInit_real_T(&y, 2);
  i7 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, i7, (int)sizeof(double));
  y->data[0] = 1.0;
  if (ndbl > 1) {
    y->data[ndbl - 1] = apnd;
    cdiff = (ndbl - 1) / 2;
    for (absb = 1; absb < cdiff; absb++) {
      y->data[absb] = 1.0 + (double)absb;
      y->data[(ndbl - absb) - 1] = apnd - absb;
    }

    if (cdiff << 1 == ndbl - 1) {
      y->data[cdiff] = (1.0 + (double)apnd) / 2.0;
    } else {
      y->data[cdiff] = 1.0 + (double)cdiff;
      y->data[cdiff + 1] = apnd - cdiff;
    }
  }

  emxInit_real_T1(&iTemp, 1);
  i7 = iTemp->size[0];
  iTemp->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)iTemp, i7, (int)sizeof(double));
  cdiff = y->size[1];
  for (i7 = 0; i7 < cdiff; i7++) {
    iTemp->data[i7] = y->data[y->size[0] * i7];
  }

  emxFree_real_T(&y);
  emxInit_boolean_T(&yFinite, 1);
  i7 = yFinite->size[0];
  yFinite->size[0] = yTemp->size[0];
  emxEnsureCapacity((emxArray__common *)yFinite, i7, (int)sizeof(boolean_T));
  cdiff = yTemp->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    yFinite->data[i7] = rtIsNaN(yTemp->data[i7]);
  }

  i7 = yFinite->size[0];
  emxEnsureCapacity((emxArray__common *)yFinite, i7, (int)sizeof(boolean_T));
  cdiff = yFinite->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    yFinite->data[i7] = !yFinite->data[i7];
  }

  emxInit_boolean_T(&x, 1);
  cdiff = yTemp->size[0] - 1;
  i7 = x->size[0];
  x->size[0] = cdiff;
  emxEnsureCapacity((emxArray__common *)x, i7, (int)sizeof(boolean_T));
  for (i7 = 0; i7 < cdiff; i7++) {
    x->data[i7] = ((yTemp->data[i7] != yTemp->data[1 + i7]) && (yFinite->data[i7]
      || yFinite->data[1 + i7]));
  }

  emxFree_boolean_T(&yFinite);
  emxInit_int32_T(&ii, 1);
  absb = x->size[0];
  ndbl = 0;
  i7 = ii->size[0];
  ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)ii, i7, (int)sizeof(int));
  cdiff = 1;
  exitg3 = false;
  while ((!exitg3) && (cdiff <= absb)) {
    guard3 = false;
    if (x->data[cdiff - 1]) {
      ndbl++;
      ii->data[ndbl - 1] = cdiff;
      if (ndbl >= absb) {
        exitg3 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3) {
      cdiff++;
    }
  }

  if (x->size[0] == 1) {
    if (ndbl == 0) {
      i7 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i7, (int)sizeof(int));
    }
  } else {
    i7 = ii->size[0];
    if (1 > ndbl) {
      ii->size[0] = 0;
    } else {
      ii->size[0] = ndbl;
    }

    emxEnsureCapacity((emxArray__common *)ii, i7, (int)sizeof(int));
  }

  emxInit_int32_T(&r13, 1);
  i7 = r13->size[0];
  r13->size[0] = 1 + ii->size[0];
  emxEnsureCapacity((emxArray__common *)r13, i7, (int)sizeof(int));
  r13->data[0] = 1;
  cdiff = ii->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    r13->data[i7 + 1] = ii->data[i7] + 1;
  }

  emxInit_real_T1(&b_iTemp, 1);
  i7 = b_iTemp->size[0];
  b_iTemp->size[0] = r13->size[0];
  emxEnsureCapacity((emxArray__common *)b_iTemp, i7, (int)sizeof(double));
  cdiff = r13->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    b_iTemp->data[i7] = iTemp->data[r13->data[i7] - 1];
  }

  emxFree_int32_T(&r13);
  i7 = iTemp->size[0];
  iTemp->size[0] = b_iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i7, (int)sizeof(double));
  cdiff = b_iTemp->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    iTemp->data[i7] = b_iTemp->data[i7];
  }

  emxFree_real_T(&b_iTemp);
  emxInit_real_T1(&b_yTemp, 1);
  i7 = b_yTemp->size[0];
  b_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)b_yTemp, i7, (int)sizeof(double));
  cdiff = iTemp->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    b_yTemp->data[i7] = yTemp->data[(int)iTemp->data[i7] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_real_T1(&r14, 1);
  diff(b_yTemp, s);
  b_sign(s);
  diff(s, r14);
  i7 = x->size[0];
  x->size[0] = r14->size[0];
  emxEnsureCapacity((emxArray__common *)x, i7, (int)sizeof(boolean_T));
  cdiff = r14->size[0];
  emxFree_real_T(&b_yTemp);
  for (i7 = 0; i7 < cdiff; i7++) {
    x->data[i7] = (r14->data[i7] < 0.0);
  }

  emxFree_real_T(&r14);
  absb = x->size[0];
  ndbl = 0;
  i7 = ii->size[0];
  ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)ii, i7, (int)sizeof(int));
  cdiff = 1;
  exitg2 = false;
  while ((!exitg2) && (cdiff <= absb)) {
    guard2 = false;
    if (x->data[cdiff - 1]) {
      ndbl++;
      ii->data[ndbl - 1] = cdiff;
      if (ndbl >= absb) {
        exitg2 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2) {
      cdiff++;
    }
  }

  if (x->size[0] == 1) {
    if (ndbl == 0) {
      i7 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i7, (int)sizeof(int));
    }
  } else {
    i7 = ii->size[0];
    if (1 > ndbl) {
      ii->size[0] = 0;
    } else {
      ii->size[0] = ndbl;
    }

    emxEnsureCapacity((emxArray__common *)ii, i7, (int)sizeof(int));
  }

  if (1.0 > (double)s->size[0] - 1.0) {
    cdiff = 0;
  } else {
    cdiff = (int)((double)s->size[0] - 1.0);
  }

  if (2 > s->size[0]) {
    i7 = 0;
  } else {
    i7 = 1;
  }

  absb = x->size[0];
  x->size[0] = cdiff;
  emxEnsureCapacity((emxArray__common *)x, absb, (int)sizeof(boolean_T));
  for (absb = 0; absb < cdiff; absb++) {
    x->data[absb] = (s->data[absb] != s->data[i7 + absb]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&b_ii, 1);
  absb = x->size[0];
  ndbl = 0;
  i7 = b_ii->size[0];
  b_ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i7, (int)sizeof(int));
  cdiff = 1;
  exitg1 = false;
  while ((!exitg1) && (cdiff <= absb)) {
    guard1 = false;
    if (x->data[cdiff - 1]) {
      ndbl++;
      b_ii->data[ndbl - 1] = cdiff;
      if (ndbl >= absb) {
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

  if (x->size[0] == 1) {
    if (ndbl == 0) {
      i7 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i7, (int)sizeof(int));
    }
  } else {
    i7 = b_ii->size[0];
    if (1 > ndbl) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = ndbl;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i7, (int)sizeof(int));
  }

  emxFree_boolean_T(&x);
  i7 = iInflect->size[0];
  iInflect->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i7, (int)sizeof(double));
  cdiff = b_ii->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    iInflect->data[i7] = iTemp->data[b_ii->data[i7]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  i7 = iPk->size[0];
  iPk->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i7, (int)sizeof(double));
  cdiff = ii->size[0];
  for (i7 = 0; i7 < cdiff; i7++) {
    iPk->data[i7] = iTemp->data[ii->data[i7]] - 1.0;
  }

  emxFree_int32_T(&ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : const emxArray_real_T *iPk
//                emxArray_real_T *idx
// Return Type  : void
//
static void c_findPeaksSeparatedByMoreThanM(const emxArray_real_T *iPk,
  emxArray_real_T *idx)
{
  int ndbl;
  int apnd;
  int cdiff;
  int absb;
  emxArray_real_T *y;
  if (iPk->size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)iPk->size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - iPk->size[0]) + 1;
    absb = iPk->size[0];
    if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
      ndbl++;
      apnd = iPk->size[0];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }
  }

  emxInit_real_T(&y, 2);
  absb = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, absb, (int)sizeof(double));
  if (ndbl > 0) {
    y->data[0] = 1.0;
    if (ndbl > 1) {
      y->data[ndbl - 1] = apnd;
      cdiff = (ndbl - 1) / 2;
      for (absb = 1; absb < cdiff; absb++) {
        y->data[absb] = 1.0 + (double)absb;
        y->data[(ndbl - absb) - 1] = apnd - absb;
      }

      if (cdiff << 1 == ndbl - 1) {
        y->data[cdiff] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[cdiff] = 1.0 + (double)cdiff;
        y->data[cdiff + 1] = apnd - cdiff;
      }
    }
  }

  absb = idx->size[0];
  idx->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)idx, absb, (int)sizeof(double));
  cdiff = y->size[1];
  for (absb = 0; absb < cdiff; absb++) {
    idx->data[absb] = y->data[y->size[0] * absb];
  }

  emxFree_real_T(&y);
}

//
// Arguments    : const emxArray_real_T *Yin
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void c_findpeaks(const emxArray_real_T *Yin, emxArray_real_T *Ypk,
  emxArray_real_T *Xpk)
{
  int nm1d2;
  int cdiff;
  int ndbl;
  int apnd;
  emxArray_real_T *y;
  emxArray_real_T *x;
  emxArray_real_T *iPk;
  emxArray_real_T *idx;
  emxArray_real_T *iInfite;
  emxArray_real_T *b_iPk;
  int b_Yin[1];
  emxArray_real_T c_Yin;
  int d_Yin[1];
  int e_Yin[1];
  int f_Yin[1];
  nm1d2 = Yin->size[1];
  if (nm1d2 < 1) {
    cdiff = 0;
    apnd = Yin->size[1];
  } else {
    nm1d2 = Yin->size[1];
    ndbl = (int)std::floor(((double)nm1d2 - 1.0) + 0.5);
    apnd = ndbl + 1;
    nm1d2 = Yin->size[1];
    cdiff = (ndbl - nm1d2) + 1;
    nm1d2 = Yin->size[1];
    if (1 >= nm1d2) {
      nm1d2 = 1;
    }

    if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)nm1d2) {
      ndbl++;
      apnd = Yin->size[1];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }

    if (ndbl >= 0) {
      cdiff = ndbl;
    } else {
      cdiff = 0;
    }
  }

  emxInit_real_T(&y, 2);
  ndbl = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = cdiff;
  emxEnsureCapacity((emxArray__common *)y, ndbl, (int)sizeof(double));
  if (cdiff > 0) {
    y->data[0] = 1.0;
    if (cdiff > 1) {
      y->data[cdiff - 1] = apnd;
      nm1d2 = (cdiff - 1) / 2;
      for (ndbl = 1; ndbl < nm1d2; ndbl++) {
        y->data[ndbl] = 1.0 + (double)ndbl;
        y->data[(cdiff - ndbl) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == cdiff - 1) {
        y->data[nm1d2] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[nm1d2] = 1.0 + (double)nm1d2;
        y->data[nm1d2 + 1] = apnd - nm1d2;
      }
    }
  }

  emxInit_real_T1(&x, 1);
  ndbl = x->size[0];
  x->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)x, ndbl, (int)sizeof(double));
  nm1d2 = y->size[1];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    x->data[ndbl] = y->data[y->size[0] * ndbl];
  }

  emxFree_real_T(&y);
  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&idx, 1);
  emxInit_real_T1(&iInfite, 1);
  emxInit_real_T1(&b_iPk, 1);
  b_Yin[0] = Yin->size[1];
  c_Yin = *Yin;
  c_Yin.size = (int *)&b_Yin;
  c_Yin.numDimensions = 1;
  getAllPeaks(&c_Yin, iPk, iInfite, idx);
  d_Yin[0] = Yin->size[1];
  c_Yin = *Yin;
  c_Yin.size = (int *)&d_Yin;
  c_Yin.numDimensions = 1;
  removePeaksBelowMinPeakHeight(&c_Yin, iPk, rtMinusInf);
  e_Yin[0] = Yin->size[1];
  c_Yin = *Yin;
  c_Yin.size = (int *)&e_Yin;
  c_Yin.numDimensions = 1;
  removePeaksBelowThreshold(&c_Yin, iPk, 0.0);
  combinePeaks(iPk, iInfite, b_iPk);
  c_findPeaksSeparatedByMoreThanM(b_iPk, idx);
  f_Yin[0] = Yin->size[1];
  c_Yin = *Yin;
  c_Yin.size = (int *)&f_Yin;
  c_Yin.numDimensions = 1;
  orderPeaks(&c_Yin, b_iPk, idx);
  nm1d2 = Yin->size[1];
  keepAtMostNpPeaks(idx, (double)nm1d2);
  ndbl = iPk->size[0];
  iPk->size[0] = idx->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, ndbl, (int)sizeof(double));
  nm1d2 = idx->size[0];
  emxFree_real_T(&iInfite);
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    iPk->data[ndbl] = b_iPk->data[(int)idx->data[ndbl] - 1];
  }

  emxFree_real_T(&b_iPk);
  emxFree_real_T(&idx);
  ndbl = Ypk->size[0] * Ypk->size[1];
  Ypk->size[0] = 1;
  Ypk->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Ypk, ndbl, (int)sizeof(double));
  nm1d2 = iPk->size[0];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    Ypk->data[Ypk->size[0] * ndbl] = Yin->data[(int)iPk->data[ndbl] - 1];
  }

  ndbl = Xpk->size[0] * Xpk->size[1];
  Xpk->size[0] = 1;
  Xpk->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Xpk, ndbl, (int)sizeof(double));
  nm1d2 = iPk->size[0];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    Xpk->data[Xpk->size[0] * ndbl] = x->data[(int)iPk->data[ndbl] - 1];
  }

  emxFree_real_T(&x);
  emxFree_real_T(&iPk);
}

//
// Arguments    : emxArray_real_T *x
//                int dim
//                emxArray_int32_T *idx
// Return Type  : void
//
static void c_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
  int i13;
  emxArray_real_T *vwork;
  int k;
  unsigned int unnamed_idx_0;
  int vstride;
  emxArray_int32_T *iidx;
  int j;
  if (dim <= 1) {
    i13 = x->size[0];
  } else {
    i13 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  k = vwork->size[0];
  vwork->size[0] = i13;
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
    for (k = 0; k + 1 <= i13; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, iidx);
    for (k = 0; k + 1 <= i13; k++) {
      x->data[j + k * vstride] = vwork->data[k];
      idx->data[j + k * vstride] = iidx->data[k];
    }
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

//
// Arguments    : const boolean_T x[16]
//                double y[4]
// Return Type  : void
//
static void c_sum(const boolean_T x[16], double y[4])
{
  int j;
  double s;
  int k;
  for (j = 0; j < 4; j++) {
    s = x[j];
    for (k = 0; k < 3; k++) {
      s += (double)x[j + ((k + 1) << 2)];
    }

    y[j] = s;
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
// Arguments    : const emxArray_real_T *Yin
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void d_findpeaks(const emxArray_real_T *Yin, emxArray_real_T *Ypk,
  emxArray_real_T *Xpk)
{
  boolean_T yIsRow;
  int nm1d2;
  int cdiff;
  int ndbl;
  int apnd;
  emxArray_real_T *y;
  emxArray_real_T *x;
  emxArray_real_T *iPk;
  emxArray_real_T *idx;
  emxArray_real_T *b_iPk;
  emxArray_real_T *b_Ypk;
  int b_Yin[1];
  emxArray_real_T c_Yin;
  int d_Yin[1];
  int e_Yin[1];
  int f_Yin[1];
  emxArray_real_T *b_Xpk;
  yIsRow = (Yin->size[0] == 1);
  nm1d2 = Yin->size[0];
  if (nm1d2 < 1) {
    cdiff = 0;
    apnd = Yin->size[0];
  } else {
    nm1d2 = Yin->size[0];
    ndbl = (int)std::floor(((double)nm1d2 - 1.0) + 0.5);
    apnd = ndbl + 1;
    nm1d2 = Yin->size[0];
    cdiff = (ndbl - nm1d2) + 1;
    nm1d2 = Yin->size[0];
    if (1 >= nm1d2) {
      nm1d2 = 1;
    }

    if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)nm1d2) {
      ndbl++;
      apnd = Yin->size[0];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }

    if (ndbl >= 0) {
      cdiff = ndbl;
    } else {
      cdiff = 0;
    }
  }

  emxInit_real_T(&y, 2);
  ndbl = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = cdiff;
  emxEnsureCapacity((emxArray__common *)y, ndbl, (int)sizeof(double));
  if (cdiff > 0) {
    y->data[0] = 1.0;
    if (cdiff > 1) {
      y->data[cdiff - 1] = apnd;
      nm1d2 = (cdiff - 1) / 2;
      for (ndbl = 1; ndbl < nm1d2; ndbl++) {
        y->data[ndbl] = 1.0 + (double)ndbl;
        y->data[(cdiff - ndbl) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == cdiff - 1) {
        y->data[nm1d2] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[nm1d2] = 1.0 + (double)nm1d2;
        y->data[nm1d2 + 1] = apnd - nm1d2;
      }
    }
  }

  emxInit_real_T1(&x, 1);
  ndbl = x->size[0];
  x->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)x, ndbl, (int)sizeof(double));
  nm1d2 = y->size[1];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    x->data[ndbl] = y->data[y->size[0] * ndbl];
  }

  emxFree_real_T(&y);
  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&idx, 1);
  emxInit_real_T1(&b_iPk, 1);
  emxInit_real_T1(&b_Ypk, 1);
  b_Yin[0] = Yin->size[0];
  c_Yin = *Yin;
  c_Yin.size = (int *)&b_Yin;
  c_Yin.numDimensions = 1;
  getAllPeaks(&c_Yin, iPk, b_Ypk, idx);
  d_Yin[0] = Yin->size[0];
  c_Yin = *Yin;
  c_Yin.size = (int *)&d_Yin;
  c_Yin.numDimensions = 1;
  removePeaksBelowMinPeakHeight(&c_Yin, iPk, rtMinusInf);
  e_Yin[0] = Yin->size[0];
  c_Yin = *Yin;
  c_Yin.size = (int *)&e_Yin;
  c_Yin.numDimensions = 1;
  removePeaksBelowThreshold(&c_Yin, iPk, 0.0);
  combinePeaks(iPk, b_Ypk, b_iPk);
  c_findPeaksSeparatedByMoreThanM(b_iPk, idx);
  f_Yin[0] = Yin->size[0];
  c_Yin = *Yin;
  c_Yin.size = (int *)&f_Yin;
  c_Yin.numDimensions = 1;
  orderPeaks(&c_Yin, b_iPk, idx);
  nm1d2 = Yin->size[0];
  keepAtMostNpPeaks(idx, (double)nm1d2);
  ndbl = iPk->size[0];
  iPk->size[0] = idx->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, ndbl, (int)sizeof(double));
  nm1d2 = idx->size[0];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    iPk->data[ndbl] = b_iPk->data[(int)idx->data[ndbl] - 1];
  }

  emxFree_real_T(&b_iPk);
  ndbl = b_Ypk->size[0];
  b_Ypk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)b_Ypk, ndbl, (int)sizeof(double));
  nm1d2 = iPk->size[0];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    b_Ypk->data[ndbl] = Yin->data[(int)iPk->data[ndbl] - 1];
  }

  emxInit_real_T1(&b_Xpk, 1);
  ndbl = b_Xpk->size[0];
  b_Xpk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)b_Xpk, ndbl, (int)sizeof(double));
  nm1d2 = iPk->size[0];
  for (ndbl = 0; ndbl < nm1d2; ndbl++) {
    b_Xpk->data[ndbl] = x->data[(int)iPk->data[ndbl] - 1];
  }

  emxFree_real_T(&x);
  emxFree_real_T(&iPk);
  if (yIsRow) {
    ndbl = Ypk->size[0] * Ypk->size[1];
    Ypk->size[0] = 1;
    Ypk->size[1] = b_Ypk->size[0];
    emxEnsureCapacity((emxArray__common *)Ypk, ndbl, (int)sizeof(double));
    nm1d2 = b_Ypk->size[0];
    for (ndbl = 0; ndbl < nm1d2; ndbl++) {
      Ypk->data[Ypk->size[0] * ndbl] = b_Ypk->data[ndbl];
    }
  } else {
    ndbl = Ypk->size[0] * Ypk->size[1];
    Ypk->size[0] = idx->size[0];
    Ypk->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)Ypk, ndbl, (int)sizeof(double));
    nm1d2 = idx->size[0];
    for (ndbl = 0; ndbl < nm1d2; ndbl++) {
      Ypk->data[ndbl] = b_Ypk->data[ndbl];
    }
  }

  emxFree_real_T(&b_Ypk);
  if (yIsRow) {
    ndbl = Xpk->size[0] * Xpk->size[1];
    Xpk->size[0] = 1;
    Xpk->size[1] = b_Xpk->size[0];
    emxEnsureCapacity((emxArray__common *)Xpk, ndbl, (int)sizeof(double));
    nm1d2 = b_Xpk->size[0];
    for (ndbl = 0; ndbl < nm1d2; ndbl++) {
      Xpk->data[Xpk->size[0] * ndbl] = b_Xpk->data[ndbl];
    }
  } else {
    ndbl = Xpk->size[0] * Xpk->size[1];
    Xpk->size[0] = idx->size[0];
    Xpk->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)Xpk, ndbl, (int)sizeof(double));
    nm1d2 = idx->size[0];
    for (ndbl = 0; ndbl < nm1d2; ndbl++) {
      Xpk->data[ndbl] = b_Xpk->data[ndbl];
    }
  }

  emxFree_real_T(&b_Xpk);
  emxFree_real_T(&idx);
}

//
// Arguments    : const boolean_T x[4]
// Return Type  : double
//
static double d_sum(const boolean_T x[4])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 3; k++) {
    y += (double)x[k + 1];
  }

  return y;
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_real_T *y
// Return Type  : void
//
static void diff(const emxArray_real_T *x, emxArray_real_T *y)
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
// Arguments    : const emxArray_real_T *x
//                int N2
//                int n1
//                const emxArray_real_T *costab
//                const emxArray_real_T *sintab
//                const emxArray_real_T *sintabinv
//                emxArray_creal_T *y
// Return Type  : void
//
static void dobluesteinfft(const emxArray_real_T *x, int N2, int n1, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, const emxArray_real_T *
  sintabinv, emxArray_creal_T *y)
{
  emxArray_creal_T *wwc;
  int minNrowsNx;
  int ix;
  int xidx;
  double r;
  double twid_im;
  emxArray_creal_T *fy;
  int istart;
  emxArray_creal_T *fv;
  int nRowsD2;
  int nRowsD4;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int j;
  double fv_re;
  double fv_im;
  int ihi;
  double wwc_im;
  double b_fv_re;
  emxInit_creal_T(&wwc, 1);
  bluestein_setup(n1, wwc);
  if (n1 <= x->size[0]) {
    minNrowsNx = n1;
  } else {
    minNrowsNx = x->size[0];
  }

  ix = y->size[0];
  y->size[0] = n1;
  emxEnsureCapacity((emxArray__common *)y, ix, (int)sizeof(creal_T));
  if (n1 > x->size[0]) {
    xidx = y->size[0];
    ix = y->size[0];
    y->size[0] = xidx;
    emxEnsureCapacity((emxArray__common *)y, ix, (int)sizeof(creal_T));
    for (ix = 0; ix < xidx; ix++) {
      y->data[ix].re = 0.0;
      y->data[ix].im = 0.0;
    }
  }

  xidx = 0;
  for (ix = 0; ix + 1 <= minNrowsNx; ix++) {
    r = wwc->data[(n1 + ix) - 1].re;
    twid_im = wwc->data[(n1 + ix) - 1].im;
    y->data[ix].re = r * x->data[xidx];
    y->data[ix].im = twid_im * -x->data[xidx];
    xidx++;
  }

  while (minNrowsNx + 1 <= n1) {
    y->data[minNrowsNx].re = 0.0;
    y->data[minNrowsNx].im = 0.0;
    minNrowsNx++;
  }

  emxInit_creal_T(&fy, 1);
  r2br_r2dit_trig_impl(y, 0, N2, costab, sintab, fy);
  if (wwc->size[0] <= N2) {
    istart = wwc->size[0];
  } else {
    istart = N2;
  }

  emxInit_creal_T(&fv, 1);
  nRowsD2 = N2 / 2;
  nRowsD4 = nRowsD2 / 2;
  ix = fv->size[0];
  fv->size[0] = N2;
  emxEnsureCapacity((emxArray__common *)fv, ix, (int)sizeof(creal_T));
  if (N2 > wwc->size[0]) {
    xidx = fv->size[0];
    ix = fv->size[0];
    fv->size[0] = xidx;
    emxEnsureCapacity((emxArray__common *)fv, ix, (int)sizeof(creal_T));
    for (ix = 0; ix < xidx; ix++) {
      fv->data[ix].re = 0.0;
      fv->data[ix].im = 0.0;
    }
  }

  ix = 0;
  minNrowsNx = 0;
  xidx = 0;
  for (i = 1; i < istart; i++) {
    fv->data[xidx] = wwc->data[ix];
    xidx = N2;
    tst = true;
    while (tst) {
      xidx >>= 1;
      minNrowsNx ^= xidx;
      tst = ((minNrowsNx & xidx) == 0);
    }

    xidx = minNrowsNx;
    ix++;
  }

  fv->data[xidx] = wwc->data[ix];
  if (N2 > 1) {
    for (i = 0; i <= N2 - 2; i += 2) {
      temp_re = fv->data[i + 1].re;
      temp_im = fv->data[i + 1].im;
      fv->data[i + 1].re = fv->data[i].re - fv->data[i + 1].re;
      fv->data[i + 1].im = fv->data[i].im - fv->data[i + 1].im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }
  }

  xidx = 2;
  minNrowsNx = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += minNrowsNx) {
      temp_re = fv->data[i + xidx].re;
      temp_im = fv->data[i + xidx].im;
      fv->data[i + xidx].re = fv->data[i].re - temp_re;
      fv->data[i + xidx].im = fv->data[i].im - temp_im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }

    istart = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      r = costab->data[j];
      twid_im = sintab->data[j];
      i = istart;
      ihi = istart + ix;
      while (i < ihi) {
        temp_re = r * fv->data[i + xidx].re - twid_im * fv->data[i + xidx].im;
        temp_im = r * fv->data[i + xidx].im + twid_im * fv->data[i + xidx].re;
        fv->data[i + xidx].re = fv->data[i].re - temp_re;
        fv->data[i + xidx].im = fv->data[i].im - temp_im;
        fv->data[i].re += temp_re;
        fv->data[i].im += temp_im;
        i += minNrowsNx;
      }

      istart++;
    }

    nRowsD4 /= 2;
    xidx = minNrowsNx;
    minNrowsNx <<= 1;
    ix -= xidx;
  }

  ix = fy->size[0];
  emxEnsureCapacity((emxArray__common *)fy, ix, (int)sizeof(creal_T));
  xidx = fy->size[0];
  for (ix = 0; ix < xidx; ix++) {
    r = fy->data[ix].re;
    twid_im = fy->data[ix].im;
    fv_re = fv->data[ix].re;
    fv_im = fv->data[ix].im;
    fy->data[ix].re = r * fv_re - twid_im * fv_im;
    fy->data[ix].im = r * fv_im + twid_im * fv_re;
  }

  if (fy->size[0] <= N2) {
    istart = fy->size[0];
  } else {
    istart = N2;
  }

  nRowsD2 = N2 / 2;
  nRowsD4 = nRowsD2 / 2;
  ix = fv->size[0];
  fv->size[0] = N2;
  emxEnsureCapacity((emxArray__common *)fv, ix, (int)sizeof(creal_T));
  if (N2 > fy->size[0]) {
    xidx = fv->size[0];
    ix = fv->size[0];
    fv->size[0] = xidx;
    emxEnsureCapacity((emxArray__common *)fv, ix, (int)sizeof(creal_T));
    for (ix = 0; ix < xidx; ix++) {
      fv->data[ix].re = 0.0;
      fv->data[ix].im = 0.0;
    }
  }

  ix = 0;
  minNrowsNx = 0;
  xidx = 0;
  for (i = 1; i < istart; i++) {
    fv->data[xidx] = fy->data[ix];
    xidx = N2;
    tst = true;
    while (tst) {
      xidx >>= 1;
      minNrowsNx ^= xidx;
      tst = ((minNrowsNx & xidx) == 0);
    }

    xidx = minNrowsNx;
    ix++;
  }

  fv->data[xidx] = fy->data[ix];
  emxFree_creal_T(&fy);
  if (N2 > 1) {
    for (i = 0; i <= N2 - 2; i += 2) {
      temp_re = fv->data[i + 1].re;
      temp_im = fv->data[i + 1].im;
      fv->data[i + 1].re = fv->data[i].re - fv->data[i + 1].re;
      fv->data[i + 1].im = fv->data[i].im - fv->data[i + 1].im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }
  }

  xidx = 2;
  minNrowsNx = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += minNrowsNx) {
      temp_re = fv->data[i + xidx].re;
      temp_im = fv->data[i + xidx].im;
      fv->data[i + xidx].re = fv->data[i].re - temp_re;
      fv->data[i + xidx].im = fv->data[i].im - temp_im;
      fv->data[i].re += temp_re;
      fv->data[i].im += temp_im;
    }

    istart = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      r = costab->data[j];
      twid_im = sintabinv->data[j];
      i = istart;
      ihi = istart + ix;
      while (i < ihi) {
        temp_re = r * fv->data[i + xidx].re - twid_im * fv->data[i + xidx].im;
        temp_im = r * fv->data[i + xidx].im + twid_im * fv->data[i + xidx].re;
        fv->data[i + xidx].re = fv->data[i].re - temp_re;
        fv->data[i + xidx].im = fv->data[i].im - temp_im;
        fv->data[i].re += temp_re;
        fv->data[i].im += temp_im;
        i += minNrowsNx;
      }

      istart++;
    }

    nRowsD4 /= 2;
    xidx = minNrowsNx;
    minNrowsNx <<= 1;
    ix -= xidx;
  }

  if (fv->size[0] > 1) {
    r = 1.0 / (double)fv->size[0];
    ix = fv->size[0];
    emxEnsureCapacity((emxArray__common *)fv, ix, (int)sizeof(creal_T));
    xidx = fv->size[0];
    for (ix = 0; ix < xidx; ix++) {
      fv->data[ix].re *= r;
      fv->data[ix].im *= r;
    }
  }

  xidx = 0;
  for (ix = n1 - 1; ix + 1 <= wwc->size[0]; ix++) {
    r = wwc->data[ix].re;
    fv_re = fv->data[ix].re;
    twid_im = wwc->data[ix].im;
    fv_im = fv->data[ix].im;
    temp_re = wwc->data[ix].re;
    temp_im = fv->data[ix].im;
    wwc_im = wwc->data[ix].im;
    b_fv_re = fv->data[ix].re;
    y->data[xidx].re = r * fv_re + twid_im * fv_im;
    y->data[xidx].im = temp_re * temp_im - wwc_im * b_fv_re;
    xidx++;
  }

  emxFree_creal_T(&fv);
  emxFree_creal_T(&wwc);
}

//
// Arguments    : const double x[16]
//                double y[4]
// Return Type  : void
//
static void e_sum(const double x[16], double y[4])
{
  int j;
  double s;
  int k;
  for (j = 0; j < 4; j++) {
    s = x[j];
    for (k = 0; k < 3; k++) {
      s += x[j + ((k + 1) << 2)];
    }

    y[j] = s;
  }
}

//
// EOGCFILT EEG filter for conversion to C.
//  Vectorize:
// Arguments    : emxArray_real_T *X
//                emxArray_real_T *Y
// Return Type  : void
//
static void eegcfilt(emxArray_real_T *X, emxArray_real_T *Y)
{
  int m;
  int i;
  emxArray_real_T *x;
  emxArray_real_T *y;
  double xtmp;
  double d2;
  int md2;
  double a[10];
  emxArray_real_T *b_y;
  static const double b_a[10] = { -2.1396152021655335E-5, -2.1396152489276133E-5,
    8.558460975207999E-5, 8.5584605288149449E-5, -0.00012837690837852629,
    -0.00012837691616921775, 8.5584610596008311E-5, 8.5584607376171939E-5,
    -2.1396151855180404E-5, -2.1396152098550849E-5 };

  emxArray_real_T *c_y;
  emxArray_int32_T *r7;
  m = X->size[0];
  i = X->size[0];
  X->size[0] = m;
  emxEnsureCapacity((emxArray__common *)X, i, (int)sizeof(double));

  //  Fs = 250, N = 5
  //  flim = [8 18], bandpass
  emxInit_real_T1(&x, 1);
  if (X->size[0] == 1) {
    i = x->size[0];
    x->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)x, i, (int)sizeof(double));
    x->data[0] = X->data[0];
  } else {
    i = x->size[0];
    x->size[0] = X->size[0];
    emxEnsureCapacity((emxArray__common *)x, i, (int)sizeof(double));
    m = X->size[0];
    for (i = 0; i < m; i++) {
      x->data[i] = X->data[i];
    }
  }

  if (x->size[0] == 0) {
    i = Y->size[0] * Y->size[1];
    Y->size[0] = 0;
    Y->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)Y, i, (int)sizeof(double));
  } else {
    emxInit_real_T1(&y, 1);
    xtmp = 2.0 * x->data[0];
    d2 = 2.0 * x->data[x->size[0] - 1];
    md2 = x->size[0] - 2;
    i = y->size[0];
    y->size[0] = 60 + x->size[0];
    emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
    for (i = 0; i < 30; i++) {
      y->data[i] = xtmp - x->data[30 - i];
    }

    m = x->size[0];
    for (i = 0; i < m; i++) {
      y->data[i + 30] = x->data[i];
    }

    for (i = 0; i < 30; i++) {
      y->data[(i + x->size[0]) + 30] = d2 - x->data[md2 - i];
    }

    xtmp = y->data[0];
    for (i = 0; i < 10; i++) {
      a[i] = b_a[i] * xtmp;
    }

    emxInit_real_T1(&b_y, 1);
    i = b_y->size[0];
    b_y->size[0] = y->size[0];
    emxEnsureCapacity((emxArray__common *)b_y, i, (int)sizeof(double));
    m = y->size[0];
    for (i = 0; i < m; i++) {
      b_y->data[i] = y->data[i];
    }

    b_filter(b_y, a, y);
    m = y->size[0];
    md2 = y->size[0] >> 1;
    i = 1;
    emxFree_real_T(&b_y);
    while (i <= md2) {
      xtmp = y->data[i - 1];
      y->data[i - 1] = y->data[m - i];
      y->data[m - i] = xtmp;
      i++;
    }

    xtmp = y->data[0];
    for (i = 0; i < 10; i++) {
      a[i] = b_a[i] * xtmp;
    }

    emxInit_real_T1(&c_y, 1);
    i = c_y->size[0];
    c_y->size[0] = y->size[0];
    emxEnsureCapacity((emxArray__common *)c_y, i, (int)sizeof(double));
    m = y->size[0];
    for (i = 0; i < m; i++) {
      c_y->data[i] = y->data[i];
    }

    b_filter(c_y, a, y);
    m = y->size[0];
    md2 = y->size[0] >> 1;
    i = 1;
    emxFree_real_T(&c_y);
    while (i <= md2) {
      xtmp = y->data[i - 1];
      y->data[i - 1] = y->data[m - i];
      y->data[m - i] = xtmp;
      i++;
    }

    if (X->size[0] == 1) {
      m = x->size[0] - 1;
      i = Y->size[0] * Y->size[1];
      Y->size[0] = 1;
      Y->size[1] = m + 1;
      emxEnsureCapacity((emxArray__common *)Y, i, (int)sizeof(double));
      for (i = 0; i <= m; i++) {
        Y->data[Y->size[0] * i] = y->data[30 + i];
      }
    } else {
      emxInit_int32_T(&r7, 1);
      m = x->size[0];
      i = Y->size[0] * Y->size[1];
      Y->size[0] = m;
      Y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)Y, i, (int)sizeof(double));
      i = r7->size[0];
      r7->size[0] = m;
      emxEnsureCapacity((emxArray__common *)r7, i, (int)sizeof(int));
      for (i = 0; i < m; i++) {
        r7->data[i] = 31 + i;
      }

      for (i = 0; i < m; i++) {
        Y->data[i] = y->data[r7->data[i] - 1];
      }

      emxFree_int32_T(&r7);
    }

    emxFree_real_T(&y);
  }

  emxFree_real_T(&x);
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
// Arguments    : emxArray_creal_T **pEmxArray
// Return Type  : void
//
static void emxFree_creal_T(emxArray_creal_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_creal_T *)NULL) {
    if (((*pEmxArray)->data != (creal_T *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_creal_T *)NULL;
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
// Arguments    : emxArray_boolean_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_boolean_T1(emxArray_boolean_T **pEmxArray, int numDimensions)
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
// Arguments    : emxArray_creal_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions)
{
  emxArray_creal_T *emxArray;
  int i;
  *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
  emxArray = *pEmxArray;
  emxArray->data = (creal_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_creal_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_creal_T1(emxArray_creal_T **pEmxArray, int numDimensions)
{
  emxArray_creal_T *emxArray;
  int i;
  *pEmxArray = (emxArray_creal_T *)malloc(sizeof(emxArray_creal_T));
  emxArray = *pEmxArray;
  emxArray->data = (creal_T *)NULL;
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
// Arguments    : emxArray_int32_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions)
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
// Feature Extraction Function for Tri-channel SSVEP Feature Extraction:
//  ----- INPUTS -----
//  fch1, fch2, fch3: Tri-channel SSVEP Samples of certain window size
//  MUST BE A VECTOR!
//  MUST BE FILTERED!!
//  Fs - Sampling Rate.
//  Using 250-sample windows, feature extraction is obtained using FFT and
//  PSD
// Arguments    : const emxArray_real_T *fch1
//                const emxArray_real_T *fch2
//                const emxArray_real_T *fch3
//                double Fs
//                emxArray_real_T *F
// Return Type  : void
//
static void fSSVEPnew(const emxArray_real_T *fch1, const emxArray_real_T *fch2,
                      const emxArray_real_T *fch3, double Fs, emxArray_real_T *F)
{
  double threshFFT[8];
  int i8;
  double threshPSD[8];
  double threshSTFT[8];
  emxArray_real_T *fchw;
  int fch1_idx_0;
  int b_fch1_idx_0;
  int loop_ub;
  int c_fch1_idx_0;
  emxArray_int32_T *r15;
  emxArray_real_T *b_fch1;
  emxArray_real_T *b_fch2;
  emxArray_real_T *b_fch3;
  double FFT[4100];
  double fft_sel_loc[16];
  double fft_sel_pks[16];
  emxArray_real_T *PSD;
  double y;
  int iv0[2];
  emxArray_boolean_T *selectPSD;
  double psd_sel_loc[16];
  double psd_sel_pks[16];
  double stft_sel_loc[4];
  double stft_sel_pks[4];
  int i;
  emxArray_real_T *hW;
  emxArray_real_T *fPSD;
  emxArray_real_T *SummedRows;
  emxArray_creal_T *S1;
  emxArray_real_T *T;
  emxArray_creal_T *S2;
  emxArray_creal_T *S3;
  emxArray_real_T *psd_sel_L;
  emxArray_boolean_T *r16;
  emxArray_boolean_T *r17;
  emxArray_int32_T *r18;
  emxArray_int32_T *r19;
  emxArray_real_T *x;
  emxArray_real_T *b_x;
  emxArray_real_T *c_x;
  emxArray_real_T *b_SummedRows;
  emxArray_real_T *r20;
  emxArray_creal_T *b_S3;
  emxArray_creal_T *b_S2;
  emxArray_creal_T *b_S1;
  emxArray_real_T *b_fchw;
  emxArray_real_T *c_fchw;
  emxArray_real_T *d_fchw;
  emxArray_real_T *e_fchw;
  emxArray_real_T *f_fchw;
  emxArray_real_T *b_FFT;
  emxArray_real_T *b_PSD;
  emxArray_real_T *c_FFT;
  emxArray_real_T *c_PSD;
  emxArray_real_T *c_SummedRows;
  emxArray_real_T *g_fchw;
  emxArray_real_T *h_fchw;
  emxArray_real_T *d_FFT;
  emxArray_real_T *d_PSD;
  emxArray_real_T *e_FFT;
  emxArray_real_T *e_PSD;
  emxArray_real_T *f_PSD;
  emxArray_int32_T *r21;
  emxArray_int32_T *r22;
  emxArray_int32_T *r23;
  emxArray_int32_T *r24;
  emxArray_real_T *g_PSD;
  emxArray_int32_T *r25;
  emxArray_int32_T *r26;
  boolean_T guard2 = false;
  boolean_T b_select[4100];
  int ch;
  double f[1025];
  double dv19[1025];
  double b_F[1025];
  boolean_T guard1 = false;
  double ft_ch[72];
  double b_ft_ch[64];
  double c_ft_ch[72];
  double d_ft_ch[64];

  // %%%% - Thresholds: - %%%%%
  // ----FFT----%
  // -% windows around certain target frequencies
  for (i8 = 0; i8 < 2; i8++) {
    threshFFT[i8 << 2] = 9.5 + 1.1300000000000008 * (double)i8;
    threshFFT[1 + (i8 << 2)] = 11.9 + 0.79999999999999893 * (double)i8;
    threshFFT[2 + (i8 << 2)] = 14.6 + 0.90000000000000036 * (double)i8;
    threshFFT[3 + (i8 << 2)] = 16.1 + 0.69999999999999929 * (double)i8;
  }

  // ---FOR >= 500 DP ---%
  //  threshFFTL = zeros(4,2);
  //  threshFFTL(1,:) = [9.76  10.14];%-% windows around certain target frequencies 
  //  threshFFTL(2,:) = [12.2  12.6];
  //  threshFFTL(3,:) = [14.76 15.27];
  //  threshFFTL(4,:) = [16.1 16.80];
  // Also use wLFFT
  // ----PSD----%
  for (i8 = 0; i8 < 2; i8++) {
    threshPSD[i8 << 2] = 9.0 + 2.0 * (double)i8;
    threshPSD[1 + (i8 << 2)] = 11.0 + 2.0 * (double)i8;
    threshPSD[2 + (i8 << 2)] = 14.0 + 1.5 * (double)i8;
    threshPSD[3 + (i8 << 2)] = 15.5 + 2.0 * (double)i8;
  }

  // ---FOR >= 500 DP ---%
  //  threshPSDL = zeros(4,2);
  //  threshPSDL(1,:) = [9.79 10.25];
  //  threshPSDL(2,:) = [12.2 12.8];
  //  threshPSDL(3,:) = [14.75 15.5];
  //  threshPSDL(4,:) = [16 17];
  // ---STFT---%
  // ---FOR >= 500 DP ---%
  for (i8 = 0; i8 < 2; i8++) {
    threshSTFT[i8 << 2] = 9.6 + 0.66000000000000014 * (double)i8;
    threshSTFT[1 + (i8 << 2)] = 12.2 + 0.51000000000000156 * (double)i8;
    threshSTFT[2 + (i8 << 2)] = 14.75 + 0.64000000000000057 * (double)i8;
    threshSTFT[3 + (i8 << 2)] = 16.11 + 0.870000000000001 * (double)i8;
  }

  emxInit_real_T(&fchw, 2);

  // ---FOR >= 1250 DP ---%
  //  threshSTFT5 = zeros(4,2);
  //  threshSTFT5(1,:) = [9.76  10.14];
  //  threshSTFT5(2,:) = [12.3  12.58];
  //  threshSTFT5(3,:) = [14.85 15.25];
  //  threshSTFT5(4,:) = [16.33 16.87];
  // ----PREALLOCATE----%
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  i8 = fchw->size[0] * fchw->size[1];
  fchw->size[0] = 3;
  emxEnsureCapacity((emxArray__common *)fchw, i8, (int)sizeof(double));
  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  i8 = fchw->size[0] * fchw->size[1];
  fchw->size[1] = b_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)fchw, i8, (int)sizeof(double));
  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  loop_ub = 3 * b_fch1_idx_0;
  for (i8 = 0; i8 < loop_ub; i8++) {
    fchw->data[i8] = 0.0;
  }

  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > b_fch1_idx_0) {
    b_fch1_idx_0 = 0;
  } else {
    b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  c_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > c_fch1_idx_0) {
    c_fch1_idx_0 = 0;
  } else {
    c_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  emxInit_int32_T(&r15, 1);
  i8 = r15->size[0];
  r15->size[0] = c_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)r15, i8, (int)sizeof(int));
  for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
    r15->data[i8] = i8;
  }

  emxInit_real_T1(&b_fch1, 1);
  i8 = b_fch1->size[0];
  b_fch1->size[0] = b_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch1, i8, (int)sizeof(double));
  for (i8 = 0; i8 < b_fch1_idx_0; i8++) {
    b_fch1->data[i8] = fch1->data[i8];
  }

  c_fch1_idx_0 = r15->size[0];
  for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
    fchw->data[fchw->size[0] * r15->data[i8]] = b_fch1->data[i8];
  }

  emxFree_real_T(&b_fch1);
  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > b_fch1_idx_0) {
    b_fch1_idx_0 = 0;
  } else {
    b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  c_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > c_fch1_idx_0) {
    c_fch1_idx_0 = 0;
  } else {
    c_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  i8 = r15->size[0];
  r15->size[0] = c_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)r15, i8, (int)sizeof(int));
  for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
    r15->data[i8] = i8;
  }

  emxInit_real_T1(&b_fch2, 1);
  i8 = b_fch2->size[0];
  b_fch2->size[0] = b_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch2, i8, (int)sizeof(double));
  for (i8 = 0; i8 < b_fch1_idx_0; i8++) {
    b_fch2->data[i8] = fch2->data[i8];
  }

  c_fch1_idx_0 = r15->size[0];
  for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
    fchw->data[1 + fchw->size[0] * r15->data[i8]] = b_fch2->data[i8];
  }

  emxFree_real_T(&b_fch2);
  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > b_fch1_idx_0) {
    b_fch1_idx_0 = 0;
  } else {
    b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  c_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > c_fch1_idx_0) {
    c_fch1_idx_0 = 0;
  } else {
    c_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  i8 = r15->size[0];
  r15->size[0] = c_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)r15, i8, (int)sizeof(int));
  for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
    r15->data[i8] = i8;
  }

  emxInit_real_T1(&b_fch3, 1);
  i8 = b_fch3->size[0];
  b_fch3->size[0] = b_fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch3, i8, (int)sizeof(double));
  for (i8 = 0; i8 < b_fch1_idx_0; i8++) {
    b_fch3->data[i8] = fch3->data[i8];
  }

  c_fch1_idx_0 = r15->size[0];
  for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
    fchw->data[2 + fchw->size[0] * r15->data[i8]] = b_fch3->data[i8];
  }

  emxFree_real_T(&b_fch3);
  emxFree_int32_T(&r15);

  //  TODO: PREALLOCATE F BASED ON WINDOWLENGTH: (IF NECESSARY):
  //  all windows should be the same size:
  memset(&FFT[0], 0, 4100U * sizeof(double));
  memset(&fft_sel_loc[0], 0, sizeof(double) << 4);
  memset(&fft_sel_pks[0], 0, sizeof(double) << 4);
  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  emxInit_real_T(&PSD, 2);
  if (b_mod((double)b_fch1_idx_0, 2.0) == 1.0) {
    b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
    y = ((double)b_fch1_idx_0 - 1.0) / 2.0;
    i8 = PSD->size[0] * PSD->size[1];
    PSD->size[0] = 4;
    PSD->size[1] = (int)y;
    emxEnsureCapacity((emxArray__common *)PSD, i8, (int)sizeof(double));
    loop_ub = (int)y << 2;
    for (i8 = 0; i8 < loop_ub; i8++) {
      PSD->data[i8] = 0.0;
    }
  } else {
    b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
    y = (double)b_fch1_idx_0 / 2.0;
    i8 = PSD->size[0] * PSD->size[1];
    PSD->size[0] = 4;
    PSD->size[1] = (int)y;
    emxEnsureCapacity((emxArray__common *)PSD, i8, (int)sizeof(double));
    loop_ub = (int)y << 2;
    for (i8 = 0; i8 < loop_ub; i8++) {
      PSD->data[i8] = 0.0;
    }
  }

  //  PSD = zeros(4,windowLength/2);
  for (i8 = 0; i8 < 2; i8++) {
    iv0[i8] = PSD->size[i8];
  }

  emxInit_boolean_T1(&selectPSD, 2);
  i8 = selectPSD->size[0] * selectPSD->size[1];
  selectPSD->size[0] = 4;
  selectPSD->size[1] = iv0[1];
  emxEnsureCapacity((emxArray__common *)selectPSD, i8, (int)sizeof(boolean_T));
  loop_ub = iv0[1] << 2;
  for (i8 = 0; i8 < loop_ub; i8++) {
    selectPSD->data[i8] = false;
  }

  memset(&psd_sel_loc[0], 0, sizeof(double) << 4);
  memset(&psd_sel_pks[0], 0, sizeof(double) << 4);
  for (i = 0; i < 4; i++) {
    stft_sel_loc[i] = 0.0;
    stft_sel_pks[i] = 0.0;
  }

  //  Data is already filtered:
  //  TODO: FILL IN LATER!!!
  // {%}
  b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
  emxInit_real_T1(&hW, 1);
  emxInit_real_T(&fPSD, 2);
  emxInit_real_T1(&SummedRows, 1);
  emxInit_creal_T1(&S1, 2);
  emxInit_real_T(&T, 2);
  emxInit_creal_T1(&S2, 2);
  emxInit_creal_T1(&S3, 2);
  emxInit_real_T(&psd_sel_L, 2);
  emxInit_boolean_T1(&r16, 2);
  emxInit_boolean_T1(&r17, 2);
  emxInit_int32_T1(&r18, 2);
  emxInit_int32_T1(&r19, 2);
  emxInit_real_T(&x, 2);
  emxInit_real_T(&b_x, 2);
  emxInit_real_T(&c_x, 2);
  emxInit_real_T1(&b_SummedRows, 1);
  emxInit_real_T(&r20, 2);
  emxInit_creal_T1(&b_S3, 2);
  emxInit_creal_T1(&b_S2, 2);
  emxInit_creal_T1(&b_S1, 2);
  emxInit_real_T(&b_fchw, 2);
  emxInit_real_T(&c_fchw, 2);
  emxInit_real_T(&d_fchw, 2);
  emxInit_real_T(&e_fchw, 2);
  emxInit_real_T(&f_fchw, 2);
  emxInit_real_T(&b_FFT, 2);
  emxInit_real_T(&b_PSD, 2);
  emxInit_real_T(&c_FFT, 2);
  emxInit_real_T(&c_PSD, 2);
  emxInit_real_T1(&c_SummedRows, 1);
  emxInit_real_T(&g_fchw, 2);
  emxInit_real_T(&h_fchw, 2);
  emxInit_real_T(&d_FFT, 2);
  emxInit_real_T(&d_PSD, 2);
  emxInit_real_T(&e_FFT, 2);
  emxInit_real_T(&e_PSD, 2);
  emxInit_real_T(&f_PSD, 2);
  emxInit_int32_T(&r21, 1);
  emxInit_int32_T(&r22, 1);
  emxInit_int32_T(&r23, 1);
  emxInit_int32_T(&r24, 1);
  emxInit_real_T(&g_PSD, 2);
  emxInit_int32_T(&r25, 1);
  emxInit_int32_T(&r26, 1);
  guard2 = false;
  if (b_fch1_idx_0 < 500) {
    b_fch1_idx_0 = fch1->size[0] * fch1->size[1];
    if (b_fch1_idx_0 >= 250) {
      for (ch = 0; ch < 3; ch++) {
        loop_ub = fchw->size[1];
        i8 = h_fchw->size[0] * h_fchw->size[1];
        h_fchw->size[0] = 1;
        h_fchw->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)h_fchw, i8, (int)sizeof(double));
        for (i8 = 0; i8 < loop_ub; i8++) {
          h_fchw->data[h_fchw->size[0] * i8] = fchw->data[ch + fchw->size[0] *
            i8];
        }

        get_nfft_data(h_fchw, Fs, f, dv19);
        for (i8 = 0; i8 < 1025; i8++) {
          FFT[ch + (i8 << 2)] = dv19[i8];
        }

        for (i = 0; i < 4; i++) {
          for (i8 = 0; i8 < 1025; i8++) {
            b_select[i + (i8 << 2)] = ((f[i8] > threshFFT[i]) && (f[i8] <
              threshFFT[4 + i]));
          }

          //
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
            if (b_select[i + (c_fch1_idx_0 << 2)]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
            if (b_select[i + (c_fch1_idx_0 << 2)]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          i8 = d_FFT->size[0] * d_FFT->size[1];
          d_FFT->size[0] = 1;
          d_FFT->size[1] = r18->size[1];
          emxEnsureCapacity((emxArray__common *)d_FFT, i8, (int)sizeof(double));
          loop_ub = r18->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            d_FFT->data[d_FFT->size[0] * i8] = FFT[ch + ((r18->data[r18->size[0]
              * i8] - 1) << 2)];
          }

          c_findpeaks(d_FFT, T, psd_sel_L);
          if (!(T->size[1] == 0)) {
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
              if (b_select[i + (c_fch1_idx_0 << 2)]) {
                b_fch1_idx_0++;
              }
            }

            i8 = r18->size[0] * r18->size[1];
            r18->size[0] = 1;
            r18->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
              if (b_select[i + (c_fch1_idx_0 << 2)]) {
                r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            fft_sel_loc[ch + (i << 2)] = f[r18->data[r18->size[0] * ((int)
              psd_sel_L->data[0] - 1)] - 1];

            // verifies that maxes are peaks. [peak must occur w/i range]
            fft_sel_pks[ch + (i << 2)] = T->data[0];
          } else {
            fft_sel_loc[ch + (i << 2)] = 0.0;
            fft_sel_pks[ch + (i << 2)] = 0.0;
          }

          //
        }

        hannWin((double)fch1_idx_0, hW);
        loop_ub = fchw->size[1];
        i8 = g_fchw->size[0] * g_fchw->size[1];
        g_fchw->size[0] = 1;
        g_fchw->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)g_fchw, i8, (int)sizeof(double));
        for (i8 = 0; i8 < loop_ub; i8++) {
          g_fchw->data[g_fchw->size[0] * i8] = fchw->data[ch + fchw->size[0] *
            i8];
        }

        welch_psd(g_fchw, Fs, hW, T, fPSD);
        loop_ub = T->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          PSD->data[ch + PSD->size[0] * i8] = T->data[T->size[0] * i8];
        }

        // fin-start
        for (i = 0; i < 4; i++) {
          loop_ub = fPSD->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            selectPSD->data[i + selectPSD->size[0] * i8] = ((fPSD->data
              [fPSD->size[0] * i8] >= threshPSD[i]) && (fPSD->data[fPSD->size[0]
              * i8] <= threshPSD[4 + i]));
          }

          loop_ub = selectPSD->size[1];
          i8 = r16->size[0] * r16->size[1];
          r16->size[0] = 1;
          r16->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
          for (i8 = 0; i8 < loop_ub; i8++) {
            r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0]
              * i8];
          }

          loop_ub = selectPSD->size[1] - 1;
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          i8 = r26->size[0];
          r26->size[0] = r18->size[1];
          emxEnsureCapacity((emxArray__common *)r26, i8, (int)sizeof(int));
          loop_ub = r18->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            r26->data[i8] = r18->data[r18->size[0] * i8];
          }

          if (r26->size[0] >= 3) {
            loop_ub = selectPSD->size[1];
            i8 = r16->size[0] * r16->size[1];
            r16->size[0] = 1;
            r16->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
            for (i8 = 0; i8 < loop_ub; i8++) {
              r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->
                size[0] * i8];
            }

            loop_ub = selectPSD->size[1] - 1;
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0]) {
                b_fch1_idx_0++;
              }
            }

            i8 = r18->size[0] * r18->size[1];
            r18->size[0] = 1;
            r18->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0]) {
                r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            i8 = d_PSD->size[0] * d_PSD->size[1];
            d_PSD->size[0] = 1;
            d_PSD->size[1] = r18->size[1];
            emxEnsureCapacity((emxArray__common *)d_PSD, i8, (int)sizeof(double));
            loop_ub = r18->size[1];
            for (i8 = 0; i8 < loop_ub; i8++) {
              d_PSD->data[d_PSD->size[0] * i8] = PSD->data[ch + PSD->size[0] *
                (r18->data[r18->size[0] * i8] - 1)];
            }

            c_findpeaks(d_PSD, T, psd_sel_L);
            if (!(T->size[1] == 0)) {
              loop_ub = selectPSD->size[1];
              i8 = r16->size[0] * r16->size[1];
              r16->size[0] = 1;
              r16->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof
                                (boolean_T));
              for (i8 = 0; i8 < loop_ub; i8++) {
                r16->data[r16->size[0] * i8] = selectPSD->data[i +
                  selectPSD->size[0] * i8];
              }

              loop_ub = selectPSD->size[1] - 1;
              b_fch1_idx_0 = 0;
              for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
                if (r16->data[c_fch1_idx_0]) {
                  b_fch1_idx_0++;
                }
              }

              i8 = r18->size[0] * r18->size[1];
              r18->size[0] = 1;
              r18->size[1] = b_fch1_idx_0;
              emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
              b_fch1_idx_0 = 0;
              for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
                if (r16->data[c_fch1_idx_0]) {
                  r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                  b_fch1_idx_0++;
                }
              }

              psd_sel_loc[ch + (i << 2)] = fPSD->data[r18->data[r18->size[0] *
                ((int)psd_sel_L->data[0] - 1)] - 1];
              psd_sel_pks[ch + (i << 2)] = T->data[0];
            } else {
              psd_sel_loc[ch + (i << 2)] = 0.0;
              psd_sel_pks[ch + (i << 2)] = 0.0;
            }
          } else {
            psd_sel_loc[ch + (i << 2)] = 0.0;
            psd_sel_pks[ch + (i << 2)] = 0.0;
          }
        }

        // TODO FIND PEAKS:
      }

      for (i8 = 0; i8 < 1025; i8++) {
        FFT[3 + (i8 << 2)] = (FFT[i8 << 2] + FFT[1 + (i8 << 2)]) + FFT[2 + (i8 <<
          2)];
      }

      for (i = 0; i < 4; i++) {
        //
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if (b_select[i + (c_fch1_idx_0 << 2)]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if (b_select[i + (c_fch1_idx_0 << 2)]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = e_FFT->size[0] * e_FFT->size[1];
        e_FFT->size[0] = 1;
        e_FFT->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)e_FFT, i8, (int)sizeof(double));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          e_FFT->data[e_FFT->size[0] * i8] = FFT[3 + ((r18->data[r18->size[0] *
            i8] - 1) << 2)];
        }

        c_findpeaks(e_FFT, T, psd_sel_L);
        if (!(T->size[1] == 0)) {
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
            if (b_select[i + (c_fch1_idx_0 << 2)]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
            if (b_select[i + (c_fch1_idx_0 << 2)]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          fft_sel_loc[3 + (i << 2)] = f[r18->data[r18->size[0] * ((int)
            psd_sel_L->data[0] - 1)] - 1];

          // verifies that maxes are peaks. [peak must occur w/i range]
          fft_sel_pks[3 + (i << 2)] = T->data[0];
        } else {
          fft_sel_loc[3 + (i << 2)] = 0.0;
          fft_sel_pks[3 + (i << 2)] = 0.0;
        }
      }

      c_fch1_idx_0 = PSD->size[1];
      i8 = g_PSD->size[0] * g_PSD->size[1];
      g_PSD->size[0] = 1;
      g_PSD->size[1] = c_fch1_idx_0;
      emxEnsureCapacity((emxArray__common *)g_PSD, i8, (int)sizeof(double));
      for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
        g_PSD->data[g_PSD->size[0] * i8] = (PSD->data[PSD->size[0] * i8] +
          PSD->data[1 + PSD->size[0] * i8]) + PSD->data[2 + PSD->size[0] * i8];
      }

      loop_ub = g_PSD->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        PSD->data[3 + PSD->size[0] * i8] = g_PSD->data[g_PSD->size[0] * i8];
      }

      for (i = 0; i < 4; i++) {
        loop_ub = fPSD->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          selectPSD->data[i + selectPSD->size[0] * i8] = ((fPSD->data[fPSD->
            size[0] * i8] >= threshPSD[i]) && (fPSD->data[fPSD->size[0] * i8] <=
            threshPSD[4 + i]));
        }

        loop_ub = selectPSD->size[1];
        i8 = r16->size[0] * r16->size[1];
        r16->size[0] = 1;
        r16->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
        for (i8 = 0; i8 < loop_ub; i8++) {
          r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0] *
            i8];
        }

        loop_ub = selectPSD->size[1] - 1;
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r25->size[0];
        r25->size[0] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r25, i8, (int)sizeof(int));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r25->data[i8] = r18->data[r18->size[0] * i8];
        }

        if (r25->size[0] >= 3) {
          loop_ub = selectPSD->size[1];
          i8 = r16->size[0] * r16->size[1];
          r16->size[0] = 1;
          r16->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
          for (i8 = 0; i8 < loop_ub; i8++) {
            r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0]
              * i8];
          }

          loop_ub = selectPSD->size[1] - 1;
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          i8 = e_PSD->size[0] * e_PSD->size[1];
          e_PSD->size[0] = 1;
          e_PSD->size[1] = r18->size[1];
          emxEnsureCapacity((emxArray__common *)e_PSD, i8, (int)sizeof(double));
          loop_ub = r18->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            e_PSD->data[e_PSD->size[0] * i8] = PSD->data[3 + PSD->size[0] *
              (r18->data[r18->size[0] * i8] - 1)];
          }

          c_findpeaks(e_PSD, T, psd_sel_L);
          if (!(T->size[1] == 0)) {
            loop_ub = selectPSD->size[1];
            i8 = r16->size[0] * r16->size[1];
            r16->size[0] = 1;
            r16->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
            for (i8 = 0; i8 < loop_ub; i8++) {
              r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->
                size[0] * i8];
            }

            loop_ub = selectPSD->size[1] - 1;
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0]) {
                b_fch1_idx_0++;
              }
            }

            i8 = r18->size[0] * r18->size[1];
            r18->size[0] = 1;
            r18->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0]) {
                r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            psd_sel_loc[3 + (i << 2)] = fPSD->data[r18->data[r18->size[0] *
              ((int)psd_sel_L->data[0] - 1)] - 1];
            psd_sel_pks[3 + (i << 2)] = T->data[0];
          } else {
            psd_sel_loc[3 + (i << 2)] = 0.0;
            psd_sel_pks[3 + (i << 2)] = 0.0;
          }
        } else {
          psd_sel_loc[3 + (i << 2)] = 0.0;
          psd_sel_pks[3 + (i << 2)] = 0.0;
        }
      }
    } else {
      guard2 = true;
    }
  } else {
    guard2 = true;
  }

  if (guard2) {
    // ---------------- Data >=500 dp -----------------------%
    for (ch = 0; ch < 3; ch++) {
      loop_ub = fchw->size[1];
      i8 = f_fchw->size[0] * f_fchw->size[1];
      f_fchw->size[0] = 1;
      f_fchw->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)f_fchw, i8, (int)sizeof(double));
      for (i8 = 0; i8 < loop_ub; i8++) {
        f_fchw->data[f_fchw->size[0] * i8] = fchw->data[ch + fchw->size[0] * i8];
      }

      get_nfft_data(f_fchw, Fs, f, dv19);
      for (i8 = 0; i8 < 1025; i8++) {
        FFT[ch + (i8 << 2)] = dv19[i8];
      }

      for (i = 0; i < 4; i++) {
        for (i8 = 0; i8 < 1025; i8++) {
          b_select[i + (i8 << 2)] = ((f[i8] > threshFFT[i]) && (f[i8] <
            threshFFT[4 + i]));
        }

        //
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if (b_select[i + (c_fch1_idx_0 << 2)]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if (b_select[i + (c_fch1_idx_0 << 2)]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = b_FFT->size[0] * b_FFT->size[1];
        b_FFT->size[0] = 1;
        b_FFT->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)b_FFT, i8, (int)sizeof(double));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          b_FFT->data[b_FFT->size[0] * i8] = FFT[ch + ((r18->data[r18->size[0] *
            i8] - 1) << 2)];
        }

        c_findpeaks(b_FFT, T, psd_sel_L);
        if (!(T->size[1] == 0)) {
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
            if (b_select[i + (c_fch1_idx_0 << 2)]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
            if (b_select[i + (c_fch1_idx_0 << 2)]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          fft_sel_loc[ch + (i << 2)] = f[r18->data[r18->size[0] * ((int)
            psd_sel_L->data[0] - 1)] - 1];

          // verifies that maxes are peaks. [peak must occur w/i range]
          fft_sel_pks[ch + (i << 2)] = T->data[0];
        } else {
          fft_sel_loc[ch + (i << 2)] = 0.0;
          fft_sel_pks[ch + (i << 2)] = 0.0;
        }

        //
      }

      hannWin((double)fch1_idx_0, hW);
      loop_ub = fchw->size[1];
      i8 = e_fchw->size[0] * e_fchw->size[1];
      e_fchw->size[0] = 1;
      e_fchw->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)e_fchw, i8, (int)sizeof(double));
      for (i8 = 0; i8 < loop_ub; i8++) {
        e_fchw->data[e_fchw->size[0] * i8] = fchw->data[ch + fchw->size[0] * i8];
      }

      welch_psd(e_fchw, Fs, hW, T, fPSD);
      loop_ub = T->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        PSD->data[ch + PSD->size[0] * i8] = T->data[T->size[0] * i8];
      }

      // fin-start
      for (i = 0; i < 4; i++) {
        //              if len<1000
        //                  selectPSD(i,:) = fPSD>=threshPSD(i,1) & fPSD<=threshPSD(i,2); 
        //              else
        //                  selectPSD(i,:) = fPSD>=threshPSDL(i,1) & fPSD<=threshPSDL(i,2); 
        //              end
        loop_ub = fPSD->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          selectPSD->data[i + selectPSD->size[0] * i8] = ((fPSD->data[fPSD->
            size[0] * i8] >= threshPSD[i]) && (fPSD->data[fPSD->size[0] * i8] <=
            threshPSD[4 + i]));
        }

        loop_ub = selectPSD->size[1];
        i8 = r16->size[0] * r16->size[1];
        r16->size[0] = 1;
        r16->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
        for (i8 = 0; i8 < loop_ub; i8++) {
          r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0] *
            i8];
        }

        loop_ub = selectPSD->size[1] - 1;
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r24->size[0];
        r24->size[0] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r24, i8, (int)sizeof(int));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r24->data[i8] = r18->data[r18->size[0] * i8];
        }

        if (r24->size[0] >= 3) {
          loop_ub = selectPSD->size[1];
          i8 = r16->size[0] * r16->size[1];
          r16->size[0] = 1;
          r16->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
          for (i8 = 0; i8 < loop_ub; i8++) {
            r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0]
              * i8];
          }

          loop_ub = selectPSD->size[1] - 1;
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          i8 = b_PSD->size[0] * b_PSD->size[1];
          b_PSD->size[0] = 1;
          b_PSD->size[1] = r18->size[1];
          emxEnsureCapacity((emxArray__common *)b_PSD, i8, (int)sizeof(double));
          loop_ub = r18->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            b_PSD->data[b_PSD->size[0] * i8] = PSD->data[ch + PSD->size[0] *
              (r18->data[r18->size[0] * i8] - 1)];
          }

          c_findpeaks(b_PSD, T, psd_sel_L);
          if (!(T->size[1] == 0)) {
            loop_ub = selectPSD->size[1];
            i8 = r16->size[0] * r16->size[1];
            r16->size[0] = 1;
            r16->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
            for (i8 = 0; i8 < loop_ub; i8++) {
              r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->
                size[0] * i8];
            }

            loop_ub = selectPSD->size[1] - 1;
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0]) {
                b_fch1_idx_0++;
              }
            }

            i8 = r18->size[0] * r18->size[1];
            r18->size[0] = 1;
            r18->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0]) {
                r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            psd_sel_loc[ch + (i << 2)] = fPSD->data[r18->data[r18->size[0] *
              ((int)psd_sel_L->data[0] - 1)] - 1];
            psd_sel_pks[ch + (i << 2)] = T->data[0];
          } else {
            psd_sel_loc[ch + (i << 2)] = 0.0;
            psd_sel_pks[ch + (i << 2)] = 0.0;
          }
        } else {
          psd_sel_loc[ch + (i << 2)] = 0.0;
          psd_sel_pks[ch + (i << 2)] = 0.0;
        }
      }

      // TODO FIND PEAKS:
    }

    for (i8 = 0; i8 < 1025; i8++) {
      FFT[3 + (i8 << 2)] = (FFT[i8 << 2] + FFT[1 + (i8 << 2)]) + FFT[2 + (i8 <<
        2)];
    }

    for (i = 0; i < 4; i++) {
      for (i8 = 0; i8 < 1025; i8++) {
        b_select[i + (i8 << 2)] = ((f[i8] > threshFFT[i]) && (f[i8] < threshFFT
          [4 + i]));
      }

      //
      b_fch1_idx_0 = 0;
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
        if (b_select[i + (c_fch1_idx_0 << 2)]) {
          b_fch1_idx_0++;
        }
      }

      i8 = r18->size[0] * r18->size[1];
      r18->size[0] = 1;
      r18->size[1] = b_fch1_idx_0;
      emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
      b_fch1_idx_0 = 0;
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
        if (b_select[i + (c_fch1_idx_0 << 2)]) {
          r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
          b_fch1_idx_0++;
        }
      }

      i8 = c_FFT->size[0] * c_FFT->size[1];
      c_FFT->size[0] = 1;
      c_FFT->size[1] = r18->size[1];
      emxEnsureCapacity((emxArray__common *)c_FFT, i8, (int)sizeof(double));
      loop_ub = r18->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        c_FFT->data[c_FFT->size[0] * i8] = FFT[3 + ((r18->data[r18->size[0] * i8]
          - 1) << 2)];
      }

      c_findpeaks(c_FFT, T, psd_sel_L);
      if (!(T->size[1] == 0)) {
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if (b_select[i + (c_fch1_idx_0 << 2)]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if (b_select[i + (c_fch1_idx_0 << 2)]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        fft_sel_loc[3 + (i << 2)] = f[r18->data[r18->size[0] * ((int)
          psd_sel_L->data[0] - 1)] - 1];

        // verifies that maxes are peaks. [peak must occur w/i range]
        fft_sel_pks[3 + (i << 2)] = T->data[0];
      } else {
        fft_sel_loc[3 + (i << 2)] = 0.0;
        fft_sel_pks[3 + (i << 2)] = 0.0;
      }
    }

    c_fch1_idx_0 = PSD->size[1];
    i8 = f_PSD->size[0] * f_PSD->size[1];
    f_PSD->size[0] = 1;
    f_PSD->size[1] = c_fch1_idx_0;
    emxEnsureCapacity((emxArray__common *)f_PSD, i8, (int)sizeof(double));
    for (i8 = 0; i8 < c_fch1_idx_0; i8++) {
      f_PSD->data[f_PSD->size[0] * i8] = (PSD->data[PSD->size[0] * i8] +
        PSD->data[1 + PSD->size[0] * i8]) + PSD->data[2 + PSD->size[0] * i8];
    }

    loop_ub = f_PSD->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      PSD->data[3 + PSD->size[0] * i8] = f_PSD->data[f_PSD->size[0] * i8];
    }

    for (i = 0; i < 4; i++) {
      loop_ub = fPSD->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        selectPSD->data[i + selectPSD->size[0] * i8] = ((fPSD->data[fPSD->size[0]
          * i8] >= threshPSD[i]) && (fPSD->data[fPSD->size[0] * i8] <=
          threshPSD[4 + i]));
      }

      //          if len<1000
      //              selectPSD(i,:) = fPSD>=threshPSD(i,1) & fPSD<=threshPSD(i,2); 
      //          else
      //              selectPSD(i,:) = fPSD>=threshPSDL(i,1) & fPSD<=threshPSDL(i,2); 
      //          end
      loop_ub = selectPSD->size[1];
      i8 = r16->size[0] * r16->size[1];
      r16->size[0] = 1;
      r16->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
      for (i8 = 0; i8 < loop_ub; i8++) {
        r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0] *
          i8];
      }

      loop_ub = selectPSD->size[1] - 1;
      b_fch1_idx_0 = 0;
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
        if (r16->data[c_fch1_idx_0]) {
          b_fch1_idx_0++;
        }
      }

      i8 = r18->size[0] * r18->size[1];
      r18->size[0] = 1;
      r18->size[1] = b_fch1_idx_0;
      emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
      b_fch1_idx_0 = 0;
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
        if (r16->data[c_fch1_idx_0]) {
          r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
          b_fch1_idx_0++;
        }
      }

      i8 = r23->size[0];
      r23->size[0] = r18->size[1];
      emxEnsureCapacity((emxArray__common *)r23, i8, (int)sizeof(int));
      loop_ub = r18->size[1];
      for (i8 = 0; i8 < loop_ub; i8++) {
        r23->data[i8] = r18->data[r18->size[0] * i8];
      }

      if (r23->size[0] >= 3) {
        loop_ub = selectPSD->size[1];
        i8 = r16->size[0] * r16->size[1];
        r16->size[0] = 1;
        r16->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
        for (i8 = 0; i8 < loop_ub; i8++) {
          r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0] *
            i8];
        }

        loop_ub = selectPSD->size[1] - 1;
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = c_PSD->size[0] * c_PSD->size[1];
        c_PSD->size[0] = 1;
        c_PSD->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)c_PSD, i8, (int)sizeof(double));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          c_PSD->data[c_PSD->size[0] * i8] = PSD->data[3 + PSD->size[0] *
            (r18->data[r18->size[0] * i8] - 1)];
        }

        c_findpeaks(c_PSD, T, psd_sel_L);
        if (!(T->size[1] == 0)) {
          loop_ub = selectPSD->size[1];
          i8 = r16->size[0] * r16->size[1];
          r16->size[0] = 1;
          r16->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
          for (i8 = 0; i8 < loop_ub; i8++) {
            r16->data[r16->size[0] * i8] = selectPSD->data[i + selectPSD->size[0]
              * i8];
          }

          loop_ub = selectPSD->size[1] - 1;
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          psd_sel_loc[3 + (i << 2)] = fPSD->data[r18->data[r18->size[0] * ((int)
            psd_sel_L->data[0] - 1)] - 1];
          psd_sel_pks[3 + (i << 2)] = T->data[0];
        } else {
          psd_sel_loc[3 + (i << 2)] = 0.0;
          psd_sel_pks[3 + (i << 2)] = 0.0;
        }
      } else {
        psd_sel_loc[3 + (i << 2)] = 0.0;
        psd_sel_pks[3 + (i << 2)] = 0.0;
      }
    }

    // Classification method #2 (w/ STFT):
    // TODO:
    loop_ub = fchw->size[1];
    i8 = d_fchw->size[0] * d_fchw->size[1];
    d_fchw->size[0] = 1;
    d_fchw->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)d_fchw, i8, (int)sizeof(double));
    for (i8 = 0; i8 < loop_ub; i8++) {
      d_fchw->data[d_fchw->size[0] * i8] = fchw->data[fchw->size[0] * i8];
    }

    stft(d_fchw, Fs, S1, b_F, T);
    loop_ub = fchw->size[1];
    i8 = c_fchw->size[0] * c_fchw->size[1];
    c_fchw->size[0] = 1;
    c_fchw->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)c_fchw, i8, (int)sizeof(double));
    for (i8 = 0; i8 < loop_ub; i8++) {
      c_fchw->data[c_fchw->size[0] * i8] = fchw->data[1 + fchw->size[0] * i8];
    }

    stft(c_fchw, Fs, S2, f, T);
    loop_ub = fchw->size[1];
    i8 = b_fchw->size[0] * b_fchw->size[1];
    b_fchw->size[0] = 1;
    b_fchw->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_fchw, i8, (int)sizeof(double));
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_fchw->data[b_fchw->size[0] * i8] = fchw->data[2 + fchw->size[0] * i8];
    }

    stft(b_fchw, Fs, S3, f, T);
    b_fch1_idx_0 = 0;
    for (i = 0; i < 1025; i++) {
      if ((b_F[i] < 17.6) && (b_F[i] > 9.0)) {
        b_fch1_idx_0++;
      }
    }

    i8 = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = b_fch1_idx_0;
    emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
    b_fch1_idx_0 = 0;
    for (i = 0; i < 1025; i++) {
      if ((b_F[i] < 17.6) && (b_F[i] > 9.0)) {
        r18->data[b_fch1_idx_0] = i + 1;
        b_fch1_idx_0++;
      }
    }

    loop_ub = S1->size[1];
    i8 = b_S1->size[0] * b_S1->size[1];
    b_S1->size[0] = r18->size[1];
    b_S1->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_S1, i8, (int)sizeof(creal_T));
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_fch1_idx_0 = r18->size[1];
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 < b_fch1_idx_0; c_fch1_idx_0++) {
        b_S1->data[c_fch1_idx_0 + b_S1->size[0] * i8] = S1->data[(r18->data
          [r18->size[0] * c_fch1_idx_0] + S1->size[0] * i8) - 1];
      }
    }

    c_abs(b_S1, x);
    i8 = x->size[0] * x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i8, (int)sizeof(double));
    c_fch1_idx_0 = x->size[0];
    b_fch1_idx_0 = x->size[1];
    loop_ub = c_fch1_idx_0 * b_fch1_idx_0;
    for (i8 = 0; i8 < loop_ub; i8++) {
      x->data[i8] = x->data[i8] / 256.0 / 0.54000000000000059 + 1.0E-6;
    }

    b_log10(x);
    b_fch1_idx_0 = 0;
    for (i = 0; i < 1025; i++) {
      if ((b_F[i] < 17.6) && (b_F[i] > 9.0)) {
        b_fch1_idx_0++;
      }
    }

    i8 = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = b_fch1_idx_0;
    emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
    b_fch1_idx_0 = 0;
    for (i = 0; i < 1025; i++) {
      if ((b_F[i] < 17.6) && (b_F[i] > 9.0)) {
        r18->data[b_fch1_idx_0] = i + 1;
        b_fch1_idx_0++;
      }
    }

    loop_ub = S2->size[1];
    i8 = b_S2->size[0] * b_S2->size[1];
    b_S2->size[0] = r18->size[1];
    b_S2->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_S2, i8, (int)sizeof(creal_T));
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_fch1_idx_0 = r18->size[1];
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 < b_fch1_idx_0; c_fch1_idx_0++) {
        b_S2->data[c_fch1_idx_0 + b_S2->size[0] * i8] = S2->data[(r18->data
          [r18->size[0] * c_fch1_idx_0] + S2->size[0] * i8) - 1];
      }
    }

    c_abs(b_S2, b_x);
    i8 = b_x->size[0] * b_x->size[1];
    emxEnsureCapacity((emxArray__common *)b_x, i8, (int)sizeof(double));
    c_fch1_idx_0 = b_x->size[0];
    b_fch1_idx_0 = b_x->size[1];
    loop_ub = c_fch1_idx_0 * b_fch1_idx_0;
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_x->data[i8] = b_x->data[i8] / 256.0 / 0.54000000000000059 + 1.0E-6;
    }

    b_log10(b_x);
    b_fch1_idx_0 = 0;
    for (i = 0; i < 1025; i++) {
      if ((b_F[i] < 17.6) && (b_F[i] > 9.0)) {
        b_fch1_idx_0++;
      }
    }

    i8 = r18->size[0] * r18->size[1];
    r18->size[0] = 1;
    r18->size[1] = b_fch1_idx_0;
    emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
    b_fch1_idx_0 = 0;
    for (i = 0; i < 1025; i++) {
      if ((b_F[i] < 17.6) && (b_F[i] > 9.0)) {
        r18->data[b_fch1_idx_0] = i + 1;
        b_fch1_idx_0++;
      }
    }

    loop_ub = S3->size[1];
    i8 = b_S3->size[0] * b_S3->size[1];
    b_S3->size[0] = r18->size[1];
    b_S3->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_S3, i8, (int)sizeof(creal_T));
    for (i8 = 0; i8 < loop_ub; i8++) {
      b_fch1_idx_0 = r18->size[1];
      for (c_fch1_idx_0 = 0; c_fch1_idx_0 < b_fch1_idx_0; c_fch1_idx_0++) {
        b_S3->data[c_fch1_idx_0 + b_S3->size[0] * i8] = S3->data[(r18->data
          [r18->size[0] * c_fch1_idx_0] + S3->size[0] * i8) - 1];
      }
    }

    c_abs(b_S3, c_x);
    i8 = c_x->size[0] * c_x->size[1];
    emxEnsureCapacity((emxArray__common *)c_x, i8, (int)sizeof(double));
    c_fch1_idx_0 = c_x->size[0];
    b_fch1_idx_0 = c_x->size[1];
    loop_ub = c_fch1_idx_0 * b_fch1_idx_0;
    for (i8 = 0; i8 < loop_ub; i8++) {
      c_x->data[i8] = c_x->data[i8] / 256.0 / 0.54000000000000059 + 1.0E-6;
    }

    b_log10(c_x);
    i8 = r20->size[0] * r20->size[1];
    r20->size[0] = x->size[0];
    r20->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)r20, i8, (int)sizeof(double));
    loop_ub = x->size[0] * x->size[1];
    for (i8 = 0; i8 < loop_ub; i8++) {
      r20->data[i8] = (20.0 * x->data[i8] + 20.0 * b_x->data[i8]) + 20.0 *
        c_x->data[i8];
    }

    sum(r20, hW);
    scaleAbs(hW, SummedRows);
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
    if (fch1_idx_0 < 1250) {
      for (i = 0; i < 4; i++) {
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r16->size[0] * r16->size[1];
        r16->size[0] = 1;
        r16->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
        loop_ub = r18->size[0] * r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r16->data[i8] = (b_F[r18->data[i8] - 1] >= threshSTFT[i]);
        }

        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r17->size[0] * r17->size[1];
        r17->size[0] = 1;
        r17->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r17, i8, (int)sizeof(boolean_T));
        loop_ub = r18->size[0] * r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r17->data[i8] = (b_F[r18->data[i8] - 1] <= threshSTFT[4 + i]);
        }

        loop_ub = r16->size[1] - 1;
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r22->size[0];
        r22->size[0] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r22, i8, (int)sizeof(int));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r22->data[i8] = r18->data[r18->size[0] * i8];
        }

        if (r22->size[0] > 3) {
          loop_ub = r16->size[1] - 1;
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          i8 = c_SummedRows->size[0];
          c_SummedRows->size[0] = r18->size[1];
          emxEnsureCapacity((emxArray__common *)c_SummedRows, i8, (int)sizeof
                            (double));
          loop_ub = r18->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            c_SummedRows->data[i8] = SummedRows->data[r18->data[r18->size[0] *
              i8] - 1];
          }

          d_findpeaks(c_SummedRows, x, b_x);
          if (!((x->size[0] == 0) || (x->size[1] == 0))) {
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
              if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
                b_fch1_idx_0++;
              }
            }

            i8 = r18->size[0] * r18->size[1];
            r18->size[0] = 1;
            r18->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
              if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
                r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            loop_ub = r16->size[1] - 1;
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
                b_fch1_idx_0++;
              }
            }

            i8 = r19->size[0] * r19->size[1];
            r19->size[0] = 1;
            r19->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r19, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
                r19->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            stft_sel_loc[i] = b_F[r18->data[r18->size[0] * (r19->data[r19->size
              [0] * ((int)b_x->data[0] - 1)] - 1)] - 1];
            stft_sel_pks[i] = x->data[0];
          } else {
            stft_sel_loc[i] = 0.0;
            stft_sel_pks[i] = 0.0;
          }
        } else {
          stft_sel_loc[i] = 0.0;
          stft_sel_pks[i] = 0.0;
        }
      }
    } else {
      // (>=1250).
      for (i = 0; i < 4; i++) {
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r16->size[0] * r16->size[1];
        r16->size[0] = 1;
        r16->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r16, i8, (int)sizeof(boolean_T));
        loop_ub = r18->size[0] * r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r16->data[i8] = (b_F[r18->data[i8] - 1] >= threshSTFT[i]);
        }

        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
          if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r17->size[0] * r17->size[1];
        r17->size[0] = 1;
        r17->size[1] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r17, i8, (int)sizeof(boolean_T));
        loop_ub = r18->size[0] * r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r17->data[i8] = (b_F[r18->data[i8] - 1] <= threshSTFT[4 + i]);
        }

        loop_ub = r16->size[1] - 1;
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
            b_fch1_idx_0++;
          }
        }

        i8 = r18->size[0] * r18->size[1];
        r18->size[0] = 1;
        r18->size[1] = b_fch1_idx_0;
        emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
        b_fch1_idx_0 = 0;
        for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
          if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
            r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
            b_fch1_idx_0++;
          }
        }

        i8 = r21->size[0];
        r21->size[0] = r18->size[1];
        emxEnsureCapacity((emxArray__common *)r21, i8, (int)sizeof(int));
        loop_ub = r18->size[1];
        for (i8 = 0; i8 < loop_ub; i8++) {
          r21->data[i8] = r18->data[r18->size[0] * i8];
        }

        if (r21->size[0] >= 3) {
          loop_ub = r16->size[1] - 1;
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
              b_fch1_idx_0++;
            }
          }

          i8 = r18->size[0] * r18->size[1];
          r18->size[0] = 1;
          r18->size[1] = b_fch1_idx_0;
          emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
          b_fch1_idx_0 = 0;
          for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
            if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
              r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
              b_fch1_idx_0++;
            }
          }

          i8 = b_SummedRows->size[0];
          b_SummedRows->size[0] = r18->size[1];
          emxEnsureCapacity((emxArray__common *)b_SummedRows, i8, (int)sizeof
                            (double));
          loop_ub = r18->size[1];
          for (i8 = 0; i8 < loop_ub; i8++) {
            b_SummedRows->data[i8] = SummedRows->data[r18->data[r18->size[0] *
              i8] - 1];
          }

          d_findpeaks(b_SummedRows, x, b_x);
          if (!((x->size[0] == 0) || (x->size[1] == 0))) {
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
              if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
                b_fch1_idx_0++;
              }
            }

            i8 = r18->size[0] * r18->size[1];
            r18->size[0] = 1;
            r18->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r18, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 < 1025; c_fch1_idx_0++) {
              if ((b_F[c_fch1_idx_0] < 17.6) && (b_F[c_fch1_idx_0] > 9.0)) {
                r18->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            loop_ub = r16->size[1] - 1;
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
                b_fch1_idx_0++;
              }
            }

            i8 = r19->size[0] * r19->size[1];
            r19->size[0] = 1;
            r19->size[1] = b_fch1_idx_0;
            emxEnsureCapacity((emxArray__common *)r19, i8, (int)sizeof(int));
            b_fch1_idx_0 = 0;
            for (c_fch1_idx_0 = 0; c_fch1_idx_0 <= loop_ub; c_fch1_idx_0++) {
              if (r16->data[c_fch1_idx_0] && r17->data[c_fch1_idx_0]) {
                r19->data[b_fch1_idx_0] = c_fch1_idx_0 + 1;
                b_fch1_idx_0++;
              }
            }

            stft_sel_loc[i] = b_F[r18->data[r18->size[0] * (r19->data[r19->size
              [0] * ((int)b_x->data[0] - 1)] - 1)] - 1];
            stft_sel_pks[i] = x->data[0];
          } else {
            stft_sel_loc[i] = 0.0;
            stft_sel_pks[i] = 0.0;
          }
        } else {
          stft_sel_loc[i] = 0.0;
          stft_sel_pks[i] = 0.0;
        }
      }
    }
  }

  emxFree_int32_T(&r26);
  emxFree_int32_T(&r25);
  emxFree_real_T(&g_PSD);
  emxFree_int32_T(&r24);
  emxFree_int32_T(&r23);
  emxFree_int32_T(&r22);
  emxFree_int32_T(&r21);
  emxFree_real_T(&f_PSD);
  emxFree_real_T(&e_PSD);
  emxFree_real_T(&e_FFT);
  emxFree_real_T(&d_PSD);
  emxFree_real_T(&d_FFT);
  emxFree_real_T(&h_fchw);
  emxFree_real_T(&g_fchw);
  emxFree_real_T(&c_SummedRows);
  emxFree_real_T(&c_PSD);
  emxFree_real_T(&c_FFT);
  emxFree_real_T(&b_PSD);
  emxFree_real_T(&b_FFT);
  emxFree_real_T(&f_fchw);
  emxFree_real_T(&e_fchw);
  emxFree_real_T(&d_fchw);
  emxFree_real_T(&c_fchw);
  emxFree_real_T(&b_fchw);
  emxFree_creal_T(&b_S1);
  emxFree_creal_T(&b_S2);
  emxFree_creal_T(&b_S3);
  emxFree_real_T(&r20);
  emxFree_real_T(&b_SummedRows);
  emxFree_real_T(&c_x);
  emxFree_real_T(&b_x);
  emxFree_real_T(&x);
  emxFree_int32_T(&r19);
  emxFree_int32_T(&r18);
  emxFree_boolean_T(&r17);
  emxFree_boolean_T(&r16);
  emxFree_real_T(&psd_sel_L);
  emxFree_creal_T(&S3);
  emxFree_creal_T(&S2);
  emxFree_real_T(&T);
  emxFree_creal_T(&S1);
  emxFree_real_T(&SummedRows);
  emxFree_real_T(&fPSD);
  emxFree_real_T(&hW);
  emxFree_boolean_T(&selectPSD);
  emxFree_real_T(&PSD);
  emxFree_real_T(&fchw);

  // % Analysis & Collection:
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  guard1 = false;
  if (fch1_idx_0 < 500) {
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
    if (fch1_idx_0 >= 250) {
      // [ LOCATION    , MAGNITUDE
      for (i = 0; i < 4; i++) {
        for (i8 = 0; i8 < 4; i8++) {
          b_ft_ch[i + (i8 << 2)] = fft_sel_loc[i8 + (i << 2)];
          b_ft_ch[i + ((i8 + 4) << 2)] = psd_sel_loc[i8 + (i << 2)];
          b_ft_ch[i + ((i8 + 8) << 2)] = fft_sel_pks[i8 + (i << 2)];
          b_ft_ch[i + ((i8 + 12) << 2)] = psd_sel_pks[i8 + (i << 2)];
        }
      }

      for (i8 = 0; i8 < 16; i8++) {
        d_ft_ch[i8] = b_ft_ch[i8 << 2];
        d_ft_ch[i8 + 16] = b_ft_ch[1 + (i8 << 2)];
        d_ft_ch[i8 + 32] = b_ft_ch[2 + (i8 << 2)];
        d_ft_ch[i8 + 48] = b_ft_ch[3 + (i8 << 2)];
      }

      i8 = F->size[0] * F->size[1];
      F->size[0] = 1;
      F->size[1] = 64;
      emxEnsureCapacity((emxArray__common *)F, i8, (int)sizeof(double));
      for (i8 = 0; i8 < 64; i8++) {
        F->data[F->size[0] * i8] = d_ft_ch[i8];
      }

      // [1x64]
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    for (i = 0; i < 4; i++) {
      for (i8 = 0; i8 < 4; i8++) {
        ft_ch[i + (i8 << 2)] = fft_sel_loc[i8 + (i << 2)];
        ft_ch[i + ((i8 + 4) << 2)] = psd_sel_loc[i8 + (i << 2)];
        ft_ch[i + ((i8 + 8) << 2)] = fft_sel_pks[i8 + (i << 2)];
        ft_ch[i + ((i8 + 12) << 2)] = psd_sel_pks[i8 + (i << 2)];
      }

      ft_ch[64 + i] = stft_sel_loc[i];
      ft_ch[68 + i] = stft_sel_pks[i];
    }

    for (i8 = 0; i8 < 18; i8++) {
      c_ft_ch[i8] = ft_ch[i8 << 2];
      c_ft_ch[i8 + 18] = ft_ch[1 + (i8 << 2)];
      c_ft_ch[i8 + 36] = ft_ch[2 + (i8 << 2)];
      c_ft_ch[i8 + 54] = ft_ch[3 + (i8 << 2)];
    }

    i8 = F->size[0] * F->size[1];
    F->size[0] = 1;
    F->size[1] = 72;
    emxEnsureCapacity((emxArray__common *)F, i8, (int)sizeof(double));
    for (i8 = 0; i8 < 72; i8++) {
      F->data[F->size[0] * i8] = c_ft_ch[i8];
    }

    // [1x72]
  }

  // END FUNCTION
}

//
// featureExtraction Summary of this function goes here if I ever feel like
// writing one up.
//  samplesX = samplesX(:);
// Arguments    : const double samplesX[250]
//                emxArray_real_T *F
// Return Type  : void
//
static void featureExtractionEOG(const double samplesX[250], emxArray_real_T *F)
{
  double y;
  int ix;
  double xbar;
  int k;
  double r;
  double b_y;
  int ixstart;
  boolean_T exitg2;
  boolean_T exitg1;
  boolean_T x[250];
  double T_countmin_1;
  double T_countmin_2;
  double T_countmax;
  double s;
  double ylast;
  emxArray_real_T *peaks;
  emxArray_real_T *loc;
  emxArray_real_T *T_findpeaks_distX;
  int T_count_findpeaks;
  y = samplesX[0];
  ix = 0;
  xbar = samplesX[0];
  for (k = 0; k < 249; k++) {
    y += samplesX[k + 1];
    ix++;
    xbar += samplesX[ix];
  }

  xbar /= 250.0;
  ix = 0;
  r = samplesX[0] - xbar;
  b_y = r * r;
  for (k = 0; k < 249; k++) {
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

  for (ixstart = 0; ixstart < 250; ixstart++) {
    x[ixstart] = ((samplesX[ixstart] < -9.9999999999999991E-6) &&
                  (samplesX[ixstart] > -0.0001));
  }

  T_countmin_1 = x[0];
  for (k = 0; k < 249; k++) {
    T_countmin_1 += (double)x[k + 1];
  }

  for (ixstart = 0; ixstart < 250; ixstart++) {
    x[ixstart] = (samplesX[ixstart] < -0.0001);
  }

  T_countmin_2 = x[0];
  for (k = 0; k < 249; k++) {
    T_countmin_2 += (double)x[k + 1];
  }

  for (ixstart = 0; ixstart < 250; ixstart++) {
    x[ixstart] = (samplesX[ixstart] > 8.4999999999999993E-5);
  }

  T_countmax = x[0];
  s = 0.0;
  ixstart = 0;
  ylast = samplesX[0];
  for (k = 0; k < 249; k++) {
    T_countmax += (double)x[k + 1];
    ixstart++;
    s += (ylast + samplesX[ixstart]) / 2.0;
    ylast = samplesX[ixstart];
  }

  emxInit_real_T(&peaks, 2);
  emxInit_real_T(&loc, 2);
  findpeaks(samplesX, peaks, loc);
  emxInit_real_T(&T_findpeaks_distX, 2);
  if (peaks->size[1] == 0) {
    T_count_findpeaks = 0;
    ixstart = T_findpeaks_distX->size[0] * T_findpeaks_distX->size[1];
    T_findpeaks_distX->size[0] = 1;
    T_findpeaks_distX->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)T_findpeaks_distX, ixstart, (int)
                      sizeof(double));
    T_findpeaks_distX->data[0] = 0.0;
  } else {
    ixstart = peaks->size[1];
    if (1 >= ixstart) {
      ixstart = 1;
    }

    if (peaks->size[1] == 0) {
      ixstart = 0;
    }

    T_count_findpeaks = ixstart;
    ixstart = peaks->size[1];
    if (1 >= ixstart) {
      ixstart = 1;
    }

    if (peaks->size[1] == 0) {
      ixstart = 0;
    }

    if (ixstart > 1) {
      ixstart = T_findpeaks_distX->size[0] * T_findpeaks_distX->size[1];
      T_findpeaks_distX->size[0] = 1;
      T_findpeaks_distX->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)T_findpeaks_distX, ixstart, (int)
                        sizeof(double));
      T_findpeaks_distX->data[0] = loc->data[loc->size[1] - 1] - loc->data[0];

      // TODO: TAKE AVG, NOT MAX-MIN
    } else {
      ixstart = T_findpeaks_distX->size[0] * T_findpeaks_distX->size[1];
      T_findpeaks_distX->size[0] = 1;
      T_findpeaks_distX->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)T_findpeaks_distX, ixstart, (int)
                        sizeof(double));
      T_findpeaks_distX->data[0] = 0.0;
    }
  }

  emxFree_real_T(&loc);
  emxFree_real_T(&peaks);

  //  F = horzcat(T_mean, T_stdv, T_max, T_min, T_countmin_1, T_countmin_2, T_countmax, T_Integrate); 
  ixstart = F->size[0] * F->size[1];
  F->size[0] = 1;
  F->size[1] = 10;
  emxEnsureCapacity((emxArray__common *)F, ixstart, (int)sizeof(double));
  F->data[0] = y / 250.0;
  F->data[F->size[0]] = std::sqrt(b_y);
  F->data[F->size[0] << 1] = xbar;
  F->data[F->size[0] * 3] = r;
  F->data[F->size[0] << 2] = T_countmin_1;
  F->data[F->size[0] * 5] = T_countmin_2;
  F->data[F->size[0] * 6] = T_countmax;
  F->data[F->size[0] * 7] = s;
  F->data[F->size[0] << 3] = T_count_findpeaks;
  for (ixstart = 0; ixstart < 1; ixstart++) {
    F->data[F->size[0] * 9] = T_findpeaks_distX->data[0];
  }

  emxFree_real_T(&T_findpeaks_distX);
}

//
// Feature Extraction Function for Tri-channel SSVEP Feature Extraction:
//  ----- INPUTS -----
//  fch1, fch2, fch3: Tri-channel SSVEP Samples of certain window size
//  MUST BE A VECTOR!
//  MUST BE FILTERED!!
//  Fs - Sampling Rate.
//  Using 250-sample windows, feature extraction is obtained using FFT and
//  PSD
// ----FFT----%
// Arguments    : const emxArray_real_T *fch1
//                const emxArray_real_T *fch2
//                const emxArray_real_T *fch3
//                double Fs
//                double F[30]
// Return Type  : void
//
static void featureExtractionSSVEP(const emxArray_real_T *fch1, const
  emxArray_real_T *fch2, const emxArray_real_T *fch3, double Fs, double F[30])
{
  double threshFFT[8];
  int i3;
  signed char wLFFT[4];
  int i;
  double threshPSD[8];
  signed char wLPSD[4];
  emxArray_real_T *fchw;
  int fch1_idx_0;
  emxArray_int32_T *r8;
  emxArray_real_T *b_fch1;
  emxArray_real_T *b_fch2;
  emxArray_real_T *b_fch3;
  double FFT_Ltop[8];
  double FFT_PkRatio[4];
  emxArray_real_T *PSD;
  double y;
  double PSD_Ltop[8];
  double PSD_PkRatio[4];
  emxArray_real_T *hW;
  emxArray_real_T *fPSD;
  double f[1025];
  emxArray_real_T *FFT_PKS;
  emxArray_real_T *FFT_L;
  emxArray_real_T *b_PSD;
  emxArray_real_T *b_fchw;
  emxArray_real_T *c_fchw;
  double FFT[4100];
  int ch;
  double b_FFT[1025];
  double dv16[1025];
  emxArray_real_T *c_PSD;
  int w;
  boolean_T exitg2;
  boolean_T exitg4;
  emxArray_real_T *d_PSD;
  boolean_T exitg1;
  double FFTPeaks1[4];
  double PSDPeaks1[4];
  int chn;
  double b_FFT_Ltop[4];
  boolean_T exitg3;
  boolean_T b0;
  boolean_T b2;
  double b_PSD_Ltop[4];
  double b_wLFFT[4];
  double c_wLFFT[4];
  boolean_T b1;
  boolean_T b4;
  double b_wLPSD[4];
  double c_wLPSD[4];

  // -% windows around certain target frequencies
  for (i3 = 0; i3 < 2; i3++) {
    threshFFT[i3 << 2] = 9.5 + 1.1300000000000008 * (double)i3;
    threshFFT[1 + (i3 << 2)] = 11.9 + 0.79999999999999893 * (double)i3;
    threshFFT[2 + (i3 << 2)] = 14.6 + 0.90000000000000036 * (double)i3;
    threshFFT[3 + (i3 << 2)] = 16.2 + 0.53999999999999915 * (double)i3;
  }

  for (i = 0; i < 4; i++) {
    wLFFT[i] = 0;
  }

  // ----PSD----%
  for (i3 = 0; i3 < 2; i3++) {
    threshPSD[i3 << 2] = 9.5 + (double)i3;
    threshPSD[1 + (i3 << 2)] = 12.0 + (double)i3;
    threshPSD[2 + (i3 << 2)] = 14.9 + 0.19999999999999929 * (double)i3;
    threshPSD[3 + (i3 << 2)] = 16.0 + (double)i3;
  }

  for (i = 0; i < 4; i++) {
    wLPSD[i] = 0;
  }

  emxInit_real_T(&fchw, 2);

  // ----PREALLOCATE----%
  i3 = fchw->size[0] * fchw->size[1];
  fchw->size[0] = 3;
  emxEnsureCapacity((emxArray__common *)fchw, i3, (int)sizeof(double));
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  i3 = fchw->size[0] * fchw->size[1];
  fchw->size[1] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)fchw, i3, (int)sizeof(double));
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  i = 3 * fch1_idx_0;
  for (i3 = 0; i3 < i; i3++) {
    fchw->data[i3] = 0.0;
  }

  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > fch1_idx_0) {
    fch1_idx_0 = 0;
  } else {
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  i = fch1->size[0] * fch1->size[1];
  if (1 > i) {
    i = 0;
  } else {
    i = fch1->size[0] * fch1->size[1];
  }

  emxInit_int32_T(&r8, 1);
  i3 = r8->size[0];
  r8->size[0] = i;
  emxEnsureCapacity((emxArray__common *)r8, i3, (int)sizeof(int));
  for (i3 = 0; i3 < i; i3++) {
    r8->data[i3] = i3;
  }

  emxInit_real_T1(&b_fch1, 1);
  i3 = b_fch1->size[0];
  b_fch1->size[0] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch1, i3, (int)sizeof(double));
  for (i3 = 0; i3 < fch1_idx_0; i3++) {
    b_fch1->data[i3] = fch1->data[i3];
  }

  i = r8->size[0];
  for (i3 = 0; i3 < i; i3++) {
    fchw->data[fchw->size[0] * r8->data[i3]] = b_fch1->data[i3];
  }

  emxFree_real_T(&b_fch1);
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > fch1_idx_0) {
    fch1_idx_0 = 0;
  } else {
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  i = fch1->size[0] * fch1->size[1];
  if (1 > i) {
    i = 0;
  } else {
    i = fch1->size[0] * fch1->size[1];
  }

  i3 = r8->size[0];
  r8->size[0] = i;
  emxEnsureCapacity((emxArray__common *)r8, i3, (int)sizeof(int));
  for (i3 = 0; i3 < i; i3++) {
    r8->data[i3] = i3;
  }

  emxInit_real_T1(&b_fch2, 1);
  i3 = b_fch2->size[0];
  b_fch2->size[0] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch2, i3, (int)sizeof(double));
  for (i3 = 0; i3 < fch1_idx_0; i3++) {
    b_fch2->data[i3] = fch2->data[i3];
  }

  i = r8->size[0];
  for (i3 = 0; i3 < i; i3++) {
    fchw->data[1 + fchw->size[0] * r8->data[i3]] = b_fch2->data[i3];
  }

  emxFree_real_T(&b_fch2);
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  if (1 > fch1_idx_0) {
    fch1_idx_0 = 0;
  } else {
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
  }

  i = fch1->size[0] * fch1->size[1];
  if (1 > i) {
    i = 0;
  } else {
    i = fch1->size[0] * fch1->size[1];
  }

  i3 = r8->size[0];
  r8->size[0] = i;
  emxEnsureCapacity((emxArray__common *)r8, i3, (int)sizeof(int));
  for (i3 = 0; i3 < i; i3++) {
    r8->data[i3] = i3;
  }

  emxInit_real_T1(&b_fch3, 1);
  i3 = b_fch3->size[0];
  b_fch3->size[0] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch3, i3, (int)sizeof(double));
  for (i3 = 0; i3 < fch1_idx_0; i3++) {
    b_fch3->data[i3] = fch3->data[i3];
  }

  i = r8->size[0];
  for (i3 = 0; i3 < i; i3++) {
    fchw->data[2 + fchw->size[0] * r8->data[i3]] = b_fch3->data[i3];
  }

  emxFree_real_T(&b_fch3);
  emxFree_int32_T(&r8);

  //  all windows should be the same size:
  memset(&FFT_Ltop[0], 0, sizeof(double) << 3);
  for (i = 0; i < 4; i++) {
    FFT_PkRatio[i] = 0.0;
  }

  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  emxInit_real_T(&PSD, 2);
  if (b_mod((double)fch1_idx_0, 2.0) == 1.0) {
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
    y = ((double)fch1_idx_0 - 1.0) / 2.0;
    i3 = PSD->size[0] * PSD->size[1];
    PSD->size[0] = 4;
    PSD->size[1] = (int)y;
    emxEnsureCapacity((emxArray__common *)PSD, i3, (int)sizeof(double));
    i = (int)y << 2;
    for (i3 = 0; i3 < i; i3++) {
      PSD->data[i3] = 0.0;
    }
  } else {
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
    y = (double)fch1_idx_0 / 2.0;
    i3 = PSD->size[0] * PSD->size[1];
    PSD->size[0] = 4;
    PSD->size[1] = (int)y;
    emxEnsureCapacity((emxArray__common *)PSD, i3, (int)sizeof(double));
    i = (int)y << 2;
    for (i3 = 0; i3 < i; i3++) {
      PSD->data[i3] = 0.0;
    }
  }

  memset(&PSD_Ltop[0], 0, sizeof(double) << 3);
  for (i = 0; i < 4; i++) {
    PSD_PkRatio[i] = 0.0;
  }

  emxInit_real_T1(&hW, 1);
  emxInit_real_T(&fPSD, 2);

  // 0?
  //  Data is already filtered:
  // between 250?500dp
  //  Preallocate for spd:
  emxInit_real_T(&FFT_PKS, 2);
  emxInit_real_T(&FFT_L, 2);
  emxInit_real_T(&b_PSD, 2);
  emxInit_real_T(&b_fchw, 2);
  emxInit_real_T(&c_fchw, 2);
  for (ch = 0; ch < 3; ch++) {
    //  #1 Take FFT:
    i = fchw->size[1];
    i3 = c_fchw->size[0] * c_fchw->size[1];
    c_fchw->size[0] = 1;
    c_fchw->size[1] = i;
    emxEnsureCapacity((emxArray__common *)c_fchw, i3, (int)sizeof(double));
    for (i3 = 0; i3 < i; i3++) {
      c_fchw->data[c_fchw->size[0] * i3] = fchw->data[ch + fchw->size[0] * i3];
    }

    get_nfft_data(c_fchw, Fs, f, dv16);

    //  #1.1 Find Peaks and M/I
    for (i3 = 0; i3 < 1025; i3++) {
      FFT[ch + (i3 << 2)] = dv16[i3];
      b_FFT[i3] = FFT[ch + (i3 << 2)];
    }

    b_findpeaks(b_FFT, FFT_PKS, FFT_L);
    if (FFT_PKS->size[1] > 1) {
      // Peak max minus min
      for (i3 = 0; i3 < 2; i3++) {
        FFT_Ltop[ch + (i3 << 2)] = f[(int)FFT_L->data[i3] - 1];
      }

      w = 0;
      exitg4 = false;
      while ((!exitg4) && (w < 4)) {
        if ((FFT_Ltop[ch] > threshFFT[w]) && (FFT_Ltop[ch] < threshFFT[4 + w]))
        {
          FFT_PkRatio[ch] = FFT_PKS->data[0] / FFT_PKS->data[1];
          wLFFT[ch] = (signed char)(1 + w);
          exitg4 = true;
        } else {
          FFT_PkRatio[ch] = 0.0;
          wLFFT[ch] = 0;
          w++;
        }
      }
    }

    //  #2 Take PSD Estimate: (Welch method)
    //  Prepare hanning window:
    fch1_idx_0 = fch1->size[0] * fch1->size[1];
    hannWin((double)fch1_idx_0, hW);
    i = fchw->size[1];
    i3 = b_fchw->size[0] * b_fchw->size[1];
    b_fchw->size[0] = 1;
    b_fchw->size[1] = i;
    emxEnsureCapacity((emxArray__common *)b_fchw, i3, (int)sizeof(double));
    for (i3 = 0; i3 < i; i3++) {
      b_fchw->data[b_fchw->size[0] * i3] = fchw->data[ch + fchw->size[0] * i3];
    }

    welch_psd(b_fchw, Fs, hW, FFT_PKS, fPSD);
    i = FFT_PKS->size[1];
    for (i3 = 0; i3 < i; i3++) {
      PSD->data[ch + PSD->size[0] * i3] = FFT_PKS->data[FFT_PKS->size[0] * i3];
    }

    // fin-start
    //  #2.2 Find Peaks and Max
    i = PSD->size[1];
    i3 = b_PSD->size[0] * b_PSD->size[1];
    b_PSD->size[0] = 1;
    b_PSD->size[1] = i;
    emxEnsureCapacity((emxArray__common *)b_PSD, i3, (int)sizeof(double));
    for (i3 = 0; i3 < i; i3++) {
      b_PSD->data[b_PSD->size[0] * i3] = PSD->data[ch + PSD->size[0] * i3];
    }

    c_findpeaks(b_PSD, FFT_PKS, FFT_L);
    if (FFT_PKS->size[1] > 1) {
      for (i3 = 0; i3 < 2; i3++) {
        PSD_Ltop[ch + (i3 << 2)] = fPSD->data[(int)FFT_L->data[i3] - 1];
      }

      w = 0;
      exitg3 = false;
      while ((!exitg3) && (w < 4)) {
        if ((PSD_Ltop[ch] >= threshPSD[w]) && (PSD_Ltop[ch] <= threshPSD[4 + w]))
        {
          PSD_PkRatio[ch] = FFT_PKS->data[0] / FFT_PKS->data[1];
          wLPSD[ch] = (signed char)(1 + w);
          exitg3 = true;
        } else {
          PSD_PkRatio[ch] = 0.0;
          wLPSD[ch] = 0;
          w++;
        }
      }
    }
  }

  emxFree_real_T(&c_fchw);
  emxFree_real_T(&b_fchw);
  emxFree_real_T(&b_PSD);
  emxFree_real_T(&hW);
  emxFree_real_T(&fchw);

  // Combine data into 'fourth' channel:
  for (i3 = 0; i3 < 1025; i3++) {
    FFT[3 + (i3 << 2)] = (FFT[i3 << 2] + FFT[1 + (i3 << 2)]) + FFT[2 + (i3 << 2)];
    b_FFT[i3] = FFT[3 + (i3 << 2)];
  }

  b_findpeaks(b_FFT, FFT_PKS, FFT_L);
  if (FFT_PKS->size[1] > 1) {
    for (i3 = 0; i3 < 2; i3++) {
      FFT_Ltop[3 + (i3 << 2)] = f[(int)FFT_L->data[i3] - 1];
    }

    w = 0;
    exitg2 = false;
    while ((!exitg2) && (w < 4)) {
      if ((FFT_Ltop[3] > threshFFT[w]) && (FFT_Ltop[3] < threshFFT[4 + w])) {
        FFT_PkRatio[3] = FFT_PKS->data[0] / FFT_PKS->data[1];
        wLFFT[3] = (signed char)(1 + w);
        exitg2 = true;
      } else {
        FFT_PkRatio[3] = 0.0;
        wLFFT[3] = 0;
        w++;
      }
    }
  }

  emxInit_real_T(&c_PSD, 2);
  i = PSD->size[1];
  i3 = c_PSD->size[0] * c_PSD->size[1];
  c_PSD->size[0] = 1;
  c_PSD->size[1] = i;
  emxEnsureCapacity((emxArray__common *)c_PSD, i3, (int)sizeof(double));
  for (i3 = 0; i3 < i; i3++) {
    c_PSD->data[c_PSD->size[0] * i3] = (PSD->data[PSD->size[0] * i3] + PSD->
      data[1 + PSD->size[0] * i3]) + PSD->data[2 + PSD->size[0] * i3];
  }

  i = c_PSD->size[1];
  for (i3 = 0; i3 < i; i3++) {
    PSD->data[3 + PSD->size[0] * i3] = c_PSD->data[c_PSD->size[0] * i3];
  }

  emxFree_real_T(&c_PSD);
  emxInit_real_T(&d_PSD, 2);
  i = PSD->size[1];
  i3 = d_PSD->size[0] * d_PSD->size[1];
  d_PSD->size[0] = 1;
  d_PSD->size[1] = i;
  emxEnsureCapacity((emxArray__common *)d_PSD, i3, (int)sizeof(double));
  for (i3 = 0; i3 < i; i3++) {
    d_PSD->data[d_PSD->size[0] * i3] = PSD->data[3 + PSD->size[0] * i3];
  }

  emxFree_real_T(&PSD);
  c_findpeaks(d_PSD, FFT_PKS, FFT_L);
  emxFree_real_T(&d_PSD);
  if (FFT_PKS->size[1] > 1) {
    for (i3 = 0; i3 < 2; i3++) {
      PSD_Ltop[3 + (i3 << 2)] = fPSD->data[(int)FFT_L->data[i3] - 1];
    }

    w = 0;
    exitg1 = false;
    while ((!exitg1) && (w < 4)) {
      if ((PSD_Ltop[3] >= threshPSD[w]) && (PSD_Ltop[3] <= threshPSD[4 + w])) {
        PSD_PkRatio[3] = FFT_PKS->data[0] / FFT_PKS->data[1];
        wLPSD[3] = (signed char)(1 + w);
        exitg1 = true;
      } else {
        PSD_PkRatio[3] = 0.0;
        wLPSD[3] = 0;
        w++;
      }
    }
  }

  emxFree_real_T(&FFT_L);
  emxFree_real_T(&FFT_PKS);
  emxFree_real_T(&fPSD);
  for (chn = 0; chn < 4; chn++) {
    FFTPeaks1[chn] = FFT_Ltop[chn];
    PSDPeaks1[chn] = PSD_Ltop[chn];
  }

  b_FFT_Ltop[0] = FFT_Ltop[0];
  b_FFT_Ltop[1] = FFT_Ltop[1];
  b_FFT_Ltop[2] = FFT_Ltop[2];
  b_FFT_Ltop[3] = FFT_Ltop[3];
  if ((wLFFT[0] != 0) && (wLFFT[1] != 0) && (wLFFT[2] != 0) && (wLFFT[3] != 0))
  {
    b0 = true;
  } else {
    b0 = false;
  }

  if (b0) {
    // if a signal was detected on each FFT
    // check that they are all the same;
    b_wLFFT[0] = wLFFT[0];
    b_wLFFT[1] = wLFFT[1];
    b_wLFFT[2] = wLFFT[2];
    b_wLFFT[3] = wLFFT[3];
    c_wLFFT[0] = wLFFT[0];
    c_wLFFT[1] = wLFFT[0];
    c_wLFFT[2] = wLFFT[0];
    c_wLFFT[3] = wLFFT[0];
    b2 = isequal(b_wLFFT, c_wLFFT);
  } else {
    b2 = false;
  }

  b_PSD_Ltop[0] = PSD_Ltop[0];
  b_PSD_Ltop[1] = PSD_Ltop[1];
  b_PSD_Ltop[2] = PSD_Ltop[2];
  b_PSD_Ltop[3] = PSD_Ltop[3];
  if ((wLPSD[0] != 0) && (wLPSD[1] != 0) && (wLPSD[2] != 0) && (wLPSD[3] != 0))
  {
    b1 = true;
  } else {
    b1 = false;
  }

  if (b1) {
    // if a signal was detected for PSD on all channels:
    // check they are equivalent:
    b_wLPSD[0] = wLPSD[0];
    b_wLPSD[1] = wLPSD[1];
    b_wLPSD[2] = wLPSD[2];
    b_wLPSD[3] = wLPSD[3];
    c_wLPSD[0] = wLPSD[0];
    c_wLPSD[1] = wLPSD[0];
    c_wLPSD[2] = wLPSD[0];
    c_wLPSD[3] = wLPSD[0];
    b4 = isequal(b_wLPSD, c_wLPSD);
  } else {
    b4 = false;
  }

  // % Collect Feature data into 'F'
  // First separate features by channel: (row vects)
  //  first to remove: *FFT_Ltop(2) ... not sure how I will use this
  //  Also remove FFTPeaks2 and averageFFTPeak2
  // WANT INFO TO PRINT IN ORDER:
  //     %% FPRINTFs:
  F[16] = mean(b_FFT_Ltop);
  F[17] = mean(b_PSD_Ltop);
  for (i3 = 0; i3 < 4; i3++) {
    F[i3] = wLFFT[i3];
    F[i3 + 4] = wLPSD[i3];
    F[i3 + 8] = FFT_PkRatio[i3];
    F[i3 + 12] = PSD_PkRatio[i3];
    F[i3 + 18] = FFTPeaks1[i3];
    F[i3 + 22] = PSDPeaks1[i3];
  }

  F[26] = b0;
  F[27] = b2;
  F[28] = b1;
  F[29] = b4;

  // END FUNCTION
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_creal_T *y
// Return Type  : void
//
static void fft(const emxArray_real_T *x, emxArray_creal_T *y)
{
  int n1;
  boolean_T useRadix2;
  int pmin;
  int pmax;
  int nn1m1;
  emxArray_real_T *costab1q;
  double e;
  int nRowsD4;
  boolean_T exitg1;
  int n;
  int pow2p;
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  emxArray_real_T *sintabinv;
  int nRowsD2;
  int i;
  double temp_re;
  double temp_im;
  double twid_im;
  n1 = x->size[0];
  if (x->size[0] == 0) {
    pmax = y->size[0];
    y->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)y, pmax, (int)sizeof(creal_T));
  } else {
    useRadix2 = ((x->size[0] & (x->size[0] - 1)) == 0);
    pmin = 1;
    if (useRadix2) {
      nn1m1 = x->size[0];
    } else {
      nn1m1 = (x->size[0] + x->size[0]) - 1;
      pmax = 31;
      if (nn1m1 > MIN_int32_T) {
        if (nn1m1 < 0) {
          nn1m1 = -nn1m1;
        }

        if (nn1m1 <= 1) {
          pmax = 0;
        } else {
          pmin = 0;
          exitg1 = false;
          while ((!exitg1) && (pmax - pmin > 1)) {
            n = (pmin + pmax) >> 1;
            pow2p = 1 << n;
            if (pow2p == nn1m1) {
              pmax = n;
              exitg1 = true;
            } else if (pow2p > nn1m1) {
              pmax = n;
            } else {
              pmin = n;
            }
          }
        }
      }

      pmin = 1 << pmax;
      nn1m1 = pmin;
    }

    emxInit_real_T(&costab1q, 2);
    e = 6.2831853071795862 / (double)nn1m1;
    nRowsD4 = nn1m1 / 2 / 2;
    pmax = costab1q->size[0] * costab1q->size[1];
    costab1q->size[0] = 1;
    costab1q->size[1] = nRowsD4 + 1;
    emxEnsureCapacity((emxArray__common *)costab1q, pmax, (int)sizeof(double));
    costab1q->data[0] = 1.0;
    nn1m1 = nRowsD4 / 2;
    for (pmax = 1; pmax <= nn1m1; pmax++) {
      costab1q->data[pmax] = std::cos(e * (double)pmax);
    }

    for (pmax = nn1m1 + 1; pmax < nRowsD4; pmax++) {
      costab1q->data[pmax] = std::sin(e * (double)(nRowsD4 - pmax));
    }

    costab1q->data[nRowsD4] = 0.0;
    emxInit_real_T(&costab, 2);
    emxInit_real_T(&sintab, 2);
    emxInit_real_T(&sintabinv, 2);
    if (!useRadix2) {
      n = costab1q->size[1] - 1;
      nn1m1 = (costab1q->size[1] - 1) << 1;
      pmax = costab->size[0] * costab->size[1];
      costab->size[0] = 1;
      costab->size[1] = nn1m1 + 1;
      emxEnsureCapacity((emxArray__common *)costab, pmax, (int)sizeof(double));
      pmax = sintab->size[0] * sintab->size[1];
      sintab->size[0] = 1;
      sintab->size[1] = nn1m1 + 1;
      emxEnsureCapacity((emxArray__common *)sintab, pmax, (int)sizeof(double));
      costab->data[0] = 1.0;
      sintab->data[0] = 0.0;
      pmax = sintabinv->size[0] * sintabinv->size[1];
      sintabinv->size[0] = 1;
      sintabinv->size[1] = nn1m1 + 1;
      emxEnsureCapacity((emxArray__common *)sintabinv, pmax, (int)sizeof(double));
      for (pmax = 1; pmax <= n; pmax++) {
        sintabinv->data[pmax] = costab1q->data[n - pmax];
      }

      for (pmax = costab1q->size[1]; pmax <= nn1m1; pmax++) {
        sintabinv->data[pmax] = costab1q->data[pmax - n];
      }

      for (pmax = 1; pmax <= n; pmax++) {
        costab->data[pmax] = costab1q->data[pmax];
        sintab->data[pmax] = -costab1q->data[n - pmax];
      }

      for (pmax = costab1q->size[1]; pmax <= nn1m1; pmax++) {
        costab->data[pmax] = -costab1q->data[nn1m1 - pmax];
        sintab->data[pmax] = -costab1q->data[pmax - n];
      }
    } else {
      n = costab1q->size[1] - 1;
      nn1m1 = (costab1q->size[1] - 1) << 1;
      pmax = costab->size[0] * costab->size[1];
      costab->size[0] = 1;
      costab->size[1] = nn1m1 + 1;
      emxEnsureCapacity((emxArray__common *)costab, pmax, (int)sizeof(double));
      pmax = sintab->size[0] * sintab->size[1];
      sintab->size[0] = 1;
      sintab->size[1] = nn1m1 + 1;
      emxEnsureCapacity((emxArray__common *)sintab, pmax, (int)sizeof(double));
      costab->data[0] = 1.0;
      sintab->data[0] = 0.0;
      for (pmax = 1; pmax <= n; pmax++) {
        costab->data[pmax] = costab1q->data[pmax];
        sintab->data[pmax] = -costab1q->data[n - pmax];
      }

      for (pmax = costab1q->size[1]; pmax <= nn1m1; pmax++) {
        costab->data[pmax] = -costab1q->data[nn1m1 - pmax];
        sintab->data[pmax] = -costab1q->data[pmax - n];
      }

      pmax = sintabinv->size[0] * sintabinv->size[1];
      sintabinv->size[0] = 1;
      sintabinv->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)sintabinv, pmax, (int)sizeof(double));
    }

    emxFree_real_T(&costab1q);
    if (useRadix2) {
      pow2p = x->size[0];
      nRowsD2 = x->size[0] / 2;
      nRowsD4 = nRowsD2 / 2;
      nn1m1 = x->size[0];
      pmax = y->size[0];
      y->size[0] = nn1m1;
      emxEnsureCapacity((emxArray__common *)y, pmax, (int)sizeof(creal_T));
      pmax = 0;
      pmin = 0;
      nn1m1 = 0;
      for (i = 1; i < pow2p; i++) {
        y->data[nn1m1].re = x->data[pmax];
        y->data[nn1m1].im = 0.0;
        n = n1;
        useRadix2 = true;
        while (useRadix2) {
          n >>= 1;
          pmin ^= n;
          useRadix2 = ((pmin & n) == 0);
        }

        nn1m1 = pmin;
        pmax++;
      }

      y->data[nn1m1].re = x->data[pmax];
      y->data[nn1m1].im = 0.0;
      if (x->size[0] > 1) {
        for (i = 0; i <= n1 - 2; i += 2) {
          temp_re = y->data[i + 1].re;
          temp_im = y->data[i + 1].im;
          y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
          y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
        }
      }

      nn1m1 = 2;
      pmax = 4;
      pmin = 1 + ((nRowsD4 - 1) << 2);
      while (nRowsD4 > 0) {
        for (i = 0; i < pmin; i += pmax) {
          temp_re = y->data[i + nn1m1].re;
          temp_im = y->data[i + nn1m1].im;
          y->data[i + nn1m1].re = y->data[i].re - temp_re;
          y->data[i + nn1m1].im = y->data[i].im - temp_im;
          y->data[i].re += temp_re;
          y->data[i].im += temp_im;
        }

        n = 1;
        for (pow2p = nRowsD4; pow2p < nRowsD2; pow2p += nRowsD4) {
          e = costab->data[pow2p];
          twid_im = sintab->data[pow2p];
          i = n;
          n1 = n + pmin;
          while (i < n1) {
            temp_re = e * y->data[i + nn1m1].re - twid_im * y->data[i + nn1m1].
              im;
            temp_im = e * y->data[i + nn1m1].im + twid_im * y->data[i + nn1m1].
              re;
            y->data[i + nn1m1].re = y->data[i].re - temp_re;
            y->data[i + nn1m1].im = y->data[i].im - temp_im;
            y->data[i].re += temp_re;
            y->data[i].im += temp_im;
            i += pmax;
          }

          n++;
        }

        nRowsD4 /= 2;
        nn1m1 = pmax;
        pmax <<= 1;
        pmin -= nn1m1;
      }
    } else {
      dobluesteinfft(x, pmin, x->size[0], costab, sintab, sintabinv, y);
    }

    emxFree_real_T(&sintabinv);
    emxFree_real_T(&sintab);
    emxFree_real_T(&costab);
  }
}

//
// Arguments    : double b[4]
//                double a[4]
//                const double x[268]
//                const double zi[3]
//                double y[268]
// Return Type  : void
//
static void filter(double b[4], double a[4], const double x[268], const double
                   zi[3], double y[268])
{
  double a1;
  int k;
  double dbuffer[4];
  int j;
  a1 = a[0];
  if ((!((!rtIsInf(a[0])) && (!rtIsNaN(a[0])))) || (a[0] == 0.0) || (!(a[0] !=
        1.0))) {
  } else {
    for (k = 0; k < 4; k++) {
      b[k] /= a1;
    }

    for (k = 0; k < 3; k++) {
      a[k + 1] /= a1;
    }

    a[0] = 1.0;
  }

  for (k = 0; k < 3; k++) {
    dbuffer[k + 1] = zi[k];
  }

  for (j = 0; j < 268; j++) {
    for (k = 0; k < 3; k++) {
      dbuffer[k] = dbuffer[k + 1];
    }

    dbuffer[3] = 0.0;
    for (k = 0; k < 4; k++) {
      dbuffer[k] += x[j] * b[k];
    }

    for (k = 0; k < 3; k++) {
      dbuffer[k + 1] -= dbuffer[0] * a[k + 1];
    }

    y[j] = dbuffer[0];
  }
}

//
// Arguments    : const double x_in[250]
//                double y_out[250]
// Return Type  : void
//
static void filtfilt(const double x_in[250], double y_out[250])
{
  double xtmp;
  double d0;
  int i;
  double y[268];
  double dv6[4];
  double dv7[4];
  double a[3];
  static const double dv8[4] = { 0.00156701035058832, 0.00470103105176495,
    0.00470103105176495, 0.00156701035058832 };

  static const double dv9[4] = { 1.0, -2.49860834469118, 2.11525412700316,
    -0.604109699507275 };

  double b_y[268];
  static const double b_a[3] = { 0.99843298964950811, -1.5048763860936776,
    0.60567670985792155 };

  double c_y[268];
  xtmp = 2.0 * x_in[0];
  d0 = 2.0 * x_in[249];
  for (i = 0; i < 9; i++) {
    y[i] = xtmp - x_in[9 - i];
  }

  memcpy(&y[9], &x_in[0], 250U * sizeof(double));
  for (i = 0; i < 9; i++) {
    y[i + 259] = d0 - x_in[248 - i];
  }

  for (i = 0; i < 4; i++) {
    dv6[i] = dv8[i];
    dv7[i] = dv9[i];
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 268U * sizeof(double));
  filter(dv6, dv7, b_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  for (i = 0; i < 4; i++) {
    dv6[i] = dv8[i];
    dv7[i] = dv9[i];
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&c_y[0], &y[0], 268U * sizeof(double));
  filter(dv6, dv7, c_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  memcpy(&y_out[0], &y[9], 250U * sizeof(double));
}

//
// Arguments    : const double yTemp[250]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void findLocalMaxima(const double yTemp[250], emxArray_real_T *iPk,
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
  emxArray_int32_T *r5;
  boolean_T guard3 = false;
  emxArray_real_T *iTemp;
  emxArray_real_T *c_yTemp;
  emxArray_real_T *s;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r6;
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

  emxInit_int32_T(&r5, 1);
  i1 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  i1 = r5->size[0];
  r5->size[0] = 1 + b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r5, i1, (int)sizeof(int));
  r5->data[0] = 1;
  ii = b_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    r5->data[i1 + 1] = b_ii->data[i1] + 1;
  }

  emxInit_real_T1(&iTemp, 1);
  i1 = iTemp->size[0];
  iTemp->size[0] = r5->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i1, (int)sizeof(double));
  ii = r5->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iTemp->data[i1] = 1.0 + (double)(r5->data[i1] - 1);
  }

  emxFree_int32_T(&r5);
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
  emxInit_real_T1(&r6, 1);
  diff(c_yTemp, s);
  b_sign(s);
  diff(s, r6);
  i1 = b_x->size[0];
  b_x->size[0] = r6->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i1, (int)sizeof(boolean_T));
  ii = r6->size[0];
  emxFree_real_T(&c_yTemp);
  for (i1 = 0; i1 < ii; i1++) {
    b_x->data[i1] = (r6->data[i1] < 0.0);
  }

  emxFree_real_T(&r6);
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
// Arguments    : double x
//                const emxArray_real_T *bin_edges
// Return Type  : int
//
static int findbin(double x, const emxArray_real_T *bin_edges)
{
  int k;
  int low_ip1;
  int high_i;
  int mid_i;
  k = 0;
  if ((!(bin_edges->size[1] == 0)) && (!rtIsNaN(x))) {
    if ((x >= bin_edges->data[0]) && (x < bin_edges->data[bin_edges->size[1] - 1]))
    {
      k = 1;
      low_ip1 = 2;
      high_i = bin_edges->size[1];
      while (high_i > low_ip1) {
        mid_i = (k >> 1) + (high_i >> 1);
        if (((k & 1) == 1) && ((high_i & 1) == 1)) {
          mid_i++;
        }

        if (x >= bin_edges->data[mid_i - 1]) {
          k = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }
    }

    if (x == bin_edges->data[bin_edges->size[1] - 1]) {
      k = bin_edges->size[1];
    }
  }

  return k;
}

//
// Arguments    : const double Yin[250]
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void findpeaks(const double Yin[250], emxArray_real_T *Ypk,
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
  findLocalMaxima(yTemp, iPk, b_idx);
  if (!(iPk->size[0] == 0)) {
    cdiff = iPk->size[0] - 1;
    ndbl = 0;
    for (idx = 0; idx <= cdiff; idx++) {
      if (Yin[(int)iPk->data[idx] - 1] > 6.5E-5) {
        ndbl++;
      }
    }

    k = 0;
    for (idx = 0; idx <= cdiff; idx++) {
      if (Yin[(int)iPk->data[idx] - 1] > 6.5E-5) {
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
// Arguments    : const emxArray_real_T *y
//                emxArray_real_T *iPk
//                emxArray_real_T *iInf
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void getAllPeaks(const emxArray_real_T *y, emxArray_real_T *iPk,
  emxArray_real_T *iInf, emxArray_real_T *iInflect)
{
  emxArray_boolean_T *x;
  int i6;
  int ii;
  emxArray_int32_T *b_ii;
  int nx;
  int idx;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxArray_real_T *yTemp;
  emxInit_boolean_T(&x, 1);
  i6 = x->size[0];
  x->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)x, i6, (int)sizeof(boolean_T));
  ii = y->size[0];
  for (i6 = 0; i6 < ii; i6++) {
    x->data[i6] = rtIsInf(y->data[i6]);
  }

  i6 = x->size[0];
  emxEnsureCapacity((emxArray__common *)x, i6, (int)sizeof(boolean_T));
  ii = x->size[0];
  for (i6 = 0; i6 < ii; i6++) {
    x->data[i6] = (x->data[i6] && (y->data[i6] > 0.0));
  }

  emxInit_int32_T(&b_ii, 1);
  nx = x->size[0];
  idx = 0;
  i6 = b_ii->size[0];
  b_ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i6, (int)sizeof(int));
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx)) {
    guard1 = false;
    if (x->data[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
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

  if (x->size[0] == 1) {
    if (idx == 0) {
      i6 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i6, (int)sizeof(int));
    }
  } else {
    i6 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i6, (int)sizeof(int));
  }

  emxFree_boolean_T(&x);
  i6 = iInf->size[0];
  iInf->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInf, i6, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i6 = 0; i6 < ii; i6++) {
    iInf->data[i6] = b_ii->data[i6];
  }

  emxInit_real_T1(&yTemp, 1);
  i6 = yTemp->size[0];
  yTemp->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)yTemp, i6, (int)sizeof(double));
  ii = y->size[0];
  for (i6 = 0; i6 < ii; i6++) {
    yTemp->data[i6] = y->data[i6];
  }

  i6 = b_ii->size[0];
  b_ii->size[0] = iInf->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i6, (int)sizeof(int));
  ii = iInf->size[0];
  for (i6 = 0; i6 < ii; i6++) {
    b_ii->data[i6] = (int)iInf->data[i6];
  }

  ii = b_ii->size[0];
  for (i6 = 0; i6 < ii; i6++) {
    yTemp->data[b_ii->data[i6] - 1] = rtNaN;
  }

  emxFree_int32_T(&b_ii);
  c_findLocalMaxima(yTemp, iPk, iInflect);
  emxFree_real_T(&yTemp);
}

//
// get_fft_data:
//  X is filtered data
//  L = size(X,1);
//  L = number of FFT points
// Arguments    : const emxArray_real_T *X
//                double Fs
//                double f[1025]
//                double C[1025]
// Return Type  : void
//
static void get_nfft_data(const emxArray_real_T *X, double Fs, double f[1025],
  double C[1025])
{
  int iy;
  emxArray_real_T *x;
  creal_T b_y1[2048];
  int iDelta2;
  int i;
  creal_T c_y1[2048];
  int k;
  double B[2048];
  int ix;
  int ju;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int j;
  double twid_re;
  static const double dv17[1025] = { 1.0, 0.99999529380957619,
    0.99998117528260111, 0.9999576445519639, 0.9999247018391445,
    0.99988234745421256, 0.9998305817958234, 0.99976940535121528,
    0.99969881869620425, 0.99961882249517864, 0.99952941750109314,
    0.99943060455546173, 0.99932238458834954, 0.99920475861836389,
    0.99907772775264536, 0.99894129318685687, 0.99879545620517241,
    0.99864021818026527, 0.99847558057329477, 0.99830154493389289,
    0.99811811290014918, 0.997925286198596, 0.99772306664419164,
    0.99751145614030345, 0.99729045667869021, 0.997060070339483,
    0.99682029929116567, 0.99657114579055484, 0.996312612182778,
    0.996044700901252, 0.99576741446765982, 0.99548075549192694,
    0.99518472667219693, 0.99487933079480562, 0.99456457073425542,
    0.9942404494531879, 0.99390697000235606, 0.9935641355205953,
    0.9932119492347945, 0.9928504144598651, 0.99247953459871, 0.9920993131421918,
    0.99170975366909953, 0.99131085984611544, 0.99090263542778,
    0.99048508425645709, 0.99005821026229712, 0.98962201746320089,
    0.989176509964781, 0.98872169196032378, 0.98825756773074946,
    0.98778414164457218, 0.98730141815785843, 0.98680940181418553,
    0.98630809724459867, 0.98579750916756748, 0.98527764238894122,
    0.98474850180190421, 0.984210092386929, 0.98366241921173025,
    0.98310548743121629, 0.98253930228744124, 0.98196386910955524,
    0.98137919331375456, 0.98078528040323043, 0.98018213596811743,
    0.97956976568544052, 0.9789481753190622, 0.97831737071962765,
    0.97767735782450993, 0.97702814265775439, 0.97636973133002114,
    0.97570213003852857, 0.97502534506699412, 0.97433938278557586,
    0.973644249650812, 0.97293995220556018, 0.97222649707893627,
    0.97150389098625178, 0.97077214072895035, 0.970031253194544,
    0.96928123535654853, 0.96852209427441727, 0.96775383709347551,
    0.96697647104485207, 0.9661900034454125, 0.9653944416976894,
    0.96458979328981276, 0.96377606579543984, 0.96295326687368388,
    0.96212140426904158, 0.96128048581132064, 0.96043051941556579,
    0.95957151308198452, 0.9587034748958716, 0.95782641302753291,
    0.95694033573220882, 0.95604525134999641, 0.95514116830577078,
    0.95422809510910567, 0.95330604035419386, 0.95237501271976588,
    0.95143502096900834, 0.9504860739494817, 0.94952818059303667,
    0.94856134991573027, 0.94758559101774109, 0.94660091308328353,
    0.94560732538052128, 0.94460483726148026, 0.94359345816196039,
    0.94257319760144687, 0.94154406518302081, 0.9405060705932683,
    0.93945922360218992, 0.93840353406310806, 0.937339011912575,
    0.93626566717027826, 0.93518350993894761, 0.93409255040425887,
    0.932992798834739, 0.93188426558166815, 0.93076696107898371,
    0.92964089584318121, 0.92850608047321559, 0.92736252565040111,
    0.92621024213831138, 0.92504924078267758, 0.92387953251128674,
    0.92270112833387863, 0.9215140393420419, 0.92031827670911059,
    0.91911385169005777, 0.9179007756213905, 0.9166790599210427,
    0.91544871608826783, 0.91420975570353069, 0.91296219042839821,
    0.91170603200542988, 0.91044129225806725, 0.90916798309052238,
    0.90788611648766626, 0.90659570451491533, 0.90529675931811882,
    0.90398929312344334, 0.90267331823725883, 0.901348847046022,
    0.90001589201616017, 0.89867446569395382, 0.89732458070541832,
    0.89596624975618522, 0.8945994856313827, 0.89322430119551532,
    0.89184070939234272, 0.89044872324475788, 0.88904835585466457,
    0.88763962040285393, 0.88622253014888064, 0.88479709843093779,
    0.88336333866573158, 0.881921264348355, 0.88047088905216075,
    0.87901222642863353, 0.87754529020726135, 0.8760700941954066,
    0.87458665227817611, 0.87309497841829009, 0.87159508665595109,
    0.87008699110871146, 0.8685707059713409, 0.86704624551569265,
    0.86551362409056909, 0.8639728561215867, 0.8624239561110405,
    0.86086693863776731, 0.85930181835700847, 0.85772861000027212,
    0.85614732837519447, 0.85455798836540053, 0.85296060493036363,
    0.8513551931052652, 0.84974176800085255, 0.84812034480329723,
    0.84649093877405213, 0.84485356524970712, 0.84320823964184544,
    0.84155497743689844, 0.83989379419599952, 0.83822470555483808,
    0.836547727223512, 0.83486287498638, 0.83317016470191319,
    0.83146961230254524, 0.829761233794523, 0.8280450452577558,
    0.82632106284566353, 0.82458930278502529, 0.82284978137582643,
    0.82110251499110465, 0.819347520076797, 0.81758481315158371,
    0.81581441080673378, 0.81403632970594841, 0.81225058658520388,
    0.81045719825259477, 0.808656181588175, 0.80684755354379933,
    0.80503133114296366, 0.80320753148064494, 0.80137617172314024,
    0.799537269107905, 0.79769084094339116, 0.79583690460888357,
    0.79397547755433717, 0.79210657730021239, 0.79023022143731,
    0.78834642762660634, 0.78645521359908577, 0.78455659715557524,
    0.78265059616657573, 0.78073722857209449, 0.778816512381476,
    0.77688846567323244, 0.77495310659487393, 0.773010453362737,
    0.77106052426181382, 0.7691033376455797, 0.7671389119358204,
    0.765167265622459, 0.76318841726338127, 0.76120238548426178,
    0.759209188978388, 0.75720884650648457, 0.75520137689653655,
    0.75318679904361252, 0.75116513190968637, 0.74913639452345937,
    0.74710060598018013, 0.745057785441466, 0.74300795213512172,
    0.74095112535495922, 0.73888732446061511, 0.73681656887736979,
    0.7347388780959635, 0.73265427167241282, 0.73056276922782759,
    0.7284643904482252, 0.726359155084346, 0.724247082951467,
    0.72212819392921535, 0.72000250796138165, 0.71787004505573171,
    0.71573082528381859, 0.71358486878079352, 0.71143219574521643,
    0.70927282643886569, 0.70710678118654757, 0.70493408037590488,
    0.7027547444572253, 0.70056879394324834, 0.69837624940897292,
    0.696177131491463, 0.69397146088965389, 0.69175925836415775,
    0.68954054473706683, 0.687315340891759, 0.68508366777270036,
    0.68284554638524808, 0.680600997795453, 0.67835004312986147,
    0.67609270357531592, 0.673829000378756, 0.67155895484701833,
    0.669282588346636, 0.66699992230363747, 0.66471097820334479,
    0.66241577759017178, 0.66011434206742048, 0.65780669329707864,
    0.65549285299961535, 0.65317284295377676, 0.650846684996381,
    0.64851440102211244, 0.64617601298331628, 0.64383154288979139,
    0.641481012808583, 0.63912444486377573, 0.6367618612362842,
    0.63439328416364549, 0.63201873593980906, 0.629638238914927,
    0.62725181549514408, 0.62485948814238634, 0.62246127937415,
    0.6200572117632891, 0.61764730793780387, 0.61523159058062682,
    0.61281008242940971, 0.61038280627630948, 0.60794978496777363,
    0.60551104140432555, 0.60306659854034816, 0.600616479383869,
    0.59816070699634238, 0.59569930449243336, 0.5932322950397998,
    0.59075970185887416, 0.58828154822264522, 0.58579785745643886,
    0.58330865293769829, 0.58081395809576453, 0.57831379641165559,
    0.57580819141784534, 0.5732971666980422, 0.57078074588696726,
    0.56825895267013149, 0.56573181078361312, 0.56319934401383409,
    0.560661576197336, 0.5581185312205561, 0.55557023301960218,
    0.55301670558002747, 0.55045797293660481, 0.54789405917310019,
    0.54532498842204646, 0.54275078486451589, 0.54017147272989285,
    0.53758707629564539, 0.53499761988709715, 0.5324031278771979,
    0.52980362468629461, 0.52719913478190128, 0.524589682678469,
    0.52197529293715439, 0.51935599016558964, 0.51673179901764987,
    0.51410274419322166, 0.5114688504379703, 0.508830142543107,
    0.50618664534515523, 0.50353838372571758, 0.50088538261124071,
    0.49822766697278181, 0.49556526182577254, 0.49289819222978404,
    0.49022648328829116, 0.487550160148436, 0.48486924800079106,
    0.48218377207912272, 0.47949375766015295, 0.47679923006332209,
    0.47410021465054997, 0.47139673682599764, 0.46868882203582796,
    0.46597649576796618, 0.46325978355186015, 0.46053871095824,
    0.45781330359887717, 0.45508358712634384, 0.45234958723377089,
    0.44961132965460654, 0.44686884016237416, 0.4441221445704292,
    0.44137126873171667, 0.43861623853852766, 0.43585707992225547,
    0.43309381885315196, 0.43032648134008261, 0.42755509343028208,
    0.42477968120910881, 0.42200027079979968, 0.41921688836322391,
    0.41642956009763715, 0.4136383122384345, 0.41084317105790391,
    0.40804416286497869, 0.40524131400498986, 0.40243465085941843,
    0.39962419984564679, 0.39680998741671031, 0.3939920400610481,
    0.39117038430225387, 0.38834504669882625, 0.38551605384391885,
    0.38268343236508978, 0.37984720892405116, 0.37700741021641826,
    0.37416406297145793, 0.37131719395183749, 0.36846682995337232,
    0.36561299780477385, 0.36275572436739723, 0.35989503653498811,
    0.35703096123343, 0.35416352542049034, 0.35129275608556709,
    0.34841868024943456, 0.34554132496398909, 0.34266071731199438,
    0.33977688440682685, 0.33688985339222005, 0.33399965144200938,
    0.33110630575987643, 0.3282098435790925, 0.32531029216226293,
    0.32240767880106985, 0.31950203081601569, 0.31659337555616585,
    0.31368174039889152, 0.31076715274961147, 0.30784964004153487,
    0.30492922973540237, 0.30200594931922808, 0.29907982630804048,
    0.29615088824362379, 0.29321916269425863, 0.29028467725446233,
    0.28734745954472951, 0.28440753721127188, 0.28146493792575794,
    0.27851968938505306, 0.27557181931095814, 0.272621355449949,
    0.26966832557291509, 0.26671275747489837, 0.26375467897483135,
    0.26079411791527551, 0.257831102162159, 0.25486565960451457,
    0.25189781815421697, 0.24892760574572015, 0.24595505033579459,
    0.24298017990326387, 0.2400030224487415, 0.2370236059943672,
    0.23404195858354343, 0.23105810828067111, 0.22807208317088573,
    0.22508391135979283, 0.22209362097320351, 0.2191012401568698,
    0.21610679707621952, 0.21311031991609136, 0.21011183688046961,
    0.20711137619221856, 0.20410896609281687, 0.2011046348420919,
    0.19809841071795356, 0.19509032201612825, 0.19208039704989244,
    0.18906866414980619, 0.18605515166344663, 0.18303988795514095,
    0.18002290140569951, 0.17700422041214875, 0.17398387338746382,
    0.17096188876030122, 0.16793829497473117, 0.16491312048996992,
    0.16188639378011183, 0.15885814333386145, 0.15582839765426523,
    0.15279718525844344, 0.14976453467732151, 0.14673047445536175,
    0.14369503315029447, 0.14065823933284921, 0.13762012158648604,
    0.13458070850712617, 0.13154002870288312, 0.12849811079379317,
    0.12545498341154623, 0.1224106751992162, 0.11936521481099135,
    0.11631863091190475, 0.11327095217756435, 0.11022220729388306,
    0.10717242495680884, 0.10412163387205459, 0.10106986275482782,
    0.0980171403295606, 0.094963495329638992, 0.091908956497132724,
    0.0888535525825246, 0.0857973123444399, 0.082740264549375692,
    0.079682437971430126, 0.076623861392031492, 0.073564563599667426,
    0.070504573389613856, 0.067443919563664051, 0.064382630929857465,
    0.061320736302208578, 0.058258264500435752, 0.055195244349689941,
    0.052131704680283324, 0.049067674327418015, 0.046003182130914623,
    0.04293825693494082, 0.039872927587739811, 0.036807222941358832,
    0.03374117185137758, 0.030674803176636626, 0.02760814577896574,
    0.024541228522912288, 0.021474080275469508, 0.01840672990580482,
    0.0153392062849881, 0.012271538285719925, 0.00920375478205982,
    0.0061358846491544753, 0.0030679567629659761, 0.0, -0.0030679567629659761,
    -0.0061358846491544753, -0.00920375478205982, -0.012271538285719925,
    -0.0153392062849881, -0.01840672990580482, -0.021474080275469508,
    -0.024541228522912288, -0.02760814577896574, -0.030674803176636626,
    -0.03374117185137758, -0.036807222941358832, -0.039872927587739811,
    -0.04293825693494082, -0.046003182130914623, -0.049067674327418015,
    -0.052131704680283324, -0.055195244349689941, -0.058258264500435752,
    -0.061320736302208578, -0.064382630929857465, -0.067443919563664051,
    -0.070504573389613856, -0.073564563599667426, -0.076623861392031492,
    -0.079682437971430126, -0.082740264549375692, -0.0857973123444399,
    -0.0888535525825246, -0.091908956497132724, -0.094963495329638992,
    -0.0980171403295606, -0.10106986275482782, -0.10412163387205459,
    -0.10717242495680884, -0.11022220729388306, -0.11327095217756435,
    -0.11631863091190475, -0.11936521481099135, -0.1224106751992162,
    -0.12545498341154623, -0.12849811079379317, -0.13154002870288312,
    -0.13458070850712617, -0.13762012158648604, -0.14065823933284921,
    -0.14369503315029447, -0.14673047445536175, -0.14976453467732151,
    -0.15279718525844344, -0.15582839765426523, -0.15885814333386145,
    -0.16188639378011183, -0.16491312048996992, -0.16793829497473117,
    -0.17096188876030122, -0.17398387338746382, -0.17700422041214875,
    -0.18002290140569951, -0.18303988795514095, -0.18605515166344663,
    -0.18906866414980619, -0.19208039704989244, -0.19509032201612825,
    -0.19809841071795356, -0.2011046348420919, -0.20410896609281687,
    -0.20711137619221856, -0.21011183688046961, -0.21311031991609136,
    -0.21610679707621952, -0.2191012401568698, -0.22209362097320351,
    -0.22508391135979283, -0.22807208317088573, -0.23105810828067111,
    -0.23404195858354343, -0.2370236059943672, -0.2400030224487415,
    -0.24298017990326387, -0.24595505033579459, -0.24892760574572015,
    -0.25189781815421697, -0.25486565960451457, -0.257831102162159,
    -0.26079411791527551, -0.26375467897483135, -0.26671275747489837,
    -0.26966832557291509, -0.272621355449949, -0.27557181931095814,
    -0.27851968938505306, -0.28146493792575794, -0.28440753721127188,
    -0.28734745954472951, -0.29028467725446233, -0.29321916269425863,
    -0.29615088824362379, -0.29907982630804048, -0.30200594931922808,
    -0.30492922973540237, -0.30784964004153487, -0.31076715274961147,
    -0.31368174039889152, -0.31659337555616585, -0.31950203081601569,
    -0.32240767880106985, -0.32531029216226293, -0.3282098435790925,
    -0.33110630575987643, -0.33399965144200938, -0.33688985339222005,
    -0.33977688440682685, -0.34266071731199438, -0.34554132496398909,
    -0.34841868024943456, -0.35129275608556709, -0.35416352542049034,
    -0.35703096123343, -0.35989503653498811, -0.36275572436739723,
    -0.36561299780477385, -0.36846682995337232, -0.37131719395183749,
    -0.37416406297145793, -0.37700741021641826, -0.37984720892405116,
    -0.38268343236508978, -0.38551605384391885, -0.38834504669882625,
    -0.39117038430225387, -0.3939920400610481, -0.39680998741671031,
    -0.39962419984564679, -0.40243465085941843, -0.40524131400498986,
    -0.40804416286497869, -0.41084317105790391, -0.4136383122384345,
    -0.41642956009763715, -0.41921688836322391, -0.42200027079979968,
    -0.42477968120910881, -0.42755509343028208, -0.43032648134008261,
    -0.43309381885315196, -0.43585707992225547, -0.43861623853852766,
    -0.44137126873171667, -0.4441221445704292, -0.44686884016237416,
    -0.44961132965460654, -0.45234958723377089, -0.45508358712634384,
    -0.45781330359887717, -0.46053871095824, -0.46325978355186015,
    -0.46597649576796618, -0.46868882203582796, -0.47139673682599764,
    -0.47410021465054997, -0.47679923006332209, -0.47949375766015295,
    -0.48218377207912272, -0.48486924800079106, -0.487550160148436,
    -0.49022648328829116, -0.49289819222978404, -0.49556526182577254,
    -0.49822766697278181, -0.50088538261124071, -0.50353838372571758,
    -0.50618664534515523, -0.508830142543107, -0.5114688504379703,
    -0.51410274419322166, -0.51673179901764987, -0.51935599016558964,
    -0.52197529293715439, -0.524589682678469, -0.52719913478190128,
    -0.52980362468629461, -0.5324031278771979, -0.53499761988709715,
    -0.53758707629564539, -0.54017147272989285, -0.54275078486451589,
    -0.54532498842204646, -0.54789405917310019, -0.55045797293660481,
    -0.55301670558002747, -0.55557023301960218, -0.5581185312205561,
    -0.560661576197336, -0.56319934401383409, -0.56573181078361312,
    -0.56825895267013149, -0.57078074588696726, -0.5732971666980422,
    -0.57580819141784534, -0.57831379641165559, -0.58081395809576453,
    -0.58330865293769829, -0.58579785745643886, -0.58828154822264522,
    -0.59075970185887416, -0.5932322950397998, -0.59569930449243336,
    -0.59816070699634238, -0.600616479383869, -0.60306659854034816,
    -0.60551104140432555, -0.60794978496777363, -0.61038280627630948,
    -0.61281008242940971, -0.61523159058062682, -0.61764730793780387,
    -0.6200572117632891, -0.62246127937415, -0.62485948814238634,
    -0.62725181549514408, -0.629638238914927, -0.63201873593980906,
    -0.63439328416364549, -0.6367618612362842, -0.63912444486377573,
    -0.641481012808583, -0.64383154288979139, -0.64617601298331628,
    -0.64851440102211244, -0.650846684996381, -0.65317284295377676,
    -0.65549285299961535, -0.65780669329707864, -0.66011434206742048,
    -0.66241577759017178, -0.66471097820334479, -0.66699992230363747,
    -0.669282588346636, -0.67155895484701833, -0.673829000378756,
    -0.67609270357531592, -0.67835004312986147, -0.680600997795453,
    -0.68284554638524808, -0.68508366777270036, -0.687315340891759,
    -0.68954054473706683, -0.69175925836415775, -0.69397146088965389,
    -0.696177131491463, -0.69837624940897292, -0.70056879394324834,
    -0.7027547444572253, -0.70493408037590488, -0.70710678118654757,
    -0.70927282643886569, -0.71143219574521643, -0.71358486878079352,
    -0.71573082528381859, -0.71787004505573171, -0.72000250796138165,
    -0.72212819392921535, -0.724247082951467, -0.726359155084346,
    -0.7284643904482252, -0.73056276922782759, -0.73265427167241282,
    -0.7347388780959635, -0.73681656887736979, -0.73888732446061511,
    -0.74095112535495922, -0.74300795213512172, -0.745057785441466,
    -0.74710060598018013, -0.74913639452345937, -0.75116513190968637,
    -0.75318679904361252, -0.75520137689653655, -0.75720884650648457,
    -0.759209188978388, -0.76120238548426178, -0.76318841726338127,
    -0.765167265622459, -0.7671389119358204, -0.7691033376455797,
    -0.77106052426181382, -0.773010453362737, -0.77495310659487393,
    -0.77688846567323244, -0.778816512381476, -0.78073722857209449,
    -0.78265059616657573, -0.78455659715557524, -0.78645521359908577,
    -0.78834642762660634, -0.79023022143731, -0.79210657730021239,
    -0.79397547755433717, -0.79583690460888357, -0.79769084094339116,
    -0.799537269107905, -0.80137617172314024, -0.80320753148064494,
    -0.80503133114296366, -0.80684755354379933, -0.808656181588175,
    -0.81045719825259477, -0.81225058658520388, -0.81403632970594841,
    -0.81581441080673378, -0.81758481315158371, -0.819347520076797,
    -0.82110251499110465, -0.82284978137582643, -0.82458930278502529,
    -0.82632106284566353, -0.8280450452577558, -0.829761233794523,
    -0.83146961230254524, -0.83317016470191319, -0.83486287498638,
    -0.836547727223512, -0.83822470555483808, -0.83989379419599952,
    -0.84155497743689844, -0.84320823964184544, -0.84485356524970712,
    -0.84649093877405213, -0.84812034480329723, -0.84974176800085255,
    -0.8513551931052652, -0.85296060493036363, -0.85455798836540053,
    -0.85614732837519447, -0.85772861000027212, -0.85930181835700847,
    -0.86086693863776731, -0.8624239561110405, -0.8639728561215867,
    -0.86551362409056909, -0.86704624551569265, -0.8685707059713409,
    -0.87008699110871146, -0.87159508665595109, -0.87309497841829009,
    -0.87458665227817611, -0.8760700941954066, -0.87754529020726135,
    -0.87901222642863353, -0.88047088905216075, -0.881921264348355,
    -0.88336333866573158, -0.88479709843093779, -0.88622253014888064,
    -0.88763962040285393, -0.88904835585466457, -0.89044872324475788,
    -0.89184070939234272, -0.89322430119551532, -0.8945994856313827,
    -0.89596624975618522, -0.89732458070541832, -0.89867446569395382,
    -0.90001589201616017, -0.901348847046022, -0.90267331823725883,
    -0.90398929312344334, -0.90529675931811882, -0.90659570451491533,
    -0.90788611648766626, -0.90916798309052238, -0.91044129225806725,
    -0.91170603200542988, -0.91296219042839821, -0.91420975570353069,
    -0.91544871608826783, -0.9166790599210427, -0.9179007756213905,
    -0.91911385169005777, -0.92031827670911059, -0.9215140393420419,
    -0.92270112833387863, -0.92387953251128674, -0.92504924078267758,
    -0.92621024213831138, -0.92736252565040111, -0.92850608047321559,
    -0.92964089584318121, -0.93076696107898371, -0.93188426558166815,
    -0.932992798834739, -0.93409255040425887, -0.93518350993894761,
    -0.93626566717027826, -0.937339011912575, -0.93840353406310806,
    -0.93945922360218992, -0.9405060705932683, -0.94154406518302081,
    -0.94257319760144687, -0.94359345816196039, -0.94460483726148026,
    -0.94560732538052128, -0.94660091308328353, -0.94758559101774109,
    -0.94856134991573027, -0.94952818059303667, -0.9504860739494817,
    -0.95143502096900834, -0.95237501271976588, -0.95330604035419386,
    -0.95422809510910567, -0.95514116830577078, -0.95604525134999641,
    -0.95694033573220882, -0.95782641302753291, -0.9587034748958716,
    -0.95957151308198452, -0.96043051941556579, -0.96128048581132064,
    -0.96212140426904158, -0.96295326687368388, -0.96377606579543984,
    -0.96458979328981276, -0.9653944416976894, -0.9661900034454125,
    -0.96697647104485207, -0.96775383709347551, -0.96852209427441727,
    -0.96928123535654853, -0.970031253194544, -0.97077214072895035,
    -0.97150389098625178, -0.97222649707893627, -0.97293995220556018,
    -0.973644249650812, -0.97433938278557586, -0.97502534506699412,
    -0.97570213003852857, -0.97636973133002114, -0.97702814265775439,
    -0.97767735782450993, -0.97831737071962765, -0.9789481753190622,
    -0.97956976568544052, -0.98018213596811743, -0.98078528040323043,
    -0.98137919331375456, -0.98196386910955524, -0.98253930228744124,
    -0.98310548743121629, -0.98366241921173025, -0.984210092386929,
    -0.98474850180190421, -0.98527764238894122, -0.98579750916756748,
    -0.98630809724459867, -0.98680940181418553, -0.98730141815785843,
    -0.98778414164457218, -0.98825756773074946, -0.98872169196032378,
    -0.989176509964781, -0.98962201746320089, -0.99005821026229712,
    -0.99048508425645709, -0.99090263542778, -0.99131085984611544,
    -0.99170975366909953, -0.9920993131421918, -0.99247953459871,
    -0.9928504144598651, -0.9932119492347945, -0.9935641355205953,
    -0.99390697000235606, -0.9942404494531879, -0.99456457073425542,
    -0.99487933079480562, -0.99518472667219693, -0.99548075549192694,
    -0.99576741446765982, -0.996044700901252, -0.996312612182778,
    -0.99657114579055484, -0.99682029929116567, -0.997060070339483,
    -0.99729045667869021, -0.99751145614030345, -0.99772306664419164,
    -0.997925286198596, -0.99811811290014918, -0.99830154493389289,
    -0.99847558057329477, -0.99864021818026527, -0.99879545620517241,
    -0.99894129318685687, -0.99907772775264536, -0.99920475861836389,
    -0.99932238458834954, -0.99943060455546173, -0.99952941750109314,
    -0.99961882249517864, -0.99969881869620425, -0.99976940535121528,
    -0.9998305817958234, -0.99988234745421256, -0.9999247018391445,
    -0.9999576445519639, -0.99998117528260111, -0.99999529380957619, -1.0 };

  double twid_im;
  static const double dv18[1025] = { 0.0, -0.0030679567629659761,
    -0.0061358846491544753, -0.00920375478205982, -0.012271538285719925,
    -0.0153392062849881, -0.01840672990580482, -0.021474080275469508,
    -0.024541228522912288, -0.02760814577896574, -0.030674803176636626,
    -0.03374117185137758, -0.036807222941358832, -0.039872927587739811,
    -0.04293825693494082, -0.046003182130914623, -0.049067674327418015,
    -0.052131704680283324, -0.055195244349689941, -0.058258264500435752,
    -0.061320736302208578, -0.064382630929857465, -0.067443919563664051,
    -0.070504573389613856, -0.073564563599667426, -0.076623861392031492,
    -0.079682437971430126, -0.082740264549375692, -0.0857973123444399,
    -0.0888535525825246, -0.091908956497132724, -0.094963495329638992,
    -0.0980171403295606, -0.10106986275482782, -0.10412163387205459,
    -0.10717242495680884, -0.11022220729388306, -0.11327095217756435,
    -0.11631863091190475, -0.11936521481099135, -0.1224106751992162,
    -0.12545498341154623, -0.12849811079379317, -0.13154002870288312,
    -0.13458070850712617, -0.13762012158648604, -0.14065823933284921,
    -0.14369503315029447, -0.14673047445536175, -0.14976453467732151,
    -0.15279718525844344, -0.15582839765426523, -0.15885814333386145,
    -0.16188639378011183, -0.16491312048996992, -0.16793829497473117,
    -0.17096188876030122, -0.17398387338746382, -0.17700422041214875,
    -0.18002290140569951, -0.18303988795514095, -0.18605515166344663,
    -0.18906866414980619, -0.19208039704989244, -0.19509032201612825,
    -0.19809841071795356, -0.2011046348420919, -0.20410896609281687,
    -0.20711137619221856, -0.21011183688046961, -0.21311031991609136,
    -0.21610679707621952, -0.2191012401568698, -0.22209362097320351,
    -0.22508391135979283, -0.22807208317088573, -0.23105810828067111,
    -0.23404195858354343, -0.2370236059943672, -0.2400030224487415,
    -0.24298017990326387, -0.24595505033579459, -0.24892760574572015,
    -0.25189781815421697, -0.25486565960451457, -0.257831102162159,
    -0.26079411791527551, -0.26375467897483135, -0.26671275747489837,
    -0.26966832557291509, -0.272621355449949, -0.27557181931095814,
    -0.27851968938505306, -0.28146493792575794, -0.28440753721127188,
    -0.28734745954472951, -0.29028467725446233, -0.29321916269425863,
    -0.29615088824362379, -0.29907982630804048, -0.30200594931922808,
    -0.30492922973540237, -0.30784964004153487, -0.31076715274961147,
    -0.31368174039889152, -0.31659337555616585, -0.31950203081601569,
    -0.32240767880106985, -0.32531029216226293, -0.3282098435790925,
    -0.33110630575987643, -0.33399965144200938, -0.33688985339222005,
    -0.33977688440682685, -0.34266071731199438, -0.34554132496398909,
    -0.34841868024943456, -0.35129275608556709, -0.35416352542049034,
    -0.35703096123343, -0.35989503653498811, -0.36275572436739723,
    -0.36561299780477385, -0.36846682995337232, -0.37131719395183749,
    -0.37416406297145793, -0.37700741021641826, -0.37984720892405116,
    -0.38268343236508978, -0.38551605384391885, -0.38834504669882625,
    -0.39117038430225387, -0.3939920400610481, -0.39680998741671031,
    -0.39962419984564679, -0.40243465085941843, -0.40524131400498986,
    -0.40804416286497869, -0.41084317105790391, -0.4136383122384345,
    -0.41642956009763715, -0.41921688836322391, -0.42200027079979968,
    -0.42477968120910881, -0.42755509343028208, -0.43032648134008261,
    -0.43309381885315196, -0.43585707992225547, -0.43861623853852766,
    -0.44137126873171667, -0.4441221445704292, -0.44686884016237416,
    -0.44961132965460654, -0.45234958723377089, -0.45508358712634384,
    -0.45781330359887717, -0.46053871095824, -0.46325978355186015,
    -0.46597649576796618, -0.46868882203582796, -0.47139673682599764,
    -0.47410021465054997, -0.47679923006332209, -0.47949375766015295,
    -0.48218377207912272, -0.48486924800079106, -0.487550160148436,
    -0.49022648328829116, -0.49289819222978404, -0.49556526182577254,
    -0.49822766697278181, -0.50088538261124071, -0.50353838372571758,
    -0.50618664534515523, -0.508830142543107, -0.5114688504379703,
    -0.51410274419322166, -0.51673179901764987, -0.51935599016558964,
    -0.52197529293715439, -0.524589682678469, -0.52719913478190128,
    -0.52980362468629461, -0.5324031278771979, -0.53499761988709715,
    -0.53758707629564539, -0.54017147272989285, -0.54275078486451589,
    -0.54532498842204646, -0.54789405917310019, -0.55045797293660481,
    -0.55301670558002747, -0.55557023301960218, -0.5581185312205561,
    -0.560661576197336, -0.56319934401383409, -0.56573181078361312,
    -0.56825895267013149, -0.57078074588696726, -0.5732971666980422,
    -0.57580819141784534, -0.57831379641165559, -0.58081395809576453,
    -0.58330865293769829, -0.58579785745643886, -0.58828154822264522,
    -0.59075970185887416, -0.5932322950397998, -0.59569930449243336,
    -0.59816070699634238, -0.600616479383869, -0.60306659854034816,
    -0.60551104140432555, -0.60794978496777363, -0.61038280627630948,
    -0.61281008242940971, -0.61523159058062682, -0.61764730793780387,
    -0.6200572117632891, -0.62246127937415, -0.62485948814238634,
    -0.62725181549514408, -0.629638238914927, -0.63201873593980906,
    -0.63439328416364549, -0.6367618612362842, -0.63912444486377573,
    -0.641481012808583, -0.64383154288979139, -0.64617601298331628,
    -0.64851440102211244, -0.650846684996381, -0.65317284295377676,
    -0.65549285299961535, -0.65780669329707864, -0.66011434206742048,
    -0.66241577759017178, -0.66471097820334479, -0.66699992230363747,
    -0.669282588346636, -0.67155895484701833, -0.673829000378756,
    -0.67609270357531592, -0.67835004312986147, -0.680600997795453,
    -0.68284554638524808, -0.68508366777270036, -0.687315340891759,
    -0.68954054473706683, -0.69175925836415775, -0.69397146088965389,
    -0.696177131491463, -0.69837624940897292, -0.70056879394324834,
    -0.7027547444572253, -0.70493408037590488, -0.70710678118654757,
    -0.70927282643886569, -0.71143219574521643, -0.71358486878079352,
    -0.71573082528381859, -0.71787004505573171, -0.72000250796138165,
    -0.72212819392921535, -0.724247082951467, -0.726359155084346,
    -0.7284643904482252, -0.73056276922782759, -0.73265427167241282,
    -0.7347388780959635, -0.73681656887736979, -0.73888732446061511,
    -0.74095112535495922, -0.74300795213512172, -0.745057785441466,
    -0.74710060598018013, -0.74913639452345937, -0.75116513190968637,
    -0.75318679904361252, -0.75520137689653655, -0.75720884650648457,
    -0.759209188978388, -0.76120238548426178, -0.76318841726338127,
    -0.765167265622459, -0.7671389119358204, -0.7691033376455797,
    -0.77106052426181382, -0.773010453362737, -0.77495310659487393,
    -0.77688846567323244, -0.778816512381476, -0.78073722857209449,
    -0.78265059616657573, -0.78455659715557524, -0.78645521359908577,
    -0.78834642762660634, -0.79023022143731, -0.79210657730021239,
    -0.79397547755433717, -0.79583690460888357, -0.79769084094339116,
    -0.799537269107905, -0.80137617172314024, -0.80320753148064494,
    -0.80503133114296366, -0.80684755354379933, -0.808656181588175,
    -0.81045719825259477, -0.81225058658520388, -0.81403632970594841,
    -0.81581441080673378, -0.81758481315158371, -0.819347520076797,
    -0.82110251499110465, -0.82284978137582643, -0.82458930278502529,
    -0.82632106284566353, -0.8280450452577558, -0.829761233794523,
    -0.83146961230254524, -0.83317016470191319, -0.83486287498638,
    -0.836547727223512, -0.83822470555483808, -0.83989379419599952,
    -0.84155497743689844, -0.84320823964184544, -0.84485356524970712,
    -0.84649093877405213, -0.84812034480329723, -0.84974176800085255,
    -0.8513551931052652, -0.85296060493036363, -0.85455798836540053,
    -0.85614732837519447, -0.85772861000027212, -0.85930181835700847,
    -0.86086693863776731, -0.8624239561110405, -0.8639728561215867,
    -0.86551362409056909, -0.86704624551569265, -0.8685707059713409,
    -0.87008699110871146, -0.87159508665595109, -0.87309497841829009,
    -0.87458665227817611, -0.8760700941954066, -0.87754529020726135,
    -0.87901222642863353, -0.88047088905216075, -0.881921264348355,
    -0.88336333866573158, -0.88479709843093779, -0.88622253014888064,
    -0.88763962040285393, -0.88904835585466457, -0.89044872324475788,
    -0.89184070939234272, -0.89322430119551532, -0.8945994856313827,
    -0.89596624975618522, -0.89732458070541832, -0.89867446569395382,
    -0.90001589201616017, -0.901348847046022, -0.90267331823725883,
    -0.90398929312344334, -0.90529675931811882, -0.90659570451491533,
    -0.90788611648766626, -0.90916798309052238, -0.91044129225806725,
    -0.91170603200542988, -0.91296219042839821, -0.91420975570353069,
    -0.91544871608826783, -0.9166790599210427, -0.9179007756213905,
    -0.91911385169005777, -0.92031827670911059, -0.9215140393420419,
    -0.92270112833387863, -0.92387953251128674, -0.92504924078267758,
    -0.92621024213831138, -0.92736252565040111, -0.92850608047321559,
    -0.92964089584318121, -0.93076696107898371, -0.93188426558166815,
    -0.932992798834739, -0.93409255040425887, -0.93518350993894761,
    -0.93626566717027826, -0.937339011912575, -0.93840353406310806,
    -0.93945922360218992, -0.9405060705932683, -0.94154406518302081,
    -0.94257319760144687, -0.94359345816196039, -0.94460483726148026,
    -0.94560732538052128, -0.94660091308328353, -0.94758559101774109,
    -0.94856134991573027, -0.94952818059303667, -0.9504860739494817,
    -0.95143502096900834, -0.95237501271976588, -0.95330604035419386,
    -0.95422809510910567, -0.95514116830577078, -0.95604525134999641,
    -0.95694033573220882, -0.95782641302753291, -0.9587034748958716,
    -0.95957151308198452, -0.96043051941556579, -0.96128048581132064,
    -0.96212140426904158, -0.96295326687368388, -0.96377606579543984,
    -0.96458979328981276, -0.9653944416976894, -0.9661900034454125,
    -0.96697647104485207, -0.96775383709347551, -0.96852209427441727,
    -0.96928123535654853, -0.970031253194544, -0.97077214072895035,
    -0.97150389098625178, -0.97222649707893627, -0.97293995220556018,
    -0.973644249650812, -0.97433938278557586, -0.97502534506699412,
    -0.97570213003852857, -0.97636973133002114, -0.97702814265775439,
    -0.97767735782450993, -0.97831737071962765, -0.9789481753190622,
    -0.97956976568544052, -0.98018213596811743, -0.98078528040323043,
    -0.98137919331375456, -0.98196386910955524, -0.98253930228744124,
    -0.98310548743121629, -0.98366241921173025, -0.984210092386929,
    -0.98474850180190421, -0.98527764238894122, -0.98579750916756748,
    -0.98630809724459867, -0.98680940181418553, -0.98730141815785843,
    -0.98778414164457218, -0.98825756773074946, -0.98872169196032378,
    -0.989176509964781, -0.98962201746320089, -0.99005821026229712,
    -0.99048508425645709, -0.99090263542778, -0.99131085984611544,
    -0.99170975366909953, -0.9920993131421918, -0.99247953459871,
    -0.9928504144598651, -0.9932119492347945, -0.9935641355205953,
    -0.99390697000235606, -0.9942404494531879, -0.99456457073425542,
    -0.99487933079480562, -0.99518472667219693, -0.99548075549192694,
    -0.99576741446765982, -0.996044700901252, -0.996312612182778,
    -0.99657114579055484, -0.99682029929116567, -0.997060070339483,
    -0.99729045667869021, -0.99751145614030345, -0.99772306664419164,
    -0.997925286198596, -0.99811811290014918, -0.99830154493389289,
    -0.99847558057329477, -0.99864021818026527, -0.99879545620517241,
    -0.99894129318685687, -0.99907772775264536, -0.99920475861836389,
    -0.99932238458834954, -0.99943060455546173, -0.99952941750109314,
    -0.99961882249517864, -0.99969881869620425, -0.99976940535121528,
    -0.9998305817958234, -0.99988234745421256, -0.9999247018391445,
    -0.9999576445519639, -0.99998117528260111, -0.99999529380957619, -1.0,
    -0.99999529380957619, -0.99998117528260111, -0.9999576445519639,
    -0.9999247018391445, -0.99988234745421256, -0.9998305817958234,
    -0.99976940535121528, -0.99969881869620425, -0.99961882249517864,
    -0.99952941750109314, -0.99943060455546173, -0.99932238458834954,
    -0.99920475861836389, -0.99907772775264536, -0.99894129318685687,
    -0.99879545620517241, -0.99864021818026527, -0.99847558057329477,
    -0.99830154493389289, -0.99811811290014918, -0.997925286198596,
    -0.99772306664419164, -0.99751145614030345, -0.99729045667869021,
    -0.997060070339483, -0.99682029929116567, -0.99657114579055484,
    -0.996312612182778, -0.996044700901252, -0.99576741446765982,
    -0.99548075549192694, -0.99518472667219693, -0.99487933079480562,
    -0.99456457073425542, -0.9942404494531879, -0.99390697000235606,
    -0.9935641355205953, -0.9932119492347945, -0.9928504144598651,
    -0.99247953459871, -0.9920993131421918, -0.99170975366909953,
    -0.99131085984611544, -0.99090263542778, -0.99048508425645709,
    -0.99005821026229712, -0.98962201746320089, -0.989176509964781,
    -0.98872169196032378, -0.98825756773074946, -0.98778414164457218,
    -0.98730141815785843, -0.98680940181418553, -0.98630809724459867,
    -0.98579750916756748, -0.98527764238894122, -0.98474850180190421,
    -0.984210092386929, -0.98366241921173025, -0.98310548743121629,
    -0.98253930228744124, -0.98196386910955524, -0.98137919331375456,
    -0.98078528040323043, -0.98018213596811743, -0.97956976568544052,
    -0.9789481753190622, -0.97831737071962765, -0.97767735782450993,
    -0.97702814265775439, -0.97636973133002114, -0.97570213003852857,
    -0.97502534506699412, -0.97433938278557586, -0.973644249650812,
    -0.97293995220556018, -0.97222649707893627, -0.97150389098625178,
    -0.97077214072895035, -0.970031253194544, -0.96928123535654853,
    -0.96852209427441727, -0.96775383709347551, -0.96697647104485207,
    -0.9661900034454125, -0.9653944416976894, -0.96458979328981276,
    -0.96377606579543984, -0.96295326687368388, -0.96212140426904158,
    -0.96128048581132064, -0.96043051941556579, -0.95957151308198452,
    -0.9587034748958716, -0.95782641302753291, -0.95694033573220882,
    -0.95604525134999641, -0.95514116830577078, -0.95422809510910567,
    -0.95330604035419386, -0.95237501271976588, -0.95143502096900834,
    -0.9504860739494817, -0.94952818059303667, -0.94856134991573027,
    -0.94758559101774109, -0.94660091308328353, -0.94560732538052128,
    -0.94460483726148026, -0.94359345816196039, -0.94257319760144687,
    -0.94154406518302081, -0.9405060705932683, -0.93945922360218992,
    -0.93840353406310806, -0.937339011912575, -0.93626566717027826,
    -0.93518350993894761, -0.93409255040425887, -0.932992798834739,
    -0.93188426558166815, -0.93076696107898371, -0.92964089584318121,
    -0.92850608047321559, -0.92736252565040111, -0.92621024213831138,
    -0.92504924078267758, -0.92387953251128674, -0.92270112833387863,
    -0.9215140393420419, -0.92031827670911059, -0.91911385169005777,
    -0.9179007756213905, -0.9166790599210427, -0.91544871608826783,
    -0.91420975570353069, -0.91296219042839821, -0.91170603200542988,
    -0.91044129225806725, -0.90916798309052238, -0.90788611648766626,
    -0.90659570451491533, -0.90529675931811882, -0.90398929312344334,
    -0.90267331823725883, -0.901348847046022, -0.90001589201616017,
    -0.89867446569395382, -0.89732458070541832, -0.89596624975618522,
    -0.8945994856313827, -0.89322430119551532, -0.89184070939234272,
    -0.89044872324475788, -0.88904835585466457, -0.88763962040285393,
    -0.88622253014888064, -0.88479709843093779, -0.88336333866573158,
    -0.881921264348355, -0.88047088905216075, -0.87901222642863353,
    -0.87754529020726135, -0.8760700941954066, -0.87458665227817611,
    -0.87309497841829009, -0.87159508665595109, -0.87008699110871146,
    -0.8685707059713409, -0.86704624551569265, -0.86551362409056909,
    -0.8639728561215867, -0.8624239561110405, -0.86086693863776731,
    -0.85930181835700847, -0.85772861000027212, -0.85614732837519447,
    -0.85455798836540053, -0.85296060493036363, -0.8513551931052652,
    -0.84974176800085255, -0.84812034480329723, -0.84649093877405213,
    -0.84485356524970712, -0.84320823964184544, -0.84155497743689844,
    -0.83989379419599952, -0.83822470555483808, -0.836547727223512,
    -0.83486287498638, -0.83317016470191319, -0.83146961230254524,
    -0.829761233794523, -0.8280450452577558, -0.82632106284566353,
    -0.82458930278502529, -0.82284978137582643, -0.82110251499110465,
    -0.819347520076797, -0.81758481315158371, -0.81581441080673378,
    -0.81403632970594841, -0.81225058658520388, -0.81045719825259477,
    -0.808656181588175, -0.80684755354379933, -0.80503133114296366,
    -0.80320753148064494, -0.80137617172314024, -0.799537269107905,
    -0.79769084094339116, -0.79583690460888357, -0.79397547755433717,
    -0.79210657730021239, -0.79023022143731, -0.78834642762660634,
    -0.78645521359908577, -0.78455659715557524, -0.78265059616657573,
    -0.78073722857209449, -0.778816512381476, -0.77688846567323244,
    -0.77495310659487393, -0.773010453362737, -0.77106052426181382,
    -0.7691033376455797, -0.7671389119358204, -0.765167265622459,
    -0.76318841726338127, -0.76120238548426178, -0.759209188978388,
    -0.75720884650648457, -0.75520137689653655, -0.75318679904361252,
    -0.75116513190968637, -0.74913639452345937, -0.74710060598018013,
    -0.745057785441466, -0.74300795213512172, -0.74095112535495922,
    -0.73888732446061511, -0.73681656887736979, -0.7347388780959635,
    -0.73265427167241282, -0.73056276922782759, -0.7284643904482252,
    -0.726359155084346, -0.724247082951467, -0.72212819392921535,
    -0.72000250796138165, -0.71787004505573171, -0.71573082528381859,
    -0.71358486878079352, -0.71143219574521643, -0.70927282643886569,
    -0.70710678118654757, -0.70493408037590488, -0.7027547444572253,
    -0.70056879394324834, -0.69837624940897292, -0.696177131491463,
    -0.69397146088965389, -0.69175925836415775, -0.68954054473706683,
    -0.687315340891759, -0.68508366777270036, -0.68284554638524808,
    -0.680600997795453, -0.67835004312986147, -0.67609270357531592,
    -0.673829000378756, -0.67155895484701833, -0.669282588346636,
    -0.66699992230363747, -0.66471097820334479, -0.66241577759017178,
    -0.66011434206742048, -0.65780669329707864, -0.65549285299961535,
    -0.65317284295377676, -0.650846684996381, -0.64851440102211244,
    -0.64617601298331628, -0.64383154288979139, -0.641481012808583,
    -0.63912444486377573, -0.6367618612362842, -0.63439328416364549,
    -0.63201873593980906, -0.629638238914927, -0.62725181549514408,
    -0.62485948814238634, -0.62246127937415, -0.6200572117632891,
    -0.61764730793780387, -0.61523159058062682, -0.61281008242940971,
    -0.61038280627630948, -0.60794978496777363, -0.60551104140432555,
    -0.60306659854034816, -0.600616479383869, -0.59816070699634238,
    -0.59569930449243336, -0.5932322950397998, -0.59075970185887416,
    -0.58828154822264522, -0.58579785745643886, -0.58330865293769829,
    -0.58081395809576453, -0.57831379641165559, -0.57580819141784534,
    -0.5732971666980422, -0.57078074588696726, -0.56825895267013149,
    -0.56573181078361312, -0.56319934401383409, -0.560661576197336,
    -0.5581185312205561, -0.55557023301960218, -0.55301670558002747,
    -0.55045797293660481, -0.54789405917310019, -0.54532498842204646,
    -0.54275078486451589, -0.54017147272989285, -0.53758707629564539,
    -0.53499761988709715, -0.5324031278771979, -0.52980362468629461,
    -0.52719913478190128, -0.524589682678469, -0.52197529293715439,
    -0.51935599016558964, -0.51673179901764987, -0.51410274419322166,
    -0.5114688504379703, -0.508830142543107, -0.50618664534515523,
    -0.50353838372571758, -0.50088538261124071, -0.49822766697278181,
    -0.49556526182577254, -0.49289819222978404, -0.49022648328829116,
    -0.487550160148436, -0.48486924800079106, -0.48218377207912272,
    -0.47949375766015295, -0.47679923006332209, -0.47410021465054997,
    -0.47139673682599764, -0.46868882203582796, -0.46597649576796618,
    -0.46325978355186015, -0.46053871095824, -0.45781330359887717,
    -0.45508358712634384, -0.45234958723377089, -0.44961132965460654,
    -0.44686884016237416, -0.4441221445704292, -0.44137126873171667,
    -0.43861623853852766, -0.43585707992225547, -0.43309381885315196,
    -0.43032648134008261, -0.42755509343028208, -0.42477968120910881,
    -0.42200027079979968, -0.41921688836322391, -0.41642956009763715,
    -0.4136383122384345, -0.41084317105790391, -0.40804416286497869,
    -0.40524131400498986, -0.40243465085941843, -0.39962419984564679,
    -0.39680998741671031, -0.3939920400610481, -0.39117038430225387,
    -0.38834504669882625, -0.38551605384391885, -0.38268343236508978,
    -0.37984720892405116, -0.37700741021641826, -0.37416406297145793,
    -0.37131719395183749, -0.36846682995337232, -0.36561299780477385,
    -0.36275572436739723, -0.35989503653498811, -0.35703096123343,
    -0.35416352542049034, -0.35129275608556709, -0.34841868024943456,
    -0.34554132496398909, -0.34266071731199438, -0.33977688440682685,
    -0.33688985339222005, -0.33399965144200938, -0.33110630575987643,
    -0.3282098435790925, -0.32531029216226293, -0.32240767880106985,
    -0.31950203081601569, -0.31659337555616585, -0.31368174039889152,
    -0.31076715274961147, -0.30784964004153487, -0.30492922973540237,
    -0.30200594931922808, -0.29907982630804048, -0.29615088824362379,
    -0.29321916269425863, -0.29028467725446233, -0.28734745954472951,
    -0.28440753721127188, -0.28146493792575794, -0.27851968938505306,
    -0.27557181931095814, -0.272621355449949, -0.26966832557291509,
    -0.26671275747489837, -0.26375467897483135, -0.26079411791527551,
    -0.257831102162159, -0.25486565960451457, -0.25189781815421697,
    -0.24892760574572015, -0.24595505033579459, -0.24298017990326387,
    -0.2400030224487415, -0.2370236059943672, -0.23404195858354343,
    -0.23105810828067111, -0.22807208317088573, -0.22508391135979283,
    -0.22209362097320351, -0.2191012401568698, -0.21610679707621952,
    -0.21311031991609136, -0.21011183688046961, -0.20711137619221856,
    -0.20410896609281687, -0.2011046348420919, -0.19809841071795356,
    -0.19509032201612825, -0.19208039704989244, -0.18906866414980619,
    -0.18605515166344663, -0.18303988795514095, -0.18002290140569951,
    -0.17700422041214875, -0.17398387338746382, -0.17096188876030122,
    -0.16793829497473117, -0.16491312048996992, -0.16188639378011183,
    -0.15885814333386145, -0.15582839765426523, -0.15279718525844344,
    -0.14976453467732151, -0.14673047445536175, -0.14369503315029447,
    -0.14065823933284921, -0.13762012158648604, -0.13458070850712617,
    -0.13154002870288312, -0.12849811079379317, -0.12545498341154623,
    -0.1224106751992162, -0.11936521481099135, -0.11631863091190475,
    -0.11327095217756435, -0.11022220729388306, -0.10717242495680884,
    -0.10412163387205459, -0.10106986275482782, -0.0980171403295606,
    -0.094963495329638992, -0.091908956497132724, -0.0888535525825246,
    -0.0857973123444399, -0.082740264549375692, -0.079682437971430126,
    -0.076623861392031492, -0.073564563599667426, -0.070504573389613856,
    -0.067443919563664051, -0.064382630929857465, -0.061320736302208578,
    -0.058258264500435752, -0.055195244349689941, -0.052131704680283324,
    -0.049067674327418015, -0.046003182130914623, -0.04293825693494082,
    -0.039872927587739811, -0.036807222941358832, -0.03374117185137758,
    -0.030674803176636626, -0.02760814577896574, -0.024541228522912288,
    -0.021474080275469508, -0.01840672990580482, -0.0153392062849881,
    -0.012271538285719925, -0.00920375478205982, -0.0061358846491544753,
    -0.0030679567629659761, -0.0 };

  int ihi;
  iy = X->size[1];
  if (iy == 0) {
    iy = X->size[1];
    if (2048 > iy) {
      for (i = 0; i < 2048; i++) {
        b_y1[i].re = 0.0;
        b_y1[i].im = 0.0;
      }
    }
  } else {
    emxInit_real_T1(&x, 1);
    iy = X->size[1];
    iDelta2 = x->size[0];
    x->size[0] = iy;
    emxEnsureCapacity((emxArray__common *)x, iDelta2, (int)sizeof(double));
    for (iDelta2 = 0; iDelta2 < iy; iDelta2++) {
      x->data[iDelta2] = X->data[iDelta2];
    }

    if (x->size[0] <= 2048) {
      k = x->size[0];
    } else {
      k = 2048;
    }

    if (2048 > x->size[0]) {
      for (i = 0; i < 2048; i++) {
        b_y1[i].re = 0.0;
        b_y1[i].im = 0.0;
      }
    }

    emxFree_real_T(&x);
    ix = 0;
    ju = 0;
    iy = 0;
    for (i = 1; i < k; i++) {
      b_y1[iy].re = X->data[ix];
      b_y1[iy].im = 0.0;
      iDelta2 = 2048;
      tst = true;
      while (tst) {
        iDelta2 >>= 1;
        ju ^= iDelta2;
        tst = ((ju & iDelta2) == 0);
      }

      iy = ju;
      ix++;
    }

    b_y1[iy].re = X->data[ix];
    b_y1[iy].im = 0.0;
    for (i = 0; i <= 2047; i += 2) {
      temp_re = b_y1[i + 1].re;
      temp_im = b_y1[i + 1].im;
      b_y1[i + 1].re = b_y1[i].re - b_y1[i + 1].re;
      b_y1[i + 1].im = b_y1[i].im - b_y1[i + 1].im;
      b_y1[i].re += temp_re;
      b_y1[i].im += temp_im;
    }

    iy = 2;
    iDelta2 = 4;
    k = 512;
    ix = 2045;
    while (k > 0) {
      for (i = 0; i < ix; i += iDelta2) {
        temp_re = b_y1[i + iy].re;
        temp_im = b_y1[i + iy].im;
        b_y1[i + iy].re = b_y1[i].re - temp_re;
        b_y1[i + iy].im = b_y1[i].im - temp_im;
        b_y1[i].re += temp_re;
        b_y1[i].im += temp_im;
      }

      ju = 1;
      for (j = k; j < 1024; j += k) {
        twid_re = dv17[j];
        twid_im = dv18[j];
        i = ju;
        ihi = ju + ix;
        while (i < ihi) {
          temp_re = twid_re * b_y1[i + iy].re - twid_im * b_y1[i + iy].im;
          temp_im = twid_re * b_y1[i + iy].im + twid_im * b_y1[i + iy].re;
          b_y1[i + iy].re = b_y1[i].re - temp_re;
          b_y1[i + iy].im = b_y1[i].im - temp_im;
          b_y1[i].re += temp_re;
          b_y1[i].im += temp_im;
          i += iDelta2;
        }

        ju++;
      }

      k /= 2;
      iy = iDelta2;
      iDelta2 <<= 1;
      ix -= iy;
    }
  }

  for (iDelta2 = 0; iDelta2 < 2048; iDelta2++) {
    if (b_y1[iDelta2].im == 0.0) {
      c_y1[iDelta2].re = b_y1[iDelta2].re / 2048.0;
      c_y1[iDelta2].im = 0.0;
    } else if (b_y1[iDelta2].re == 0.0) {
      c_y1[iDelta2].re = 0.0;
      c_y1[iDelta2].im = b_y1[iDelta2].im / 2048.0;
    } else {
      c_y1[iDelta2].re = b_y1[iDelta2].re / 2048.0;
      c_y1[iDelta2].im = b_y1[iDelta2].im / 2048.0;
    }
  }

  b_abs(c_y1, B);
  memcpy(&C[0], &B[0], 1025U * sizeof(double));
  for (iDelta2 = 0; iDelta2 < 1023; iDelta2++) {
    C[1 + iDelta2] = 2.0 * B[1 + iDelta2];
  }

  for (iDelta2 = 0; iDelta2 < 1025; iDelta2++) {
    f[iDelta2] = Fs * (double)iDelta2 / 2048.0;
  }
}

//
// Calculate the generalized cosine window samples
//  x is the length of the window
// Arguments    : double x
//                emxArray_real_T *w
// Return Type  : void
//
static void hannWin(double x, emxArray_real_T *w)
{
  int n;
  double anew;
  double apnd;
  double ndbl;
  emxArray_real_T *y;
  double cdiff;
  int k;
  emxArray_real_T *b_y;
  int nm1d2;
  if (rtIsNaN(x - 1.0)) {
    n = 1;
    anew = rtNaN;
    apnd = x - 1.0;
  } else if (x - 1.0 < 0.0) {
    n = 0;
    anew = 0.0;
    apnd = x - 1.0;
  } else if (rtIsInf(x - 1.0)) {
    n = 1;
    anew = rtNaN;
    apnd = x - 1.0;
  } else {
    anew = 0.0;
    ndbl = std::floor((x - 1.0) + 0.5);
    apnd = ndbl;
    cdiff = ndbl - (x - 1.0);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(x - 1.0)) {
      ndbl++;
      apnd = x - 1.0;
    } else if (cdiff > 0.0) {
      apnd = ndbl - 1.0;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&y, 2);
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = n;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  if (n > 0) {
    y->data[0] = anew;
    if (n > 1) {
      y->data[n - 1] = apnd;
      nm1d2 = (n - 1) / 2;
      for (k = 1; k < nm1d2; k++) {
        y->data[k] = anew + (double)k;
        y->data[(n - k) - 1] = apnd - (double)k;
      }

      if (nm1d2 << 1 == n - 1) {
        y->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        y->data[nm1d2] = anew + (double)nm1d2;
        y->data[nm1d2 + 1] = apnd - (double)nm1d2;
      }
    }
  }

  emxInit_real_T1(&b_y, 1);
  k = b_y->size[0];
  b_y->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)b_y, k, (int)sizeof(double));
  nm1d2 = y->size[1];
  for (k = 0; k < nm1d2; k++) {
    b_y->data[k] = 6.2831853071795862 * y->data[y->size[0] * k] / (x - 1.0);
  }

  emxFree_real_T(&y);
  k = w->size[0];
  w->size[0] = b_y->size[0];
  emxEnsureCapacity((emxArray__common *)w, k, (int)sizeof(double));
  nm1d2 = b_y->size[0];
  for (k = 0; k < nm1d2; k++) {
    w->data[k] = b_y->data[k];
  }

  for (k = 0; k + 1 <= b_y->size[0]; k++) {
    w->data[k] = std::cos(w->data[k]);
  }

  emxFree_real_T(&b_y);
  k = w->size[0];
  emxEnsureCapacity((emxArray__common *)w, k, (int)sizeof(double));
  nm1d2 = w->size[0];
  for (k = 0; k < nm1d2; k++) {
    w->data[k] = 0.5 * (1.0 - w->data[k]);
  }
}

//
// Arguments    : const double varargin_1[4]
//                const double varargin_2[4]
// Return Type  : boolean_T
//
static boolean_T isequal(const double varargin_1[4], const double varargin_2[4])
{
  boolean_T p;
  boolean_T b_p;
  int k;
  boolean_T exitg1;
  p = false;
  b_p = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 4)) {
    if (!(varargin_1[k] == varargin_2[k])) {
      b_p = false;
      exitg1 = true;
    } else {
      k++;
    }
  }

  if (!b_p) {
  } else {
    p = true;
  }

  return p;
}

//
// Arguments    : emxArray_real_T *idx
//                double Np
// Return Type  : void
//
static void keepAtMostNpPeaks(emxArray_real_T *idx, double Np)
{
  int loop_ub;
  emxArray_real_T *b_idx;
  int i15;
  if (idx->size[0] > Np) {
    if (1.0 > Np) {
      loop_ub = 0;
    } else {
      loop_ub = (int)Np;
    }

    emxInit_real_T1(&b_idx, 1);
    i15 = b_idx->size[0];
    b_idx->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_idx, i15, (int)sizeof(double));
    for (i15 = 0; i15 < loop_ub; i15++) {
      b_idx->data[i15] = idx->data[i15];
    }

    i15 = idx->size[0];
    idx->size[0] = b_idx->size[0];
    emxEnsureCapacity((emxArray__common *)idx, i15, (int)sizeof(double));
    loop_ub = b_idx->size[0];
    for (i15 = 0; i15 < loop_ub; i15++) {
      idx->data[i15] = b_idx->data[i15];
    }

    emxFree_real_T(&b_idx);
  }
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
// Arguments    : const double tsX[40]
// Return Type  : double
//
static double knn(const double tsX[40])
{
  double yfit;
  int iwork[712];
  int idx[712];
  int k;
  int i;
  static const signed char a[712] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  emxArray_real_T *Uc;
  int i2;
  int nb;
  int pEnd;
  int p;
  int q;
  int qEnd;
  int kEnd;
  int b_i2;
  int exitg5;
  boolean_T eok;
  static double x[28480];
  static double y[28480];
  int b_k;
  static const double tX[28480] = { 3.41758581347478E-7, -3.73853011436922E-7,
    7.94193180133536E-7, -1.51695397945949E-7, -5.88638007205537E-7,
    2.10615347916054E-8, 9.31928748239796E-7, 4.52103848979293E-7,
    -3.45154501461168E-7, 7.22698117305282E-8, -3.19617643504927E-8,
    -1.87303081596653E-7, 1.7342460386657E-7, -3.51378036654497E-8,
    -4.87088717602582E-7, 6.84055595064394E-7, -3.51459257589008E-7,
    2.98643033674179E-7, -3.83325491011277E-7, -1.73462308656191E-7,
    9.83829948551746E-8, 6.01254033212384E-7, 6.34344223551183E-7,
    -1.19071901576878E-7, -1.75660173403836E-6, 9.31781697109708E-8,
    4.11175339736361E-7, -6.29434986673917E-7, -1.02089267972595E-7,
    -3.88185779687104E-8, -1.9435083641075E-6, -1.79300493373602E-7,
    1.12228930687985E-6, -5.65434350525707E-8, -2.55012532255475E-7,
    -5.92733339645342E-8, -7.96957054840871E-8, -1.24164183602236E-7,
    1.8865584392115E-7, 3.31643411151701E-8, 1.84776044301267E-7,
    -1.4031129996089E-6, -7.43661522958991E-7, -1.39160025992671E-6,
    -6.47815892992094E-6, -2.08586061072724E-6, -1.75318725722575E-6,
    5.72968789237729E-7, -8.77932680210488E-7, -4.62869019008671E-7,
    -5.0882380385589E-7, -1.34696311197378E-6, -3.60164363478411E-7,
    1.5695952790503E-7, -4.77642741616376E-7, -4.68553277552346E-6,
    -1.23351452726702E-6, 4.6809774603368E-7, 2.59531711328329E-7,
    -1.25831120588678E-8, -6.55646052993802E-7, -1.00327575349692E-6,
    5.75484636202391E-7, 2.83607000586568E-7, 1.19285145212695E-6,
    4.94456815892879E-7, 1.57330582458045E-6, -1.18427519002042E-6,
    1.10975119744787E-5, 1.08259092090671E-5, -3.26530762055362E-7,
    6.81610270997693E-7, -5.73849911324069E-7, -2.2779217405836E-7,
    1.36974485581363E-6, -1.31128464772478E-7, -5.05627892815765E-7,
    2.9131837729525E-7, 1.21210791705137E-6, 2.24714814972526E-7,
    5.73146544831913E-7, -1.48885094370922E-6, -1.9610550772654E-6,
    1.54589857496828E-7, -5.19938139601191E-7, -2.80060754823305E-7,
    -4.4808527178394E-7, -7.80871273744878E-7, 3.20888458097897E-7,
    -5.76224295167178E-7, 4.73222993849056E-7, -4.20169208081224E-7,
    -1.45126647526452E-6, -2.94766646197402E-7, -2.85246771522821E-6,
    -1.06693161057714E-5, -9.62421279222713E-7, 1.04192837678945E-6,
    1.24636409747449E-7, 1.21275185871191E-6, 2.40996756781962E-7,
    -1.15929998999623E-6, -3.55191390186998E-7, -5.17652489909098E-7,
    9.39523399008572E-7, 2.75630912001392E-7, 1.29265751420238E-6,
    -8.60146661783756E-7, -3.05423209714767E-6, 1.23714326839866E-5,
    -6.00820170379874E-7, -5.85173348664024E-7, -7.82161643175577E-7,
    -4.17253747809525E-7, 1.95694773391337E-7, -6.82088262283156E-7,
    2.83952893646375E-8, 7.79469915722601E-7, -5.85472968550083E-7,
    -1.45626043952407E-6, 5.41519076960662E-6, 9.67411218276663E-6,
    -7.40988774514953E-7, 2.65489657794818E-6, -3.85377557284484E-6,
    -9.46439222973075E-8, -3.59944966212447E-7, -1.22513852956506E-6,
    -8.62051940814332E-8, -8.07124163572434E-8, -9.18805325244614E-7,
    6.03689685453994E-7, 8.4223115890572E-8, 5.70375748014523E-7,
    -6.15476324305069E-10, -8.56291118797572E-6, 2.93251824008453E-6,
    -1.58779446240779E-6, 1.29855816460075E-6, 1.1878773102139E-6,
    -2.50692494596309E-8, -1.6615827193305E-7, 9.31195552167648E-7,
    7.10782959951766E-7, 1.05800322420443E-6, -2.26140347134434E-6,
    1.123219840185E-5, -1.77749994484697E-6, -1.57699724770092E-6,
    -1.50510272046927E-7, 1.30213492825046E-7, -7.64310160568104E-7,
    2.10150249933307E-7, 6.48186174931065E-7, 2.21965867209897E-7,
    3.15909722606167E-7, 4.68370027990502E-7, -1.33496590611818E-6,
    1.13485380717683E-5, -1.01405898959236E-6, 8.18522829435867E-8,
    6.77169998651615E-7, -4.27161522987947E-7, -1.42032690138228E-7,
    -9.45303792356822E-8, 7.60397868098717E-7, 1.33187485300896E-6,
    -3.47657021721662E-7, -1.23107310057852E-6, 1.4252235962216E-6,
    -1.14904029721572E-6, -4.71532654340569E-6, 9.48349922037925E-8,
    2.75719274973276E-7, -2.72125739422531E-7, -1.5911346180864E-6,
    -5.82178171951542E-7, 3.07134639277691E-7, 1.87410709433356E-8,
    -2.004896948091E-8, -4.01032607373481E-7, -2.3468427343032E-7,
    -2.4107027762391E-7, -5.15817127879003E-6, 2.35694555223733E-6,
    -2.26743940676022E-7, -5.98386471320756E-7, 3.2649402594531E-8,
    9.54001397786338E-8, 3.642982100504E-7, -6.26820995172936E-8,
    1.49494411017911E-6, -1.81219863280028E-6, -1.512290351011E-6,
    9.96814569607451E-7, 1.51646242341356E-7, 6.92001655914916E-8,
    4.75280896453075E-7, -1.13789086189474E-7, -4.31634941613977E-7,
    8.22088196630166E-7, 5.88775511249256E-7, -2.01755486649016E-6,
    2.40183750982504E-6, -1.17092382896941E-7, -8.63294275961495E-6,
    -3.48341416863646E-7, -2.94824229324413E-7, 2.33492407609746E-7,
    -3.93491433520184E-7, -5.15716433486713E-7, 1.360191319277E-7,
    -3.26310662409092E-7, 4.63171651422925E-7, 7.10702067040629E-7,
    -2.2581583426703E-6, 7.82366095667484E-6, -5.29971261320825E-7,
    -6.34935961420051E-6, -3.5517218214688E-6, -4.46872140971938E-7,
    -7.39794016109189E-7, 7.74505474029717E-7, -1.58250952521956E-7,
    4.23328976233843E-7, 3.52521892148858E-7, 5.71796440444647E-7,
    -4.56269708722291E-8, 5.7679416688755E-8, -1.39783556868361E-6,
    -8.17375374675165E-6, 2.84155582004756E-6, -1.00654670421315E-7,
    2.89808083864586E-8, -8.67570492156805E-7, -2.34848358508759E-7,
    1.87476834434387E-6, -1.10488616210101E-6, 1.50722378127978E-7,
    -1.12541992316202E-6, -1.20070938518628E-6, 7.99869903873654E-6,
    -1.95900887440595E-6, 5.82143717010784E-7, -4.41464856619368E-7,
    -8.72822885047768E-7, 1.10194978627705E-6, 6.41393522175969E-8,
    6.11263857091478E-7, -4.76396505842436E-7, 2.19928136899103E-7,
    -6.32776326957134E-7, 2.24987525289836E-7, -1.89561680573892E-6,
    -6.10079970788766E-6, -8.57018335730253E-9, 1.30521868477E-7,
    -3.71839910469275E-7, -3.9983384945859E-8, -6.79893581755988E-7,
    -7.95678301591032E-7, 6.02283878286187E-7, -1.21393092056851E-6,
    -2.3233584799416E-6, 1.70814821657539E-6, 6.14935202444982E-7,
    -7.79206656337671E-7, 1.84154164899236E-8, 7.46404280490722E-7,
    -3.52227612305713E-7, -3.59723527423396E-7, 8.62404477198034E-7,
    -1.38560780655557E-6, 5.46846191720777E-6, -8.83487566452224E-7,
    -1.00304379985896E-6, 1.24823621903476E-7, 4.87873190522374E-7,
    1.86012374231358E-7, 9.42374682877385E-7, 2.80245001848939E-7,
    -2.22275830101048E-7, 4.29158252304581E-6, 4.41224648508458E-6,
    -3.97558378342121E-6, -2.02220215831996E-7, -2.86015053641937E-7,
    6.43775164077289E-7, -2.04137961472868E-7, 9.68208480039883E-7,
    -9.44937468702337E-8, -1.43643067126962E-6, 1.24520738957912E-6,
    -2.88141284152792E-6, -4.50954735658379E-6, -1.30834659451147E-6,
    -6.69270879786552E-7, 3.39474752099534E-9, -7.69103951599065E-7,
    9.21445961606906E-7, 2.93236964815356E-7, 6.11618456371054E-7,
    -6.19274981822427E-7, -7.84064304128755E-7, 4.8346432243001E-6,
    -1.61082769944491E-6, -6.90872150185503E-6, 9.76386545072582E-8,
    -9.89305474009422E-7, -6.51295108409342E-7, 1.11232859937237E-6,
    -9.21325976930584E-7, -9.58756208969248E-8, -4.44415948111284E-7,
    4.43472869942043E-7, -1.2479922093223E-6, -4.53050128290494E-6,
    8.37250975506191E-6, 8.77572831078434E-7, 3.90114109445702E-7,
    -5.18502676853484E-7, -6.60690147762972E-7, -7.54892109957456E-7,
    9.81662203923645E-7, -1.18680306941435E-6, -4.43827892533626E-7,
    -6.24145172674901E-7, -5.4092084460273E-6, 7.95220635249705E-6,
    -4.21075834787185E-7, -1.58446823316492E-7, 1.4087352766369E-6,
    6.24866204835764E-7, 2.54524408676655E-7, 2.44893539622806E-7,
    2.06975639743121E-7, 1.15657225420601E-6, -2.06563828474468E-6,
    7.85439103432026E-6, 2.76814745501405E-7, -4.84557223107628E-6,
    4.60996110503665E-7, -1.03894507511968E-6, 5.08681437844455E-7,
    7.35499901082259E-7, 9.89779168134081E-7, 7.91359757626456E-7,
    2.89691693356478E-7, -3.65824962618441E-7, -1.55167211716771E-6,
    -4.03209946201928E-6, 9.62661164284273E-6, -2.46617977317214E-6,
    -1.41362581953715E-6, 1.82081743758296E-6, 1.24826324874688E-6,
    -1.62737265114843E-8, -1.51556884573084E-8, -4.5530957780161E-7,
    1.53903018211316E-7, 7.14969712696585E-7, -8.18102675925234E-7,
    7.59288950069738E-6, 7.6948656194552E-6, -3.73913653303511E-6,
    1.31644072013248E-6, -1.70687811978242E-6, -1.02010149226415E-6,
    1.80273653197096E-6, -6.63455708921713E-7, 1.57656111869046E-6,
    -3.33147346738586E-7, 1.10805009098998E-6, -1.82594880511295E-6,
    2.0173874435858E-6, -3.82815646207549E-6, -2.64181435404068E-6,
    -1.95976551306109E-6, 5.37711200020275E-7, 8.50432116883041E-7,
    4.67310124127081E-8, -3.90949956685601E-7, -4.57820256928373E-7,
    -4.84029415898358E-7, 3.25437744260796E-7, -8.00453641604182E-7,
    8.4052884422919E-6, 1.89842055261312E-6, -2.92876752271134E-6,
    1.88024882963422E-7, -5.5087381219557E-7, 3.23414094695325E-7,
    -4.75129331737217E-8, -3.82061923462785E-7, -5.62384194794364E-7,
    -2.66238493733549E-9, 1.69196358154865E-7, 1.47719817151961E-8,
    9.8262568557099E-7, -9.95596984715836E-7, -1.71845208666931E-7,
    2.36126520504465E-7, -8.67646026402237E-7, -5.39400290746964E-7,
    -2.55563682898892E-7, -1.0742180223866E-6, -3.55740213264254E-7,
    1.07665361945784E-7, 2.36290385042989E-7, 5.8691767959808E-7,
    -9.66577595526007E-8, -6.35329734435515E-9, -2.80122457945033E-7,
    -4.61892938222625E-8, -2.64360577499178E-7, -4.45668832038028E-7,
    -6.92894649493627E-8, 4.83429879627724E-7, 2.97895340390396E-7,
    7.77540405790845E-8, 3.74314950615489E-7, 4.90459302692555E-7,
    -4.49765052182671E-7, 3.0528260357924E-8, -8.31763939754096E-7,
    -2.78315956004066E-7, -2.25310999379682E-6, 9.88220718599157E-6,
    2.88152568929944E-6, -3.18168533698817E-6, -2.5585887719694E-6,
    -9.31629789042777E-7, 9.83310824590118E-7, -1.16901224750623E-6,
    -1.97902640083012E-7, -1.40980230378951E-6, 2.45656145380085E-7,
    1.07930687232373E-6, 3.35785318861931E-7, 9.75140360515448E-7,
    1.16737706566678E-6, 7.06069245451397E-7, 4.35949727069427E-8,
    -8.30455340813526E-7, -1.07795603910637E-6, 1.78411761028131E-7,
    -1.06173248662461E-6, 8.41564707487816E-6, 1.45820512883933E-6,
    1.11412446056484E-6, -3.7053743139656E-7, 1.70587757923938E-6,
    9.94204061472848E-8, -8.85352955134848E-7, -6.67131359278815E-7,
    1.09435644390501E-6, -2.12847236124602E-6, 1.07484533394236E-6,
    2.10597203000018E-6, -2.97982470827037E-6, 2.76765488782477E-7,
    -3.25569875395796E-7, -5.30711124437643E-8, 4.75128142491469E-7,
    -3.8840734962032E-7, 1.92721831691382E-7, 1.51449985482401E-6,
    -2.95672554590052E-7, -5.05562497913281E-9, -5.93886165071791E-7,
    2.00042067459989E-6, 6.17384641858189E-7, -5.43247476146684E-8,
    1.26913371963216E-7, 7.34666943785848E-8, -2.5965863057066E-7,
    1.52129364613048E-6, -1.11283056348493E-6, 3.97983146147256E-7,
    3.26108085780431E-7, 1.25748062354302E-6, -8.65830864087449E-7,
    7.68367645444263E-6, 1.01611048529938E-5, 9.33709109241616E-7,
    -9.15141488078784E-8, 4.39516301676725E-7, 5.73935854371017E-7,
    1.73885112894941E-6, 9.79258347767084E-8, -1.32460929253173E-6,
    2.38444880514433E-7, 2.6345669340132E-7, -8.19587115830498E-7,
    8.09309973809569E-6, 5.57823006901977E-6, -5.51841196552641E-7,
    8.36154427591295E-7, -5.95487482870482E-7, 2.24949731007072E-7,
    -5.38288447693826E-7, 1.01586535270111E-6, -3.01065024887249E-7,
    3.89651373424864E-7, -1.43396081102421E-7, -5.92026638787574E-7,
    -9.37113160380664E-6, -5.13569717375334E-6, -7.80369517167185E-7,
    1.27764170345973E-6, 4.20582171071963E-7, -3.40950747425915E-7,
    -3.51015316744511E-8, 1.53218242886678E-6, -7.25846070991076E-7,
    1.46906466845823E-6, -4.48057828087328E-7, -1.82706543261079E-6,
    1.90877942022816E-6, 5.88488936179829E-7, 7.95875965096526E-7,
    8.25860742397859E-7, -3.40553860382465E-7, -2.20300978678151E-8,
    -2.05724371880376E-7, -8.85110379902947E-8, -8.71346329617082E-9,
    -9.93065612442474E-7, -3.03125700759114E-7, -6.80080025071354E-7,
    -1.46855525546586E-7, 7.86194430172217E-8, -1.27315578527639E-6,
    6.87636409507246E-6, 4.60121683294168E-7, 8.87346721366158E-7,
    1.69592559825775E-7, -1.06898633760851E-6, 3.20625936836356E-7,
    -6.00788584305324E-7, 5.40656281251196E-6, -2.76817779731583E-7,
    6.26074590809651E-7, 2.52130424984408E-7, 2.97706261064433E-7,
    1.45716390770917E-6, 3.08152141330984E-7, -4.9347696609788E-7,
    3.88410597391895E-7, -1.49706362013683E-7, -4.15961517603735E-6,
    9.91299737226995E-6, -1.00664465464865E-6, 3.48718726508172E-7,
    -8.87682923406452E-7, 3.7308810614572E-7, -6.01506800721516E-7,
    -5.55838700318703E-7, 4.79085534015419E-7, -8.09127624033625E-7,
    -1.87120567036787E-6, 1.30388413368236E-5, 8.78217464443972E-7,
    1.5566931669705E-6, 1.24386216133855E-8, -2.85731164165451E-7,
    -4.69093321305006E-7, 7.18718200077887E-7, 7.51243436313539E-8,
    1.34046620496463E-6, -1.77800034236625E-6, 5.67217640994241E-6,
    1.66545541657326E-6, -4.37831335299469E-7, 5.73126203066308E-7,
    1.28341365458239E-7, 6.7333395971399E-7, 1.07029880560163E-7,
    1.3946107919239E-6, 5.10653783205173E-7, -7.59367457587263E-7,
    5.36169241204201E-6, -1.00351104926258E-6, -4.4412705153131E-6,
    2.23768729890495E-7, 3.38294908648856E-7, -9.29245456250725E-7,
    1.14159512770609E-6, 1.3268001582744E-6, -1.3968824849299E-7,
    1.63913856925057E-6, 2.39724908755935E-6, -5.953968135332E-6,
    1.37848416667371E-7, -1.18014651462101E-6, 8.26125031640555E-7,
    -7.67431551904669E-7, -4.30640329287084E-7, -1.16972542572264E-6,
    1.1774747218925E-6, -4.14668971637033E-6, -2.98365376024127E-6,
    5.99735147183616E-7, -2.37869990369863E-6, 4.05859720910376E-7,
    3.40530025200837E-7, -5.1949770631191E-9, 1.15310398829893E-6,
    1.97159832997263E-7, -6.86769783804525E-7, -2.00625295680939E-7,
    2.57729784017981E-7, -8.1668250631855E-7, -1.27309608458157E-6,
    1.48247971266042E-6, -2.07215021161651E-6, -8.78440147860788E-7,
    -1.17714380217921E-6, -2.03159740122242E-7, 4.7744599168633E-7,
    -5.44780287299377E-7, -6.67616475614027E-6, -8.37511482425393E-7,
    1.18426864726926E-6, 1.85678861718497E-7, 6.26272405341613E-7,
    -7.00807736146094E-7, -9.56160269584605E-8, 2.55889201934774E-7,
    5.81319997032456E-8, -8.30721193819276E-8, -1.2570485264891E-7,
    9.61994916678507E-7, 7.81637570948292E-7, -1.12230338559124E-6,
    -3.62625229662508E-7, -1.00647466852256E-6, -8.92951770422288E-7,
    1.08086186499048E-7, -2.69936060141736E-7, -1.10203194488824E-7,
    2.21693456403357E-7, -5.07020102265263E-10, 5.91371013743667E-8,
    -7.22964722124555E-8, -8.73283017578545E-7, -7.91628442631194E-9,
    4.31210128170662E-7, -9.62960371841711E-7, -4.91415173333036E-7,
    -2.89164625569235E-6, -6.65665189047362E-7, -6.72806234244423E-8,
    -2.13096719500749E-7, 1.59426664627891E-7, -1.17075998572271E-6,
    -8.98577788235011E-6, -9.28134615504654E-6, -7.01337719997126E-7,
    4.10927191765328E-8, 1.10101403955647E-7, -7.8814542867644E-7,
    -1.167598646461E-7, 4.45100492009447E-7, 9.82487047900416E-8,
    -1.56066235953631E-7, 1.88786167775912E-7, -6.86778483402234E-7,
    4.49820726393985E-7, -1.35858475114561E-6, -5.93122860540669E-6,
    -2.60296568648478E-6, -3.59442178518445E-7, -1.40865489686478E-6,
    3.0123236396972E-7, 7.80970223250501E-8, -4.89693271027284E-7,
    5.09510150165489E-7, 4.81194416995403E-7, -3.82603930748131E-7,
    -6.07909506948315E-7, -6.67177717469766E-7, -4.55981428613748E-7,
    -1.74587626438279E-6, -6.41014030338652E-6, 5.84936354642605E-6,
    -4.10517026882251E-7, 4.93386557074611E-7, 3.31784192812869E-7,
    -4.34436139806086E-7, -5.24223588575897E-8, -4.86909206286241E-7,
    -6.06706890423337E-7, -1.08000229203287E-7, -4.36110487415131E-8,
    3.27997835472806E-7, -7.82417584499185E-7, 6.24578071257475E-7,
    -2.16635777376267E-6, 5.7332230074683E-7, -8.28365388237837E-7,
    -5.95137454506616E-7, 2.14998792696261E-7, -5.32494928411371E-7,
    -6.05373670276305E-7, -6.27253643018388E-8, -1.32495430944984E-6,
    -6.76034591375343E-6, -6.47400265818212E-7, 7.91843470154805E-7,
    -1.20306196816929E-7, 3.91583606855573E-9, 3.37163525349492E-6,
    3.01581675280723E-6, 4.45303619898528E-6, 3.8746583894074E-6,
    4.27128303279192E-6, 3.71202644717309E-6, 4.72364216504853E-6,
    3.07891817753881E-6, 3.0540469180255E-6, 2.75197215453608E-6,
    3.03528448873142E-6, 3.33562746372562E-6, 3.1391175702738E-6,
    3.04845802569086E-6, 3.83286254508844E-6, 4.23384825777839E-6,
    3.25024107281396E-6, 3.13735095748544E-6, 2.65820705050366E-6,
    3.48763508949945E-6, 8.18440605094934E-6, 9.29196749843436E-6,
    9.71921590996643E-6, 9.96798280995384E-6, 7.43482368632031E-6,
    3.48027969836685E-6, 4.23020878915733E-6, 1.85676048126266E-5,
    5.58594305739654E-5, 5.84037832093669E-5, 5.84785074333519E-5,
    2.86765371377512E-5, 7.43594887577444E-6, 6.58513631255477E-6,
    5.00028977026814E-6, 4.25512998294096E-6, 3.18712980463183E-6,
    2.86539015512444E-6, 2.99800002136049E-6, 2.83855820003086E-6,
    3.92061016480683E-5, 6.34489869638103E-5, 6.58589141604876E-5,
    6.6267317102111E-5, 2.63988154487295E-5, 9.85834496015902E-6,
    8.30112895139651E-6, 5.07701522963831E-6, 6.38332065096207E-6,
    5.33750173418382E-6, 5.24001519998768E-6, 4.62317521658038E-5,
    6.0753761772772E-5, 6.27308743408122E-5, 6.17090594688315E-5,
    1.91706412481446E-5, 7.53787042371478E-6, 4.66579429898472E-6,
    4.46540794123315E-6, 4.49600961261206E-6, 5.49601221731849E-6,
    5.01867794189886E-6, 4.24709121146971E-6, 3.65966525578086E-6,
    8.58627515878571E-6, 5.92387936200976E-5, 9.47070878675335E-5,
    9.46282816691573E-5, 9.6541294950205E-5, 6.62639016354164E-5,
    7.80165011809295E-6, 7.27878924494435E-6, 5.01537791352341E-6,
    3.3195183789682E-6, 7.82123939753714E-6, 3.3349950319776E-6,
    4.56231751620005E-6, 4.56207165576638E-6, 2.1715575835302E-5,
    7.22312566552709E-5, 9.13426704883417E-5, 9.19323824141281E-5,
    8.04203433047752E-5, 2.93587956758375E-5, 4.74909600550067E-6,
    3.70426932999475E-6, 4.22699558214942E-6, 4.43503459377251E-6,
    4.80908700497788E-6, 5.57724400217702E-6, 4.94519392035696E-5,
    8.70713328000233E-5, 0.000108893252033073, 0.000110762236882951,
    8.31114313774344E-5, 4.34272230792535E-5, 6.37769318494896E-6,
    7.17337250914136E-6, 4.98721619275946E-6, 6.34851067682111E-6,
    3.52704060566423E-6, 6.81210911749086E-6, 4.33531802262118E-6,
    4.07047977105945E-6, 1.76790711103559E-5, 5.22886709819686E-5,
    9.22072598461378E-5, 9.74901993755876E-5, 8.25473781322544E-5,
    7.33153910468797E-5, 4.58663799694062E-6, 3.7669960035485E-6,
    4.21701918703968E-6, 3.56061344007076E-6, 3.19092101559162E-6,
    4.98005734115965E-6, 1.67338857109079E-5, 5.62290337079925E-5,
    7.86212191217394E-5, 7.85740717831668E-5, 7.81951425483408E-5,
    7.39224185919122E-5, 5.38958717789279E-5, 4.96683220882142E-5,
    1.46779064107276E-5, 3.24034781319626E-6, 3.63506936458615E-6,
    6.15152510886289E-6, 3.59477974718692E-6, 3.70253501897745E-6,
    5.68602293425476E-6, 5.12005400182739E-5, 7.49489013983178E-5,
    7.79741681904094E-5, 7.49848295831071E-5, 3.23406236593305E-5,
    1.30291577925706E-5, 7.53260514227871E-6, 7.94765109607022E-6,
    5.94701834633236E-6, 4.12293992368427E-6, 3.5434292100504E-6,
    1.51839919800075E-5, 5.96007292179343E-5, 7.11446622760433E-5,
    7.18705776063672E-5, 6.29937041141338E-5, 7.72008894947328E-6,
    7.37900149535339E-6, 4.07911365750321E-6, 4.96261656702992E-6,
    5.67228758671878E-6, 4.51960040474236E-6, 5.05420730257276E-6,
    5.75388427611822E-6, 6.09252746379312E-5, 6.4453519599096E-5,
    6.45595699594184E-5, 6.39718875069051E-5, 6.37583977777968E-6,
    3.93449657351325E-6, 4.91368023142723E-6, 4.30537741586828E-6,
    3.97063117563471E-6, 4.02848078725021E-6, 5.3596951993463E-6,
    4.74867439985735E-5, 8.19712058406747E-5, 9.52668077460041E-5,
    9.66580917160873E-5, 6.89437605324935E-5, 2.17559380159031E-5,
    6.55589086849598E-6, 5.7967969375906E-6, 4.47342663234447E-6,
    8.73183897739809E-6, 4.79351329721744E-6, 4.25599693212427E-6,
    3.61869019735306E-6, 3.91542451449723E-5, 5.66602479505313E-5,
    8.52008660589415E-5, 8.79721128675155E-5, 6.76172273073046E-5,
    3.49187670661104E-5, 4.32786338935498E-6, 4.93572636393766E-6,
    3.94122145722346E-6, 4.03622932949329E-6, 1.74197237691461E-5,
    5.26507039299363E-5, 8.05828223168701E-5, 8.83516026455662E-5,
    6.90917689702969E-5, 6.12308773887523E-5, 6.90847420853792E-6,
    3.86239026398502E-6, 5.20937074614986E-6, 3.7662925654842E-6,
    1.99680711679466E-5, 5.85512744209781E-5, 8.03790300101342E-5,
    0.000100943337428896, 8.5436103235991E-5, 7.04453849940902E-5,
    3.33804763450414E-5, 2.8863228868149E-6, 3.80291264456092E-6,
    3.04673971082865E-6, 3.50060382041115E-6, 3.90908466544871E-6,
    2.51559462999409E-6, 1.62484227224573E-5, 5.13581677782302E-5,
    5.63814820306373E-5, 6.87396008360067E-5, 8.79538536317799E-5,
    8.53587072470711E-5, 9.08377644800055E-5, 3.2196500758819E-5,
    5.40355853399902E-6, 5.53354325195075E-6, 5.81910968066078E-6,
    3.99006840774014E-6, 3.92386122355911E-6, 3.34487583242198E-6,
    3.23674858236943E-5, 6.66964587564481E-5, 9.32684368929753E-5,
    9.55520539292451E-5, 8.05267361352616E-5, 3.21823278723367E-5,
    6.74769272715504E-6, 4.08484073562941E-6, 5.5527380045684E-6,
    3.92520247309651E-6, 8.1976298078072E-6, 5.67621274229331E-6,
    5.10561315079101E-5, 6.56149052992259E-5, 6.54056785538196E-5,
    5.83711730070115E-5, 8.90510044995789E-6, 5.19690720596883E-6,
    4.63501503057897E-6, 4.11116258481624E-6, 5.04483327328424E-6,
    3.12727417646414E-6, 3.71069018139323E-6, 4.34296624503952E-6,
    2.45765436596385E-5, 5.74355284733192E-5, 5.96718851879302E-5,
    5.93898919252666E-5, 2.5901583760866E-5, 3.67722842962691E-6,
    3.74402976876796E-6, 4.05231050199874E-6, 4.50204205811439E-6,
    4.72804494363647E-5, 5.54660384162055E-5, 5.72936516404795E-5,
    5.59964620508903E-5, 1.35058607448242E-5, 8.35860454498854E-6,
    4.62817951456306E-6, 4.94590669335324E-6, 3.32998593096712E-6,
    4.56748502745941E-6, 1.54741382924574E-5, 5.12966437842638E-5,
    5.36049999607409E-5, 5.33793793665971E-5, 3.97912190433198E-5,
    4.97824350111803E-6, 5.95121686241188E-6, 4.18918193194976E-6,
    4.83557196589726E-6, 4.15599383195075E-6, 4.1280999555787E-5,
    6.31194906731089E-5, 7.78517844583054E-5, 7.96713287470517E-5,
    5.64078571582069E-5, 1.59177361279461E-5, 5.28678481073176E-6,
    4.59754344423909E-6, 4.44264977914422E-6, 4.44957218108953E-6,
    3.71413608937797E-5, 4.83108936262014E-5, 7.91928802966943E-5,
    7.96845899521593E-5, 6.23646659468745E-5, 2.30882111350936E-5,
    7.2199140549518E-6, 4.93532279477704E-6, 4.99372026840007E-6,
    5.39779020700708E-6, 5.82358352062526E-6, 3.43408822532877E-6,
    4.37028088405794E-5, 6.86478408721837E-5, 8.05530972789871E-5,
    8.20611214546861E-5, 6.04638734198067E-5, 2.64426890564119E-5,
    4.39639062565601E-6, 6.45651035250395E-6, 6.34365219810398E-6,
    7.80615380737666E-6, 7.47369094390601E-6, 4.21886003028141E-5,
    5.20919931921198E-5, 8.59854686557E-5, 9.42416133812019E-5,
    7.52829220074292E-5, 5.35879940297277E-5, 7.54023523605286E-6,
    3.72474822292296E-6, 3.93018940884495E-6, 4.37843157503942E-6,
    4.39158453386594E-6, 3.85177847849336E-5, 6.41125307235727E-5,
    8.74051850628904E-5, 9.29979129330508E-5, 7.34512428020644E-5,
    5.5598794212907E-5, 6.69789021907029E-6, 5.37183662919365E-6,
    6.174647675832E-6, 4.39170219867229E-6, 2.2879961839034E-6,
    1.57642972776475E-5, 6.03097241097034E-5, 8.54836102999401E-5,
    0.000102819904476094, 9.82148434579696E-5, 6.77093737041269E-5,
    2.08517082724527E-5, 5.59297205302429E-6, 5.24799108963945E-6,
    3.26236502728331E-6, 5.27310540821475E-6, 5.7390908991061E-6,
    1.94648536805718E-5, 5.21267867949912E-5, 7.83674884094672E-5,
    8.86793677861688E-5, 7.44864938193589E-5, 5.78210515560198E-5,
    1.20683491801492E-5, 6.89095768005079E-6, 8.75785792381548E-6,
    6.91431268925217E-6, 4.7283512116209E-6, 4.08502828457546E-6,
    7.94310499164186E-6, 5.93510808735495E-5, 8.32935205323863E-5,
    8.83469341463616E-5, 9.47249643086349E-5, 5.45620864262069E-5,
    1.75105865761416E-5, 9.76451085848904E-6, 8.57710632570575E-6,
    6.07597076039409E-6, 9.1294998224395E-6, 4.46379502449207E-6,
    2.01107671643799E-5, 5.67629902267206E-5, 7.63861711151511E-5,
    8.45587247279736E-5, 7.57612922666506E-5, 6.79830051443676E-5,
    1.16935975625059E-5, 8.97434901172726E-6, 4.07108639145029E-6,
    4.20240090960925E-6, 2.59657305378068E-6, 2.85716429119251E-6,
    7.01417541100103E-6, 4.76342004085071E-5, 7.71246173934482E-5,
    9.03605980766655E-5, 9.41437838816041E-5, 6.26921413109758E-5,
    1.12527240295339E-5, 3.47568901632636E-6, 3.18566526959147E-6,
    2.68380977984381E-6, 2.2954450791234E-6, 2.62348932973685E-6,
    4.96409987616636E-6, 6.40024144279805E-5, 4.46672496926269E-6,
    3.94530982396823E-6, 5.65804052831081E-6, 5.37671680619382E-6,
    4.20703925258726E-6, 3.89670982445196E-6, 5.27523416818513E-6,
    3.90217243495334E-6, 3.75552150163045E-6, 5.6420289779641E-6,
    3.68968340541753E-6, 3.52796375690149E-6, 4.43128813029266E-6,
    5.82006281062052E-6, 5.17291818609342E-6, 5.08683400772124E-6,
    3.781159405475E-6, 2.96724683821109E-6, 3.35381742751056E-6,
    3.58668992816509E-6, 3.56038035708736E-6, 4.05128614018993E-6,
    3.70850274505659E-6, 9.93161575263053E-6, 8.41415243937358E-6,
    8.97310128704988E-6, 8.72518935379435E-6, 1.10668911647379E-5,
    6.0895973709392E-5, 6.33288484562869E-5, 6.41438227522182E-5,
    5.82901096374134E-5, 1.24409366298288E-5, 1.32765144967705E-5,
    1.08601447505259E-5, 4.52228845923032E-6, 4.31897511655581E-6,
    5.36603790151625E-6, 3.23761846274879E-6, 7.43956775911787E-6,
    5.71536846987773E-6, 6.24223210013693E-6, 4.90157474172399E-6,
    5.09031257012283E-6, 6.99391586882541E-6, 4.7097198543778E-6,
    3.21588671268631E-6, 1.55227210854236E-5, 5.40870325606285E-5,
    5.98998118690263E-5, 6.04871800682206E-5, 5.63372601418881E-5,
    8.35304815684776E-6, 6.35589383357957E-6, 4.44700586484314E-6,
    8.23570004665506E-6, 5.11909339668164E-6, 5.83459668321736E-6,
    4.61090473316341E-6, 4.90349059266461E-5, 6.2027307614954E-5,
    6.26895264760114E-5, 6.03262622737667E-5, 1.30253816913883E-5,
    4.26344090477247E-6, 4.67554607351318E-6, 3.69605592668675E-6,
    4.82989342713673E-6, 3.47734458359234E-6, 4.80852169095003E-6,
    7.75018969101321E-6, 3.5639537365009E-5, 4.33739844343395E-5,
    4.34624830395193E-5, 4.20073791283726E-5, 8.3302224529316E-6,
    4.2391085893354E-6, 4.88601979733043E-6, 4.7758155767138E-6,
    4.24812782469994E-6, 7.28804213562255E-6, 6.13743054822089E-6,
    4.98842518468393E-6, 5.90324152786321E-5, 9.32899520541292E-5,
    0.000100284059607571, 0.000102468126823007, 6.97911898828033E-5,
    5.4381118118542E-6, 2.90584119377996E-6, 4.48068496156864E-6,
    4.7849803243442E-6, 7.99181484852237E-6, 3.29132597667133E-6,
    8.10341046358841E-6, 5.12957964222148E-5, 8.6800087334088E-5,
    8.89778823397188E-5, 9.23238250918876E-5, 6.313506760734E-5,
    3.76583750708782E-6, 4.61937989742428E-6, 3.7161609612337E-6,
    4.67759908589466E-6, 4.66518521062677E-6, 5.62742492555296E-6,
    2.59365265891999E-5, 5.16324202426295E-5, 7.60426386013154E-5,
    7.62504524394262E-5, 7.12889033635945E-5, 3.1440265282605E-5,
    4.52388730026132E-6, 7.28639121638776E-6, 4.28748967832917E-6,
    4.9652773042126E-6, 4.62337996543307E-5, 8.10007355890244E-5,
    8.97500094045253E-5, 9.01173915932007E-5, 6.84211001725633E-5,
    8.55076747402762E-6, 1.0924209482417E-5, 6.86593984112228E-6,
    6.82572146986911E-6, 5.0715193627809E-6, 4.12522651950866E-5,
    5.17342381402919E-5, 5.24615372313896E-5, 5.23283891071414E-5,
    3.30410306458216E-6, 6.17276803826912E-6, 4.04792563674224E-6,
    1.19398254214043E-5, 5.70579091829734E-5, 5.93508607156272E-5,
    5.91865203598832E-5, 5.28328454510861E-5, 4.01180669013824E-6,
    5.66421242022991E-6, 4.76570596068168E-6, 4.54232875771246E-5,
    4.70013897577599E-5, 4.76865111722796E-5, 4.67014431437146E-5,
    5.31029254566936E-6, 4.7473881094547E-6, 3.53360228023984E-6,
    6.44535042409817E-6, 8.34139370167633E-6, 1.76536613578107E-5,
    5.57157289404749E-5, 9.35979386507648E-5, 9.55929451676886E-5,
    8.63399073417559E-5, 6.13066051788653E-5, 4.5593171811243E-6,
    3.52106545486869E-6, 6.08504573263962E-6, 5.726241563008E-6,
    1.7140558950509E-5, 4.66464812227281E-5, 9.02109966445692E-5,
    9.33731217822867E-5, 8.55991200982777E-5, 7.4608973857491E-5,
    8.67395416293994E-6, 6.79708803086494E-6, 4.89017317133035E-6,
    6.10604098699005E-6, 5.18245992849081E-6, 5.72616576970111E-6,
    4.19368950223514E-5, 7.47935349563318E-5, 7.87564679806408E-5,
    7.92141206044978E-5, 5.91888340929517E-5, 5.15075183480319E-6,
    4.67708602414239E-6, 5.59610441218786E-6, 5.61015232977352E-6,
    4.73455514657458E-6, 6.33587666491987E-5, 0.000102041528138027,
    0.00010960801940981, 0.00010865994381664, 7.03867648312873E-5,
    1.94092926757091E-5, 7.97626776873351E-6, 5.81646008262735E-6,
    7.61165635633457E-6, 5.86280052426037E-5, 0.000102334174488706,
    0.000114630836888341, 0.00011421465086033, 7.72649846562451E-5,
    2.34229848275438E-5, 5.11257204988104E-6, 6.85703510934611E-6,
    6.04709372176617E-6, 4.78314463464263E-5, 6.3005781812902E-5,
    9.64883805916075E-5, 9.90344605441979E-5, 7.95140163718153E-5,
    3.33144380602005E-5, 3.24636291466372E-6, 1.00412514276337E-5,
    3.33623389079258E-6, 3.250091783573E-6, 3.60003183200257E-6,
    5.67546337082403E-6, 4.10280733394053E-6, 4.53819279412314E-5,
    5.43972635039153E-5, 5.63894368523585E-5, 5.62218842160451E-5,
    1.18422589386688E-5, 1.10516090213679E-5, 1.32885084662112E-5,
    5.81083017583013E-6, 3.84053478874889E-5, 5.50567189300989E-5,
    5.63230417845088E-5, 5.54351950610155E-5, 2.63918098206866E-5,
    6.72127286140669E-6, 7.28548111088211E-6, 3.51823242013472E-6,
    4.21596144992426E-6, 5.10161579159281E-6, 4.92895625705124E-6,
    4.28552190886687E-6, 7.93526279386593E-6, 9.68461237229735E-6,
    9.64665132378902E-6, 9.83936885700595E-6, 7.73155301610167E-6,
    8.4200512175813E-6, 5.51169982822563E-6, 7.39175511063597E-6,
    5.54898767651852E-6, 4.36579873237283E-6, 4.50333124189084E-6,
    4.40253796055544E-6, 5.06249894649884E-6, 6.60440650793578E-6,
    9.62866224617853E-6, 1.18646919988236E-5, 1.13433939740239E-5,
    1.55294235049822E-5, 2.41123193777001E-5, 2.50062613008584E-5,
    2.61516820626351E-5, 1.16750224947607E-5, 5.73861827348179E-6,
    3.59641973263232E-5, 8.36333681571591E-5, 9.65712359207826E-5,
    9.78489191855004E-5, 9.56324677956578E-5, 3.7604795441721E-5,
    5.28285273882613E-6, 4.5078324293255E-6, 4.57660567871555E-6,
    4.91328465619337E-6, 3.11311897547725E-6, 3.20679080555035E-6,
    3.71088307115406E-6, 2.5011929716079E-6, 2.11444457574763E-5,
    7.64999801464471E-5, 9.79578581191299E-5, 9.92447931006451E-5,
    9.00654913488614E-5, 3.06194162323421E-5, 3.40290482750228E-6,
    6.38769135809911E-6, 4.04580921411238E-6, 4.02909933803313E-6,
    3.73196826755433E-6, 3.66261259620457E-6, 3.65243896571208E-6,
    3.02325696497964E-6, 1.30098559038408E-5, 4.46822461850944E-5,
    7.14041121304018E-5, 7.3350683535116E-5, 6.6383184362028E-5,
    4.12297362147972E-5, 4.88027299084886E-6, 4.41296727452144E-6,
    4.73795961022973E-6, 4.24920359963017E-6, 4.29116739476273E-6,
    4.34579623465107E-6, 4.50929097487941E-6, 2.9708220229542E-6,
    4.53703097761524E-5, 5.79480858685112E-5, 5.88961971168723E-5,
    5.84605409156671E-5, 8.88290218850578E-6, 3.64101297840114E-6,
    4.78008802488361E-6, 4.22638544300582E-6, 4.11850554533426E-6,
    2.74627378037872E-5, 5.43710193016393E-5, 5.59315770379511E-5,
    5.59419902104292E-5, 2.72243375269983E-5, 5.02659696855622E-6,
    4.77969530302537E-6, 4.6331790863727E-6, 5.54303554578287E-6,
    1.0687503051201E-5, 6.75946706233107E-6, 1.52886508628176E-5,
    7.81677883437481E-6, 8.30834314372461E-6, 8.10738269299842E-6,
    1.51969827077963E-5, 1.09438558801172E-5, 6.0384743384092E-6,
    5.94819722791181E-6, 5.99176451380661E-6, 7.62231692622485E-6,
    7.77791144401441E-6, 7.20798838849436E-6, 7.4076656059198E-6,
    1.18279923551036E-5, 7.02828513612093E-6, 8.17272614856922E-6,
    7.32655336687172E-6, 6.87613229711617E-6, 2.25566355431345E-5,
    3.11701747771784E-5, 3.25566773520444E-5, 3.47120986257364E-5,
    8.3396293416293E-6, 8.05693145627725E-6, 9.68367964356002E-6,
    4.65594685243503E-5, 0.000182392908026845, 0.000190677027887626,
    0.000190990821370872, 8.67959077922145E-5, 1.99536314470984E-5,
    1.27075296371139E-5, 1.1052875151547E-5, 7.7119790837905E-6,
    7.595209194367E-6, 8.44703857012239E-6, 9.74059409556715E-6,
    7.53014433277723E-6, 0.000131054272537454, 0.000194987693089166,
    0.000201942552937496, 0.000202156495307148, 2.20799456629282E-5,
    8.38082016027723E-6, 7.23070875168714E-6, 1.29078167186318E-5,
    1.25475706847326E-5, 1.05738968810019E-5, 1.3444774737508E-5,
    0.000125057380981563, 0.000202751780428698, 0.000208318028814146,
    0.000207304069209772, 1.9035909490035E-5, 1.31003700374231E-5,
    9.03708718716824E-6, 7.21183545513936E-6, 8.10941888855177E-6,
    8.36452098003217E-6, 1.05239934089333E-5, 9.8290133968549E-6,
    7.8702560529379E-6, 2.1351091043954E-5, 0.000142993067246535,
    0.00019250147515911, 0.00019628980255219, 0.00019891257098437,
    0.000244860148805642, 1.60952263346345E-5, 1.93861813614291E-5,
    8.89424251407755E-6, 6.87819059589556E-6, 2.29991638659667E-5,
    6.57412673821209E-6, 9.04701920886069E-6, 7.50475330964605E-6,
    4.4469576223189E-5, 0.000240937469446067, 0.000213617050992454,
    0.000213548321820482, 0.000219592968211454, 8.96649596658214E-5,
    9.12224082960632E-6, 9.01370763182398E-6, 9.78608468278451E-6,
    6.32024286833998E-6, 1.0953175503588E-5, 1.04092965850766E-5,
    0.000136149303320483, 0.000240289993463775, 0.000227711495114367,
    0.000225181883230187, 0.00021628993054232, 6.48646034848174E-5,
    1.07599703767004E-5, 2.48480112866816E-5, 1.15250754501056E-5,
    2.54773221191791E-5, 7.35124346803246E-6, 8.54843425273317E-6,
    6.68355320629884E-6, 7.65327893424043E-6, 4.04795671593529E-5,
    0.00018051533216497, 0.000211513016854531, 0.000224057154114223,
    0.00022995832884821, 0.000236404501775648, 9.45496107437806E-6,
    5.40461443856234E-6, 6.92141727193684E-6, 8.62081445780773E-6,
    8.1076954951324E-6, 9.99901895935313E-6, 3.3799181047316E-5,
    0.000199093906412319, 0.000201931795500936, 0.000203415882056365,
    0.000191182356689207, 0.000214089556321459, 0.000181352912549726,
    0.000177672642946819, 2.0519202804441E-5, 7.52634585510506E-6,
    8.93530832036498E-6, 7.34268083645732E-6, 8.33033824673544E-6,
    8.63746003959862E-6, 1.07266093796647E-5, 0.000123082217527794,
    0.000199401403662162, 0.000206004660027085, 0.00020575782975844,
    3.46578754478116E-5, 5.2487605804462E-5, 1.29214554223584E-5,
    2.91216096895185E-5, 2.23288688447963E-5, 1.11129579666957E-5,
    1.08433709653944E-5, 3.46873214069178E-5, 0.000130227737367415,
    0.000162184396898947, 0.0001649954732425, 0.000168727325785244,
    7.66554257018134E-6, 7.09976651194544E-6, 6.97040828938989E-6,
    1.00038337037274E-5, 9.1507158180621E-6, 1.02399654255842E-5,
    1.00463272867812E-5, 1.19109888853156E-5, 0.000194220612286393,
    0.000204265748628445, 0.000204095964076271, 0.000227587809712162,
    1.11719320857264E-5, 1.1073404143643E-5, 1.6291316171754E-5,
    7.19320522524966E-6, 7.97745370596497E-6, 7.88900171865842E-6,
    1.51544184695889E-5, 0.000164612889585702, 0.000232988304528318,
    0.000220608248417606, 0.000219779521761084, 0.000209034101205281,
    3.55658034344122E-5, 1.65173843506373E-5, 1.45077489301168E-5,
    1.06780421407878E-5, 1.33706640832657E-5, 7.56944798374833E-6,
    6.98635108287269E-6, 6.27098135503542E-6, 0.000126981142318357,
    0.000196434544623138, 0.000205438622424672, 0.000204643075873289,
    0.000172465689225266, 0.000110145586326676, 9.61618552396232E-6,
    8.45858402159569E-6, 8.76194018171591E-6, 9.92406719685199E-6,
    3.32451383720426E-5, 0.000185413151392429, 0.0002120466566674,
    0.000222362397272023, 0.000226244169996065, 0.000213606449357302,
    1.8320786355966E-5, 1.03526848470488E-5, 1.68331546727362E-5,
    6.31495861119419E-6, 4.19888925336015E-5, 0.000203944891138257,
    0.000234834767635248, 0.00023699489772734, 0.000154848210603212,
    0.000146513154898452, 4.33193098770713E-5, 6.63429040193692E-6,
    8.59909699237623E-6, 7.39696323561077E-6, 7.68562745163928E-6,
    8.33182681644324E-6, 5.35817582167207E-6, 3.69365092386796E-5,
    0.000125881608861494, 0.000134884418626897, 0.000134272893266581,
    0.00019714650229541, 0.000214231581492537, 0.000220689526059781,
    7.51238104987704E-5, 9.73449259581761E-6, 1.02161185384595E-5,
    1.64831945598401E-5, 9.34991008302638E-6, 1.01676833986492E-5,
    6.90615455599214E-6, 0.00010039543171366, 0.000220888612812708,
    0.00022241621238541, 0.000220398601831555, 0.00021041499890405,
    0.000107475911032138, 1.93227172309166E-5, 7.84159957541047E-6,
    7.45081471450551E-6, 8.22577915328775E-6, 3.04399037452039E-5,
    7.89520015601316E-6, 0.000130336286753719, 0.000162960210836434,
    0.000165270563454948, 0.000175103020113829, 1.40071470419429E-5,
    1.12112624866002E-5, 9.02588987047006E-6, 5.44200225917512E-6,
    1.86330806652456E-5, 6.02335722814405E-6, 9.85889747986619E-6,
    9.48078439982385E-6, 6.635532576405E-5, 0.000195983109270294,
    0.000203752035344902, 0.000203279307171529, 4.53274078393988E-5,
    1.18758675461237E-5, 8.99279801474159E-6, 9.69786626871734E-6,
    9.84814215147898E-6, 0.000138789692129707, 0.000192623966803457,
    0.000198603550471606, 0.000197621722472302, 2.81929080365041E-5,
    2.92214269847574E-5, 8.40469259650139E-6, 8.8729633127233E-6,
    8.97997313638162E-6, 9.58259973845458E-6, 3.4349308022849E-5,
    0.000180262962078701, 0.000187519854057759, 0.000187235705407889,
    0.000163043711797016, 6.09663119044916E-6, 1.12378910943718E-5,
    9.41498545669321E-6, 1.29236349894246E-5, 9.69928085832781E-6,
    0.000141055920113184, 0.000165486238445531, 0.000175052405105736,
    0.00017793196538435, 0.000171971239793913, 2.0027642096376E-5,
    9.32094059313049E-6, 1.03129640510443E-5, 1.10872435520174E-5,
    1.25716487041965E-5, 0.000126747262082789, 0.000155404296556452,
    0.000208299656134655, 0.000213266697567903, 0.000207596118743351,
    5.02949726179994E-5, 1.11607994185781E-5, 1.16441224499857E-5,
    1.27452731219501E-5, 1.24018951452923E-5, 1.38426041642613E-5,
    7.49663938559334E-6, 0.000154796548056764, 0.000192072751896236,
    0.000184714800196626, 0.000192454131570645, 0.000183744340640712,
    3.10601218557711E-5, 1.01195189974224E-5, 9.72076967036776E-6,
    1.16617519379042E-5, 2.17754419078708E-5, 1.360690426929E-5,
    0.000121722884091643, 0.000177463728388293, 0.000202327040849214,
    0.000200308976226997, 0.000169634408571692, 0.000153044194884499,
    2.65522417860475E-5, 9.71309031572238E-6, 7.9143156958644E-6,
    7.16466538605128E-6, 8.31740320284423E-6, 0.0001160411124222,
    0.000199066789372968, 0.000203863729894216, 0.000202939452958022,
    0.000161811912770319, 0.000176857554268389, 1.67468569592311E-5,
    1.43668184500497E-5, 2.49685913396858E-5, 1.24254523519368E-5,
    5.88323474455508E-6, 3.61850594558556E-5, 0.000192136281901599,
    0.000230746717096186, 0.000227442331782329, 0.000211327177741283,
    0.000172933966685, 2.98922498453373E-5, 9.79959830209659E-6,
    8.02323795117609E-6, 8.17344495921623E-6, 1.09798384161832E-5,
    1.75280555320196E-5, 3.83585527839719E-5, 0.000182864855972276,
    0.000209667774154251, 0.000205595735357034, 0.000172825972520345,
    0.000165086773687226, 2.21424056849083E-5, 8.77299903877754E-6,
    3.17268223458418E-5, 2.46975359778342E-5, 9.77662972577318E-6,
    7.57558236967069E-6, 1.58823812612895E-5, 0.000203699833733067,
    0.000186954920851678, 0.0001844960848329, 0.000222854169355616,
    0.000178078038537919, 2.61142835407292E-5, 3.35910299462911E-5,
    9.12382004157381E-6, 1.17131384786661E-5, 3.47708992107623E-5,
    8.14927470269412E-6, 3.95886499246976E-5, 0.000184969420336649,
    0.000202522425874849, 0.000199631235010628, 0.000184159127782214,
    0.000186258027621732, 1.84376170416272E-5, 1.35544797876524E-5,
    1.29812585651548E-5, 1.75614136108169E-5, 4.67112858805875E-6,
    4.29204337766929E-6, 1.66064021257502E-5, 0.000164002915008637,
    0.000205227376974565, 0.000193713500047691, 0.000221552759551078,
    0.000150383372543337, 1.93546434676282E-5, 8.65351229169799E-6,
    5.05651788716869E-6, 4.92492701710892E-6, 4.8239164782765E-6,
    5.43777912741768E-6, 1.33147861102081E-5, 0.000222567398339372,
    1.01460844498245E-5, 1.01072337658764E-5, 1.33100372659702E-5,
    7.60735643036845E-6, 8.76310197686004E-6, 8.58766798133775E-6,
    7.81375005111288E-6, 8.05347323587743E-6, 8.67921239510477E-6,
    7.53483328827971E-6, 8.34632617504363E-6, 7.8102464074801E-6,
    9.21754177614702E-6, 1.59496673507755E-5, 1.46330576247228E-5,
    1.42056652630232E-5, 9.16810843922698E-6, 6.90561114590368E-6,
    7.15990610681954E-6, 5.98251867871766E-6, 7.31205377286595E-6,
    1.18901342450306E-5, 8.22157109933577E-6, 2.79178959307808E-5,
    2.59643094529404E-5, 2.54867954929405E-5, 2.62499595630359E-5,
    1.6320937964143E-5, 0.000177124213181553, 0.000183166538157828,
    0.000186343027176331, 0.000201854921404998, 4.79927948549815E-5,
    1.05181380201246E-5, 1.06831031827466E-5, 7.17191622775398E-6,
    1.69843856869611E-5, 7.20305990148614E-6, 7.2501541822744E-6,
    8.43554950387992E-6, 1.2141027221209E-5, 2.52458539968755E-5,
    9.2506059325947E-6, 1.98266015355231E-5, 2.46301593498925E-5,
    1.58994081023E-5, 7.53441517082579E-6, 3.32061444557485E-5,
    0.00016044729435023, 0.000178509652150884, 0.000180038968614073,
    0.000195144126440294, 2.78233943402766E-5, 2.33747848484418E-5,
    1.20345804572114E-5, 2.90515506593109E-5, 1.24092874850007E-5,
    8.95745537762341E-6, 8.72029619224173E-6, 0.000160665316426585,
    0.000202294278761182, 0.000203251916804623, 0.000205414593739316,
    1.54626445462589E-5, 1.04433321026038E-5, 8.30185471386404E-6,
    7.85536090292521E-6, 8.29386579156038E-6, 6.97068293950418E-6,
    8.39933382352594E-6, 2.70928625744859E-5, 0.000120317588108555,
    0.000157553721737705, 0.000158855112970551, 0.00015890904125777,
    2.1490858007152E-5, 8.01712357669074E-6, 9.17223508754399E-6,
    1.17290122324033E-5, 1.17983520436549E-5, 2.54977373504287E-5,
    7.83237900509189E-6, 1.01907289927933E-5, 0.000193841983030029,
    0.000207561910329852, 0.000200375288509805, 0.000218049990159131,
    0.000191882567908045, 1.96680216995463E-5, 5.94043382199198E-6,
    9.74984212525311E-6, 9.97510197431104E-6, 3.1195958792736E-5,
    6.61280690405296E-6, 1.55903855564446E-5, 0.000177666508965983,
    0.000211078673258689, 0.000221366808410471, 0.00022036788762769,
    0.000220030582077982, 7.22785699491629E-6, 1.45890424410825E-5,
    5.78642474870604E-6, 1.1324044234211E-5, 1.07707212816266E-5,
    1.89764805402985E-5, 6.96638709965909E-5, 0.000187887483478018,
    0.000208436782820388, 0.000209932915855354, 0.000229229080035883,
    6.05300493036053E-5, 7.86663035227478E-6, 2.59690893920873E-5,
    8.44919043280399E-6, 1.21673830027306E-5, 0.000162299163248839,
    0.000212507803842311, 0.000218801647837539, 0.000222256273471257,
    0.000218463951839776, 1.21033106898664E-5, 3.63751540851115E-5,
    1.40117043001233E-5, 1.78530684629984E-5, 1.61614078073428E-5,
    0.000143933294354087, 0.000192776579935712, 0.000195373680160098,
    0.000195254504124328, 8.85326004420613E-6, 9.24289078962152E-6,
    8.97381106656166E-6, 2.41491787220401E-5, 0.000203731278054671,
    0.000209809040271665, 0.000209258260943535, 0.000203474251807436,
    8.81928224418335E-6, 1.4490282820835E-5, 1.22885155862071E-5,
    0.00017668354441379, 0.000183890289482804, 0.000184830037867932,
    0.00019291859023317, 1.08621684175062E-5, 1.2443763085315E-5,
    6.15255128826237E-6, 1.30279431515782E-5, 3.16428072735162E-5,
    3.63065070531686E-5, 0.000206851040715549, 0.000209368163930491,
    0.000207237768254406, 0.000189466387045622, 0.000170220955034366,
    9.02878855482882E-6, 7.10391626415519E-6, 1.429289862135E-5,
    1.35877765235687E-5, 3.32196291777519E-5, 0.000181932721539842,
    0.000206791623153745, 0.000218305165477053, 0.000226427752117102,
    0.00023971675286429, 2.71229719197448E-5, 1.63187796840985E-5,
    9.66569166882634E-6, 1.26719215583954E-5, 1.26267056514271E-5,
    1.34206651807906E-5, 0.000155479453138653, 0.000202350514930538,
    0.000223483041433876, 0.000225120525038446, 0.000220405994737832,
    1.10224169552461E-5, 1.06245359111933E-5, 1.14078343930062E-5,
    1.13018978635804E-5, 8.03952033761384E-6, 0.000196473935796256,
    0.000261999938822453, 0.000242781708397672, 0.000244857697547946,
    0.000206647590202368, 2.41707406055077E-5, 1.34845217063989E-5,
    1.4746315097426E-5, 1.59082798005799E-5, 0.000184807169300158,
    0.000262584395908765, 0.00023937971800293, 0.000237387393390147,
    0.000182452975960247, 1.45402949664645E-5, 1.27927908880284E-5,
    1.14460627591939E-5, 1.85580714336511E-5, 0.000175568774902735,
    0.000213568326876014, 0.000225304765495157, 0.000237161681814445,
    0.000231855023811663, 7.57861048392616E-5, 8.10561956959648E-6,
    7.03147604693226E-6, 7.87841486510365E-6, 9.14558804903927E-6,
    7.4797797628481E-6, 2.03165780664564E-5, 6.72907998077055E-6,
    0.000133478642044135, 0.00020000413091289, 0.000204998344482471,
    0.00020512430221262, 2.39088709675666E-5, 2.54596557500542E-5,
    2.85641406256451E-5, 1.32436050220331E-5, 0.000108251616162259,
    0.000200999312279417, 0.000205988782789052, 0.000204573981208438,
    1.6914428870612E-5, 1.82089830966955E-5, 2.78892339078611E-5,
    5.69358945400568E-6, 1.06112992317037E-5, 8.62735432506446E-6,
    9.85863800091246E-6, 7.21636309759946E-6, 1.73957805618332E-5,
    2.729882758574E-5, 2.59932654110132E-5, 2.55470111195723E-5,
    1.75991495136028E-5, 1.85776169412593E-5, 9.83336379927673E-6,
    1.24194440574556E-5, 9.98045079028148E-6, 1.06456730678079E-5,
    1.24484198437086E-5, 1.23661635894197E-5, 1.21699587850492E-5,
    1.29048617671093E-5, 2.39257887072541E-5, 2.83726955302265E-5,
    2.7730742994153E-5, 4.12791094878234E-5, 8.32672060527911E-5,
    8.59617325130236E-5, 8.81187316947805E-5, 1.15552106087484E-5,
    1.1521244150066E-5, 0.000106929648539941, 0.000220662573455911,
    0.000207899500754917, 0.000215711891876928, 0.000241820383667611,
    5.04292066535223E-5, 9.77944636078748E-6, 1.03103383377605E-5,
    1.03117160982628E-5, 8.06638968793809E-6, 7.35528781740395E-6,
    7.7035964929166E-6, 1.03888862031045E-5, 5.00235245313151E-6,
    4.75533056122428E-5, 0.000232093871793515, 0.000220041087702034,
    0.000229417514269565, 0.000254248462658143, 8.29811319358672E-5,
    7.03555857840813E-6, 7.43007833115385E-6, 7.94038165977698E-6,
    7.59121359418337E-6, 7.48263221576611E-6, 8.72531274217694E-6,
    1.0694557043587E-5, 8.61714928894396E-6, 2.44785132830888E-5,
    0.000166403521948993, 0.000173778158173859, 0.00018033681347421,
    0.000195926534766487, 0.000154720071863859, 1.07555450722517E-5,
    1.05079367474144E-5, 1.17923090472505E-5, 9.83414854687996E-6,
    1.0601669694311E-5, 7.71161459876915E-6, 8.38884960490581E-6,
    5.9032276550995E-6, 0.00014632711444459, 0.000196312974735306,
    0.000200142405061184, 0.00020015242602187, 1.0816970580286E-5,
    8.40595495050319E-6, 5.66412679189899E-6, 8.95874783112723E-6,
    1.2056555538032E-5, 8.46891507236601E-5, 0.000201691895715619,
    0.000207083200518918, 0.000207365338195632, 3.48498367198272E-5,
    1.06426209250593E-5, 1.16228357522886E-5, 7.45451630415748E-6,
    1.25294229610682E-5, -6.02556852690611E-6, -7.08546275363364E-6,
    -7.92407040942471E-6, -7.92869260821444E-6, -1.01890272798788E-5,
    -7.90329866269586E-6, -7.16386785365877E-6, -4.81646931921721E-6,
    -6.27412655764235E-6, -6.32032161856858E-6, -6.4277263730826E-6,
    -7.1997134574812E-6, -7.09983221044316E-6, -7.14535501001467E-6,
    -9.2178301524846E-6, -8.89724885893466E-6, -9.31535965079218E-6,
    -4.95109618594073E-6, -4.62418809355966E-6, -9.09266796072458E-6,
    -2.30114128975628E-5, -1.57357960340977E-5, -1.55638303036247E-5,
    -1.38026113468313E-5, -2.92893026450334E-5, -7.19365183100657E-6,
    -7.34748267764355E-6, -6.1532349252882E-5, -0.000105983249276932,
    -0.000100781661270265, -9.94967554025955E-5, -7.26402335110431E-5,
    -1.17904315974394E-5, -1.26708378683127E-5, -1.17643416486314E-5,
    -8.86435437090101E-6, -7.12792624940644E-6, -5.90948057827549E-6,
    -6.08659547832763E-6, -6.43217258211434E-6, -0.000108831677439631,
    -0.000123077022091815, -0.000115662140419308, -0.000123363073659752,
    -0.000100833332819802, -4.04011258142228E-5, -3.3444828990155E-5,
    -1.2983946339104E-5, -2.28851974351959E-5, -1.22964861121435E-5,
    -8.79972212172406E-6, -0.000139227949085837, -0.000111675811395968,
    -0.000106952765925552, -0.000103035852145246, -6.70219002494882E-5,
    -2.40401373810073E-5, -8.82973061066254E-6, -9.16552974758949E-6,
    -9.2630468968964E-6, -1.60445270750024E-5, -1.30747280426073E-5,
    -8.60611907134829E-6, -8.2098238782417E-6, -1.93129826556684E-5,
    -0.000132092711252589, -0.000247334783401648, -0.000243759759206527,
    -0.000253728472608469, -8.47286935260469E-5, -1.21493332373809E-5,
    -1.22906176849695E-5, -1.90105511404098E-5, -7.4426836914305E-6,
    -1.62265679382215E-5, -7.2140585955139E-6, -1.19065487016938E-5,
    -1.29315891661446E-5, -7.53630481593429E-5, -0.00011906489878459,
    -0.000199774262221528, -0.000192483777059585, -0.000184010038098542,
    -6.52282964808453E-5, -9.24338472684946E-6, -7.1630893295294E-6,
    -1.13729891604613E-5, -1.46623490589868E-5, -8.50197973418695E-6,
    -1.4834225942083E-5, -0.000143169763561452, -0.000189150329435444,
    -0.00021756907935317, -0.00021293695349822, -0.000171743715912411,
    -0.000144815987376321, -1.78837593238031E-5, -1.37203780359655E-5,
    -7.83892545719453E-6, -1.01891643058264E-5, -9.64273831177852E-6,
    -1.91989624673716E-5, -9.13927209228733E-6, -8.73639450599239E-6,
    -6.23066133932774E-5, -9.46101341531574E-5, -0.000151644491798601,
    -0.000151247718824871, -0.000131573449297755, -0.000107443667066393,
    -1.15802409072821E-5, -9.71145384366107E-6, -1.19853858019107E-5,
    -9.7649760459759E-6, -8.29554802420788E-6, -1.04392390685457E-5,
    -5.77459802216847E-5, -0.000104927678927213, -0.000142282813310025,
    -0.000137838536995665, -0.000137814921921353, -8.29030792027732E-5,
    -9.03189230538363E-5, -7.69927025308173E-5, -4.70787371695573E-5,
    -6.71250973552667E-6, -7.83175056216015E-6, -2.06013389628086E-5,
    -5.54028584187758E-6, -8.98015158479247E-6, -1.87782538577182E-5,
    -0.000151786763051843, -0.000153430342636528, -0.000154876189612503,
    -0.000148180979320421, -0.000118847406126866, -1.55389333698992E-5,
    -3.07615899351031E-5, -1.16492666134645E-5, -7.55341808518132E-6,
    -7.36086707032974E-6, -7.8845242604491E-6, -5.02675913295017E-5,
    -0.000147344372998698, -0.000148756968703973, -0.000147800596329585,
    -0.000107693493126941, -2.26255838679738E-5, -3.19447859461986E-5,
    -6.63541192742118E-6, -8.21232994039529E-6, -1.43520430182336E-5,
    -1.10327010694137E-5, -1.3439271043007E-5, -1.63347003309234E-5,
    -0.000124293259600942, -0.000121617188856484, -0.000122671976177924,
    -7.41356544528088E-5, -1.58036830457477E-5, -8.05554195661706E-6,
    -7.84910186643456E-6, -9.4889946875425E-6, -9.27698668451761E-6,
    -9.4614706826192E-6, -1.09433362437876E-5, -0.000117917873591544,
    -0.000162711720650718, -0.000158974795087298, -0.000154785361562016,
    -0.000127005972596361, -6.57942591661849E-5, -9.44048294248277E-6,
    -8.79271581636417E-6, -6.54418249287983E-6, -2.94577178926006E-5,
    -1.23147357099566E-5, -1.17144521031596E-5, -7.64484749097131E-6,
    -0.000107553228915546, -7.97985295084535E-5, -0.000128918855010882,
    -0.000124789368442979, -0.00010675846752747, -7.2875370238731E-5,
    -1.043938271986E-5, -1.02781562532686E-5, -7.01038957873558E-6,
    -6.36953730991268E-6, -6.10623686878162E-5, -0.000105528219178754,
    -0.000134417692923569, -0.000124041018257892, -0.000108381428421387,
    -0.000109354132191514, -1.56279044876087E-5, -7.83508507943124E-6,
    -8.7571432301877E-6, -7.66435682931902E-6, -6.58242670875512E-5,
    -0.000118622810972416, -0.000149425442494876, -0.000177372757942621,
    -0.000164086473071829, -0.000160904897743612, -0.000125670981253556,
    -5.56765473331136E-6, -8.57577211605969E-6, -5.01405648571673E-6,
    -9.17712149147302E-6, -1.07694642448434E-5, -4.4793105219457E-6,
    -5.5553503207889E-5, -0.000132126625772427, -0.000133109500938134,
    -0.000132648981662207, -0.00017132958104495, -0.000166272152489918,
    -0.000177049399661493, -0.000112106206095441, -1.31483781313437E-5,
    -1.50917616103854E-5, -1.04862398841128E-5, -1.1077265076536E-5,
    -7.4695601516355E-6, -8.16848558335354E-6, -0.000104456591616292,
    -0.000115453770172418, -0.000193558471103268, -0.000189427398046921,
    -0.000168862605066673, -5.82667309775874E-5, -1.47263519202479E-5,
    -9.8648216705797E-6, -1.65789333226946E-5, -1.00896163363329E-5,
    -1.41299472946306E-5, -2.17062926416481E-5, -0.00013717847886131,
    -0.000140956442143946, -0.000140539994104315, -8.92606958492249E-5,
    -2.68471257577439E-5, -1.26423251296378E-5, -1.04673663041728E-5,
    -1.64511162784549E-5, -6.34478854749066E-6, -7.25327743672786E-6,
    -6.47602547900727E-6, -8.39251696927845E-6, -8.63666392123707E-5,
    -0.000112616579194632, -0.000108691204450629, -0.000108302541167977,
    -8.46814493183758E-5, -7.06645230402581E-6, -9.42216328772719E-6,
    -1.03670353139644E-5, -1.08283646154916E-5, -0.000136658514052204,
    -0.000103788300123494, -0.000101407250117456, -9.81445698013978E-5,
    -4.24968183825116E-5, -1.01935554113222E-5, -9.72958494031206E-6,
    -1.67604718361732E-5, -7.77658413366376E-6, -9.9668034685757E-6,
    -5.23005011851107E-5, -0.000105131084846968, -9.95786366746511E-5,
    -9.93414680603111E-5, -5.10697440092769E-5, -1.91638227262249E-5,
    -2.17068522208394E-5, -9.71604603348211E-6, -1.00985994475843E-5,
    -1.06509817119488E-5, -7.3825162140996E-5, -0.000125765100958807,
    -0.000138528638382313, -0.000144913217887655, -9.65286897053719E-5,
    -5.77816402785356E-5, -1.08206980202327E-5, -9.30314741995128E-6,
    -8.70390561855816E-6, -9.21097707212864E-6, -8.21181728630324E-5,
    -8.15786787294555E-5, -0.000129104577960112, -0.000123290046328727,
    -0.000122191366910927, -8.01332659193145E-5, -2.53646820753033E-5,
    -1.1445160823266E-5, -9.64662360824063E-6, -1.67765276569339E-5,
    -1.09073625522794E-5, -6.19147861336273E-6, -8.84824331764742E-5,
    -0.000142663313661106, -0.000138012337854178, -0.000136636552759179,
    -0.000122880387882389, -9.23051066941433E-5, -8.03204326064869E-6,
    -1.87448465700877E-5, -1.56246728031957E-5, -1.45232560711771E-5,
    -2.03823973408069E-5, -0.000120087897855194, -9.47754096112311E-5,
    -0.00016808578850942, -0.000167598653588042, -0.000156443316063371,
    -9.79150569928182E-5, -1.44706762133445E-5, -8.02723331614643E-6,
    -1.09763393851365E-5, -1.01283887750177E-5, -1.12448925668894E-5,
    -0.000114309286033622, -0.000128413456498907, -0.000175478291436554,
    -0.000172805419516137, -0.000153208728724503, -8.98891438149597E-5,
    -1.39364298278648E-5, -1.44674936371662E-5, -7.77162590564877E-6,
    -7.03255030744902E-6, -4.67062267402925E-6, -5.35566542758341E-5,
    -0.000125079301335257, -0.000147999731388913, -0.000170605997361751,
    -0.000166728028397105, -0.000139895030996547, -7.03919350070539E-5,
    -1.44186956268206E-5, -1.80400448143068E-5, -4.95958386279141E-6,
    -8.88237809913433E-6, -9.54670224674149E-6, -6.48277960672085E-5,
    -9.9775959679868E-5, -0.000137003617050119, -0.000148505563913953,
    -0.000128802627273635, -8.06098745996011E-5, -3.40551144303451E-5,
    -2.74031514911324E-5, -1.05312570004273E-5, -1.02260896640896E-5,
    -9.38721840943605E-6, -8.31378976137132E-6, -2.09432048958447E-5,
    -0.000109545178646695, -0.00020077326567707, -0.000204053787433355,
    -0.000212198471985957, -7.17771608227149E-5, -5.74772676717559E-5,
    -1.58941811421365E-5, -3.19804037942891E-5, -1.61802137906234E-5,
    -1.18715709506737E-5, -8.66607655486836E-6, -6.59414506133157E-5,
    -0.00013917065742757, -0.000171910068627247, -0.000171746421263465,
    -0.000160989938347199, -0.000181462843762779, -3.36157725920385E-5,
    -3.90188717454086E-5, -8.60526682680336E-6, -5.17818790050863E-6,
    -4.75967657331188E-6, -7.13864955359335E-6, -2.07493219819882E-5,
    -9.17173380566854E-5, -0.000167491863741817, -0.000170002014801742,
    -0.000174590459012953, -0.000135597652108803, -2.99060248271595E-5,
    -7.12895062936994E-6, -9.7414093547956E-6, -6.20946032432039E-6,
    -5.64644886332094E-6, -8.51766411736862E-6, -1.25351570346014E-5,
    -0.000117877290832724, -9.00962769799962E-6, -8.90828149768766E-6,
    -1.09754113732237E-5, -1.86534602399242E-5, -7.81772807573526E-6,
    -7.84465385743436E-6, -1.743887678613E-5, -7.54111620061341E-6,
    -7.57006425388482E-6, -2.0757123807056E-5, -7.46610257666996E-6,
    -6.06448111715701E-6, -8.56373667446394E-6, -6.93571313349603E-6,
    -8.70609481639531E-6, -8.65204835073109E-6, -6.9571136576241E-6,
    -7.65161338875964E-6, -8.03836208137425E-6, -6.90134988109705E-6,
    -6.91558729105932E-6, -7.55198619766136E-6, -7.19730363505836E-6,
    -2.53297423427984E-5, -1.64192126847277E-5, -1.77279608324483E-5,
    -1.67248144117318E-5, -3.60238458326389E-5, -0.000120508648769839,
    -0.000113353479678742, -0.000111891480693582, -8.9923096610291E-5,
    -1.09570301543002E-5, -5.79959605580285E-5, -4.12784423709585E-5,
    -1.10047195508233E-5, -5.03188380082114E-6, -2.06206513992994E-5,
    -6.71859600183061E-6, -2.9588904246708E-5, -1.06723774723492E-5,
    -7.76767805464598E-6, -9.36390394898337E-6, -7.78649836479191E-6,
    -7.7023488499523E-6, -6.90364591721322E-6, -6.14042806605345E-6,
    -5.00283285978798E-5, -0.00010760547835701, -0.000101695005976776,
    -0.000100732555256275, -9.62323838134865E-5, -1.44220321209722E-5,
    -9.97884441345415E-6, -8.83084946856155E-6, -1.22297141952209E-5,
    -1.27164368058032E-5, -1.92922529296537E-5, -1.13118126138021E-5,
    -0.000107955600148724, -9.89796647717639E-5, -0.000101159738309494,
    -8.62029822803755E-5, -4.97293533904436E-5, -8.42884395625992E-6,
    -1.26604673988516E-5, -7.44068100137793E-6, -1.01589750886908E-5,
    -8.77000773137561E-6, -1.33677241925534E-5, -1.02835793517134E-5,
    -7.7783042954993E-5, -6.71661701292855E-5, -6.47346754865072E-5,
    -6.89665582960629E-5, -1.37722664684918E-5, -9.59287188227208E-6,
    -9.51176140636695E-6, -9.26710157312329E-6, -8.7984284075384E-6,
    -1.10317701027363E-5, -2.17813442918134E-5, -9.48889880899632E-6,
    -0.000102594780300151, -0.000202092948791066, -0.000199641219281779,
    -0.000208045019527247, -0.000134156877707748, -8.87402029863197E-6,
    -5.53759444364415E-6, -9.9395350544737E-6, -9.45345607902287E-6,
    -8.04436398829416E-6, -7.34998608037067E-6, -1.88798516308301E-5,
    -9.17699377324073E-5, -0.00013369915768439, -0.000128065821357687,
    -0.000133698123893735, -0.00010725400649588, -1.02789912009329E-5,
    -7.02668359183318E-6, -1.18683719166875E-5, -9.19865745646827E-6,
    -9.73056881734642E-6, -7.83751537821542E-6, -8.78789020851384E-5,
    -8.28358246583593E-5, -0.000103706778992302, -0.000101837954572223,
    -0.000155005410576701, -0.000123815104614886, -1.34191176382055E-5,
    -1.32029164118823E-5, -1.03687389437742E-5, -9.59880126627542E-6,
    -8.35264868693903E-5, -0.000173625056360981, -0.00014670920608473,
    -0.000145425104717722, -0.000120615805202605, -2.50003431326831E-5,
    -1.21467846701004E-5, -1.30582469927957E-5, -1.47937307709011E-5,
    -9.12612049606345E-6, -0.000101425611422617, -8.02424841902751E-5,
    -7.83285740730998E-5, -8.1938831816519E-5, -1.06598926922069E-5,
    -1.78595441618198E-5, -1.27944940340356E-5, -3.9929000624151E-5,
    -9.5709309775219E-5, -8.84450582457204E-5, -8.82070953383576E-5,
    -8.55830446141189E-5, -8.53579100471376E-6, -9.35749189340938E-6,
    -8.24463140246975E-6, -7.17541472327732E-5, -6.7124758830406E-5,
    -6.73953821925198E-5, -6.25857476634165E-5, -1.06777342511327E-5,
    -1.0026233472256E-5, -1.00614823280587E-5, -1.60176533321316E-5,
    -1.40072392819721E-5, -5.58975303336541E-5, -8.63626730276465E-5,
    -0.000173501040489822, -0.000167267078599289, -0.000150554387524124,
    -0.000125472093715376, -1.03395472181311E-5, -6.92465826560909E-6,
    -1.48166057361839E-5, -1.61937189328531E-5, -5.74229597692279E-5,
    -7.54797439566358E-5, -0.000168107777783916, -0.000161856547277654,
    -0.00014800785089406, -0.00012054079894264, -1.92611446576112E-5,
    -1.23610285355771E-5, -1.13377534254899E-5, -1.18979220989464E-5,
    -1.08051341957663E-5, -1.21362991639651E-5, -7.71424976131795E-5,
    -0.000124310053881357, -0.000118665288512468, -0.000122136773344633,
    -9.39746938811192E-5, -1.22500399552478E-5, -8.17862926965515E-6,
    -1.39454827054923E-5, -1.22786952802016E-5, -1.246412654168E-5,
    -0.00011697356590065, -0.000196554962772372, -0.000182993897226551,
    -0.000183039990582428, -0.000136628964872889, -6.6824151101677E-5,
    -1.98896420549044E-5, -1.0826469717835E-5, -1.88739928774254E-5,
    -0.000128339004639744, -0.000265418759863629, -0.000267014980556481,
    -0.00026618134095067, -0.00016234065264131, -9.94203443273271E-5,
    -1.25082484003978E-5, -2.74996973019716E-5, -1.10053888162049E-5,
    -9.78138843589459E-5, -9.74094809632773E-5, -0.000145972435103968,
    -0.000140577990870888, -0.000139600379512852, -0.00010706073651693,
    -6.78747339078645E-6, -4.09913638513127E-5, -5.3042988100509E-6,
    -4.67604473097306E-6, -9.28458706837542E-6, -9.28602685533791E-6,
    -8.88481008381012E-6, -0.000129636554571217, -9.05421773508301E-5,
    -8.75212077740225E-5, -8.67014655404683E-5, -2.21641702996357E-5,
    -2.08266935821295E-5, -3.95741430044587E-5, -1.291062026496E-5,
    -0.000120057117960042, -9.43689589861217E-5, -9.08042440945085E-5,
    -8.68427113697677E-5, -0.000106585553280915, -1.95263359502321E-5,
    -8.42984081100369E-6, -7.93060149905834E-6, -8.97873706364797E-6,
    -1.45454898019641E-5, -7.79806009757867E-6, -8.83186264143596E-6,
    -1.61768511578471E-5, -1.83845431443332E-5, -1.7524673564865E-5,
    -1.98864009543984E-5, -1.49265231918684E-5, -2.16410465102776E-5,
    -1.04600963960715E-5, -2.6343766042657E-5, -1.4849394610138E-5,
    -5.78797673861621E-6, -6.19490575861097E-6, -6.62156729191772E-6,
    -1.11163628566891E-5, -1.28061524113779E-5, -1.47998967635133E-5,
    -2.62193074830547E-5, -2.76312391270592E-5, -3.20769744951696E-5,
    -3.10975501730575E-5, -3.25953479707073E-5, -3.65785582595031E-5,
    -4.31484883193128E-5, -1.18913354319119E-5, -0.000111434514757185,
    -0.000229255345562047, -0.000232568001496741, -0.000229848533189152,
    -0.000244241911906071, -0.00013058291097234, -1.80015475361894E-5,
    -9.51847284777669E-6, -9.04853811240937E-6, -1.67706571093064E-5,
    -8.28952309168169E-6, -6.54788029810853E-6, -1.08953647961458E-5,
    -6.813304261215E-6, -7.35729849884408E-5, -0.000163040622565436,
    -0.000220026205095511, -0.000216857878960276, -0.000212622518805025,
    -9.81061299938716E-5, -6.51316667724344E-6, -2.10463006742076E-5,
    -7.36784134220351E-6, -7.95297994715902E-6, -7.9689161489853E-6,
    -5.11832642471466E-6, -5.91493731755451E-6, -6.01909070581253E-6,
    -4.81976934165497E-5, -7.10181224052238E-5, -0.00013665539585833,
    -0.000132174516780171, -0.00011425743972447, -6.30219479969292E-5,
    -7.86566035936102E-6, -7.63884774927333E-6, -8.59729558614299E-6,
    -8.70452372145372E-6, -9.14068361003862E-6, -1.44164739143129E-5,
    -1.57433266453625E-5, -7.24426225636412E-6, -0.000110723989720744,
    -9.67733545683003E-5, -9.48919603491559E-5, -9.93794113423063E-5,
    -2.68420494295977E-5, -8.48694281868123E-6, -1.52411213949675E-5,
    -9.25580363265509E-6, -7.41400582177384E-6, -9.54389500940979E-5,
    -0.000102463351546602, -9.79993698824363E-5, -9.85971755874148E-5,
    -9.48335806898746E-5, -1.12847202377968E-5, -9.46067708850015E-6,
    -1.26473916823556E-5, -1.23822452350256E-5, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0,
    23.0, 25.0, 25.0, 22.0, 0.0, 0.0, 44.0, 63.0, 76.0, 91.0, 35.0, 6.0, 8.0,
    6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 39.0, 90.0, 92.0, 63.0, 49.0, 31.0, 29.0, 9.0,
    15.0, 9.0, 0.0, 43.0, 76.0, 83.0, 72.0, 38.0, 27.0, 0.0, 0.0, 0.0, 14.0, 9.0,
    0.0, 0.0, 20.0, 47.0, 67.0, 77.0, 56.0, 59.0, 31.0, 9.0, 13.0, 0.0, 17.0,
    0.0, 7.0, 11.0, 39.0, 76.0, 86.0, 88.0, 53.0, 40.0, 0.0, 0.0, 6.0, 14.0, 0.0,
    13.0, 39.0, 57.0, 95.0, 73.0, 52.0, 35.0, 24.0, 10.0, 0.0, 2.0, 0.0, 24.0,
    0.0, 0.0, 22.0, 78.0, 94.0, 103.0, 90.0, 44.0, 9.0, 0.0, 8.0, 0.0, 0.0, 6.0,
    36.0, 74.0, 99.0, 97.0, 94.0, 107.0, 82.0, 61.0, 37.0, 0.0, 0.0, 21.0, 0.0,
    0.0, 16.0, 39.0, 64.0, 89.0, 73.0, 58.0, 23.0, 22.0, 9.0, 0.0, 0.0, 0.0,
    19.0, 41.0, 62.0, 64.0, 53.0, 33.0, 19.0, 0.0, 0.0, 20.0, 5.0, 10.0, 12.0,
    55.0, 84.0, 84.0, 57.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 44.0, 72.0,
    107.0, 85.0, 90.0, 57.0, 0.0, 0.0, 0.0, 32.0, 8.0, 6.0, 0.0, 37.0, 89.0,
    97.0, 98.0, 79.0, 54.0, 5.0, 2.0, 0.0, 0.0, 42.0, 69.0, 102.0, 136.0, 111.0,
    68.0, 26.0, 0.0, 0.0, 0.0, 56.0, 65.0, 71.0, 93.0, 73.0, 44.0, 32.0, 0.0,
    0.0, 0.0, 0.0, 7.0, 0.0, 40.0, 44.0, 63.0, 82.0, 81.0, 57.0, 61.0, 38.0,
    13.0, 18.0, 5.0, 7.0, 0.0, 0.0, 35.0, 66.0, 89.0, 110.0, 71.0, 42.0, 15.0,
    0.0, 18.0, 2.0, 13.0, 15.0, 38.0, 85.0, 68.0, 64.0, 36.0, 8.0, 7.0, 14.0,
    0.0, 0.0, 0.0, 0.0, 34.0, 81.0, 89.0, 89.0, 41.0, 0.0, 0.0, 4.0, 5.0, 27.0,
    82.0, 93.0, 95.0, 31.0, 6.0, 0.0, 13.0, 0.0, 0.0, 35.0, 69.0, 88.0, 89.0,
    48.0, 16.0, 14.0, 0.0, 2.0, 5.0, 61.0, 67.0, 93.0, 66.0, 61.0, 33.0, 8.0,
    0.0, 0.0, 0.0, 63.0, 90.0, 135.0, 111.0, 86.0, 40.0, 19.0, 9.0, 0.0, 15.0,
    5.0, 0.0, 58.0, 76.0, 113.0, 94.0, 75.0, 45.0, 0.0, 22.0, 27.0, 14.0, 30.0,
    45.0, 90.0, 93.0, 107.0, 72.0, 51.0, 11.0, 0.0, 5.0, 4.0, 7.0, 29.0, 84.0,
    91.0, 91.0, 79.0, 63.0, 28.0, 17.0, 0.0, 0.0, 0.0, 39.0, 58.0, 66.0, 93.0,
    80.0, 64.0, 56.0, 16.0, 20.0, 0.0, 0.0, 0.0, 53.0, 73.0, 91.0, 104.0, 92.0,
    63.0, 34.0, 14.0, 4.0, 3.0, 0.0, 0.0, 19.0, 65.0, 68.0, 95.0, 81.0, 80.0,
    58.0, 16.0, 23.0, 21.0, 9.0, 0.0, 50.0, 57.0, 84.0, 108.0, 89.0, 70.0, 49.0,
    21.0, 0.0, 0.0, 0.0, 0.0, 14.0, 73.0, 69.0, 101.0, 72.0, 73.0, 53.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 11.0, 65.0, 0.0, 0.0, 10.0, 18.0, 0.0, 0.0, 15.0, 0.0,
    0.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    26.0, 30.0, 38.0, 36.0, 35.0, 90.0, 74.0, 94.0, 41.0, 5.0, 26.0, 29.0, 9.0,
    0.0, 15.0, 0.0, 20.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 43.0, 75.0, 76.0,
    98.0, 48.0, 19.0, 0.0, 0.0, 9.0, 7.0, 21.0, 10.0, 46.0, 110.0, 94.0, 73.0,
    32.0, 0.0, 8.0, 0.0, 2.0, 0.0, 13.0, 5.0, 75.0, 88.0, 101.0, 83.0, 22.0, 0.0,
    0.0, 0.0, 0.0, 8.0, 18.0, 0.0, 64.0, 64.0, 79.0, 51.0, 35.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 41.0, 68.0, 68.0, 89.0, 56.0, 50.0, 3.0, 0.0, 7.0, 0.0, 0.0,
    0.0, 30.0, 93.0, 106.0, 114.0, 71.0, 22.0, 14.0, 10.0, 4.0, 0.0, 66.0, 57.0,
    90.0, 85.0, 62.0, 45.0, 26.0, 21.0, 13.0, 0.0, 56.0, 80.0, 83.0, 67.0, 4.0,
    26.0, 8.0, 24.0, 77.0, 85.0, 85.0, 42.0, 0.0, 0.0, 0.0, 94.0, 88.0, 99.0,
    56.0, 10.0, 1.0, 2.0, 29.0, 16.0, 56.0, 82.0, 62.0, 77.0, 45.0, 33.0, 5.0,
    0.0, 25.0, 15.0, 42.0, 96.0, 81.0, 98.0, 64.0, 50.0, 22.0, 12.0, 8.0, 16.0,
    6.0, 7.0, 67.0, 93.0, 135.0, 74.0, 61.0, 15.0, 0.0, 11.0, 9.0, 9.0, 45.0,
    55.0, 91.0, 74.0, 86.0, 60.0, 21.0, 8.0, 25.0, 52.0, 56.0, 74.0, 60.0, 32.0,
    32.0, 10.0, 23.0, 6.0, 83.0, 99.0, 101.0, 93.0, 72.0, 36.0, 0.0, 27.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 25.0, 81.0, 109.0, 110.0, 72.0, 43.0, 30.0, 28.0, 41.0,
    72.0, 78.0, 81.0, 28.0, 19.0, 0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 29.0, 40.0,
    42.0, 27.0, 21.0, 28.0, 5.0, 17.0, 16.0, 0.0, 0.0, 0.0, 7.0, 21.0, 40.0,
    56.0, 47.0, 66.0, 85.0, 88.0, 76.0, 37.0, 15.0, 30.0, 46.0, 82.0, 86.0, 56.0,
    32.0, 12.0, 0.0, 0.0, 15.0, 0.0, 0.0, 7.0, 0.0, 34.0, 47.0, 75.0, 92.0, 73.0,
    38.0, 0.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 29.0, 87.0, 103.0, 115.0,
    84.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 13.0, 0.0, 51.0, 82.0, 97.0, 80.0,
    44.0, 0.0, 13.0, 0.0, 0.0, 26.0, 88.0, 83.0, 85.0, 42.0, 6.0, 0.0, 12.0,
    12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0,
    2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 12.0, 10.0, 12.0,
    2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 8.0, 6.0, 4.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 29.0, 30.0, 22.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 24.0, 24.0, 25.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 15.0, 26.0, 30.0, 38.0, 31.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 27.0, 29.0, 22.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 6.0, 22.0, 21.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 17.0, 18.0, 19.0, 17.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0,
    24.0, 25.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 12.0, 12.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 24.0, 28.0, 29.0, 12.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 28.0, 24.0, 13.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.0, 16.0, 15.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0,
    27.0, 40.0, 31.0, 26.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 14.0,
    25.0, 28.0, 31.0, 40.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 15.0, 26.0,
    25.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 16.0, 16.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 7.0, 7.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 15.0, 5.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 22.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 14.0, 12.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    16.0, 17.0, 17.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 0.0, 20.0, 28.0,
    26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 16.0, 26.0, 32.0, 29.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 34.0, 39.0, 31.0, 15.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 24.0, 27.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.0, 22.0, 23.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 18.0,
    22.0, 21.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 25.0, 27.0,
    12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 10.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 3.0, 3.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 4.0, 27.0, 37.0, 36.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 33.0, 33.0, 33.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0,
    4.0, 19.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 30.0, 31.0, 11.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 39.0,
    43.0, 46.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 32.0, 30.0, 13.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 12.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    16.0, 34.0, 33.0, 33.0, 14.0, 0.0, 0.0, 0.0, 0.0, 13.0, 26.0, 38.0, 39.0,
    29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 28.0, 30.0, 28.0, 7.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0,
    0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 8.0, 22.0, 22.0, 22.0, 26.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 23.0, 23.0, 23.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 18.0, 17.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 20.0, 21.0, 21.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    11.0, 23.0, 23.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 21.0, 21.0,
    21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 27.0, 47.0, 47.0,
    51.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 41.0, 41.0,
    29.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 35.0, 51.0, 52.0, 30.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 46.0, 47.0, 27.0, 30.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 37.0, 37.0, 42.0, 40.0, 18.0, 18.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 27.0, 28.0, 28.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 28.0, 31.0, 32.0, 33.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 22.0, 22.0, 22.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    16.0, 36.0, 46.0, 46.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0,
    19.0, 45.0, 46.0, 27.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 40.0, 43.0,
    22.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 31.0, 59.0, 53.0, 35.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 24.0, 30.0, 51.0, 31.0, 32.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 21.0, 44.0, 46.0, 27.0, 12.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 18.0, 27.0, 29.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 20.0, 21.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 18.0, 20.0, 20.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 19.0, 19.0, 17.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 16.0, 31.0, 42.0, 44.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 17.0,
    42.0, 41.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 31.0, 42.0, 44.0,
    23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 19.0, 49.0, 53.0, 33.0, 28.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 10.0, 20.0, 45.0, 49.0, 30.0, 24.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 21.0, 39.0, 54.0, 55.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    18.0, 41.0, 49.0, 31.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 40.0,
    43.0, 48.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 36.0, 40.0, 35.0,
    23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 23.0, 51.0, 54.0, 30.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 24.0, 24.0, 25.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 23.0, 23.0,
    23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 21.0, 21.0, 21.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 17.0, 17.0, 17.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 21.0, 47.0, 51.0, 52.0, 31.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 19.0, 42.0, 42.0, 45.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    19.0, 35.0, 35.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 37.0, 42.0, 42.0,
    24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 18.0, 18.0, 18.0, 0.0, 0.0, 0.0, 0.0,
    19.0, 20.0, 20.0, 20.0, 0.0, 0.0, 0.0, 16.0, 16.0, 17.0, 17.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 19.0, 49.0, 50.0, 38.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0,
    44.0, 45.0, 37.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 36.0, 37.0, 39.0,
    19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0, 43.0, 50.0, 51.0, 25.0, 0.0, 0.0, 0.0,
    0.0, 19.0, 41.0, 54.0, 54.0, 32.0, 0.0, 0.0, 0.0, 0.0, 15.0, 24.0, 46.0,
    47.0, 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 18.0, 18.0, 18.0,
    0.0, 0.0, 0.0, 0.0, 8.0, 19.0, 19.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 0.0, 0.0, 8.0, 31.0, 42.0, 43.0, 28.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 43.0, 43.0, 25.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 36.0, 37.0, 23.0,
    18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 21.0, 21.0, 21.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 19.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.50888644931942E-5, -9.43876234554427E-5, 0.000198407554538651,
    -3.60874643256829E-5, -0.000147672656135961, 2.9351977146719E-6,
    0.000233882398514276, 0.000112930643812503, -8.58202816129892E-5,
    1.8508715217904E-5, -9.81560324680459E-6, -4.56978015210399E-5,
    4.45865905309572E-5, -9.31364601542153E-6, -0.000121400478246564,
    0.000170037650604301, -8.70247829411185E-5, 7.55696641879192E-5,
    -9.64606632225748E-5, -4.274383783286E-5, 2.45121721498313E-5,
    0.00014940416906151, 0.0001589800920091, -2.76680439649985E-5,
    -0.000440656510722474, 2.41796770685842E-5, 0.000101854610103881,
    -0.000171068990047538, -2.88322079568207E-5, -7.96021541428457E-6,
    -0.000490814302736042, -7.29061965862717E-5, 0.000280821077895723,
    -1.46522557909077E-5, -6.49527131491182E-5, -1.3967643083328E-5,
    -1.9197233340786E-5, -3.17444859717527E-5, 4.77029034667073E-5,
    9.35283359406864E-6, 5.60755295063414E-5, -0.000358615485482464,
    -0.000186289690751449, -0.000324615840651801, -0.00162858346467851,
    -0.000521193387960573, -0.000437739682423508, 0.000142489880057252,
    -0.000216768956072802, -0.000116703674943502, -0.000128871133645297,
    -0.000334324109296193, -9.68125430029642E-5, 3.57059796854502E-5,
    -0.000103008005252445, -0.0011788588279222, -0.000306232493057427,
    0.000117598393539886, 6.64156459870002E-5, -1.90560345553017E-6,
    -0.00016235862749814, -0.000252182558683025, 0.000145601564577245,
    6.97613431611573E-5, 0.000299589368597845, 0.000154095497600908,
    0.000403434669072769, -0.000308112144782484, 0.00278696683917776,
    0.00272798035486278, -8.54331676418329E-5, 0.000167118718544605,
    -0.000141750787298336, -5.6841397865357E-5, 0.000341669594550739,
    -3.39314934540978E-5, -0.00012597836799012, 7.21496432940867E-5,
    0.000298970963503133, 4.30365136920759E-5, 0.000145382216529623,
    -0.000375341799457467, -0.000523567737125819, 9.44217870091806E-6,
    -0.00013103848062749, -7.11997912538908E-5, -0.000110579227816858,
    -0.000195114912292052, 8.09975842477245E-5, -0.000141214862316573,
    0.00012830648652324, -0.000118428613019054, -0.000374368417832723,
    -5.41816705383669E-5, -0.000700706954899709, -0.00269897721513246,
    -0.000241540054044862, 0.000259364964478112, 3.46540458090349E-5,
    0.000302445888685831, 5.87942145545425E-5, -0.000289198019310389,
    -8.82271468547814E-5, -0.000129791514431528, 0.000221122778021746,
    6.41699960933665E-5, 0.000347764652649102, -0.000219603565133459,
    -0.000797107756354032, 0.00312343005016168, -0.000152391685902439,
    -0.000146869647105491, -0.000197102442816365, -0.00010385957207254,
    4.87306511922865E-5, -0.000170430007676199, 8.38272661647116E-6,
    0.000197083328774574, -0.000139573657961739, -0.000369612812259037,
    0.00132722602839412, 0.0024026649907277, -0.000191357419779478,
    0.000694639289076205, -0.000973994113686831, -2.33336092695282E-5,
    -9.01044996888989E-5, -0.000307248498068664, -2.12808776121752E-5,
    -2.09276049336136E-5, -0.000228489152563045, 0.000164418916363754,
    1.14291859870426E-5, 0.000138377656191904, 4.34029284826488E-5,
    -0.00215455019265942, 0.00073319684584764, -0.000397272709587966,
    0.000321723055536446, 0.000296519713273183, -8.11841949090814E-6,
    -3.99876095728728E-5, 0.000233144860292511, 0.000204610684081284,
    0.000258252839462928, -0.000574171469646782, 0.00282650575509483,
    -0.000448547539980737, -0.000393628184638039, -3.62984421501849E-5,
    3.29853811795057E-5, -0.000189937322354122, 5.23134417771982E-5,
    0.000162207437256194, 5.79380732766379E-5, 9.77611458602741E-5,
    0.000110428152444845, -0.000339808798815767, 0.00285626541165701,
    -0.00025734561644113, 2.0456616148399E-5, 0.000168630528009716,
    -0.000106251767483579, -3.54649213724878E-5, -2.26766311733149E-5,
    0.000189126230895933, 0.000352336341682484, -9.9337496204495E-5,
    -0.000320851837684632, 0.0003905942833694, -0.000273994448320557,
    -0.00119755272867529, 2.32993563623116E-5, 6.83616031526228E-5,
    -6.78740651839237E-5, -0.000398174627944868, -0.000147681793711544,
    7.76640773956901E-5, 6.75441529982978E-6, 4.46066122828435E-6,
    -0.000118416673567408, -5.62288051208541E-5, -4.68061042854534E-5,
    -0.0012886136067024, 0.00056626563994556, -5.79188234671488E-5,
    -0.000149218729963813, 1.06433910509764E-5, 2.45518754170067E-5,
    8.32786027930891E-5, -1.60199392685151E-5, 0.000395237932614138,
    -0.000459800676102003, -0.000408585551390964, 0.000283762296209686,
    3.54207096431495E-5, 1.63411780773142E-5, 0.000116465685937306,
    -2.67227385701247E-5, -0.000113014502853411, 0.000208325830363955,
    0.000155387643846861, -0.000518619511182046, 0.000566227927447602,
    2.23121705468783E-5, -0.00217974847883269, -8.70177848581607E-5,
    -7.34550671406058E-5, 5.90226053626232E-5, -9.81109249022705E-5,
    -0.00012828316907182, 3.38127028811566E-5, -8.38353587695484E-5,
    0.000135557428251192, 0.00017304309202429, -0.000595065885277083,
    0.0019593494922815, -0.000139880682131855, -0.00158063396766287,
    -0.000926984499192456, -0.000114348554458223, -0.000184933048449578,
    0.00019351005262666, -4.17363517761438E-5, 0.000105983990493261,
    8.86165262374128E-5, 0.000121005360080875, -2.86913977451697E-5,
    2.04388939513123E-5, -0.000349984888899772, -0.00206618735379033,
    0.000691188709278719, -2.81841937603048E-5, 6.6421627995158E-6,
    -0.000217042775618243, -5.84329302423954E-5, 0.00046955027087653,
    -0.000275280879432315, 6.14070235479073E-5, -0.000287518125736764,
    -0.000308841800024798, 0.00204498530943444, -0.000497676440984645,
    0.00014656070954004, -0.000111341558654839, -0.000218216636686357,
    0.000276062454755423, 1.46552015977377E-5, 0.000152888604608565,
    -0.000120026723117859, 3.63306459478489E-5, -0.000160539602675026,
    5.54732101168536E-5, -0.000471769339338555, -0.00154772245920525,
    -2.11836059441561E-6, 3.25470715428421E-5, -9.21399181551613E-5,
    -9.74222054753897E-6, -0.000188911983734762, -0.00020313543425885,
    0.00015168157950355, -0.000293945885794398, -0.000595438288791347,
    0.00042675453092512, 0.000155521871380628, -0.000194073779670286,
    4.94988131262408E-6, 0.000187295255004392, -8.87621594292717E-5,
    -8.64036898479744E-5, 0.000216288900950906, -0.000349721861344779,
    0.00134319098504848, -0.000220046805065901, -0.000248215013993009,
    3.17496584594653E-5, 0.000120401579796149, 4.74353106008882E-5,
    0.000246333190546768, 8.8214146638144E-5, -6.72657742500745E-5,
    0.0010846589930856, 0.00114015970444026, -0.00100018269554128,
    -5.35711566030109E-5, -7.11179698975106E-5, 0.000160973341736203,
    -5.18742323953499E-5, 0.00025864667583992, -4.23906985579449E-5,
    -0.000361704540334321, 0.000337079833857545, -0.000720327421830912,
    -0.00115186939673097, -0.000328098016134501, -0.000168138534350703,
    2.0666420331434E-6, -0.000191366858238157, 0.000230661390013281,
    7.37393331774522E-5, 0.00016800713969414, -0.000164400039205156,
    -0.0002063813140458, 0.00123645015728566, -0.000392430064874454,
    -0.00174168016296564, 2.58217074803637E-5, -0.000248338347970333,
    -0.000160611486875985, 0.000277217940003375, -0.00022963733675461,
    -3.7279964726353E-5, -0.000119051132529562, 0.000132632802064184,
    -0.000311476019375974, -0.00114674931768006, 0.00210286908825173,
    0.00021625992066285, 9.62146551614637E-5, -0.000129092418136637,
    -0.000165469763774828, -0.000188483098576839, 0.000250335694964084,
    -0.000316728740904915, -0.000103133810559973, -0.000150556087432719,
    -0.0013555375184874, 0.00197582937093681, -0.000107404186979075,
    -4.42794771834527E-5, 0.000352670331911605, 0.00015559310896501,
    6.39106395966667E-5, 6.06704228365352E-5, 6.20662920151805E-5,
    0.000299700488965141, -0.000531587424442102, 0.00194536059793007,
    0.000114374763510366, -0.00122621032035696, 0.000114036645126978,
    -0.000259090236618164, 0.000127946428288287, 0.000185560406048283,
    0.000247444833868829, 0.000190692156750652, 7.01195685551403E-5,
    -6.43740877463368E-5, -0.000397795688046377, -0.00103907571942233,
    0.00244156389237935, -0.000625205443894144, -0.000353592567017047,
    0.000455329082847387, 0.00031248419691262, -1.52777990545896E-6,
    -4.29420752224277E-6, -0.000113270284028559, 4.70311522823962E-5,
    0.000201254601180505, -0.000218464703626706, 0.00188616125508606,
    0.00194822608242301, -0.000942078233469396, 0.000325816792090679,
    -0.000427538710809352, -0.000254254237188298, 0.00044912270190139,
    -0.000166325527176743, 0.000392028965646242, -9.28424098261618E-5,
    0.000288034162895632, -0.000465530981041058, 0.000478166935327759,
    -0.000925443573110168, -0.000669901416882753, -0.000489035863409368,
    0.000137248079215125, 0.000211752381576215, 1.09570268110336E-5,
    -9.79712873733096E-5, -0.000114870300606488, -0.00011494721134383,
    9.23982984836132E-5, -0.000211868286046468, 0.00208055433587822,
    0.000523426634429245, -0.000740971798338648, 4.89789654903793E-5,
    -0.000138642713102876, 8.10327132422668E-5, -1.26568891916501E-5,
    -9.5268778143428E-5, -0.000140196430297211, 6.74125491437282E-6,
    4.36088698523562E-5, 3.09713861259615E-6, 0.000248615469177045,
    -0.000248424826015013, -4.24509844482228E-5, 5.87086195554031E-5,
    -0.00021579413498434, -0.000137019929825246, -6.33731388047677E-5,
    -0.000267789688117202, -8.87718300714023E-5, 2.62289912939323E-5,
    6.05153083909103E-5, 0.000144939452535777, -2.45628122698846E-5,
    -1.19915569213814E-6, -7.08397327552464E-5, -1.23383817521986E-5,
    -6.571446270267E-5, -0.000110846510292058, -1.55310011466612E-5,
    0.000119632271582958, 7.6461386227924E-5, 1.55960869079604E-5,
    9.48461391413284E-5, 0.000122164632365981, -0.000113066162427331,
    6.80341369308643E-6, -0.000197999457847049, -7.05764717794772E-5,
    -0.000566478059906133, 0.00246106400447606, 0.00072060897068911,
    -0.000795059912125399, -0.00063994574709871, -0.000233969095086059,
    0.00024695485433792, -0.000291298151931173, -5.15505027660404E-5,
    -0.000350549229284078, 6.06220348856246E-5, 0.000267628836166733,
    8.34132522676299E-5, 0.000243194320892389, 0.000291313685792087,
    0.000175569375780023, 9.33774170943749E-6, -0.000209195837310753,
    -0.000267604700303752, 4.27329873063275E-5, -0.000270697545593582,
    0.00208802040044538, 0.000365376583103408, 0.000277883919824547,
    -9.2676918531675E-5, 0.00042566390572774, 2.55079268472054E-5,
    -0.000220468727574558, -0.00016705757985911, 0.00029993952982988,
    -0.000537073422437029, 0.000265820693699281, 0.000560278723123479,
    -0.000747188388989525, 6.82848360838933E-5, -8.03037387306608E-5,
    -1.38521890565994E-5, 0.000120489770888975, -9.63297121475058E-5,
    4.82211611334622E-5, 0.000378304469294538, -5.80021497958781E-5,
    -4.35750454811867E-6, -0.000151503597199052, 0.000520990900882907,
    0.000149576540181633, -1.41609323865025E-5, 3.47949152774257E-5,
    1.79019931216754E-5, -6.50416475693755E-5, 0.000380588738912725,
    -0.000277160304889764, 9.91876202737505E-5, 9.55990766156955E-5,
    0.00034241106556094, -0.000227260330023735, 0.0019283846084353,
    0.00257024695535629, 0.000232203252471987, -2.30991747253813E-5,
    0.000110907085877252, 0.000142673928597064, 0.000435230516995978,
    2.4673162553296E-5, -0.000331481013607464, 6.2551515264075E-5,
    8.69987328714956E-5, -0.000210359043116291, 0.00200603862415106,
    0.00142886592094927, -0.000136777860432722, 0.000209383336238765,
    -0.000148237979902632, 5.53920550889957E-5, -0.000135361138715078,
    0.000254412741069393, -8.92425618484115E-5, 8.52402062424838E-5,
    -3.85169489695463E-5, -0.000143254802580257, -0.0023486068622735,
    -0.00131386029367809, -0.000195666496082755, 0.000317317852624448,
    0.000106023108896835, -8.59816551902227E-5, 1.23356527268941E-5,
    0.000392159915868357, -0.00019416785639724, 0.000392021939182035,
    -8.44281178606376E-5, -0.000458819019121561, 0.000477509577521345,
    0.000146116858931175, 0.000202147743370888, 0.000207934534686625,
    -6.76727018031502E-5, -9.81122931373626E-6, -5.45423164534861E-5,
    -2.33400671503786E-6, -4.36531754139439E-6, -0.000249191144298029,
    -7.55743402057275E-5, -0.000169931463674858, -3.57054877203841E-5,
    2.01725408521446E-5, -0.000321797904394989, 0.00169390934976745,
    0.00011494113258861, 0.000222585453647843, 4.25015455199071E-5,
    -0.000258464933377956, 7.62313939664163E-5, -0.000153351563152863,
    0.00137063836075979, -7.0743329507972E-5, 0.000157298302112518,
    6.54492243064502E-5, 7.43852468764875E-5, 0.000362539442392363,
    7.41508030663428E-5, -0.000124763039352996, 0.000122411032159456,
    -4.28313378614654E-5, -0.00107918333151371, 0.00248374331812316,
    -0.000254258065775405, 8.66132145901504E-5, -0.000222531180964825,
    9.45702782273545E-5, -0.000150898612548402, -0.000141390877910151,
    0.000139200624695495, -0.000206356357815882, -0.000502077080404517,
    0.00326831362495612, 0.000215794423507677, 0.000392548546150106,
    4.62219584670517E-6, -7.2567909306885E-5, -0.000117564570001889,
    0.000183734972119936, 2.21473350552762E-5, 0.000359498788990941,
    -0.000452077412531898, 0.00142040324627255, 0.000449144024313527,
    -0.000110366547212224, 0.000145761717040863, 2.98412365985118E-5,
    0.000168776648579291, 2.59483433470022E-5, 0.000378958282981415,
    0.000132019348146719, -0.00020444644484244, 0.00138379088128464,
    -0.000222527425230313, -0.00111684288849793, 5.49657489821049E-5,
    8.64113628993043E-5, -0.000230821951801066, 0.000312799895956184,
    0.000341600026665021, -4.88826704141365E-5, 0.000445222371697805,
    0.000657056719981532, -0.00149300099398625, 3.41664575722111E-5,
    -0.000294189690095152, 0.000205304733095572, -0.000170742854026098,
    -0.000134664469084628, -0.000302634357357811, 0.000321544614574326,
    -0.00103702018050024, -0.000783216851996058, 0.000150524017931236,
    -0.00059418346073352, 0.00010251132963681, 8.50750301958185E-5,
    -9.28688643199403E-7, 0.000287952453268857, 4.95010161922663E-5,
    -0.00017107573406808, -5.7533641891428E-5, 6.68684282979708E-5,
    -0.000191268058556646, -0.000324104271594639, 0.000373080170518091,
    -0.000517039114425079, -0.000222676713083622, -0.000306522768302406,
    -5.49712772050893E-5, 0.000116884234098593, -0.000127463066670356,
    -0.00167847834368394, -0.000207183952106143, 0.000294663153004674,
    4.71913953714152E-5, 0.000155949184760889, -0.000175213841091472,
    -2.4082704846513E-5, 6.31145914478798E-5, 1.62476458147861E-5,
    -2.02165320599614E-5, -3.55973229615221E-5, 0.000245251048453192,
    0.000196858160540659, -0.000279508670077919, -9.02571946873005E-5,
    -0.000249199542997161, -0.000224411507342804, 2.78083005234914E-5,
    -6.74672555732132E-5, -2.75344684814766E-5, 5.81225244286229E-5,
    -1.03097705762847E-6, 1.74673505698843E-5, -1.47704807122269E-5,
    -0.000219781835235625, -6.37214727303428E-6, 0.000105904505061425,
    -0.000243649865544985, -0.000110910312734075, -0.000727153234142347,
    -0.00016813155944918, -3.72394408127978E-5, -6.97379233274955E-5,
    3.19201170798084E-5, -0.000290313740157946, -0.00225036061983439,
    -0.00234688703059834, -0.000174457981789554, 1.09912209252244E-5,
    2.66556808161908E-5, -0.000195871858116598, -2.98059527137526E-5,
    0.000110877581223143, 2.33900880052803E-5, -3.95380060010992E-5,
    3.25281853818578E-5, -0.000198743223891721, 0.000108568966765552,
    -0.000343387136986799, -0.00151244303626248, -0.000686555414419321,
    -9.06035315232723E-5, -0.000353468972323282, 7.59777165807192E-5,
    1.9620796712034E-5, -0.000123626187011062, 0.000128026951133158,
    0.000119258725220096, -9.65694321334514E-5, -0.000158350050401769,
    -0.000171867067271506, -0.000109488936034705, -0.000435166987464053,
    -0.00162374792982749, 0.00144502813398221, -0.000105802743544681,
    0.00012511761260521, 8.27244548705559E-5, -0.000108285741399216,
    -1.35491463253935E-5, -0.000120079639468215, -0.000150493209943579,
    -2.73688252109492E-5, 1.07018894332507E-5, 7.50413865456129E-5,
    -0.000198426528000136, 0.000186605665048775, -0.000545158658399539,
    0.000145980754212293, -0.000206879868365668, -0.000148606294210189,
    5.43370973796797E-5, -0.000153647107227017, -0.00015552186844642,
    -1.74264867819771E-5, -0.000327668331426884, -0.00170811146230917,
    -0.000161511433706866, 0.000198146537116708, -3.11954495008673E-5,
    1.77720141798136E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 2.0, 2.0, 3.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
    2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 3.0, 3.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 64.0, 64.0, 67.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 62.0, 62.0, 61.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 69.0, 76.0,
    76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 86.0, 87.0,
    82.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 68.0, 68.0, 225.0, 160.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    80.0, 82.0, 82.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    91.0, 90.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 94.0, 94.0, 93.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 88.0, 108.0, 107.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 166.0, 173.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 75.0, 75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 79.0, 79.0, 79.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 91.0, 91.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 84.0, 86.0, 86.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 92.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 71.0, 94.0, 93.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 88.0, 94.0, 94.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 87.0, 89.0, 90.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 74.0, 74.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    97.0, 98.0, 96.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 83.0, 89.0,
    91.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 73.0, 73.0, 75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 75.0, 75.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 73.0,
    73.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 71.0, 71.0, 71.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 80.0,
    79.0, 79.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 76.0, 75.0, 75.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 77.0, 77.0, 78.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 75.0, 76.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 72.0, 73.0, 73.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 85.0, 91.0, 91.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 57.0, 60.0, 60.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    56.0, 62.0, 62.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    72.0, 72.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.49410487919147E-9, -1.77399548310145E-7, 3.70261522386153E-7,
    1.10933684827816E-7, -2.23970826470141E-7, -2.41312192520006E-7,
    2.24593014832439E-7, 3.67002642052142E-7, -2.67115989496057E-7,
    -7.53401128432776E-8, -2.39786048268382E-7, 3.95259015343077E-8,
    1.40128265784856E-7, -5.0449417818721E-8, -2.5561857042407E-7,
    5.41861644732514E-7, -2.67296852459316E-7, -2.17153396868031E-8,
    -4.87081314470822E-7, -2.53297972569263E-7, 1.93342470680598E-7,
    3.91290785074195E-7, 3.29281684864193E-7, 2.43562761256044E-7,
    -1.30878005772663E-6, -9.92809973183704E-8, -7.6716243595483E-8,
    -4.77535975513675E-7, -1.67241094293936E-7, -8.64169096053968E-8,
    -1.79323590933717E-6, -1.17580928592029E-6, 2.95668600096595E-7,
    2.66040432754138E-9, -1.00205926635528E-6, 7.47826391424462E-7,
    -7.60822757389474E-7, 5.29608899067815E-7, 1.31497898336725E-6,
    -7.07655409173017E-8, 1.31416352784861E-7, -3.49100123995101E-7,
    -8.0763710325037E-7, -1.67714535125581E-7, -4.1876262746589E-6,
    -1.51907874493529E-6, -1.35419828418886E-6, 6.33047160317903E-7,
    -1.2575912577529E-6, -9.11616588911865E-7, 6.1979699991727E-7,
    -1.21489578694005E-6, -6.64377211508108E-7, -4.16365643831453E-7,
    -1.1816994226851E-6, -4.7686658735571E-6, 1.942216041931E-7,
    1.22306642656837E-6, 2.55909279579005E-8, 2.92846199080872E-7,
    -7.56900934282463E-7, -4.23108599877844E-7, 9.63649097015651E-7,
    5.35272244309114E-8, 7.51958718182502E-8, 1.05173693020968E-6,
    3.90131979742416E-7, -1.96457574696962E-6, 1.10063329100953E-5,
    1.18032941203291E-5, -1.06487709920148E-6, 1.58950612789816E-7,
    -2.00331917156241E-7, -4.27043048118589E-7, 1.05319681120588E-6,
    -6.32832610023559E-7, -4.60176855196956E-7, -7.80724994437348E-7,
    6.19052339701094E-7, 1.4975545285561E-7, 8.57103133260947E-7,
    -2.40196538705468E-6, -3.56199956553292E-6, -1.36567991797411E-7,
    -6.48114870760792E-7, -9.93667621422381E-8, -7.22595316818012E-7,
    -4.15451829657658E-7, 1.14169013737035E-6, -6.70130728852505E-7,
    4.13756557080286E-7, 3.44837999076188E-7, -1.42667935239251E-6,
    1.00889691313801E-6, -5.1448230158353E-6, -1.24460954932044E-5,
    -1.41868772688421E-6, 1.0805176458439E-6, -4.99142160376811E-7,
    8.79404532026942E-7, -1.30281400387913E-7, -3.37506590334608E-7,
    -8.80977624333446E-8, -8.7119658145969E-8, 1.75251710236245E-6,
    1.12374685080497E-7, 1.00640413455996E-6, -1.95647034304539E-6,
    -2.80867325738065E-6, 1.39895419194693E-5, -9.48706319521135E-7,
    -1.02983828370773E-6, -9.36126686065056E-7, -4.70276678586733E-7,
    4.98884227049495E-7, -8.75448635707654E-7, -7.19608130823918E-7,
    -3.66275256972699E-8, -3.13354263470957E-7, -1.94826606822422E-6,
    4.82290527180209E-6, 1.06339581927186E-5, -5.78354484932417E-7,
    2.75087753327894E-6, -3.78199265612641E-6, -1.36573701895719E-7,
    -2.72010344343454E-7, -7.82912052029208E-7, 2.88172537934638E-7,
    -1.29150545889164E-7, -7.86987473913468E-7, 4.32395316376573E-7,
    -1.45583864850988E-7, 5.60740838095341E-7, -6.41356403242106E-7,
    -6.37092000482258E-6, 3.27652141824388E-6, -9.50321388854807E-7,
    1.1402935594566E-6, 1.47488986734301E-6, 2.97540050523437E-7,
    -4.38842022503323E-7, 3.6977468850577E-7, 7.69084187691554E-7,
    3.59733453385681E-7, -2.82678165407473E-6, 1.31497157260033E-5,
    -1.65167090028126E-6, -1.00317190046979E-6, -1.69386442028517E-7,
    2.20967278315658E-8, -1.11177605642733E-6, 8.59184635034552E-9,
    2.5640981086724E-7, -5.14649484068615E-7, 9.12579117069153E-7,
    6.59058417424665E-7, -1.73691095014144E-6, 1.2393718224917E-5,
    1.31377932901256E-6, -4.70010576219794E-7, -8.69850298894143E-8,
    -5.05771220650826E-7, -5.00691700456803E-8, -1.34211014002985E-7,
    1.08590108107097E-6, 1.49445077570049E-6, -4.81178114946739E-7,
    -1.03561674564225E-6, 2.86612721613332E-6, -3.54532997470486E-7,
    -1.10412631470571E-6, -7.58736073926876E-8, -1.12392092133602E-7,
    6.24004984654116E-7, -1.68593694487237E-6, -2.55409316076376E-7,
    -2.01018544428269E-7, 2.42499061466872E-7, 6.7092192903751E-7,
    -5.53735938578005E-7, -3.71709361995714E-9, -9.17057721506511E-8,
    -5.17724663200877E-6, 2.12270232247335E-6, 3.6050310238019E-7,
    -1.46097840558611E-6, -1.76612145294918E-7, -1.46910091856262E-7,
    -4.21493002531954E-7, -1.0819756691491E-7, 1.66013415486594E-6,
    -1.60314781982502E-6, -2.28838531175407E-6, 1.31855797331778E-6,
    -4.4648773768429E-7, 8.34564393236865E-7, 1.78953510716252E-7,
    -1.53143902170136E-7, -2.46333429285991E-7, 2.23170276394322E-7,
    1.27988568747181E-6, -1.95608564892315E-6, 2.0120293325487E-6,
    4.12027818739383E-7, -8.32113660803551E-6, -5.11673215256132E-7,
    -6.66924802696801E-7, 2.32214295303133E-7, 1.58146575580485E-7,
    -2.73159267833711E-7, -2.5885785382138E-7, -3.40058389074717E-7,
    5.17532061297101E-7, -1.34483141313484E-7, -2.01906683348746E-6,
    8.29634086673676E-6, -8.55106275997853E-7, -7.05688052401037E-6,
    -5.331056217565E-6, -1.96106204922012E-7, -7.5133488368258E-7,
    4.78966699905011E-7, -2.53772858813166E-7, -2.23173544828722E-7,
    -1.06904899866864E-7, -5.03096223290511E-7, -1.07691758767323E-6,
    3.38216887159359E-7, -1.68705881726538E-6, -1.15121185832892E-5,
    2.67312379689747E-6, 4.70872077960056E-7, -3.21125246051653E-7,
    -3.60084171713845E-7, -9.40860549651979E-7, 1.67248819978484E-6,
    -9.72905998441655E-8, 7.88040897166919E-7, -7.61206299976333E-7,
    -1.3753082899678E-6, 9.51696558552981E-6, -2.7553326196332E-6,
    -5.30951160533263E-8, -5.55782971890339E-7, -8.23538894016302E-7,
    4.62820539733961E-7, -1.8614470679338E-7, 3.6456300477927E-7,
    8.60963167632938E-8, 2.16183671348629E-7, -3.52438034098849E-7,
    7.49438764085889E-8, -2.31081583716467E-6, -7.09175600509408E-6,
    5.92519605702017E-8, -2.54033563784497E-7, -6.28406914536686E-8,
    -7.33812148290726E-8, -4.4543452794574E-7, -5.37513360808425E-7,
    2.38412182507591E-7, -1.95102105630781E-6, -2.41146789337226E-6,
    6.78808801985317E-7, 9.66094307498806E-7, -2.06622352742104E-7,
    -2.45860825581953E-8, 8.28202062602034E-7, -2.79799726573193E-7,
    2.53270097402823E-7, 5.21959840559191E-7, -1.57042720896721E-6,
    6.28400128794134E-6, -3.96892208368286E-7, -6.67880201384665E-7,
    1.71291376079959E-8, -3.30571308428118E-9, 5.71001784737183E-7,
    5.82623373654881E-7, 5.85463006639939E-7, -7.33929214802639E-7,
    5.21166151190657E-6, 4.52485320743296E-6, -4.6621215502651E-6,
    -4.06907834475514E-7, -3.77198538359274E-7, 2.19107384362711E-7,
    4.4613309538667E-7, 6.78871592938739E-7, -1.81571441900757E-7,
    -9.79372839501478E-7, 2.23890421825013E-6, -3.39614014095902E-6,
    -6.50858123064973E-6, -5.61771539394895E-7, -6.51105477069143E-7,
    5.04764027809947E-7, -4.7302548742423E-7, 1.18492304918153E-7,
    6.2593890699723E-7, 2.25154554909133E-7, -4.76833024864913E-8,
    -3.88388779145771E-7, 5.78021391421226E-6, -2.73545264356802E-6,
    -7.36162836505729E-6, -1.25854709586891E-7, 1.01311418732427E-7,
    -2.4493926382967E-7, 8.68522491510016E-7, -4.49035460647981E-7,
    -3.59942605087129E-7, -4.57327990864501E-7, 9.7630691705417E-7,
    -1.15741571118937E-6, -4.28724714774336E-6, 9.69546435448436E-6,
    -3.65055029882358E-7, 5.84917514407975E-7, -1.1437231197523E-7,
    -5.16831136071503E-7, -1.25920216062376E-6, 4.40168796942987E-7,
    -1.58452807850124E-7, -6.01540877474423E-7, -1.29745133069474E-6,
    -6.91429891442259E-6, 7.73709453329524E-6, -8.24428882608626E-7,
    -1.60433899682309E-7, 3.85139482359275E-7, 3.87498661625482E-7,
    1.55648140260511E-7, -2.25426914707579E-7, 1.3330360732476E-7,
    8.47938485481661E-7, -2.25516569540719E-6, 1.06445084276273E-5,
    -6.51275584577723E-8, -3.90755370909557E-6, -3.16030720574179E-7,
    -2.18566298159962E-7, 1.83337954361846E-7, -3.85050626936049E-7,
    2.15967622951252E-7, 3.20340474446912E-7, -3.34198990315077E-8,
    5.80891409135394E-7, -1.96278950645229E-6, -6.06437286563132E-6,
    1.02671818940343E-5, -3.30270159718521E-6, -8.45984639565835E-7,
    6.8997796765221E-7, 5.03024061132106E-7, 4.63210769365495E-7,
    -1.11055366701321E-6, -2.82053157645784E-7, -2.85046802575159E-7,
    1.24525848346231E-6, -1.76919769174963E-7, 8.51156004952796E-6,
    7.84587895871366E-6, -5.18545663835074E-6, -2.83085068217527E-7,
    1.2633467445623E-6, -1.34809021167392E-7, 7.74813061411302E-7,
    -4.92765882378587E-7, -9.77813618103785E-9, -6.9076253358336E-7,
    1.38260368796646E-6, -1.5112932513661E-6, 4.14507184010007E-7,
    -3.00956252286062E-6, -2.80020484578445E-6, -1.63023191620609E-6,
    1.01395021483522E-6, 4.31429038920023E-7, -1.09356838917398E-7,
    -5.00453881534323E-7, -7.22933620564296E-7, 1.99410833858678E-7,
    2.63193794129055E-7, -1.86393646700719E-6, 1.10059448599331E-5,
    2.99246837865195E-6, -3.45299804279721E-6, 9.72719906863088E-8,
    -5.63606183826856E-7, 6.4209953365635E-7, 4.07006920346603E-7,
    -1.66362302590796E-8, -1.57907563319543E-7, 5.11448237595762E-7,
    3.66828296165416E-8, -6.02690110017071E-7, 6.24014822816583E-7,
    -2.56332859640482E-7, 1.98602706039377E-9, 1.56852667009787E-7,
    -1.86620664557656E-7, -6.99758459286613E-7, -1.48544287786505E-7,
    -5.67085043365324E-7, -4.84464244805607E-7, 2.82785374620417E-7,
    1.72688664882575E-7, 3.16418655629067E-7, 8.99786519834002E-8,
    -2.12385514674543E-7, 6.66699566282875E-8, 1.20294113504488E-7,
    -5.08061011739271E-7, -3.1085147068144E-7, -2.23131670235809E-8,
    1.72342416881501E-7, 1.78550333052874E-7, 1.79074199742459E-7,
    3.31115484199973E-7, 5.17584975000614E-7, -3.25614031207397E-7,
    1.7198905061677E-7, -4.04839932977245E-7, -3.86949632447566E-7,
    -2.19733714550109E-6, 1.01129898163263E-5, 2.47106481644178E-6,
    -4.2877838065396E-6, -1.53170081448786E-6, -1.99875135398831E-7,
    1.40704623670148E-7, -1.10150253536715E-6, 1.96656091650658E-7,
    6.91854982229552E-9, 2.89734745451057E-7, 1.09574437895136E-6,
    9.05237968344502E-8, 1.16490773968049E-6, 1.02557681813486E-6,
    1.01893533733138E-7, -1.01458836093806E-8, -8.66516317518478E-7,
    -1.73786396252731E-7, 2.73517866319532E-7, -1.26674909764026E-6,
    7.63749507561332E-6, 3.10572183880318E-7, -1.27288592046781E-6,
    -7.58936314055265E-7, 1.05644221149292E-8, 1.47944984624301E-7,
    -6.63307241570243E-7, -4.46486528392907E-7, 9.5299644483624E-8,
    -6.34403791190321E-7, -5.04903966478115E-7, 2.42488655719327E-6,
    -2.83118256450797E-6, 1.15049759892018E-6, -3.53840025184861E-7,
    3.86129299530184E-7, 6.84028822860344E-7, -7.26284747882822E-7,
    -7.7037525333298E-8, 8.95232273863473E-7, 9.79544166467801E-7,
    5.551100143567E-7, -3.07123083170378E-7, 2.48714243400388E-6,
    2.46017589360714E-7, -6.20026748449746E-10, 1.53801784867172E-6,
    3.39377061383626E-7, -6.22121755798326E-7, 2.77551200162598E-7,
    1.08220158583694E-7, -5.93594622008764E-7, 7.14561909611901E-7,
    1.11555083757145E-6, -5.8995379823872E-7, 8.08807516969814E-6,
    9.78363405899369E-6, 7.55911908789913E-7, -3.61675424130769E-7,
    4.38585976874804E-7, 2.80510170173829E-7, 7.82707019074568E-7,
    -5.91522876286206E-7, -1.58504456731264E-7, -4.80443514068553E-7,
    -4.84309641990382E-8, -1.58991161854962E-6, 9.67180606467222E-6,
    6.28003330052786E-6, -8.20641460157837E-7, 9.67276678101713E-7,
    1.38570946334123E-7, -5.61565509102009E-7, 9.83571656379307E-7,
    -2.00058042014324E-8, -2.39515336714945E-7, -5.40408642493893E-8,
    -5.72813040413377E-7, -8.96104318854474E-7, -8.21255669679632E-6,
    -5.97396351427621E-6, -8.92614716526964E-7, 6.04447326201256E-7,
    -2.60916059192648E-7, -8.81841553400511E-7, -9.05124006309929E-8,
    7.0240331129622E-7, -3.29363559014481E-7, 1.94309674571266E-6,
    -3.62548826591351E-7, -8.27046055873861E-7, 1.20563632883395E-6,
    3.07939706101253E-7, 6.80060353436396E-7, 4.67508979860223E-7,
    6.797211978546E-7, -3.5839505515165E-7, -1.54333995473555E-7,
    1.09017927800487E-6, -8.66118545399356E-7, -7.736142069989E-7,
    -3.63238794406747E-7, -5.56025287693279E-7, -4.81467682612302E-7,
    6.24694451644623E-7, -1.74387847021742E-6, 7.58942723534603E-6,
    7.62990955359378E-7, 3.66326720969741E-7, -1.86458897269905E-7,
    -1.18327746803715E-7, 8.52073260658279E-7, -3.4167207092944E-7,
    6.66807235950896E-6, -4.01666208320355E-7, 8.18766375815269E-8,
    6.90295092024862E-7, -1.10433527871343E-7, 2.66591506575947E-7,
    2.93045417930519E-7, -8.87584561669597E-7, 6.41767978526817E-7,
    -9.39409204985468E-7, -5.66518584404062E-6, 1.02159485362362E-5,
    -1.04373349064348E-6, -2.28272963852339E-7, -7.77302723057036E-7,
    6.63216151740139E-7, -6.84182977518784E-7, -1.58964689788325E-7,
    1.31707894153467E-6, -1.18254729613008E-6, -3.63561853786407E-6,
    1.27004956949194E-5, 1.11844694158112E-6, 8.40109077941752E-7,
    -1.67472905665631E-7, -4.2697563376895E-7, 1.72876048178001E-8,
    2.40531112099481E-7, -2.14684641054511E-7, 5.58887370599701E-7,
    -5.3703171917786E-7, 6.74811106567392E-6, 2.00232111751399E-6,
    -5.77549953547124E-7, 8.41509586900815E-7, 3.58500292338836E-8,
    5.40470639583405E-7, 1.04666942402068E-7, 1.20528700645903E-6,
    8.16970757684991E-7, -5.78356510984577E-7, 6.22010095893319E-6,
    -1.5177651823554E-6, -5.16131346307554E-6, 7.64073366388552E-8,
    -7.88425460280528E-8, -6.2273483155405E-7, 1.36139802060154E-6,
    8.08505378722481E-7, -5.31484320445565E-7, 2.24218056574842E-6,
    2.11382684820126E-6, -6.14648504207771E-6, -4.4959492961693E-7,
    -1.13937330016903E-6, -7.46675023225105E-7, -5.46724115801493E-7,
    2.05840630430627E-7, -1.23824677260111E-6, 1.69867948962687E-6,
    -4.55444242788072E-6, -5.0027068884175E-6, -3.43284260703973E-7,
    -6.22516037125179E-7, -2.64336051068854E-7, 1.28037078687261E-6,
    -4.13325993805052E-8, 7.46599210023199E-7, -5.98184931209229E-7,
    -2.69181279163052E-7, -6.27082554422842E-7, 8.92613028186136E-8,
    -4.7070387986662E-7, -1.78294975207625E-6, 4.86514978945004E-7,
    -1.35164800734605E-6, -7.75811655113712E-7, -6.98090271198565E-7,
    -4.4346319844085E-7, 3.34216499435752E-7, -1.24338629882085E-6,
    -3.45237024375332E-6, -7.89917749094707E-8, 4.85459411145583E-7,
    -9.100217462444E-8, 5.27763077345424E-7, -7.40017911997938E-7,
    1.11769408009368E-6, -8.83735759719813E-7, 3.27451900424231E-8,
    -6.92934517127778E-8, -8.14322211198948E-7, 5.61302040266285E-7,
    5.03447892729415E-7, -5.05818067617038E-7, 3.73946243817475E-7,
    -6.05618752283844E-7, -5.44882124669933E-7, -2.07945308787178E-7,
    1.59935883422849E-9, -1.17882781194147E-7, 2.67443946076025E-7,
    -2.53614326349177E-7, -9.30026330558698E-8, 1.96156862598169E-7,
    -1.21728640365122E-6, 7.3304234593606E-8, 4.21664942744072E-7,
    -8.47165824669116E-7, -3.5299513308455E-7, -3.22563261786472E-6,
    -5.15476106497661E-7, -1.52229329313349E-7, -1.56097254174572E-7,
    -9.27227074317848E-8, -1.25626534204116E-6, -9.61605755843194E-6,
    -9.04051615877529E-6, -2.02991824879962E-7, -2.15364208951968E-7,
    -1.8734957313029E-10, -3.96591040722627E-7, -4.6975796146588E-7,
    2.51914966501862E-7, 7.31200096743119E-8, 8.12373194945429E-7,
    -1.784791280317E-7, -6.4795554857411E-7, 8.30583740319009E-7,
    -1.66703543017537E-6, -6.70085385604237E-6, -4.1520475407399E-6,
    -1.08799492303505E-6, -1.04411885365705E-6, 1.67850557828373E-7,
    1.55858326322445E-7, -2.81900880918289E-7, 8.82645990556158E-7,
    1.00402528408928E-7, -6.50376370128825E-7, -9.38302541711888E-7,
    -7.86626174414669E-7, -3.1544476547014E-7, -1.12032643946279E-6,
    -8.0555288397156E-6, 4.57778332922177E-6, -1.01064597111581E-7,
    2.08333516282326E-7, -7.58565631121415E-8, -4.43965430020097E-7,
    -2.05461354369822E-7, -7.40070109677944E-7, -4.06199671048471E-7,
    -5.09881182479505E-7, -8.48736838302862E-9, 8.79518811022071E-8,
    -8.03505560928851E-7, 6.9485172852719E-7, -1.98119003275563E-6,
    -1.07258308239692E-7, -2.60370547151883E-7, -7.38988276093041E-7,
    1.08969236511184E-7, -5.37870452122213E-7, -3.61036552549866E-7,
    2.57298408391653E-7, -1.55793982263258E-6, -6.80106339416733E-6,
    -3.35929808850143E-7, 1.07086205244994E-7, -2.20914138550445E-7,
    3.06974573795734E-7, 2.39040881390268E-6, 2.1877967536632E-6,
    2.64321697044229E-6, 2.74717412858624E-6, 3.11108641891506E-6,
    3.13526479591749E-6, 2.85165828953434E-6, 2.64286677825311E-6,
    2.63621323651625E-6, 2.42162357297142E-6, 2.39058897351681E-6,
    2.8018236220105E-6, 2.89844217401683E-6, 2.93817886147459E-6,
    3.11071894490607E-6, 3.80447585457131E-6, 3.33903573333192E-6,
    3.43989247491351E-6, 3.19652461368132E-6, 3.53682655550244E-6,
    7.11361830721374E-6, 8.69698457431248E-6, 9.22545434339934E-6,
    9.50045646603944E-6, 5.7321821858405E-6, 3.18216531008014E-6,
    3.60442569305928E-6, 1.82664444780842E-5, 5.71189671348819E-5,
    5.94021160574881E-5, 5.93869799624486E-5, 2.72799273563802E-5,
    3.810219910598E-6, 3.87989727672066E-6, 5.51640279440955E-6,
    4.72478320015473E-6, 4.61328552415324E-6, 3.81368222524027E-6,
    6.11312590622409E-6, 2.43639162699009E-6, 4.42421821297301E-5,
    6.72004061197198E-5, 6.98395202743371E-5, 6.84429434824554E-5,
    1.82889700786285E-5, 7.861950574149E-6, 7.77092922675049E-6,
    5.57742613831297E-6, 7.91531843522668E-6, 6.09002951711802E-6,
    5.43663637075566E-6, 4.88893058730682E-5, 6.521886527364E-5,
    6.62955539306596E-5, 6.65575286954908E-5, 1.92844287157005E-5,
    4.57195863576563E-6, 6.5005734962793E-6, 4.23114877222317E-6,
    4.71928867668403E-6, 5.84988395509514E-6, 3.82136623083065E-6,
    5.44922865216159E-6, 3.5697078526991E-6, 7.42481301484944E-6,
    6.03687467805805E-5, 9.74543372611766E-5, 9.75855928529458E-5,
    9.81361202284441E-5, 7.05694783120525E-5, 7.28450606198253E-6,
    5.07490954022495E-6, 3.48485336065433E-6, 2.78839083957566E-6,
    5.6719382147238E-6, 3.6228397478166E-6, 3.83639301858771E-6,
    4.49108266777185E-6, 2.1451293228392E-5, 7.43629109985054E-5,
    9.37707679177969E-5, 9.47580091546443E-5, 8.21292937920961E-5,
    2.64981128797952E-5, 5.21724554650799E-6, 3.12661407593749E-6,
    4.95929136418215E-6, 2.9219653237184E-6, 5.91788092317344E-6,
    5.81067686351294E-6, 4.95184318610343E-5, 9.28088942819589E-5,
    0.000114493554139766, 0.000116356699086926, 9.36549755342024E-5,
    4.73559825535589E-5, 7.62389009494836E-6, 6.74406754115049E-6,
    5.06971684674428E-6, 4.49397069354233E-6, 2.522835679735E-6,
    4.00908955823791E-6, 4.25940869063024E-6, 4.10573248439071E-6,
    2.17665054931236E-5, 5.4467360129703E-5, 9.6484379420124E-5,
    0.000103740841124755, 8.67435898641273E-5, 7.92318617741986E-5,
    4.58347062978903E-6, 4.83193638961499E-6, 4.44610040029361E-6,
    3.77449859153658E-6, 3.72247815512597E-6, 5.21101401219164E-6,
    2.00447820696483E-5, 5.98433506335309E-5, 8.40554664232941E-5,
    8.45914628811997E-5, 8.36352323010352E-5, 8.00520247147969E-5,
    5.90227541146379E-5, 5.4425096700558E-5, 1.45917423249731E-5,
    2.67583876190344E-6, 3.1907567423465E-6, 4.82469107444696E-6,
    4.02153481593531E-6, 3.8876161033736E-6, 5.30020642515062E-6,
    5.83291168205931E-5, 8.02797267938678E-5, 8.34602033641605E-5,
    8.16214338128515E-5, 2.42425443127868E-5, 1.40575208332694E-5,
    5.60377142453254E-6, 6.9285806265451E-6, 6.99989906011299E-6,
    4.88888408812709E-6, 4.12922433925892E-6, 1.66493160672535E-5,
    7.07159762191663E-5, 8.11260284703162E-5, 8.15171993275835E-5,
    7.11834690919941E-5, 7.17503297760532E-6, 5.64800765978293E-6,
    3.92091852537777E-6, 4.07150305895025E-6, 6.0567994868363E-6,
    3.12792566516204E-6, 3.41279591268577E-6, 6.98710952847597E-6,
    6.64577458816057E-5, 7.12119373398235E-5, 7.12967195389555E-5,
    6.9973532813827E-5, 9.44154233464657E-6, 4.65280691654789E-6,
    2.95049734912365E-6, 4.07872623615248E-6, 3.36646791791692E-6,
    3.29959440775226E-6, 5.98005085372387E-6, 5.1132883111193E-5,
    9.31466555822772E-5, 0.000104055532070124, 0.000105469836633711,
    7.39401041146433E-5, 1.76014976501851E-5, 5.80813350268318E-6,
    4.84875461926249E-6, 4.87607832924183E-6, 8.8313635224114E-6,
    4.04313851796136E-6, 3.42643003423944E-6, 2.86451086161188E-6,
    4.26727198722884E-5, 6.20454799300768E-5, 9.08452150898216E-5,
    9.43848378034554E-5, 7.08876923478078E-5, 3.44752128440845E-5,
    4.86236089513007E-6, 7.08427754994256E-6, 3.47094937232034E-6,
    3.65429402002591E-6, 1.92890152222464E-5, 5.87494344013384E-5,
    8.90348030438984E-5, 9.73315079137697E-5, 7.52508972606571E-5,
    6.66011861691507E-5, 5.10542570491394E-6, 4.63042670795263E-6,
    3.33166326871799E-6, 3.09474396262906E-6, 2.31583324926192E-5,
    6.42097913427002E-5, 8.98082669763975E-5, 0.000108711535685639,
    9.00139998171709E-5, 7.30493191064925E-5, 3.03802388139782E-5,
    2.97748737603446E-6, 4.82220285708784E-6, 3.10163359654589E-6,
    3.02891878739277E-6, 3.45494860859665E-6, 2.86074125227078E-6,
    1.72404887640438E-5, 5.61435498026696E-5, 6.014179721411E-5,
    7.67054582862145E-5, 9.71554015943556E-5, 9.55495911403471E-5,
    0.000101303900024976, 3.32639990769386E-5, 3.76995183333288E-6,
    4.10946298741298E-6, 4.2808029404153E-6, 3.21382561628909E-6,
    3.62375492496178E-6, 2.42160249238637E-6, 3.43304613154578E-5,
    7.29545607740517E-5, 0.000100961236531297, 0.000102760674794344,
    8.93736176883342E-5, 3.28711412322088E-5, 5.9247718269757E-6,
    3.61584952530683E-6, 3.70941475806615E-6, 4.80869541513546E-6,
    7.43824552587013E-6, 5.00504932320708E-6, 6.19545804256471E-5,
    7.52959161561384E-5, 7.61255401404417E-5, 6.74589726670609E-5,
    1.11936250484687E-5, 3.11292485036554E-6, 3.5044431051498E-6,
    3.77953298497483E-6, 3.01305079272393E-6, 2.11193339000006E-6,
    3.97923573981974E-6, 4.30981646153727E-6, 3.05394664017907E-5,
    6.31737557281446E-5, 6.58851348905156E-5, 6.59527089579451E-5,
    2.70837844962547E-5, 3.32199134078801E-6, 2.70500203778135E-6,
    2.57205695107415E-6, 3.60247589845647E-6, 6.37575924075456E-5,
    6.55507249531898E-5, 6.75140961296578E-5, 6.71966055166708E-5,
    1.44214135509876E-5, 5.02680392394262E-6, 5.16752579518979E-6,
    3.58557527053553E-6, 3.27691027593669E-6, 4.77057647885974E-6,
    1.48766303284767E-5, 5.57269176841959E-5, 5.85048608148025E-5,
    5.85073032686988E-5, 4.48062629624209E-5, 3.3730006785485E-6,
    4.48800118255683E-6, 3.53812104917684E-6, 3.26145524522379E-6,
    4.60776682234047E-6, 4.51524343932621E-5, 7.04722587569931E-5,
    8.41946899312608E-5, 8.70237299400592E-5, 6.07693628537915E-5,
    1.81255793292198E-5, 5.06243739409674E-6, 4.44047585840595E-6,
    3.92296426374457E-6, 3.95882672908101E-6, 4.16556515775368E-5,
    5.50505955318498E-5, 8.79928072660817E-5, 8.87959021108583E-5,
    6.87011121037674E-5, 2.61886471071455E-5, 5.48057586001438E-6,
    4.83322279198482E-6, 4.56879298530239E-6, 4.29566461185433E-6,
    4.53171760204955E-6, 4.65123162439393E-6, 4.80860356310314E-5,
    7.57165765995109E-5, 8.83724397919879E-5, 8.99619306363473E-5,
    6.80079013774773E-5, 2.77783783619659E-5, 4.77835784075599E-6,
    4.56657626554301E-6, 4.88744753465119E-6, 5.24142374686901E-6,
    4.61892762919582E-6, 4.39894195469134E-5, 5.55387280478217E-5,
    9.46195711669138E-5, 0.000102311079184063, 8.2113085354902E-5,
    5.78775138489047E-5, 6.41303007087149E-6, 3.99159918666199E-6,
    2.3008359016659E-6, 3.20324865970858E-6, 5.8859058063666E-6,
    4.01390770266955E-5, 7.25405549501591E-5, 9.42603443153434E-5,
    0.000100821042650271, 8.00446807109749E-5, 5.31067070946147E-5,
    6.6278637794763E-6, 4.76004229769002E-6, 2.725800423591E-6,
    3.13202874456816E-6, 1.17337155259713E-6, 1.79315255164596E-5,
    6.94692802692367E-5, 9.81914863419444E-5, 0.000114801035966721,
    0.000112674842627399, 7.47859493966386E-5, 1.83012719278307E-5,
    3.70962908384962E-6, 2.83418660903019E-6, 2.07787900912405E-6,
    4.12874595223668E-6, 3.72439315429689E-6, 2.07255410649778E-5,
    5.86622881695595E-5, 9.02311050804934E-5, 9.9294957806358E-5,
    8.3022767185864E-5, 6.42453198626671E-5, 1.53011795208218E-5,
    5.57756965011147E-6, 4.95504141797356E-6, 4.62605314026116E-6,
    3.92704627219181E-6, 5.20521281945633E-6, 8.60120273795676E-6,
    6.49458811864753E-5, 9.40616601555425E-5, 9.67172743740651E-5,
    0.000103877707249693, 5.58223569784318E-5, 2.28803828726773E-5,
    5.04560955680037E-6, 6.16165513321026E-6, 3.83132552491554E-6,
    5.75196561835411E-6, 4.0836201271829E-6, 2.12922715418838E-5,
    6.2428882216181E-5, 8.25097262045854E-5, 9.02531692654633E-5,
    7.70917881193944E-5, 6.91088390183048E-5, 1.28152831599406E-5,
    7.37397623606084E-6, 5.13707767905972E-6, 2.95399039196918E-6,
    2.62758220256501E-6, 2.99639614447782E-6, 7.27144283553492E-6,
    5.39777415469871E-5, 8.72821662795839E-5, 0.000102716365726279,
    0.000111217173957904, 6.92188504462448E-5, 1.35578786986814E-5,
    3.1860719098798E-6, 3.40578907383871E-6, 3.36596355752642E-6,
    2.97780033401281E-6, 2.09257425654942E-6, 8.53554050939821E-6,
    7.28902056258754E-5, 4.37289664282821E-6, 4.01151936701142E-6,
    3.69878982568611E-6, 2.47232968202935E-6, 3.48660217694154E-6,
    3.01208133507764E-6, 3.028382392078E-6, 3.55610659597454E-6,
    2.3884877575317E-6, 3.96903234402125E-6, 3.9205184389925E-6,
    4.11554078496497E-6, 4.22576909716903E-6, 4.21967218737642E-6,
    3.86697984089161E-6, 4.08038883082562E-6, 3.06297945180986E-6,
    3.07949654161874E-6, 3.54426974820329E-6, 2.56695987287902E-6,
    2.83824096229164E-6, 2.67174094394124E-6, 2.83566947052313E-6,
    8.90166793623077E-6, 7.72261073956285E-6, 8.67389599720929E-6,
    8.28788758283413E-6, 1.26920439394346E-5, 6.26778976288583E-5,
    6.69975193610545E-5, 6.78810488342444E-5, 6.06451262024109E-5,
    1.12157172714345E-5, 1.76407353657194E-5, 6.99190449919092E-6,
    2.70137232396738E-6, 2.0993881017729E-6, 4.86642743204881E-6,
    3.37813961223675E-6, 3.97780910426916E-6, 4.87493456567912E-6,
    6.24870447683507E-6, 4.12991195968019E-6, 5.78254589374941E-6,
    5.75377578722769E-6, 3.16358979702257E-6, 3.10373376799361E-6,
    1.63453516502298E-5, 5.70273260299197E-5, 6.09738600092587E-5,
    6.10552940765585E-5, 5.30200943285293E-5, 4.96996419162653E-6,
    7.71908390765174E-6, 4.73129238278121E-6, 4.96516459016135E-6,
    3.25311589409923E-6, 4.37693643447189E-6, 3.16863471717889E-6,
    5.06585707834711E-5, 6.3286609719412E-5, 6.38843234496874E-5,
    6.20637257315669E-5, 1.22691256474674E-5, 6.41531779583696E-6,
    3.98234110006159E-6, 3.58556699365089E-6, 4.23128315861282E-6,
    4.22051568109264E-6, 5.62468363133469E-6, 6.70751335431149E-6,
    3.28147622607241E-5, 4.05987722851446E-5, 4.06458116123722E-5,
    3.86527544406804E-5, 6.07037258235291E-6, 3.7406063091226E-6,
    7.04053125189026E-6, 3.89601331842904E-6, 4.00819997726635E-6,
    3.57822631124383E-6, 2.95068451204384E-6, 5.56623031867729E-6,
    6.42576765630111E-5, 9.8759847987182E-5, 0.000105704307397404,
    0.00010737188067604, 7.17464404687464E-5, 5.27777763188703E-6,
    3.30055125833526E-6, 3.60081661838205E-6, 3.30198111662893E-6,
    4.48795552643811E-6, 3.58856397223941E-6, 7.47931034405415E-6,
    5.60183357337011E-5, 9.29037636987214E-5, 9.50522184041274E-5,
    9.99473637343184E-5, 6.57706172973844E-5, 4.91807923268559E-6,
    5.75247169572778E-6, 3.83098284473306E-6, 4.72521892025695E-6,
    5.91594028120179E-6, 4.35433574922292E-6, 2.88902410715252E-5,
    5.53217689269608E-5, 8.07236789869238E-5, 8.12251894019065E-5,
    7.09776082802378E-5, 3.22684892107032E-5, 4.79221163220949E-6,
    4.22733665737522E-6, 3.45686816872897E-6, 4.93701772468862E-6,
    4.75963681922878E-5, 8.57480773669545E-5, 9.58619807414148E-5,
    9.59360155470924E-5, 7.36094534672519E-5, 5.63211711229479E-6,
    7.76744444300724E-6, 5.96333876737277E-6, 5.78964787371808E-6,
    4.0551366391796E-6, 4.46101961910311E-5, 5.52615239974492E-5,
    5.62715042880741E-5, 5.47101152736644E-5, 4.59509507983742E-6,
    4.9784704321796E-6, 3.5138089242664E-6, 1.1738255358365E-5,
    5.94363070631394E-5, 6.21690986090775E-5, 6.2157204891823E-5,
    5.58014966365797E-5, 4.43614185203746E-6, 4.18389964007543E-6,
    4.5308501983272E-6, 4.64677629280254E-5, 4.91077548132988E-5,
    4.97787542439848E-5, 5.03902476945079E-5, 5.36200999923805E-6,
    4.31865670437779E-6, 4.84768337073974E-6, 5.38482194659233E-6,
    5.16324694789728E-6, 1.90984097450707E-5, 5.90132885048258E-5,
    9.82742409565403E-5, 0.000101056474938251, 9.07639684076738E-5,
    6.27793412436279E-5, 5.17708856447711E-6, 4.31544629182851E-6,
    6.28835105929732E-6, 6.30108385686107E-6, 1.80965398979339E-5,
    5.16260049309334E-5, 9.51656608139782E-5, 9.94256769120367E-5,
    8.99334094659817E-5, 7.42599305955909E-5, 9.30442811429702E-6,
    5.63111207845004E-6, 4.24091016805485E-6, 5.33072335717281E-6,
    4.05057810794366E-6, 5.41219850849692E-6, 4.24442302924147E-5,
    7.56916602689507E-5, 8.00604312681692E-5, 8.23631009244041E-5,
    6.02672877898356E-5, 5.02711205468087E-6, 4.69288338327329E-6,
    4.45292971673053E-6, 5.2506526630888E-6, 5.20115726795699E-6,
    6.89228134088893E-5, 0.000112299127466728, 0.000120827544158559,
    0.000119576993224155, 7.77913358143456E-5, 2.17541622505794E-5,
    7.00096168196962E-6, 5.55381388987179E-6, 6.32980421663699E-6,
    6.32378767108573E-5, 0.00011054227535894, 0.000123635504937794,
    0.000123304165113158, 8.38470267283971E-5, 2.33834223157649E-5,
    6.16209233940652E-6, 7.3346914867639E-6, 6.26689743533629E-6,
    5.1916851334366E-5, 6.84247359358734E-5, 0.000106067494544123,
    0.000108308048314746, 8.73184583137394E-5, 3.27250084200902E-5,
    3.45319490045117E-6, 4.27310927544787E-6, 3.0734421981992E-6,
    6.5445957412986E-6, 3.14825883399938E-6, 4.60609154323265E-6,
    4.12305124410584E-6, 4.83920036101422E-5, 5.79324245027515E-5,
    6.00975629853881E-5, 5.95274390434125E-5, 1.16787129091408E-5,
    8.51514662029265E-6, 1.02769231981693E-5, 4.74422389693661E-6,
    4.4517881073929E-5, 5.8445968149216E-5, 5.95626745706248E-5,
    5.96596647563805E-5, 1.48334387377484E-5, 4.29242017426348E-6,
    4.61754440734072E-6, 3.69929649167543E-6, 3.48980816708564E-6,
    5.04971301460107E-6, 6.52285149150768E-6, 5.86371420722567E-6,
    7.8605238720337E-6, 9.46594219620062E-6, 9.1748734375373E-6,
    8.9558575953418E-6, 6.09864320819487E-6, 5.77195377183024E-6,
    4.92961357555968E-6, 5.30451680528198E-6, 3.73011260282318E-6,
    3.52572307528065E-6, 3.40960380186356E-6, 3.66287811653018E-6,
    4.45359677428307E-6, 6.23424838801031E-6, 1.15460511375179E-5,
    1.38073786519245E-5, 1.40539347210455E-5, 1.80412194578567E-5,
    2.85704919961542E-5, 2.98958685746815E-5, 3.044447778027E-5,
    1.32455922311301E-5, 4.61975426783239E-6, 3.79345519070433E-5,
    8.96756698616027E-5, 0.000101861438000085, 0.000103288552382586,
    0.000100766584215778, 3.58285681053934E-5, 3.80362629393251E-6,
    3.80477280286408E-6, 3.37430358234177E-6, 3.41354743415773E-6,
    3.49721975497695E-6, 2.85624436421603E-6, 4.01066519167502E-6,
    4.94589903126949E-6, 2.22697928302512E-5, 8.01437751553286E-5,
    0.000102856028663355, 0.000103977024304048, 9.43111997169689E-5,
    3.00340888666615E-5, 5.58828003776478E-6, 6.2794929232604E-6,
    4.08588390483696E-6, 3.90679002509024E-6, 3.64635510576236E-6,
    4.31304057834735E-6, 2.86009495767779E-6, 3.60514416402886E-6,
    1.38543021121258E-5, 4.85638652198351E-5, 7.27115009425769E-5,
    7.42458090363112E-5, 6.86989708767907E-5, 3.56002542462639E-5,
    4.40852227803682E-6, 3.83497488939386E-6, 4.212078540889E-6,
    3.70321762447527E-6, 3.7732081363474E-6, 4.39630396399895E-6,
    3.86356083580552E-6, 3.5459841312687E-6, 4.98203012374737E-5,
    6.27866914308701E-5, 6.39147773849797E-5, 6.33608715085605E-5,
    8.24037181009855E-6, 3.36288582242499E-6, 2.22756929690863E-6,
    3.60053095848774E-6, 3.21796773051262E-6, 3.00919143966675E-5,
    5.85274010675089E-5, 6.05252935176465E-5, 6.02342469426079E-5,
    2.79431397335538E-5, 4.1564416813742E-6, 3.4115033042747E-6,
    3.78065923705723E-6, 4.62586749020353E-6, 5.80223556721483E-6,
    4.21702194475772E-6, 7.40948133560412E-6, 5.80645656799796E-6,
    7.09106421235151E-6, 7.24064944890498E-6, 6.07322906991207E-6,
    9.14219305425637E-6, 4.7962691654774E-6, 5.04218019938409E-6,
    5.46019607839739E-6, 7.08749048453329E-6, 7.22207596751114E-6,
    6.75200861793186E-6, 4.82992945915562E-6, 9.81250798748027E-6,
    6.63300861956864E-6, 7.85686873201941E-6, 9.05264324409014E-6,
    8.28599948784329E-6, 2.10611801956147E-5, 3.24635007962563E-5,
    3.38031729201974E-5, 3.53394931843541E-5, 6.85378378211606E-6,
    6.00053253211181E-6, 8.45728197764643E-6, 4.52743267088724E-5,
    0.000194098089238665, 0.000200765744394538, 0.000200971122292239,
    7.60074488861986E-5, 1.00890112210488E-5, 9.32377687795456E-6,
    8.07149214734864E-6, 1.18916703501204E-5, 1.02949673757791E-5,
    1.29298592005337E-5, 2.22574940600175E-5, 6.59159543722674E-6,
    0.00014214583667982, 0.000208501907417226, 0.000215778238317319,
    0.000215486580846835, 2.32161589938821E-5, 8.65094902578923E-6,
    1.13475969287394E-5, 1.45016891221867E-5, 1.01126768143806E-5,
    1.15715616591222E-5, 1.42134863231135E-5, 0.00014074321951669,
    0.000219094785434623, 0.0002245430754572, 0.000224536198414772,
    1.90974800514884E-5, 9.6715995836E-6, 2.29757740134982E-5,
    7.88900455255683E-6, 7.79302182612247E-6, 9.80490779408211E-6,
    1.03356561371464E-5, 1.4021083990197E-5, 8.04227897436819E-6,
    1.49999526575126E-5, 0.000147742710182712, 0.000204925290496809,
    0.000209134286864534, 0.000205554492625152, 0.000265912868528718,
    1.49965633173676E-5, 1.15672899357582E-5, 9.82346071949936E-6,
    3.59222298193669E-6, 1.92189974699081E-5, 6.83855696700954E-6,
    8.25780445257258E-6, 8.43211374427361E-6, 4.17938886013235E-5,
    0.000244601487780211, 0.00021964288644738, 0.000220975795628403,
    0.000233665576783134, 8.64002519352095E-5, 8.36822306803259E-6,
    9.17438727264871E-6, 1.1122551770085E-5, 5.75860171287827E-6,
    2.16401460393841E-5, 1.05805907595079E-5, 0.000158838453768717,
    0.000257870611005722, 0.000243492629215998, 0.000246762063615371,
    0.000247683504129472, 6.31173863878714E-5, 1.43446064892888E-5,
    2.66323444101845E-5, 7.81965954855742E-6, 1.65232508551754E-5,
    6.23086741303243E-6, 9.76186046953144E-6, 8.91001746709873E-6,
    8.77364339809793E-6, 6.05962794469404E-5, 0.000186645614167384,
    0.000221468797594136, 0.00023943163319176, 0.000244801096006459,
    0.000259116174214362, 6.38369449076707E-6, 8.68367713658798E-6,
    6.31498100571182E-6, 6.62380273959037E-6, 8.89901134147388E-6,
    8.45494994987449E-6, 3.97235077399736E-5, 0.00020968109989087,
    0.000212989728755133, 0.000214814993872846, 0.000190626467303655,
    0.000232647637882973, 0.000206862000412825, 0.000203093306939271,
    2.29411167141407E-5, 5.19220325095649E-6, 7.26295526984367E-6,
    8.76054639472664E-6, 8.6205307816178E-6, 8.22011012065994E-6,
    8.89361627719901E-6, 0.000151882967350967, 0.000225328570778079,
    0.000231854649782746, 0.000232199409772332, 3.03032841558307E-5,
    5.72095273195436E-5, 1.34546490383165E-5, 2.71080977239073E-5,
    2.65127571202225E-5, 1.20138144491681E-5, 1.18799930543998E-5,
    4.00031436081657E-5, 0.000159577070244878, 0.000188928960068708,
    0.000191442198331944, 0.00019751854561994, 9.53875767412204E-6,
    7.6873901075178E-6, 7.55137657866865E-6, 7.80370765601199E-6,
    7.26108393861147E-6, 5.76858796195746E-6, 5.5269604430511E-6,
    1.44161227503752E-5, 0.000216845472545256, 0.000230729478762451,
    0.000230655289486617, 0.000255434383263672, 4.0231766133707E-5,
    1.07962223568505E-5, 4.29194171779307E-6, 6.30679795229133E-6,
    5.39323713940714E-6, 5.93218630863612E-6, 2.02198401971373E-5,
    0.00017846527233731, 0.000261962830972567, 0.000245450674690088,
    0.000247304090900767, 0.00023408597798559, 4.87822993968253E-5,
    1.356691256158E-5, 1.19539116564626E-5, 1.32092624159245E-5,
    1.32808143133881E-5, 8.01853751481072E-6, 6.83613823060833E-6,
    3.82409668251749E-6, 0.000139495882555727, 0.000214564544599606,
    0.000223029278997471, 0.000222356857795338, 0.000189821340288827,
    0.000113947021416947, 1.16188886593943E-5, 9.98534830221226E-6,
    7.74966879933428E-6, 7.68403519104296E-6, 3.80371388352345E-5,
    0.000208272633938651, 0.000236425902142884, 0.000247690865547902,
    0.000251715225051182, 0.000237268050453691, 1.06336056042131E-5,
    1.60920071461222E-5, 7.87708504425558E-6, 5.65209246572275E-6,
    4.99721888451914E-5, 0.000229377105834955, 0.000266182041179706,
    0.000265092701163078, 0.000169656568208938, 0.000157992494473419,
    3.9681649992502E-5, 7.68662839043061E-6, 8.8048434631187E-6,
    8.11884989405442E-6, 8.42178106839613E-6, 8.6403656743576E-6,
    6.93942467571474E-6, 4.1647349662671E-5, 0.000142667937616043,
    0.000149897621800291, 0.000148327736656219, 0.000233873579554392,
    0.00025242959373661, 0.000259047979015071, 8.39237551769859E-5,
    6.45538376331651E-6, 7.01168837824624E-6, 1.00708656381057E-5,
    7.9682944145102E-6, 8.14584899029088E-6, 5.12800088639959E-6,
    0.000100278424357691, 0.000254783803911971, 0.000256535104825764,
    0.000255377535125339, 0.000235329572370351, 0.000116301398006982,
    1.39958979708735E-5, 5.75235991366703E-6, 6.84355121227466E-6,
    7.20754969460098E-6, 2.83242042174438E-5, 1.31102237979715E-5,
    0.000175594438736936, 0.000206676759540012, 0.000210062749394374,
    0.00022125848409852, 1.63378357760047E-5, 7.05161765316479E-6,
    6.30748213394436E-6, 7.05338104506321E-6, 9.17357816608359E-6,
    4.62886905046212E-6, 8.2160460318826E-6, 1.00907364775472E-5,
    8.70336278934774E-5, 0.000221542870870228, 0.000230287552180588,
    0.000230671496153883, 4.33009221901416E-5, 1.07282516235854E-5,
    6.80251333877556E-6, 7.46448408665416E-6, 7.55123354533848E-6,
    0.000192672404251372, 0.000235791570717948, 0.000242767272301547,
    0.000242865853605913, 3.13008504184542E-5, 1.10043927120023E-5,
    1.28784426271111E-5, 7.38707140312372E-6, 5.46632780219793E-6,
    1.30725647718093E-5, 2.94243649984876E-5, 0.000200807760898222,
    0.000209177167271413, 0.000209446609281519, 0.000189804103110968,
    4.74131075279192E-6, 9.67321016610304E-6, 7.69505995229364E-6,
    8.64209870502699E-6, 9.22356012974214E-6, 0.000155336847517985,
    0.00017821820883883, 0.000206012510340349, 0.00020943526859338,
    0.000202720144076206, 2.14293922686029E-5, 9.14760375833713E-6,
    9.43263398459369E-6, 7.95413286781186E-6, 1.03342149432421E-5,
    0.00014665783407799, 0.000181960401894995, 0.000236580193393092,
    0.000243852007279116, 0.000236669006579994, 5.02044884161848E-5,
    1.21407772580859E-5, 1.20992992820259E-5, 1.17321815969122E-5,
    9.83222375670139E-6, 1.02645414503946E-5, 1.29457802846442E-5,
    0.000169302581236668, 0.000208030855729391, 0.000211684692061741,
    0.000219777146193999, 0.000212131618574468, 3.34916946881897E-5,
    1.12478565612694E-5, 1.11442016348267E-5, 1.15377878480729E-5,
    1.53832514140181E-5, 8.51924801111077E-6, 0.000125737477399005,
    0.000182505113392914, 0.000210579527751452, 0.000209076118404657,
    0.000184455588707939, 0.000168047549316161, 1.10496817391584E-5,
    1.31624455702371E-5, 6.13557158096866E-6, 4.8129736005578E-6,
    1.13796133510836E-5, 0.000124296475573439, 0.000216917309638302,
    0.000221489735403122, 0.000221292944223275, 0.000182738778812745,
    0.000179338885058238, 1.65116493835482E-5, 1.32535197400243E-5,
    6.35688914466095E-6, 1.02485838250455E-5, 3.10327484144078E-6,
    4.15619901565002E-5, 0.000235871609581453, 0.000277682143234576,
    0.000273555870006473, 0.000281933575768288, 0.000200607265696385,
    3.41434100851697E-5, 9.00280824586835E-6, 5.78810868995409E-6,
    4.5445912888846E-6, 8.68067556260797E-6, 7.48923502912004E-6,
    4.01350242313085E-5, 0.0002106735152003, 0.000239535665597438,
    0.000236829434753224, 0.000211423981924952, 0.000201247721884313,
    2.39012913998873E-5, 8.20553041822229E-6, 1.14520896734734E-5,
    1.38496567601757E-5, 1.0102078590991E-5, 6.15363096165643E-6,
    1.85310400931181E-5, 0.000231772263670244, 0.000215097973190853,
    0.00021253528240535, 0.000260372529167864, 0.000191051244137624,
    2.46474882171426E-5, 8.66082615993072E-6, 2.51984105346397E-5,
    9.09219964273121E-6, 1.36712308153747E-5, 8.41348788023417E-6,
    4.44229482872253E-5, 0.000216281645881108, 0.000235336232887369,
    0.000231718928360835, 0.00020546932724882, 0.000205270501995181,
    2.27378149981816E-5, 1.1620886130957E-5, 1.36078253473298E-5,
    7.28709189768988E-6, 5.40756591879675E-6, 5.36604980317322E-6,
    1.64534635439646E-5, 0.000190252234555733, 0.000234850519803366,
    0.000227254264009385, 0.000289614054764386, 0.000169973251546798,
    2.6213188535409E-5, 6.46392883130162E-6, 7.32324330335927E-6,
    9.99775502160871E-6, 9.13453484065325E-6, 5.45080410127442E-6,
    1.81424283383324E-5, 0.000257205906518799, 1.09596902192399E-5,
    1.09642737150098E-5, 1.20206453890105E-5, 4.50667204950299E-6,
    7.09207957957558E-6, 5.17433020569278E-6, 5.55784484778798E-6,
    5.23280123416358E-6, 5.82381476277214E-6, 8.40512474753852E-6,
    6.65865253717089E-6, 7.36330472858258E-6, 7.0806914528484E-6,
    9.16693248826382E-6, 9.61733179322206E-6, 9.23022805722015E-6,
    7.94831485694321E-6, 8.12848400661433E-6, 8.74393930604796E-6,
    5.07819663976887E-6, 5.93969870093067E-6, 6.16346218323057E-6,
    5.4586828462232E-6, 2.50249475832939E-5, 2.33821421074437E-5,
    2.241003352767E-5, 2.28434263849694E-5, 2.03448143493785E-5,
    0.000186297319777801, 0.000196579189119066, 0.000200130847147958,
    0.000213749653495412, 4.32129570877736E-5, 1.17675364530474E-5,
    9.32420983859863E-6, 6.73148739051182E-6, 5.62204453801782E-6,
    5.31865278916219E-6, 7.27520142762917E-6, 8.09623839788057E-6,
    8.4105775716837E-6, 2.7052153959257E-5, 8.35078951625562E-6,
    2.28097675371975E-5, 2.26849387929978E-5, 6.009654394285E-6,
    6.51572302767534E-6, 3.49753591304099E-5, 0.000165010935032712,
    0.000175760542919402, 0.000175621610929302, 0.000176273712898155,
    8.32172742556533E-6, 1.09469858616017E-5, 1.18486276198579E-5,
    1.11495603994183E-5, 5.75222411415273E-6, 7.81266474688219E-6,
    5.37850593866092E-6, 0.000164779272023816, 0.000206919767695325,
    0.000208919654511058, 0.000209659578630415, 1.17001584734468E-5,
    2.2538255863618E-5, 7.30064469525475E-6, 9.43592465140707E-6,
    1.16915660709663E-5, 5.86858496019193E-6, 8.47614955870911E-6,
    1.4378554081706E-5, 0.000109098429626353, 0.000142300970003995,
    0.000143049886538384, 0.000142881390913662, 1.53814672897848E-5,
    6.79835194495823E-6, 2.24377929578268E-5, 9.86833430802189E-6,
    8.63572411583917E-6, 8.18697835083573E-6, 7.00304629446347E-6,
    1.21171965817088E-5, 0.000206133253683722, 0.00021723607485042,
    0.00020920576021068, 0.000224250192721086, 0.000204489048775246,
    1.8324945800982E-5, 6.84753478690411E-6, 7.98402281675406E-6,
    8.73239483520915E-6, 1.33608283436909E-5, 6.07850233835504E-6,
    1.81203699256052E-5, 0.000191381370890958, 0.000223486781796715,
    0.000233489040225088, 0.000251316449304874, 0.000233060212296912,
    8.16819814345163E-6, 1.79045817372367E-5, 7.00121859356785E-6,
    9.72918714133106E-6, 1.99386485667393E-5, 9.1542023834071E-6,
    7.81866116991131E-5, 0.000198665753010216, 0.000217007150612706,
    0.000219879329681146, 0.000235192865158078, 6.08932654003616E-5,
    7.08151112874999E-6, 1.4857645456211E-5, 8.99627504869157E-6,
    9.0884035679161E-6, 0.000162365940303282, 0.000215932321030369,
    0.000239672393180456, 0.000242982588243558, 0.00023933204974138,
    1.04543494445945E-5, 2.13506573379237E-5, 1.24454209177528E-5,
    1.47081075898458E-5, 1.18379640302222E-5, 0.000158860391456002,
    0.000206823482761058, 0.000210025634380247, 0.000209939470499135,
    5.94830660094743E-6, 7.21179242643411E-6, 5.83236386407401E-6,
    2.29112821473664E-5, 0.000213981045975408, 0.000221813060369434,
    0.000221920921394449, 0.000221059101910895, 1.50712274521268E-5,
    1.04262855645871E-5, 9.99343416026336E-6, 0.000180552308684926,
    0.000189937871409115, 0.000190802778360276, 0.000206081073413429,
    9.06208068207312E-6, 7.3180904242741E-6, 1.09306141758193E-5,
    1.23956413891565E-5, 1.17173397218265E-5, 4.14577768068239E-5,
    0.000213110875215499, 0.000217332820298518, 0.000216144802336477,
    0.000202586868030604, 0.000177720518207558, 1.07687307459785E-5,
    1.05464876067824E-5, 1.14228407516015E-5, 1.15982751751429E-5,
    3.53107262483565E-5, 0.000200626481750737, 0.000216980014041915,
    0.000232975719279315, 0.000243365675084885, 0.000242476167602054,
    3.20308921067568E-5, 1.4714736254249E-5, 8.27045064635257E-6,
    7.574017137144E-6, 6.62856218825104E-6, 1.03017553844639E-5,
    0.000156499801700265, 0.000206955489730536, 0.000230568993001765,
    0.000231356067842675, 0.000226699897789953, 1.26648942165882E-5,
    1.34274593110033E-5, 7.8497532365845E-6, 1.09755653599022E-5,
    9.42962204891836E-6, 0.000216218543641904, 0.000290760004066407,
    0.000270434096657594, 0.000273072531325484, 0.000221422709986779,
    2.38860142440389E-5, 1.08518538879296E-5, 1.23834296628986E-5,
    1.35800926166342E-5, 0.000207569402992533, 0.000288519155132951,
    0.000265190522576387, 0.000263971448288524, 0.000202793570517056,
    1.42189664460602E-5, 1.41798015835404E-5, 1.15071035726365E-5,
    1.26368477622886E-5, 0.000194457615666765, 0.00023231706652901,
    0.000257593706854882, 0.000269710893210582, 0.000262821248395729,
    7.62160582835757E-5, 4.640330501613E-6, 5.68163027775729E-6,
    5.7556927050958E-6, 2.41068455100703E-5, 7.85014679108249E-6,
    1.59545517033479E-5, 7.17797204362512E-6, 0.000147620330646747,
    0.00021552330765132, 0.000221541024838604, 0.000221334879337183,
    2.1326160010381E-5, 2.26641030087163E-5, 2.49609923927452E-5,
    8.05173883160311E-6, 0.000130442929293477, 0.000214673128229075,
    0.000219634259310163, 0.000220250317076931, 2.04491566009238E-5,
    1.49609912734452E-5, 1.56656940588024E-5, 7.32146117822463E-6,
    8.96469992577839E-6, 8.15898813380496E-6, 1.73134924438589E-5,
    9.18418333937525E-6, 1.85049418176102E-5, 2.71036101903619E-5,
    2.62573618142145E-5, 2.5032593507078E-5, 1.3412562569395E-5,
    1.30059679358963E-5, 9.9933888574957E-6, 8.76036342936558E-6,
    7.40195546568819E-6, 1.02564228094586E-5, 1.18834710586776E-5,
    1.24231026694366E-5, 1.25521759492692E-5, 1.15900973457643E-5,
    3.25732281749054E-5, 3.75885122246021E-5, 3.69461113817438E-5,
    4.54786489914812E-5, 9.59721614444912E-5, 9.95355066665189E-5,
    0.000100693811091316, 1.58824587435852E-5, 1.07706272265075E-5,
    0.000112128541866211, 0.000238354779404608, 0.000225413827953219,
    0.000233885468246058, 0.000261951724976752, 5.11894024268106E-5,
    7.19195112001675E-6, 7.10401821116832E-6, 7.74878296393154E-6,
    5.68560569515975E-6, 4.71904410747539E-6, 4.71362117354849E-6,
    1.02149859009447E-5, 1.75987219262929E-5, 5.34024730796451E-5,
    0.000245080325641367, 0.000239558019364973, 0.000247992116036288,
    0.000274419344305198, 7.29891204270367E-5, 9.90068739701881E-6,
    8.08128910974578E-6, 9.04245116264978E-6, 8.03164258493281E-6,
    4.89209194623596E-6, 1.49239878225634E-5, 5.62486395862316E-6,
    6.97242628129481E-6, 3.05627307978072E-5, 0.000180022754866146,
    0.000175850468693063, 0.000180100386796363, 0.000199429754185523,
    0.000136074926224514, 1.03246994696719E-5, 9.0639024851237E-6,
    1.15020807169558E-5, 9.61247188848824E-6, 9.82654126951091E-6,
    6.79375642479543E-6, 7.39549168984213E-6, 7.23952594770451E-6,
    0.00016417005892018, 0.00021637043334069, 0.000220858101321891,
    0.000220839198201028, 1.04756688863393E-5, 6.46294894613708E-6,
    4.09019412316226E-6, 7.95326473171179E-6, 1.04890371016734E-5,
    9.39677052234668E-5, 0.000214278897507895, 0.000219855323861365,
    0.000220117612684773, 3.85022773946537E-5, 1.07740532143626E-5,
    9.97002480581788E-6, 6.5674284852007E-6, 1.11261980300687E-5,
    -5.81417802082694E-6, -4.57410496243388E-6, -4.19090733802833E-6,
    -5.80219817934621E-6, -5.80395125224818E-6, -5.50159945415242E-6,
    -6.20571354119405E-6, -3.77221993778238E-6, -7.25602930950814E-6,
    -4.91119437927555E-6, -4.04958330779896E-6, -4.06112877940922E-6,
    -4.73715389959926E-6, -5.43702790964034E-6, -6.53321927533671E-6,
    -7.13321354688316E-6, -8.56695553467089E-6, -7.38735529967747E-6,
    -6.55025017850873E-6, -9.48102587479701E-6, -2.10238957618044E-5,
    -1.53716557585198E-5, -1.48327457358272E-5, -1.16901378813156E-5,
    -2.22415191475587E-5, -7.45844158595069E-6, -8.26035397809089E-6,
    -6.14659415865872E-5, -0.000105624148556416, -9.94119047785598E-5,
    -9.85233507346804E-5, -7.31215135780557E-5, -6.69574548431137E-6,
    -5.93255008951999E-6, -1.98656482289177E-5, -5.91085477761447E-6,
    -1.52541093759167E-5, -7.37025372455572E-6, -4.6560875414869E-6,
    -4.8101744243707E-6, -0.000125702935393041, -0.000132780481392352,
    -0.000125827308275393, -0.000120130309854022, -6.8452200751655E-5,
    -3.19789525533879E-5, -2.70213454112894E-5, -1.09428978312E-5,
    -3.02686762146738E-5, -1.33006535466075E-5, -1.31882372051547E-5,
    -0.000143537249307763, -0.000120430794475779, -0.000116735870784881,
    -0.000116926601379545, -6.67497677617501E-5, -1.02914613607147E-5,
    -7.77045821254273E-6, -8.57826927105897E-6, -9.14421528792033E-6,
    -1.881739253282E-5, -7.12162076123646E-6, -1.0939734851169E-5,
    -1.08904903250066E-5, -2.11831283799499E-5, -0.000138782726643405,
    -0.000255599778921426, -0.000252206241496083, -0.000260326524190299,
    -8.58871742948096E-5, -1.23304144855135E-5, -1.31772819836357E-5,
    -1.03041565509613E-5, -1.10861199549767E-5, -9.37854126250337E-6,
    -9.9290526499856E-6, -1.05336703211876E-5, -1.29528974249045E-5,
    -7.65245543447913E-5, -0.000127952766812558, -0.000209905791859966,
    -0.000203108448234362, -0.000192482769116773, -5.55531374598857E-5,
    -1.5028300841081E-5, -8.14855116035095E-6, -1.51389085019946E-5,
    -6.90103364460653E-6, -8.76545624732409E-6, -1.54899218969218E-5,
    -0.000134067470811521, -0.000208900822861356, -0.000234773152915712,
    -0.00023092006424846, -0.000200498744613946, -0.000158691012411059,
    -2.61383673427166E-5, -1.04659919199766E-5, -1.74144014086928E-5,
    -6.28672360096458E-6, -4.70628897360989E-6, -7.61429171037324E-6,
    -1.06017311440432E-5, -9.63648283611238E-6, -7.0931843502794E-5,
    -0.000101819859107129, -0.000164695255732649, -0.000164910164078544,
    -0.000146298624615335, -0.000112092983067495, -1.15240710738432E-5,
    -1.85125338248256E-5, -1.36638293059254E-5, -1.02751224494691E-5,
    -1.01760557672737E-5, -1.26607083375791E-5, -6.63438770579026E-5,
    -0.000113116750211601, -0.000157905414544727, -0.000153702362071324,
    -0.000151370544919665, -8.743327621823E-5, -9.64765719579326E-5,
    -8.2464755968278E-5, -5.03982720911113E-5, -5.60967092055387E-6,
    -5.62825280346759E-6, -1.21150948642209E-5, -7.85493763129598E-6,
    -7.33829288360258E-6, -1.45673860208864E-5, -0.000167498505301328,
    -0.000170774862036058, -0.000172232786432008, -0.000177226081422134,
    -7.87942843590486E-5, -1.62587759466679E-5, -1.70690734452856E-5,
    -9.82832702731731E-6, -9.67806086783772E-6, -9.16765556083215E-6,
    -9.23629239861476E-6, -5.56057288557296E-5, -0.000167874403644134,
    -0.000169306934316817, -0.00016830750096515, -0.000104393400817937,
    -1.93433518719372E-5, -2.20231731843491E-5, -6.52106036494431E-6,
    -6.46854342385146E-6, -2.14659516461578E-5, -7.97066900222301E-6,
    -9.10850698866902E-6, -2.22288247962421E-5, -0.000136780633395742,
    -0.000133497197001001, -0.000134215259620234, -8.47018635526069E-5,
    -1.30286045790474E-5, -1.16496511902565E-5, -7.5232369660538E-6,
    -1.24060508265481E-5, -8.2960411973769E-6, -8.21535005108182E-6,
    -1.06771325818634E-5, -0.00012859092428572, -0.00018961710271484,
    -0.000176091789764176, -0.000172677349514767, -0.000136332951351528,
    -4.48943308689077E-5, -1.09763125243936E-5, -7.90485588297396E-6,
    -8.68734094505806E-6, -3.02194208966878E-5, -8.74419284665922E-6,
    -7.28402519241588E-6, -7.40457937547272E-6, -0.000116152018650708,
    -8.54048647005788E-5, -0.000146517260675238, -0.000143573653628228,
    -0.000120566910985593, -6.63726376833597E-5, -9.97805875942432E-6,
    -2.63899125329157E-5, -9.63290053941565E-6, -7.96153924877404E-6,
    -6.60454613762568E-5, -0.000116802140421558, -0.000145512930531551,
    -0.00013560563675535, -0.000117789771305079, -0.000118909775539482,
    -1.12610315301677E-5, -5.57949290495262E-6, -6.04810679330525E-6,
    -7.76598847250333E-6, -7.61115030834805E-5, -0.00013140917583298,
    -0.000169932633405193, -0.0001927589156251, -0.000177191706852337,
    -0.000166775501198536, -0.000110390510136975, -8.38661319214754E-6,
    -1.63445704880961E-5, -5.70556982589265E-6, -5.07318377841507E-6,
    -6.70316461704801E-6, -6.6484582660322E-6, -5.73654046429135E-5,
    -0.000142705586531326, -0.0001436530193905, -0.000155207272650468,
    -0.000189881797781869, -0.000184246989085916, -0.000196767025140944,
    -0.00011152103097824, -1.06037947588458E-5, -1.40551047069593E-5,
    -8.86583148907711E-6, -8.33596390340824E-6, -9.91631620135037E-6,
    -6.49169614092165E-6, -0.000113704915557407, -0.000125128274028103,
    -0.000204496802248263, -0.000200017608805741, -0.00017822566198763,
    -4.99534489073849E-5, -1.20458128621537E-5, -1.15670737742211E-5,
    -8.72472268083471E-6, -1.83191026021962E-5, -1.18586230645515E-5,
    -1.43328928199102E-5, -0.000157613019576642, -0.000161913058503873,
    -0.000161289792669172, -0.00010033307977965, -3.5599860706493E-5,
    -7.2615382229817E-6, -1.23031683796063E-5, -1.32319208159505E-5,
    -5.53261782483923E-6, -6.15032190286214E-6, -8.58884722817049E-6,
    -9.05279171934421E-6, -0.000103658615408108, -0.000128807168561671,
    -0.00012387572765566, -0.000123573636290012, -8.28917437811648E-5,
    -4.5007500915339E-6, -4.38634272776751E-6, -6.14076751355716E-6,
    -8.15412453941435E-6, -0.000177765637556398, -0.000123479515393462,
    -0.000120239784402895, -0.000118823624753994, -4.88371210183146E-5,
    -9.18968783588726E-6, -9.03213727918035E-6, -8.84813363632271E-6,
    -9.08726160054474E-6, -1.11549108502938E-5, -5.09368950176721E-5,
    -0.000112693265829226, -0.000106760884851536, -0.000106250481579905,
    -5.51010467013737E-5, -1.1046537265038E-5, -1.49440623578564E-5,
    -7.71869491163913E-6, -7.59387948885352E-6, -9.0701411975117E-6,
    -8.09484784388537E-5, -0.000139269622245629, -0.00015294850922656,
    -0.000160341916017552, -8.82841821570198E-5, -6.41366139107793E-5,
    -9.15450478037172E-6, -9.18122989131312E-6, -8.23749571549455E-6,
    -7.58888031463243E-6, -9.38419464374454E-5, -9.07694183771614E-5,
    -0.000134872635367845, -0.000129966273582087, -0.000128678753050767,
    -8.0245907259936E-5, -1.21914798586079E-5, -1.13442089462156E-5,
    -9.80241605085448E-6, -9.16175783513117E-6, -9.07475826233103E-6,
    -1.04145555583366E-5, -9.80641887216871E-5, -0.000157323500725904,
    -0.000150943610814883, -0.000150752613785059, -0.000139824643174416,
    -9.46135266714556E-5, -1.03650063187625E-5, -1.06344234388798E-5,
    -1.09993840827385E-5, -9.80918634544873E-6, -1.18216680747354E-5,
    -0.000127442790034481, -0.000103020941096374, -0.000184408669789255,
    -0.000184548955903519, -0.000174488834048393, -9.77754669302826E-5,
    -1.75405526699597E-5, -8.97197841651975E-6, -4.84699017935973E-6,
    -9.07011398258407E-6, -1.90280489152836E-5, -0.000117312428857652,
    -0.000148327079475084, -0.000190396709108275, -0.000188985710755011,
    -0.000171263204232439, -7.79190414806197E-5, -1.69773283579407E-5,
    -1.29230602042652E-5, -7.1145255282087E-6, -6.89895365429885E-6,
    -2.7876894444875E-6, -5.9192475845926E-5, -0.000141372685811291,
    -0.000174229045695759, -0.000193491720150584, -0.000194147834293775,
    -0.000167880246329518, -5.64114978356616E-5, -9.51511183524426E-6,
    -6.11566684684795E-6, -4.63651990908311E-6, -1.03362111034725E-5,
    -6.74129927998092E-6, -7.0224135333394E-5, -0.000110077593346609,
    -0.000152558378665613, -0.000164304710836649, -0.000140868281608613,
    -8.3286360472578E-5, -4.55260467796109E-5, -1.5514995817799E-5,
    -9.47529427449462E-6, -7.2767840781641E-6, -6.24187597296297E-6,
    -2.17178144005681E-5, -2.79592037778289E-5, -0.000116398501518218,
    -0.000223736967315427, -0.000226877911672668, -0.000236193022482179,
    -7.10720245275073E-5, -7.82122290728334E-5, -1.6736077978166E-5,
    -9.82698368873746E-6, -9.49592401240821E-6, -1.05080166080649E-5,
    -7.93237821589824E-6, -7.42879969192932E-5, -0.000138645268130188,
    -0.00016727643646666, -0.000170007828419646, -0.000155129743912739,
    -0.000175620983410261, -3.60181026186464E-5, -3.12551842730173E-5,
    -8.56406179029016E-6, -4.39273968608281E-6, -4.49907570844699E-6,
    -9.36981862088738E-6, -2.07728374223961E-5, -0.000104107073492174,
    -0.000186949316826489, -0.000195954182807205, -0.000204786809403208,
    -0.00014647046312111, -3.87419629214469E-5, -7.04778462713218E-6,
    -9.00781899105237E-6, -6.31480000037082E-6, -5.97457629107038E-6,
    -4.55131408435362E-6, -2.90060606221672E-5, -0.000134758628631468,
    -8.88848664649415E-6, -9.66762712412739E-6, -9.10621526974487E-6,
    -7.72183872356522E-6, -8.17228613323658E-6, -6.61228453209804E-6,
    -7.00651913099572E-6, -1.10997532517877E-5, -4.48458503750047E-6,
    -1.11290369165819E-5, -1.05065490132372E-5, -5.9392562337396E-6,
    -7.33064464331615E-6, -6.48276177117091E-6, -6.75815402261873E-6,
    -7.93194056369333E-6, -6.10440114700798E-6, -5.91595487849781E-6,
    -1.03419533086625E-5, -6.36885768199041E-6, -4.4587131470596E-6,
    -5.16725743729813E-6, -5.33076469308502E-6, -2.4728706047477E-5,
    -1.69016467118379E-5, -1.82586647576243E-5, -1.67622820516346E-5,
    -4.04199150436132E-5, -0.000126001233157507, -0.000119685310294842,
    -0.000118476788292662, -8.35877085831509E-5, -1.49131226298595E-5,
    -7.53866001438795E-5, -2.43240680630421E-5, -7.011886030424E-6,
    -5.19594422790116E-6, -2.19509687165127E-5, -9.06934999467458E-6,
    -9.85770911194175E-6, -1.20654336575566E-5, -7.83915621840403E-6,
    -7.06430369024077E-6, -7.98504420287804E-6, -5.86468562658409E-6,
    -5.50399379682946E-6, -7.32184975772988E-6, -5.32339268940002E-5,
    -0.000112264463776969, -0.000104832727739611, -0.000104135884963682,
    -9.70982116745427E-5, -1.21163937372525E-5, -2.05293720983417E-5,
    -1.40725175912413E-5, -8.96188369635217E-6, -8.46901173003073E-6,
    -1.51616350726867E-5, -9.32878764393974E-6, -0.00010769482474605,
    -9.9435591388887E-5, -9.8346695552415E-5, -8.99858869237922E-5,
    -4.68726773393452E-5, -9.97553868702592E-6, -1.21129144062506E-5,
    -6.96848894577882E-6, -6.37590275587772E-6, -1.5035836467368E-5,
    -1.73264568247804E-5, -1.09454309698979E-5, -6.62346814153326E-5,
    -6.34427627941546E-5, -5.98868968661015E-5, -6.54640079222726E-5,
    -1.42436569615283E-5, -1.0430957200357E-5, -1.02553955460948E-5,
    -6.9605765019099E-6, -1.18490784587516E-5, -6.94058646108589E-6,
    -6.98734384208268E-6, -1.16309113494351E-5, -0.000106466892588552,
    -0.000214321150940487, -0.000213178757575661, -0.000220663665730863,
    -0.000134369501959848, -8.77351602837678E-6, -6.96175184958102E-6,
    -7.00978025516524E-6, -7.13674659810867E-6, -7.95688407022563E-6,
    -1.19683614641229E-5, -2.47264115336503E-5, -0.000101278437241446,
    -0.000142678911749731, -0.000137926231476783, -0.000145723023191903,
    -0.00010819194659338, -1.19918904038833E-5, -7.34147749845107E-6,
    -1.03725059467566E-5, -9.01042971663008E-6, -9.4875772177575E-6,
    -8.0315705548071E-6, -9.4911739143774E-5, -8.93553357763155E-5,
    -0.000111246451139282, -0.000108631072154566, -0.00013445804399219,
    -0.000124399346832737, -1.84444225484239E-5, -1.10901134192176E-5,
    -1.02614563412836E-5, -1.46354451933101E-5, -9.17084006264266E-5,
    -0.000184673361220414, -0.000159450778740429, -0.000159527920163272,
    -0.000130214604156413, -1.07528035859931E-5, -1.20287252628066E-5,
    -1.20212249510578E-5, -1.35007676087732E-5, -9.13951952470167E-6,
    -0.000105345396329678, -8.77500284456775E-5, -8.52032382279723E-5,
    -7.55873436418766E-5, -1.98536401016967E-5, -1.32807659211731E-5,
    -1.06843716060577E-5, -4.00433112347339E-5, -0.000102227242942546,
    -9.47247404378414E-5, -9.38929045417912E-5, -8.22598660177022E-5,
    -7.92905982076316E-6, -8.12385532891306E-6, -8.34710671815879E-6,
    -7.8501535168606E-5, -7.29370417826699E-5, -7.38012698470006E-5,
    -6.20471809409293E-5, -9.71844375499481E-6, -1.10389241468684E-5,
    -1.07598212234009E-5, -1.17509874127541E-5, -1.04945289152084E-5,
    -6.00682717668959E-5, -9.53636553894692E-5, -0.000189562117153352,
    -0.00018427214429803, -0.000164942032276026, -0.000122681350392739,
    -1.34554060110225E-5, -8.40658601096132E-6, -1.53239310630455E-5,
    -1.59884111706281E-5, -6.09343778027071E-5, -8.2508907705765E-5,
    -0.000175916445244115, -0.000170111318981523, -0.000151891735510357,
    -0.000114801347397399, -1.6001063684143E-5, -1.03019158733504E-5,
    -9.42099489718723E-6, -1.13630952430203E-5, -1.00769720680715E-5,
    -1.52371912063014E-5, -8.38587075972128E-5, -0.000123389436560836,
    -0.000116977260869652, -0.000122860255267328, -9.4781459434386E-5,
    -1.25445967593553E-5, -7.97703143861707E-6, -1.19467934976958E-5,
    -1.19570386707541E-5, -1.24633330408001E-5, -0.000131941471829384,
    -0.000220962571735873, -0.000209699996557113, -0.00021032405516764,
    -0.000160466802440738, -7.63348675582132E-5, -1.86561511627317E-5,
    -1.14406588504485E-5, -1.35885812533652E-5, -0.000136882016855614,
    -0.00028685291681349, -0.000289272691528409, -0.000289536091360151,
    -0.000187475220762035, -9.68817640006474E-5, -1.41521590195203E-5,
    -2.4032052140654E-5, -1.3249629218567E-5, -0.00010617292818913,
    -0.000109997107251859, -0.000167537755874944, -0.000160462921815358,
    -0.000158955771402222, -0.000105847879657498, -1.01126495108705E-5,
    -1.39308304825036E-5, -7.08473829437018E-6, -7.02295026587967E-6,
    -6.70639922672723E-6, -7.84260803569458E-6, -9.25126657151726E-6,
    -0.000135070876898362, -9.78487644247501E-5, -9.3341101644596E-5,
    -9.02703216493307E-5, -1.94994177843486E-5, -1.66230881365404E-5,
    -2.63880613033673E-5, -1.07931616141554E-5, -0.000134175695324004,
    -9.77662651259379E-5, -9.50350847182784E-5, -9.3493876991659E-5,
    -6.01199558296694E-5, -7.02324429716871E-6, -7.44985104951059E-6,
    -6.99222553212976E-6, -7.04757344774036E-6, -1.24976348082499E-5,
    -8.09609927600924E-6, -1.47211790036485E-5, -1.56710535279468E-5,
    -1.50160665415126E-5, -1.43268591662523E-5, -1.50981662343124E-5,
    -1.23417353871628E-5, -1.20069952237878E-5, -1.05403439153204E-5,
    -1.88372090568534E-5, -1.0251454464012E-5, -5.04144921822426E-6,
    -3.98843050355779E-6, -4.40767119535042E-6, -8.9462964456696E-6,
    -1.07910903823339E-5, -1.73182676952325E-5, -2.61631424144494E-5,
    -2.82713415149332E-5, -3.38294323661282E-5, -3.66006014401637E-5,
    -3.92743200807778E-5, -3.97467033455829E-5, -4.54767145477611E-5,
    -8.05928376403149E-6, -0.000119776201416713, -0.000245197314019147,
    -0.000243872368352079, -0.000240838521021773, -0.000255646694395525,
    -0.000120552750509985, -8.00521241964021E-6, -9.11894011173878E-6,
    -6.47460927269326E-6, -1.02031170895455E-5, -9.93887183513E-6,
    -5.98988728947729E-6, -1.11917435852053E-5, -7.4669389144203E-6,
    -7.74241719846991E-5, -0.000175701477130939, -0.000236094679260908,
    -0.000230723593418985, -0.000223321053229229, -9.25418572083832E-5,
    -1.92205665988439E-5, -1.66911687793397E-5, -9.57623485041263E-6,
    -1.04519469236287E-5, -9.99072621088906E-6, -7.85563622809201E-6,
    -5.98844698608733E-6, -1.17022860684803E-5, -5.10511276877264E-5,
    -7.52500833988002E-5, -0.000137408143675861, -0.00013345591155873,
    -0.000114717906411846, -5.31937089480141E-5, -9.06010198629713E-6,
    -8.87753064011961E-6, -8.47276302990229E-6, -5.97054212344256E-6,
    -6.48979365444272E-6, -1.7361760824415E-5, -9.72764326723978E-6,
    -7.54316965776224E-6, -0.00012014649799911, -0.000106597591301689,
    -0.000104355613608076, -0.000107948680532507, -2.18637621865687E-5,
    -9.84774175909065E-6, -5.78490912254631E-6, -8.97622692851363E-6,
    -4.66065870010025E-6, -0.00010009431996042, -0.000106792234111124,
    -0.000101818882636135, -0.000101661396736081, -9.23719073806131E-5,
    -8.10144734159555E-6, -5.32981604832421E-6, -9.01667539043E-6,
    -8.34840169450386E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 21.0, 23.0, 15.0, 17.0,
    0.0, 0.0, 41.0, 67.0, 79.0, 80.0, 36.0, 0.0, 0.0, 17.0, 0.0, 13.0, 0.0, 0.0,
    0.0, 34.0, 69.0, 88.0, 73.0, 39.0, 24.0, 28.0, 12.0, 26.0, 36.0, 9.0, 41.0,
    71.0, 78.0, 66.0, 39.0, 2.0, 0.0, 0.0, 0.0, 16.0, 0.0, 6.0, 6.0, 16.0, 45.0,
    63.0, 81.0, 57.0, 63.0, 35.0, 10.0, 3.0, 6.0, 0.0, 0.0, 4.0, 16.0, 51.0,
    70.0, 85.0, 98.0, 55.0, 42.0, 15.0, 0.0, 13.0, 0.0, 0.0, 17.0, 40.0, 56.0,
    89.0, 78.0, 48.0, 26.0, 22.0, 4.0, 13.0, 0.0, 0.0, 0.0, 5.0, 0.0, 33.0, 79.0,
    86.0, 109.0, 87.0, 42.0, 24.0, 15.0, 13.0, 3.0, 2.0, 15.0, 62.0, 67.0, 93.0,
    101.0, 95.0, 112.0, 86.0, 64.0, 38.0, 0.0, 0.0, 16.0, 0.0, 0.0, 14.0, 38.0,
    62.0, 89.0, 67.0, 68.0, 22.0, 15.0, 0.0, 0.0, 0.0, 0.0, 29.0, 38.0, 65.0,
    66.0, 58.0, 33.0, 14.0, 0.0, 0.0, 19.0, 0.0, 0.0, 23.0, 52.0, 83.0, 84.0,
    58.0, 15.0, 8.0, 0.0, 10.0, 0.0, 0.0, 6.0, 44.0, 80.0, 102.0, 83.0, 102.0,
    39.0, 8.0, 0.0, 0.0, 25.0, 0.0, 0.0, 0.0, 37.0, 89.0, 100.0, 98.0, 77.0,
    54.0, 0.0, 21.0, 0.0, 0.0, 58.0, 67.0, 90.0, 119.0, 106.0, 62.0, 22.0, 0.0,
    0.0, 0.0, 61.0, 62.0, 61.0, 93.0, 82.0, 55.0, 41.0, 0.0, 14.0, 0.0, 0.0, 0.0,
    0.0, 42.0, 44.0, 61.0, 77.0, 80.0, 56.0, 51.0, 36.0, 6.0, 12.0, 0.0, 0.0,
    0.0, 0.0, 27.0, 60.0, 93.0, 112.0, 64.0, 47.0, 12.0, 6.0, 0.0, 13.0, 10.0,
    10.0, 35.0, 81.0, 70.0, 58.0, 42.0, 0.0, 10.0, 14.0, 0.0, 0.0, 0.0, 0.0,
    27.0, 66.0, 83.0, 85.0, 46.0, 0.0, 0.0, 0.0, 0.0, 18.0, 66.0, 82.0, 83.0,
    31.0, 0.0, 0.0, 0.0, 0.0, 8.0, 35.0, 69.0, 85.0, 85.0, 49.0, 8.0, 10.0, 0.0,
    0.0, 0.0, 62.0, 68.0, 95.0, 69.0, 66.0, 39.0, 0.0, 0.0, 0.0, 0.0, 57.0, 98.0,
    122.0, 113.0, 87.0, 50.0, 10.0, 8.0, 0.0, 0.0, 0.0, 3.0, 59.0, 72.0, 108.0,
    86.0, 69.0, 45.0, 5.0, 6.0, 8.0, 0.0, 9.0, 34.0, 71.0, 92.0, 100.0, 66.0,
    55.0, 18.0, 0.0, 0.0, 0.0, 25.0, 30.0, 65.0, 91.0, 95.0, 77.0, 61.0, 22.0,
    13.0, 0.0, 0.0, 0.0, 43.0, 57.0, 59.0, 90.0, 69.0, 58.0, 49.0, 0.0, 0.0, 0.0,
    3.0, 0.0, 56.0, 66.0, 90.0, 109.0, 91.0, 71.0, 46.0, 10.0, 0.0, 0.0, 0.0,
    15.0, 16.0, 65.0, 69.0, 98.0, 82.0, 70.0, 57.0, 16.0, 0.0, 0.0, 5.0, 0.0,
    45.0, 58.0, 88.0, 111.0, 93.0, 73.0, 47.0, 19.0, 0.0, 0.0, 0.0, 0.0, 25.0,
    71.0, 71.0, 105.0, 70.0, 59.0, 48.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 60.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 25.0, 29.0, 36.0, 36.0, 36.0, 75.0,
    78.0, 89.0, 44.0, 10.0, 32.0, 28.0, 0.0, 0.0, 15.0, 0.0, 0.0, 9.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 48.0, 61.0, 79.0, 84.0, 43.0, 8.0, 39.0, 10.0, 0.0, 0.0,
    12.0, 0.0, 46.0, 93.0, 96.0, 67.0, 30.0, 0.0, 9.0, 0.0, 0.0, 11.0, 18.0, 5.0,
    68.0, 83.0, 94.0, 75.0, 13.0, 3.0, 4.0, 0.0, 7.0, 0.0, 0.0, 11.0, 59.0, 63.0,
    77.0, 44.0, 37.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 15.0, 70.0, 76.0, 104.0,
    56.0, 51.0, 15.0, 0.0, 4.0, 0.0, 0.0, 0.0, 32.0, 96.0, 96.0, 104.0, 74.0,
    24.0, 16.0, 5.0, 3.0, 11.0, 64.0, 57.0, 88.0, 71.0, 54.0, 14.0, 8.0, 13.0,
    10.0, 0.0, 51.0, 83.0, 86.0, 69.0, 15.0, 22.0, 4.0, 18.0, 82.0, 87.0, 89.0,
    43.0, 0.0, 0.0, 0.0, 71.0, 88.0, 94.0, 58.0, 0.0, 6.0, 5.0, 18.0, 9.0, 41.0,
    87.0, 63.0, 81.0, 47.0, 43.0, 13.0, 0.0, 32.0, 15.0, 37.0, 98.0, 80.0, 101.0,
    55.0, 52.0, 19.0, 4.0, 0.0, 7.0, 2.0, 9.0, 73.0, 90.0, 113.0, 75.0, 62.0,
    12.0, 0.0, 8.0, 8.0, 15.0, 43.0, 47.0, 80.0, 75.0, 80.0, 59.0, 18.0, 9.0,
    20.0, 56.0, 57.0, 77.0, 67.0, 35.0, 41.0, 13.0, 29.0, 28.0, 78.0, 87.0,
    105.0, 89.0, 60.0, 32.0, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 28.0, 81.0,
    107.0, 107.0, 72.0, 16.0, 30.0, 10.0, 35.0, 73.0, 80.0, 80.0, 29.0, 0.0, 0.0,
    0.0, 0.0, 18.0, 0.0, 13.0, 28.0, 46.0, 32.0, 19.0, 17.0, 12.0, 4.0, 14.0,
    4.0, 0.0, 0.0, 0.0, 0.0, 15.0, 46.0, 60.0, 65.0, 73.0, 85.0, 88.0, 77.0,
    39.0, 0.0, 24.0, 46.0, 82.0, 86.0, 56.0, 41.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0,
    8.0, 0.0, 38.0, 50.0, 83.0, 92.0, 66.0, 48.0, 19.0, 31.0, 0.0, 4.0, 0.0, 0.0,
    0.0, 8.0, 33.0, 83.0, 100.0, 113.0, 67.0, 41.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    12.0, 0.0, 0.0, 51.0, 75.0, 93.0, 74.0, 49.0, 0.0, 0.0, 0.0, 0.0, 28.0, 73.0,
    81.0, 81.0, 43.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 11.0, 14.0, 12.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0,
    12.0, 11.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0,
    32.0, 32.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0,
    25.0, 24.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 26.0, 31.0, 35.0,
    36.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 29.0, 31.0,
    25.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 25.0, 23.0, 23.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 21.0, 21.0, 21.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 23.0, 26.0, 10.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 28.0, 32.0, 32.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    9.0, 0.0, 29.0, 27.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 26.0, 24.0,
    11.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 34.0, 40.0, 24.0, 17.0, 8.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 17.0, 32.0, 31.0, 36.0, 43.0, 9.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 21.0, 26.0, 26.0, 41.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 19.0, 19.0, 19.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    5.0, 13.0, 11.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0, 13.0, 13.0, 11.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 6.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 21.0, 23.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 14.0,
    12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 21.0, 24.0, 16.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 5.0, 24.0, 34.0, 30.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 10.0, 22.0, 30.0, 30.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    17.0, 42.0, 46.0, 33.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 28.0,
    30.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 23.0, 23.0, 24.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 18.0, 22.0, 23.0, 18.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 25.0, 30.0, 31.0, 13.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 15.0, 13.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 5.0, 5.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 28.0,
    37.0, 37.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 33.0, 33.0, 33.0,
    7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 19.0, 17.0, 11.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 31.0, 31.0, 31.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 38.0, 43.0, 47.0, 12.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0, 35.0, 42.0, 12.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 14.0, 13.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 42.0,
    44.0, 35.0, 17.0, 0.0, 0.0, 0.0, 0.0, 15.0, 26.0, 38.0, 39.0, 31.0, 0.0, 0.0,
    0.0, 0.0, 7.0, 12.0, 26.0, 36.0, 34.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0,
    24.0, 23.0, 23.0, 28.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    25.0, 24.0, 24.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 20.0, 20.0, 33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0,
    7.0, 6.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 6.0, 3.0, 3.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    20.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 23.0,
    24.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 22.0, 22.0, 22.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 28.0, 47.0, 47.0, 50.0, 25.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 40.0, 41.0, 22.0, 3.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 37.0, 50.0, 51.0, 29.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 47.0, 49.0, 27.0, 31.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 38.0, 38.0, 47.0, 41.0, 19.0, 18.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 26.0, 27.0, 27.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 29.0, 32.0, 32.0, 33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 22.0, 23.0, 23.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 40.0,
    46.0, 46.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 20.0, 47.0,
    48.0, 26.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 41.0, 43.0, 23.0, 22.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 19.0, 34.0, 58.0, 51.0, 35.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 21.0, 26.0, 36.0, 52.0, 31.0, 31.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 7.0, 21.0, 44.0, 46.0, 27.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    21.0, 27.0, 28.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 20.0,
    21.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 20.0, 20.0, 20.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 18.0, 20.0, 20.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0,
    36.0, 41.0, 43.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 18.0, 41.0, 43.0,
    21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 33.0, 43.0, 44.0, 23.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 19.0, 52.0, 55.0, 33.0, 29.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 10.0, 30.0, 50.0, 53.0, 31.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    22.0, 44.0, 54.0, 54.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 44.0,
    48.0, 30.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 44.0, 44.0, 49.0,
    24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 34.0, 39.0, 34.0, 22.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 43.0, 53.0, 56.0, 32.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 24.0, 25.0, 25.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 24.0, 24.0, 23.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 22.0, 22.0, 22.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 11.0, 16.0, 16.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 24.0, 49.0, 51.0, 53.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    20.0, 43.0, 44.0, 46.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 37.0,
    38.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 39.0, 44.0, 44.0, 23.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 14.0, 18.0, 19.0, 18.0, 0.0, 0.0, 0.0, 0.0, 20.0, 20.0, 21.0,
    20.0, 0.0, 0.0, 0.0, 16.0, 16.0, 17.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    20.0, 50.0, 51.0, 32.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 45.0, 46.0,
    35.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 36.0, 37.0, 40.0, 20.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 23.0, 46.0, 52.0, 53.0, 27.0, 0.0, 0.0, 0.0, 0.0, 19.0,
    45.0, 55.0, 55.0, 33.0, 0.0, 0.0, 0.0, 0.0, 15.0, 26.0, 46.0, 47.0, 26.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 18.0, 19.0, 19.0, 0.0, 0.0,
    0.0, 0.0, 11.0, 19.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 9.0, 11.0, 11.0, 0.0, 0.0, 9.0, 35.0, 44.0, 44.0, 27.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 42.0, 42.0, 24.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 37.0, 38.0, 23.0, 16.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 21.0, 21.0, 21.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 5.0, 19.0, 19.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    4.68941801192126E-7, -4.47520889997594E-5, 9.26531952805842E-5,
    2.89024877067483E-5, -5.65347706349665E-5, -6.17421710958912E-5,
    5.72292592010413E-5, 9.14012997597007E-5, -6.57554108071996E-5,
    -1.79525734332504E-5, -6.16315032723875E-5, 1.06149692900965E-5,
    3.60906591152276E-5, -1.23618020628261E-5, -6.42055126721976E-5,
    0.000134602252039803, -6.5171878399102E-5, -4.13345496177953E-6,
    -0.000122040799237274, -6.27629393369019E-5, 4.88530132341822E-5,
    9.73473434430994E-5, 8.31817937262999E-5, 6.33115353789544E-5,
    -0.000328248893054477, -2.42946777218574E-5, -1.99969026638784E-5,
    -0.000132943964539858, -4.48039472817843E-5, -2.04888602208624E-5,
    -0.000451769006612743, -0.00032341444515133, 7.40265613433918E-5,
    5.99941632424491E-7, -0.000250532548082012, 0.000187487611189527,
    -0.000189593111221929, 0.000130399376520488, 0.000329434287165668,
    -1.66655259033594E-5, 4.35743599432142E-5, -9.55153121195523E-5,
    -0.000203392027741188, -1.64182383899623E-5, -0.00105782541035405,
    -0.000378669162701965, -0.000339269886130007, 0.000157628977493603,
    -0.000311920239359903, -0.000228795902709589, 0.000153550433005536,
    -0.000296709920183479, -0.000174222173386368, -0.000108135038088201,
    -0.000275612836140447, -0.00120095146045776, 5.02510524787133E-5,
    0.000306009099695254, 7.47654695456179E-6, 7.40899560330009E-5,
    -0.000187704253605732, -0.000106927852228361, 0.000242247753879295,
    1.2749080149719E-5, 2.03284969386499E-5, 0.000293358930148539,
    0.000108133684421675, -0.000503442274106349, 0.00276065911070455,
    0.00297125313904727, -0.000269550923796835, 3.70073277786773E-5,
    -4.85747275149755E-5, -0.000106128318528738, 0.00026226236713001,
    -0.000158487495361269, -0.000114717328279972, -0.000194685274376856,
    0.000148786225634442, 2.18062814143358E-5, 0.000216385493215861,
    -0.000604116555136793, -0.000923391026720637, -6.37362920298888E-5,
    -0.000163780488069602, -2.6611700346526E-5, -0.00017938527575009,
    -0.000103981984716734, 0.000285950799594575, -0.000164960381755693,
    0.000117705420806098, 7.21354711696419E-5, -0.000368828034480384,
    0.000273943141460437, -0.00127002221249867, -0.00314234677778485,
    -0.000354789757420769, 0.000268156530796646, -0.000122340188712194,
    0.000219963690153347, -3.42445913633844E-5, -8.4168417307605E-5,
    -2.26578586822362E-5, -2.28861254238982E-5, 0.000420013572032383,
    2.27900082471464E-5, 0.00027538336068384, -0.000493346753397949,
    -0.000734424452417519, 0.00352822310695463, -0.000239805070537469,
    -0.000257600755729144, -0.000235621357362965, -0.000117849301666226,
    0.000124365706826271, -0.000219036313889358, -0.000178589013986116,
    -7.16475488799993E-6, -7.09505439121573E-5, -0.000492650463840751,
    0.00117651676049044, 0.00263993092023711, -0.000152279772707651,
    0.000720129804026525, -0.00095711708898327, -3.37363014094619E-5,
    -6.87060046917378E-5, -0.000196551712952328, 7.29960454329905E-5,
    -3.25311410009533E-5, -0.00019543267731095, 0.000127407122737846,
    -4.53486059875381E-5, 0.000134963645551456, -0.000110767888331843,
    -0.00160792289360275, 0.000819321970256896, -0.000238647441568062,
    0.000282554716931077, 0.000368540177470259, 7.25096919498153E-5,
    -0.000108157928432637, 9.12846367680967E-5, 0.000220163425900721,
    8.27143803833466E-5, -0.000717491167778807, 0.0033074133484154,
    -0.000418042124624598, -0.000250675915893097, -4.11531132560076E-5,
    6.85280924719484E-6, -0.000276757839389285, 2.55083597235917E-6,
    6.40623881100293E-5, -0.00012673123222291, 0.000248638242005523,
    0.000156997115760359, -0.000441240221664005, 0.00311827901475958,
    0.000323260229373212, -0.000117594293394437, -2.2669877617147E-5,
    -0.000125668421148087, -1.27270584891606E-5, -3.33948842393938E-5,
    0.000270510204278634, 0.000392865152047128, -0.000132207319247568,
    -0.000272472741585842, 0.000754526451911103, -7.3020005942351E-5,
    -0.000295332034305219, -1.94216807413586E-5, -2.73852585423854E-5,
    0.000155484066621936, -0.000422120942206268, -6.6611607056804E-5,
    -4.96613106585621E-5, 6.16635120285066E-5, 0.000179328145890317,
    -0.000157140778032979, 4.39661499908365E-7, -7.18356564718619E-6,
    -0.00129174927567381, 0.000503550250702193, 8.85567851571596E-5,
    -0.000365055269617655, -4.19485350990292E-5, -3.6092670888902E-5,
    -0.00011485083493004, -2.70560269179202E-5, 0.000437834523952433,
    -0.000408286698919442, -0.000605224348565255, 0.000368981238270964,
    -0.000114558581670228, 0.000208036666589239, 4.32600028613926E-5,
    -3.65571108450642E-5, -6.73197412223731E-5, 5.82767250822228E-5,
    0.000327769171425686, -0.000503275291491759, 0.000462803360007771,
    0.000159004950008473, -0.00209982821582186, -0.00012919172048751,
    -0.000166339053330558, 5.85369084817973E-5, 3.93820125223685E-5,
    -6.78469863564424E-5, -6.47366793224519E-5, -8.88619126599427E-5,
    0.000148704467495678, -3.74742384588617E-5, -0.000540229358182179,
    0.00207574217921535, -0.00022264731237661, -0.00175693918367932,
    -0.00137532612064172, -5.15641497676031E-5, -0.00018852330826162,
    0.000119149576969543, -6.49326729608811E-5, -5.59525687322679E-5,
    -2.66899687199667E-5, -0.000149925019649744, -0.000288358586107309,
    9.03714407344019E-5, -0.000421540985647363, -0.00289651771473142,
    0.000645076718136908, 0.000115169113322664, -8.08958314580895E-5,
    -8.94888310763245E-5, -0.000234639802690413, 0.000418736278474952,
    -2.26877389375728E-5, 0.000223664965869309, -0.000197186299961408,
    -0.000353128690121778, 0.0024296520575011, -0.000697547793649279,
    -1.24534122324577E-5, -0.000139453385341166, -0.000206444854299617,
    0.000115718325140615, -4.65611306503721E-5, 9.08867280662389E-5,
    2.07764477399576E-5, 3.20696075823701E-5, -9.01435484680505E-5,
    1.84106696013576E-5, -0.000575281502610568, -0.00179411642255289,
    1.55949303761195E-5, -6.33876520630389E-5, -1.52700777659318E-5,
    -1.81962955980185E-5, -0.000134158164977549, -0.000139735724774658,
    6.01774334001349E-5, -0.000476791596414418, -0.000618685519573416,
    0.000170181360663934, 0.000243131809312348, -5.11486588535495E-5,
    -5.95420592576925E-6, 0.000206848109220938, -7.1626816914734E-5,
    6.68350682632735E-5, 0.00013049034241873, -0.000396621416485141,
    0.0015433051947039, -9.8299584320565E-5, -0.000165068457337302,
    5.12992954753627E-6, -1.00191559452127E-6, 0.000144236806819069,
    0.000157070840563658, 0.000170377198581386, -0.000196042860333643,
    0.00131317938193165, 0.00117102036815193, -0.00117315611994805,
    -0.000103322083895893, -9.38490859909085E-5, 5.40753463218356E-5,
    0.000111154057539205, 0.000187746261210808, -6.79316860080723E-5,
    -0.000249920595441214, 0.000589432401730316, -0.000847374326313941,
    -0.001651297959371, -0.000141183006716675, -0.000164264743102582,
    0.000127446672624103, -0.000117923360316198, 3.06042936477954E-5,
    0.000156283197898147, 7.34108734302145E-5, -1.91395798278088E-5,
    -0.000108235179894576, 0.00147398504602809, -0.000669778960476077,
    -0.00185594542322475, -2.99032765313878E-5, 2.49526093480546E-5,
    -5.86439428482879E-5, 0.000217332215945855, -0.000112340460982966,
    -9.98705158193614E-5, -0.000124326000147707, 0.000266381419394045,
    -0.000288300302619933, -0.00108918524513002, 0.00243203962692672,
    -9.52725086252156E-5, 0.000145537673650656, -2.74520647333761E-5,
    -0.000129127052175092, -0.000314815688789246, 0.000117518053380441,
    -6.29325444632407E-5, -0.000141480776828415, -0.000317360471638659,
    -0.00173049801734065, 0.00192036720005255, -0.00021008563912804,
    -4.42399852138152E-5, 9.68369779881528E-5, 9.66659495634183E-5,
    3.90570329499046E-5, -5.80873965106313E-5, 4.35110067492228E-5,
    0.000226588613919741, -0.000580341492911316, 0.0026334231237239,
    3.55268480899049E-5, -0.0009936809801404, -8.0077525774055E-5,
    -5.43450376287896E-5, 4.6573986909773E-5, -9.43919899576337E-5,
    5.43678026292614E-5, 7.04755186484432E-5, -1.2788989322195E-5,
    0.000175951817392615, -0.000501878796071592, -0.00155105467356554,
    0.00260611227240493, -0.00083781323249798, -0.000211825059882332,
    0.000172829046443882, 0.000125437856282936, 0.000117417885440644,
    -0.000277463929863093, -7.04701491963915E-5, -6.43581860078477E-5,
    0.000335200481930144, -5.92699939711276E-5, 0.00210881860096867,
    0.00198698478486069, -0.00130628851319717, -7.38793145563119E-5,
    0.000315106868036464, -3.3595457357394E-5, 0.000192137572766898,
    -0.000123581047246153, -5.37686126512168E-6, -0.000180896061427182,
    0.000362065799901334, -0.000387090813576246, 7.21581878890499E-5,
    -0.000717375337309478, -0.000711652562489547, -0.000407478825707658,
    0.000256550850225918, 0.000107985158295337, -2.73569051430595E-5,
    -0.000125095752409932, -0.000183237942750769, 5.60288574892103E-5,
    8.21017553613096E-5, -0.000480243284069348, 0.00272285275753677,
    0.000805551408512664, -0.000876082244764847, 2.57984615474135E-5,
    -0.000141807900254671, 0.00016092500254481, 0.000101049152276143,
    -4.08029415019531E-6, -3.85032449694062E-5, 0.000137367048044443,
    9.52603468007452E-6, -0.000150459858962436, 0.000158052417731176,
    -6.49432452787033E-5, 1.09081906292166E-6, 3.90519448506868E-5,
    -4.59592744581275E-5, -0.000176188265498516, -3.67647089324323E-5,
    -0.000141121141478148, -0.000120502121024817, 7.06623324945868E-5,
    4.49778171662749E-5, 7.81568350973655E-5, 2.20905724896623E-5,
    -5.26375233808804E-5, 1.63209871991769E-5, 2.92125711659395E-5,
    -0.000127155826509518, -7.69024447096207E-5, -4.53037902073279E-6,
    4.21073495829392E-5, 4.57571480302356E-5, 4.11755040565207E-5,
    8.47371316737943E-5, 0.000127969829680396, -8.24488707988234E-5,
    4.23399085801879E-5, -9.20458226237113E-5, -9.72989720761144E-5,
    -0.000553504086578028, 0.00251508334468853, 0.000617062524567034,
    -0.00107234868270509, -0.000383006108368999, -5.10338446373996E-5,
    3.59296070220254E-5, -0.000273998312436988, 4.82052539908757E-5,
    3.06833700184359E-6, 7.28275224884327E-5, 0.00027156298612847,
    2.27291510892733E-5, 0.000291295481440078, 0.000255384850589694,
    2.50242850778375E-5, -4.30339230004385E-6, -0.000218693001937336,
    -3.97512195434707E-5, 6.55282028923202E-5, -0.000321196171461838,
    0.00189398994123257, 7.83553395437744E-5, -0.000318540868084679,
    -0.000189798455034479, 1.94810993063173E-6, 3.69031847206475E-5,
    -0.000165389002228479, -0.00011139154270362, 5.01262620938345E-5,
    -0.000164110719525142, -0.000129473168114597, 0.000642927538967625,
    -0.000710458504492001, 0.000288022593928342, -8.68215747926866E-5,
    9.62101937764481E-5, 0.000171518113969588, -0.00018084981194626,
    -1.89261931402749E-5, 0.000223192693931495, 0.000260185807116011,
    0.000136311587077155, -8.14263951589235E-5, 0.000642959661316545,
    5.82709235598242E-5, -1.55064891860491E-6, 0.00038544517016325,
    8.39452348934916E-5, -0.000155224276843714, 6.93696317203058E-5,
    2.72966741508898E-5, -0.000148492161494204, 0.000194528636541459,
    0.000308021830533714, -0.000159004208968599, 0.00202661738458157,
    0.0024805536146669, 0.000187182446153363, -8.99840752097087E-5,
    0.000110376201375793, 6.906916482341E-5, 0.000196883957011701,
    -0.000147640581974238, -3.97953869920424E-5, -0.000116760942360435,
    9.66340015090154E-6, -0.000403980609273964, 0.00239750065895876,
    0.00160556317232445, -0.000205334843069109, 0.000242131267393751,
    3.48255932113791E-5, -0.000141591585381148, 0.00024584228391745,
    -5.03821011960081E-6, -7.50178973496976E-5, -2.68134606814004E-5,
    -0.000146034077740607, -0.000218288573806682, -0.00206047865152875,
    -0.00152319696997034, -0.000223017686134476, 0.000149878561599496,
    -6.55047622784313E-5, -0.000219962029337108, -6.31667890502034E-7,
    0.000187796747557117, -9.52316206106119E-5, 0.00051256612970811,
    -6.01231655390386E-5, -0.000208565674508949, 0.000301932965007424,
    7.57965758237215E-5, 0.000172328267929379, 0.000118086921404815,
    0.000189058603857151, -9.43777790825905E-5, -4.18252162376106E-5,
    0.00029529380938875, -0.000219601134365553, -0.000194595197884616,
    -9.08844038375001E-5, -0.000139725723501701, -0.000119238905293035,
    0.000156647355911637, -0.000439714753291358, 0.00186888731198835,
    0.000190169241984307, 9.23216530530984E-5, -4.64479571333379E-5,
    -2.09863982336705E-5, 0.000209129935990503, -8.92197399728918E-5,
    0.00168457694495842, -0.000101757445086858, 2.01718899633295E-5,
    0.000174265447782458, -2.78563302151498E-5, 6.57713270275463E-5,
    6.88869711210107E-5, -0.000224310633158445, 0.000186357262877436,
    -0.000240803707839816, -0.00145668210836423, 0.00255890409388853,
    -0.000264099486721156, -5.78317809436399E-5, -0.000196210710639618,
    0.000168092371846118, -0.000172946268298242, -4.26780415779591E-5,
    0.000349310720594424, -0.000301136865230547, -0.000944687409571683,
    0.00318344003612424, 0.000275360692137981, 0.000214535854525235,
    -4.06366976853139E-5, -0.000107351362640147, 4.29284713609145E-6,
    6.37602298258214E-5, -5.04720266385232E-5, 0.000164826073965974,
    -0.000142032056690597, 0.00168361314675674, 0.000535028662568885,
    -0.000145732917382276, 0.000213057896892891, 7.2603511480077E-6,
    0.000134984880768442, 2.55367318409854E-5, 0.000333288014407944,
    0.00020823744888523, -0.000160457354682064, 0.00160219779398887,
    -0.000346923848491622, -0.00129956814222784, 1.84971004995804E-5,
    -1.76664337674503E-5, -0.000154264666361138, 0.000369184280907131,
    0.000215163559584663, -0.000148990698175536, 0.000598888376723531,
    0.000593333657374474, -0.00154352799387609, -0.000113969608040541,
    -0.000284264808015407, -0.000185965272977229, -0.000115851865860314,
    2.06742867798154E-5, -0.000319366738175955, 0.000455109448453227,
    -0.00113964628357746, -0.00128983965905868, -8.48681983419496E-5,
    -0.00015432680886412, -6.50506913122872E-5, 0.000319624890806118,
    -9.78103819985648E-6, 0.00018499524816312, -0.000149844870093361,
    -6.38943022811028E-5, -0.000163639716357602, 2.41083313003435E-5,
    -0.00010250582792823, -0.000451290554824196, 0.000124990338879453,
    -0.000336628512993755, -0.000196508907832225, -0.000185689385586427,
    -0.00011562042588975, 8.07681129380619E-5, -0.000299847320795887,
    -0.000874499969449713, -1.88039307839833E-5, 0.000121347665409026,
    -2.12371464910018E-5, 0.000131840239790144, -0.000186276137785101,
    0.000278961423522441, -0.000221541426625495, 1.06269871005322E-5,
    -1.69274263321749E-5, -0.000206506721034909, 0.00014525746567906,
    0.000126617861892308, -0.000125384815121471, 9.35121747552548E-5,
    -0.000149244518728397, -0.000136952955480185, -5.04892594999974E-5,
    -1.25123900214882E-7, -2.90605522062075E-5, 6.97018687300297E-5,
    -6.43773418485511E-5, -2.14265273684975E-5, 5.16036097562677E-5,
    -0.00030490831613934, 1.31601956734885E-5, 0.000103148295527111,
    -0.000214459793250805, -7.52516860906172E-5, -0.000813404143848123,
    -0.000130193401762272, -6.04957081017617E-5, -5.72963950607197E-5,
    -3.22268175066348E-5, -0.000311677721928083, -0.00240755989288043,
    -0.00228634694763538, -4.96386573561953E-5, -5.39066384651254E-5,
    -1.69215404165407E-7, -9.79617507237069E-5, -0.000117846735727783,
    6.27956176059119E-5, 1.77016993051885E-5, 0.000202276417018978,
    -6.12088707822611E-5, -0.000191145264469984, 0.000204613593607701,
    -0.000419281075470489, -0.00170342687017127, -0.0010740999960372,
    -0.000272682276480384, -0.0002631469002715, 4.28293641502678E-5,
    3.92064532990767E-5, -7.15439968724529E-5, 0.00022140380106428,
    2.39509901915298E-5, -0.000163109978861303, -0.000243189113941892,
    -0.000202764079465496, -7.47433447211477E-5, -0.000278862453864805,
    -0.00203402111847207, 0.00112525715721016, -2.786563091121E-5,
    5.36130668999387E-5, -1.86968081421565E-5, -0.000111877365070131,
    -5.15696082704026E-5, -0.000183412047559974, -0.000100870218287374,
    -0.00012844860582449, 2.16283617861772E-5, 1.50137682289112E-5,
    -0.00020466057707876, 0.000207399083251492, -0.000500454379703544,
    -2.3624925479789E-5, -6.53742479008787E-5, -0.000184973390331862,
    2.79404444202174E-5, -0.000156577306225598, -9.53849126382047E-5,
    6.30663647602874E-5, -0.000384767171324047, -0.0017195730777127,
    -8.43349940527921E-5, 2.74517174975806E-5, -5.60266341309705E-5,
    7.7368970430019E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 2.0, 2.0, 3.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
    2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 3.0, 2.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 65.0, 65.0, 67.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 62.0, 62.0, 62.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 71.0, 75.0,
    75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 86.0, 87.0,
    81.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 68.0, 68.0, 225.0, 160.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    81.0, 82.0, 82.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 78.0,
    89.0, 88.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 93.0, 93.0, 93.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 89.0, 106.0, 104.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 167.0, 175.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 75.0, 75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 80.0, 80.0, 80.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 89.0, 92.0, 92.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 83.0, 85.0, 85.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 92.0, 94.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 72.0, 93.0, 93.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 90.0, 94.0, 95.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 89.0, 90.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 75.0, 75.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    97.0, 99.0, 97.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 82.0, 92.0,
    94.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 73.0, 74.0, 74.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 76.0, 76.0, 77.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 73.0,
    73.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 71.0, 71.0, 71.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 80.0,
    80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 77.0, 77.0, 77.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 78.0, 78.0, 79.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 75.0, 76.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 72.0, 72.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 86.0, 91.0, 91.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 57.0, 59.0, 59.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    57.0, 63.0, 63.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    72.0, 72.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.38609430295611E-8, -1.90065368002879E-7, 4.4499870093796E-7,
    -6.10366317101492E-8, -2.0361316729231E-7, -2.89918941041277E-7,
    4.75976584602191E-7, 3.60775009760028E-7, -2.44106248686244E-7,
    -2.26570425080015E-7, -6.65482222156752E-8, 7.35054750663825E-8,
    2.45836352367109E-7, -2.19577683930374E-7, -2.16062061040051E-7,
    5.89192753960608E-7, -1.55526958556169E-7, 2.28996856804694E-7,
    -3.24000046152895E-7, -2.53940415276264E-7, 1.63036309032582E-7,
    4.66660478371486E-7, 5.30022872802066E-7, 7.60215659946777E-8,
    -1.47269545149564E-6, -8.33566426511535E-8, 1.56358044459513E-7,
    -5.26610554787197E-7, -1.92022721853301E-7, -1.22011971076587E-7,
    -1.67828110052491E-6, -1.62710438101492E-6, -5.04419169569953E-8,
    -2.07528324575027E-8, -3.49214044332391E-7, 2.25444301860948E-7,
    -6.79128966136803E-7, -1.5330669918234E-7, 5.39832430412364E-7,
    -5.02070797583782E-8, 5.52770190018364E-7, -8.30690355762635E-7,
    -6.06732162592585E-7, -7.8159987499796E-7, -2.02347855616198E-6,
    -7.35531232902179E-7, -7.38920647532639E-7, 2.52150781649078E-7,
    -5.96349039412702E-7, -8.52772414666493E-7, 4.30688963661426E-7,
    -7.77293164771209E-7, -5.95366146770899E-7, -2.95675971532358E-8,
    -7.25459963886975E-7, -3.62290731858061E-6, 1.10825468808058E-8,
    6.27937425798512E-7, 7.81456582333952E-8, 2.48629473434464E-7,
    -2.87139290496756E-7, -3.00635375346683E-7, 6.42853505007129E-7,
    2.03936695259321E-7, 5.55055998062734E-7, 5.3038712839182E-7,
    8.37271026615758E-7, -1.42626101128362E-6, 9.49073908081478E-6,
    9.99287601121475E-6, -2.42599463988583E-7, 3.88655478778442E-7,
    -1.99973991278545E-7, -3.14261994493223E-7, 6.8855170014609E-7,
    -8.64997684745066E-8, -4.00324004156285E-7, -1.36436300183419E-7,
    5.97236536456704E-7, 1.47268982687286E-7, 5.49061475597967E-7,
    -1.40123425429514E-6, -2.59426867500574E-6, -4.81851662883646E-7,
    -6.2265987819738E-7, 1.63357385628188E-7, -3.92122358128571E-7,
    -3.54498478687107E-7, 4.88536327900885E-7, -4.80999113413271E-7,
    6.34182760194271E-7, -1.08219057563835E-7, -9.04243436562913E-7,
    1.27109716933012E-7, -2.84233129561759E-6, -7.35806246545524E-6,
    -9.51031549342034E-7, 1.19096657449647E-6, 2.37434436961483E-7,
    8.36124205691084E-7, 2.25657871556309E-7, -7.06554930190212E-7,
    -2.13218594605133E-7, -2.96736076245562E-7, 6.3716380410588E-7,
    3.69261456134313E-8, 8.72155460456674E-7, -1.27773997204808E-6,
    -4.47078131497551E-6, 1.01024018499125E-5, -2.28593903688432E-7,
    -6.23820111542759E-7, -6.12741618341166E-7, 3.72847793814934E-8,
    -1.15294484424379E-7, -5.78913008926242E-7, -4.12498222379079E-7,
    3.07015540831661E-7, -4.40728308784009E-7, -1.14232511878322E-6,
    4.60530814845035E-6, 8.97164857924289E-6, -4.93354113585733E-7,
    2.50886711428265E-6, -2.94221361541117E-6, -2.40493496791571E-7,
    -3.30126421704585E-7, -6.08588453553501E-7, 1.09886566286425E-7,
    -2.24647796131903E-7, -5.67321470205334E-7, 6.28589609139987E-7,
    -2.15970276588242E-7, 3.22418871273614E-7, 4.53400415088877E-8,
    -4.66520855440701E-6, 1.3200361671223E-6, -8.76328454828469E-7,
    9.15715981669718E-7, 1.21750940654623E-6, 1.55989256094216E-7,
    -3.69298118463465E-7, 5.54690666393976E-7, 5.6923248621919E-7,
    9.95217925389516E-7, -1.84329831961806E-6, 1.00943954199519E-5,
    -9.39045824358774E-7, -5.86725060306271E-7, -3.67641719352781E-7,
    -9.81699240127375E-8, -1.12120146676738E-6, -3.20125484365916E-8,
    2.59059776881293E-7, -5.42140032733613E-8, 8.05040263285108E-7,
    2.72091406982215E-7, -1.68441095241457E-6, 1.01681302530904E-5,
    -2.26449239145627E-7, -5.42229338340088E-7, 6.4259157595407E-7,
    -3.06934823707665E-7, -3.43839460612378E-7, 9.60390168238064E-8,
    5.03886264871699E-7, 7.99061414532927E-7, -7.252563294148E-9,
    -1.06020442457302E-6, 1.38397659419379E-6, -1.50692893480563E-6,
    -4.1975334692135E-6, 1.45411268669989E-7, 1.0279616971295E-7,
    8.41106030623035E-8, -1.48253733913851E-6, -1.28966038645984E-7,
    -1.40093091836292E-7, 2.2242706359191E-7, 4.06231215452068E-7,
    -1.67427902132321E-7, -1.33074983004701E-7, -1.47756899865572E-7,
    -4.43205217891177E-6, 1.44572378591852E-6, -1.78054176942962E-7,
    -7.14792657724712E-7, 2.15724584415322E-8, 3.8854117622093E-7,
    1.54098297607272E-7, -2.6329438393943E-7, 1.23692583193082E-6,
    -1.20255244299074E-6, -1.87774860418208E-6, 6.07300502525351E-7,
    1.15910578047933E-8, 3.75125607681113E-7, 5.96680473691959E-8,
    1.23576672480355E-7, -8.56855645312959E-8, 2.73883429072323E-7,
    7.18462635259117E-7, -1.46844738137315E-6, 1.59879037900761E-6,
    4.53620489813996E-7, -6.95330026262405E-6, -3.66284757480322E-7,
    -5.158796852978E-7, 6.92201865587204E-8, 8.92931277103046E-8,
    -1.63814820325163E-7, 1.9011731682313E-7, -4.44790262585369E-7,
    2.02204686926017E-7, 2.92953621969543E-7, -1.79202541433277E-6,
    6.71839301457512E-6, -2.72268721374305E-7, -4.81025354567133E-6,
    -4.25540891329078E-6, -2.42396137578068E-7, -4.29148921835649E-7,
    7.61497363953489E-7, 1.90571907264289E-7, 1.03972096359817E-7,
    6.24784404388995E-8, -2.61672604396101E-7, -3.66977301862851E-7,
    1.11385051753039E-7, -1.43764290126944E-6, -8.10842376586061E-6,
    2.20541089391727E-6, 3.80901303425028E-7, -5.22099918017611E-7,
    -3.33392573496095E-7, -5.35320150237223E-7, 1.51880998226525E-6,
    -3.50790413083779E-7, 4.53659867433827E-7, -8.69746132780962E-7,
    -1.16518685075289E-6, 8.0161102846636E-6, -1.67920495118585E-6,
    -5.50347934223449E-8, -2.7279557594899E-7, -6.7903623327663E-7,
    5.17512225765096E-7, 2.96992684558808E-7, 1.79845721469417E-7,
    5.16590013576955E-8, -1.2273015254292E-7, -6.55910306551882E-7,
    -2.54464884741151E-7, -1.9595833237436E-6, -5.18965591552023E-6,
    1.97870360889406E-7, 3.65023914172544E-8, -2.74262301520385E-7,
    1.28053625954437E-7, -4.81200425800158E-7, -4.95245941191749E-7,
    3.35469203703913E-7, -1.53728491450804E-6, -2.87448532498829E-6,
    5.04709508048634E-7, 9.48427502339515E-7, -2.41875428930724E-7,
    2.58962641477257E-8, 3.22100924459214E-7, -2.30542609036753E-8,
    5.00705777611197E-7, 5.74242224874926E-7, -1.25294883417973E-6,
    4.38140067918634E-6, -5.24683781666892E-7, -8.83201716489544E-7,
    -9.66576975627607E-8, 3.11453955023109E-7, 2.84800685404223E-7,
    7.76781926902685E-7, 5.18999574819696E-7, -3.07671061862331E-7,
    4.15783864142925E-6, 3.89887815218788E-6, -2.54909413057525E-6,
    -3.13351361308993E-7, -2.20622947985284E-7, 3.10616246371789E-7,
    3.00431749903832E-7, 5.68983541800155E-7, -7.82768555001951E-8,
    -1.16169620998889E-6, 1.39986112937007E-6, -2.64443333000987E-6,
    -4.81757588588731E-6, -8.69414568905576E-7, -3.70270112025309E-7,
    2.58099412147547E-7, -2.99571029280572E-7, 3.67290370436976E-7,
    4.30370945522282E-7, 3.50732044904965E-7, -4.17011706053686E-7,
    -3.92686848370678E-7, 4.31089327222902E-6, -2.26529367403005E-6,
    -5.42646227990925E-6, 1.34844822874288E-7, -2.53363306332696E-7,
    -3.05001105997824E-7, 7.16131385322909E-7, -4.95000054862816E-7,
    -2.43069477486665E-7, -1.26331937680738E-7, 8.06935269364336E-7,
    -1.38221893885227E-6, -4.89612733256722E-6, 7.20925779858538E-6,
    -5.70316163664181E-7, 5.34038741807381E-7, -1.59512293056528E-7,
    -4.97796713333462E-7, -8.01353583414593E-7, 6.54304431145212E-7,
    -3.83926706546174E-7, -5.86904928276593E-7, -1.08492922906668E-6,
    -5.29245403405115E-6, 5.93179785525979E-6, -9.25709170302369E-7,
    -5.7020192602449E-7, 4.95720582327738E-7, 7.50611931382776E-7,
    -4.9515915647208E-8, -4.6650866640745E-7, 1.38201086307861E-7,
    9.30600656055539E-7, -1.62991949299549E-6, 7.37700415852725E-6,
    4.16026154151144E-7, -5.85288442673585E-7, -1.02143958925685E-7,
    -2.37895417789923E-7, 3.01056136465939E-7, 3.44296495002264E-7,
    3.9162261191937E-7, 5.99220570595601E-7, 1.6464690717559E-7,
    1.62253775316959E-7, -1.52589805482389E-6, -4.72318268769218E-6,
    8.24289337338292E-6, -2.68985897083455E-6, -7.05034189966573E-7,
    8.24457665984796E-7, 5.50804486372279E-7, 3.53959636592812E-8,
    -2.94861445294333E-7, -3.34748742613728E-7, 1.89317890903065E-7,
    4.21175728220871E-7, -8.2560930473296E-7, 6.71051989891374E-6,
    5.78118004622968E-6, -3.64280986082077E-6, 2.50182427050697E-7,
    -2.71421282294039E-8, 8.77724567356324E-8, 3.7403567281373E-7,
    -8.77390804736175E-8, 5.53285329631232E-7, -5.41275930485371E-7,
    1.5186271608337E-6, -1.32493719745453E-6, 8.53444883185383E-7,
    -2.80632034692174E-6, -1.17705969599355E-6, -6.42880690216192E-8,
    7.02126944065689E-7, 5.00088175040451E-7, -1.26731263115082E-7,
    -1.0752031956964E-9, -5.44093127982374E-8, -1.32151082850169E-7,
    7.32694933707744E-8, -1.09980519393034E-6, 8.54058164659511E-6,
    2.10026252777839E-6, -1.80938093961433E-6, 1.5519244896472E-7,
    -7.84405763977925E-7, 2.95752518406434E-7, 1.4739978841828E-7,
    -4.47268406648451E-7, -1.67239078982978E-7, 2.13729604019403E-7,
    -6.49877688404553E-8, -7.11884773062845E-7, 5.6796343166749E-7,
    -4.7089970112147E-7, -1.39499741574066E-7, 2.69895763681944E-7,
    -6.46359232659546E-7, -3.65218942242194E-7, -4.79694727059533E-8,
    -5.54901964992889E-7, -3.30172052161744E-7, 5.78779978202528E-8,
    1.02403129850991E-7, 3.56679500076638E-7, 1.08600303493301E-7,
    -9.09334712663715E-8, -8.05264583316424E-8, 3.04412690213134E-8,
    -2.43410933919341E-7, -4.41921928193421E-7, -1.08809300709344E-7,
    3.68923189988064E-7, 3.19708794606374E-7, 2.32740445246296E-7,
    4.34578244716437E-7, 5.69933867575287E-7, -3.59416489075378E-7,
    3.59192435326826E-7, -4.59942073991057E-7, -2.88398947841912E-7,
    -1.94938103980061E-6, 9.01274903472196E-6, 1.21897313329079E-6,
    -2.66167521677804E-6, -9.32368488084013E-7, -4.02119316099051E-7,
    -2.78829599192693E-7, -1.21714342499607E-6, -3.13379284427001E-7,
    1.00260754861957E-7, 1.6173347175456E-7, 4.65170527834064E-7,
    4.58261788607301E-8, 4.35095554206286E-7, 1.1039359287303E-6,
    6.61659387313187E-8, 5.77016358423767E-8, -8.45375053178543E-7,
    -4.76327184034328E-7, 2.43550637090977E-7, -9.08188304791491E-7,
    7.39885830743949E-6, -9.74360312287936E-7, -6.86262841673852E-9,
    -9.88297729387503E-8, 4.0759981461498E-7, 7.91608551664421E-8,
    -5.39990493380848E-7, -4.27014189960787E-7, 3.07025148349034E-7,
    -8.29638242091119E-7, 6.857236875758E-8, 1.43140132275767E-6,
    -2.12071457847958E-6, 2.08298279645839E-7, 4.38434885624327E-7,
    3.24005084766601E-7, 1.93438190091785E-7, -6.48143733944862E-7,
    6.6442371448708E-8, 7.36600678210769E-7, 2.01551180877062E-7,
    -1.86260322299232E-7, -3.99137598952127E-7, 1.47277357590159E-6,
    2.22697142465229E-7, -3.74567002988658E-8, 7.88005738806102E-7,
    5.04030828997628E-7, -7.26061159588667E-8, 1.46316393817601E-7,
    -4.0465281651915E-7, -7.34714284016414E-8, 3.67806904762475E-7,
    8.21626498853118E-7, -5.23883936509977E-7, 7.41292918363709E-6,
    8.91593662706984E-6, 2.99295812417775E-7, 1.94918213499824E-7,
    3.72302853166659E-7, 7.26545523996677E-7, 6.6350776421884E-7,
    -3.53100894520831E-7, -6.2990284824114E-7, -1.63213525052159E-8,
    2.48955862246157E-7, -9.68044008742802E-7, 8.52835201826781E-6,
    5.31310511406796E-6, -3.08057487653019E-7, 5.76567840641953E-7,
    -6.97422564412399E-7, -2.25910167452437E-7, -4.94519785558511E-7,
    2.93681960303791E-7, -3.56335724120823E-7, -3.25536690908933E-7,
    -2.77270327180361E-7, -8.19834356187605E-7, -6.98577354078522E-6,
    -4.57701579582458E-6, -7.04472393845303E-7, 2.46830806708165E-7,
    8.56868858377918E-8, -3.04548131064788E-7, -1.31259804771641E-7,
    9.62759790001995E-7, 6.15263161609532E-8, 1.63347527092684E-6,
    -5.46494029229836E-7, -1.12285810373572E-6, 7.32259413422856E-7,
    6.21414339599224E-7, 7.50631399103618E-7, 5.14984773260652E-7,
    5.63230952842554E-8, -2.1575807630182E-7, -2.10556808225799E-7,
    1.38662156261106E-7, -1.14600610565074E-7, -8.76709832653743E-7,
    -2.27182390743367E-7, -3.63172692397033E-7, -4.60352024061005E-7,
    3.28767030692365E-7, -1.56458591549597E-6, 5.88495837461689E-6,
    -1.27371172079559E-7, 6.9255432552134E-7, -1.66661631572967E-7,
    -5.61722128308469E-7, 3.1995768040921E-7, -3.32267119836717E-7,
    5.38570798177314E-6, -3.80915162477671E-7, 4.17167002502058E-7,
    5.70840804331297E-7, -3.18434251516756E-7, 6.51470431800968E-7,
    1.06679283263009E-7, -7.77170211209132E-7, 5.43433309161039E-7,
    -4.62181403709485E-7, -4.41902490959503E-6, 8.41480998414467E-6,
    -5.80258597270991E-7, -6.55980042055472E-8, -8.0261069584768E-7,
    3.13806029950866E-7, -6.53171627039138E-7, -2.32446895634561E-7,
    8.78903046890248E-7, -7.10710076080408E-7, -2.22829521446025E-6,
    1.11773366091782E-5, 3.94868923943946E-7, 1.04092507043933E-6,
    3.85641247491641E-9, -6.09020617823204E-7, 1.32584520014147E-7,
    4.37901156701094E-7, 1.22020665098517E-7, 5.4088834098976E-7,
    -1.29087964257764E-6, 5.98330860176741E-6, 1.5541153553927E-6,
    -4.13570663619231E-7, 7.39706763034032E-8, 8.80571540080679E-8,
    4.81922443648394E-7, -4.35085211547331E-7, 1.05445761617979E-6,
    7.71689850085654E-7, -3.77681226026062E-7, 5.08901476582619E-6,
    -1.0092336624876E-6, -3.78554191024591E-6, 3.02670859269351E-7,
    1.32146974583893E-7, -7.01726565706499E-7, 1.30839528937705E-6,
    1.04597580199473E-6, -3.29600451312187E-7, 1.52654027005908E-6,
    2.52791386355691E-6, -3.83370591419447E-6, -2.78994590968961E-7,
    -6.71877403285712E-7, -1.88942226807853E-7, -2.35232666763798E-7,
    -1.48431952509483E-9, -1.14684689944402E-6, 1.08834420009321E-6,
    -3.65914816542362E-6, -3.35588623342772E-6, 2.62886391458325E-7,
    -1.24086435955484E-6, 3.03769738550538E-7, 3.13687603061245E-7,
    -9.51941095648412E-8, 6.1786747431574E-7, -4.5779888129067E-7,
    -6.42726743398469E-7, -3.75603323726508E-7, 9.58690844655312E-8,
    -9.91276189504886E-7, -8.90043782372324E-7, 7.36840042888446E-7,
    -1.32779014532245E-6, -9.39778084724582E-7, -5.30561104664465E-7,
    -6.11076852050883E-7, -3.82350097515175E-8, -9.79184827624799E-7,
    -4.1597250324E-6, -6.08408706728088E-7, 3.86094455782139E-7,
    2.26553950277178E-7, 2.68958686106123E-7, -6.59641042656877E-7,
    4.3461033154215E-7, -6.45663129904225E-8, -1.02303801699555E-7,
    3.36038227735756E-8, -7.61907328105703E-7, 7.51176193433239E-7,
    6.41204979693598E-7, -7.84777599729511E-7, -2.89082383526522E-7,
    -6.83283616918552E-7, -4.85448549320076E-7, 1.57983872391968E-7,
    -3.96181813320725E-8, 5.4265228791329E-8, 4.41652887882065E-7,
    -1.92872286238472E-7, 3.80138249442212E-8, -4.27484283924544E-8,
    -7.34661328861748E-7, 2.84856481710732E-7, 4.6542294502749E-7,
    -9.21003861248246E-7, -2.60328992698652E-7, -2.61728718665754E-6,
    -7.26304212029859E-7, 5.28457658520818E-8, -2.34582955034869E-7,
    3.42709302073351E-8, -1.19803472875379E-6, -7.09565313146418E-6,
    -7.54006019323857E-6, -1.5209562498293E-7, -4.97408074462914E-8,
    1.07079078145286E-7, -7.35502622739587E-7, -4.35550486102784E-7,
    2.89537414855328E-7, -3.60544129822769E-9, 3.59812703807974E-7,
    -3.07737444243747E-7, -3.1442201785191E-7, 4.97361116404047E-7,
    -1.6593816175637E-6, -5.22589296128378E-6, -3.33528882178279E-6,
    -6.1100574433531E-7, -1.01323953111579E-6, 2.25364635629454E-7,
    3.00205168998514E-7, -4.87329354099516E-7, 3.65882416416079E-7,
    1.42225347024376E-7, -3.51015192795452E-7, -3.63704296096916E-7,
    -6.23606995774812E-7, -4.98653170676253E-7, -9.96767921292078E-7,
    -6.42539166884778E-6, 4.34532644300995E-6, -4.91416970253217E-7,
    3.46891921867206E-7, -1.8224364046532E-7, -4.25080975275969E-7,
    2.94578699673017E-8, -4.60493741175789E-7, -7.05114374483299E-7,
    -3.51968747890157E-7, 9.51053466595897E-8, 2.31289869484194E-7,
    -8.69731324462925E-7, 8.4034213771355E-7, -1.87291030025115E-6,
    4.58854334871803E-8, -3.42371203907856E-7, -6.92211462110594E-7,
    2.62697541689396E-7, -5.22227666335688E-7, -3.82433944994586E-7,
    -2.50738579379026E-8, -1.2141148158592E-6, -5.6221093060761E-6,
    -5.69534277277005E-7, 6.33208515146002E-7, -3.31206973942456E-7,
    2.83924980141854E-7, 2.99613661166535E-6, 2.72512114841905E-6,
    3.34925011149379E-6, 3.5704714915038E-6, 3.57228598531562E-6,
    3.57746411468121E-6, 3.40143051612268E-6, 2.64823500162601E-6,
    2.56700369426273E-6, 2.58815103929614E-6, 2.40655880268736E-6,
    2.57247725895056E-6, 2.79448851549388E-6, 3.04518801013605E-6,
    3.24799496046383E-6, 3.89613993280173E-6, 2.92485106798821E-6,
    3.51268834578631E-6, 2.85718439941956E-6, 3.64070893537906E-6,
    7.00329340206195E-6, 8.25957392751978E-6, 8.76060509191027E-6,
    8.90623675963764E-6, 6.2957163005847E-6, 2.96519261686137E-6,
    3.55135872258598E-6, 1.71743574759111E-5, 5.06621875971417E-5,
    5.24006602551024E-5, 5.21867239080609E-5, 2.20045789596414E-5,
    3.51481300631157E-6, 3.36442240320249E-6, 3.20682142165378E-6,
    3.48928254341875E-6, 4.49831598213105E-6, 2.84729439533592E-6,
    3.55723243852392E-6, 2.46952705260973E-6, 3.8488766786047E-5,
    5.41496173165739E-5, 5.56357875516694E-5, 5.50750147838114E-5,
    1.07950882522955E-5, 5.25676274932707E-6, 5.56729015767863E-6,
    5.43327479184165E-6, 5.94301695352909E-6, 5.80266664181815E-6,
    5.04632895178671E-6, 4.17315880171973E-5, 5.43653549318046E-5,
    5.55414441386979E-5, 5.50762275353044E-5, 1.51967543035494E-5,
    3.85366348652627E-6, 4.27465661023675E-6, 3.86878731237079E-6,
    4.579580146173E-6, 4.59065715728161E-6, 3.94436933823454E-6,
    4.50237745878029E-6, 3.93155660164645E-6, 7.15147774285243E-6,
    5.07372339001201E-5, 8.03562756832869E-5, 8.06074793559013E-5,
    8.14239225433925E-5, 5.9634384850526E-5, 6.29334062221197E-6,
    5.38848970008698E-6, 3.40749602188955E-6, 3.04199391979043E-6,
    4.68023765652132E-6, 3.03276065200105E-6, 4.05777065308215E-6,
    3.83348165300744E-6, 1.91803683518778E-5, 6.22445391903588E-5,
    7.67906448023826E-5, 7.71593910671325E-5, 6.60996967133665E-5,
    2.14965778809776E-5, 4.02836295503206E-6, 3.07757656064632E-6,
    3.95307359608428E-6, 2.95238601132187E-6, 3.86298796689172E-6,
    4.98068266375018E-6, 4.33783719084648E-5, 7.557505847518E-5,
    9.09450673872869E-5, 9.21309581293712E-5, 7.08893780173147E-5,
    2.89610663161122E-5, 6.22650343507993E-6, 7.35719716840623E-6,
    3.96878254336469E-6, 4.54298553668143E-6, 3.12202979883363E-6,
    4.90406639195788E-6, 3.83796285175518E-6, 3.69078634042373E-6,
    1.66646059230832E-5, 4.58589501479257E-5, 7.95081714947771E-5,
    8.45553171374169E-5, 7.00591615800341E-5, 6.13652352136007E-5,
    3.81595355860208E-6, 3.2978534477026E-6, 3.57827939057729E-6,
    3.12288444798097E-6, 3.18538696420482E-6, 4.19886617779842E-6,
    1.67640803662011E-5, 5.08297684487257E-5, 7.07951005957984E-5,
    7.1077282584265E-5, 7.04576033103496E-5, 6.70592158318959E-5,
    4.88746432097901E-5, 4.44928869597392E-5, 1.15406194441524E-5,
    3.3075214685241E-6, 3.5296058676165E-6, 4.52294072102914E-6,
    3.79355027306584E-6, 3.64987828872993E-6, 4.65781505278911E-6,
    4.56214041619982E-5, 6.28090610069223E-5, 6.48045010408393E-5,
    6.25020034476303E-5, 1.78291908716102E-5, 7.04340109221355E-6,
    5.54165762368774E-6, 6.89508769181105E-6, 6.26657618376098E-6,
    4.62473508884876E-6, 4.25196329151035E-6, 1.45252771259011E-5,
    5.49182740765128E-5, 6.18685355142691E-5, 6.20993217678111E-5,
    5.41161905138408E-5, 5.17697391245454E-6, 4.05832080933587E-6,
    4.21407886467268E-6, 4.68347521533661E-6, 6.65038257077073E-6,
    3.77202319862239E-6, 3.81355236998549E-6, 4.89999575045033E-6,
    5.33526757117311E-5, 5.71507585887401E-5, 5.73950433376445E-5,
    5.63069246394626E-5, 5.98493611091363E-6, 4.7998563422057E-6,
    4.6396049673082E-6, 3.83944980477667E-6, 3.85014779410577E-6,
    3.52954414233701E-6, 4.22319797681436E-6, 4.36144462990196E-5,
    7.35462709179535E-5, 8.30752620339018E-5, 8.4208134543368E-5,
    6.0794833089958E-5, 1.80171260827097E-5, 6.01344004352862E-6,
    5.19288681307519E-6, 4.33146393723511E-6, 8.22391192299635E-6,
    4.02489831458962E-6, 3.6831360531888E-6, 3.31967810128712E-6,
    3.6008374870457E-5, 5.19156950247573E-5, 7.5124081692792E-5,
    7.74912442588214E-5, 5.78533999236029E-5, 2.68373204592357E-5,
    4.36813007104099E-6, 4.85933479717689E-6, 3.89735350172782E-6,
    4.04727131749125E-6, 1.64206441778875E-5, 4.76591414960361E-5,
    7.36672002198084E-5, 7.98001582830388E-5, 6.18152660617796E-5,
    5.54309578784142E-5, 5.78979427988296E-6, 3.92736485725817E-6,
    3.83062530195886E-6, 3.31751912112027E-6, 1.88263015869346E-5,
    5.23827344029105E-5, 7.32369336096163E-5, 8.8694379546014E-5,
    7.29948361388686E-5, 5.93499101295961E-5, 2.57506963519426E-5,
    2.37382360774849E-6, 3.9725056106799E-6, 2.80746663929212E-6,
    3.10755663511536E-6, 3.09417906296206E-6, 2.39591086698252E-6,
    1.51909727644306E-5, 4.54556620152453E-5, 4.9207809598499E-5,
    6.26655287318988E-5, 7.86646767920064E-5, 7.6872426946558E-5,
    7.99127029039421E-5, 2.78512266694424E-5, 3.61526592805267E-6,
    3.76871988493329E-6, 4.96181746529575E-6, 4.03516579788589E-6,
    3.52938016602973E-6, 2.84288304087068E-6, 2.98572652163986E-5,
    5.90200768448167E-5, 7.99953008056045E-5, 8.17391123659828E-5,
    6.92518920251532E-5, 2.69467036894147E-5, 5.15999347740376E-6,
    4.04661597748433E-6, 4.12387705914259E-6, 4.06386934425898E-6,
    7.02833078242966E-6, 3.72560109342825E-6, 4.87955634258403E-5,
    5.84682854266633E-5, 5.8840358896956E-5, 5.25872326082846E-5,
    7.8164644216153E-6, 3.20235390022909E-6, 3.30225593716532E-6,
    3.22038851507482E-6, 2.82346349259091E-6, 2.84285200695726E-6,
    3.27341074942223E-6, 4.33986515236066E-6, 2.5341086260604E-5,
    5.1634921891107E-5, 5.33810468754525E-5, 5.34076578019261E-5,
    2.06713480112751E-5, 3.39285730327574E-6, 2.67774316417198E-6,
    3.33055523965079E-6, 3.74945867477131E-6, 4.48478995935622E-5,
    5.15483152709466E-5, 5.31043465271284E-5, 5.27061836620365E-5,
    1.26517616732112E-5, 4.83549220918311E-6, 5.34520038538825E-6,
    3.58589298462521E-6, 3.16761831584766E-6, 3.57247689953592E-6,
    1.35335704994892E-5, 4.52240520397644E-5, 4.73367545247638E-5,
    4.72468121351614E-5, 3.26968015270013E-5, 3.39501138520501E-6,
    5.52479506670903E-6, 4.04915853382241E-6, 4.19516833170521E-6,
    4.21541037618314E-6, 3.70478002750911E-5, 5.58958616900624E-5,
    6.70505481148631E-5, 6.91720807307066E-5, 4.98151574761449E-5,
    1.06026151097828E-5, 4.71591191385882E-6, 4.29338144417878E-6,
    3.61454343610556E-6, 3.92259616064565E-6, 3.46489521413602E-5,
    4.41711466127222E-5, 6.95613172277547E-5, 7.00894358274408E-5,
    5.44527600577126E-5, 2.02139349584067E-5, 5.40609880257563E-6,
    4.25091228324292E-6, 4.6753946937027E-6, 4.46458482481067E-6,
    4.92415192663346E-6, 4.01770189804234E-6, 3.97792160489681E-5,
    6.28480218285486E-5, 7.22208322260172E-5, 7.30782748676854E-5,
    5.53354133510519E-5, 2.0334521810369E-5, 4.04451166632331E-6,
    4.41667643872101E-6, 4.98526592503216E-6, 5.828373706353E-6,
    5.56400164818217E-6, 4.05477105219206E-5, 4.80077770562574E-5,
    7.76251160426226E-5, 8.28409633977413E-5, 6.59210582923723E-5,
    4.37236748286054E-5, 5.77648233240165E-6, 3.92458898103038E-6,
    2.71442347768323E-6, 3.48916416407267E-6, 4.14781169531991E-6,
    3.43795598809134E-5, 5.90941181031342E-5, 7.5468090042502E-5,
    7.98629744284133E-5, 6.19151133113117E-5, 4.14904221773706E-5,
    6.88285281492575E-6, 5.01912225956404E-6, 3.07957438103992E-6,
    4.34959496476632E-6, 1.53497154902289E-6, 1.49403712547175E-5,
    5.45305435196273E-5, 7.76485462200722E-5, 8.98191400118416E-5,
    8.59523610288434E-5, 5.86746705425902E-5, 1.0764548522833E-5,
    3.35015458242979E-6, 2.48655100136955E-6, 2.76426469084559E-6,
    4.24953937866972E-6, 4.21461931861959E-6, 1.74057803941963E-5,
    4.67142801633418E-5, 6.9290242904035E-5, 7.47845217012391E-5,
    6.12337941034843E-5, 4.67149925518762E-5, 1.25994496345924E-5,
    4.86900490169015E-6, 5.44518711955964E-6, 4.90281334717307E-6,
    4.27535888715256E-6, 3.57054348504071E-6, 7.31006221364636E-6,
    5.10649691880101E-5, 7.09043914078212E-5, 7.22619214009202E-5,
    7.68847941524371E-5, 4.10528288493962E-5, 1.66013558283645E-5,
    4.77469640936989E-6, 3.61093249727001E-6, 4.24630800541541E-6,
    5.27946606653651E-6, 4.05745563599282E-6, 1.80790910465968E-5,
    5.1374682323938E-5, 6.75473956945901E-5, 7.26044207993991E-5,
    6.2674984747668E-5, 5.56351535904691E-5, 8.98871646134995E-6,
    4.78515749343296E-6, 3.93206648362298E-6, 3.00058750338106E-6,
    2.46886051162581E-6, 2.07974906026565E-6, 5.64284915511134E-6,
    4.30111774784606E-5, 6.86918613319276E-5, 7.95332280465413E-5,
    8.49424239542522E-5, 5.32610466664646E-5, 8.7940450131724E-6,
    3.34918460919935E-6, 3.8416152289891E-6, 2.20132377032709E-6,
    2.29087289074503E-6, 2.85626323311929E-6, 6.44756728199207E-6,
    5.87734413476308E-5, 3.99256368709442E-6, 4.01078537640644E-6,
    4.13549064795811E-6, 3.16211010195541E-6, 3.91681530820563E-6,
    3.6206250644165E-6, 4.40989946660087E-6, 3.13697693253009E-6,
    2.81982270160358E-6, 3.69924153956583E-6, 3.62070619739507E-6,
    3.68688194666972E-6, 4.15777895445071E-6, 4.63326851113088E-6,
    4.2440229288666E-6, 4.12609394567279E-6, 3.47736436587233E-6,
    3.23525695610736E-6, 3.23559684302811E-6, 3.19132711657608E-6,
    3.24800347639462E-6, 3.38831270887235E-6, 3.23909852613738E-6,
    8.99110623619116E-6, 7.66392132381244E-6, 8.55462742188169E-6,
    8.12454069673811E-6, 1.12163347541792E-5, 5.10031851241997E-5,
    5.58705046197747E-5, 5.67178202092221E-5, 5.05248624649765E-5,
    6.01293079731308E-6, 1.07313987954831E-5, 4.5619410931017E-6,
    3.278352729884E-6, 2.44218972306449E-6, 5.2102406850991E-6,
    2.86467765254225E-6, 3.15078455304892E-6, 4.48837316237124E-6,
    4.69159141819496E-6, 4.22387455485434E-6, 3.84994728259045E-6,
    6.37886650964164E-6, 2.84966101949718E-6, 3.06611569763497E-6,
    1.49404752239308E-5, 4.99980595184818E-5, 5.26682975538561E-5,
    5.27232341599314E-5, 4.76504990520461E-5, 5.76731544412367E-6,
    3.99353119120367E-6, 3.47531054501939E-6, 5.18270013822645E-6,
    3.48960582556012E-6, 4.37366443057932E-6, 3.58030525712915E-6,
    4.41310750464293E-5, 5.46108512549037E-5, 5.49754900143981E-5,
    5.35143724260435E-5, 9.21269819209655E-6, 3.82931681537725E-6,
    3.67220902125653E-6, 3.42093820565461E-6, 3.41414557106369E-6,
    3.7445606127682E-6, 4.50886289526162E-6, 5.45577892125924E-6,
    3.11327375069928E-5, 3.71039803674289E-5, 3.70332782602805E-5,
    3.54424072560757E-5, 5.23756150339856E-6, 3.43622817538838E-6,
    4.76064348572383E-6, 4.43547879114636E-6, 3.60406856212476E-6,
    4.02497781666058E-6, 3.64425801682281E-6, 4.61260807928957E-6,
    5.40743357548039E-5, 8.30818865632846E-5, 8.89440072619853E-5,
    9.11196078998744E-5, 6.13624771595352E-5, 4.16182133737307E-6,
    2.81750524782855E-6, 3.51994531339411E-6, 4.40447743512413E-6,
    4.04568995442457E-6, 3.1496559658583E-6, 7.36452826572787E-6,
    4.79252350427228E-5, 7.86728851795971E-5, 8.055120388146E-5,
    8.53317790598875E-5, 5.65689128300917E-5, 3.55013998090827E-6,
    4.11407028712133E-6, 4.52625676361886E-6, 3.84244317882215E-6,
    3.94358745289516E-6, 3.89840667348299E-6, 2.99427950989079E-5,
    4.87866056769304E-5, 6.95785510302311E-5, 6.98816003115208E-5,
    6.04549122398419E-5, 2.64513379989706E-5, 3.97657430528837E-6,
    3.77326051134017E-6, 3.49693269834797E-6, 4.19045259405404E-6,
    4.22756202910173E-5, 7.33714384858474E-5, 8.0786269442203E-5,
    8.08987252225124E-5, 6.17220184570366E-5, 6.65154078083684E-6,
    7.51989728211001E-6, 6.8629327575393E-6, 6.19757385189194E-6,
    4.03551106525793E-6, 3.78683923249074E-5, 4.6781869248819E-5,
    4.7625668262401E-5, 4.7221195988602E-5, 3.6161687311511E-6,
    5.31167536871377E-6, 3.47513580814862E-6, 1.12197725518482E-5,
    5.23574226118336E-5, 5.42458518983452E-5, 5.41985328838756E-5,
    4.56943536554903E-5, 2.76408361975219E-6, 4.73096646036018E-6,
    4.10896983440831E-6, 4.12422250787783E-5, 4.32128418200475E-5,
    4.39736249829466E-5, 4.35062929557783E-5, 5.44057528163483E-6,
    4.51348759855923E-6, 4.20337274124438E-6, 5.74641576459508E-6,
    5.75081036202445E-6, 1.73184500571364E-5, 5.21484029773843E-5,
    8.32834232602847E-5, 8.53740092731975E-5, 7.57784402516597E-5,
    5.15779933072613E-5, 4.24880071664898E-6, 3.82028410063992E-6,
    5.54894958033031E-6, 5.7335106914637E-6, 1.59039285879874E-5,
    4.3786391936912E-5, 7.93846826058298E-5, 8.24568896962937E-5,
    7.45159202412449E-5, 6.32053498177827E-5, 7.66896936506875E-6,
    5.74566996254739E-6, 4.60956953454359E-6, 5.95350103230669E-6,
    4.89466462052587E-6, 5.27240448047774E-6, 3.80838006812154E-5,
    6.58844025810201E-5, 6.88667980583428E-5, 7.02315478096755E-5,
    5.08538475283915E-5, 5.1339897976474E-6, 4.43224492140464E-6,
    5.06421622891044E-6, 5.19710622348903E-6, 5.58987857523346E-6,
    5.74942989034147E-5, 9.08464256170885E-5, 9.71640719524622E-5,
    9.60399515499583E-5, 6.19024307694663E-5, 1.67330389439975E-5,
    7.895039775664E-6, 5.97083273290715E-6, 7.23766707165698E-6,
    5.35547390859609E-5, 8.89823632029374E-5, 9.8436954781646E-5,
    9.80077602472795E-5, 6.50640948835854E-5, 1.46331762566584E-5,
    5.38599950637454E-6, 5.85195775708312E-6, 4.8127778148817E-6,
    4.45445391213327E-5, 5.7590450380748E-5, 8.6607491614969E-5,
    8.84707498672395E-5, 7.07516054081547E-5, 2.70988348617217E-5,
    3.19032465494327E-6, 6.29270230634534E-6, 3.22378056810727E-6,
    3.54630750583374E-6, 3.84035093636784E-6, 4.50425108089968E-6,
    4.45883879898895E-6, 4.14220997125688E-5, 4.8611084645869E-5,
    5.05703981946576E-5, 5.0703951909432E-5, 1.11931365518207E-5,
    9.36534612097564E-6, 1.02950380932755E-5, 5.35512723899651E-6,
    3.9328502771927E-5, 4.94541855530759E-5, 5.02439545890061E-5,
    5.03198627284666E-5, 1.61275100701913E-5, 5.71538445851476E-6,
    4.41108105782289E-6, 3.31320478263367E-6, 3.06630201814251E-6,
    4.78445466470281E-6, 4.44106604169995E-6, 3.85166328050224E-6,
    8.14522117747109E-6, 9.80123170797475E-6, 9.73198078093842E-6,
    9.80161635665862E-6, 7.38735521714992E-6, 7.52458328964049E-6,
    5.07289294292643E-6, 6.08842168119275E-6, 4.40488078799189E-6,
    4.3502381024663E-6, 4.36462198573762E-6, 4.57641517530983E-6,
    4.70370942663583E-6, 6.25309150932274E-6, 1.06320838214236E-5,
    1.30150088594815E-5, 1.24140788704103E-5, 1.58647212710922E-5,
    2.36051060385816E-5, 2.46195454060085E-5, 2.5170536002842E-5,
    1.0697839951512E-5, 5.53679297190484E-6, 3.39064337945428E-5,
    7.33383821282215E-5, 8.49377733836002E-5, 8.61169185661822E-5,
    8.12830128367936E-5, 2.95091884181501E-5, 3.90118185989933E-6,
    4.012000470505E-6, 4.11621646343797E-6, 4.33159292065012E-6,
    3.08219820744002E-6, 2.71445856364749E-6, 3.91415660315981E-6,
    3.32401542931464E-6, 1.94293422056092E-5, 6.68285438746371E-5,
    8.38023961340823E-5, 8.48589013143557E-5, 7.62960435093884E-5,
    2.48729224413987E-5, 4.42123413481215E-6, 5.57269277580341E-6,
    4.08182610741938E-6, 3.87093109344697E-6, 3.56311417010694E-6,
    3.2126090410882E-6, 2.97367418636939E-6, 2.71881627140405E-6,
    1.15550750178655E-5, 4.08143616568466E-5, 6.19004372702329E-5,
    6.30941734274797E-5, 5.77878400099285E-5, 3.14320538409802E-5,
    4.59543412049927E-6, 3.72891144490448E-6, 4.2808838968956E-6,
    3.9203309015539E-6, 4.00409085876673E-6, 4.41921470339603E-6,
    4.8198166480823E-6, 3.42665626428859E-6, 4.14651461226287E-5,
    5.22617930297481E-5, 5.31539791456535E-5, 5.24921496403447E-5,
    8.10330858372169E-6, 3.86483501945309E-6, 3.13460894755591E-6,
    4.03625710763347E-6, 3.8450747426427E-6, 2.58281892681049E-5,
    4.82634254045784E-5, 4.96715451922572E-5, 4.96534001784986E-5,
    2.30489267248033E-5, 4.71138413077017E-6, 4.44840652228877E-6,
    4.10311508581094E-6, 4.78917350937397E-6, 7.71780188389997E-6,
    7.03922968001312E-6, 8.77099318004451E-6, 7.72221599322951E-6,
    8.72240933510269E-6, 9.01216028580313E-6, 7.91705933129176E-6,
    9.71012744418553E-6, 3.63927732931736E-6, 4.85973962303358E-6,
    5.08337649721472E-6, 5.83448538304837E-6, 5.97319326936747E-6,
    6.19575813131141E-6, 6.09054737270892E-6, 9.43751565354675E-6,
    6.16249307906997E-6, 7.62411478417764E-6, 7.60476752921986E-6,
    7.19936903478282E-6, 1.92229843655058E-5, 2.89044596427211E-5,
    3.02411952762425E-5, 3.20524372655847E-5, 6.83340014928872E-6,
    5.90747885365704E-6, 7.93604613933683E-6, 3.62993761798427E-5,
    0.000171214111782448, 0.000176834290255517, 0.000176696016629046,
    5.20488363876914E-5, 9.29513524696439E-6, 8.93394657704911E-6,
    7.34067571367851E-6, 8.45813705247355E-6, 7.7828986594271E-6,
    8.60810641232185E-6, 1.2786837427223E-5, 6.75054707794176E-6,
    0.000126682113117918, 0.000175030828878265, 0.000180558640980506,
    0.000180142700573427, 1.57852810721031E-5, 9.61121532452453E-6,
    8.96612298207062E-6, 8.46584715070011E-6, 9.05888515395467E-6,
    9.99200838014242E-6, 1.09934518511056E-5, 0.000119427949590588,
    0.000184563767581232, 0.000189197242217301, 0.000188552024764323,
    1.76162257037378E-5, 8.73236139943908E-6, 1.06342750728569E-5,
    8.33868257144363E-6, 8.91065768894465E-6, 8.78503157083306E-6,
    8.75864879407302E-6, 9.98812777013474E-6, 7.75360200630577E-6,
    1.2800402068905E-5, 0.000124001527669781, 0.000170005137625753,
    0.000175869172036123, 0.000172568795029189, 0.000225356977840176,
    1.25019848114729E-5, 1.32617775098681E-5, 8.58753816373763E-6,
    6.5174151521272E-6, 1.29634549436694E-5, 7.06254141334501E-6,
    8.33982880350906E-6, 6.76513617004409E-6, 3.93542063232601E-5,
    0.000211708753272182, 0.000190588094976856, 0.000189977956637937,
    0.000187727232241485, 6.5232615428176E-5, 6.96014370213152E-6,
    8.62022150320157E-6, 9.3346409937334E-6, 6.18468936225166E-6,
    8.01927331966168E-6, 1.0637343597694E-5, 0.000142465982654588,
    0.000219103319702104, 0.000207023406123406, 0.000204699504413481,
    0.00019337699840022, 4.71815343547189E-5, 1.11184059380836E-5,
    2.93387134388429E-5, 9.95028885387643E-6, 1.71388640337323E-5,
    7.3217322843828E-6, 6.80256468280376E-6, 5.92447203953564E-6,
    6.67336345921464E-6, 4.43052953723228E-5, 0.000163171433592816,
    0.000189222045592125, 0.000197175526343604, 0.000204619190051856,
    0.000205565615210374, 6.07573818688862E-6, 5.83846859082489E-6,
    6.29135289625969E-6, 7.88156827580864E-6, 7.23285164766811E-6,
    9.28412507457357E-6, 3.25551236636028E-5, 0.000180738755266141,
    0.000184419916593445, 0.000184697782607153, 0.000169689834766563,
    0.000195635833315542, 0.000168128869053387, 0.000164347522716442,
    1.94487193505547E-5, 5.90974382734468E-6, 8.65821320797801E-6,
    8.48624216565845E-6, 8.61646136158153E-6, 7.42478037187819E-6,
    8.68790212425792E-6, 0.00011508679441164, 0.000170982006747438,
    0.000175686437904501, 0.000175781437545955, 2.05276987042079E-5,
    2.0091162255163E-5, 1.34192604584476E-5, 2.37196540762592E-5,
    2.17203904051435E-5, 1.21883754524987E-5, 1.23701011149091E-5,
    3.41339966918524E-5, 0.000120406418171415, 0.000140291522140902,
    0.000141890285593998, 0.000146366619876853, 8.2246089001806E-6,
    6.91394515030639E-6, 6.38196955921151E-6, 7.55619567200587E-6,
    8.08004745642704E-6, 7.13832934987173E-6, 6.93180112048461E-6,
    1.02552252301126E-5, 0.000173928480985996, 0.00018467360517765,
    0.000185448651134412, 0.000205319653365707, 1.27501270760362E-5,
    1.13473361150481E-5, 1.53238163343737E-5, 6.46275428577644E-6,
    6.49443937927581E-6, 6.63720111376938E-6, 1.09857464250675E-5,
    0.00015960058188021, 0.000213916001547755, 0.000202217732819723,
    0.000199791636224781, 0.000185786056759335, 2.8250403845985E-5,
    1.45313857670097E-5, 1.27962528422405E-5, 1.00272548368087E-5,
    1.33963606530824E-5, 6.37789950351962E-6, 5.65228816810167E-6,
    5.54049698115387E-6, 0.000118947466122299, 0.000180681956149894,
    0.00018765777767706, 0.000187292728656142, 0.000149506859885304,
    8.43937466506638E-5, 9.42563325412779E-6, 9.81694954797317E-6,
    8.70633922378357E-6, 1.06356989222813E-5, 3.29159529302199E-5,
    0.000167588654884104, 0.00019106337827393, 0.000202988445653005,
    0.00020658702534677, 0.000194167119818647, 1.32110353445923E-5,
    9.1135782996977E-6, 8.25590135222275E-6, 6.04651993921227E-6,
    3.91768738336766E-5, 0.000183712940138297, 0.000212521315350599,
    0.000212242107664415, 0.000131480582584359, 0.00012185716341675,
    3.19013802673913E-5, 6.01710029326363E-6, 7.56531635553744E-6,
    7.46404899742896E-6, 7.74941607403038E-6, 7.70048219319752E-6,
    5.50149030894963E-6, 3.53652155874619E-5, 0.000116558430293055,
    0.000123889834239001, 0.000122240881305039, 0.000186230013966893,
    0.000199720885774587, 0.000202950361756535, 6.81759198574489E-5,
    5.0480674759336E-6, 5.55623087538805E-6, 1.34778373324442E-5,
    7.70078923157955E-6, 6.55988906023309E-6, 6.55024863009883E-6,
    9.10018728711259E-5, 0.00020442648435568, 0.000204863340071687,
    0.000203465007276813, 0.000183932514686886, 9.46916374924833E-5,
    1.06476361480539E-5, 6.13533146408052E-6, 7.46785731214856E-6,
    7.66930123487177E-6, 2.42564130633543E-5, 7.42421604005616E-6,
    0.000135055492763724, 0.000157091919169652, 0.000159823368764802,
    0.000172055190695344, 1.47755429128494E-5, 7.06852172271378E-6,
    6.80682061459082E-6, 4.03421095794842E-6, 8.76552671598768E-6,
    7.65198701493202E-6, 5.7984028412257E-6, 1.02173211356063E-5,
    7.32553110804496E-5, 0.000180386530004898, 0.000186916004017518,
    0.000187242551630581, 3.40075198234931E-5, 1.19276111778562E-5,
    7.41630582281803E-6, 8.08103157946906E-6, 7.79875432954762E-6,
    0.000130879812221589, 0.000182731518997476, 0.000188255842365212,
    0.000188114395036704, 2.2527637230838E-5, 1.00965292469557E-5,
    1.2445067837784E-5, 8.04783213770746E-6, 7.99604800162854E-6,
    7.8325399307789E-6, 2.87894665629589E-5, 0.000160711481984668,
    0.000167537373700869, 0.000167596865683523, 0.000137117311338016,
    4.33972992519375E-6, 1.11499696260935E-5, 8.85640237524758E-6,
    9.68183103551247E-6, 9.62581340438162E-6, 0.000132653143169043,
    0.000150734114648484, 0.000163719703076256, 0.000166598905491913,
    0.000163909669201214, 1.55912767001839E-5, 7.92986146534875E-6,
    1.00780947507522E-5, 8.44370761894035E-6, 1.05357938917927E-5,
    0.000123493265036258, 0.000147432408122694, 0.000189025565194555,
    0.000194526757675381, 0.000190143715160877, 4.1199962285906E-5,
    1.04335845937258E-5, 1.09968310494907E-5, 1.15771317729658E-5,
    1.14129441738122E-5, 1.13310198717225E-5, 1.09679779254515E-5,
    0.000142353708796986, 0.000174643983507921, 0.000168230485136382,
    0.000175039474592756, 0.000169286514481601, 2.54021888323176E-5,
    9.98936386042353E-6, 9.41855015183268E-6, 9.74165993326622E-6,
    1.3258570003152E-5, 1.10705237392547E-5, 0.000117333722140816,
    0.000162720500118738, 0.000185043972553572, 0.000184174535777075,
    0.000146373456879942, 0.000127402231299929, 1.1854953188348E-5,
    1.11414958806797E-5, 5.71361756276303E-6, 5.69725210286728E-6,
    7.57947321435211E-6, 0.000106239140673664, 0.000181401084172403,
    0.000183763576612252, 0.000183547828211831, 0.000134499714419757,
    0.000141198789986318, 1.70175992228618E-5, 8.86908046213204E-6,
    8.92348003907378E-6, 1.52652245399514E-5, 3.84850635510477E-6,
    3.46418366401988E-5, 0.000182556362434547, 0.000216240916109494,
    0.000211232567178526, 0.000202844121706527, 0.000151826237066571,
    3.01110153466592E-5, 5.57854863597186E-6, 5.10970960807692E-6,
    6.19830753219746E-6, 9.02961133921801E-6, 8.66788795231223E-6,
    3.42487736428767E-5, 0.00016757733702273, 0.000186500989637556,
    0.000185129589539519, 0.000144098185555516, 0.000136982507256584,
    1.88120543492135E-5, 8.27185877820829E-6, 1.42856846707926E-5,
    1.31346379465996E-5, 9.33928544135933E-6, 6.75251701177791E-6,
    1.47999030924232E-5, 0.000186314244799197, 0.000167930908485708,
    0.000168038594903829, 0.000200542142054648, 0.000136212493005796,
    2.2276059893933E-5, 1.14172560464683E-5, 8.11431047472767E-6,
    8.13350824826142E-6, 1.02112104168833E-5, 7.32715089651084E-6,
    3.73417360063217E-5, 0.00017541390247907, 0.000191646428989796,
    0.000187150573227188, 0.00015954163918969, 0.000160618486418427,
    1.67109821720494E-5, 1.09763071341104E-5, 1.27696290690234E-5,
    9.23225565099539E-6, 5.45348019026755E-6, 5.63447761477831E-6,
    1.3915039249446E-5, 0.000154704710996008, 0.000190433311252732,
    0.000182312975750952, 0.00022483293795375, 0.000128056682607901,
    1.79623958291601E-5, 5.96164414575867E-6, 5.65775193490814E-6,
    4.87427033621036E-6, 5.08360828684861E-6, 5.64507569105814E-6,
    1.49656693785739E-5, 0.000207305174786015, 8.30963566020589E-6,
    8.83240241601226E-6, 9.79218066306115E-6, 6.93257319393915E-6,
    7.86607856654563E-6, 7.40632235748521E-6, 8.23877637987571E-6,
    8.28246713102136E-6, 8.42092819274844E-6, 8.22706798210133E-6,
    7.12679420606118E-6, 7.77608893440445E-6, 7.23155893739801E-6,
    1.07623403386682E-5, 1.08883387330942E-5, 1.02848156138217E-5,
    7.64911462227293E-6, 7.86970041641992E-6, 8.03295886851414E-6,
    5.80905794591016E-6, 6.73050384498216E-6, 9.89459913185093E-6,
    8.35874723483246E-6, 2.42956791569815E-5, 2.1511031630104E-5,
    2.06181417617334E-5, 2.0947654298088E-5, 1.81996272848841E-5,
    0.000154596236971603, 0.000167435381171399, 0.000170938973600765,
    0.000183471065687664, 2.20299664768237E-5, 8.18320130568198E-6,
    5.89166015422496E-6, 6.11135264003396E-6, 4.26469443746957E-6,
    5.42503414453847E-6, 6.68032472606053E-6, 7.46352992823557E-6,
    8.5097908535504E-6, 1.46332360016514E-5, 8.75614837465149E-6,
    1.06058888347418E-5, 2.33248599894066E-5, 6.83773026411856E-6,
    6.9360875211809E-6, 2.96675059911513E-5, 0.000148971898989248,
    0.000157039186347777, 0.000156940588527483, 0.000167588982694493,
    1.14511370440481E-5, 9.01961772540708E-6, 1.03101085410389E-5,
    1.28700950202054E-5, 7.1952987323348E-6, 6.90430222649351E-6,
    6.354060563884E-6, 0.000144531909527263, 0.000180275576809142,
    0.000181540349614601, 0.000181479634301162, 9.26624668742471E-6,
    8.7373580050192E-6, 9.00470955730525E-6, 6.92836665588625E-6,
    6.21923162617571E-6, 7.80625742713701E-6, 9.43550654816546E-6,
    1.25308166335016E-5, 0.000108566958847887, 0.000135072528970758,
    0.0001360656495883, 0.000135735808987537, 1.19347722619963E-5,
    6.78724920620006E-6, 9.98359276740667E-6, 1.12035913125319E-5,
    9.43790162484453E-6, 8.9546802390677E-6, 7.30449635827978E-6,
    9.75469031622197E-6, 0.000174772338209225, 0.000186223621508292,
    0.000179588002590334, 0.000197476285587807, 0.000176285569498328,
    1.15839381117546E-5, 6.94980617640822E-6, 7.41269288912933E-6,
    1.41459016993733E-5, 1.15829928418905E-5, 6.72849053424152E-6,
    1.72800779056569E-5, 0.000170387994969771, 0.000188213970477409,
    0.000197620314319065, 0.000221964589500811, 0.000199175117063518,
    7.53525052880131E-6, 9.84484729256211E-6, 7.58838445651814E-6,
    1.04735339127106E-5, 1.02324027750564E-5, 9.28195309809859E-6,
    8.67466806832367E-5, 0.000179728367052525, 0.000188502567358496,
    0.000190611367965602, 0.00020454262792609, 5.14376638803093E-5,
    5.99661622080278E-6, 8.1976923063347E-6, 8.94735249144547E-6,
    9.81825895765932E-6, 0.000150668881063169, 0.000194929913028259,
    0.000198791583459457, 0.000201702634400614, 0.000199602317642013,
    9.75270784632434E-6, 1.47935759364581E-5, 1.33718238344293E-5,
    1.67605099706368E-5, 1.24669296535988E-5, 0.000135169222685992,
    0.000175701009989375, 0.000178481447744706, 0.000178350784971035,
    9.41108260528611E-6, 8.0407602860787E-6, 7.1844660007466E-6,
    2.25415348884264E-5, 0.000189226556875085, 0.000194757790516201,
    0.000194793713282904, 0.000179509097404988, 5.84138741296378E-6,
    1.18972881517363E-5, 1.08120042583952E-5, 0.000161075518292425,
    0.000169169475223437, 0.00016989449995699, 0.00018018699461457,
    1.10006800414135E-5, 9.48646514087132E-6, 9.42790123393488E-6,
    1.18686536565591E-5, 1.65751955457006E-5, 3.76946055406422E-5,
    0.000191602521827494, 0.000194047824688985, 0.000192435077622202,
    0.000164150558123319, 0.000143431027308591, 1.0196312496284E-5,
    8.59480222723045E-6, 1.20713535065826E-5, 1.13573612992558E-5,
    3.01012856302691E-5, 0.000172041531982016, 0.000184140137845603,
    0.000190505529578442, 0.000198791202102316, 0.000206537533626424,
    1.88011889255614E-5, 1.53664440186827E-5, 9.91809855191449E-6,
    1.05413586103614E-5, 1.04304069097568E-5, 1.12040570150435E-5,
    0.000137371573433187, 0.000179431367951825, 0.000195447190902955,
    0.000196329184685834, 0.000192879432742, 1.14810068328191E-5,
    8.85134564506755E-6, 1.09654720707843E-5, 9.98084054733988E-6,
    8.61043849195495E-6, 0.000179696104847472, 0.000236006991677063,
    0.000219757614613027, 0.000222115092468106, 0.000179166550858339,
    2.14497090892979E-5, 1.39233476248226E-5, 1.55646531643257E-5,
    1.66714562555803E-5, 0.000171807719652834, 0.000232627670196006,
    0.000213607272114998, 0.000212181269943503, 0.000151823068479357,
    1.26957757111968E-5, 1.43425023386619E-5, 1.12292025717348E-5,
    1.25692944555698E-5, 0.000167548251205821, 0.000196320379303976,
    0.000206237452678408, 0.000213791641747558, 0.000209076459370999,
    6.39205457517384E-5, 5.87937332866721E-6, 7.12227083774528E-6,
    7.84954886524982E-6, 9.29530915324456E-6, 8.3968051344022E-6,
    1.27878228947974E-5, 7.82381515055282E-6, 0.000125377426254418,
    0.000181755597639771, 0.000187041276320814, 0.000187062741824012,
    2.3136004696335E-5, 2.45222054540092E-5, 2.61448866754171E-5,
    9.4985996761883E-6, 0.000114909217786928, 0.000183632013683129,
    0.000188025835007587, 0.000188299128339105, 1.35168464101074E-5,
    1.59199581832883E-5, 1.42547886043713E-5, 7.00156376243729E-6,
    6.62558537140024E-6, 8.06262037092864E-6, 8.62852125672651E-6,
    7.77353658195917E-6, 1.99193304088725E-5, 2.90709619117735E-5,
    2.82335438133689E-5, 2.71301161404601E-5, 1.90267276526926E-5,
    1.92940278041694E-5, 8.18819174997747E-6, 9.67953414979765E-6,
    9.85594944551089E-6, 1.14464181513026E-5, 1.30551523973067E-5,
    1.35921926326413E-5, 1.35151223139563E-5, 1.28639613116965E-5,
    2.83797388577453E-5, 3.32239488541743E-5, 3.22191320888227E-5,
    3.88081368217141E-5, 7.99155731062778E-5, 8.28536937741968E-5,
    8.42449430034073E-5, 1.22406632608382E-5, 1.19244390656088E-5,
    0.000101815171324154, 0.000204276148330749, 0.000185087398273051,
    0.000192142277056291, 0.00021402948692242, 4.04839976844212E-5,
    9.40313576516357E-6, 9.5133490905717E-6, 9.72666287291884E-6,
    5.89453681612172E-6, 4.81055140455391E-6, 5.92942448390965E-6,
    1.04915288755555E-5, 8.81912287399776E-6, 4.40796817048327E-5,
    0.000208791599192561, 0.000191409011246054, 0.000198896386358573,
    0.000221789841266304, 6.06232139986211E-5, 7.2247772930078E-6,
    9.33499554288948E-6, 9.67390448019106E-6, 8.3229416211266E-6,
    6.43604669326029E-6, 6.9661491566804E-6, 6.06833247404437E-6,
    6.77806109548577E-6, 2.37014517694929E-5, 0.000151355067540881,
    0.000148160581536294, 0.000152824491742515, 0.000168295327564102,
    0.00012017003073071, 1.05328678905131E-5, 9.43309166799686E-6,
    1.19835681382049E-5, 1.0339337642445E-5, 1.0615529186209E-5,
    8.45142166751238E-6, 9.34266329491774E-6, 7.12239938561719E-6,
    0.000135650255717206, 0.000179258133771805, 0.000183097382279143,
    0.000183333042752223, 1.0606537771759E-5, 8.57821721489041E-6,
    5.745876766734E-6, 8.50104575140536E-6, 1.18840636148347E-5,
    8.08636825249652E-5, 0.000178335601265803, 0.000182814943647756,
    0.000183171807907016, 3.0402235564564E-5, 9.09108309441164E-6,
    9.72815786687408E-6, 5.30534602232229E-6, 1.17164744356953E-5,
    -6.38340899015181E-6, -6.18543976795429E-6, -6.24580561290306E-6,
    -8.01896004457985E-6, -7.96393256714011E-6, -7.46077503222289E-6,
    -6.88869147946523E-6, -4.69672065694966E-6, -6.07272134554605E-6,
    -4.81074587109591E-6, -4.63922110913348E-6, -4.74874525032186E-6,
    -6.11554361864991E-6, -6.65395771262599E-6, -8.17164409631411E-6,
    -8.90237066714313E-6, -8.15612596102311E-6, -7.15228608118149E-6,
    -5.08414616408279E-6, -9.70911069473729E-6, -1.92374149494576E-5,
    -1.3635245067784E-5, -1.32829752963314E-5, -1.12419486114878E-5,
    -2.52563156443E-5, -6.3676676537403E-6, -8.29042649594303E-6,
    -5.74985219451162E-5, -9.83146440272397E-5, -9.32593527336332E-5,
    -9.25189625514671E-5, -6.45954155865403E-5, -5.17015739972714E-6,
    -5.52896790890385E-6, -7.08270538980533E-6, -8.29313359837472E-6,
    -1.36405256147616E-5, -7.26426223197179E-6, -5.37436966033258E-6,
    -4.73961827315292E-6, -0.000107602525696006, -0.000113109667915497,
    -0.000107599639391951, -0.000108183270614034, -3.84580610887323E-5,
    -1.52435578020956E-5, -1.44401199411066E-5, -1.26620141995657E-5,
    -1.88510525216585E-5, -1.33027254316073E-5, -8.7262240708872E-6,
    -0.00012578940950297, -0.000103499437879966, -0.000100495213306833,
    -9.88009037761125E-5, -5.09911647080099E-5, -7.13685742156097E-6,
    -6.83273319250493E-6, -7.2316034292789E-6, -7.87008605009384E-6,
    -9.46657094149083E-6, -6.9619231179589E-6, -1.04378934353006E-5,
    -1.00691202830119E-5, -1.73422298811448E-5, -0.000119132583851459,
    -0.000207360936403746, -0.000203705538132031, -0.000211384141096332,
    -7.16497791799887E-5, -1.19440574966607E-5, -1.02574936890057E-5,
    -9.37850583700642E-6, -7.68460121170897E-6, -1.10763840728035E-5,
    -6.46131450246317E-6, -9.25000786477894E-6, -1.21751790679318E-5,
    -6.93350714455283E-5, -0.000102036475324896, -0.000165230134493758,
    -0.000159688836386879, -0.000150196546876913, -4.59348089521248E-5,
    -1.1288384765936E-5, -5.86616964991636E-6, -9.95116449556288E-6,
    -7.22398114735333E-6, -8.67288193956946E-6, -1.26818767585621E-5,
    -0.000116817881614321, -0.000167902416616617, -0.000182987116705657,
    -0.000179771993474156, -0.000150492896239045, -9.51045715892242E-5,
    -1.95632466308315E-5, -1.20362165073256E-5, -5.4692286604702E-6,
    -8.13590590061654E-6, -7.70575872462703E-6, -1.31446256167786E-5,
    -9.3083201482102E-6, -9.10594029288963E-6, -6.01239985383504E-5,
    -8.86115082640276E-5, -0.00013475823522146, -0.000134518776767965,
    -0.000112664248834932, -8.59146202317202E-5, -9.40044389570666E-6,
    -1.02558927964578E-5, -8.23021923498398E-6, -9.13782574641912E-6,
    -9.70275936152425E-6, -8.87996705201286E-6, -5.72019285031065E-5,
    -9.90818908820685E-5, -0.000131477552284168, -0.000127307640674484,
    -0.000126500157047996, -7.7765789724233E-5, -8.5287957709401E-5,
    -6.59802286872896E-5, -3.82253246575156E-5, -6.53642690152356E-6,
    -6.37733646681391E-6, -9.78260414459615E-6, -7.29256027931285E-6,
    -6.86080801648244E-6, -1.10161816832764E-5, -0.000132548158665864,
    -0.000135470626698316, -0.000136803992251013, -0.000132844710199401,
    -5.5324463119071E-5, -1.13323467410252E-5, -1.52910840907905E-5,
    -1.16680303957447E-5, -7.83584739267695E-6, -7.96668835038872E-6,
    -8.81855770201449E-6, -4.89139562604141E-5, -0.000133034122386663,
    -0.000132884262976013, -0.000132721577995584, -8.11057439261302E-5,
    -1.69388129947534E-5, -1.33292747631475E-5, -9.00616883152697E-6,
    -8.30400847733418E-6, -2.05646233282631E-5, -9.61352743356155E-6,
    -1.05277983912367E-5, -1.36646737309232E-5, -0.000111479982673311,
    -0.000109606735439164, -0.000109435779568991, -6.23586290579086E-5,
    -1.53518346097812E-5, -1.22710922793822E-5, -6.50432165571121E-6,
    -8.85570948934099E-6, -8.4502360525971E-6, -8.40155108020576E-6,
    -9.51981014334514E-6, -0.000102050373845928, -0.000149197012326862,
    -0.000141487382559787, -0.000137800920390336, -0.000116701845731597,
    -4.94959884937465E-5, -8.64181310499091E-6, -7.91115302959908E-6,
    -8.2438135869922E-6, -2.74596430130876E-5, -9.76639770453178E-6,
    -8.5245547921219E-6, -7.93096000754115E-6, -9.89204723873253E-5,
    -7.4082874066958E-5, -0.000118949341681179, -0.000115742975450116,
    -9.84239993219516E-5, -5.27512835865653E-5, -1.03002701748716E-5,
    -1.12861882534127E-5, -7.02838227183154E-6, -5.89403541929232E-6,
    -5.70480891343103E-5, -9.83117382536026E-5, -0.000124225315300904,
    -0.000116425915119603, -0.000101525859974714, -0.000106541385130338,
    -1.18861718085244E-5, -5.93872493581622E-6, -8.06227570411751E-6,
    -7.8042027276739E-6, -6.21346346432807E-5, -0.000110314893969004,
    -0.000139179169550037, -0.000159871923906587, -0.000146740570101285,
    -0.000135334625413979, -9.70803276453698E-5, -5.26490193860778E-6,
    -1.24020493086722E-5, -5.25552655968676E-6, -4.48146756683935E-6,
    -5.1972888946854E-6, -3.8046919678161E-6, -5.11728182328601E-5,
    -0.000118851068188529, -0.000119421569049432, -0.000128421252602376,
    -0.000154087014107822, -0.000148809711162014, -0.000156587985449356,
    -9.83196186933653E-5, -9.4056859696581E-6, -1.03846459934005E-5,
    -1.04372438083658E-5, -1.16109761857082E-5, -9.21366997578013E-6,
    -6.93845260990026E-6, -9.75086307447044E-5, -0.000103742287991364,
    -0.000161584114373589, -0.00015903614999152, -0.000139927901063684,
    -4.9137867770653E-5, -9.75075580938731E-6, -1.10035298173292E-5,
    -9.03089616668909E-6, -1.20932260553404E-5, -1.33927105927424E-5,
    -9.13990125124233E-6, -0.000127587055738161, -0.000131464134568945,
    -0.000130875116961978, -7.86326699378058E-5, -2.26899222277067E-5,
    -7.12873327143919E-6, -8.28524952220046E-6, -1.2099815014312E-5,
    -4.47488400537254E-6, -6.23798016607255E-6, -8.60660427641757E-6,
    -9.41302233945687E-6, -8.63655464188522E-5, -0.000105193229960805,
    -0.000101873749720304, -0.000101790627401466, -6.66445013947076E-5,
    -6.02835441795078E-6, -5.32730610880539E-6, -7.15445542276358E-6,
    -7.61646231127312E-6, -0.000130053349001777, -9.90229581352028E-5,
    -9.66618343417993E-5, -9.50456969728395E-5, -3.72682440856382E-5,
    -9.42817200027075E-6, -9.32948837665549E-6, -9.10598004701104E-6,
    -9.18229106137903E-6, -1.0106496944033E-5, -4.59179238851403E-5,
    -9.44344079925657E-5, -9.01360498785687E-5, -8.98305303695483E-5,
    -3.92827424021999E-5, -1.15301907470917E-5, -1.92146989893495E-5,
    -9.58515634922948E-6, -9.4447043946406E-6, -1.02085503676262E-5,
    -7.00839585946823E-5, -0.000103537557790755, -0.000113169167000681,
    -0.00011951318895003, -7.50949655307418E-5, -3.45693801730495E-5,
    -1.07952006081798E-5, -7.33245488962192E-6, -8.12186499259137E-6,
    -6.75100441504174E-6, -7.63964718923457E-5, -7.61266179920219E-5,
    -0.000106795762137278, -0.000102787805214802, -0.000100863417681924,
    -6.2795520981109E-5, -1.57471302603946E-5, -8.21383375507586E-6,
    -1.01105776009431E-5, -9.23071917005134E-6, -9.50673389865762E-6,
    -7.3775327977167E-6, -8.40780659398485E-5, -0.000133063196674287,
    -0.000126719002308489, -0.000126306602136303, -0.000116791354024818,
    -6.81624615794423E-5, -7.47842749364729E-6, -9.07096995048303E-6,
    -1.04236196925696E-5, -1.13375356801328E-5, -1.3035672037833E-5,
    -0.000117240349115531, -9.23604406255084E-5, -0.00015148620017929,
    -0.000151194769096583, -0.000139281769446007, -7.20058096836945E-5,
    -1.60790672813129E-5, -8.18785660718519E-6, -5.77655690359957E-6,
    -6.89494317023971E-6, -1.09214217864509E-5, -0.000102383483441712,
    -0.000119215231019789, -0.000152983968005343, -0.000151300616465195,
    -0.000136188139771676, -5.3454354238843E-5, -1.74112082420183E-5,
    -1.39723159407891E-5, -5.94493368004468E-6, -6.52463353139518E-6,
    -3.23180936078005E-6, -4.99663183620293E-5, -0.000113682808279283,
    -0.000139848645901661, -0.000151517786343366, -0.000149053848600058,
    -0.000127979561586396, -2.57251035445744E-5, -9.41077672851492E-6,
    -4.99120617115597E-6, -6.00108641591208E-6, -8.21112352080357E-6,
    -7.98034898923314E-6, -5.81416627729947E-5, -8.99466654291615E-5,
    -0.00011959214524599, -0.000126796905942358, -0.000108051899878355,
    -5.85900116721383E-5, -3.8995960588854E-5, -1.3008940835938E-5,
    -8.14509205377657E-6, -7.55462445867775E-6, -7.11276145493855E-6,
    -7.04819228156059E-6, -2.13622179671612E-5, -8.78057924416942E-5,
    -0.000163392903143893, -0.000163162975760932, -0.000170134819603656,
    -5.89718107817614E-5, -5.67406251816368E-5, -1.13906468033177E-5,
    -7.44198253832803E-6, -7.38874311263026E-6, -1.0832970498744E-5,
    -7.27150532216456E-6, -6.2768431589933E-5, -0.000117940830649288,
    -0.000139958138299287, -0.000135095924770699, -0.000125377818480477,
    -0.00014678126623579, -2.07570262794706E-5, -1.22143805269204E-5,
    -5.25265875473229E-6, -4.84956190411096E-6, -4.78356976632831E-6,
    -4.57799952547848E-6, -1.77953828036531E-5, -8.46937234496623E-5,
    -0.000148939854143091, -0.00015253130136727, -0.000159137936057821,
    -0.000116000000994039, -2.115470171797E-5, -7.80294910845377E-6,
    -1.23716533820859E-5, -5.68849260577528E-6, -5.32896444492429E-6,
    -9.94613037897698E-6, -2.09268520123957E-5, -0.000111016517947474,
    -8.40712340210654E-6, -1.07874493013384E-5, -9.0440720053169E-6,
    -8.88984106968718E-6, -7.78570101511379E-6, -6.97817502580117E-6,
    -1.30112990315916E-5, -4.67016730143145E-6, -4.28270797100859E-6,
    -1.03394125645678E-5, -7.23020270989166E-6, -5.36100131188597E-6,
    -7.50358869710114E-6, -5.89392189821691E-6, -7.59113914848407E-6,
    -6.96164239591637E-6, -7.09498348150116E-6, -7.44031359786709E-6,
    -7.34656555870977E-6, -8.000346298158E-6, -5.70044517609097E-6,
    -5.69430339959273E-6, -5.52597463805485E-6, -2.40776635208592E-5,
    -1.58118644542712E-5, -1.74294555120647E-5, -1.61606790152087E-5,
    -3.53399075015325E-5, -0.000112171628715413, -0.000108324490697361,
    -0.000106989623120982, -6.71856375338191E-5, -7.93196697847358E-6,
    -4.55781376812244E-5, -1.35746522805372E-5, -9.49705913466203E-6,
    -7.47678083621094E-6, -2.19279525175762E-5, -6.40712350505209E-6,
    -6.96178843489187E-6, -8.80738591950372E-6, -6.49133318950107E-6,
    -6.0811120528284E-6, -6.3986964192441E-6, -6.80348701557454E-6,
    -5.24983652370387E-6, -7.04916474779538E-6, -4.83182198849207E-5,
    -0.000101107091736156, -9.50515619108406E-5, -9.45325065244011E-5,
    -8.00857756221376E-5, -1.64929843246208E-5, -7.27291179912561E-6,
    -6.76048182244631E-6, -1.11211738087871E-5, -7.3505737243689E-6,
    -1.30408567418807E-5, -7.35753077940265E-6, -9.67209122186429E-5,
    -8.87683508890333E-5, -8.90074127423184E-5, -8.85248087865261E-5,
    -3.36943738941594E-5, -8.68409263815628E-6, -6.62690832290232E-6,
    -5.70123374962376E-6, -6.64563685963496E-6, -1.36656830070401E-5,
    -1.36298504583948E-5, -1.12809642917767E-5, -6.11238312806364E-5,
    -5.24093463070458E-5, -5.18099581207047E-5, -5.3099403778227E-5,
    -1.03146488999336E-5, -7.77058535499485E-6, -8.52184286018954E-6,
    -8.16413361971398E-6, -6.81508875658748E-6, -9.70069700490432E-6,
    -9.30004874062175E-6, -1.17056668767971E-5, -9.95780659735836E-5,
    -0.000180111295872069, -0.000178313601343483, -0.000186375411690305,
    -0.000113348914935778, -7.47961074636149E-6, -5.84899879755561E-6,
    -7.050088419815E-6, -6.3538629180251E-6, -6.25915244346789E-6,
    -8.08323896911386E-6, -2.23995034730879E-5, -9.0233420583228E-5,
    -0.000123555210868439, -0.00011865410923019, -0.000125790849895555,
    -9.3917909241736E-5, -5.774253855824E-6, -5.53864877947572E-6,
    -1.49578683801481E-5, -6.56991257583894E-6, -8.36120325040224E-6,
    -7.97387218015783E-6, -9.58475540276453E-5, -8.19992573350198E-5,
    -9.86971361560109E-5, -9.48668365713517E-5, -0.000113602051369762,
    -0.000103567058808864, -1.35273157663208E-5, -1.11860147449598E-5,
    -8.79583882099343E-6, -8.56414916164345E-6, -8.26989804544984E-5,
    -0.00016122824633582, -0.000137882678734717, -0.000137961477897818,
    -0.000115946190847001, -1.35063512024158E-5, -1.26703621927736E-5,
    -1.35786095082467E-5, -1.46994096356889E-5, -8.20072549464703E-6,
    -8.91587468174406E-5, -7.37182080698916E-5, -7.16370295273752E-5,
    -7.40114658565836E-5, -1.27700604243985E-5, -1.49925910574514E-5,
    -9.80155512124325E-6, -3.90554219238789E-5, -9.00670230693954E-5,
    -8.35783426107729E-5, -8.27766698221549E-5, -7.02790081376336E-5,
    -6.25129922714006E-6, -6.472434987707E-6, -7.32205881728722E-6,
    -6.93563558815604E-5, -6.5005475881964E-5, -6.5860250764773E-5,
    -5.3571902734093E-5, -1.05236464302889E-5, -1.02999321070529E-5,
    -1.03027283153391E-5, -1.36282042803548E-5, -1.20559673168855E-5,
    -5.48737351284592E-5, -8.51776067374289E-5, -0.000161335674627783,
    -0.000156666741705255, -0.000140173307359071, -0.00010313491452512,
    -8.04699482685978E-6, -7.15885565160342E-6, -1.25380499981326E-5,
    -1.49160214725769E-5, -5.35136057149177E-5, -7.07728046561004E-5,
    -0.000150136979917189, -0.000145463756096896, -0.000131279697529669,
    -9.22954022067666E-5, -1.58924601135135E-5, -1.07652105585337E-5,
    -1.04175175445922E-5, -1.0932665476848E-5, -1.0433640078416E-5,
    -1.28355025543536E-5, -7.71293705302382E-5, -0.000106451611427554,
    -0.00010241451360783, -0.000107000695469219, -8.13013295635425E-5,
    -1.23120440787114E-5, -9.13857086651006E-6, -1.1736342313031E-5,
    -1.0523589733393E-5, -1.06236034730625E-5, -0.000106992556779914,
    -0.000178170088275536, -0.000167307249519804, -0.000167843845191329,
    -0.000126546021844419, -5.47583688059934E-5, -2.00666943903624E-5,
    -1.06923920021332E-5, -1.40478845216997E-5, -0.000114286709010509,
    -0.000227274389751229, -0.00022751907798943, -0.000227195425336159,
    -0.000138422950778249, -6.02421494262996E-5, -1.2656508823119E-5,
    -1.85836819487128E-5, -1.03451671089862E-5, -9.05744766734938E-5,
    -9.36553831565711E-5, -0.00013856851932578, -0.000133176690271422,
    -0.000131300311704605, -8.59981869706198E-5, -7.31929343070065E-6,
    -2.3546360585123E-5, -6.65900832130499E-6, -5.9513535604983E-6,
    -8.79410258011095E-6, -9.41825977214702E-6, -8.79428502324887E-6,
    -0.00011758738332233, -8.62150013422318E-5, -8.26174465877172E-5,
    -8.29372520123201E-5, -1.85617675164594E-5, -1.67678682617389E-5,
    -2.62321167139069E-5, -1.31840072020076E-5, -0.000119960374736014,
    -8.75208419145712E-5, -8.43413919764454E-5, -8.26550353970565E-5,
    -6.4663499281838E-5, -1.52998907492649E-5, -7.08975714082045E-6,
    -7.14542835483761E-6, -7.38036183029689E-6, -1.25640802255692E-5,
    -7.60370074149366E-6, -8.64400446713304E-6, -1.70690335035989E-5,
    -1.65514847217619E-5, -1.57661322182435E-5, -1.72598565377992E-5,
    -1.38003187942966E-5, -1.62698622211159E-5, -1.0816150968738E-5,
    -2.05217130571832E-5, -9.17844924990449E-6, -5.98601565197344E-6,
    -6.76328624505604E-6, -6.66515698374444E-6, -1.0384358667011E-5,
    -1.18953708459441E-5, -1.60033245556614E-5, -2.65402456497354E-5,
    -2.80920778478361E-5, -2.73606002925807E-5, -3.02827480845473E-5,
    -3.17346993025537E-5, -3.50223008396849E-5, -3.81943577444453E-5,
    -1.32805215510047E-5, -0.000105418908755065, -0.000196324857889353,
    -0.00020131063518464, -0.00019881695135449, -0.000204446027510452,
    -9.7899954048346E-5, -8.43382309822106E-6, -7.79749028591011E-6,
    -7.34899995269359E-6, -1.54682303584911E-5, -9.25374941501581E-6,
    -6.27249329196924E-6, -1.14917990279835E-5, -7.37711613922583E-6,
    -6.80089910640623E-5, -0.000140168972063061, -0.00018419557780116,
    -0.000180437376046296, -0.000175530975128285, -8.28932998480502E-5,
    -1.06304542866216E-5, -1.34327491536807E-5, -7.67164425510918E-6,
    -8.66643160544043E-6, -9.84235026849154E-6, -5.08622166400611E-6,
    -4.7745695756829E-6, -6.09015457797284E-6, -4.39739069990154E-5,
    -6.39397274321716E-5, -0.000120674507718564, -0.00011724138414223,
    -0.00010074141114113, -4.51693782793887E-5, -8.62313145188628E-6,
    -7.20042310411897E-6, -7.80246377580908E-6, -7.64451225875337E-6,
    -8.64643349864813E-6, -1.46329155812771E-5, -1.63009161968269E-5,
    -8.01250220710031E-6, -0.00010149952249082, -8.96959822940957E-5,
    -8.76335655665426E-5, -8.63930010376524E-5, -2.08598733770774E-5,
    -8.07691818197857E-6, -7.07121792789965E-6, -1.01035343950162E-5,
    -6.35521284323931E-6, -8.82851737406733E-5, -9.20604219129436E-5,
    -8.80200887452686E-5, -8.83380897665321E-5, -7.72315303097289E-5,
    -8.93081367927715E-6, -8.39761414370676E-6, -1.08336446535989E-5,
    -1.03663704462573E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 20.0, 22.0, 13.0, 18.0,
    0.0, 0.0, 43.0, 71.0, 75.0, 76.0, 35.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0,
    0.0, 39.0, 86.0, 87.0, 68.0, 41.0, 21.0, 21.0, 18.0, 19.0, 13.0, 0.0, 41.0,
    78.0, 85.0, 76.0, 38.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 1.0, 22.0, 49.0,
    73.0, 85.0, 55.0, 56.0, 13.0, 4.0, 0.0, 0.0, 8.0, 0.0, 0.0, 9.0, 40.0, 88.0,
    91.0, 91.0, 59.0, 41.0, 9.0, 0.0, 0.0, 0.0, 0.0, 10.0, 39.0, 60.0, 92.0,
    88.0, 66.0, 53.0, 16.0, 10.0, 0.0, 0.0, 0.0, 15.0, 0.0, 0.0, 22.0, 80.0,
    95.0, 113.0, 95.0, 55.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 38.0, 76.0, 100.0,
    101.0, 97.0, 105.0, 81.0, 61.0, 34.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 39.0,
    66.0, 84.0, 66.0, 64.0, 7.0, 10.0, 8.0, 0.0, 0.0, 0.0, 21.0, 43.0, 72.0,
    74.0, 60.0, 18.0, 10.0, 0.0, 0.0, 27.0, 0.0, 4.0, 10.0, 57.0, 88.0, 91.0,
    58.0, 16.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 53.0, 75.0, 111.0, 95.0, 87.0,
    45.0, 0.0, 0.0, 0.0, 24.0, 0.0, 0.0, 0.0, 45.0, 88.0, 106.0, 106.0, 91.0,
    54.0, 5.0, 7.0, 0.0, 0.0, 28.0, 76.0, 102.0, 130.0, 116.0, 68.0, 15.0, 0.0,
    0.0, 0.0, 52.0, 67.0, 73.0, 107.0, 84.0, 57.0, 43.0, 0.0, 10.0, 0.0, 0.0,
    0.0, 0.0, 40.0, 48.0, 64.0, 81.0, 90.0, 69.0, 68.0, 43.0, 0.0, 3.0, 4.0, 8.0,
    0.0, 0.0, 36.0, 72.0, 96.0, 109.0, 78.0, 45.0, 0.0, 7.0, 0.0, 7.0, 12.0, 0.0,
    40.0, 85.0, 73.0, 58.0, 41.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 30.0, 86.0,
    91.0, 93.0, 41.0, 0.0, 0.0, 0.0, 0.0, 24.0, 77.0, 95.0, 95.0, 36.0, 0.0, 0.0,
    0.0, 0.0, 2.0, 32.0, 74.0, 87.0, 86.0, 47.0, 9.0, 13.0, 0.0, 0.0, 3.0, 61.0,
    82.0, 102.0, 74.0, 61.0, 31.0, 9.0, 0.0, 0.0, 0.0, 56.0, 98.0, 130.0, 120.0,
    95.0, 42.0, 16.0, 0.0, 3.0, 0.0, 0.0, 0.0, 58.0, 73.0, 111.0, 97.0, 76.0,
    46.0, 0.0, 0.0, 8.0, 7.0, 16.0, 42.0, 78.0, 96.0, 113.0, 83.0, 44.0, 17.0,
    0.0, 0.0, 0.0, 6.0, 35.0, 69.0, 95.0, 100.0, 89.0, 63.0, 28.0, 22.0, 0.0,
    0.0, 0.0, 19.0, 62.0, 75.0, 106.0, 76.0, 62.0, 35.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    50.0, 72.0, 94.0, 113.0, 102.0, 64.0, 35.0, 7.0, 0.0, 0.0, 0.0, 0.0, 14.0,
    70.0, 72.0, 102.0, 77.0, 63.0, 47.0, 5.0, 0.0, 0.0, 6.0, 0.0, 35.0, 52.0,
    90.0, 117.0, 99.0, 75.0, 44.0, 9.0, 0.0, 0.0, 0.0, 0.0, 13.0, 74.0, 74.0,
    106.0, 80.0, 79.0, 52.0, 0.0, 15.0, 0.0, 0.0, 0.0, 13.0, 65.0, 0.0, 4.0, 0.0,
    0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 29.0, 40.0, 36.0, 28.0, 81.0, 76.0, 89.0,
    42.0, 0.0, 26.0, 14.0, 0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 43.0, 68.0, 81.0, 81.0, 43.0, 20.0, 0.0, 0.0, 8.0, 0.0, 10.0, 0.0,
    54.0, 97.0, 85.0, 69.0, 27.0, 0.0, 0.0, 0.0, 0.0, 9.0, 13.0, 8.0, 67.0, 87.0,
    87.0, 72.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 67.0, 67.0, 88.0, 53.0,
    40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 68.0, 86.0, 110.0, 74.0, 58.0, 0.0,
    0.0, 11.0, 0.0, 0.0, 0.0, 32.0, 91.0, 115.0, 120.0, 82.0, 29.0, 13.0, 6.0,
    0.0, 0.0, 63.0, 65.0, 95.0, 82.0, 64.0, 37.0, 21.0, 25.0, 12.0, 0.0, 50.0,
    83.0, 86.0, 69.0, 8.0, 19.0, 0.0, 17.0, 81.0, 86.0, 88.0, 42.0, 0.0, 0.0,
    0.0, 82.0, 88.0, 95.0, 60.0, 7.0, 3.0, 3.0, 22.0, 11.0, 40.0, 80.0, 73.0,
    92.0, 57.0, 38.0, 0.0, 0.0, 22.0, 14.0, 31.0, 94.0, 89.0, 110.0, 80.0, 59.0,
    18.0, 6.0, 4.0, 16.0, 5.0, 8.0, 64.0, 94.0, 126.0, 79.0, 61.0, 13.0, 0.0,
    10.0, 5.0, 8.0, 54.0, 57.0, 90.0, 77.0, 85.0, 59.0, 20.0, 6.0, 22.0, 59.0,
    57.0, 85.0, 64.0, 36.0, 29.0, 11.0, 19.0, 6.0, 78.0, 97.0, 110.0, 107.0,
    75.0, 38.0, 0.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 28.0, 81.0, 108.0, 110.0,
    65.0, 22.0, 29.0, 21.0, 28.0, 72.0, 82.0, 93.0, 31.0, 14.0, 0.0, 0.0, 0.0,
    12.0, 0.0, 0.0, 30.0, 43.0, 40.0, 27.0, 23.0, 29.0, 6.0, 15.0, 0.0, 0.0, 0.0,
    0.0, 4.0, 17.0, 52.0, 65.0, 49.0, 73.0, 83.0, 87.0, 76.0, 37.0, 11.0, 29.0,
    50.0, 81.0, 88.0, 60.0, 43.0, 0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 9.0, 0.0, 34.0,
    56.0, 79.0, 96.0, 77.0, 40.0, 5.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0,
    75.0, 115.0, 122.0, 98.0, 39.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 14.0, 0.0,
    52.0, 80.0, 97.0, 81.0, 43.0, 0.0, 0.0, 2.0, 0.0, 27.0, 75.0, 80.0, 81.0,
    42.0, 0.0, 0.0, 7.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.0, 9.0, 6.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 5.0, 2.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 23.0, 23.0, 20.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 21.0, 21.0, 22.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 21.0, 26.0, 25.0, 17.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 21.0, 12.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 16.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 14.0, 15.0, 16.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 15.0, 14.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 8.0,
    8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 20.0, 22.0, 21.0, 11.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 17.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 12.0, 3.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0,
    22.0, 24.0, 19.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 11.0,
    23.0, 18.0, 21.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 21.0,
    21.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 14.0, 14.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 3.0, 3.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 13.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 7.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0,
    15.0, 14.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 0.0, 19.0, 20.0, 16.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 14.0, 23.0, 23.0, 19.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 10.0, 24.0, 25.0, 25.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 18.0, 21.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    18.0, 18.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 15.0, 14.0,
    12.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 22.0, 23.0, 9.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 24.0, 25.0, 25.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    14.0, 13.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    10.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0, 23.0, 24.0, 9.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 26.0,
    35.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 18.0, 15.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 5.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0,
    31.0, 32.0, 32.0, 12.0, 0.0, 0.0, 0.0, 0.0, 9.0, 25.0, 25.0, 25.0, 22.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 17.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
    20.0, 22.0, 21.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    15.0, 21.0, 21.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 12.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0,
    19.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 20.0, 20.0,
    20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 20.0, 20.0, 20.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 22.0, 41.0, 42.0, 45.0, 23.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 36.0, 37.0, 21.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 30.0, 45.0, 45.0, 27.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 41.0, 43.0, 25.0, 27.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 18.0, 34.0, 34.0, 39.0, 37.0, 17.0, 17.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 13.0, 25.0, 25.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 27.0, 30.0, 30.0, 31.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    20.0, 21.0, 21.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 32.0, 42.0,
    42.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 18.0, 43.0, 44.0,
    25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 38.0, 40.0, 21.0, 20.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 17.0, 25.0, 54.0, 47.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 15.0, 18.0, 25.0, 42.0, 30.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 5.0, 20.0, 40.0, 41.0, 25.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0,
    23.0, 24.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 19.0,
    19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 18.0, 19.0, 19.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 16.0, 18.0, 18.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 28.0,
    36.0, 38.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 16.0, 36.0, 37.0, 19.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 28.0, 38.0, 40.0, 21.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 10.0, 17.0, 46.0, 49.0, 30.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    8.0, 19.0, 39.0, 43.0, 25.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 35.0,
    50.0, 50.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 38.0, 44.0, 27.0,
    26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 38.0, 39.0, 42.0, 19.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 29.0, 36.0, 29.0, 19.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 15.0, 20.0, 47.0, 50.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 21.0, 22.0, 22.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 22.0, 22.0, 22.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 16.0, 20.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 10.0, 14.0, 14.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    21.0, 44.0, 48.0, 49.0, 29.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 40.0,
    40.0, 43.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 17.0, 33.0, 33.0, 19.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 35.0, 40.0, 40.0, 22.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 17.0, 17.0, 17.0, 0.0, 0.0, 0.0, 0.0, 19.0, 19.0, 19.0, 19.0, 0.0,
    0.0, 0.0, 15.0, 16.0, 16.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 47.0,
    47.0, 30.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 41.0, 42.0, 31.0, 28.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 14.0, 34.0, 34.0, 37.0, 18.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 21.0, 41.0, 48.0, 49.0, 24.0, 0.0, 0.0, 0.0, 0.0, 19.0, 24.0, 52.0,
    52.0, 31.0, 0.0, 0.0, 0.0, 0.0, 14.0, 17.0, 43.0, 43.0, 25.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 17.0, 17.0, 17.0, 0.0, 0.0, 0.0, 0.0, 10.0,
    17.0, 17.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 7.0, 27.0, 39.0, 40.0, 25.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 21.0, 39.0, 39.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 32.0, 33.0, 22.0, 14.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 13.0, 19.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    17.0, 18.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.99661184214039E-5,
    -4.83883903823533E-5, 0.000111250382886763, -1.38552786537679E-5,
    -5.16296999956353E-5, -7.40426697644171E-5, 0.000119886559501568,
    8.96236912826785E-5, -6.0369263345886E-5, -5.63553171532901E-5,
    -1.8260517791725E-5, 1.92756614186633E-5, 6.25483062035625E-5,
    -5.48237564485524E-5, -5.37067074220348E-5, 0.000146745259458647,
    -3.76673572457482E-5, 5.8789122922614E-5, -8.13989733515558E-5,
    -6.29443670838187E-5, 4.14719822517947E-5, 0.000116024655986205,
    0.000133075140261335, 2.12381210716772E-5, -0.00036925404056195,
    -2.0318185994889E-5, 3.84534300583143E-5, -0.000143542088617853,
    -5.02488851527272E-5, -2.9113764593172E-5, -0.000424225002995508,
    -0.000429881201458318, -1.23175444475402E-5, -5.43750956644272E-6,
    -8.74986554676105E-5, 5.65510199964231E-5, -0.000168975876244208,
    -3.90080583770316E-5, 0.000135775369029094, -1.21492657559227E-5,
    0.000148130116486164, -0.000213907805885652, -0.000152987674717695,
    -0.000175162241275013, -0.00051380306136938, -0.00018383199781581,
    -0.000185548985474961, 6.21647934901723E-5, -0.000146509561850953,
    -0.000214894099919053, 0.000106266959201533, -0.000189411718298229,
    -0.000155729772077741, -1.07355220302602E-5, -0.0001659989270138,
    -0.000914012429685098, 4.09168247459765E-6, 0.000157653141087482,
    2.01958109095454E-5, 6.3501064692521E-5, -7.04563741843465E-5,
    -7.66119681902442E-5, 0.000162398081524288, 4.97868317207596E-5,
    0.000140282541099157, 0.000157826141440218, 0.000217357416509922,
    -0.000367792426931194, 0.00238214334419193, 0.00251408761224393,
    -6.42658076040324E-5, 9.46787970235643E-5, -4.92416690343127E-5,
    -7.89829017748042E-5, 0.000171491983454404, -2.25306686529296E-5,
    -0.000100015255591063, -3.4011024954528E-5, 0.00014489189701326,
    2.4719944465674E-5, 0.000139253276011336, -0.000353606598035064,
    -0.000675411132420654, -0.000145038388238393, -0.000157198743142563,
    3.96750800744602E-5, -9.66300677017648E-5, -8.90294722167508E-5,
    0.000123054219385912, -0.00011738335956293, 0.000170628863719401,
    -3.7195690247264E-5, -0.000235268682423777, 4.89548570034522E-5,
    -0.00069827907696995, -0.00186238116060834, -0.000237930253051631,
    0.000295609899901909, 6.19793484839828E-5, 0.000208706146011918,
    5.4736815159999E-5, -0.000176413087582748, -5.34729200342001E-5,
    -7.48622686531363E-5, 0.00014537483376658, 6.03162713004863E-6,
    0.000237226268858162, -0.000323751998368016, -0.00114223129818175,
    0.00255112829862107, -5.95663021151186E-5, -0.000156542618122559,
    -0.000154735174838469, 9.72828547874436E-6, -2.86988128251813E-5,
    -0.00014499827574117, -0.000101770766298747, 7.9100400454904E-5,
    -0.000104169661300832, -0.000290921802552061, 0.00112720375114643,
    0.00222790549688484, -0.000129804873490489, 0.00065509651680939,
    -0.00074556844219838, -6.00554060002804E-5, -8.32561555751829E-5,
    -0.000153402204107869, 2.797201062408E-5, -5.62066556571904E-5,
    -0.000140661951204936, 0.000172440596331288, -6.13813314025921E-5,
    7.65737586494643E-5, 5.12506869382042E-5, -0.00117651219020627,
    0.000330941440822158, -0.000220195215711117, 0.00022583074338696,
    0.000304156173726459, 3.70920560755515E-5, -9.04847729993007E-5,
    0.000138953518935461, 0.000163858386422068, 0.000243403958957621,
    -0.000468351054653535, 0.00253712212928492, -0.000238984804599238,
    -0.000146666212780839, -9.05723935967939E-5, -2.39510236557474E-5,
    -0.000279886993133393, -7.87743107085691E-6, 6.50707985515324E-5,
    -1.20271927620357E-5, 0.000217276396271568, 6.09852645749566E-5,
    -0.000426921130088973, 0.00255786663237921, -6.08663561453341E-5,
    -0.000135919762543708, 0.00015978557875916, -7.61029171816642E-5,
    -8.61564789108961E-5, 2.46949486097114E-5, 0.000125201726946759,
    0.000216341232568039, -1.1109812194549E-5, -0.000276550226350726,
    0.000377459824284565, -0.000364381669681732, -0.00106445977306311,
    3.55125531840176E-5, 2.60078046787412E-5, 2.11096542197098E-5,
    -0.000371123019074766, -3.48243017940637E-5, -3.44146801948779E-5,
    5.69873365865647E-5, 0.000111108646488518, -5.76987380213348E-5,
    -3.1552075985861E-5, -2.44951596436551E-5, -0.00110631527075522,
    0.000340932565326582, -4.60563866970758E-5, -0.000178932478785907,
    7.93024627713332E-6, 9.76376808135516E-5, 3.08182223431743E-5,
    -6.5590790775827E-5, 0.000328687037271832, -0.000306637649786796,
    -0.000495825872436841, 0.000184028160416892, 2.10979218914289E-7,
    9.2643722690494E-5, 1.29601782636632E-5, 3.23444409510625E-5,
    -2.68075967578758E-5, 7.12522801598802E-5, 0.00018749698340599,
    -0.000378944021064318, 0.000369005584589379, 0.000159727109408348,
    -0.00175426841665234, -9.2314281571934E-5, -0.000128727032838198,
    1.77957949708725E-5, 2.23268844089127E-5, -4.04119200064145E-5,
    4.73756566194408E-5, -0.000113485149968275, 6.61905624312549E-5,
    6.90280152354599E-5, -0.000476013946499906, 0.00168036000258004,
    -7.55051574509867E-5, -0.0011960836635833, -0.0010985651354727,
    -6.26218441123006E-5, -0.000107888472701622, 0.000190211682700487,
    4.53336836382378E-5, 2.58453521908546E-5, 1.62112657523091E-5,
    -8.60034466909409E-5, -0.000106955933270376, 3.15450281876396E-5,
    -0.000359251638334574, -0.00204139385365953, 0.000532758134878425,
    9.29898923342644E-5, -0.000130178257603278, -8.28261726445385E-5,
    -0.000133316222981108, 0.000380746763692719, -8.63731972673901E-5,
    0.000134128356516981, -0.000222786713629402, -0.000298826518452565,
    0.00204372705256637, -0.000427889095700262, -1.30936141322072E-5,
    -6.89425252204435E-5, -0.0001697947593588, 0.000129928919061717,
    7.3157006437499E-5, 4.51530869450085E-5, 1.21151081856888E-5,
    -4.9046000793027E-5, -0.000164728126173668, -6.43142972373946E-5,
    -0.000488176616350041, -0.00131429468701181, 4.95474535028276E-5,
    9.11410192459967E-6, -6.76866912413478E-5, 3.25274057415789E-5,
    -0.000137537842608277, -0.000127364759518887, 8.49123597398495E-5,
    -0.000375741375288804, -0.00073043682810548, 0.000126552907517905,
    0.00023844800733505, -6.02925815223142E-5, 6.47321147780342E-6,
    8.07332977373374E-5, -6.65173509204399E-6, 0.000128251530438613,
    0.000143848736180973, -0.000316652331403533, 0.00107408865080654,
    -0.000130573785673232, -0.00021839175114402, -2.36282797153674E-5,
    7.7295034121368E-5, 7.250927185224E-5, 0.00020236402438748,
    0.000149317198555055, -8.74516946853483E-5, 0.00105083930759351,
    0.00100706421743763, -0.000643812498323467, -8.01407955663684E-5,
    -5.45760274970834E-5, 7.73074355265911E-5, 7.45374788043339E-5,
    0.000157235224822037, -3.70425666092791E-5, -0.000294009256240681,
    0.000373762351031016, -0.000660738469124902, -0.00122429923372736,
    -0.000218320904572873, -9.39582233920309E-5, 6.52320795410603E-5,
    -7.41785783415947E-5, 9.22224069652566E-5, 0.000107374921498922,
    0.000100968101278929, -0.000112153784963324, -0.000107333677685605,
    0.00110364637728335, -0.000556549543163517, -0.00136838790426048,
    3.45687168398042E-5, -6.43348995842805E-5, -7.42646097305854E-5,
    0.000178639993974416, -0.000123311746679814, -7.09551308531805E-5,
    -3.92938333555321E-5, 0.000218621173390307, -0.000344724457858257,
    -0.00123625637639267, 0.0018079284032108, -0.000146736584145176,
    0.000132915991727913, -3.90320332768504E-5, -0.000124808880574441,
    -0.000200571760342878, 0.000169386916774506, -0.000113807998056499,
    -0.000142316596918384, -0.000265436138249938, -0.00132482696191435,
    0.00147083166811172, -0.000233838171442275, -0.000146663288237695,
    0.000124122669889836, 0.000187021143237408, -1.21406978386868E-5,
    -0.000117473517972892, 4.2672583359238E-5, 0.000243893125795331,
    -0.000420125380830518, 0.00182389395569224, 0.000145215221979295,
    -0.000158022696214918, -2.68062536317087E-5, -5.9325580771086E-5,
    7.63467008338353E-5, 8.81196086447828E-5, 9.8033582829148E-5,
    0.00014229582493884, 3.86162571735068E-5, 6.38348488575017E-5,
    -0.000390586767039346, -0.00120706543419969, 0.00208659103073508,
    -0.000681740132883372, -0.00017670149454041, 0.00020625962307432,
    0.000137498287168978, 1.11690244167932E-5, -7.40512455645516E-5,
    -8.31297833371222E-5, 5.29440873899594E-5, 0.000121298499108464,
    -0.000218302232765283, 0.00166472782099895, 0.00145905963563188,
    -0.000917303918996183, 6.01959745695163E-5, -7.40395368502035E-6,
    2.26791190986166E-5, 9.2192679648074E-5, -2.25922396366081E-5,
    0.000136384758638213, -0.000141864997792187, 0.00039226162065277,
    -0.000339039958328376, 0.000188376224063555, -0.000673922132775752,
    -0.000302766382320815, -1.5508495337715E-5, 0.000178379962622114,
    0.00012477314201759, -3.21477234220859E-5, -5.36793231222627E-7,
    -1.51124595652423E-5, -2.75942650915729E-5, 3.00503669319718E-5,
    -0.000285625346261483, 0.00211350481628527, 0.000568482873225118,
    -0.000461126212060126, 4.09855605421938E-5, -0.000197142474434416,
    7.42338389661123E-5, 3.64688095496847E-5, -0.000111369918952771,
    -4.09170631692389E-5, 6.09568359964349E-5, -1.54725114559023E-5,
    -0.000178212771822675, 0.000144521021902355, -0.000118212594195449,
    -3.40893455307623E-5, 6.69716040200858E-5, -0.000161055484339949,
    -9.29060300216524E-5, -1.1793441639194E-5, -0.000138500103454577,
    -8.22655382011021E-5, 1.39263832491066E-5, 2.65883364657087E-5,
    8.79324604629103E-5, 2.66535176873677E-5, -2.26216670347648E-5,
    -2.07093593160861E-5, 6.47438399886312E-6, -6.07485604243762E-5,
    -0.000109646468995589, -2.60131271001756E-5, 9.11271442716157E-5,
    8.1281848446639E-5, 5.48654005796258E-5, 0.000110336574263014,
    0.000141196923930406, -9.10251890827426E-5, 9.04563400857603E-5,
    -0.000106208283900409, -7.3057987210218E-5, -0.000491283470653045,
    0.00224527153065518, 0.000304606445203462, -0.000666190700652125,
    -0.000233556045816482, -0.000101514070062547, -6.84819621126517E-5,
    -0.000303641489809699, -7.98110139356428E-5, 2.62818390677418E-5,
    4.03905891696722E-5, 0.000114758263216579, 1.14223967904315E-5,
    0.000108445123417078, 0.000275377647260453, 1.62089775839415E-5,
    1.26075566218426E-5, -0.000212685444426407, -0.000115769502488818,
    5.78661530300788E-5, -0.000231739372859844, 0.00183619763975945,
    -0.000242765991355045, -2.45057143367774E-6, -2.46894848252109E-5,
    0.000101600721555261, 2.00451672285389E-5, -0.000134550855116562,
    -0.000106400048829932, 9.89753547288184E-5, -0.000212069230102595,
    1.41535908503909E-5, 0.000389661402961303, -0.000532245517949337,
    5.2100913674531E-5, 0.000110944675638728, 8.05263376039559E-5,
    4.96516106366493E-5, -0.0001614181740211, 1.7150603692645E-5,
    0.000183930729009619, 6.42307219319906E-5, -4.85403890349205E-5,
    -0.000103461807435336, 0.000387720207758573, 5.26090738835075E-5,
    -1.03796403949571E-5, 0.000198935010378375, 0.000125798624367615,
    -1.82449296054633E-5, 3.69261211115981E-5, -0.00010055243100584,
    -1.77130996236656E-5, 0.000105575486505621, 0.000229912747059791,
    -0.00014139947349795, 0.0018596001093483, 0.00225653224600177,
    7.31598763551953E-5, 4.88488538311812E-5, 9.34114217936503E-5,
    0.00018079610347821, 0.000166631893366037, -8.82097919481664E-5,
    -0.000157866649535857, -9.47842424325911E-7, 8.04888127161595E-5,
    -0.00024781360280172, 0.00211656119583303, 0.00135834863723899,
    -7.69928245435499E-5, 0.000144971919564722, -0.000173708001172715,
    -5.74610513481574E-5, -0.000123806296190745, 7.36024484229414E-5,
    -0.000103794534759513, -9.26359039495072E-5, -7.13007666784366E-5,
    -0.000200821486000497, -0.00175104298583566, -0.00116969969846181,
    -0.00017623384701319, 6.03029222855941E-5, 2.15579448421089E-5,
    -7.71216162398885E-5, -1.47192665459261E-5, 0.000250653979067571,
    3.47745772232969E-6, 0.000432286011793628, -0.000112098141709971,
    -0.000283192831641065, 0.000183331288786308, 0.000154485437900542,
    0.000190724794439264, 0.000129992602804042, 3.07613104998158E-5,
    -5.79479357144414E-5, -5.58536251785097E-5, 5.40343465516897E-5,
    -3.21451778894633E-5, -0.000220464075302394, -5.67784751545527E-5,
    -9.09128597502381E-5, -0.000113592257874702, 8.21486411445089E-5,
    -0.00039411306576463, 0.00144724196088172, -3.22076419909975E-5,
    0.000173774616193367, -4.12972738229067E-5, -0.000132829205944464,
    7.61591680255883E-5, -8.57771414003449E-5, 0.00136327090795169,
    -9.65833719038963E-5, 0.000104229423459235, 0.000144059067673365,
    -7.94916981134552E-5, 0.000161915192757918, 2.34410810082524E-5,
    -0.000195337959609628, 0.000157518841582731, -0.000120781500215565,
    -0.00113987854596986, 0.00210736075496591, -0.000148772496738032,
    -1.75235273607846E-5, -0.000201887284214114, 8.05115220355425E-5,
    -0.000164292696698201, -5.98967831887334E-5, 0.000235213405047762,
    -0.000182271332021751, -0.000587347583015183, 0.00280026363538977,
    9.45524480810821E-5, 0.00026464187740693, 2.25833821800251E-6,
    -0.000152749563972936, 3.30414931511923E-5, 0.000113779594014048,
    3.35945463885263E-5, 0.000156077711352856, -0.00032953746310281,
    0.00149441410645624, 0.000417256604892638, -0.000104634572186881,
    2.17707173145016E-5, 1.99967748025424E-5, 0.000120315124892753,
    -0.00010906515280456, 0.000289814380438369, 0.000197169647041821,
    -0.000107478955192141, 0.00131220613065154, -0.000226654021217725,
    -0.000953667544288956, 7.49208558159692E-5, 3.50863412274145E-5,
    -0.000174197016213529, 0.000351192037406786, 0.000271713840521435,
    -9.38821295829389E-5, 0.000413441673037327, 0.000683217316204949,
    -0.000963785822837837, -7.07771871688743E-5, -0.00016686805041901,
    -4.7748615778156E-5, -4.0827444250253E-5, -2.55338041194153E-5,
    -0.00029501985172592, 0.0002976001976277, -0.000915834745299217,
    -0.000871374296330663, 6.6198117373428E-5, -0.000308791069963571,
    7.65067813775794E-5, 7.79646255457184E-5, -2.34598794979023E-5,
    0.000153308897442116, -0.000113582817772663, -0.000158493833394834,
    -0.000100180940474921, 2.61685618029979E-5, -0.000235935606355943,
    -0.000228689577850124, 0.000187201667515263, -0.000330798348053006,
    -0.000237437889755027, -0.000144217570365617, -0.000155777127801809,
    -1.19627668381397E-5, -0.000236666415754381, -0.00104794378543355,
    -0.000150002655608424, 9.56903352855338E-5, 5.7527887035811E-5,
    6.7366447917474E-5, -0.000165641985275158, 0.000108859109481904,
    -1.66602988674266E-5, -2.35567372052554E-5, 8.82202599698765E-6,
    -0.000194282135361685, 0.000192652865340413, 0.000161727497951426,
    -0.000195045776130903, -7.28655789009784E-5, -0.000168445523286836,
    -0.000122140074031956, 4.08198788064556E-5, -1.02131022602421E-5,
    1.41319163338427E-5, 0.000112971910145257, -4.97778142210353E-5,
    1.21638235255996E-5, -7.97426819856026E-6, -0.000184848946306788,
    6.71272421981683E-5, 0.000115024160914686, -0.000233110279212054,
    -5.25194262698599E-5, -0.000659373717823057, -0.000183261315678956,
    -6.04613256345022E-6, -7.2161252269819E-5, 1.60637896900606E-6,
    -0.000297317302189945, -0.00177585578392532, -0.00190632372727807,
    -3.71961861296862E-5, -1.20406369473212E-5, 2.64872431360678E-5,
    -0.00018292622724658, -0.000109542507965167, 7.2444655343621E-5,
    -2.15466653133488E-6, 8.93700507132949E-5, -9.02424681272061E-5,
    -0.000102232426398393, 0.0001216436976532, -0.000418585640767557,
    -0.00132994484183482, -0.000863665251998653, -0.000153858278241726,
    -0.000254901563823361, 5.72362133667335E-5, 7.48512376310948E-5,
    -0.000123344878972768, 9.1928388660846E-5, 3.4897134007845E-5,
    -8.85143901379312E-5, -9.75225151353729E-5, -0.000160696079881299,
    -0.000120149838233437, -0.000248002387200356, -0.00162407627756439,
    0.00107100724663451, -0.000125409510851037, 8.81162679413792E-5,
    -4.53696861049128E-5, -0.000106413966036138, 7.12840251648E-6,
    -0.000112991439885765, -0.000175152095440573, -8.91231225644013E-5,
    4.33715148515628E-5, 5.16924811762388E-5, -0.000220242031255919,
    0.000237515303770877, -0.000472619665200404, 1.46954171547805E-5,
    -8.56473433728604E-5, -0.000173009139583444, 6.60891610331126E-5,
    -0.000149390164826322, -9.99739352622395E-5, -7.76141719137235E-6,
    -0.000300147525526068, -0.00142101972614065, -0.000142485322469937,
    0.000158636057070292, -8.36025692294827E-5, 7.1746976459127E-5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
    2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
    2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 3.0, 2.0,
    1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0,
    2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0, 2.0, 3.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 3.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 65.0, 65.0, 67.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 62.0, 62.0, 62.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 75.0, 75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 86.0, 87.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 68.0, 68.0, 225.0, 159.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 81.0, 82.0, 82.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 89.0, 89.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 93.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 88.0, 105.0,
    103.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 166.0, 173.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 75.0, 75.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 79.0, 79.0, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    91.0, 91.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 84.0, 86.0, 86.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 91.0, 93.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 71.0, 93.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 88.0, 93.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 88.0,
    89.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 74.0, 74.0, 75.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 97.0, 97.0, 96.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 82.0, 91.0, 93.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 73.0, 73.0, 74.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 75.0, 75.0, 77.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 73.0, 73.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    71.0, 71.0, 71.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 81.0, 80.0, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    76.0, 76.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 78.0, 78.0, 78.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 75.0, 75.0, 76.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 71.0, 72.0, 72.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 86.0, 91.0, 91.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 56.0, 59.0, 59.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 57.0, 61.0, 61.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 72.0, 72.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.21387221225495E-8, -4.34147993818704E-7,
    1.82536849903729E-7, 1.32342420837361E-8, -8.67910904622389E-8,
    -6.287912598747E-8, 1.33644147926427E-7, 6.17104980887688E-8,
    -1.67817947671719E-7, -2.00669886719305E-7, -2.93542397805959E-7,
    9.4587248481987E-8, 3.60003371680568E-8, 1.28494650955011E-7,
    -2.25170125206385E-7, 2.77535397047812E-7, -3.454019246998E-7,
    -5.5517914229208E-8, -2.80206133464837E-7, -1.44147904968319E-6,
    -6.03784154321408E-7, 6.1360757963723E-7, 5.24659684920928E-8,
    7.46844810415829E-7, -1.50751304621711E-6, 7.5387730427209E-7,
    -1.1192026099793E-6, -2.59786034483552E-7, 2.02282866434643E-7,
    1.54610875319253E-7, -3.0231757605572E-8, -6.83963512091439E-7,
    3.38027014495573E-7, -8.19367208682695E-8, -3.73787120235786E-7,
    3.47577067881697E-7, -4.13300860940485E-9, 4.34510223781827E-7,
    5.54654912321001E-7, 1.13990066080849E-7, 3.65349728742261E-9,
    -3.43552596365223E-7, -2.31569567587481E-7, -3.08815271979719E-7,
    -3.73729773591444E-7, -7.09079284276833E-7, -4.29800246053968E-7,
    4.15762759040251E-7, -4.92871211091781E-8, -1.78280736455164E-7,
    1.45142979296428E-7, -1.45040429777783E-7, 9.01040990517429E-8,
    -2.14457685703128E-7, -4.12799431442163E-7, -2.15383423657993E-7,
    -1.68176907593573E-7, -6.37202495185698E-8, 1.27675082200064E-7,
    1.02513874888799E-7, -1.22399613088269E-7, -8.49582772211381E-8,
    3.05408542167219E-7, 2.00280778397579E-8, 9.35409415354554E-8,
    -2.53479281232369E-7, -1.03762163254356E-7, -4.98728500295364E-7,
    1.21015813928551E-7, 1.1248752834222E-6, -2.57532603977376E-7,
    2.89409070217071E-7, 2.07329294326247E-8, 2.04599287245018E-7,
    -1.5842699229976E-9, -2.4245776206802E-7, -1.51586701085658E-7,
    -4.25537820126081E-7, 5.17877676767292E-8, 1.67312978332866E-7,
    4.312131276383E-7, -6.09823205350397E-7, -1.26019748324702E-6,
    1.14918848124383E-6, -5.4507879799898E-7, 3.28724165681415E-8,
    -1.69769468793486E-7, -1.31104443792093E-7, 1.73004648089608E-7,
    -5.59989025443241E-8, 2.13064769183182E-7, 5.5559549721084E-8,
    -6.51273404608334E-8, 2.39956849613628E-7, -1.65767786701502E-6,
    -4.27484020144242E-6, -4.59312030303499E-7, -1.80839364798371E-7,
    -4.58620340535408E-8, 3.88314839933596E-7, 2.36665149749838E-8,
    1.36587999055394E-7, -1.89010075416978E-7, 1.35915410443041E-7,
    2.14822132040819E-7, -3.85552869822426E-7, 2.16146119272309E-7,
    -4.13637751569626E-7, -8.3051346281607E-7, 1.14581746866472E-6,
    8.51397652693388E-8, -1.79677052730138E-7, -2.26118309176084E-7,
    -6.30235682695462E-8, 1.24417383949267E-7, -4.6502038174354E-8,
    -3.84221605354552E-7, -3.96863075911149E-7, -1.72280653207355E-7,
    -7.45674449800134E-7, -1.59140476623444E-6, 5.80747334348844E-7,
    2.35619899354363E-7, 5.06761699822533E-7, -8.50922024234922E-8,
    -2.21349783586895E-7, -1.00712876386901E-7, -2.55904472673006E-7,
    2.50050981996248E-7, 4.80816276953756E-8, -2.5516434927496E-7,
    1.21565272944946E-7, -1.77713312950963E-7, 3.2047430228777E-7,
    -9.38413874013483E-7, -2.46535212942031E-6, 4.72416718782355E-7,
    -2.31319549031047E-7, 1.42700652693703E-7, 2.92199802456638E-7,
    2.40489302571021E-7, -2.24554668321188E-8, 1.15002509773962E-7,
    7.8440146009424E-8, -1.11875041772199E-7, -4.95695359173662E-7,
    1.51225147748733E-6, -4.24598141006425E-7, -4.61584193094808E-7,
    1.42432841922905E-7, 1.2278804249658E-7, -3.32004793129954E-7,
    2.88561224196928E-8, -1.1847720001321E-7, 2.35134454458572E-7,
    4.26636519666656E-7, 2.26068450725345E-7, -3.76080546092257E-7,
    7.53379515590179E-7, 1.92105285157686E-7, -2.72505747093668E-7,
    7.29515394687025E-8, -3.07202880681882E-7, -4.63198331585441E-7,
    -6.05472019208017E-8, 2.2196469908005E-7, 1.86260694958549E-7,
    -1.06042151601055E-7, -4.84438265134206E-7, 6.49565964851388E-7,
    -1.20611019045248E-7, -9.50114187911376E-7, -1.0014217435E-7,
    2.07550608655405E-7, 2.4685402736002E-7, -4.81230144383421E-7,
    3.88054509568579E-7, -2.21755399143209E-7, 3.16899895816293E-7,
    -6.81171259526476E-8, 1.0605127089907E-7, -1.38039803211878E-7,
    -4.35549847371745E-8, -5.45819194132947E-7, -2.7495717625224E-7,
    -2.0367425771068E-8, -9.12277437916706E-7, 8.97648264677918E-8,
    -6.00563910807296E-8, 2.07605178397181E-8, -1.34913532152634E-7,
    2.40097054172149E-7, -1.71458891171781E-7, -1.01963906274961E-6,
    1.6171408828731E-7, -1.35556911803458E-8, 2.70782107492764E-7,
    -1.59809764479738E-7, 1.32106338289337E-7, -6.64359965929797E-8,
    7.84538004036237E-8, 1.00353186362396E-7, -5.34165334871469E-7,
    -6.13002939260199E-7, -4.08063940192893E-7, -1.0590087128706E-6,
    -1.55421879926806E-7, -4.15228987163271E-7, 1.3146569558674E-7,
    3.66740532421404E-8, -9.92003437296354E-8, -9.0781859298512E-8,
    2.37469722115559E-7, 2.73019684098682E-8, 1.02334964893977E-8,
    4.48648565692982E-8, 5.32816366233888E-8, 5.95983039796086E-7,
    -2.54686974318462E-6, 8.05744773248914E-6, -5.97049138040323E-7,
    -1.28194717786657E-7, 5.68426234381676E-8, 3.81021065334576E-7,
    1.13011521184933E-7, 1.64546127445175E-9, -2.12857693969063E-7,
    -2.47348173123397E-7, 3.87553149591797E-7, -6.5089032700213E-7,
    -3.24550368789062E-6, 3.65761432019731E-6, -3.5825395938042E-7,
    -2.26653078681352E-7, -6.50948100009558E-8, 5.90524325418533E-7,
    2.86640084374724E-7, 5.75554311086342E-7, 4.97728977117464E-7,
    -4.16708292390348E-8, -3.21182531246976E-7, 2.21421770198098E-6,
    -3.56888721321938E-7, 1.76524713499491E-7, -1.27231077328156E-7,
    -2.21536196397904E-8, 1.69994123677722E-7, -3.77546023510618E-8,
    -1.23010128742808E-7, -1.53565439434989E-8, 9.98281702140296E-8,
    -2.21961444013728E-7, -2.68353141612493E-7, -5.20386465207446E-7,
    -2.44233356496474E-6, 3.65544305998612E-7, -1.97002896249716E-7,
    7.89824335136602E-8, 1.77923975760872E-7, -4.37525034855908E-7,
    -3.06911693219353E-7, 6.6657292471533E-8, -6.37701334486052E-7,
    -1.9568423452181E-6, 3.85092377396597E-7, 2.24826798288039E-7,
    -2.27016164213588E-7, -2.63241798567719E-7, 3.07473125848956E-8,
    -1.03871226567447E-7, 2.93949752322307E-7, -1.73650176584285E-7,
    -1.50964491235963E-7, 4.85366136114469E-6, -1.14509252442502E-7,
    -1.5937827665034E-7, 1.24187709633056E-7, 4.77795606004144E-8,
    -3.52556727708943E-8, 2.477770497437E-7, 6.93279774549559E-8,
    -4.41103304618483E-7, 1.72576984840495E-6, -1.01740974268839E-6,
    -1.30042175239087E-6, -3.3340792162475E-7, -4.35710460957934E-7,
    -1.53526679816506E-9, 1.09329323906299E-8, 2.01885060521446E-7,
    1.44293735775098E-7, -4.79874634144717E-7, 3.94487994317977E-7,
    -7.73269681227645E-7, 9.34798855344636E-8, -2.35069059532046E-7,
    7.76637527591855E-9, 3.96450089800592E-7, -1.9628943144157E-7,
    3.39364862983266E-8, 6.60026686028629E-7, 3.05437328437232E-8,
    2.03624853984596E-8, 3.10344795518969E-7, 1.3004566948991E-6,
    -1.01010094232944E-6, -1.37354081955839E-6, 6.90543257452361E-8,
    -2.22403998116776E-7, -1.40768870422153E-7, -2.52509316183282E-7,
    -1.47645071901609E-7, -4.17435152267063E-7, -9.69088912351049E-8,
    6.02613498810367E-7, -4.66706231760182E-7, -1.40697524960799E-6,
    2.94053988465562E-6, -1.86524538926467E-7, 6.59454138776638E-8,
    -3.35953715064532E-7, -1.67566578007276E-7, -2.73978021756104E-7,
    1.02155775276153E-7, -3.43841570271836E-8, 6.34732739745947E-8,
    2.21784905716785E-7, -4.01411653276297E-6, 1.45808063682311E-5,
    -6.22547637641174E-7, -2.08619074568861E-7, -2.60525112931793E-7,
    9.91999994461817E-8, -3.41579383877827E-7, -5.40461434845203E-8,
    4.53998180012116E-8, 1.26008483100116E-7, -7.30465211892496E-7,
    3.68593181691068E-6, -1.59685495346548E-6, -3.32592069054025E-6,
    -2.95183939575432E-7, 3.37780462420602E-8, -1.10946216541855E-7,
    2.48801986443879E-8, 3.50010036407317E-8, 3.57014033553113E-7,
    1.88893718203271E-7, 8.57122265843088E-7, -8.56339646365797E-7,
    -1.72418263790451E-6, 3.66822667112338E-6, -2.01607209003086E-6,
    -5.82104060357884E-7, 2.9776890521792E-7, 1.17430632038249E-7,
    1.35120163850123E-7, -3.4736879412187E-7, -1.75784074750505E-7,
    -1.13057223983298E-7, 1.02033132935168E-6, -1.93818687654298E-7,
    1.52937281807613E-6, 9.87796184570069E-6, -2.37794758724378E-6,
    -1.02744019695713E-6, -2.89868784272856E-7, 5.118137578632E-7,
    1.76441011872612E-7, -1.49747440734065E-8, 1.72954304498351E-7,
    -6.82228603887448E-7, 6.73363382362979E-7, -9.08284482820289E-7,
    1.88671260095249E-6, -1.54475846635435E-6, -2.0835885523508E-6,
    -2.26649602490173E-7, -2.60452100415509E-7, 3.97474152635774E-8,
    -3.90793697983679E-7, -1.08456556218963E-7, -2.33269512568103E-7,
    2.57392282948277E-7, 2.72845222116628E-7, -4.76852264104231E-7,
    4.742674478868E-6, -3.32966229901254E-7, -1.35413992355682E-6,
    3.71965803201859E-8, -3.10822387577835E-7, 1.1658272297991E-7,
    1.30500485848977E-7, 9.1701178293048E-8, 3.64899518089503E-7,
    4.10433267180787E-8, 2.07276147749452E-7, -7.03587100977338E-7,
    -4.05613695594246E-8, -3.46179770136107E-7, -1.70502170004645E-8,
    2.66511895109686E-7, 9.88617583370393E-8, -4.92869698613758E-7,
    -3.09132527905215E-7, -7.55499238021323E-8, -3.40992705244802E-7,
    -4.48241692236504E-8, 1.42767278974067E-8, 5.83672349176384E-8,
    7.60845147452632E-8, -1.85780080244746E-7, 7.78585135610372E-8,
    1.21300168297654E-7, -1.88916789957541E-7, -3.23813932475989E-7,
    -1.20665485972965E-7, 2.90674102360731E-7, -6.66570177276908E-8,
    1.22027788426204E-7, -4.20312990588783E-8, 3.98748729597136E-7,
    -4.5893436749014E-7, 1.33843043268908E-6, -2.069682741682E-7,
    -4.6238297292115E-7, -1.65695571928222E-7, 3.70179411852453E-6,
    -6.58320628392751E-7, -2.53657617337266E-7, 5.10762872885862E-8,
    1.7140255238847E-7, 7.09837950926684E-8, -1.6599920470024E-7,
    -3.4773227433878E-8, 4.06220354620026E-8, 2.29181714909359E-7,
    -6.90481724342654E-8, 1.76415561100688E-7, -1.05942183500573E-7,
    -2.09424917464419E-7, -1.19058028881863E-9, 1.01136082405995E-7,
    5.72545944465657E-8, -9.28381502110859E-8, 9.4014659095447E-8,
    3.39191975025297E-7, 7.74524951539058E-8, -3.69747491220668E-7,
    -4.95973712285286E-7, -2.19500850590034E-7, 7.05101413043694E-8,
    1.2301259827391E-7, 1.18351552223698E-7, -2.38870418866121E-7,
    -7.98291018350435E-8, 7.98642728195048E-9, 6.361139205349E-8,
    7.22411610250224E-8, -1.60503841648241E-7, 4.9501772938554E-7,
    5.5281893235679E-7, 5.24896075175304E-8, 1.32530470289077E-8,
    1.92831755629139E-7, 2.66356430168674E-7, -2.65661028879006E-7,
    7.49367294758279E-8, -1.50129659030212E-7, -5.82212584677291E-8,
    3.19448813945543E-7, -1.53667723057383E-7, -1.54625601812867E-7,
    2.80166483528883E-7, -4.20632347924949E-8, 3.0522506812769E-7,
    7.8526388798524E-8, -8.93295060311113E-8, -2.88127450023908E-7,
    1.64967027877003E-7, 1.17928341398848E-7, 1.42180647623447E-7,
    -2.48804755631885E-7, -7.69073799334133E-8, 5.79417432574061E-7,
    4.52794376032709E-8, -1.30414477932812E-7, 1.31887923043464E-7,
    2.82649551980925E-7, -2.09067300840612E-7, 9.57906246994279E-8,
    2.05356902278808E-9, -2.62592674547183E-7, 6.67150556115416E-8,
    3.01414702502433E-8, 3.45795930414852E-7, 1.14065402282684E-7,
    -1.54847174450485E-7, -1.38802845338772E-7, -1.77292863480284E-7,
    1.45505339858222E-7, 6.00082487864117E-8, -2.97703848161412E-7,
    -8.87230716229942E-8, 1.72104487593203E-7, -4.91452829510269E-8,
    -3.53398302625921E-7, -9.22203694590175E-7, -6.4757821127882E-8,
    2.67971800307498E-9, 3.3012337996502E-8, 5.72839275799989E-9,
    -2.63298061646363E-7, 1.54487781346343E-7, -1.80621870557968E-7,
    5.07666808216369E-7, 7.14754620554708E-8, 9.96262694816564E-8,
    -1.75210312059961E-7, 2.2106811671848E-7, 3.30254641856145E-7,
    3.34529055767098E-7, 7.11845159794045E-8, -1.89734954727591E-7,
    -4.73091055655546E-8, 1.0115490515198E-7, -1.86422286269249E-7,
    -8.30763824446997E-7, -3.70777511518896E-7, -1.45241448466763E-7,
    -3.57371930149554E-7, -8.64619281286212E-9, -2.93837198524022E-7,
    1.49915751770992E-6, 3.0513375199555E-7, 1.92710716647462E-7,
    -1.31670000960177E-7, -1.09988392760662E-7, 2.19770261808478E-7,
    1.24838020752354E-7, 4.62041645797303E-7, -2.18268482446889E-7,
    2.4543797067881E-7, 2.31148134989515E-7, -2.25886522796555E-7,
    -5.51818674398657E-8, 2.33236544868659E-7, -6.4359918418424E-8,
    -2.65495120537919E-7, -9.01938930217087E-8, -1.5845109479757E-6,
    1.82887082696289E-6, 8.95605583021985E-8, -3.12381511804967E-7,
    -5.23761833070722E-7, -1.11246109517977E-7, -7.54182613341989E-8,
    -1.13202423379598E-7, 2.84061193715022E-7, -1.21760622494617E-7,
    -7.22313853925491E-7, 4.07785181851684E-7, 2.47476260923004E-7,
    7.20796042340423E-7, 3.23574034888103E-8, 2.58333711839562E-7,
    2.14552306047613E-7, 3.76968279999818E-8, 3.21482543775645E-7,
    2.99542155810474E-7, 2.55728899557814E-7, 6.06555129707328E-7,
    -1.33750840225235E-7, -9.12705065247556E-8, -3.03180498204053E-7,
    -4.70771963059736E-8, 8.10998515877123E-8, -3.91935987997029E-7,
    4.15371010964197E-7, -3.94948231343251E-9, -4.30994511294868E-8,
    1.14170827226919E-6, -6.36911699408703E-7, -1.84964775422925E-6,
    6.10492595259479E-8, 2.67614303931156E-7, -5.32773567391948E-7,
    1.58958990497916E-7, -1.8793250051729E-7, -1.86998199096988E-7,
    1.24810800599679E-6, -4.6436735288878E-7, -1.31566072199545E-6,
    -7.19491868783991E-7, -1.27981984521813E-7, -9.93738883299089E-8,
    -5.7258312713062E-7, 6.11107924659968E-8, -6.14696789033879E-8,
    3.53793020895042E-7, -9.4450879780214E-7, -1.01505333821183E-6,
    7.42719516846614E-8, -1.004198027383E-7, 1.003767182288E-7,
    4.84893777569306E-7, -2.36526078962919E-7, 1.57544312444847E-7,
    -7.22192791234162E-7, -3.93983260488023E-8, -5.38781706241147E-7,
    5.23577651680159E-8, -1.56659655232594E-7, -1.67101809817036E-7,
    4.21962014529032E-7, -1.08375402715852E-8, -2.99022459453455E-8,
    -3.33525298248466E-7, -1.65297192790947E-7, 4.46601396833959E-7,
    -1.25159770898812E-7, -1.92009938477642E-6, 1.77399524725958E-6,
    -3.74054983855853E-7, 2.73216770396306E-7, -6.44378117523976E-8,
    -4.11353973372133E-8, 2.34010829436904E-8, -4.60867482596274E-7,
    -1.209011992261E-7, 2.29442956379382E-7, -6.87477863064042E-7,
    -5.15432726953189E-8, -4.46280570609723E-7, -2.5490719817955E-7,
    7.38090560795402E-7, 5.33737318650733E-7, -6.19639887511499E-7,
    7.95424228837466E-8, 5.85840641311656E-7, 1.25227484305807E-7,
    5.50221968932148E-7, -1.40332660668747E-7, -1.72320879192877E-7,
    1.91886206384262E-7, -1.55009930920035E-6, 3.47198828730379E-7,
    -3.37465534137556E-7, -1.0335050020754E-6, 1.02123100930929E-7,
    5.05533920117405E-7, -4.19630445421544E-7, -3.36779164181297E-7,
    -1.2401454201482E-7, -1.56670835716133E-7, -1.01790981958289E-7,
    -8.02424947676753E-7, -2.36508285500548E-6, -1.69729400370019E-7,
    -1.7673313385326E-7, 2.36286824411956E-7, -5.21403151664751E-7,
    -3.86968879425773E-7, 1.44168736258205E-7, -1.94001270310786E-7,
    1.92705647301597E-7, -9.31087572585292E-8, -1.58452291202876E-7,
    3.11205716628185E-7, 3.81538498848545E-8, -9.24245028603215E-7,
    2.28829867959898E-6, -2.26654179296542E-7, -2.14505782320171E-7,
    -1.80047596439112E-7, -2.57106541211492E-7, 2.76125025663579E-9,
    1.32501843390742E-7, -1.45197324791844E-7, -4.37392699098805E-7,
    -1.84336567383689E-7, -8.01491873748554E-8, 2.05557396801388E-7,
    -3.51395938493151E-7, -9.27732478727331E-7, 3.43862247768876E-7,
    -3.31360431285629E-7, 2.80556734548765E-7, 1.71923410501001E-7,
    -2.97107292607669E-7, 4.42233049737186E-8, -3.31650411062493E-7,
    -2.28408114261145E-7, 2.17250209719656E-7, -1.39700585873548E-7,
    -9.54933769263372E-8, -1.77019419734345E-7, 9.68224174048893E-8,
    -3.76589088037005E-7, -1.61389574713946E-6, -2.23690854601178E-8,
    -5.58151153109008E-7, -6.06840356485864E-8, -3.80946040480979E-7,
    -9.03343118798607E-8, -1.97307788379552E-8, -3.1390092698135E-7,
    -6.35120573063135E-7, 5.73113324617419E-7, 1.25176154771964E-7,
    -4.00736132611301E-7, 4.88791627140089E-7, 1.69917464022095E-6,
    2.20357265280238E-6, 1.83038172682403E-6, 2.29084302700059E-6,
    2.33891399441625E-6, 2.19399889748294E-6, 2.12796509729049E-6,
    1.2102473399736E-6, 1.52649025088967E-6, 1.89337397021136E-6,
    1.7591324247284E-6, 2.20323616858915E-6, 2.25646254297948E-6,
    2.1993191798843E-6, 1.99399104667711E-6, 3.05417857512349E-6,
    3.49648662025727E-6, 2.54227200277779E-6, 2.66832145668082E-6,
    6.4498137254616E-6, 4.90146851561292E-6, 4.53350258042618E-6,
    3.2400756075939E-6, 4.20799378255727E-6, 6.79726723472891E-6,
    4.28490942042981E-6, 5.0181390138552E-6, 5.40905184114923E-6,
    7.33896276944054E-6, 7.34127287582937E-6, 7.2947073188812E-6,
    3.44076728165768E-6, 2.7625187268969E-6, 2.64905247103082E-6,
    2.71125276747894E-6, 2.5953809757241E-6, 1.74337302059205E-6,
    2.80049660435682E-6, 2.94932064897425E-6, 1.72413581840104E-6,
    1.21055284213686E-5, 1.46195523471776E-5, 1.48643830627129E-5,
    1.49194150683693E-5, 3.32406382284551E-6, 4.34054650880219E-6,
    3.26854454671677E-6, 3.41091178965359E-6, 2.88947355345667E-6,
    2.74729368434204E-6, 2.88497067596578E-6, 1.09194480916923E-5,
    1.2210069492663E-5, 1.23058161417373E-5, 1.26547078297678E-5,
    2.415150621842E-6, 2.57336872693777E-6, 2.30562629392596E-6,
    2.34845371206356E-6, 2.58298355383248E-6, 2.59000867324391E-6,
    2.24176173811383E-6, 3.41578923262424E-6, 3.12874680577012E-6,
    5.71870023117546E-6, 7.98143764187031E-6, 9.46327698814463E-6,
    9.73391914552469E-6, 7.75459515305937E-6, 6.26734607894293E-6,
    2.8366269315593E-6, 2.71523416196276E-6, 1.93220363230183E-6,
    1.66274541844814E-6, 1.5922649765145E-6, 1.90086067540581E-6,
    2.08853501548801E-6, 2.59417090005043E-6, 4.26078273568225E-6,
    1.38103422575341E-5, 1.74227688195778E-5, 1.73470788665962E-5,
    1.47031762858632E-5, 7.22964372821043E-6, 3.5953691857071E-6,
    2.2976490616978E-6, 2.04047713502774E-6, 1.79265026484253E-6,
    2.11539455219557E-6, 2.08976579778499E-6, 1.38528614386035E-5,
    2.50181880209467E-5, 3.54296698956108E-5, 3.58077154907821E-5,
    3.01817009740791E-5, 1.65462578564086E-5, 4.25948325474797E-6,
    3.37194794830021E-6, 3.00763749791901E-6, 2.55293447130965E-6,
    1.79441060989495E-6, 1.88232994584582E-6, 2.0713503673261E-6,
    2.51009949687178E-6, 3.35727631658371E-6, 1.04497891601004E-5,
    1.32417930967041E-5, 1.33653334397669E-5, 7.9224281534023E-6,
    6.35244360291047E-6, 2.18253503101316E-6, 1.7637275340335E-6,
    1.86087146469749E-6, 1.62779228064574E-6, 1.96743910800404E-6,
    2.74115715293099E-6, 4.74774496796328E-6, 8.28249174805473E-6,
    8.63715254400202E-6, 8.91908818783249E-6, 8.22127177159759E-6,
    6.15024008134294E-6, 5.7250435246084E-6, 4.54399153556724E-6,
    1.87951695195063E-6, 1.8014543385134E-6, 1.82859271803069E-6,
    2.40045133409318E-6, 2.53999298595791E-6, 2.33304826961147E-6,
    2.82835173574863E-6, 3.4639088714679E-5, 4.54290695864468E-5,
    4.69402772138048E-5, 4.61014861879997E-5, 9.24393427342232E-6,
    3.35780902449746E-6, 2.86369864648282E-6, 2.79780909206032E-6,
    2.42423563694989E-6, 2.47492864856783E-6, 1.57968230406577E-6,
    3.3944607181953E-6, 4.18070336585956E-5, 4.52746879783779E-5,
    4.5440623282344E-5, 3.96629847804472E-5, 4.13182197527691E-6,
    3.81188384884144E-6, 2.67683966364038E-6, 3.45338070307124E-6,
    3.29673243981258E-6, 2.33038971198353E-6, 1.96721769363732E-6,
    2.53194214023291E-6, 9.47504708859868E-6, 1.01649207301537E-5,
    1.01906076604288E-5, 6.64503424554834E-6, 3.71288806137322E-6,
    2.8268688350947E-6, 2.6217004114962E-6, 3.6368443962353E-6,
    3.81487381452619E-6, 2.81032441905988E-6, 3.23438686842423E-6,
    1.06293213126084E-5, 1.72527596183714E-5, 2.0725169203219E-5,
    2.10838992876109E-5, 1.52187705011147E-5, 5.26671290772697E-6,
    3.23652771598086E-6, 2.93996292576001E-6, 3.09385200609859E-6,
    3.41889827090511E-6, 3.5834837522879E-6, 2.75463679495666E-6,
    2.25919691218505E-6, 7.20536089924793E-6, 9.42548192481193E-6,
    1.28072949491695E-5, 1.30270810563481E-5, 8.63307679339572E-6,
    4.64434576600273E-6, 3.08789588609844E-6, 4.39547695884575E-6,
    2.34965566832504E-6, 2.33212859734073E-6, 3.64128920427264E-6,
    7.96298391081464E-6, 1.09282390881211E-5, 1.15085132997393E-5,
    8.02353676173625E-6, 6.69108474734081E-6, 1.77381148016694E-6,
    2.47309754370427E-6, 2.2716892506198E-6, 2.10646877669932E-6,
    6.10747462789183E-6, 1.04561233827552E-5, 1.61090337382805E-5,
    1.86858793180402E-5, 1.40986296884473E-5, 1.38760886650207E-5,
    5.19996722630617E-6, 1.64311051159894E-6, 3.10438942906373E-6,
    1.86475739955495E-6, 2.08419573389662E-6, 2.51756992204708E-6,
    1.93631969626569E-6, 3.95699483532312E-6, 8.79767200201716E-6,
    1.03640598043787E-5, 1.29072500274219E-5, 5.94779872923924E-5,
    6.09590863003083E-5, 6.10417307860555E-5, 5.02587927533147E-5,
    3.27625174956534E-6, 2.44887833523735E-6, 2.39351933662443E-6,
    3.65414539131668E-6, 1.47920299631106E-6, 1.42032614085561E-6,
    6.01762493071713E-6, 1.83336656538042E-5, 3.75014440921699E-5,
    3.9008222535163E-5, 3.58846504559342E-5, 2.57283315436722E-5,
    2.6407022587959E-6, 2.45742449186519E-6, 1.88742695767816E-6,
    3.14727069846264E-6, 1.99054397405851E-6, 3.29835419975568E-6,
    2.309420268169E-5, 2.46446311104181E-5, 2.47563399713788E-5,
    2.27105537678598E-5, 2.38605338499359E-6, 1.81243699213032E-6,
    1.5663222213785E-6, 1.89495665937391E-6, 1.70961665010115E-6,
    1.59605870569096E-6, 2.47623187132655E-6, 2.33332098104334E-6,
    7.38407456064418E-6, 1.4840046881136E-5, 1.53871958095778E-5,
    1.53810354351071E-5, 9.79639881983203E-6, 3.02796905147582E-6,
    2.14470474626664E-6, 1.98560737035027E-6, 2.00337368517486E-6,
    1.28382455442784E-5, 1.79645771367026E-5, 1.87253893237741E-5,
    1.858134581108E-5, 7.82284162207559E-6, 3.65880966934651E-6,
    2.69800839597032E-6, 2.88838107240395E-6, 2.75563737376148E-6,
    2.835210512859E-6, 4.2255495861706E-6, 1.50670568051526E-5,
    1.58137839170983E-5, 1.58563297679928E-5, 2.44746192706363E-5,
    1.65913371103053E-6, 1.77636049767866E-6, 1.94812613646415E-6,
    1.88953955655969E-6, 2.02306453950496E-6, 1.19722242756406E-5,
    2.46451292326117E-5, 2.75524979527313E-5, 2.81495959092809E-5,
    2.46430124619075E-5, 5.29640757158139E-6, 2.66224029908469E-6,
    2.89535720985467E-6, 2.06476871450037E-6, 1.98394991705246E-6,
    1.01504065459146E-5, 1.3168202066903E-5, 2.85063757439549E-5,
    2.90190711111449E-5, 2.52078971147193E-5, 9.96558991294234E-6,
    3.25557253671738E-6, 2.84297489501194E-6, 2.5529489255382E-6,
    2.38114354068005E-6, 3.48997564182048E-6, 4.56350614353128E-6,
    1.01670179598023E-5, 1.55806357913136E-5, 1.72628445949319E-5,
    1.71604510962204E-5, 1.37491963391874E-5, 5.36753769969395E-6,
    3.07275785074561E-6, 2.94587560787288E-6, 2.89685742427236E-6,
    2.24155151060863E-6, 2.54787666469963E-6, 1.09004435521947E-5,
    1.55675375883903E-5, 3.19348082195718E-5, 3.41935931379485E-5,
    3.06268351758597E-5, 2.48022850236006E-5, 3.35908995696608E-6,
    2.26224844326647E-6, 2.4210977653519E-6, 1.88408749558188E-6,
    1.96178706448992E-6, 6.78012917743034E-6, 1.6703714651835E-5,
    8.46909797931996E-5, 8.89184350651558E-5, 8.72143005527649E-5,
    8.11300391070702E-5, 3.81392661655647E-6, 2.40427976888495E-6,
    1.90125004075791E-6, 1.79161361209457E-6, 1.71338730892153E-6,
    3.78926022019861E-6, 1.97170837567881E-5, 4.76334550368638E-5,
    5.86116374180261E-5, 5.97834021852006E-5, 5.21973449668761E-5,
    1.30173828930714E-5, 2.8724363808694E-6, 1.6895336975259E-6,
    1.37658357088792E-6, 1.70409694345794E-6, 1.91826730246299E-6,
    5.01021427510081E-6, 1.97244460541087E-5, 8.17286867216795E-5,
    9.89402604315148E-5, 9.66178936069005E-5, 8.97916099562152E-5,
    7.73001179502812E-6, 4.11748742087268E-6, 2.89013197248178E-6,
    2.91527279551229E-6, 2.59292246216655E-6, 2.64174218491239E-6,
    4.4902804566734E-6, 2.24661180722198E-5, 8.81709714664856E-5,
    0.000101661198231997, 0.000103112767474352, 8.61750940905506E-5,
    9.89220004815011E-6, 7.12604121635449E-6, 2.7735514937439E-6,
    4.34333421415385E-6, 3.68425231819438E-6, 3.09631052554358E-6,
    1.06340437942183E-5, 2.41616382408746E-5, 4.46410390533726E-5,
    5.2712680003312E-5, 4.98668132755613E-5, 4.65486479885535E-5,
    8.07804229769316E-6, 2.74199404612672E-6, 2.41979990523643E-6,
    1.81265814611458E-6, 2.57523070370417E-6, 1.87485287784845E-6,
    3.85075426473067E-6, 1.68590637169689E-5, 3.35130905591592E-5,
    3.98813804519928E-5, 4.49695408204927E-5, 3.28997515076946E-5,
    5.41571928929532E-6, 2.04991681057722E-6, 2.23080631081952E-6,
    1.47299118480746E-6, 1.89191117130486E-6, 1.05346829958611E-6,
    5.65507836356329E-6, 1.43626356267299E-5, 2.6703132637041E-6,
    3.38030210565072E-6, 2.17035282782434E-6, 1.88039644267192E-6,
    1.53954143891083E-6, 2.29460556016449E-6, 2.2043906139281E-6,
    2.75626581511297E-6, 2.25682223607411E-6, 2.06525855418931E-6,
    2.76651092745725E-6, 2.8329583079877E-6, 2.78720959577209E-6,
    2.85263416631476E-6, 2.78489778865095E-6, 2.748778356149E-6,
    1.98812866912214E-6, 2.32935334366113E-6, 2.45428025910659E-6,
    2.54741235110322E-6, 2.37159104612046E-6, 2.11524353631226E-6,
    1.91907226394418E-6, 1.95606096361336E-6, 8.79804007184064E-6,
    7.94648908643481E-6, 7.72189837959884E-6, 6.03665392897761E-6,
    5.77242673455945E-6, 6.1183119582187E-6, 6.04451950381288E-6,
    1.73360471255544E-5, 2.69746204637626E-6, 2.36135723923603E-6,
    1.47308796123788E-6, 2.04544835761576E-6, 1.28481683774533E-6,
    1.49019640968903E-6, 1.87230926143453E-6, 1.84079817418567E-6,
    2.64381753047857E-6, 2.38646062566769E-6, 2.16023098910329E-6,
    2.00053953223986E-6, 2.37739305730902E-6, 2.35246287383827E-6,
    2.51567896008896E-6, 2.620024230913E-6, 2.86444589540218E-6,
    3.20780592718113E-6, 3.34973502090201E-6, 3.20407763574375E-6,
    3.4261247084153E-6, 3.24800518332321E-6, 2.66423511605264E-6,
    2.54083469505129E-6, 2.03743021672869E-6, 2.05987981016763E-6,
    1.80108643070455E-6, 3.65347457924338E-6, 3.34938014192492E-6,
    3.57503623072545E-6, 3.98632115257973E-6, 2.78674614220201E-6,
    3.93200816655508E-6, 3.30785499840875E-6, 2.6029749553654E-6,
    2.22948463392361E-6, 2.24086807568748E-6, 3.28251359927246E-6,
    2.91132171153715E-6, 2.98937395157915E-6, 3.40834406758377E-6,
    3.54969152191917E-6, 3.40204670072107E-6, 2.41012216426309E-6,
    1.62373659558586E-6, 1.71762502622387E-6, 1.4291098947411E-6,
    2.11041008749009E-6, 1.50875261525518E-6, 1.29690462972252E-6,
    2.19208377524785E-6, 4.19217249350823E-6, 5.35633566710596E-6,
    5.55587185563649E-6, 4.36188894174623E-6, 3.66158360294308E-6,
    4.04378840370112E-6, 1.80108498514442E-6, 2.0106264910633E-6,
    2.39928574474579E-6, 2.26683421075126E-6, 1.92258533567581E-6,
    3.33462386036009E-6, 2.6848477528153E-6, 4.62422249743775E-6,
    4.7525801364159E-6, 4.57248902175055E-6, 4.20971936189489E-6,
    2.56702569900114E-6, 2.2657413182539E-6, 2.22709235244141E-6,
    2.98397317135631E-6, 3.04035834774704E-6, 3.0707500843955E-6,
    3.9301729497091E-6, 5.45108447865708E-6, 5.96474503344837E-6,
    5.9162964908761E-6, 3.30575072923888E-6, 4.31046925706741E-6,
    1.54004112034711E-6, 2.26406956060344E-6, 2.40829448842333E-6,
    2.18090038435304E-6, 5.44689046316766E-6, 7.72858586382364E-6,
    7.99943736350214E-6, 8.00369245200577E-6, 5.93053629819045E-6,
    2.0973204523194E-6, 2.75909958903236E-6, 2.93756639370372E-6,
    2.8082729761367E-6, 2.86298618955252E-6, 5.89668919427287E-6,
    6.48632782989695E-6, 6.64673723243134E-6, 6.43345350184602E-6,
    2.61631968382151E-6, 4.24775916673217E-6, 2.82623405753764E-6,
    2.84221333291443E-6, 5.18042138596141E-6, 4.91727865844577E-6,
    4.9853727310287E-6, 8.03533357600893E-6, 2.14157653936554E-6,
    2.7170783083722E-6, 2.88156361195859E-6, 3.96112234852728E-6,
    4.54948743720181E-6, 4.54551790191018E-6, 4.70845083025563E-6,
    2.38273495326574E-6, 2.40767640473199E-6, 2.77932729487705E-6,
    2.49983400889361E-6, 2.76683048483914E-6, 3.1061365688251E-6,
    6.04797898634501E-6, 9.21471961394426E-6, 9.20105485371253E-6,
    8.82500120858895E-6, 8.24679704366118E-6, 3.93650007101678E-6,
    3.66804069279849E-6, 4.32302210567006E-6, 4.50474033573741E-6,
    4.51127439274362E-6, 7.19204644275096E-6, 1.72898015359726E-5,
    2.02736428987894E-5, 1.75340847040288E-5, 1.51112130067923E-5,
    6.434834505982E-6, 4.93925239431223E-6, 1.83165226810988E-6,
    2.57944370588958E-6, 3.03063003142828E-6, 3.57765505289657E-6,
    8.22574566069632E-6, 1.9705212497918E-5, 2.14820988916525E-5,
    2.18034029163354E-5, 1.93156978928301E-5, 3.26030649522624E-6,
    3.0150141666075E-6, 1.92151311128738E-6, 3.36593597612278E-6,
    4.63947622024273E-6, 1.36234793967722E-5, 2.35139468786882E-5,
    2.63255281486665E-5, 2.62851397391829E-5, 1.81722396005852E-5,
    9.27917037612926E-6, 5.16521305992622E-6, 4.2934991758809E-6,
    5.0206722974244E-6, 1.32099546892382E-5, 4.22546081257112E-5,
    4.82313974794809E-5, 4.91840153026967E-5, 4.23206757256908E-5,
    5.24910960186829E-6, 5.40527510591431E-6, 5.03297731479864E-6,
    5.04893416803668E-6, 1.78769663626489E-5, 1.99155086605319E-5,
    2.62165217968964E-5, 2.57436171026296E-5, 1.68019550700417E-5,
    6.8072911447162E-6, 3.76497307706352E-6, 3.9536935426878E-6,
    2.52469946440819E-6, 3.45700422240666E-6, 2.60803949422759E-6,
    3.38295934331604E-6, 3.95860007638396E-6, 1.25713511926258E-5,
    1.43101456373531E-5, 1.50582053274483E-5, 1.58005265866164E-5,
    7.10332547981235E-6, 6.47979074145994E-6, 5.4089630183411E-6,
    2.41731736696558E-6, 7.44537912876952E-6, 1.23232121882883E-5,
    1.18032590933668E-5, 1.16940260269531E-5, 9.20098318138249E-6,
    8.46466377903573E-6, 2.75364707711935E-6, 2.54151423578185E-6,
    2.28609599507081E-6, 2.49661871260943E-6, 4.0157073856288E-6,
    3.94217587863838E-6, 5.95252984616858E-6, 7.3607819004567E-6,
    7.45084221241554E-6, 7.30831543744514E-6, 5.019678561542E-6,
    3.94599586239302E-6, 4.53595364611592E-6, 3.15436721633163E-6,
    3.14063700669077E-6, 2.07254506496741E-6, 3.37414796140676E-6,
    2.89390761108874E-6, 4.1927930911195E-6, 5.72785107713753E-6,
    1.09327871267255E-5, 1.13301194330376E-5, 1.28698683316394E-5,
    1.11832364344021E-5, 8.84345150160976E-6, 8.12833351935908E-6,
    7.35849133925843E-6, 3.40731285257518E-6, 2.9342528933812E-6,
    8.55259783775368E-6, 1.34538841261291E-5, 1.97519142736575E-5,
    2.02108898099328E-5, 1.7739191970555E-5, 9.1965734795521E-6,
    2.3586450480869E-6, 2.31348899633416E-6, 2.10146130470924E-6,
    2.68175869974829E-6, 2.56544297916926E-6, 2.04747060937373E-6,
    2.0513678627121E-6, 2.43245768009547E-6, 3.05035287182558E-6,
    9.06102258642226E-6, 2.23118807457135E-5, 2.32291757831723E-5,
    2.26992312500823E-5, 1.55704297787567E-5, 3.11039589126316E-6,
    2.51625948934958E-6, 2.39983458534733E-6, 3.05486125617132E-6,
    2.46176327353772E-6, 2.37555925102595E-6, 2.44123299751192E-6,
    2.67726931007984E-6, 2.18836627187516E-6, 6.98643293367651E-6,
    8.42642292838454E-6, 8.9745899206758E-6, 8.26191830701239E-6,
    2.9021901581308E-6, 2.72472701618836E-6, 2.59100875620985E-6,
    2.00034314546292E-6, 2.52890600765142E-6, 2.28277200787814E-6,
    2.22293181792091E-6, 1.83374914635594E-6, 2.65162054851593E-6,
    5.15239849481277E-6, 5.43727795366418E-6, 5.79105827902518E-6,
    5.58749168720752E-6, 2.87254836165249E-6, 7.8055604481712E-6,
    1.64724424743214E-6, 2.68087612012021E-6, 2.10922364612082E-6,
    3.02380244686545E-6, 7.2271985413888E-6, 8.24866906346741E-6,
    8.71452201362507E-6, 6.06619195415513E-6, 5.17708143345852E-6,
    2.68692151576354E-6, 2.55230905399963E-6, 2.93291266385739E-6,
    3.59757108269936E-6, 4.29749360351553E-6, 4.85665768790896E-6,
    5.0204214322455E-6, 4.87767193492799E-6, 4.84843480964456E-6,
    4.75110123899422E-6, 2.24893049692103E-6, 2.38805518568181E-6,
    3.24871560619148E-6, 3.56943619414011E-6, 4.81785746863827E-6,
    4.95914464880755E-6, 4.40730429713608E-6, 4.057553659123E-6,
    6.49724098230692E-6, 5.57531325893907E-6, 4.93501915022181E-6,
    5.35861722812005E-6, 6.19709122964692E-6, 6.52093363489576E-6,
    1.29222737431597E-5, 6.19116963522851E-6, 1.21781554864597E-5,
    7.74117645627206E-6, 1.43565980807126E-5, 7.30314723404081E-6,
    9.63362033284774E-6, 2.33616271999951E-5, 2.43540366197727E-5,
    2.39001612169181E-5, 6.93796267805559E-6, 7.24117792026927E-6,
    7.24691107335949E-6, 3.85637753932707E-6, 4.96836709319677E-6,
    4.88759335209377E-6, 1.06516122260278E-5, 9.40239583503366E-6,
    2.74798887655952E-6, 3.81625450499504E-5, 5.3367379557213E-5,
    5.46082946310087E-5, 5.49298405125502E-5, 7.98803044882501E-6,
    8.52165135911324E-6, 5.2484178970597E-6, 9.91922098926479E-6,
    4.79698616995163E-6, 5.75162319917013E-6, 6.04183355031744E-6,
    2.79415101259386E-5, 4.09809463050913E-5, 4.16792894882756E-5,
    4.18146187931366E-5, 7.26774281212742E-6, 6.26205464639477E-6,
    5.26856489099579E-6, 4.55216823147672E-6, 4.97899984997823E-6,
    4.78381461004219E-6, 4.51966844415477E-6, 6.4276274122182E-6,
    5.89255350376837E-6, 8.53132002871114E-6, 2.59542703905714E-5,
    2.59349511025224E-5, 2.66206069438013E-5, 2.03763577376454E-5,
    2.38800709671214E-5, 6.88672971264075E-6, 7.77089831958184E-6,
    5.16855273131121E-6, 4.57365575619857E-6, 3.90967207842104E-6,
    3.65629543501935E-6, 5.236480023763E-6, 4.27711374115347E-6,
    8.47048163883327E-6, 5.19309980785037E-5, 5.0246682805287E-5,
    5.02324775156243E-5, 3.81828528010698E-5, 2.73424916065309E-5,
    4.84292415874782E-6, 5.22770063456556E-6, 5.62621647592243E-6,
    3.41119054679705E-6, 4.1576509080615E-6, 3.95420361036366E-6,
    4.56442382845235E-5, 7.54824715928262E-5, 0.00010123253868761,
    0.000104538338471362, 0.000104206073796819, 1.21485757635351E-5,
    9.67627208666509E-6, 6.72655064196589E-6, 6.174370650798E-6,
    7.7572109671655E-6, 4.55994585130607E-6, 5.46711108693508E-6,
    3.07568062573429E-6, 5.06570704969015E-6, 7.64392974906049E-6,
    3.71650317894715E-5, 4.12367343878451E-5, 4.10202982382243E-5,
    2.25767028115493E-5, 2.19956152696231E-5, 3.82145282880158E-6,
    4.24445799073634E-6, 3.11679992435969E-6, 3.35998441915744E-6,
    3.99416503167664E-6, 4.21639372621868E-6, 7.45812057636487E-6,
    2.81142233524613E-5, 2.71955007625006E-5, 2.8645360083078E-5,
    1.97110004807364E-5, 1.71581799826111E-5, 1.81822845716247E-5,
    1.7171702338535E-5, 3.64035095955178E-6, 3.27962463591771E-6,
    4.64214177545632E-6, 5.08958384835415E-6, 5.34889995842117E-6,
    4.58781168706062E-6, 5.20672554247118E-6, 9.23676405454544E-5,
    0.000152113429821276, 0.00015658871712336, 0.00015659580675269,
    1.18805605443372E-5, 9.47108910285353E-6, 6.71627557109795E-6,
    6.76903088326515E-6, 5.90594828271037E-6, 8.08014566810009E-6,
    4.85461364118103E-6, 1.41166177064312E-5, 0.000138996098071969,
    0.000150613175662492, 0.000151572482285941, 0.000143340601539063,
    9.95083129028804E-6, 5.36997583634072E-6, 4.75277255813184E-6,
    5.56909647435659E-6, 5.52838201968956E-6, 5.82456909312089E-6,
    4.58795806039094E-6, 6.16751886348198E-6, 3.34024998608034E-5,
    3.49914715122679E-5, 3.48272611295227E-5, 2.6665489550309E-5,
    1.13063705627368E-5, 5.6850614659476E-6, 5.8966265720548E-6,
    6.79289357246483E-6, 7.14725466052425E-6, 6.53019089477715E-6,
    7.04340792240698E-6, 3.99877826769221E-5, 5.20957946452661E-5,
    5.01505336042303E-5, 5.20564121088345E-5, 4.88946104626163E-5,
    9.27311494290596E-6, 7.44190468129479E-6, 6.69618879753395E-6,
    5.75341415954831E-6, 6.05103057651948E-6, 1.07510857020707E-5,
    5.18598164588266E-6, 4.91165839301956E-6, 2.26872346069572E-5,
    3.36292605785285E-5, 3.49445876674337E-5, 3.49714534380113E-5,
    2.77623513072805E-5, 8.29990868742412E-6, 6.34580983058899E-6,
    5.21763981083155E-6, 5.66875948255217E-6, 6.06583194263955E-6,
    6.17682183680051E-6, 2.70862209720474E-5, 3.02841538977845E-5,
    2.96422332958949E-5, 2.53209359714757E-5, 2.30794197594517E-5,
    4.12876818673174E-6, 6.20269275267156E-6, 6.42662630893302E-6,
    6.07183544684181E-6, 1.37629920950919E-5, 3.55053478373647E-5,
    4.15063723604297E-5, 5.40494905962939E-5, 5.38555242311796E-5,
    5.27615141853472E-5, 6.51188824699612E-6, 3.41280544301843E-6,
    6.71831715540151E-6, 5.6188015450176E-6, 5.90720844081941E-6,
    5.81291424258095E-6, 5.59354859739537E-6, 1.1171964339657E-5,
    3.09396998888494E-5, 2.98810252621286E-5, 3.01632346085162E-5,
    0.000209136424369127, 0.000216296599187877, 0.00021671411302036,
    0.000211095451605146, 4.34574817232293E-6, 5.09500159018561E-6,
    4.91539903118372E-6, 1.16994514124199E-5, 3.45861863455921E-6,
    2.74865757691666E-6, 1.46142710255463E-5, 6.51148864456919E-5,
    0.000108994052903273, 0.000114533494705246, 0.000119706471118365,
    0.000103907901712145, 4.61060489546027E-6, 4.73336902977617E-6,
    4.26113138719364E-6, 1.19585358964516E-5, 3.91049164295347E-6,
    9.35920300603747E-6, 7.83072171955969E-5, 8.35626481978095E-5,
    8.42130779337982E-5, 8.55640437405677E-5, 4.89397619204077E-6,
    4.5194732524579E-6, 2.72544806398556E-6, 4.02019126329054E-6,
    4.23086470792102E-6, 3.46236143897508E-6, 6.26347477656949E-6,
    5.84442813237076E-6, 2.01039556729732E-5, 5.58161166555034E-5,
    5.76467466239449E-5, 5.76798860066296E-5, 1.07289619132944E-5,
    7.66357733404598E-6, 4.02260964048006E-6, 4.43046786342527E-6,
    4.13401858524485E-6, 3.88330109958483E-5, 6.41932456611498E-5,
    6.66282280937767E-5, 6.6502924590802E-5, 8.6967584961857E-6,
    7.31918670928619E-6, 6.65087316833222E-6, 5.37694055147469E-6,
    4.45823081228274E-6, 4.60559719414139E-6, 8.57560859123639E-6,
    5.75219452107916E-5, 6.01129052333837E-5, 5.99416119085187E-5,
    0.000109299376138281, 5.49311302877608E-6, 4.93384625435662E-6,
    5.09878860409546E-6, 5.19549559957886E-6, 4.6373233320983E-6,
    4.22990923883767E-5, 5.80022602183109E-5, 8.57297285290993E-5,
    8.77586553548433E-5, 8.74669761944412E-5, 7.59719695711196E-6,
    5.16483561374695E-6, 6.03635740551193E-6, 4.59802567588235E-6,
    5.09795919948752E-6, 3.48744705783104E-5, 4.07268866513615E-5,
    8.75847403984222E-5, 9.05544841424338E-5, 8.82942987844395E-5,
    3.54664273927025E-5, 8.21918833276392E-6, 5.92697020018169E-6,
    6.04237865826674E-6, 4.51893412307595E-6, 6.82563549155264E-6,
    1.34806129043064E-5, 3.46283353628072E-5, 4.20564799802605E-5,
    4.12847669522544E-5, 4.31703158774049E-5, 4.25534649832476E-5,
    8.51384708382411E-6, 7.00290615662967E-6, 7.21684600705256E-6,
    6.86063031304797E-6, 3.40361406778059E-6, 4.64994589394835E-6,
    3.18435261997312E-5, 4.71551769690662E-5, 9.1006784022897E-5,
    9.90218341022881E-5, 9.93138146531009E-5, 9.4668070486284E-5,
    6.5919173726693E-6, 5.01894017649586E-6, 4.75447105850568E-6,
    3.61021420341598E-6, 4.17649548316822E-6, 2.04758519753245E-5,
    3.62521990576789E-5, 0.000271427315626173, 0.000284645127306594,
    0.000285034970977107, 0.000310964497069726, 9.83460805971317E-6,
    4.05312632425645E-6, 3.93814798964658E-6, 4.11476061586633E-6,
    2.07975877824903E-6, 1.00302773120724E-5, 6.5720616854555E-5,
    0.000106882605216936, 0.000173420553545146, 0.000179009646654362,
    0.000173888111861698, 1.79670694391353E-5, 5.7141638127032E-6,
    3.86164085436333E-6, 3.50481277254149E-6, 4.02969601054744E-6,
    4.63532565996357E-6, 9.41252048195846E-6, 6.72459977581509E-5,
    0.000227909602033918, 0.000291033071739385, 0.000296686039097583,
    0.000291226252649975, 1.60887440309365E-5, 6.95772530161442E-6,
    6.01115933000589E-6, 7.55050937469584E-6, 4.80147670732729E-6,
    3.45681775369496E-6, 8.7242786143204E-6, 8.78091026846871E-5,
    0.000245440179118192, 0.000292695304682217, 0.000295180713348594,
    0.000295300991091275, 1.87318825223258E-5, 1.1655201433804E-5,
    8.14902502712329E-6, 1.37530189782982E-5, 8.26231357983106E-6,
    6.9287679844347E-6, 2.32149428000404E-5, 8.97845513085354E-5,
    0.000106439471581566, 0.000157641545169013, 0.000161363271768833,
    0.000159530171241765, 1.39732726261122E-5, 7.94860782241163E-6,
    5.2586770684825E-6, 3.95693852170061E-6, 4.64421201657597E-6,
    4.26256123818673E-6, 8.49729344929414E-6, 5.8882710295439E-5,
    7.81025394392583E-5, 0.000111977553599711, 0.000115001479287373,
    0.000111334638570585, 8.261665386934E-6, 3.99141962580005E-6,
    4.31208037318594E-6, 3.75944658816157E-6, 4.77763386345566E-6,
    2.71362256088164E-6, 2.13076404434629E-5, 5.19747729418592E-5,
    6.40540620816021E-6, 6.74524749589366E-6, 5.82378155210042E-6,
    3.71918275437391E-6, 2.89536417869592E-6, 5.8014463169827E-6,
    4.83725815887117E-6, 4.58442423649749E-6, 5.78666480723449E-6,
    5.14186704087984E-6, 4.61476485500233E-6, 5.83062194330215E-6,
    5.50424165898677E-6, 6.85368855442333E-6, 7.48427448966708E-6,
    7.15727715353816E-6, 4.79799974860167E-6, 4.97115793884218E-6,
    4.95093561121635E-6, 4.85178278495927E-6, 5.12314649264284E-6,
    6.97000582782137E-6, 3.7364550539616E-6, 3.33758672942517E-6,
    2.35995695911071E-5, 1.84623008530142E-5, 1.86141120783453E-5,
    1.87175635920262E-5, 1.35842333250438E-5, 1.24390489210531E-5,
    1.26525536260806E-5, 7.7579404081639E-5, 5.64491993584525E-6,
    5.08947529085113E-6, 3.70668307034924E-6, 6.23306414216628E-6,
    3.62508013575445E-6, 3.25446577012214E-6, 3.47901302669878E-6,
    3.53755303431861E-6, 4.83992849894642E-6, 4.48855552408016E-6,
    4.43233036190311E-6, 4.58009143243129E-6, 6.13816749573013E-6,
    6.10062922451893E-6, 5.77079234809645E-6, 5.90641378922079E-6,
    8.3069206908447E-6, 7.9452289269128E-6, 8.59834190679807E-6,
    7.037075240102E-6, 6.33860488839392E-6, 7.04533064351905E-6,
    7.45690831360455E-6, 7.35036555131913E-6, 3.73785753508774E-6,
    3.56399491332087E-6, 2.6785705739429E-6, 1.41976272064573E-5,
    1.17594261538911E-5, 1.18593964627991E-5, 1.19108163913777E-5,
    4.71222354517315E-6, 1.03379823833552E-5, 5.98373440812713E-6,
    5.29386783200563E-6, 4.75347868651074E-6, 5.2120773897769E-6,
    5.72839001588184E-6, 6.53615867079585E-6, 6.54247949474019E-6,
    8.62850913765622E-6, 9.55454503408907E-6, 8.42502582911889E-6,
    4.99407458965782E-6, 3.41030217653342E-6, 4.27586149960863E-6,
    2.844977124132E-6, 4.30076646179627E-6, 2.75218961382979E-6,
    2.59444228028763E-6, 3.7983288287496E-6, 1.20501681200268E-5,
    1.38653909480483E-5, 1.24763509959662E-5, 8.93286553825126E-6,
    7.56599025983936E-6, 1.55150893387625E-5, 3.56883818432559E-6,
    4.01089980696473E-6, 4.42142680772314E-6, 4.11491850165092E-6,
    3.87940870821523E-6, 8.91520336630122E-6, 9.7794700516187E-6,
    1.08577243948464E-5, 1.13345268425251E-5, 1.20391112235549E-5,
    1.24234231892071E-5, 4.48542047341693E-6, 3.66606830443288E-6,
    3.03760018726393E-6, 5.36073767416486E-6, 5.36887619046861E-6,
    4.91152260079736E-6, 9.48575006578446E-6, 2.10811054094021E-5,
    2.12315713590178E-5, 2.11095490665643E-5, 1.05558173691531E-5,
    4.70721251767142E-6, 3.24552235614733E-6, 4.64652692355505E-6,
    6.00521719095922E-6, 5.55785657692419E-6, 2.10273838048052E-5,
    2.37173390561546E-5, 2.24377217592883E-5, 2.25024139267199E-5,
    2.24371087915406E-5, 5.38740562235213E-6, 5.51711273162555E-6,
    6.79778383915401E-6, 6.14522573825807E-6, 5.59642042207807E-6,
    2.48399005477907E-5, 2.68837645667894E-5, 2.80169789762662E-5,
    2.79924547165423E-5, 5.27607574359271E-6, 4.19451165973981E-6,
    3.88004189137328E-6, 5.75356127725779E-6, 2.01339175909416E-5,
    2.09858359287173E-5, 2.12718567468943E-5, 3.70915550365206E-5,
    7.17280501875678E-6, 7.78507096837385E-6, 7.26107020120995E-6,
    1.2220426613538E-5, 1.42834134872691E-5, 1.4584423660749E-5,
    1.69693397464921E-5, 4.00804947101037E-6, 4.42669500580785E-6,
    6.91039954822871E-6, 4.8274893604517E-6, 5.32035873367924E-6,
    7.42952991527384E-6, 1.67108943162674E-5, 2.39886795268205E-5,
    2.35716423184707E-5, 1.77341502827901E-5, 2.76222567801137E-5,
    9.53372815441649E-6, 8.88442185794868E-6, 8.69357653753683E-6,
    7.42740734628362E-6, 7.99573399174439E-6, 2.55048856620017E-5,
    5.07477181183261E-5, 5.93711265782484E-5, 5.81789535838185E-5,
    5.14710891052974E-5, 1.38494596854003E-5, 1.36875092852911E-5,
    3.37908179004377E-6, 6.18780154970008E-6, 5.98604585371575E-6,
    7.15496127503976E-6, 3.22432897455592E-5, 6.02755781655008E-5,
    7.22592634083997E-5, 7.34519981256829E-5, 7.31482364653597E-5,
    5.32694628010503E-6, 4.90812273110531E-6, 4.27287734465254E-6,
    6.32594324690121E-6, 6.49514300712212E-6, 4.26045617208248E-5,
    6.45355872351466E-5, 6.10255809004766E-5, 6.21369058661038E-5,
    5.35179995799992E-5, 1.43275565627489E-5, 8.58881022678234E-6,
    1.11077002072676E-5, 1.18701004942607E-5, 4.78179088015198E-5,
    9.95740185723038E-5, 0.000151064233915962, 0.000154064795129491,
    0.00015074489020773, 5.33313455589433E-6, 8.85669555274405E-6,
    9.44391781291524E-6, 1.03281335366453E-5, 6.24114901119111E-5,
    7.47036637714579E-5, 7.96885312068913E-5, 7.83003796154191E-5,
    5.66842403558565E-5, 1.47594710620327E-5, 8.11183979672871E-6,
    8.49174360070361E-6, 4.97785659802044E-6, 1.08264857751101E-5,
    5.56483301392658E-6, 7.83288700641976E-6, 6.35717140674608E-6,
    3.77318368121163E-5, 5.52198389408646E-5, 5.79223430360716E-5,
    5.80764953314816E-5, 1.76091561014473E-5, 1.83941516006741E-5,
    1.83323198353897E-5, 5.16761656639552E-6, 2.14900944673738E-5,
    4.29031652799175E-5, 4.34932775303028E-5, 4.38970919037267E-5,
    1.19726407611593E-5, 3.65502804497321E-5, 5.87606252709577E-6,
    5.86820456773064E-6, 5.95211872752937E-6, 5.32479182582175E-6,
    1.07734467053646E-5, 7.25120763155724E-6, 1.66618462138486E-5,
    2.07991753361014E-5, 1.98361412686559E-5, 1.83199788374481E-5,
    7.28526925848006E-6, 7.54258541952236E-6, 1.15433891933024E-5,
    8.15283434311577E-6, 5.90196637248849E-6, 6.82383804763112E-6,
    9.59012929936676E-6, 8.41848077491784E-6, 8.75525201262037E-6,
    1.20683584125403E-5, 2.98986889023054E-5, 3.02736341131153E-5,
    2.95549396786816E-5, 2.39585482153831E-5, 1.49903375699576E-5,
    1.63881710126907E-5, 1.65991241012376E-5, 1.03005265522255E-5,
    5.92070498044276E-6, 2.95605247344092E-5, 4.76432914098672E-5,
    5.87026424060971E-5, 6.10240199251903E-5, 6.26316255027551E-5,
    1.27996504766276E-5, 5.41603359846054E-6, 5.60250546505074E-6,
    5.45659157787847E-6, 2.49286274256702E-6, 3.49082595504641E-6,
    3.55396095321406E-6, 3.57191997554519E-6, 4.84107925199944E-6,
    6.83244051572996E-6, 3.21912083886969E-5, 7.19213598494889E-5,
    7.57308302901345E-5, 7.70722834806546E-5, 6.51848045348662E-5,
    9.49417171727564E-6, 5.21368620483016E-6, 5.34105335681219E-6,
    5.77707676361767E-6, 4.44773862504892E-6, 3.72951299737431E-6,
    3.85209294432279E-6, 4.41929029894167E-6, 2.9060277742712E-6,
    1.8808616081887E-5, 2.08880588165288E-5, 2.12422244467588E-5,
    2.33120654281007E-5, 6.9969711516192E-6, 6.37730749679586E-6,
    6.77150700089842E-6, 4.78300897748106E-6, 5.2598780909945E-6,
    5.63912885832038E-6, 4.83363064667272E-6, 5.24263561684768E-6,
    8.35576603909928E-6, 1.60972412139457E-5, 1.58187486894878E-5,
    1.76798206592251E-5, 1.71754481379855E-5, 3.53360978225421E-6,
    8.32429976116981E-6, 4.66336721391995E-6, 4.75014122538689E-6,
    5.67825503694859E-6, 6.97529842423174E-6, 2.3672005923804E-5,
    2.42599935801173E-5, 2.45090130271272E-5, 9.90404454983964E-6,
    1.45857690098089E-5, 9.28691451181283E-6, 4.90103260712302E-6,
    7.86668965375052E-6, -3.62832732425267E-6, -6.42819312073003E-6,
    -2.86149474174078E-6, -3.25816709300678E-6, -4.84988957960099E-6,
    -5.03746168576385E-6, -5.51439225150272E-6, -2.98943217560174E-6,
    -3.88321231455986E-6, -3.77323911313249E-6, -3.85253593142384E-6,
    -4.27273368616464E-6, -4.61595708819832E-6, -4.38926493693197E-6,
    -4.09144493559169E-6, -7.46918923526722E-6, -8.39588604047186E-6,
    -6.98657466572268E-6, -7.06918136965189E-6, -2.5595186781851E-5,
    -1.60047225114217E-5, -9.87823469233368E-6, -9.56576278284489E-6,
    -9.71293052106286E-6, -2.70754873047073E-5, -5.97900998570491E-6,
    -1.64984430196592E-5, -1.45311508089764E-5, -1.48991797473175E-5,
    -1.449465268277E-5, -1.31753095244828E-5, -8.63978046504566E-6,
    -4.35159395400282E-6, -5.10564070683507E-6, -6.85876289755758E-6,
    -5.385410677679E-6, -3.26679606416878E-6, -4.73515242514011E-6,
    -4.50979423534506E-6, -5.0677348460427E-6, -3.15641397507283E-5,
    -2.72173302911467E-5, -2.58674084719785E-5, -2.56594144700795E-5,
    -7.03223125577791E-6, -9.74026081882554E-6, -8.10227885299325E-6,
    -7.77199596276859E-6, -7.24533850665864E-6, -7.03765424576698E-6,
    -6.90238045243041E-6, -3.0363551787491E-5, -2.46720571818027E-5,
    -2.41601251394264E-5, -2.61087735864917E-5, -4.33703686810859E-6,
    -5.63488092999498E-6, -4.79998250767011E-6, -5.49621981269767E-6,
    -5.51527310769045E-6, -6.2853183820357E-6, -4.36503626731494E-6,
    -8.41907078729051E-6, -8.51711597440183E-6, -1.77907640311405E-5,
    -1.08927131796539E-5, -1.82969901529359E-5, -1.93086416814054E-5,
    -1.67982590942064E-5, -5.92521049620824E-6, -5.15178159077114E-6,
    -5.14522431154462E-6, -3.5827556896531E-6, -2.47142328742924E-6,
    -3.3246795925271E-6, -4.99046079832306E-6, -3.7751386999198E-6,
    -8.34937541800191E-6, -1.36044647961622E-5, -1.82144744695253E-5,
    -3.1067803242604E-5, -3.06454760873117E-5, -2.9355803436052E-5,
    -8.02357981102048E-6, -9.29233133343923E-6, -5.2645182464816E-6,
    -5.00277835402328E-6, -3.90201603138007E-6, -3.92508279977265E-6,
    -4.11480721847605E-6, -3.17349553216898E-5, -5.02618667790068E-5,
    -6.4041949469592E-5, -6.29829574353098E-5, -5.62017310414454E-5,
    -6.71612066392021E-5, -1.0342803808047E-5, -6.16619066803454E-6,
    -5.45029448070157E-6, -3.74891488335804E-6, -3.54923085826011E-6,
    -3.75762130948033E-6, -5.36685765767304E-6, -6.23803415476915E-6,
    -9.28437635459391E-6, -1.84270355178567E-5, -1.8524689298604E-5,
    -1.93950527363177E-5, -1.44407710852813E-5, -7.91669105178246E-6,
    -4.24202144806498E-6, -3.57766598020336E-6, -3.74035679968519E-6,
    -3.57190753709615E-6, -4.38738012672199E-6, -6.10998031149234E-6,
    -1.27649660291261E-5, -1.66700932300955E-5, -1.48474246546664E-5,
    -1.4174645976722E-5, -2.17112839359661E-5, -1.03767158545415E-5,
    -1.15875732107332E-5, -6.2018203042808E-6, -5.67753901547992E-6,
    -4.94108740172065E-6, -3.49458942163743E-6, -3.81130912424387E-6,
    -4.86720978400996E-6, -4.3690693778676E-6, -6.43521083460499E-6,
    -9.56335187616894E-5, -7.98753655842576E-5, -7.87567885879765E-5,
    -7.69221513535945E-5, -3.25339899817636E-5, -5.90532001350327E-6,
    -6.28239143883191E-6, -5.0705028704912E-6, -3.58525586413114E-6,
    -3.83766084269595E-6, -3.35285178336817E-6, -9.47552954125881E-6,
    -7.00825538875914E-5, -6.91858623194512E-5, -6.96045955607108E-5,
    -8.02547846387755E-5, -7.82319894252675E-6, -1.52980999930562E-5,
    -4.63678260658578E-6, -6.91028313183311E-6, -6.88548929385536E-6,
    -6.1462838891098E-6, -3.96689561578243E-6, -6.18787298096445E-6,
    -1.83483695918658E-5, -1.86050667501932E-5, -1.8784312391683E-5,
    -9.38308168189519E-6, -7.50276446442025E-6, -6.33850359651578E-6,
    -5.16344347579058E-6, -7.88376289262701E-6, -7.32627691570818E-6,
    -7.35141341254875E-6, -7.88604089224489E-6, -1.90242627397191E-5,
    -3.17500305064217E-5, -3.21955384846548E-5, -3.14504491132172E-5,
    -2.63308006556365E-5, -1.57114043339199E-5, -5.46582814486547E-6,
    -5.62100740855546E-6, -6.40554057170124E-6, -9.35954317170767E-6,
    -6.31321044801533E-6, -6.30542402697101E-6, -5.8277112851586E-6,
    -1.89539845356724E-5, -1.42119289093029E-5, -1.98623584266375E-5,
    -1.88001225461147E-5, -1.51534649921099E-5, -1.35908796638177E-5,
    -7.77691320213969E-6, -1.5840639732593E-5, -5.67560589652612E-6,
    -5.6425293792633E-6, -8.1893290496328E-6, -1.70406539530606E-5,
    -1.8245803732269E-5, -1.79106697628158E-5, -1.4106909370771E-5,
    -1.27380825389328E-5, -3.21611779009471E-6, -4.79937016818529E-6,
    -4.86110013772063E-6, -4.93069173188488E-6, -1.81841067910503E-5,
    -2.37963208470058E-5, -2.39114217954195E-5, -2.30893928056557E-5,
    -1.79419220609946E-5, -1.994457687007E-5, -2.11791586719551E-5,
    -4.56237186947544E-6, -9.31369641955931E-6, -4.27244600175162E-6,
    -4.2953533120344E-6, -4.6214236909305E-6, -3.60227254796995E-6,
    -9.36449712989079E-6, -1.44682356707173E-5, -1.83720128790332E-5,
    -2.78992784447297E-5, -9.74601849456878E-5, -9.27459594539994E-5,
    -9.34447862292546E-5, -5.19958774893848E-5, -8.18516887079217E-6,
    -6.65797434356959E-6, -6.51223815040645E-6, -6.91507383419934E-6,
    -2.72476889934474E-6, -2.71346461459722E-6, -1.85159191128217E-5,
    -2.93187467615432E-5, -6.197551473969E-5, -6.1901301833719E-5,
    -5.35714339440928E-5, -2.92739559833883E-5, -5.61886000439358E-6,
    -6.15755995951777E-6, -3.72172554979659E-6, -4.49530871879121E-6,
    -4.37400460289205E-6, -5.54415474587099E-6, -4.30909078167793E-5,
    -4.353377620435E-5, -4.38783875543233E-5, -3.46678212577011E-5,
    -5.54237735887952E-6, -3.99879395667307E-6, -3.89193510763911E-6,
    -3.88839929435232E-6, -2.94512882499994E-6, -3.50209579255611E-6,
    -4.62453733302628E-6, -4.8602590773286E-6, -2.45251577436259E-5,
    -3.06975745660398E-5, -2.98266511511597E-5, -2.97925921957928E-5,
    -3.70896208922688E-5, -6.32101217543599E-6, -3.90734996190923E-6,
    -4.53436489857856E-6, -4.72848139163666E-6, -3.46686242717017E-5,
    -3.40545762960916E-5, -3.31619087186669E-5, -3.25655362415915E-5,
    -2.6456187045577E-5, -6.83528958238384E-6, -3.96447282002102E-6,
    -9.13122198478781E-6, -9.00110834742093E-6, -8.98649982219857E-6,
    -1.14068491577988E-5, -2.89790023344912E-5, -2.79606408954246E-5,
    -2.80946136796404E-5, -2.03738956610972E-5, -3.53663720789395E-6,
    -3.29425013294828E-6, -2.8230948724434E-6, -3.23536721729547E-6,
    -3.15747106075555E-6, -2.47264528522648E-5, -5.54564524287627E-5,
    -4.11847341994079E-5, -4.11109076599655E-5, -4.13287019049577E-5,
    -1.87606004661996E-5, -7.45912005863449E-6, -8.26143429479406E-6,
    -6.10747558104988E-6, -3.53602990555796E-6, -1.93972670686244E-5,
    -2.4337132940865E-5, -4.9175156435505E-5, -4.81246455610124E-5,
    -4.82032111787324E-5, -2.19363095404323E-5, -5.00995898251699E-6,
    -4.53600790350429E-6, -4.03629254678019E-6, -5.6040538872538E-6,
    -8.59434301076824E-6, -7.30853476705367E-6, -2.05887738945217E-5,
    -3.13345897059636E-5, -2.93741905799273E-5, -2.94907079491736E-5,
    -2.86202759051839E-5, -1.81978002879399E-5, -5.37635876287457E-6,
    -5.01208068602773E-6, -6.12799124853519E-6, -5.86459788719991E-6,
    -5.22277184150103E-6, -2.86359966951539E-5, -3.38515579359103E-5,
    -5.75950468473174E-5, -5.81078692017515E-5, -5.50279091079637E-5,
    -3.23226799663165E-5, -8.42259041708598E-6, -5.19117981236969E-6,
    -7.01523405134472E-6, -3.88088102061525E-6, -5.34369943808769E-6,
    -1.99597991702382E-5, -4.7631959831116E-5, -0.000151120976567729,
    -0.000146942119146512, -0.000145376747460635, -8.368866602911E-5,
    -1.02398120090432E-5, -6.51673472679738E-6, -4.66486720076451E-6,
    -5.10194107744454E-6, -6.11731854013032E-6, -1.13241758549179E-5,
    -4.00730623427879E-5, -0.000110397750380373, -9.80017616852466E-5,
    -9.61733279548753E-5, -9.19720854267592E-5, -4.42668093442439E-5,
    -7.39947302318695E-6, -3.12725682971212E-6, -2.68986866034648E-6,
    -3.28906520629112E-6, -4.11523970123288E-6, -1.53572067049241E-5,
    -3.2437983540888E-5, -0.000172707738260233, -0.000170787342464723,
    -0.000162511449121592, -0.000175091381193721, -2.07235526467114E-5,
    -1.3136526643676E-5, -6.29642714766699E-6, -5.42795742965756E-6,
    -5.00308096225941E-6, -9.54512097003233E-6, -1.57967449938616E-5,
    -2.94753596816717E-5, -0.00018180517229045, -0.000187216056573951,
    -0.000189056850635704, -0.000125357790917227, -2.83088313360213E-5,
    -2.77627054487049E-5, -5.87197410049373E-6, -6.576060414206E-6,
    -8.19947218239726E-6, -6.27482389653755E-6, -3.74262598181813E-5,
    -4.32078362573408E-5, -9.89364176394123E-5, -8.36646646798506E-5,
    -7.9706501320143E-5, -8.34037451531321E-5, -2.80089248646535E-5,
    -7.70773732728221E-6, -6.75149607796964E-6, -4.08276045012039E-6,
    -7.37701943610176E-6, -4.41450911493339E-6, -1.17761766215501E-5,
    -3.08365696473898E-5, -6.64047471029104E-5, -6.27442695886414E-5,
    -6.76954247454917E-5, -5.50163882742936E-5, -1.87933338798594E-5,
    -3.78993687971273E-6, -5.26868173535933E-6, -3.1320441979689E-6,
    -3.49562808200577E-6, -1.84751652735224E-6, -1.23849237919771E-5,
    -2.4619810928113E-5, -5.17089439325705E-6, -1.12011948721114E-5,
    -5.23207957063204E-6, -6.04003898337484E-6, -3.84250648502937E-6,
    -5.65382599218594E-6, -5.1358940807609E-6, -7.85273125855203E-6,
    -5.4261573911015E-6, -4.07941944466396E-6, -9.06379506062008E-6,
    -5.36922574029695E-6, -6.06102087778365E-6, -6.52140322980285E-6,
    -6.76867078735266E-6, -5.27459837257351E-6, -4.63846185389806E-6,
    -4.52500527799303E-6, -4.57810732709765E-6, -8.22196939505578E-6,
    -5.60291477324712E-6, -4.06835480388925E-6, -4.65159778271642E-6,
    -3.98164875887441E-6, -1.90165041820792E-5, -2.15140179654201E-5,
    -2.12414378277917E-5, -1.22005967829719E-5, -1.84151198606202E-5,
    -1.89052741571227E-5, -1.8624054179956E-5, -1.17106733118352E-5,
    -7.0190657676318E-6, -6.21760232879406E-6, -3.53375217862702E-6,
    -4.33458540484055E-6, -2.19262401211155E-6, -4.08409276825997E-6,
    -4.37937311897479E-6, -4.61190997805534E-6, -5.7326011723676E-6,
    -6.79645884016865E-6, -6.14351067490928E-6, -4.01265608866492E-6,
    -4.33284306258094E-6, -4.35068390866589E-6, -5.27899859025014E-6,
    -6.06255933380222E-6, -5.0292460326516E-6, -4.94026553792926E-6,
    -4.4584399843337E-6, -5.15657076988699E-6, -7.47796638830846E-6,
    -7.49067284446495E-6, -8.1862824458509E-6, -5.02332187399664E-6,
    -4.54013548581205E-6, -4.59125443956743E-6, -5.62513517738948E-6,
    -6.91800310643333E-6, -7.7145814685458E-6, -7.24491829288724E-6,
    -7.01090274466686E-6, -7.47326753447054E-6, -8.9295492109895E-6,
    -5.44628308284414E-6, -5.2312040310983E-6, -4.77279108116082E-6,
    -4.18143968615649E-6, -9.32492345381329E-6, -8.27623397163682E-6,
    -8.55718832674517E-6, -7.21130022074481E-6, -8.29991964941042E-6,
    -8.48200439094225E-6, -7.28938403352932E-6, -5.40842780051245E-6,
    -2.53932082079705E-6, -2.95249768773367E-6, -4.5385419757237E-6,
    -3.5125205916571E-6, -3.12070605253877E-6, -5.29328213772263E-6,
    -9.87186398875761E-6, -1.02524522120549E-5, -1.06987354102612E-5,
    -8.81105116121902E-6, -7.65886309523605E-6, -5.31937585719603E-6,
    -3.3611384043556E-6, -3.77696294818744E-6, -4.13602259404235E-6,
    -4.25811260783071E-6, -4.61700433204001E-6, -9.56745810558212E-6,
    -4.17678700030306E-6, -9.73818133174744E-6, -9.49069439271623E-6,
    -8.31564675816148E-6, -7.09357904728939E-6, -5.7169376699121E-6,
    -4.73219306734979E-6, -5.74556736050426E-6, -6.33137439182145E-6,
    -6.87939538820691E-6, -7.5643231199347E-6, -1.02901949755952E-5,
    -9.02641056770153E-6, -8.57250378792546E-6, -8.60460030872288E-6,
    -7.07550373843216E-6, -1.47906017638034E-5, -3.50925094199102E-6,
    -6.16328030722339E-6, -4.17015591548634E-6, -4.45310971476395E-6,
    -7.47411293317457E-6, -1.84050536891022E-5, -1.63817010640176E-5,
    -1.66652820905262E-5, -1.3884553839997E-5, -4.45874213218176E-6,
    -4.5328448859695E-6, -5.6700211001568E-6, -5.87188324245722E-6,
    -6.59982455163526E-6, -1.08673315249989E-5, -1.12098021421461E-5,
    -9.72140437728763E-6, -9.45413554041113E-6, -7.66682471919212E-6,
    -1.38615661090651E-5, -6.87795980573505E-6, -6.05450302873474E-6,
    -7.95985070862359E-6, -7.26241169714473E-6, -7.1214434886503E-6,
    -7.2313320381201E-6, -3.48844033505743E-6, -3.61988701998848E-6,
    -5.50244327316812E-6, -7.72999937363397E-6, -5.72790422092546E-6,
    -7.09843654530386E-6, -7.36744793874095E-6, -6.8325708075304E-6,
    -7.03910514137076E-6, -7.16812613433661E-6, -7.07160469823664E-6,
    -7.77860095343374E-6, -7.76209691707997E-6, -1.33270413026214E-5,
    -1.954166766158E-5, -1.96041428088502E-5, -1.83724496279249E-5,
    -1.22507189017692E-5, -7.02506609070943E-6, -7.32285594770707E-6,
    -9.11750236741813E-6, -1.11427981636017E-5, -1.09315546113375E-5,
    -1.11067278966177E-5, -3.54160478556506E-5, -3.38353896043705E-5,
    -3.18972998855099E-5, -3.7633711887074E-5, -1.26126371716298E-5,
    -8.69051874101983E-6, -4.97518229246094E-6, -4.73073355157871E-6,
    -6.68953704843519E-6, -9.20457799694775E-6, -1.20743974865094E-5,
    -4.58434978681065E-5, -4.16349205851129E-5, -4.18249744608868E-5,
    -3.96093309495226E-5, -8.21367976996022E-6, -1.04996947472116E-5,
    -3.55211574037956E-6, -7.11886271923962E-6, -1.01190240375585E-5,
    -2.12650474796383E-5, -4.6536913410563E-5, -4.83691002011198E-5,
    -4.75557257465743E-5, -3.7908236561602E-5, -2.86642480936511E-5,
    -1.44745575668027E-5, -7.63412841968139E-6, -1.15290604403657E-5,
    -2.29788639667883E-5, -0.000101738100989733, -8.48084301670209E-5,
    -8.47411568616962E-5, -6.77341175798612E-5, -1.93576338692719E-5,
    -1.25117828950966E-5, -1.25541189393058E-5, -1.16576981480663E-5,
    -3.75623808991638E-5, -3.47570627984982E-5, -3.19963400740322E-5,
    -3.04637648989914E-5, -2.86424317921381E-5, -1.91667792733217E-5,
    -7.17884906401209E-6, -6.9735663539915E-6, -6.03758292254913E-6,
    -6.25711947041277E-6, -6.74220715121918E-6, -7.3036458094739E-6,
    -1.12860494358863E-5, -3.53108657331139E-5, -2.57002375409129E-5,
    -2.21463712666932E-5, -2.10712535179088E-5, -1.18527170433721E-5,
    -9.48886092055916E-6, -8.38328197018201E-6, -4.90037459080157E-6,
    -2.13159923211614E-5, -2.57706160444269E-5, -2.07919439113322E-5,
    -2.04663363269162E-5, -3.52563356013857E-5, -7.48903700357838E-6,
    -8.43221295075208E-6, -4.53363457637147E-6, -4.34960157300642E-6,
    -5.5968325432403E-6, -7.44703865124906E-6, -8.21261182344985E-6,
    -9.82359009297809E-6, -1.0572433948886E-5, -1.1343794207173E-5,
    -1.54704136102395E-5, -1.16431938620353E-5, -1.1405746618623E-5,
    -1.13519366890112E-5, -5.20805454453657E-6, -1.09420743112205E-5,
    -3.47752171319599E-6, -4.85982144356475E-6, -5.22287776137299E-6,
    -8.76181636579607E-6, -1.06834875652884E-5, -1.71255340789903E-5,
    -1.68062840078815E-5, -2.35416446083005E-5, -1.91158198315596E-5,
    -1.62318116958033E-5, -1.52432408374042E-5, -2.23072320596203E-5,
    -6.2851825493599E-6, -6.564039431088E-6, -2.51129122632114E-5,
    -2.06710733140536E-5, -3.33702535591516E-5, -3.29105709835637E-5,
    -2.9705185966356E-5, -2.81604281399298E-5, -5.0506729676205E-6,
    -5.04861474800631E-6, -3.99163362052132E-6, -1.02313651513283E-5,
    -8.04631234204354E-6, -3.7711887092825E-6, -5.43199463527512E-6,
    -6.94565123346414E-6, -7.92345402668373E-6, -1.41047400085224E-5,
    -4.24723359801255E-5, -4.13498449829165E-5, -3.911524181983E-5,
    -1.74353833661171E-5, -8.65776899648507E-6, -5.18765105045317E-6,
    -5.32110746333797E-6, -5.64980195365149E-6, -5.32059068397447E-6,
    -5.38005342736191E-6, -5.45988753433982E-6, -7.49514924502177E-6,
    -5.2936849339314E-6, -2.04202895201545E-5, -2.36165468844441E-5,
    -2.42318248981833E-5, -2.14665843676285E-5, -4.71960735948972E-6,
    -4.86389612218165E-6, -5.20756577652859E-6, -3.60059549585979E-6,
    -4.46176496987432E-6, -3.71379579244588E-6, -7.50723097568673E-6,
    -4.40132783196195E-6, -5.76939506914906E-6, -1.10999808958132E-5,
    -1.30183563825762E-5, -1.22487276674299E-5, -1.24114423988509E-5,
    -9.83995424540267E-6, -3.24969360766381E-5, -3.1271606966572E-6,
    -8.00357737049397E-6, -2.89699683520839E-6, -6.56198756264668E-6,
    -1.53608654756005E-5, -2.03962636604068E-5, -1.90177392285211E-5,
    -1.67502533706942E-5, -1.10095761008048E-5, -3.88675981683517E-6,
    -4.87492859172088E-6, -4.19669433578091E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 17.0,
    0.0, 0.0, 0.0, 20.0, 0.0, 23.0, 14.0, 17.0, 17.0, 12.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 29.0, 38.0, 48.0, 48.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 28.0, 37.0, 35.0, 35.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 12.0, 20.0, 25.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 11.0, 62.0, 49.0, 51.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 37.0,
    75.0, 108.0, 113.0, 87.0, 33.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    23.0, 54.0, 58.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 18.0, 26.0,
    26.0, 26.0, 4.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 39.0, 76.0,
    90.0, 84.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60.0, 86.0, 88.0, 60.0,
    0.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 23.0, 31.0, 31.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 33.0, 62.0, 81.0, 83.0, 49.0, 19.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 21.0, 26.0, 46.0, 54.0, 22.0, 10.0, 0.0, 17.0, 0.0, 0.0,
    0.0, 20.0, 42.0, 41.0, 24.0, 14.0, 0.0, 0.0, 0.0, 0.0, 16.0, 24.0, 70.0,
    90.0, 53.0, 51.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 27.0, 39.0, 62.0,
    113.0, 95.0, 98.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 61.0, 104.0,
    123.0, 101.0, 46.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 53.0, 71.0, 72.0, 52.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18.0, 36.0, 27.0, 27.0, 30.0, 0.0,
    0.0, 0.0, 0.0, 30.0, 54.0, 71.0, 73.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0,
    30.0, 49.0, 53.0, 37.0, 0.0, 0.0, 0.0, 0.0, 0.0, 37.0, 68.0, 98.0, 88.0,
    74.0, 18.0, 0.0, 0.0, 0.0, 0.0, 38.0, 56.0, 107.0, 111.0, 83.0, 25.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 26.0, 60.0, 62.0, 39.0, 43.0, 18.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 30.0, 50.0, 85.0, 103.0, 76.0, 44.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0,
    50.0, 59.0, 97.0, 81.0, 55.0, 4.0, 0.0, 0.0, 0.0, 0.0, 7.0, 55.0, 71.0,
    124.0, 102.0, 84.0, 39.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 56.0, 71.0, 108.0,
    75.0, 57.0, 35.0, 10.0, 0.0, 0.0, 0.0, 0.0, 11.0, 65.0, 29.0, 67.0, 70.0,
    48.0, 47.0, 20.0, 0.0, 0.0, 0.0, 0.0, 18.0, 54.0, 81.0, 127.0, 101.0, 81.0,
    33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 56.0, 82.0, 120.0, 88.0, 66.0, 21.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 13.0, 50.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    30.0, 26.0, 26.0, 8.0, 18.0, 13.0, 13.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 5.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 14.0, 13.0, 13.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 4.0,
    0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 23.0, 30.0, 46.0, 9.0, 0.0,
    0.0, 0.0, 10.0, 10.0, 9.0, 48.0, 76.0, 52.0, 14.0, 13.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 18.0, 42.0, 61.0, 64.0, 43.0, 0.0, 4.0, 0.0, 0.0, 2.0, 54.0, 84.0, 92.0,
    98.0, 84.0, 42.0, 14.0, 0.0, 10.0, 49.0, 49.0, 83.0, 90.0, 78.0, 18.0, 11.0,
    11.0, 8.0, 66.0, 51.0, 96.0, 107.0, 69.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    7.0, 31.0, 35.0, 50.0, 68.0, 17.0, 0.0, 0.0, 0.0, 21.0, 38.0, 39.0, 34.0,
    25.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 14.0, 21.0, 11.0, 8.0,
    7.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 6.0, 32.0, 42.0, 82.0, 64.0, 47.0, 39.0,
    16.0, 0.0, 0.0, 20.0, 47.0, 72.0, 76.0, 61.0, 36.0, 0.0, 0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 26.0, 58.0, 74.0, 76.0, 34.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 12.0, 18.0, 18.0, 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 8.0, 14.0, 12.0, 12.0, 0.0, 22.0, 0.0, 0.0, 0.0, 0.0, 12.0,
    18.0, 17.0, 14.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 23.0, 22.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 29.0, 29.0, 27.0,
    15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 33.0, 27.0, 28.0, 17.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 18.0,
    18.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.0, 18.0, 18.0, 16.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 19.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 12.0, 13.0, 9.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 6.0, 4.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.0, 10.0, 10.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.0,
    25.0, 25.0, 26.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 17.0, 23.0, 25.0, 20.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 33.0, 36.0, 27.0, 27.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 37.0, 38.0, 40.0, 28.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 4.0, 17.0, 26.0, 18.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 12.0, 25.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 17.0, 17.0, 17.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.55985115349262E-5, -0.000109009369310881, 4.57065639346384E-5,
    4.59826788243839E-6, -2.22437285674695E-5, -1.60700841818448E-5,
    3.376095614036E-5, 1.52092871263923E-5, -4.17062850445728E-5,
    -4.96307036057834E-5, -7.44698055543479E-5, 2.39349448651928E-5,
    9.22783828319296E-6, 3.20499947328342E-5, -5.74335417868234E-5,
    6.92487588673511E-5, -8.56081600161013E-5, -1.24912922115002E-5,
    -7.06111102575528E-5, -0.000360805311217588, -0.000147735312460892,
    0.000152988732924093, 1.42399623439779E-5, 0.000187236325619947,
    -0.000376576367034171, 0.000187775525228067, -0.000281173674226119,
    -6.60017316727674E-5, 4.99853036292593E-5, 3.91320773780938E-5,
    -8.7507247814593E-6, -0.000172950335857869, 8.50796308482241E-5,
    -2.00083058084694E-5, -9.37158226835222E-5, 8.71657764944583E-5,
    -1.31633363066082E-7, 0.000107373816892528, 0.00013912560502043,
    2.92278096671567E-5, 2.70421414080196E-6, -8.69795111351953E-5,
    -5.99339684585455E-5, -7.14044143333748E-5, -9.46761572769038E-5,
    -0.000178339830259628, -0.000107569858956417, 0.000102920766222964,
    -1.13799353987682E-5, -4.47861141321267E-5, 3.58815481344256E-5,
    -3.39122739865786E-5, 2.11379319384506E-5, -5.41879947642099E-5,
    -9.75869694358711E-5, -5.46955798068776E-5, -4.09801229685324E-5,
    -1.52982108083519E-5, 3.20216987801497E-5, 2.56471947617937E-5,
    -2.94484343215426E-5, -2.14568218829804E-5, 7.66890807368119E-5,
    4.46194593840378E-6, 2.48271303055878E-5, -6.05882557245734E-5,
    -2.60654744045065E-5, -0.0001258383783393, 2.68557705641615E-5,
    0.000281791644936045, -6.64061490201921E-5, 7.18517229257593E-5,
    5.95239018838072E-6, 5.07343042801081E-5, -7.37345071163475E-7,
    -6.03161286840239E-5, -3.76979796297485E-5, -0.00010591286489778,
    1.21859522089919E-5, 4.07175399984001E-5, 0.000109427676101342,
    -0.00015375384912604, -0.000322371613715936, 0.000284457984632843,
    -0.000138068238374067, 7.68027948060255E-6, -4.18977341688417E-5,
    -3.34425016865999E-5, 4.30306419068798E-5, -1.30327194531882E-5,
    5.61335269203582E-5, 4.08812266485077E-6, -1.70591852888957E-5,
    6.47292761985292E-5, -0.000413719291166533, -0.00107401354405059,
    -0.000114536031751878, -4.72904052320401E-5, -1.01313051304795E-5,
    9.73981153647395E-5, 4.32626714306753E-6, 3.36511823807333E-5,
    -4.74778206491932E-5, 3.32749233880402E-5, 5.1108950162142E-5,
    -9.69792370977189E-5, 5.42387021902441E-5, -0.000103068566909393,
    -0.000212539785635063, 0.000290231885465444, 2.05699684837761E-5,
    -4.5656897760919E-5, -5.73939037494708E-5, -1.59166896213087E-5,
    3.07770624785563E-5, -1.13206496539605E-5, -9.55667168124465E-5,
    -9.91313638841464E-5, -4.26921104437983E-5, -0.000187023527100891,
    -0.000401599822788601, 0.000143335600378957, 5.74623190034711E-5,
    0.000128961547872232, -2.27828615500877E-5, -5.51175950715078E-5,
    -2.56882974944718E-5, -6.4269121142214E-5, 6.31510607887944E-5,
    1.22743152797635E-5, -6.40945556274128E-5, 2.9035425643814E-5,
    -4.77699683101777E-5, 7.63763383042864E-5, -0.000217827148516344,
    -0.000622104057248757, 0.000117662962381402, -5.84047081281385E-5,
    3.45243786192757E-5, 7.28212284935219E-5, 5.98055905344255E-5,
    -5.1744891061358E-6, 2.72676273163541E-5, 3.70413459883458E-5,
    -3.14562263380898E-5, -0.000128680310703135, 0.000404320996846256,
    -0.000109507064697142, -0.000114170796972528, 3.64965944969471E-5,
    3.13001341313522E-5, -8.3093472246317E-5, 6.99809500668954E-6,
    -2.96475929647875E-5, 5.88561395958966E-5, 0.000109180737569836,
    5.46145753265728E-5, -9.44439602390571E-5, 0.000190571518152833,
    4.59840527868465E-5, -6.87937634811285E-5, 1.77636309019932E-5,
    -7.64197916751059E-5, -0.000116646633843294, -1.55561169990572E-5,
    5.47284225921176E-5, 5.09066999465082E-5, -3.03851810815431E-5,
    -0.000123430957602276, 0.000169202921207164, -2.72120587197033E-5,
    -0.000242916337369051, -2.5386012927943E-5, 5.28249214496962E-5,
    6.13228018981227E-5, -0.000120889062416887, 9.50635775629634E-5,
    -5.49207795198914E-5, 7.95132850617036E-5, -1.59402617093077E-5,
    2.52085937151432E-5, -3.49308914964985E-5, -7.16125494310056E-6,
    -0.000136429538719547, -7.30877527141592E-5, -5.53943326391684E-6,
    -0.000228771293702302, 2.32436555758277E-5, -1.51437599059712E-5,
    3.12668022916221E-6, -3.39365969932603E-5, 6.1891341219838E-5,
    -4.37124006735178E-5, -0.000258987092095069, 4.4926885912451E-5,
    -2.95785976498677E-6, 6.74031115938194E-5, -4.04639073196617E-5,
    3.3598781620633E-5, -1.73420090398852E-5, 2.00706739455979E-5,
    2.18882043909985E-5, -0.000135095312305049, -0.000159120138504474,
    -9.69053887271975E-5, -0.000264186470548089, -3.96576022595033E-5,
    -0.000103977574896918, 3.26625336191766E-5, 9.04436193440572E-6,
    -2.43980136158688E-5, -2.23589256043436E-5, 5.79140480262693E-5,
    6.67662826290824E-6, 3.30974589467753E-6, 4.85643202928362E-6,
    1.30462417971257E-5, 0.000145798325481526, -0.000641459841746467,
    0.0019933933647124, -0.000152161334457356, -3.26682645819439E-5,
    1.36697414288214E-5, 9.38137457579203E-5, 2.86797579914266E-5,
    -2.52282074935531E-7, -5.75691212881682E-5, -6.48200020092495E-5,
    0.000102816587363928, -0.000164830898574337, -0.000820323602362805,
    0.000921182955337841, -9.12874045459993E-5, -5.68312064062876E-5,
    -1.62875129367265E-5, 0.000147354163615051, 7.24189943592972E-5,
    0.000144819758664092, 0.000133654273987425, -1.17187908933591E-5,
    -8.22428933613079E-5, 0.000570540226559682, -9.08457430269171E-5,
    4.4039326832628E-5, -3.21567619033615E-5, -5.02586515396484E-6,
    4.27136459594776E-5, -1.0259104672821E-5, -3.08539178017436E-5,
    -4.22861673152845E-6, 1.94079212657426E-5, -5.52646052911176E-5,
    -6.65983266107484E-5, -0.0001299733759545, -0.000613064051691269,
    9.1238042588384E-5, -4.90768956591499E-5, 2.05017078648248E-5,
    4.43652565862443E-5, -0.000116816186444249, -7.66257396681945E-5,
    1.77794011889245E-5, -0.000157288734087305, -0.000493222401046109,
    9.61919135653039E-5, 5.71666440287865E-5, -5.69322121869818E-5,
    -6.54123243551359E-5, 7.40098753368048E-6, -2.52700521680312E-5,
    7.49414297571794E-5, -4.39629557606876E-5, -3.95716572346806E-5,
    0.00120535700104355, -2.79540260692688E-5, -3.92907797586989E-5,
    3.12471657148818E-5, 1.19328638730274E-5, -8.17774997207662E-6,
    6.51934249555938E-5, 1.72137828696465E-5, -0.00011321705565882,
    0.00043780901476233, -0.000247382952000251, -0.000327741367278609,
    -8.29054862281073E-5, -0.000109448568000276, -9.29590477996963E-7,
    2.66256500954763E-6, 5.49503299096092E-5, 3.14354758048658E-5,
    -0.00012037984629443, 0.000105271148409204, -0.000195661279265197,
    1.31385648186623E-5, -5.93795931521326E-5, 4.8144327539852E-7,
    9.99967130968912E-5, -4.88476991945361E-5, 8.80919271141632E-6,
    0.000164905744580893, 1.13199354426257E-5, 1.6625079147935E-6,
    7.6375560137502E-5, 0.000331369888229607, -0.000250363129414236,
    -0.000345244598769424, 1.75137933154058E-5, -5.54727095562416E-5,
    -3.37682975605581E-5, -6.27338037486978E-5, -3.75274595767462E-5,
    -0.000108261799124404, -2.37996193703062E-5, 0.000159970707502927,
    -0.000117947485163548, -0.000356693642759777, 0.000749336894981465,
    -4.93760536068905E-5, 1.62603318471301E-5, -8.33704855505988E-5,
    -4.18456794911411E-5, -6.94776291237592E-5, 2.75228112549634E-5,
    -1.76129033036673E-5, 2.91161919316464E-5, 5.45509673005289E-5,
    -0.00101523259878682, 0.00364342984855039, -0.000160846760408255,
    -5.38472708877987E-5, -6.54283171120827E-5, 2.43948800387349E-5,
    -8.56270063714942E-5, -1.38250658271661E-5, 1.58511710165826E-5,
    2.24175902572988E-5, -0.000188802059350943, 0.000920043632422046,
    -0.000384432192282491, -0.000840648210936845, -7.57535905429291E-5,
    9.13707291863339E-6, -2.78527121613176E-5, 6.74981177014312E-6,
    8.90520183960782E-6, 8.70729390525621E-5, 4.74168838226886E-5,
    0.000253589563232489, -0.000222173208512235, -0.00045120741651554,
    0.000986234540867152, -0.000512518126262467, -0.000146312555353678,
    7.49726930891263E-5, 2.78827146599621E-5, 3.40322583854725E-5,
    -8.59005866130937E-5, -4.42570221761604E-5, -2.62677979225623E-5,
    0.000293591888574106, -5.85515750177421E-5, 0.000369037183357296,
    0.00253179183136765, -0.000602623962272384, -0.00025638403468601,
    -7.33355683266472E-5, 0.00012725042464764, 4.26595408409134E-5,
    -3.65804731962651E-6, 4.36296042941502E-5, -0.000170972748499493,
    0.000164160863523444, -0.000232855539284059, 0.000458341217291923,
    -0.000371561912354231, -0.000527962181015112, -5.70285351014475E-5,
    -6.41913763932463E-5, 1.01069762688032E-5, -9.74092606882373E-5,
    -2.69215076490265E-5, -5.96267840959306E-5, 6.68676236662453E-5,
    7.05145147951093E-5, -0.000124351218311792, 0.00118129463179898,
    -6.7404298682073E-5, -0.000342272989601104, 8.77197674162789E-6,
    -7.82144531222921E-5, 2.94728361922114E-5, 3.16743495010559E-5,
    2.29791247435891E-5, 8.91781400060151E-5, 1.12699555738382E-5,
    5.20028760912573E-5, -0.000175961210573015, -8.90670994584749E-6,
    -8.67660321589039E-5, -4.53893627908293E-6, 6.63415743019405E-5,
    2.53493652596609E-5, -0.000123410521066378, -7.70705088130965E-5,
    -1.83673778395099E-5, -8.44021027331352E-5, -1.12275263861787E-5,
    4.62240643705099E-6, 1.39884493536065E-5, 1.85095415928486E-5,
    -4.57991442944805E-5, 1.97999570554609E-5, 3.01020675793577E-5,
    -4.69261470656003E-5, -8.00946465176629E-5, -2.95730069036935E-5,
    7.21678956947901E-5, -1.58337863454024E-5, 3.01020826152532E-5,
    -7.38170019547137E-6, 9.80410565182062E-5, -0.000116970110588088,
    0.00033827155177363, -5.46658618107226E-5, -0.000116594584186546,
    -4.10576364701508E-5, 0.00092165198759865, -0.000167265840682314,
    -6.55700456130861E-5, 1.21593287878327E-5, 4.17353465147218E-5,
    1.84657704224687E-5, -4.09363437240981E-5, -8.58556575907527E-6,
    1.10531813962953E-5, 5.86257927322211E-5, -1.76605514697537E-5,
    4.45091942282788E-5, -2.58811964993925E-5, -5.27625504672161E-5,
    5.18264634694427E-7, 2.4327362502162E-5, 1.40903830080749E-5,
    -2.28173235969643E-5, 2.26997139370869E-5, 8.39467458728117E-5,
    1.95594569743015E-5, -9.1714618114299E-5, -0.000124912038803736,
    -5.35884305714225E-5, 1.74636809757899E-5, 3.05869379312561E-5,
    3.02870464433977E-5, -5.94174834506495E-5, -1.88143559137068E-5,
    2.18156055584044E-6, 1.65070357190288E-5, 2.04978443239143E-5,
    -4.19276386420468E-5, 0.000123893308291604, 0.000139774003334635,
    1.20781825964099E-5, 2.68471027471587E-6, 4.9104542208558E-5,
    6.62196786200816E-5, -6.70897341132847E-5, 1.8238833124441E-5,
    -3.74657085486297E-5, -1.58267083140583E-5, 8.07521912794163E-5,
    -3.98846905316652E-5, -3.90486252928186E-5, 7.04113443293353E-5,
    -1.08994490273814E-5, 7.6530588774686E-5, 2.01069772001305E-5,
    -2.23574473970756E-5, -7.28489220987955E-5, 4.25492900515457E-5,
    2.96082938933284E-5, 3.44989065506725E-5, -6.19775592406192E-5,
    -1.718180794906E-5, 0.000143535827694031, 1.17796739389225E-5,
    -3.23619519430589E-5, 3.21775988305877E-5, 7.09304414750715E-5,
    -5.23247145730137E-5, 2.41766450757279E-5, 1.4016591518748E-6,
    -6.53066783011751E-5, 1.73970681455207E-5, 7.92853148474428E-6,
    8.82597678775161E-5, 2.81046817151577E-5, -3.86870381449958E-5,
    -3.46718337326534E-5, -4.54781240261082E-5, 3.60482682014554E-5,
    1.41140744478245E-5, -7.76479537348258E-5, -2.19261108392175E-5,
    4.32602998775663E-5, -1.22768564592497E-5, -8.92510762353532E-5,
    -0.000230960973797202, -1.62428831187702E-5, 2.4591947338344E-7,
    7.16632004516289E-6, 1.172709176049E-6, -6.41880496653572E-5,
    4.00673659264039E-5, -4.68450541063858E-5, 0.000128238355021182,
    1.87449699188693E-5, 2.51730375839311E-5, -4.39898004760286E-5,
    5.45530010051699E-5, 8.3708312477079E-5, 8.40310817426451E-5,
    2.0757545585588E-5, -4.72253639995745E-5, -1.21931867367455E-5,
    2.77082708592233E-5, -4.86497877556122E-5, -0.000208638819955353,
    -9.27783328104835E-5, -3.69502588825128E-5, -8.89813836536687E-5,
    -2.26761479230465E-6, -7.35275376427778E-5, 0.000371241085018761,
    7.5596174782627E-5, 4.82189409757652E-5, -3.27019468693628E-5,
    -2.72821876505565E-5, 5.30628874285511E-5, 3.10718390181898E-5,
    0.000117044956790562, -5.51991721076217E-5, 6.0693698573143E-5,
    5.72763261917656E-5, -5.72569787669031E-5, -1.53015841375449E-5,
    5.82756135942191E-5, -1.79013417007232E-5, -6.3785820365108E-5,
    -2.41659934732377E-5, -0.000402136338594117, 0.000458629107653646,
    2.07468128644345E-5, -7.85175400341918E-5, -0.000132598964135397,
    -2.61531070953716E-5, -1.84845606725637E-5, -2.82497231269891E-5,
    7.58584288341571E-5, -3.37525840690407E-5, -0.000184555731060176,
    0.000111601575809037, 5.92633525459095E-5, 0.000184023503933112,
    8.21712113247176E-6, 6.43815370903126E-5, 5.40027030456665E-5,
    1.12331077358394E-5, 8.10908033302855E-5, 8.21967380602687E-5,
    6.22464305190958E-5, 0.000154713248443489, -2.54854194268543E-5,
    -2.38233611623088E-5, -7.44022626595037E-5, -1.24482815585029E-5,
    1.89938547375786E-5, -9.85393995040764E-5, 0.000110478981036936,
    -5.12749034221657E-6, -1.38700692342889E-5, 0.000295963768593391,
    -0.000155984896971662, -0.000465063898900171, 1.41530235936513E-5,
    6.80208392212862E-5, -0.00013288975833002, 4.61984414080804E-5,
    -4.78430660071671E-5, -5.03210807368182E-5, 0.000315270945681791,
    -9.94817974280249E-5, -0.000331965605220653, -0.00018193487960156,
    -3.19914223579307E-5, -2.41825818925781E-5, -0.000137132962470787,
    8.95311045519511E-6, -1.64413617686559E-5, 9.6031435272398E-5,
    -0.000238778029214358, -0.000261875929718427, 1.97308245648021E-5,
    -2.41599570789707E-5, 2.50687293250363E-5, 0.000119984504905377,
    -6.0487219364524E-5, 3.71130091666366E-5, -0.00018122289836131,
    -1.29931325715102E-5, -0.000136275539941182, 1.48610908481787E-5,
    -3.75333645331361E-5, -4.43298591604322E-5, 0.000108292487685947,
    -1.40159092143085E-6, -8.73471847062137E-6, -8.55917658204204E-5,
    -4.28288493029599E-5, 0.000111663072184186, -2.92561696181375E-5,
    -0.000481509804574914, 0.00044149152569559, -9.37237645644068E-5,
    6.92347894279247E-5, -1.58244470673278E-5, -1.11321006168758E-5,
    5.65828872237397E-6, -0.000116968476801453, -2.75557556856896E-5,
    5.82498385512582E-5, -0.000174865928323874, -7.77130219154503E-6,
    -0.000111394566507705, -6.31605862740359E-5, 0.000185173795442734,
    0.000133637147717297, -0.00015511143384555, 2.05552558571651E-5,
    0.000146669360137554, 3.14491232929518E-5, 0.000139727844114733,
    -3.60430706464162E-5, -4.15605570383211E-5, 5.09969089738745E-5,
    -0.000388261286930688, 8.98456824332055E-5, -8.29829124348053E-5,
    -0.000260713813201554, 3.33056259761169E-5, 0.000125835755428848,
    -0.000105876464729781, -8.93544839082206E-5, -3.55294128420902E-5,
    -3.98003256977639E-5, -2.56172270669912E-5, -0.000201558512037932,
    -0.000597875097073922, -4.21106466799166E-5, -4.48392484959181E-5,
    6.01245433938899E-5, -0.000129967632515684, -9.68295783900478E-5,
    3.58671149315304E-5, -4.80617750382445E-5, 4.83225280080969E-5,
    -2.40116492673398E-5, -4.22682858259973E-5, 7.97791545236589E-5,
    8.61947189205942E-6, -0.000237954367422859, 0.000564070186938103,
    -5.77489032947224E-5, -5.45395427583604E-5, -4.46067553393868E-5,
    -6.50001835312645E-5, -1.26828102264218E-7, 3.28620751929251E-5,
    -3.7359201230177E-5, -0.000109615013785388, -4.63470590925585E-5,
    -1.87135128821438E-5, 5.27589383764098E-5, -8.74679681225696E-5,
    -0.000233737895015852, 8.60215074683868E-5, -8.40759052659422E-5,
    7.06726196871885E-5, 4.28912305435538E-5, -7.49611752959564E-5,
    1.13412804322399E-5, -8.23080170140794E-5, -5.79422183861728E-5,
    5.29021205583756E-5, -3.26745919640207E-5, -2.39767533187131E-5,
    -4.56962543417932E-5, 2.87156793375714E-5, -9.56648509499969E-5,
    -0.000401899378474376, -6.17468132262393E-6, -0.000139686645514084,
    -1.47005759454119E-5, -9.72826084620796E-5, -2.41734689192378E-5,
    -5.13724386656498E-6, -7.6961906262739E-5, -0.000161454176384144,
    0.00014229147205248, 3.22389651929362E-5, -0.000100808727017233,
    0.000122689588610696, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0,
    2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 78.0, 78.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 82.0,
    82.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 95.0, 97.0, 97.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    90.0, 91.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 77.0, 78.0,
    78.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 102.0, 102.0, 102.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 93.0, 93.0, 95.0, 0.0, 0.0, 0.0, 0.0,
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
    0.0, 0.0, 73.0, 73.0, 73.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  double dist[712];
  int c_k;
  double s;
  signed char Y[3];
  emxArray_real_T *edges;
  int exitg4;
  int exitg3;
  emxArray_real_T *nn;
  short outsize_idx_0;
  boolean_T guard1 = false;
  int exitg2;
  boolean_T exitg1;
  for (k = 0; k <= 711; k += 2) {
    if (a[k] <= a[k + 1]) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  i = 2;
  while (i < 712) {
    i2 = i << 1;
    nb = 1;
    for (pEnd = 1 + i; pEnd < 713; pEnd = qEnd + i) {
      p = nb;
      q = pEnd;
      qEnd = nb + i2;
      if (qEnd > 713) {
        qEnd = 713;
      }

      k = 0;
      kEnd = qEnd - nb;
      while (k + 1 <= kEnd) {
        if (a[idx[p - 1] - 1] <= a[idx[q - 1] - 1]) {
          iwork[k] = idx[p - 1];
          p++;
          if (p == pEnd) {
            while (q < qEnd) {
              k++;
              iwork[k] = idx[q - 1];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q - 1];
          q++;
          if (q == qEnd) {
            while (p < pEnd) {
              k++;
              iwork[k] = idx[p - 1];
              p++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx[(nb + k) - 1] = iwork[k];
      }

      nb = qEnd;
    }

    i = i2;
  }

  emxInit_real_T1(&Uc, 1);
  i2 = Uc->size[0];
  Uc->size[0] = 712;
  emxEnsureCapacity((emxArray__common *)Uc, i2, (int)sizeof(double));
  for (k = 0; k < 712; k++) {
    Uc->data[k] = a[idx[k] - 1];
  }

  nb = -1;
  k = 0;
  while (k + 1 <= 712) {
    i = (int)Uc->data[k];
    do {
      exitg5 = 0;
      k++;
      if (k + 1 > 712) {
        exitg5 = 1;
      } else {
        eok = (fabs((double)i - Uc->data[k]) < eps((double)i / 2.0));
        if (!eok) {
          exitg5 = 1;
        }
      }
    } while (exitg5 == 0);

    nb++;
    Uc->data[nb] = i;
  }

  i2 = Uc->size[0];
  if (1 > nb + 1) {
    b_i2 = -1;
  } else {
    b_i2 = nb;
  }

  Uc->size[0] = b_i2 + 1;
  emxEnsureCapacity((emxArray__common *)Uc, i2, (int)sizeof(double));
  for (i2 = 0; i2 < 712; i2++) {
    for (i = 0; i < 40; i++) {
      x[i2 + 712 * i] = tX[i2 + 712 * i] - tsX[i];
    }
  }

//#pragma omp parallel for \
// num_threads(omp_get_max_threads()) \
// private(c_k)

  for (b_k = 1; b_k < 28481; b_k++) {
    c_k = b_k;
    y[c_k - 1] = x[c_k - 1] * x[c_k - 1];
  }

  for (nb = 0; nb < 712; nb++) {
    s = y[nb];
    for (k = 0; k < 39; k++) {
      s += y[nb + (k + 1) * 712];
    }

    dist[nb] = s;
  }

  sort(dist, idx);
  for (i = 0; i < 3; i++) {
    Y[i] = a[idx[i] - 1];
  }

  emxInit_real_T(&edges, 2);
  i = Uc->size[0];
  i2 = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = (short)(i + 1);
  emxEnsureCapacity((emxArray__common *)edges, i2, (int)sizeof(double));
  k = 0;
  do {
    exitg4 = 0;
    i = Uc->size[0];
    if (k <= i - 2) {
      edges->data[1 + k] = Uc->data[k] + (Uc->data[1 + k] - Uc->data[k]) / 2.0;
      k++;
    } else {
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  edges->data[0] = rtMinusInf;
  edges->data[edges->size[1] - 1] = rtInf;
  k = 1;
  do {
    exitg3 = 0;
    i = Uc->size[0];
    if (k - 1 <= i - 2) {
      edges->data[k] += eps(edges->data[k]);
      k++;
    } else {
      exitg3 = 1;
    }
  } while (exitg3 == 0);

  emxInit_real_T1(&nn, 1);
  outsize_idx_0 = (short)edges->size[1];
  i2 = nn->size[0];
  nn->size[0] = outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)nn, i2, (int)sizeof(double));
  i = outsize_idx_0;
  for (i2 = 0; i2 < i; i2++) {
    nn->data[i2] = 0.0;
  }

  i = edges->size[1];
  guard1 = false;
  if (i > 1) {
    nb = 1;
    do {
      exitg2 = 0;
      if (nb + 1 <= i) {
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
    i2 = nn->size[0];
    nn->size[0] = outsize_idx_0;
    emxEnsureCapacity((emxArray__common *)nn, i2, (int)sizeof(double));
    i = outsize_idx_0;
    for (i2 = 0; i2 < i; i2++) {
      nn->data[i2] = rtNaN;
    }
  } else {
    i = 0;
    for (k = 0; k < 3; k++) {
      nb = findbin((double)Y[i], edges);
      if (nb > 0) {
        nn->data[nb - 1]++;
      }

      i++;
    }
  }

  i2 = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = nn->size[0] - 1;
  emxEnsureCapacity((emxArray__common *)edges, i2, (int)sizeof(double));
  for (k = 0; k <= nn->size[0] - 2; k++) {
    edges->data[k] = nn->data[k];
  }

  if (nn->size[0] - 1 > 0) {
    edges->data[edges->size[1] - 1] += nn->data[nn->size[0] - 1];
  }

  emxFree_real_T(&nn);
  i = 1;
  nb = edges->size[1];
  s = edges->data[0];
  i2 = 0;
  if (edges->size[1] > 1) {
    if (rtIsNaN(edges->data[0])) {
      pEnd = 2;
      exitg1 = false;
      while ((!exitg1) && (pEnd <= nb)) {
        i = pEnd;
        if (!rtIsNaN(edges->data[pEnd - 1])) {
          s = edges->data[pEnd - 1];
          i2 = pEnd - 1;
          exitg1 = true;
        } else {
          pEnd++;
        }
      }
    }

    if (i < edges->size[1]) {
      while (i + 1 <= nb) {
        if (edges->data[i] > s) {
          s = edges->data[i];
          i2 = i;
        }

        i++;
      }
    }
  }

  emxFree_real_T(&edges);
  yfit = Uc->data[i2];
  emxFree_real_T(&Uc);
  return yfit;
}

//
// Arguments    : const double x[4]
// Return Type  : double
//
static double mean(const double x[4])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 3; k++) {
    y += x[k + 1];
  }

  y /= 4.0;
  return y;
}

//
// Arguments    : int idx[712]
//                double x[712]
//                int offset
//                int np
//                int nq
//                int iwork[712]
//                double xwork[712]
// Return Type  : void
//
static void merge(int idx[712], double x[712], int offset, int np, int nq, int
                  iwork[712], double xwork[712])
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
// Arguments    : int idx[712]
//                double x[712]
//                int offset
//                int n
//                int preSortLevel
//                int iwork[712]
//                double xwork[712]
// Return Type  : void
//
static void merge_block(int idx[712], double x[712], int offset, int n, int
  preSortLevel, int iwork[712], double xwork[712])
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
// Arguments    : int idx[712]
//                double x[712]
//                int offset
// Return Type  : void
//
static void merge_pow2_block(int idx[712], double x[712], int offset)
{
  int iwork[256];
  double xwork[256];
  int b;
  int bLen;
  int bLen2;
  int nPairs;
  int k;
  int blockOffset;
  int q;
  int p;
  int exitg1;
  for (b = 0; b < 6; b++) {
    bLen = 1 << (b + 2);
    bLen2 = bLen << 1;
    nPairs = 256 >> (b + 3);
    for (k = 1; k <= nPairs; k++) {
      blockOffset = (offset + (k - 1) * bLen2) - 1;
      for (q = 1; q <= bLen2; q++) {
        iwork[q - 1] = idx[blockOffset + q];
        xwork[q - 1] = x[blockOffset + q];
      }

      p = 0;
      q = bLen;
      do {
        exitg1 = 0;
        blockOffset++;
        if (xwork[p] <= xwork[q]) {
          idx[blockOffset] = iwork[p];
          x[blockOffset] = xwork[p];
          if (p + 1 < bLen) {
            p++;
          } else {
            exitg1 = 1;
          }
        } else {
          idx[blockOffset] = iwork[q];
          x[blockOffset] = xwork[q];
          if (q + 1 < bLen2) {
            q++;
          } else {
            q = blockOffset - p;
            while (p + 1 <= bLen) {
              idx[(q + p) + 1] = iwork[p];
              x[(q + p) + 1] = xwork[p];
              p++;
            }

            exitg1 = 1;
          }
        }
      } while (exitg1 == 0);
    }
  }
}

//
// Arguments    : const emxArray_real_T *Y
//                const emxArray_real_T *iPk
//                emxArray_real_T *idx
// Return Type  : void
//
static void orderPeaks(const emxArray_real_T *Y, const emxArray_real_T *iPk,
  emxArray_real_T *idx)
{
  emxArray_real_T *x;
  int i14;
  int loop_ub;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_idx;
  if (idx->size[0] == 0) {
  } else {
    emxInit_real_T1(&x, 1);
    i14 = x->size[0];
    x->size[0] = idx->size[0];
    emxEnsureCapacity((emxArray__common *)x, i14, (int)sizeof(double));
    loop_ub = idx->size[0];
    for (i14 = 0; i14 < loop_ub; i14++) {
      x->data[i14] = Y->data[(int)iPk->data[(int)idx->data[i14] - 1] - 1];
    }

    emxInit_int32_T(&iidx, 1);
    emxInit_real_T1(&b_idx, 1);
    b_sort(x, iidx);
    i14 = b_idx->size[0];
    b_idx->size[0] = iidx->size[0];
    emxEnsureCapacity((emxArray__common *)b_idx, i14, (int)sizeof(double));
    loop_ub = iidx->size[0];
    emxFree_real_T(&x);
    for (i14 = 0; i14 < loop_ub; i14++) {
      b_idx->data[i14] = idx->data[iidx->data[i14] - 1];
    }

    emxFree_int32_T(&iidx);
    i14 = idx->size[0];
    idx->size[0] = b_idx->size[0];
    emxEnsureCapacity((emxArray__common *)idx, i14, (int)sizeof(double));
    loop_ub = b_idx->size[0];
    for (i14 = 0; i14 < loop_ub; i14++) {
      idx->data[i14] = b_idx->data[i14];
    }

    emxFree_real_T(&b_idx);
  }
}

//
// Arguments    : const emxArray_creal_T *x
//                int xoffInit
//                int unsigned_nRows
//                const emxArray_real_T *costab
//                const emxArray_real_T *sintab
//                emxArray_creal_T *y
// Return Type  : void
//
static void r2br_r2dit_trig_impl(const emxArray_creal_T *x, int xoffInit, int
  unsigned_nRows, const emxArray_real_T *costab, const emxArray_real_T *sintab,
  emxArray_creal_T *y)
{
  int j;
  int nRowsD2;
  int nRowsD4;
  int iy;
  int iDelta;
  int ix;
  int ju;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  double twid_re;
  double twid_im;
  int ihi;
  if (x->size[0] <= unsigned_nRows) {
    j = x->size[0];
  } else {
    j = unsigned_nRows;
  }

  nRowsD2 = unsigned_nRows / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
  iy = x->size[0];
  if (unsigned_nRows > iy) {
    iDelta = y->size[0];
    iy = y->size[0];
    y->size[0] = iDelta;
    emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
    for (iy = 0; iy < iDelta; iy++) {
      y->data[iy].re = 0.0;
      y->data[iy].im = 0.0;
    }
  }

  ix = xoffInit;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy] = x->data[ix];
    iDelta = unsigned_nRows;
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy] = x->data[ix];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= unsigned_nRows - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  iDelta = 2;
  iy = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += iy) {
      temp_re = y->data[i + iDelta].re;
      temp_im = y->data[i + iDelta].im;
      y->data[i + iDelta].re = y->data[i].re - temp_re;
      y->data[i + iDelta].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    ju = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      twid_re = costab->data[j];
      twid_im = sintab->data[j];
      i = ju;
      ihi = ju + ix;
      while (i < ihi) {
        temp_re = twid_re * y->data[i + iDelta].re - twid_im * y->data[i +
          iDelta].im;
        temp_im = twid_re * y->data[i + iDelta].im + twid_im * y->data[i +
          iDelta].re;
        y->data[i + iDelta].re = y->data[i].re - temp_re;
        y->data[i + iDelta].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += iy;
      }

      ju++;
    }

    nRowsD4 /= 2;
    iDelta = iy;
    iy <<= 1;
    ix -= iDelta;
  }
}

//
// Arguments    : const emxArray_real_T *Y
//                emxArray_real_T *iPk
//                double Ph
// Return Type  : void
//
static void removePeaksBelowMinPeakHeight(const emxArray_real_T *Y,
  emxArray_real_T *iPk, double Ph)
{
  int end;
  int trueCount;
  int i;
  int partialTrueCount;
  if (!(iPk->size[0] == 0)) {
    end = iPk->size[0] - 1;
    trueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y->data[(int)iPk->data[i] - 1] > Ph) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y->data[(int)iPk->data[i] - 1] > Ph) {
        iPk->data[partialTrueCount] = iPk->data[i];
        partialTrueCount++;
      }
    }

    end = iPk->size[0];
    iPk->size[0] = trueCount;
    emxEnsureCapacity((emxArray__common *)iPk, end, (int)sizeof(double));
  }
}

//
// Arguments    : const emxArray_real_T *Y
//                emxArray_real_T *iPk
//                double Th
// Return Type  : void
//
static void removePeaksBelowThreshold(const emxArray_real_T *Y, emxArray_real_T *
  iPk, double Th)
{
  int c;
  emxArray_real_T *base;
  int k;
  int trueCount;
  double extremum;
  int partialTrueCount;
  c = iPk->size[0];
  emxInit_real_T1(&base, 1);
  k = base->size[0];
  base->size[0] = c;
  emxEnsureCapacity((emxArray__common *)base, k, (int)sizeof(double));
  for (k = 0; k + 1 <= c; k++) {
    if ((Y->data[(int)(iPk->data[k] - 1.0) - 1] >= Y->data[(int)(iPk->data[k] +
          1.0) - 1]) || rtIsNaN(Y->data[(int)(iPk->data[k] + 1.0) - 1])) {
      extremum = Y->data[(int)(iPk->data[k] - 1.0) - 1];
    } else {
      extremum = Y->data[(int)(iPk->data[k] + 1.0) - 1];
    }

    base->data[k] = extremum;
  }

  k = iPk->size[0] - 1;
  trueCount = 0;
  for (c = 0; c <= k; c++) {
    if (Y->data[(int)iPk->data[c] - 1] - base->data[c] >= Th) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (c = 0; c <= k; c++) {
    if (Y->data[(int)iPk->data[c] - 1] - base->data[c] >= Th) {
      iPk->data[partialTrueCount] = iPk->data[c];
      partialTrueCount++;
    }
  }

  emxFree_real_T(&base);
  k = iPk->size[0];
  iPk->size[0] = trueCount;
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
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
// Arguments    : double u
// Return Type  : double
//
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

//
// scaleAbs :: scales to the maximum of the absolute value of the signal ...
//  (For a ratio of 1).
// Arguments    : emxArray_real_T *X
//                emxArray_real_T *Y
// Return Type  : void
//
static void scaleAbs(emxArray_real_T *X, emxArray_real_T *Y)
{
  int ixstart;
  int n;
  emxArray_real_T *y;
  unsigned int unnamed_idx_0;
  double mtmp;
  int ix;
  boolean_T exitg1;
  ixstart = X->size[0];
  n = X->size[0];
  X->size[0] = ixstart;
  emxEnsureCapacity((emxArray__common *)X, n, (int)sizeof(double));
  emxInit_real_T1(&y, 1);
  unnamed_idx_0 = (unsigned int)X->size[0];
  n = y->size[0];
  y->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)y, n, (int)sizeof(double));
  for (ixstart = 0; ixstart + 1 <= X->size[0]; ixstart++) {
    y->data[ixstart] = fabs(X->data[ixstart]);
  }

  ixstart = 1;
  n = y->size[0];
  mtmp = y->data[0];
  if (y->size[0] > 1) {
    if (rtIsNaN(y->data[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(y->data[ix - 1])) {
          mtmp = y->data[ix - 1];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < y->size[0]) {
      while (ixstart + 1 <= n) {
        if (y->data[ixstart] > mtmp) {
          mtmp = y->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxFree_real_T(&y);
  n = Y->size[0];
  Y->size[0] = X->size[0];
  emxEnsureCapacity((emxArray__common *)Y, n, (int)sizeof(double));
  ixstart = X->size[0];
  for (n = 0; n < ixstart; n++) {
    Y->data[n] = X->data[n] / mtmp;
  }
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
// Arguments    : double x[712]
//                int idx[712]
// Return Type  : void
//
static void sort(double x[712], int idx[712])
{
  double x4[4];
  short idx4[4];
  int m;
  double xwork[712];
  int nNaNs;
  int ib;
  int k;
  signed char perm[4];
  int iwork[712];
  int i2;
  int i3;
  int i4;
  memset(&idx[0], 0, 712U * sizeof(int));
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  memset(&xwork[0], 0, 712U * sizeof(double));
  nNaNs = -711;
  ib = 0;
  for (k = 0; k < 712; k++) {
    if (rtIsNaN(x[k])) {
      idx[-nNaNs] = k + 1;
      xwork[-nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = (short)(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        ib = (k - nNaNs) - 714;
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

  m = (nNaNs + 711) >> 1;
  for (k = 1; k <= m; k++) {
    ib = idx[k - nNaNs];
    idx[k - nNaNs] = idx[712 - k];
    idx[712 - k] = ib;
    x[k - nNaNs] = xwork[712 - k];
    x[712 - k] = xwork[k - nNaNs];
  }

  if (((nNaNs + 711) & 1) != 0) {
    x[(m - nNaNs) + 1] = xwork[(m - nNaNs) + 1];
  }

  memset(&iwork[0], 0, 712U * sizeof(int));
  m = 2;
  if (1 - nNaNs > 1) {
    ib = (1 - nNaNs) >> 8;
    if (ib > 0) {
      for (m = 1; m <= ib; m++) {
        merge_pow2_block(idx, x, (m - 1) << 8);
      }

      m = ib << 8;
      ib = 1 - (nNaNs + m);
      if (ib > 0) {
        merge_block(idx, x, m, ib, 2, iwork, xwork);
      }

      m = 8;
    }

    merge_block(idx, x, 0, 1 - nNaNs, m, iwork, xwork);
  }
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
          b_merge_block(idx, b_x, m, ib, 2, iwork, xwork);
        }

        m = 8;
      }
    }

    b_merge_block(idx, b_x, 0, nNonNaN, m, iwork, xwork);
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
// function: [s, f, t] = stft(x, wlen, h, nfft, fs)
//  x - signal in the time domain
//  wlen - length of the hamming window
//  h - hop size
//  nfft - number of FFT points
//  fs - sampling frequency, Hz
//  f - frequency vector, Hz
//  t - time vector, s
//  s - STFT matrix (only unique points, time across columns, freq across rows)
// Arguments    : const emxArray_real_T *x
//                double fs
//                emxArray_creal_T *s
//                double f[1025]
//                emxArray_real_T *t
// Return Type  : void
//
static void stft(const emxArray_real_T *x, double fs, emxArray_creal_T *s,
                 double f[1025], emxArray_real_T *t)
{
  int nm1d2;
  double twid_re;
  int b_x;
  int cdiff;
  double indx;
  int col;
  int exitg1;
  double xw[256];
  static const double win[256] = { 0.080000000000000016, 0.080138543399746076,
    0.0805540901456207, 0.081246389927802531, 0.082215025730789426,
    0.083459414084593453, 0.084978805416200731, 0.086772284501087038,
    0.08883877101451404, 0.091177020182276858, 0.093785623530509787,
    0.0966630097340977, 0.099807445563183939, 0.10321703692720313,
    0.10688973001581042, 0.11082331253602012, 0.11501541504480811,
    0.11946351237637592, 0.12416492516321609, 0.12911682145006298,
    0.13431621839975671, 0.13975998408999274, 0.14544483939987485,
    0.15136735998513473, 0.15752397834082921, 0.16391098595027154,
    0.17052453551890334, 0.1773606432917611, 0.184415191453141,
    0.19168393060701711, 0.19916248233671879, 0.2068463418423252,
    0.21473088065418816, 0.22281134942094921, 0.23108288077037153,
    0.23954049224126267, 0.2481790892847231, 0.2569934683329117,
    0.26597831993348064, 0.27512823194779118, 0.28443769281098297,
    0.29390109485193527, 0.30351273767111808, 0.31326683157429935,
    0.32315750106004104, 0.333178788358881, 0.34332465702207021,
    0.35358899555770473, 0.36396562111205871, 0.37444828319390544,
    0.38503066743957881, 0.39570639941650987, 0.40646904846294735,
    0.41731213156154678, 0.42822911724449858, 0.43921342952783993,
    0.450258451872581, 0.4613575311702614, 0.47250398175053365,
    0.48369108940836053, 0.49491211544840208, 0.506160300744153,
    0.51742886980938774, 0.52871103487946036, 0.54, 0.5512889651205396,
    0.56257113019061233, 0.573839699255847, 0.58508788455159788,
    0.59630891059163948, 0.60749601824946642, 0.61864246882973861,
    0.629741548127419, 0.64078657047216, 0.65177088275550144,
    0.66268786843845329, 0.67353095153705267, 0.68429360058349009,
    0.69496933256042115, 0.70555171680609463, 0.7160343788879413,
    0.72641100444229534, 0.73667534297792969, 0.74682121164111914,
    0.75684249893995892, 0.76673316842570072, 0.776487262328882,
    0.78609890514806469, 0.795562307189017, 0.80487176805220884,
    0.81402168006651943, 0.82300653166708826, 0.831820910715277,
    0.84045950775873723, 0.84891711922962854, 0.85718865057905091,
    0.86526911934581185, 0.87315365815767476, 0.88083751766328111,
    0.888316069392983, 0.89558480854685907, 0.90263935670823892,
    0.90947546448109673, 0.91608901404972853, 0.922476021659171,
    0.92863264001486534, 0.93455516060012522, 0.94024001591000728,
    0.94568378160024325, 0.950883178549937, 0.955835074836784,
    0.9605364876236242, 0.9649845849551919, 0.96917668746398,
    0.97311026998418959, 0.97678296307279688, 0.98019255443681619,
    0.98333699026590238, 0.98621437646949028, 0.98882297981772316,
    0.99116122898548609, 0.99322771549891309, 0.99502119458379934,
    0.99654058591540662, 0.99778497426921064, 0.99875361007219754,
    0.99944590985437931, 0.99986145660025394, 1.0, 0.99986145660025394,
    0.99944590985437931, 0.99875361007219754, 0.99778497426921064,
    0.99654058591540662, 0.99502119458379934, 0.99322771549891309,
    0.99116122898548609, 0.98882297981772316, 0.98621437646949028,
    0.98333699026590238, 0.98019255443681619, 0.97678296307279688,
    0.97311026998418959, 0.96917668746398, 0.9649845849551919,
    0.9605364876236242, 0.955835074836784, 0.950883178549937,
    0.94568378160024325, 0.94024001591000728, 0.93455516060012522,
    0.92863264001486534, 0.922476021659171, 0.91608901404972853,
    0.90947546448109673, 0.90263935670823892, 0.89558480854685907,
    0.888316069392983, 0.88083751766328111, 0.87315365815767476,
    0.86526911934581185, 0.85718865057905091, 0.84891711922962854,
    0.84045950775873723, 0.831820910715277, 0.82300653166708826,
    0.81402168006651943, 0.80487176805220884, 0.795562307189017,
    0.78609890514806469, 0.776487262328882, 0.76673316842570072,
    0.75684249893995892, 0.74682121164111914, 0.73667534297792969,
    0.72641100444229534, 0.7160343788879413, 0.70555171680609463,
    0.69496933256042115, 0.68429360058349009, 0.67353095153705267,
    0.66268786843845329, 0.65177088275550144, 0.64078657047216,
    0.629741548127419, 0.61864246882973861, 0.60749601824946642,
    0.59630891059163948, 0.58508788455159788, 0.573839699255847,
    0.56257113019061233, 0.5512889651205396, 0.54, 0.52871103487946036,
    0.51742886980938774, 0.506160300744153, 0.49491211544840208,
    0.48369108940836053, 0.47250398175053365, 0.4613575311702614,
    0.450258451872581, 0.43921342952783993, 0.42822911724449858,
    0.41731213156154678, 0.40646904846294735, 0.39570639941650987,
    0.38503066743957881, 0.37444828319390544, 0.36396562111205871,
    0.35358899555770473, 0.34332465702207021, 0.333178788358881,
    0.32315750106004104, 0.31326683157429935, 0.30351273767111808,
    0.29390109485193527, 0.28443769281098297, 0.27512823194779118,
    0.26597831993348064, 0.2569934683329117, 0.2481790892847231,
    0.23954049224126267, 0.23108288077037153, 0.22281134942094921,
    0.21473088065418816, 0.2068463418423252, 0.19916248233671879,
    0.19168393060701711, 0.184415191453141, 0.1773606432917611,
    0.17052453551890334, 0.16391098595027154, 0.15752397834082921,
    0.15136735998513473, 0.14544483939987485, 0.13975998408999274,
    0.13431621839975671, 0.12911682145006298, 0.12416492516321609,
    0.11946351237637592, 0.11501541504480811, 0.11082331253602012,
    0.10688973001581042, 0.10321703692720313, 0.099807445563183939,
    0.0966630097340977, 0.093785623530509787, 0.091177020182276858,
    0.08883877101451404, 0.086772284501087038, 0.084978805416200731,
    0.083459414084593453, 0.082215025730789426, 0.081246389927802531,
    0.0805540901456207, 0.080138543399746076 };

  creal_T X[2048];
  int ndbl;
  int i;
  int apnd;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int k;
  int j;
  static const double dv20[1025] = { 1.0, 0.99999529380957619,
    0.99998117528260111, 0.9999576445519639, 0.9999247018391445,
    0.99988234745421256, 0.9998305817958234, 0.99976940535121528,
    0.99969881869620425, 0.99961882249517864, 0.99952941750109314,
    0.99943060455546173, 0.99932238458834954, 0.99920475861836389,
    0.99907772775264536, 0.99894129318685687, 0.99879545620517241,
    0.99864021818026527, 0.99847558057329477, 0.99830154493389289,
    0.99811811290014918, 0.997925286198596, 0.99772306664419164,
    0.99751145614030345, 0.99729045667869021, 0.997060070339483,
    0.99682029929116567, 0.99657114579055484, 0.996312612182778,
    0.996044700901252, 0.99576741446765982, 0.99548075549192694,
    0.99518472667219693, 0.99487933079480562, 0.99456457073425542,
    0.9942404494531879, 0.99390697000235606, 0.9935641355205953,
    0.9932119492347945, 0.9928504144598651, 0.99247953459871, 0.9920993131421918,
    0.99170975366909953, 0.99131085984611544, 0.99090263542778,
    0.99048508425645709, 0.99005821026229712, 0.98962201746320089,
    0.989176509964781, 0.98872169196032378, 0.98825756773074946,
    0.98778414164457218, 0.98730141815785843, 0.98680940181418553,
    0.98630809724459867, 0.98579750916756748, 0.98527764238894122,
    0.98474850180190421, 0.984210092386929, 0.98366241921173025,
    0.98310548743121629, 0.98253930228744124, 0.98196386910955524,
    0.98137919331375456, 0.98078528040323043, 0.98018213596811743,
    0.97956976568544052, 0.9789481753190622, 0.97831737071962765,
    0.97767735782450993, 0.97702814265775439, 0.97636973133002114,
    0.97570213003852857, 0.97502534506699412, 0.97433938278557586,
    0.973644249650812, 0.97293995220556018, 0.97222649707893627,
    0.97150389098625178, 0.97077214072895035, 0.970031253194544,
    0.96928123535654853, 0.96852209427441727, 0.96775383709347551,
    0.96697647104485207, 0.9661900034454125, 0.9653944416976894,
    0.96458979328981276, 0.96377606579543984, 0.96295326687368388,
    0.96212140426904158, 0.96128048581132064, 0.96043051941556579,
    0.95957151308198452, 0.9587034748958716, 0.95782641302753291,
    0.95694033573220882, 0.95604525134999641, 0.95514116830577078,
    0.95422809510910567, 0.95330604035419386, 0.95237501271976588,
    0.95143502096900834, 0.9504860739494817, 0.94952818059303667,
    0.94856134991573027, 0.94758559101774109, 0.94660091308328353,
    0.94560732538052128, 0.94460483726148026, 0.94359345816196039,
    0.94257319760144687, 0.94154406518302081, 0.9405060705932683,
    0.93945922360218992, 0.93840353406310806, 0.937339011912575,
    0.93626566717027826, 0.93518350993894761, 0.93409255040425887,
    0.932992798834739, 0.93188426558166815, 0.93076696107898371,
    0.92964089584318121, 0.92850608047321559, 0.92736252565040111,
    0.92621024213831138, 0.92504924078267758, 0.92387953251128674,
    0.92270112833387863, 0.9215140393420419, 0.92031827670911059,
    0.91911385169005777, 0.9179007756213905, 0.9166790599210427,
    0.91544871608826783, 0.91420975570353069, 0.91296219042839821,
    0.91170603200542988, 0.91044129225806725, 0.90916798309052238,
    0.90788611648766626, 0.90659570451491533, 0.90529675931811882,
    0.90398929312344334, 0.90267331823725883, 0.901348847046022,
    0.90001589201616017, 0.89867446569395382, 0.89732458070541832,
    0.89596624975618522, 0.8945994856313827, 0.89322430119551532,
    0.89184070939234272, 0.89044872324475788, 0.88904835585466457,
    0.88763962040285393, 0.88622253014888064, 0.88479709843093779,
    0.88336333866573158, 0.881921264348355, 0.88047088905216075,
    0.87901222642863353, 0.87754529020726135, 0.8760700941954066,
    0.87458665227817611, 0.87309497841829009, 0.87159508665595109,
    0.87008699110871146, 0.8685707059713409, 0.86704624551569265,
    0.86551362409056909, 0.8639728561215867, 0.8624239561110405,
    0.86086693863776731, 0.85930181835700847, 0.85772861000027212,
    0.85614732837519447, 0.85455798836540053, 0.85296060493036363,
    0.8513551931052652, 0.84974176800085255, 0.84812034480329723,
    0.84649093877405213, 0.84485356524970712, 0.84320823964184544,
    0.84155497743689844, 0.83989379419599952, 0.83822470555483808,
    0.836547727223512, 0.83486287498638, 0.83317016470191319,
    0.83146961230254524, 0.829761233794523, 0.8280450452577558,
    0.82632106284566353, 0.82458930278502529, 0.82284978137582643,
    0.82110251499110465, 0.819347520076797, 0.81758481315158371,
    0.81581441080673378, 0.81403632970594841, 0.81225058658520388,
    0.81045719825259477, 0.808656181588175, 0.80684755354379933,
    0.80503133114296366, 0.80320753148064494, 0.80137617172314024,
    0.799537269107905, 0.79769084094339116, 0.79583690460888357,
    0.79397547755433717, 0.79210657730021239, 0.79023022143731,
    0.78834642762660634, 0.78645521359908577, 0.78455659715557524,
    0.78265059616657573, 0.78073722857209449, 0.778816512381476,
    0.77688846567323244, 0.77495310659487393, 0.773010453362737,
    0.77106052426181382, 0.7691033376455797, 0.7671389119358204,
    0.765167265622459, 0.76318841726338127, 0.76120238548426178,
    0.759209188978388, 0.75720884650648457, 0.75520137689653655,
    0.75318679904361252, 0.75116513190968637, 0.74913639452345937,
    0.74710060598018013, 0.745057785441466, 0.74300795213512172,
    0.74095112535495922, 0.73888732446061511, 0.73681656887736979,
    0.7347388780959635, 0.73265427167241282, 0.73056276922782759,
    0.7284643904482252, 0.726359155084346, 0.724247082951467,
    0.72212819392921535, 0.72000250796138165, 0.71787004505573171,
    0.71573082528381859, 0.71358486878079352, 0.71143219574521643,
    0.70927282643886569, 0.70710678118654757, 0.70493408037590488,
    0.7027547444572253, 0.70056879394324834, 0.69837624940897292,
    0.696177131491463, 0.69397146088965389, 0.69175925836415775,
    0.68954054473706683, 0.687315340891759, 0.68508366777270036,
    0.68284554638524808, 0.680600997795453, 0.67835004312986147,
    0.67609270357531592, 0.673829000378756, 0.67155895484701833,
    0.669282588346636, 0.66699992230363747, 0.66471097820334479,
    0.66241577759017178, 0.66011434206742048, 0.65780669329707864,
    0.65549285299961535, 0.65317284295377676, 0.650846684996381,
    0.64851440102211244, 0.64617601298331628, 0.64383154288979139,
    0.641481012808583, 0.63912444486377573, 0.6367618612362842,
    0.63439328416364549, 0.63201873593980906, 0.629638238914927,
    0.62725181549514408, 0.62485948814238634, 0.62246127937415,
    0.6200572117632891, 0.61764730793780387, 0.61523159058062682,
    0.61281008242940971, 0.61038280627630948, 0.60794978496777363,
    0.60551104140432555, 0.60306659854034816, 0.600616479383869,
    0.59816070699634238, 0.59569930449243336, 0.5932322950397998,
    0.59075970185887416, 0.58828154822264522, 0.58579785745643886,
    0.58330865293769829, 0.58081395809576453, 0.57831379641165559,
    0.57580819141784534, 0.5732971666980422, 0.57078074588696726,
    0.56825895267013149, 0.56573181078361312, 0.56319934401383409,
    0.560661576197336, 0.5581185312205561, 0.55557023301960218,
    0.55301670558002747, 0.55045797293660481, 0.54789405917310019,
    0.54532498842204646, 0.54275078486451589, 0.54017147272989285,
    0.53758707629564539, 0.53499761988709715, 0.5324031278771979,
    0.52980362468629461, 0.52719913478190128, 0.524589682678469,
    0.52197529293715439, 0.51935599016558964, 0.51673179901764987,
    0.51410274419322166, 0.5114688504379703, 0.508830142543107,
    0.50618664534515523, 0.50353838372571758, 0.50088538261124071,
    0.49822766697278181, 0.49556526182577254, 0.49289819222978404,
    0.49022648328829116, 0.487550160148436, 0.48486924800079106,
    0.48218377207912272, 0.47949375766015295, 0.47679923006332209,
    0.47410021465054997, 0.47139673682599764, 0.46868882203582796,
    0.46597649576796618, 0.46325978355186015, 0.46053871095824,
    0.45781330359887717, 0.45508358712634384, 0.45234958723377089,
    0.44961132965460654, 0.44686884016237416, 0.4441221445704292,
    0.44137126873171667, 0.43861623853852766, 0.43585707992225547,
    0.43309381885315196, 0.43032648134008261, 0.42755509343028208,
    0.42477968120910881, 0.42200027079979968, 0.41921688836322391,
    0.41642956009763715, 0.4136383122384345, 0.41084317105790391,
    0.40804416286497869, 0.40524131400498986, 0.40243465085941843,
    0.39962419984564679, 0.39680998741671031, 0.3939920400610481,
    0.39117038430225387, 0.38834504669882625, 0.38551605384391885,
    0.38268343236508978, 0.37984720892405116, 0.37700741021641826,
    0.37416406297145793, 0.37131719395183749, 0.36846682995337232,
    0.36561299780477385, 0.36275572436739723, 0.35989503653498811,
    0.35703096123343, 0.35416352542049034, 0.35129275608556709,
    0.34841868024943456, 0.34554132496398909, 0.34266071731199438,
    0.33977688440682685, 0.33688985339222005, 0.33399965144200938,
    0.33110630575987643, 0.3282098435790925, 0.32531029216226293,
    0.32240767880106985, 0.31950203081601569, 0.31659337555616585,
    0.31368174039889152, 0.31076715274961147, 0.30784964004153487,
    0.30492922973540237, 0.30200594931922808, 0.29907982630804048,
    0.29615088824362379, 0.29321916269425863, 0.29028467725446233,
    0.28734745954472951, 0.28440753721127188, 0.28146493792575794,
    0.27851968938505306, 0.27557181931095814, 0.272621355449949,
    0.26966832557291509, 0.26671275747489837, 0.26375467897483135,
    0.26079411791527551, 0.257831102162159, 0.25486565960451457,
    0.25189781815421697, 0.24892760574572015, 0.24595505033579459,
    0.24298017990326387, 0.2400030224487415, 0.2370236059943672,
    0.23404195858354343, 0.23105810828067111, 0.22807208317088573,
    0.22508391135979283, 0.22209362097320351, 0.2191012401568698,
    0.21610679707621952, 0.21311031991609136, 0.21011183688046961,
    0.20711137619221856, 0.20410896609281687, 0.2011046348420919,
    0.19809841071795356, 0.19509032201612825, 0.19208039704989244,
    0.18906866414980619, 0.18605515166344663, 0.18303988795514095,
    0.18002290140569951, 0.17700422041214875, 0.17398387338746382,
    0.17096188876030122, 0.16793829497473117, 0.16491312048996992,
    0.16188639378011183, 0.15885814333386145, 0.15582839765426523,
    0.15279718525844344, 0.14976453467732151, 0.14673047445536175,
    0.14369503315029447, 0.14065823933284921, 0.13762012158648604,
    0.13458070850712617, 0.13154002870288312, 0.12849811079379317,
    0.12545498341154623, 0.1224106751992162, 0.11936521481099135,
    0.11631863091190475, 0.11327095217756435, 0.11022220729388306,
    0.10717242495680884, 0.10412163387205459, 0.10106986275482782,
    0.0980171403295606, 0.094963495329638992, 0.091908956497132724,
    0.0888535525825246, 0.0857973123444399, 0.082740264549375692,
    0.079682437971430126, 0.076623861392031492, 0.073564563599667426,
    0.070504573389613856, 0.067443919563664051, 0.064382630929857465,
    0.061320736302208578, 0.058258264500435752, 0.055195244349689941,
    0.052131704680283324, 0.049067674327418015, 0.046003182130914623,
    0.04293825693494082, 0.039872927587739811, 0.036807222941358832,
    0.03374117185137758, 0.030674803176636626, 0.02760814577896574,
    0.024541228522912288, 0.021474080275469508, 0.01840672990580482,
    0.0153392062849881, 0.012271538285719925, 0.00920375478205982,
    0.0061358846491544753, 0.0030679567629659761, 0.0, -0.0030679567629659761,
    -0.0061358846491544753, -0.00920375478205982, -0.012271538285719925,
    -0.0153392062849881, -0.01840672990580482, -0.021474080275469508,
    -0.024541228522912288, -0.02760814577896574, -0.030674803176636626,
    -0.03374117185137758, -0.036807222941358832, -0.039872927587739811,
    -0.04293825693494082, -0.046003182130914623, -0.049067674327418015,
    -0.052131704680283324, -0.055195244349689941, -0.058258264500435752,
    -0.061320736302208578, -0.064382630929857465, -0.067443919563664051,
    -0.070504573389613856, -0.073564563599667426, -0.076623861392031492,
    -0.079682437971430126, -0.082740264549375692, -0.0857973123444399,
    -0.0888535525825246, -0.091908956497132724, -0.094963495329638992,
    -0.0980171403295606, -0.10106986275482782, -0.10412163387205459,
    -0.10717242495680884, -0.11022220729388306, -0.11327095217756435,
    -0.11631863091190475, -0.11936521481099135, -0.1224106751992162,
    -0.12545498341154623, -0.12849811079379317, -0.13154002870288312,
    -0.13458070850712617, -0.13762012158648604, -0.14065823933284921,
    -0.14369503315029447, -0.14673047445536175, -0.14976453467732151,
    -0.15279718525844344, -0.15582839765426523, -0.15885814333386145,
    -0.16188639378011183, -0.16491312048996992, -0.16793829497473117,
    -0.17096188876030122, -0.17398387338746382, -0.17700422041214875,
    -0.18002290140569951, -0.18303988795514095, -0.18605515166344663,
    -0.18906866414980619, -0.19208039704989244, -0.19509032201612825,
    -0.19809841071795356, -0.2011046348420919, -0.20410896609281687,
    -0.20711137619221856, -0.21011183688046961, -0.21311031991609136,
    -0.21610679707621952, -0.2191012401568698, -0.22209362097320351,
    -0.22508391135979283, -0.22807208317088573, -0.23105810828067111,
    -0.23404195858354343, -0.2370236059943672, -0.2400030224487415,
    -0.24298017990326387, -0.24595505033579459, -0.24892760574572015,
    -0.25189781815421697, -0.25486565960451457, -0.257831102162159,
    -0.26079411791527551, -0.26375467897483135, -0.26671275747489837,
    -0.26966832557291509, -0.272621355449949, -0.27557181931095814,
    -0.27851968938505306, -0.28146493792575794, -0.28440753721127188,
    -0.28734745954472951, -0.29028467725446233, -0.29321916269425863,
    -0.29615088824362379, -0.29907982630804048, -0.30200594931922808,
    -0.30492922973540237, -0.30784964004153487, -0.31076715274961147,
    -0.31368174039889152, -0.31659337555616585, -0.31950203081601569,
    -0.32240767880106985, -0.32531029216226293, -0.3282098435790925,
    -0.33110630575987643, -0.33399965144200938, -0.33688985339222005,
    -0.33977688440682685, -0.34266071731199438, -0.34554132496398909,
    -0.34841868024943456, -0.35129275608556709, -0.35416352542049034,
    -0.35703096123343, -0.35989503653498811, -0.36275572436739723,
    -0.36561299780477385, -0.36846682995337232, -0.37131719395183749,
    -0.37416406297145793, -0.37700741021641826, -0.37984720892405116,
    -0.38268343236508978, -0.38551605384391885, -0.38834504669882625,
    -0.39117038430225387, -0.3939920400610481, -0.39680998741671031,
    -0.39962419984564679, -0.40243465085941843, -0.40524131400498986,
    -0.40804416286497869, -0.41084317105790391, -0.4136383122384345,
    -0.41642956009763715, -0.41921688836322391, -0.42200027079979968,
    -0.42477968120910881, -0.42755509343028208, -0.43032648134008261,
    -0.43309381885315196, -0.43585707992225547, -0.43861623853852766,
    -0.44137126873171667, -0.4441221445704292, -0.44686884016237416,
    -0.44961132965460654, -0.45234958723377089, -0.45508358712634384,
    -0.45781330359887717, -0.46053871095824, -0.46325978355186015,
    -0.46597649576796618, -0.46868882203582796, -0.47139673682599764,
    -0.47410021465054997, -0.47679923006332209, -0.47949375766015295,
    -0.48218377207912272, -0.48486924800079106, -0.487550160148436,
    -0.49022648328829116, -0.49289819222978404, -0.49556526182577254,
    -0.49822766697278181, -0.50088538261124071, -0.50353838372571758,
    -0.50618664534515523, -0.508830142543107, -0.5114688504379703,
    -0.51410274419322166, -0.51673179901764987, -0.51935599016558964,
    -0.52197529293715439, -0.524589682678469, -0.52719913478190128,
    -0.52980362468629461, -0.5324031278771979, -0.53499761988709715,
    -0.53758707629564539, -0.54017147272989285, -0.54275078486451589,
    -0.54532498842204646, -0.54789405917310019, -0.55045797293660481,
    -0.55301670558002747, -0.55557023301960218, -0.5581185312205561,
    -0.560661576197336, -0.56319934401383409, -0.56573181078361312,
    -0.56825895267013149, -0.57078074588696726, -0.5732971666980422,
    -0.57580819141784534, -0.57831379641165559, -0.58081395809576453,
    -0.58330865293769829, -0.58579785745643886, -0.58828154822264522,
    -0.59075970185887416, -0.5932322950397998, -0.59569930449243336,
    -0.59816070699634238, -0.600616479383869, -0.60306659854034816,
    -0.60551104140432555, -0.60794978496777363, -0.61038280627630948,
    -0.61281008242940971, -0.61523159058062682, -0.61764730793780387,
    -0.6200572117632891, -0.62246127937415, -0.62485948814238634,
    -0.62725181549514408, -0.629638238914927, -0.63201873593980906,
    -0.63439328416364549, -0.6367618612362842, -0.63912444486377573,
    -0.641481012808583, -0.64383154288979139, -0.64617601298331628,
    -0.64851440102211244, -0.650846684996381, -0.65317284295377676,
    -0.65549285299961535, -0.65780669329707864, -0.66011434206742048,
    -0.66241577759017178, -0.66471097820334479, -0.66699992230363747,
    -0.669282588346636, -0.67155895484701833, -0.673829000378756,
    -0.67609270357531592, -0.67835004312986147, -0.680600997795453,
    -0.68284554638524808, -0.68508366777270036, -0.687315340891759,
    -0.68954054473706683, -0.69175925836415775, -0.69397146088965389,
    -0.696177131491463, -0.69837624940897292, -0.70056879394324834,
    -0.7027547444572253, -0.70493408037590488, -0.70710678118654757,
    -0.70927282643886569, -0.71143219574521643, -0.71358486878079352,
    -0.71573082528381859, -0.71787004505573171, -0.72000250796138165,
    -0.72212819392921535, -0.724247082951467, -0.726359155084346,
    -0.7284643904482252, -0.73056276922782759, -0.73265427167241282,
    -0.7347388780959635, -0.73681656887736979, -0.73888732446061511,
    -0.74095112535495922, -0.74300795213512172, -0.745057785441466,
    -0.74710060598018013, -0.74913639452345937, -0.75116513190968637,
    -0.75318679904361252, -0.75520137689653655, -0.75720884650648457,
    -0.759209188978388, -0.76120238548426178, -0.76318841726338127,
    -0.765167265622459, -0.7671389119358204, -0.7691033376455797,
    -0.77106052426181382, -0.773010453362737, -0.77495310659487393,
    -0.77688846567323244, -0.778816512381476, -0.78073722857209449,
    -0.78265059616657573, -0.78455659715557524, -0.78645521359908577,
    -0.78834642762660634, -0.79023022143731, -0.79210657730021239,
    -0.79397547755433717, -0.79583690460888357, -0.79769084094339116,
    -0.799537269107905, -0.80137617172314024, -0.80320753148064494,
    -0.80503133114296366, -0.80684755354379933, -0.808656181588175,
    -0.81045719825259477, -0.81225058658520388, -0.81403632970594841,
    -0.81581441080673378, -0.81758481315158371, -0.819347520076797,
    -0.82110251499110465, -0.82284978137582643, -0.82458930278502529,
    -0.82632106284566353, -0.8280450452577558, -0.829761233794523,
    -0.83146961230254524, -0.83317016470191319, -0.83486287498638,
    -0.836547727223512, -0.83822470555483808, -0.83989379419599952,
    -0.84155497743689844, -0.84320823964184544, -0.84485356524970712,
    -0.84649093877405213, -0.84812034480329723, -0.84974176800085255,
    -0.8513551931052652, -0.85296060493036363, -0.85455798836540053,
    -0.85614732837519447, -0.85772861000027212, -0.85930181835700847,
    -0.86086693863776731, -0.8624239561110405, -0.8639728561215867,
    -0.86551362409056909, -0.86704624551569265, -0.8685707059713409,
    -0.87008699110871146, -0.87159508665595109, -0.87309497841829009,
    -0.87458665227817611, -0.8760700941954066, -0.87754529020726135,
    -0.87901222642863353, -0.88047088905216075, -0.881921264348355,
    -0.88336333866573158, -0.88479709843093779, -0.88622253014888064,
    -0.88763962040285393, -0.88904835585466457, -0.89044872324475788,
    -0.89184070939234272, -0.89322430119551532, -0.8945994856313827,
    -0.89596624975618522, -0.89732458070541832, -0.89867446569395382,
    -0.90001589201616017, -0.901348847046022, -0.90267331823725883,
    -0.90398929312344334, -0.90529675931811882, -0.90659570451491533,
    -0.90788611648766626, -0.90916798309052238, -0.91044129225806725,
    -0.91170603200542988, -0.91296219042839821, -0.91420975570353069,
    -0.91544871608826783, -0.9166790599210427, -0.9179007756213905,
    -0.91911385169005777, -0.92031827670911059, -0.9215140393420419,
    -0.92270112833387863, -0.92387953251128674, -0.92504924078267758,
    -0.92621024213831138, -0.92736252565040111, -0.92850608047321559,
    -0.92964089584318121, -0.93076696107898371, -0.93188426558166815,
    -0.932992798834739, -0.93409255040425887, -0.93518350993894761,
    -0.93626566717027826, -0.937339011912575, -0.93840353406310806,
    -0.93945922360218992, -0.9405060705932683, -0.94154406518302081,
    -0.94257319760144687, -0.94359345816196039, -0.94460483726148026,
    -0.94560732538052128, -0.94660091308328353, -0.94758559101774109,
    -0.94856134991573027, -0.94952818059303667, -0.9504860739494817,
    -0.95143502096900834, -0.95237501271976588, -0.95330604035419386,
    -0.95422809510910567, -0.95514116830577078, -0.95604525134999641,
    -0.95694033573220882, -0.95782641302753291, -0.9587034748958716,
    -0.95957151308198452, -0.96043051941556579, -0.96128048581132064,
    -0.96212140426904158, -0.96295326687368388, -0.96377606579543984,
    -0.96458979328981276, -0.9653944416976894, -0.9661900034454125,
    -0.96697647104485207, -0.96775383709347551, -0.96852209427441727,
    -0.96928123535654853, -0.970031253194544, -0.97077214072895035,
    -0.97150389098625178, -0.97222649707893627, -0.97293995220556018,
    -0.973644249650812, -0.97433938278557586, -0.97502534506699412,
    -0.97570213003852857, -0.97636973133002114, -0.97702814265775439,
    -0.97767735782450993, -0.97831737071962765, -0.9789481753190622,
    -0.97956976568544052, -0.98018213596811743, -0.98078528040323043,
    -0.98137919331375456, -0.98196386910955524, -0.98253930228744124,
    -0.98310548743121629, -0.98366241921173025, -0.984210092386929,
    -0.98474850180190421, -0.98527764238894122, -0.98579750916756748,
    -0.98630809724459867, -0.98680940181418553, -0.98730141815785843,
    -0.98778414164457218, -0.98825756773074946, -0.98872169196032378,
    -0.989176509964781, -0.98962201746320089, -0.99005821026229712,
    -0.99048508425645709, -0.99090263542778, -0.99131085984611544,
    -0.99170975366909953, -0.9920993131421918, -0.99247953459871,
    -0.9928504144598651, -0.9932119492347945, -0.9935641355205953,
    -0.99390697000235606, -0.9942404494531879, -0.99456457073425542,
    -0.99487933079480562, -0.99518472667219693, -0.99548075549192694,
    -0.99576741446765982, -0.996044700901252, -0.996312612182778,
    -0.99657114579055484, -0.99682029929116567, -0.997060070339483,
    -0.99729045667869021, -0.99751145614030345, -0.99772306664419164,
    -0.997925286198596, -0.99811811290014918, -0.99830154493389289,
    -0.99847558057329477, -0.99864021818026527, -0.99879545620517241,
    -0.99894129318685687, -0.99907772775264536, -0.99920475861836389,
    -0.99932238458834954, -0.99943060455546173, -0.99952941750109314,
    -0.99961882249517864, -0.99969881869620425, -0.99976940535121528,
    -0.9998305817958234, -0.99988234745421256, -0.9999247018391445,
    -0.9999576445519639, -0.99998117528260111, -0.99999529380957619, -1.0 };

  double twid_im;
  static const double dv21[1025] = { 0.0, -0.0030679567629659761,
    -0.0061358846491544753, -0.00920375478205982, -0.012271538285719925,
    -0.0153392062849881, -0.01840672990580482, -0.021474080275469508,
    -0.024541228522912288, -0.02760814577896574, -0.030674803176636626,
    -0.03374117185137758, -0.036807222941358832, -0.039872927587739811,
    -0.04293825693494082, -0.046003182130914623, -0.049067674327418015,
    -0.052131704680283324, -0.055195244349689941, -0.058258264500435752,
    -0.061320736302208578, -0.064382630929857465, -0.067443919563664051,
    -0.070504573389613856, -0.073564563599667426, -0.076623861392031492,
    -0.079682437971430126, -0.082740264549375692, -0.0857973123444399,
    -0.0888535525825246, -0.091908956497132724, -0.094963495329638992,
    -0.0980171403295606, -0.10106986275482782, -0.10412163387205459,
    -0.10717242495680884, -0.11022220729388306, -0.11327095217756435,
    -0.11631863091190475, -0.11936521481099135, -0.1224106751992162,
    -0.12545498341154623, -0.12849811079379317, -0.13154002870288312,
    -0.13458070850712617, -0.13762012158648604, -0.14065823933284921,
    -0.14369503315029447, -0.14673047445536175, -0.14976453467732151,
    -0.15279718525844344, -0.15582839765426523, -0.15885814333386145,
    -0.16188639378011183, -0.16491312048996992, -0.16793829497473117,
    -0.17096188876030122, -0.17398387338746382, -0.17700422041214875,
    -0.18002290140569951, -0.18303988795514095, -0.18605515166344663,
    -0.18906866414980619, -0.19208039704989244, -0.19509032201612825,
    -0.19809841071795356, -0.2011046348420919, -0.20410896609281687,
    -0.20711137619221856, -0.21011183688046961, -0.21311031991609136,
    -0.21610679707621952, -0.2191012401568698, -0.22209362097320351,
    -0.22508391135979283, -0.22807208317088573, -0.23105810828067111,
    -0.23404195858354343, -0.2370236059943672, -0.2400030224487415,
    -0.24298017990326387, -0.24595505033579459, -0.24892760574572015,
    -0.25189781815421697, -0.25486565960451457, -0.257831102162159,
    -0.26079411791527551, -0.26375467897483135, -0.26671275747489837,
    -0.26966832557291509, -0.272621355449949, -0.27557181931095814,
    -0.27851968938505306, -0.28146493792575794, -0.28440753721127188,
    -0.28734745954472951, -0.29028467725446233, -0.29321916269425863,
    -0.29615088824362379, -0.29907982630804048, -0.30200594931922808,
    -0.30492922973540237, -0.30784964004153487, -0.31076715274961147,
    -0.31368174039889152, -0.31659337555616585, -0.31950203081601569,
    -0.32240767880106985, -0.32531029216226293, -0.3282098435790925,
    -0.33110630575987643, -0.33399965144200938, -0.33688985339222005,
    -0.33977688440682685, -0.34266071731199438, -0.34554132496398909,
    -0.34841868024943456, -0.35129275608556709, -0.35416352542049034,
    -0.35703096123343, -0.35989503653498811, -0.36275572436739723,
    -0.36561299780477385, -0.36846682995337232, -0.37131719395183749,
    -0.37416406297145793, -0.37700741021641826, -0.37984720892405116,
    -0.38268343236508978, -0.38551605384391885, -0.38834504669882625,
    -0.39117038430225387, -0.3939920400610481, -0.39680998741671031,
    -0.39962419984564679, -0.40243465085941843, -0.40524131400498986,
    -0.40804416286497869, -0.41084317105790391, -0.4136383122384345,
    -0.41642956009763715, -0.41921688836322391, -0.42200027079979968,
    -0.42477968120910881, -0.42755509343028208, -0.43032648134008261,
    -0.43309381885315196, -0.43585707992225547, -0.43861623853852766,
    -0.44137126873171667, -0.4441221445704292, -0.44686884016237416,
    -0.44961132965460654, -0.45234958723377089, -0.45508358712634384,
    -0.45781330359887717, -0.46053871095824, -0.46325978355186015,
    -0.46597649576796618, -0.46868882203582796, -0.47139673682599764,
    -0.47410021465054997, -0.47679923006332209, -0.47949375766015295,
    -0.48218377207912272, -0.48486924800079106, -0.487550160148436,
    -0.49022648328829116, -0.49289819222978404, -0.49556526182577254,
    -0.49822766697278181, -0.50088538261124071, -0.50353838372571758,
    -0.50618664534515523, -0.508830142543107, -0.5114688504379703,
    -0.51410274419322166, -0.51673179901764987, -0.51935599016558964,
    -0.52197529293715439, -0.524589682678469, -0.52719913478190128,
    -0.52980362468629461, -0.5324031278771979, -0.53499761988709715,
    -0.53758707629564539, -0.54017147272989285, -0.54275078486451589,
    -0.54532498842204646, -0.54789405917310019, -0.55045797293660481,
    -0.55301670558002747, -0.55557023301960218, -0.5581185312205561,
    -0.560661576197336, -0.56319934401383409, -0.56573181078361312,
    -0.56825895267013149, -0.57078074588696726, -0.5732971666980422,
    -0.57580819141784534, -0.57831379641165559, -0.58081395809576453,
    -0.58330865293769829, -0.58579785745643886, -0.58828154822264522,
    -0.59075970185887416, -0.5932322950397998, -0.59569930449243336,
    -0.59816070699634238, -0.600616479383869, -0.60306659854034816,
    -0.60551104140432555, -0.60794978496777363, -0.61038280627630948,
    -0.61281008242940971, -0.61523159058062682, -0.61764730793780387,
    -0.6200572117632891, -0.62246127937415, -0.62485948814238634,
    -0.62725181549514408, -0.629638238914927, -0.63201873593980906,
    -0.63439328416364549, -0.6367618612362842, -0.63912444486377573,
    -0.641481012808583, -0.64383154288979139, -0.64617601298331628,
    -0.64851440102211244, -0.650846684996381, -0.65317284295377676,
    -0.65549285299961535, -0.65780669329707864, -0.66011434206742048,
    -0.66241577759017178, -0.66471097820334479, -0.66699992230363747,
    -0.669282588346636, -0.67155895484701833, -0.673829000378756,
    -0.67609270357531592, -0.67835004312986147, -0.680600997795453,
    -0.68284554638524808, -0.68508366777270036, -0.687315340891759,
    -0.68954054473706683, -0.69175925836415775, -0.69397146088965389,
    -0.696177131491463, -0.69837624940897292, -0.70056879394324834,
    -0.7027547444572253, -0.70493408037590488, -0.70710678118654757,
    -0.70927282643886569, -0.71143219574521643, -0.71358486878079352,
    -0.71573082528381859, -0.71787004505573171, -0.72000250796138165,
    -0.72212819392921535, -0.724247082951467, -0.726359155084346,
    -0.7284643904482252, -0.73056276922782759, -0.73265427167241282,
    -0.7347388780959635, -0.73681656887736979, -0.73888732446061511,
    -0.74095112535495922, -0.74300795213512172, -0.745057785441466,
    -0.74710060598018013, -0.74913639452345937, -0.75116513190968637,
    -0.75318679904361252, -0.75520137689653655, -0.75720884650648457,
    -0.759209188978388, -0.76120238548426178, -0.76318841726338127,
    -0.765167265622459, -0.7671389119358204, -0.7691033376455797,
    -0.77106052426181382, -0.773010453362737, -0.77495310659487393,
    -0.77688846567323244, -0.778816512381476, -0.78073722857209449,
    -0.78265059616657573, -0.78455659715557524, -0.78645521359908577,
    -0.78834642762660634, -0.79023022143731, -0.79210657730021239,
    -0.79397547755433717, -0.79583690460888357, -0.79769084094339116,
    -0.799537269107905, -0.80137617172314024, -0.80320753148064494,
    -0.80503133114296366, -0.80684755354379933, -0.808656181588175,
    -0.81045719825259477, -0.81225058658520388, -0.81403632970594841,
    -0.81581441080673378, -0.81758481315158371, -0.819347520076797,
    -0.82110251499110465, -0.82284978137582643, -0.82458930278502529,
    -0.82632106284566353, -0.8280450452577558, -0.829761233794523,
    -0.83146961230254524, -0.83317016470191319, -0.83486287498638,
    -0.836547727223512, -0.83822470555483808, -0.83989379419599952,
    -0.84155497743689844, -0.84320823964184544, -0.84485356524970712,
    -0.84649093877405213, -0.84812034480329723, -0.84974176800085255,
    -0.8513551931052652, -0.85296060493036363, -0.85455798836540053,
    -0.85614732837519447, -0.85772861000027212, -0.85930181835700847,
    -0.86086693863776731, -0.8624239561110405, -0.8639728561215867,
    -0.86551362409056909, -0.86704624551569265, -0.8685707059713409,
    -0.87008699110871146, -0.87159508665595109, -0.87309497841829009,
    -0.87458665227817611, -0.8760700941954066, -0.87754529020726135,
    -0.87901222642863353, -0.88047088905216075, -0.881921264348355,
    -0.88336333866573158, -0.88479709843093779, -0.88622253014888064,
    -0.88763962040285393, -0.88904835585466457, -0.89044872324475788,
    -0.89184070939234272, -0.89322430119551532, -0.8945994856313827,
    -0.89596624975618522, -0.89732458070541832, -0.89867446569395382,
    -0.90001589201616017, -0.901348847046022, -0.90267331823725883,
    -0.90398929312344334, -0.90529675931811882, -0.90659570451491533,
    -0.90788611648766626, -0.90916798309052238, -0.91044129225806725,
    -0.91170603200542988, -0.91296219042839821, -0.91420975570353069,
    -0.91544871608826783, -0.9166790599210427, -0.9179007756213905,
    -0.91911385169005777, -0.92031827670911059, -0.9215140393420419,
    -0.92270112833387863, -0.92387953251128674, -0.92504924078267758,
    -0.92621024213831138, -0.92736252565040111, -0.92850608047321559,
    -0.92964089584318121, -0.93076696107898371, -0.93188426558166815,
    -0.932992798834739, -0.93409255040425887, -0.93518350993894761,
    -0.93626566717027826, -0.937339011912575, -0.93840353406310806,
    -0.93945922360218992, -0.9405060705932683, -0.94154406518302081,
    -0.94257319760144687, -0.94359345816196039, -0.94460483726148026,
    -0.94560732538052128, -0.94660091308328353, -0.94758559101774109,
    -0.94856134991573027, -0.94952818059303667, -0.9504860739494817,
    -0.95143502096900834, -0.95237501271976588, -0.95330604035419386,
    -0.95422809510910567, -0.95514116830577078, -0.95604525134999641,
    -0.95694033573220882, -0.95782641302753291, -0.9587034748958716,
    -0.95957151308198452, -0.96043051941556579, -0.96128048581132064,
    -0.96212140426904158, -0.96295326687368388, -0.96377606579543984,
    -0.96458979328981276, -0.9653944416976894, -0.9661900034454125,
    -0.96697647104485207, -0.96775383709347551, -0.96852209427441727,
    -0.96928123535654853, -0.970031253194544, -0.97077214072895035,
    -0.97150389098625178, -0.97222649707893627, -0.97293995220556018,
    -0.973644249650812, -0.97433938278557586, -0.97502534506699412,
    -0.97570213003852857, -0.97636973133002114, -0.97702814265775439,
    -0.97767735782450993, -0.97831737071962765, -0.9789481753190622,
    -0.97956976568544052, -0.98018213596811743, -0.98078528040323043,
    -0.98137919331375456, -0.98196386910955524, -0.98253930228744124,
    -0.98310548743121629, -0.98366241921173025, -0.984210092386929,
    -0.98474850180190421, -0.98527764238894122, -0.98579750916756748,
    -0.98630809724459867, -0.98680940181418553, -0.98730141815785843,
    -0.98778414164457218, -0.98825756773074946, -0.98872169196032378,
    -0.989176509964781, -0.98962201746320089, -0.99005821026229712,
    -0.99048508425645709, -0.99090263542778, -0.99131085984611544,
    -0.99170975366909953, -0.9920993131421918, -0.99247953459871,
    -0.9928504144598651, -0.9932119492347945, -0.9935641355205953,
    -0.99390697000235606, -0.9942404494531879, -0.99456457073425542,
    -0.99487933079480562, -0.99518472667219693, -0.99548075549192694,
    -0.99576741446765982, -0.996044700901252, -0.996312612182778,
    -0.99657114579055484, -0.99682029929116567, -0.997060070339483,
    -0.99729045667869021, -0.99751145614030345, -0.99772306664419164,
    -0.997925286198596, -0.99811811290014918, -0.99830154493389289,
    -0.99847558057329477, -0.99864021818026527, -0.99879545620517241,
    -0.99894129318685687, -0.99907772775264536, -0.99920475861836389,
    -0.99932238458834954, -0.99943060455546173, -0.99952941750109314,
    -0.99961882249517864, -0.99969881869620425, -0.99976940535121528,
    -0.9998305817958234, -0.99988234745421256, -0.9999247018391445,
    -0.9999576445519639, -0.99998117528260111, -0.99999529380957619, -1.0,
    -0.99999529380957619, -0.99998117528260111, -0.9999576445519639,
    -0.9999247018391445, -0.99988234745421256, -0.9998305817958234,
    -0.99976940535121528, -0.99969881869620425, -0.99961882249517864,
    -0.99952941750109314, -0.99943060455546173, -0.99932238458834954,
    -0.99920475861836389, -0.99907772775264536, -0.99894129318685687,
    -0.99879545620517241, -0.99864021818026527, -0.99847558057329477,
    -0.99830154493389289, -0.99811811290014918, -0.997925286198596,
    -0.99772306664419164, -0.99751145614030345, -0.99729045667869021,
    -0.997060070339483, -0.99682029929116567, -0.99657114579055484,
    -0.996312612182778, -0.996044700901252, -0.99576741446765982,
    -0.99548075549192694, -0.99518472667219693, -0.99487933079480562,
    -0.99456457073425542, -0.9942404494531879, -0.99390697000235606,
    -0.9935641355205953, -0.9932119492347945, -0.9928504144598651,
    -0.99247953459871, -0.9920993131421918, -0.99170975366909953,
    -0.99131085984611544, -0.99090263542778, -0.99048508425645709,
    -0.99005821026229712, -0.98962201746320089, -0.989176509964781,
    -0.98872169196032378, -0.98825756773074946, -0.98778414164457218,
    -0.98730141815785843, -0.98680940181418553, -0.98630809724459867,
    -0.98579750916756748, -0.98527764238894122, -0.98474850180190421,
    -0.984210092386929, -0.98366241921173025, -0.98310548743121629,
    -0.98253930228744124, -0.98196386910955524, -0.98137919331375456,
    -0.98078528040323043, -0.98018213596811743, -0.97956976568544052,
    -0.9789481753190622, -0.97831737071962765, -0.97767735782450993,
    -0.97702814265775439, -0.97636973133002114, -0.97570213003852857,
    -0.97502534506699412, -0.97433938278557586, -0.973644249650812,
    -0.97293995220556018, -0.97222649707893627, -0.97150389098625178,
    -0.97077214072895035, -0.970031253194544, -0.96928123535654853,
    -0.96852209427441727, -0.96775383709347551, -0.96697647104485207,
    -0.9661900034454125, -0.9653944416976894, -0.96458979328981276,
    -0.96377606579543984, -0.96295326687368388, -0.96212140426904158,
    -0.96128048581132064, -0.96043051941556579, -0.95957151308198452,
    -0.9587034748958716, -0.95782641302753291, -0.95694033573220882,
    -0.95604525134999641, -0.95514116830577078, -0.95422809510910567,
    -0.95330604035419386, -0.95237501271976588, -0.95143502096900834,
    -0.9504860739494817, -0.94952818059303667, -0.94856134991573027,
    -0.94758559101774109, -0.94660091308328353, -0.94560732538052128,
    -0.94460483726148026, -0.94359345816196039, -0.94257319760144687,
    -0.94154406518302081, -0.9405060705932683, -0.93945922360218992,
    -0.93840353406310806, -0.937339011912575, -0.93626566717027826,
    -0.93518350993894761, -0.93409255040425887, -0.932992798834739,
    -0.93188426558166815, -0.93076696107898371, -0.92964089584318121,
    -0.92850608047321559, -0.92736252565040111, -0.92621024213831138,
    -0.92504924078267758, -0.92387953251128674, -0.92270112833387863,
    -0.9215140393420419, -0.92031827670911059, -0.91911385169005777,
    -0.9179007756213905, -0.9166790599210427, -0.91544871608826783,
    -0.91420975570353069, -0.91296219042839821, -0.91170603200542988,
    -0.91044129225806725, -0.90916798309052238, -0.90788611648766626,
    -0.90659570451491533, -0.90529675931811882, -0.90398929312344334,
    -0.90267331823725883, -0.901348847046022, -0.90001589201616017,
    -0.89867446569395382, -0.89732458070541832, -0.89596624975618522,
    -0.8945994856313827, -0.89322430119551532, -0.89184070939234272,
    -0.89044872324475788, -0.88904835585466457, -0.88763962040285393,
    -0.88622253014888064, -0.88479709843093779, -0.88336333866573158,
    -0.881921264348355, -0.88047088905216075, -0.87901222642863353,
    -0.87754529020726135, -0.8760700941954066, -0.87458665227817611,
    -0.87309497841829009, -0.87159508665595109, -0.87008699110871146,
    -0.8685707059713409, -0.86704624551569265, -0.86551362409056909,
    -0.8639728561215867, -0.8624239561110405, -0.86086693863776731,
    -0.85930181835700847, -0.85772861000027212, -0.85614732837519447,
    -0.85455798836540053, -0.85296060493036363, -0.8513551931052652,
    -0.84974176800085255, -0.84812034480329723, -0.84649093877405213,
    -0.84485356524970712, -0.84320823964184544, -0.84155497743689844,
    -0.83989379419599952, -0.83822470555483808, -0.836547727223512,
    -0.83486287498638, -0.83317016470191319, -0.83146961230254524,
    -0.829761233794523, -0.8280450452577558, -0.82632106284566353,
    -0.82458930278502529, -0.82284978137582643, -0.82110251499110465,
    -0.819347520076797, -0.81758481315158371, -0.81581441080673378,
    -0.81403632970594841, -0.81225058658520388, -0.81045719825259477,
    -0.808656181588175, -0.80684755354379933, -0.80503133114296366,
    -0.80320753148064494, -0.80137617172314024, -0.799537269107905,
    -0.79769084094339116, -0.79583690460888357, -0.79397547755433717,
    -0.79210657730021239, -0.79023022143731, -0.78834642762660634,
    -0.78645521359908577, -0.78455659715557524, -0.78265059616657573,
    -0.78073722857209449, -0.778816512381476, -0.77688846567323244,
    -0.77495310659487393, -0.773010453362737, -0.77106052426181382,
    -0.7691033376455797, -0.7671389119358204, -0.765167265622459,
    -0.76318841726338127, -0.76120238548426178, -0.759209188978388,
    -0.75720884650648457, -0.75520137689653655, -0.75318679904361252,
    -0.75116513190968637, -0.74913639452345937, -0.74710060598018013,
    -0.745057785441466, -0.74300795213512172, -0.74095112535495922,
    -0.73888732446061511, -0.73681656887736979, -0.7347388780959635,
    -0.73265427167241282, -0.73056276922782759, -0.7284643904482252,
    -0.726359155084346, -0.724247082951467, -0.72212819392921535,
    -0.72000250796138165, -0.71787004505573171, -0.71573082528381859,
    -0.71358486878079352, -0.71143219574521643, -0.70927282643886569,
    -0.70710678118654757, -0.70493408037590488, -0.7027547444572253,
    -0.70056879394324834, -0.69837624940897292, -0.696177131491463,
    -0.69397146088965389, -0.69175925836415775, -0.68954054473706683,
    -0.687315340891759, -0.68508366777270036, -0.68284554638524808,
    -0.680600997795453, -0.67835004312986147, -0.67609270357531592,
    -0.673829000378756, -0.67155895484701833, -0.669282588346636,
    -0.66699992230363747, -0.66471097820334479, -0.66241577759017178,
    -0.66011434206742048, -0.65780669329707864, -0.65549285299961535,
    -0.65317284295377676, -0.650846684996381, -0.64851440102211244,
    -0.64617601298331628, -0.64383154288979139, -0.641481012808583,
    -0.63912444486377573, -0.6367618612362842, -0.63439328416364549,
    -0.63201873593980906, -0.629638238914927, -0.62725181549514408,
    -0.62485948814238634, -0.62246127937415, -0.6200572117632891,
    -0.61764730793780387, -0.61523159058062682, -0.61281008242940971,
    -0.61038280627630948, -0.60794978496777363, -0.60551104140432555,
    -0.60306659854034816, -0.600616479383869, -0.59816070699634238,
    -0.59569930449243336, -0.5932322950397998, -0.59075970185887416,
    -0.58828154822264522, -0.58579785745643886, -0.58330865293769829,
    -0.58081395809576453, -0.57831379641165559, -0.57580819141784534,
    -0.5732971666980422, -0.57078074588696726, -0.56825895267013149,
    -0.56573181078361312, -0.56319934401383409, -0.560661576197336,
    -0.5581185312205561, -0.55557023301960218, -0.55301670558002747,
    -0.55045797293660481, -0.54789405917310019, -0.54532498842204646,
    -0.54275078486451589, -0.54017147272989285, -0.53758707629564539,
    -0.53499761988709715, -0.5324031278771979, -0.52980362468629461,
    -0.52719913478190128, -0.524589682678469, -0.52197529293715439,
    -0.51935599016558964, -0.51673179901764987, -0.51410274419322166,
    -0.5114688504379703, -0.508830142543107, -0.50618664534515523,
    -0.50353838372571758, -0.50088538261124071, -0.49822766697278181,
    -0.49556526182577254, -0.49289819222978404, -0.49022648328829116,
    -0.487550160148436, -0.48486924800079106, -0.48218377207912272,
    -0.47949375766015295, -0.47679923006332209, -0.47410021465054997,
    -0.47139673682599764, -0.46868882203582796, -0.46597649576796618,
    -0.46325978355186015, -0.46053871095824, -0.45781330359887717,
    -0.45508358712634384, -0.45234958723377089, -0.44961132965460654,
    -0.44686884016237416, -0.4441221445704292, -0.44137126873171667,
    -0.43861623853852766, -0.43585707992225547, -0.43309381885315196,
    -0.43032648134008261, -0.42755509343028208, -0.42477968120910881,
    -0.42200027079979968, -0.41921688836322391, -0.41642956009763715,
    -0.4136383122384345, -0.41084317105790391, -0.40804416286497869,
    -0.40524131400498986, -0.40243465085941843, -0.39962419984564679,
    -0.39680998741671031, -0.3939920400610481, -0.39117038430225387,
    -0.38834504669882625, -0.38551605384391885, -0.38268343236508978,
    -0.37984720892405116, -0.37700741021641826, -0.37416406297145793,
    -0.37131719395183749, -0.36846682995337232, -0.36561299780477385,
    -0.36275572436739723, -0.35989503653498811, -0.35703096123343,
    -0.35416352542049034, -0.35129275608556709, -0.34841868024943456,
    -0.34554132496398909, -0.34266071731199438, -0.33977688440682685,
    -0.33688985339222005, -0.33399965144200938, -0.33110630575987643,
    -0.3282098435790925, -0.32531029216226293, -0.32240767880106985,
    -0.31950203081601569, -0.31659337555616585, -0.31368174039889152,
    -0.31076715274961147, -0.30784964004153487, -0.30492922973540237,
    -0.30200594931922808, -0.29907982630804048, -0.29615088824362379,
    -0.29321916269425863, -0.29028467725446233, -0.28734745954472951,
    -0.28440753721127188, -0.28146493792575794, -0.27851968938505306,
    -0.27557181931095814, -0.272621355449949, -0.26966832557291509,
    -0.26671275747489837, -0.26375467897483135, -0.26079411791527551,
    -0.257831102162159, -0.25486565960451457, -0.25189781815421697,
    -0.24892760574572015, -0.24595505033579459, -0.24298017990326387,
    -0.2400030224487415, -0.2370236059943672, -0.23404195858354343,
    -0.23105810828067111, -0.22807208317088573, -0.22508391135979283,
    -0.22209362097320351, -0.2191012401568698, -0.21610679707621952,
    -0.21311031991609136, -0.21011183688046961, -0.20711137619221856,
    -0.20410896609281687, -0.2011046348420919, -0.19809841071795356,
    -0.19509032201612825, -0.19208039704989244, -0.18906866414980619,
    -0.18605515166344663, -0.18303988795514095, -0.18002290140569951,
    -0.17700422041214875, -0.17398387338746382, -0.17096188876030122,
    -0.16793829497473117, -0.16491312048996992, -0.16188639378011183,
    -0.15885814333386145, -0.15582839765426523, -0.15279718525844344,
    -0.14976453467732151, -0.14673047445536175, -0.14369503315029447,
    -0.14065823933284921, -0.13762012158648604, -0.13458070850712617,
    -0.13154002870288312, -0.12849811079379317, -0.12545498341154623,
    -0.1224106751992162, -0.11936521481099135, -0.11631863091190475,
    -0.11327095217756435, -0.11022220729388306, -0.10717242495680884,
    -0.10412163387205459, -0.10106986275482782, -0.0980171403295606,
    -0.094963495329638992, -0.091908956497132724, -0.0888535525825246,
    -0.0857973123444399, -0.082740264549375692, -0.079682437971430126,
    -0.076623861392031492, -0.073564563599667426, -0.070504573389613856,
    -0.067443919563664051, -0.064382630929857465, -0.061320736302208578,
    -0.058258264500435752, -0.055195244349689941, -0.052131704680283324,
    -0.049067674327418015, -0.046003182130914623, -0.04293825693494082,
    -0.039872927587739811, -0.036807222941358832, -0.03374117185137758,
    -0.030674803176636626, -0.02760814577896574, -0.024541228522912288,
    -0.021474080275469508, -0.01840672990580482, -0.0153392062849881,
    -0.012271538285719925, -0.00920375478205982, -0.0061358846491544753,
    -0.0030679567629659761, -0.0 };

  int ihi;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //               Short-Time Fourier Transform            %
  //                with MATLAB Implementation             %
  //                                                       %
  //  Author: M.Sc. Eng. Hristo Zhivomirov       12/21/13  %
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //  represent x as column-vector if it is not x = x(:);
  //  if size(x, 2) > 1
  //      x = x';
  //  end
  //  length of the signal
  //  wlen = 2^nextpow2(fs);
  //  form a periodic hamming window (!REPLACED EXTERNALLY)
  //  win = hamming(256, 'periodic');
  // WINDOW LENGTH:
  //  form the stft matrix
  //  calculate the total number of rows
  nm1d2 = x->size[1];
  twid_re = ((double)nm1d2 - 256.0) / 64.0;
  if (twid_re < 0.0) {
    b_x = (int)std::ceil(twid_re);
  } else {
    b_x = (int)std::floor(twid_re);
  }

  //  calculate the total number of columns
  nm1d2 = s->size[0] * s->size[1];
  s->size[0] = 1025;
  s->size[1] = 1 + b_x;
  emxEnsureCapacity((emxArray__common *)s, nm1d2, (int)sizeof(creal_T));
  cdiff = 1025 * (1 + b_x);
  for (nm1d2 = 0; nm1d2 < cdiff; nm1d2++) {
    s->data[nm1d2].re = 0.0;
    s->data[nm1d2].im = 0.0;
  }

  //  form the stft matrix
  //  for r = 1:rown
  //      for c = 1:coln
  //          s(r,c) = complex(0,0);
  //      end
  //  end
  //  initialize the indexes
  indx = 0.0;
  col = 0;

  //  perform STFT
  do {
    exitg1 = 0;
    nm1d2 = x->size[1];
    if (indx + 256.0 <= nm1d2) {
      //  windowing
      for (nm1d2 = 0; nm1d2 < 256; nm1d2++) {
        xw[nm1d2] = x->data[(int)(indx + (1.0 + (double)nm1d2)) - 1] * win[nm1d2];
      }

      //  FFT
      for (i = 0; i < 2048; i++) {
        X[i].re = 0.0;
        X[i].im = 0.0;
      }

      nm1d2 = 0;
      ndbl = 0;
      cdiff = 0;
      for (i = 0; i < 255; i++) {
        X[cdiff].re = xw[nm1d2];
        X[cdiff].im = 0.0;
        cdiff = 2048;
        tst = true;
        while (tst) {
          cdiff >>= 1;
          ndbl ^= cdiff;
          tst = ((ndbl & cdiff) == 0);
        }

        cdiff = ndbl;
        nm1d2++;
      }

      X[cdiff].re = xw[nm1d2];
      X[cdiff].im = 0.0;
      for (i = 0; i <= 2047; i += 2) {
        temp_re = X[i + 1].re;
        temp_im = X[i + 1].im;
        X[i + 1].re = X[i].re - X[i + 1].re;
        X[i + 1].im = X[i].im - X[i + 1].im;
        X[i].re += temp_re;
        X[i].im += temp_im;
      }

      cdiff = 2;
      nm1d2 = 4;
      k = 512;
      ndbl = 2045;
      while (k > 0) {
        for (i = 0; i < ndbl; i += nm1d2) {
          temp_re = X[i + cdiff].re;
          temp_im = X[i + cdiff].im;
          X[i + cdiff].re = X[i].re - temp_re;
          X[i + cdiff].im = X[i].im - temp_im;
          X[i].re += temp_re;
          X[i].im += temp_im;
        }

        apnd = 1;
        for (j = k; j < 1024; j += k) {
          twid_re = dv20[j];
          twid_im = dv21[j];
          i = apnd;
          ihi = apnd + ndbl;
          while (i < ihi) {
            temp_re = twid_re * X[i + cdiff].re - twid_im * X[i + cdiff].im;
            temp_im = twid_re * X[i + cdiff].im + twid_im * X[i + cdiff].re;
            X[i + cdiff].re = X[i].re - temp_re;
            X[i + cdiff].im = X[i].im - temp_im;
            X[i].re += temp_re;
            X[i].im += temp_im;
            i += nm1d2;
          }

          apnd++;
        }

        k /= 2;
        cdiff = nm1d2;
        nm1d2 <<= 1;
        ndbl -= cdiff;
      }

      //  update the stft matrix
      //      s(:, col) = X(1:rown); %%%% OLD
      for (nm1d2 = 0; nm1d2 < 1025; nm1d2++) {
        s->data[nm1d2 + s->size[0] * col] = X[nm1d2];
      }

      //  update the indexes
      indx += 64.0;
      col++;
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  //  calculate the time and frequency vectors
  nm1d2 = 128 + (b_x << 6);
  if (nm1d2 < 128) {
    ndbl = 0;
    apnd = nm1d2;
  } else {
    ndbl = (int)std::floor(((double)nm1d2 - 128.0) / 64.0 + 0.5);
    apnd = 128 + (ndbl << 6);
    cdiff = apnd - nm1d2;
    if (fabs((double)cdiff) < 4.4408920985006262E-16 * (double)nm1d2) {
      ndbl++;
      apnd = nm1d2;
    } else if (cdiff > 0) {
      apnd = 128 + ((ndbl - 1) << 6);
    } else {
      ndbl++;
    }
  }

  nm1d2 = t->size[0] * t->size[1];
  t->size[0] = 1;
  t->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)t, nm1d2, (int)sizeof(double));
  if (ndbl > 0) {
    t->data[0] = 128.0;
    if (ndbl > 1) {
      t->data[ndbl - 1] = apnd;
      nm1d2 = (ndbl - 1) / 2;
      for (k = 1; k < nm1d2; k++) {
        cdiff = k << 6;
        t->data[k] = 128.0 + (double)cdiff;
        t->data[(ndbl - k) - 1] = apnd - cdiff;
      }

      if (nm1d2 << 1 == ndbl - 1) {
        t->data[nm1d2] = (128.0 + (double)apnd) / 2.0;
      } else {
        cdiff = nm1d2 << 6;
        t->data[nm1d2] = 128.0 + (double)cdiff;
        t->data[nm1d2 + 1] = apnd - cdiff;
      }
    }
  }

  nm1d2 = t->size[0] * t->size[1];
  t->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)t, nm1d2, (int)sizeof(double));
  nm1d2 = t->size[0];
  cdiff = t->size[1];
  cdiff *= nm1d2;
  for (nm1d2 = 0; nm1d2 < cdiff; nm1d2++) {
    t->data[nm1d2] /= fs;
  }

  for (nm1d2 = 0; nm1d2 < 1025; nm1d2++) {
    f[nm1d2] = (double)nm1d2 * fs / 2048.0;
  }
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_real_T *y
// Return Type  : void
//
static void sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int sz[2];
  int vstride;
  int j;
  double s;
  int k;
  for (vstride = 0; vstride < 2; vstride++) {
    sz[vstride] = (unsigned int)x->size[vstride];
  }

  vstride = y->size[0];
  y->size[0] = (int)sz[0];
  emxEnsureCapacity((emxArray__common *)y, vstride, (int)sizeof(double));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    j = y->size[0];
    vstride = y->size[0];
    y->size[0] = j;
    emxEnsureCapacity((emxArray__common *)y, vstride, (int)sizeof(double));
    for (vstride = 0; vstride < j; vstride++) {
      y->data[vstride] = 0.0;
    }
  } else {
    vstride = x->size[0];
    for (j = 0; j + 1 <= vstride; j++) {
      s = x->data[j];
      for (k = 2; k <= x->size[1]; k++) {
        s += x->data[j + (k - 1) * vstride];
      }

      y->data[j] = s;
    }
  }
}

//
// treeClassifier - For SSVEP Classification.
// Arguments    : const double F[30]
//                const emxArray_real_T *F2
//                double Y[7]
// Return Type  : void
//
static void treeClassifier(const double F[30], const emxArray_real_T *F2, double
  Y[7])
{
  int i;
  emxArray_real_T *A2;
  int pEnd;
  double dv22[4];
  double b_F[4];
  static const signed char iv1[4] = { 26, 27, 28, 29 };

  static const signed char iv2[4] = { 10, 12, 15, 16 };

  boolean_T p;
  int iwork[4];
  int idx[4];
  int k;
  static const signed char iv3[4] = { 0, 1, 2, 3 };

  emxArray_real_T *unqwLFFT;
  int khi;
  int j;
  static const signed char iv4[2] = { 1, 4 };

  int b_p;
  int b_k;
  int nb;
  int kEnd;
  double x;
  int exitg4;
  int i9;
  static const signed char iv5[4] = { 4, 5, 6, 7 };

  emxArray_real_T *unqwPSD;
  int exitg3;
  int i10;
  double stft_sel_loc[4];
  double fft_sel_loc[16];
  double psd_sel_loc[16];
  double fft_sel_pks[16];
  double psd_sel_pks[16];
  boolean_T b_fft_sel_loc[16];
  double sumB1[4];
  double sumB2[4];
  boolean_T B1_1[4];
  boolean_T b_B1_1[4];
  boolean_T B2_1[4];
  emxArray_int32_T *ndx;
  boolean_T c_B1_1;
  boolean_T b_B2_1;
  boolean_T exitg2;
  boolean_T exitg1;

  // if L = 72
  // % F1 %# OF FEATURES IS CONSTANT:
  for (i = 0; i < 7; i++) {
    Y[i] = 0.0;
  }

  emxInit_real_T(&A2, 2);

  // short classifier
  //  unpack variables:
  //  averageFFTL = F(17); %May not be very useful.
  //  averagePSDL = F(18);
  //  Locations of major peaks
  pEnd = A2->size[0] * A2->size[1];
  A2->size[0] = 1;
  A2->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)A2, pEnd, (int)sizeof(double));
  A2->data[0] = 0.0;
  if (F[26] + F[27] == 2.0) {
    // FFTs match
    for (i = 0; i < 4; i++) {
      dv22[0] = 1.0 + (double)i;
      dv22[1] = 1.0 + (double)i;
      dv22[2] = 1.0 + (double)i;
      dv22[3] = 1.0 + (double)i;
      if (isequal(*(double (*)[4])&F[0], dv22)) {
        Y[1] = iv2[i];
      }
    }
  }

  if (F[28] + F[29] == 2.0) {
    for (i = 0; i < 4; i++) {
      b_F[0] = 1.0 + (double)i;
      b_F[1] = 1.0 + (double)i;
      b_F[2] = 1.0 + (double)i;
      b_F[3] = 1.0 + (double)i;
      if (isequal(*(double (*)[4])&F[4], b_F)) {
        Y[2] = iv2[i];
      }
    }
  } else {
    Y[2] = 0.0;
  }

  for (pEnd = 0; pEnd < 4; pEnd++) {
    b_F[pEnd] = F[iv1[pEnd]];
  }

  if (b_sum(b_F) == 4.0) {
    if ((Y[1] == Y[2]) && (Y[1] != 0.0)) {
      p = true;
    } else {
      p = false;
    }

    Y[3] = p;

    //  IF Frequencies match up entirely.
  }

  //      FFTMatch
  for (k = 0; k <= 3; k += 2) {
    if ((F[iv3[k]] <= F[iv3[k + 1]]) || rtIsNaN(F[k + 1])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  i = 2;
  while (i < 4) {
    khi = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 5; pEnd = nb + i) {
      b_p = j;
      b_k = pEnd - 1;
      nb = j + khi;
      k = 0;
      kEnd = nb - j;
      while (k + 1 <= kEnd) {
        if ((F[iv3[idx[b_p - 1] - 1]] <= F[iv3[idx[b_k] - 1]]) || rtIsNaN
            (F[idx[b_k] - 1])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (b_k + 1 < nb) {
              k++;
              iwork[k] = idx[b_k];
              b_k++;
            }
          }
        } else {
          iwork[k] = idx[b_k];
          b_k++;
          if (b_k + 1 == nb) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }

      j = nb;
    }

    i = khi;
  }

  emxInit_real_T(&unqwLFFT, 2);
  for (pEnd = 0; pEnd < 2; pEnd++) {
    khi = unqwLFFT->size[0] * unqwLFFT->size[1];
    unqwLFFT->size[pEnd] = iv4[pEnd];
    emxEnsureCapacity((emxArray__common *)unqwLFFT, khi, (int)sizeof(double));
  }

  for (k = 0; k < 4; k++) {
    unqwLFFT->data[k] = F[idx[k] - 1];
  }

  k = 0;
  while ((k + 1 <= 4) && rtIsInf(unqwLFFT->data[k]) && (unqwLFFT->data[k] < 0.0))
  {
    k++;
  }

  b_k = k;
  k = 4;
  while ((k >= 1) && rtIsNaN(unqwLFFT->data[k - 1])) {
    k--;
  }

  b_p = 4 - k;
  while ((k >= 1) && rtIsInf(unqwLFFT->data[k - 1]) && (unqwLFFT->data[k - 1] >
          0.0)) {
    k--;
  }

  pEnd = 4 - (k + b_p);
  nb = -1;
  if (b_k > 0) {
    nb = 0;
  }

  khi = (b_k + k) - b_k;
  while (b_k + 1 <= khi) {
    x = unqwLFFT->data[b_k];
    do {
      exitg4 = 0;
      b_k++;
      if (b_k + 1 > khi) {
        exitg4 = 1;
      } else {
        if ((fabs(x - unqwLFFT->data[b_k]) < eps(x / 2.0)) || (rtIsInf
             (unqwLFFT->data[b_k]) && rtIsInf(x) && ((unqwLFFT->data[b_k] > 0.0)
              == (x > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (!p) {
          exitg4 = 1;
        }
      }
    } while (exitg4 == 0);

    nb++;
    unqwLFFT->data[nb] = x;
  }

  if (pEnd > 0) {
    nb++;
    unqwLFFT->data[nb] = unqwLFFT->data[khi];
  }

  b_k = khi + pEnd;
  for (j = 1; j <= b_p; j++) {
    nb++;
    unqwLFFT->data[nb] = unqwLFFT->data[(b_k + j) - 1];
  }

  pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
  if (1 > nb + 1) {
    i9 = -1;
  } else {
    i9 = nb;
  }

  unqwLFFT->size[1] = i9 + 1;
  emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
  for (k = 0; k <= 3; k += 2) {
    if ((F[iv5[k]] <= F[iv5[k + 1]]) || rtIsNaN(F[k + 5])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  i = 2;
  while (i < 4) {
    khi = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 5; pEnd = nb + i) {
      b_p = j;
      b_k = pEnd - 1;
      nb = j + khi;
      k = 0;
      kEnd = nb - j;
      while (k + 1 <= kEnd) {
        if ((F[iv5[idx[b_p - 1] - 1]] <= F[iv5[idx[b_k] - 1]]) || rtIsNaN
            (F[idx[b_k] + 3])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (b_k + 1 < nb) {
              k++;
              iwork[k] = idx[b_k];
              b_k++;
            }
          }
        } else {
          iwork[k] = idx[b_k];
          b_k++;
          if (b_k + 1 == nb) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }

      j = nb;
    }

    i = khi;
  }

  emxInit_real_T(&unqwPSD, 2);
  for (pEnd = 0; pEnd < 2; pEnd++) {
    khi = unqwPSD->size[0] * unqwPSD->size[1];
    unqwPSD->size[pEnd] = iv4[pEnd];
    emxEnsureCapacity((emxArray__common *)unqwPSD, khi, (int)sizeof(double));
  }

  for (k = 0; k < 4; k++) {
    unqwPSD->data[k] = F[idx[k] + 3];
  }

  k = 0;
  while ((k + 1 <= 4) && rtIsInf(unqwPSD->data[k]) && (unqwPSD->data[k] < 0.0))
  {
    k++;
  }

  b_k = k;
  k = 4;
  while ((k >= 1) && rtIsNaN(unqwPSD->data[k - 1])) {
    k--;
  }

  b_p = 4 - k;
  while ((k >= 1) && rtIsInf(unqwPSD->data[k - 1]) && (unqwPSD->data[k - 1] >
          0.0)) {
    k--;
  }

  pEnd = 4 - (k + b_p);
  nb = -1;
  if (b_k > 0) {
    nb = 0;
  }

  khi = (b_k + k) - b_k;
  while (b_k + 1 <= khi) {
    x = unqwPSD->data[b_k];
    do {
      exitg3 = 0;
      b_k++;
      if (b_k + 1 > khi) {
        exitg3 = 1;
      } else {
        if ((fabs(x - unqwPSD->data[b_k]) < eps(x / 2.0)) || (rtIsInf
             (unqwPSD->data[b_k]) && rtIsInf(x) && ((unqwPSD->data[b_k] > 0.0) ==
              (x > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (!p) {
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);

    nb++;
    unqwPSD->data[nb] = x;
  }

  if (pEnd > 0) {
    nb++;
    unqwPSD->data[nb] = unqwPSD->data[khi];
  }

  b_k = khi + pEnd;
  for (j = 1; j <= b_p; j++) {
    nb++;
    unqwPSD->data[nb] = unqwPSD->data[(b_k + j) - 1];
  }

  pEnd = unqwPSD->size[0] * unqwPSD->size[1];
  if (1 > nb + 1) {
    i10 = -1;
  } else {
    i10 = nb;
  }

  unqwPSD->size[1] = i10 + 1;
  emxEnsureCapacity((emxArray__common *)unqwPSD, pEnd, (int)sizeof(double));
  if (b_isequal(unqwLFFT, unqwPSD) && (Y[3] != 0.0)) {
    Y[4] = 1.0;
  }

  emxFree_real_T(&unqwPSD);

  // % F2
  // Unpack shared data:
  for (pEnd = 0; pEnd < 4; pEnd++) {
    stft_sel_loc[pEnd] = 0.0;
  }

  // COMPARE PEAKS
  if (F2->size[1] == 64) {
    for (i = 0; i < 4; i++) {
      // LOC BLOCKS:
      khi = i << 4;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        fft_sel_loc[i + (pEnd << 2)] = F2->data[pEnd + khi];
      }

      khi = i << 4;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        psd_sel_loc[i + (pEnd << 2)] = F2->data[(pEnd + khi) + 4];
      }

      khi = i << 4;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        fft_sel_pks[i + (pEnd << 2)] = F2->data[(pEnd + khi) + 8];
      }

      khi = i << 4;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        psd_sel_pks[i + (pEnd << 2)] = F2->data[(pEnd + khi) + 12];
      }
    }
  } else {
    // numFeatures = 72
    for (i = 0; i < 4; i++) {
      // LOC BLOCKS:
      khi = i * 18;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        fft_sel_loc[i + (pEnd << 2)] = F2->data[pEnd + khi];
      }

      khi = i * 18;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        psd_sel_loc[i + (pEnd << 2)] = F2->data[(pEnd + khi) + 4];
      }

      khi = i * 18;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        fft_sel_pks[i + (pEnd << 2)] = F2->data[(pEnd + khi) + 8];
      }

      khi = i * 18;
      for (pEnd = 0; pEnd < 4; pEnd++) {
        psd_sel_pks[i + (pEnd << 2)] = F2->data[(pEnd + khi) + 12];
      }

      stft_sel_loc[i] = F2->data[i * 18 + 16];
    }
  }

  // % Analysis F2
  for (pEnd = 0; pEnd < 16; pEnd++) {
    b_fft_sel_loc[pEnd] = (fft_sel_loc[pEnd] != 0.0);
  }

  c_sum(b_fft_sel_loc, sumB1);
  for (pEnd = 0; pEnd < 16; pEnd++) {
    b_fft_sel_loc[pEnd] = (psd_sel_loc[pEnd] != 0.0);
  }

  c_sum(b_fft_sel_loc, sumB2);

  // B1_1 & B2_1; %IF ROWS IN FFT AND PSD ARE COMPLETE
  for (i = 0; i < 4; i++) {
    c_B1_1 = (sumB1[i] == 4.0);
    b_B2_1 = (sumB2[i] == 4.0);
    B1_1[i] = (c_B1_1 && b_B2_1);
    b_B1_1[i] = c_B1_1;
    B2_1[i] = b_B2_1;
  }

  emxInit_int32_T(&ndx, 1);
  if (d_sum(B1_1) > 1.0) {
    // COMPARE PEAKS (SUM)
    //      PkValuesToCompare = sum(fft_sel_pks(B_Short,:),2);
    e_sum(fft_sel_pks, sumB1);
    e_sum(psd_sel_pks, sumB2);
    khi = 1;
    x = sumB1[0];
    b_k = 0;
    if (rtIsNaN(sumB1[0])) {
      pEnd = 2;
      exitg2 = false;
      while ((!exitg2) && (pEnd < 5)) {
        khi = pEnd;
        if (!rtIsNaN(sumB1[pEnd - 1])) {
          x = sumB1[pEnd - 1];
          b_k = pEnd - 1;
          exitg2 = true;
        } else {
          pEnd++;
        }
      }
    }

    if (khi < 4) {
      while (khi + 1 < 5) {
        if (sumB1[khi] > x) {
          x = sumB1[khi];
          b_k = khi;
        }

        khi++;
      }
    }

    khi = 1;
    x = sumB2[0];
    b_p = 1;
    if (rtIsNaN(sumB2[0])) {
      pEnd = 2;
      exitg1 = false;
      while ((!exitg1) && (pEnd < 5)) {
        khi = pEnd;
        if (!rtIsNaN(sumB2[pEnd - 1])) {
          x = sumB2[pEnd - 1];
          b_p = pEnd;
          exitg1 = true;
        } else {
          pEnd++;
        }
      }
    }

    if (khi < 4) {
      while (khi + 1 < 5) {
        if (sumB2[khi] > x) {
          x = sumB2[khi];
          b_p = khi + 1;
        }

        khi++;
      }
    }

    if (b_k + 1 == b_p) {
      pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
      unqwLFFT->size[0] = 1;
      unqwLFFT->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
      unqwLFFT->data[0] = iv2[b_k];
    } else {
      pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
      unqwLFFT->size[0] = 1;
      unqwLFFT->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
      unqwLFFT->data[0] = 0.0;
    }
  } else {
    for (i = 0; i < 4; i++) {
      B1_1[i] = (b_B1_1[i] && B2_1[i]);
    }

    x = B1_1[0];
    for (k = 0; k < 3; k++) {
      x += (double)B1_1[k + 1];
    }

    if (x == 1.0) {
      khi = 0;
      for (i = 0; i < 4; i++) {
        if (b_B1_1[i] && B2_1[i]) {
          khi++;
        }
      }

      pEnd = ndx->size[0];
      ndx->size[0] = khi;
      emxEnsureCapacity((emxArray__common *)ndx, pEnd, (int)sizeof(int));
      khi = 0;
      for (i = 0; i < 4; i++) {
        if (b_B1_1[i] && B2_1[i]) {
          ndx->data[khi] = i + 1;
          khi++;
        }
      }

      pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
      unqwLFFT->size[0] = 1;
      unqwLFFT->size[1] = ndx->size[0];
      emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
      khi = ndx->size[0];
      for (pEnd = 0; pEnd < khi; pEnd++) {
        unqwLFFT->data[unqwLFFT->size[0] * pEnd] = iv2[ndx->data[pEnd] - 1];
      }
    } else {
      pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
      unqwLFFT->size[0] = 1;
      unqwLFFT->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
      unqwLFFT->data[0] = 0.0;
    }
  }

  if (unqwLFFT->size[1] == 0) {
    x = 0.0;
  } else {
    x = unqwLFFT->data[0];
    for (k = 2; k <= unqwLFFT->size[1]; k++) {
      x += unqwLFFT->data[k - 1];
    }
  }

  if (x == 1.0) {
    pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
    unqwLFFT->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
    khi = unqwLFFT->size[1];
    for (pEnd = 0; pEnd < khi; pEnd++) {
      unqwLFFT->data[unqwLFFT->size[0] * pEnd] = iv2[(int)unqwLFFT->
        data[unqwLFFT->size[0] * pEnd] - 1];
    }
  }

  if (!(unqwLFFT->size[1] == 0)) {
    if ((unqwLFFT->data[0] == Y[1]) && (Y[3] == 1.0)) {
      Y[0] = 1.0;
    }
  } else {
    Y[0] = 0.0;
  }

  if (F2->size[1] == 72) {
    for (i = 0; i < 4; i++) {
      c_B1_1 = (b_B1_1[i] && B2_1[i]);
      b_B2_1 = (stft_sel_loc[i] != 0.0);
      B1_1[i] = (c_B1_1 && b_B2_1);
      b_B1_1[i] = c_B1_1;
      B2_1[i] = b_B2_1;
    }

    x = B1_1[0];
    for (k = 0; k < 3; k++) {
      x += (double)B1_1[k + 1];
    }

    if (x == 1.0) {
      khi = 0;
      for (i = 0; i < 4; i++) {
        if (b_B1_1[i] && B2_1[i]) {
          khi++;
        }
      }

      pEnd = ndx->size[0];
      ndx->size[0] = khi;
      emxEnsureCapacity((emxArray__common *)ndx, pEnd, (int)sizeof(int));
      khi = 0;
      for (i = 0; i < 4; i++) {
        if (b_B1_1[i] && B2_1[i]) {
          ndx->data[khi] = i + 1;
          khi++;
        }
      }

      pEnd = A2->size[0] * A2->size[1];
      A2->size[0] = 1;
      A2->size[1] = ndx->size[0];
      emxEnsureCapacity((emxArray__common *)A2, pEnd, (int)sizeof(double));
      khi = ndx->size[0];
      for (pEnd = 0; pEnd < khi; pEnd++) {
        A2->data[A2->size[0] * pEnd] = iv2[ndx->data[pEnd] - 1];
      }
    }

    if (unqwLFFT->size[1] == 0) {
      pEnd = unqwLFFT->size[0] * unqwLFFT->size[1];
      unqwLFFT->size[0] = 1;
      unqwLFFT->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)unqwLFFT, pEnd, (int)sizeof(double));
      unqwLFFT->data[0] = 0.0;
    }

    if (A2->size[1] == 0) {
      pEnd = A2->size[0] * A2->size[1];
      A2->size[0] = 1;
      A2->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)A2, pEnd, (int)sizeof(double));
      A2->data[0] = 0.0;
    }

    //      if (A1 == A2) && (Y(4)==1) && (A1 == Y(2))
    if (((int)unqwLFFT->data[0] == (int)A2->data[0]) && ((int)unqwLFFT->data[0]
         != 0)) {
      Y[0] = 1.0;
    } else {
      Y[0] = 0.0;
    }
  }

  emxFree_int32_T(&ndx);
  Y[5] = unqwLFFT->data[0];
  Y[6] = A2->data[0];

  //  Y = [Y];
  // TreeClassifier Function
  emxFree_real_T(&unqwLFFT);
  emxFree_real_T(&A2);
}

//
// if size(signals,2) > size(signals,1)
//      signals = signals.';
//  end
// Arguments    : const emxArray_real_T *signals
//                double fs
//                emxArray_real_T *window
//                emxArray_real_T *CSM
//                emxArray_real_T *frequencies
// Return Type  : void
//
static void welch_psd(const emxArray_real_T *signals, double fs, emxArray_real_T
                      *window, emxArray_real_T *CSM, emxArray_real_T
                      *frequencies)
{
  int k;
  int i5;
  int varargin_1;
  double y;
  double cdiff;
  int ndbl;
  int apnd;
  int loop_ub;
  emxArray_real_T *data_taper;
  double back_shift;
  double number_of_blocks;
  emxArray_real_T *S;
  int a;
  emxArray_real_T *Data_Block;
  emxArray_real_T *P;
  emxArray_creal_T *b_Data_Block;
  emxArray_creal_T *c_Data_Block;
  unsigned int unnamed_idx_0;
  int b_k;
  int c_k;
  emxArray_int32_T *r11;
  emxArray_real_T *b_S;

  // % Function for spectra estimation by Welch's method
  //  Developed by Luiz A. Baccala, Fl?vio Caduda and Luciano Caldas, all from
  //  Escola Polit?cnica - Poli-USP, with cooperation of Carlos Pagani and Felipe 
  //  Amaral from Escola de Engenharia de S?o Carlos - EESC-USP.
  //
  //  Cross-spectra matrix are estimated by Welch's method with 50% overlap and
  //  the window energy loss are compasated by a factor of 1/sum(Wi.^2) where
  //  Wi are the elements of the window [1]. Then, the spectra becomes:
  //  Sxy = fft(x)*conj(fft(y))/sum(Wi.^2)
  //
  //  Code was tested with a known- spectra signal from a white noise filtered
  //  by a filter. The variance (power) of the signal checks with the integral
  //  of the PSD estimated.
  //
  //  INPUT:
  //  -- signals: matrix of signals to perform the spectra estimatino. Size is
  //  [Samples x number of sensors];
  //  -- fs: samplerate in Hertz;
  //  -- window: data taper desired. Must be a vector. For best performance it
  //  should be a power of 2. For general applications do: window=hanning(1024); 
  //
  //  OUTPUT:
  //  -- CSM: Cross Spectral Matrix: Unilateral (0:fs/2) spectra. Welch's
  //  method is used with 50% overlap. Matrix size: sensors x sensors x
  //  windowsize/2
  //  -- frequencies: vector with all frequencies corresponding to each layer
  //  (3rd layer in depth) of CSM.
  //
  //  LAST REVISION: Aug - 18 - 2016
  //  ADDED 'fs' missing term in line 82, for calibration factor
  //  [1] Trobs,M.; Heinzel,G. "Improved spectrum estimation from digitized
  //  time series on a logarithmic frequency axis"
  //  doi:10.1016/j.measurement.2005.10.010
  k = window->size[0];
  i5 = window->size[0];
  window->size[0] = k;
  emxEnsureCapacity((emxArray__common *)window, i5, (int)sizeof(double));
  varargin_1 = window->size[0];
  y = (double)window->size[0] / 2.0;
  cdiff = (double)window->size[0] / 2.0 - 1.0;
  if (y - 1.0 < 0.0) {
    apnd = 0;
    cdiff = y - 1.0;
  } else {
    ndbl = (int)std::floor(cdiff + 0.5);
    apnd = ndbl;
    cdiff = (double)ndbl - (y - 1.0);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(y - 1.0)) {
      ndbl++;
      cdiff = y - 1.0;
    } else if (cdiff > 0.0) {
      cdiff = (double)ndbl - 1.0;
    } else {
      ndbl++;
      cdiff = apnd;
    }

    if (ndbl >= 0) {
      apnd = ndbl;
    } else {
      apnd = 0;
    }
  }

  i5 = frequencies->size[0] * frequencies->size[1];
  frequencies->size[0] = 1;
  frequencies->size[1] = apnd;
  emxEnsureCapacity((emxArray__common *)frequencies, i5, (int)sizeof(double));
  if (apnd > 0) {
    frequencies->data[0] = 0.0;
    if (apnd > 1) {
      frequencies->data[apnd - 1] = cdiff;
      ndbl = (apnd - 1) / 2;
      for (k = 1; k < ndbl; k++) {
        frequencies->data[k] = k;
        frequencies->data[(apnd - k) - 1] = cdiff - (double)k;
      }

      if (ndbl << 1 == apnd - 1) {
        frequencies->data[ndbl] = cdiff / 2.0;
      } else {
        frequencies->data[ndbl] = ndbl;
        frequencies->data[ndbl + 1] = cdiff - (double)ndbl;
      }
    }
  }

  i5 = frequencies->size[0] * frequencies->size[1];
  frequencies->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)frequencies, i5, (int)sizeof(double));
  ndbl = frequencies->size[0];
  apnd = frequencies->size[1];
  k = window->size[0];
  loop_ub = ndbl * apnd;
  for (i5 = 0; i5 < loop_ub; i5++) {
    frequencies->data[i5] = frequencies->data[i5] * fs / (double)k;
  }

  emxInit_real_T1(&data_taper, 1);

  // must be even, best if 2^n
  back_shift = (double)window->size[0] / 2.0;

  // ORIGINAL;
  ndbl = signals->size[1];
  number_of_blocks = std::floor(2.0 * (double)ndbl / (double)window->size[0]) -
    1.0;
  ndbl = window->size[0];
  i5 = data_taper->size[0];
  data_taper->size[0] = ndbl;
  emxEnsureCapacity((emxArray__common *)data_taper, i5, (int)sizeof(double));
  if ((!(window->size[0] == 0)) && (!(ndbl == 0))) {
    for (k = 0; k + 1 <= window->size[0]; k++) {
      data_taper->data[k] = window->data[k];
    }
  }

  emxInit_real_T1(&S, 1);

  //  Data segmentation into blocks of size block_samples:
  y = (double)window->size[0] / 2.0;
  i5 = S->size[0];
  S->size[0] = (int)y;
  emxEnsureCapacity((emxArray__common *)S, i5, (int)sizeof(double));
  loop_ub = (int)y;
  for (i5 = 0; i5 < loop_ub; i5++) {
    S->data[i5] = 0.0;
  }

  // ORIGINAL
  //  S = zeros(ceil(block_samples/2),number_of_signals.^2);
  a = 0;
  emxInit_real_T1(&Data_Block, 1);
  emxInit_real_T1(&P, 1);
  emxInit_creal_T(&b_Data_Block, 1);
  emxInit_creal_T(&c_Data_Block, 1);
  while (a <= (int)number_of_blocks - 1) {
    //  Retrieve current data block
    cdiff = ((1.0 + (double)a) - 1.0) * back_shift + 1.0;
    y = (double)varargin_1 + ((1.0 + (double)a) - 1.0) * back_shift;
    if (cdiff > y) {
      i5 = 0;
      ndbl = 0;
    } else {
      i5 = (int)cdiff - 1;
      ndbl = (int)y;
    }

    apnd = Data_Block->size[0];
    Data_Block->size[0] = ndbl - i5;
    emxEnsureCapacity((emxArray__common *)Data_Block, apnd, (int)sizeof(double));
    loop_ub = ndbl - i5;
    for (apnd = 0; apnd < loop_ub; apnd++) {
      Data_Block->data[apnd] = signals->data[i5 + apnd];
    }

    if (ndbl - i5 == 0) {
      y = 0.0;
    } else {
      y = signals->data[i5];
      for (k = 0; k + 2 <= Data_Block->size[0]; k++) {
        y += signals->data[(i5 + k) + 1];
      }
    }

    y /= (double)Data_Block->size[0];
    i5 = Data_Block->size[0];
    emxEnsureCapacity((emxArray__common *)Data_Block, i5, (int)sizeof(double));
    loop_ub = Data_Block->size[0];
    for (i5 = 0; i5 < loop_ub; i5++) {
      Data_Block->data[i5] -= y;
    }

    i5 = Data_Block->size[0];
    emxEnsureCapacity((emxArray__common *)Data_Block, i5, (int)sizeof(double));
    loop_ub = Data_Block->size[0];
    for (i5 = 0; i5 < loop_ub; i5++) {
      Data_Block->data[i5] *= data_taper->data[i5];
    }

    // Taper it
    fft(Data_Block, b_Data_Block);

    // FFT it,
    //  bilateral DFT
    //  viii
    cdiff = (double)varargin_1 / 2.0;
    if (1.0 > cdiff) {
      loop_ub = 0;
    } else {
      loop_ub = (int)cdiff;
    }

    i5 = c_Data_Block->size[0];
    c_Data_Block->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)c_Data_Block, i5, (int)sizeof(creal_T));
    for (i5 = 0; i5 < loop_ub; i5++) {
      c_Data_Block->data[i5] = b_Data_Block->data[i5];
    }

    i5 = b_Data_Block->size[0];
    b_Data_Block->size[0] = c_Data_Block->size[0];
    emxEnsureCapacity((emxArray__common *)b_Data_Block, i5, (int)sizeof(creal_T));
    loop_ub = c_Data_Block->size[0];
    for (i5 = 0; i5 < loop_ub; i5++) {
      b_Data_Block->data[i5] = c_Data_Block->data[i5];
    }

    // ORIGINAL
    //  Data_Block = Data_Block(1:ceil(block_samples/2),:);
    // All spectral combinations:
    y = (double)varargin_1 / 2.0;
    i5 = P->size[0];
    P->size[0] = (int)y;
    emxEnsureCapacity((emxArray__common *)P, i5, (int)sizeof(double));
    loop_ub = (int)y;
    for (i5 = 0; i5 < loop_ub; i5++) {
      P->data[i5] = 0.0;
    }

    // ORIGINAL
    //  P = zeros(ceil(block_samples/2)/2,number_of_signals.^2);
    //  P(:,c) = Data_Block(:,b).*conj(Data_Block(:,aa)); % THIS
    //  IS FOR WIND TUNNEL EESC-USP BEAMFORMING CODE
    //  P(:,c) = Data_Block(:,aa).*conj(Data_Block(:,b)); % THIS IS THE ORIGINAL 
    //  LINE
    loop_ub = b_Data_Block->size[0] - 1;
    for (i5 = 0; i5 <= loop_ub; i5++) {
      cdiff = b_Data_Block->data[i5].re;
      y = -b_Data_Block->data[i5].im;
      cdiff = b_Data_Block->data[i5].re * cdiff - b_Data_Block->data[i5].im * y;
      P->data[i5] = cdiff;
    }

    //  IS FOR FAN RIG BEAMFORMING CODE
    //  Sum the spectrums up ...
    i5 = S->size[0];
    emxEnsureCapacity((emxArray__common *)S, i5, (int)sizeof(double));
    loop_ub = S->size[0];
    for (i5 = 0; i5 < loop_ub; i5++) {
      S->data[i5] += P->data[i5];
    }

    a++;
  }

  emxFree_creal_T(&c_Data_Block);
  emxFree_creal_T(&b_Data_Block);
  emxFree_real_T(&data_taper);
  i5 = S->size[0];
  emxEnsureCapacity((emxArray__common *)S, i5, (int)sizeof(double));
  loop_ub = S->size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    S->data[i5] *= 2.0;
  }

  i5 = Data_Block->size[0];
  Data_Block->size[0] = window->size[0];
  emxEnsureCapacity((emxArray__common *)Data_Block, i5, (int)sizeof(double));
  loop_ub = window->size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    Data_Block->data[i5] = window->data[i5];
  }

  unnamed_idx_0 = (unsigned int)window->size[0];
  i5 = P->size[0];
  P->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)P, i5, (int)sizeof(double));
  ndbl = window->size[0];

//#pragma omp parallel for \
// num_threads(omp_get_max_threads()) \
// private(c_k)

  for (b_k = 1; b_k <= ndbl; b_k++) {
    c_k = b_k;
    P->data[c_k - 1] = Data_Block->data[c_k - 1] * Data_Block->data[c_k - 1];
  }

  emxFree_real_T(&Data_Block);
  if (P->size[0] == 0) {
    y = 0.0;
  } else {
    y = P->data[0];
    for (k = 2; k <= P->size[0]; k++) {
      y += P->data[k - 1];
    }
  }

  emxFree_real_T(&P);
  cdiff = y * fs * number_of_blocks;
  i5 = S->size[0];
  emxEnsureCapacity((emxArray__common *)S, i5, (int)sizeof(double));
  loop_ub = S->size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    S->data[i5] /= cdiff;
  }

  //  Average them out
  i5 = CSM->size[0] * CSM->size[1];
  CSM->size[0] = 1;
  CSM->size[1] = S->size[0];
  emxEnsureCapacity((emxArray__common *)CSM, i5, (int)sizeof(double));
  loop_ub = S->size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    CSM->data[i5] = 0.0;
  }

  emxInit_int32_T(&r11, 1);

  //  for a = 1:sensors
  ndbl = S->size[0];
  i5 = r11->size[0];
  r11->size[0] = ndbl;
  emxEnsureCapacity((emxArray__common *)r11, i5, (int)sizeof(int));
  for (i5 = 0; i5 < ndbl; i5++) {
    r11->data[i5] = i5;
  }

  emxInit_real_T1(&b_S, 1);
  loop_ub = S->size[0];
  i5 = b_S->size[0];
  b_S->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)b_S, i5, (int)sizeof(double));
  for (i5 = 0; i5 < loop_ub; i5++) {
    b_S->data[i5] = S->data[i5];
  }

  emxFree_real_T(&S);
  ndbl = r11->size[0];
  for (i5 = 0; i5 < ndbl; i5++) {
    CSM->data[CSM->size[0] * r11->data[i5]] = b_S->data[i5];
  }

  emxFree_real_T(&b_S);
  emxFree_int32_T(&r11);

  //  end
  //  clear S
  CSM->data[0] = (CSM->data[0] + CSM->data[0]) - CSM->data[0];
}

//
// Arguments    : int numDimensions
//                int *size
// Return Type  : emxArray_real_T *
//
emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T1(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

//
// Arguments    : double *data
//                int numDimensions
//                int *size
// Return Type  : emxArray_real_T *
//
emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
  size)
{
  emxArray_real_T *emx;
  int numEl;
  int i;
  emxInit_real_T1(&emx, numDimensions);
  numEl = 1;
  for (i = 0; i < numDimensions; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

//
// Arguments    : double *data
//                int rows
//                int cols
// Return Type  : emxArray_real_T *
//
emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T1(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = false;
  return emx;
}

//
// Arguments    : int rows
//                int cols
// Return Type  : emxArray_real_T *
//
emxArray_real_T *emxCreate_real_T(int rows, int cols)
{
  emxArray_real_T *emx;
  int size[2];
  int numEl;
  int i;
  size[0] = rows;
  size[1] = cols;
  emxInit_real_T1(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

//
// Arguments    : emxArray_real_T *emxArray
// Return Type  : void
//
void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxInit_real_T1(pEmxArray, numDimensions);
}

//
// Full hybrid EOG/EEG classifier
//  Ch1 = Fp1
//  Ch2 = Fp2
//  Ch3 = Fpz
//  Ch4 = Right Eye
//
//  tXEOG = Training Data for EOG
//  tYEOG = Classes for EOG
//
//  EOGOnly is a boolean value that is used only if EOG function is to be
//  called.
//
//  The first part of this classifier determines if an emergency stop is
//  requested in the form of a double-blink by the subject:
//  A '1' is a double blink.
//  Any other result is a pass, and classification continues.
//  The SSVEP Portion will output one of the following corresponding to its
//  frequency makeup:
//  CLASS :: CORRESPONDING FREQ (ACTUAL)
//  10    :: 10.00Hz
//  12    :: 12.50Hz
//  15    :: 15.15Hz
//  16    :: 16.66Hz
// Arguments    : const emxArray_real_T *ch1
//                const emxArray_real_T *ch2
//                const emxArray_real_T *ch3
//                const emxArray_real_T *ch4
//                double Fs
//                boolean_T EOGOnly
//                double Y[7]
// Return Type  : void
//
void fullHybridClassifier(const emxArray_real_T *ch1, const emxArray_real_T *ch2,
  const emxArray_real_T *ch3, const emxArray_real_T *ch4, double Fs, boolean_T
  EOGOnly, double Y[7])
{
  int k;
  emxArray_real_T *F2;
  int i0;
  double ch4f[250];
  double dv0[250];
  double ch1f[250];
  double dv1[250];
  double ch2f[250];
  double dv2[250];
  double ch3f[250];
  emxArray_real_T *r0;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *r3;
  double dv3[250];
  double dv4[250];
  double dv5[40];
  boolean_T DB;
  boolean_T exitg1;
  emxArray_real_T *sch1f;
  emxArray_real_T *sch2f;
  emxArray_real_T *sch3f;
  emxArray_real_T *r4;
  emxArray_real_T *b_ch3;
  emxArray_real_T *b_ch2;
  emxArray_real_T *b_ch1;
  double F[30];

  //  Window length??
  //  Y = 0;
  for (k = 0; k < 7; k++) {
    Y[k] = 0.0;
  }

  // Default Value.
  //  DB = false;
  // Load Training Data:
  //  numberSSVEPFeaturesLong = ??
  k = ch1->size[0];
  emxInit_real_T(&F2, 2);
  if (k < 500) {
    i0 = F2->size[0] * F2->size[1];
    F2->size[0] = 1;
    F2->size[1] = 64;
    emxEnsureCapacity((emxArray__common *)F2, i0, (int)sizeof(double));
    for (i0 = 0; i0 < 64; i0++) {
      F2->data[i0] = 0.0;
    }
  } else {
    i0 = F2->size[0] * F2->size[1];
    F2->size[0] = 1;
    F2->size[1] = 72;
    emxEnsureCapacity((emxArray__common *)F2, i0, (int)sizeof(double));
    for (i0 = 0; i0 < 72; i0++) {
      F2->data[i0] = 0.0;
    }
  }

  // Temporary:
  //  F = zeros(1,100);
  k = ch1->size[0];
  if (k >= 250) {
    //  Filter using optimized EOG filter:
    //  (take last second, regardless of actual window length):
    k = ch1->size[0];
    for (i0 = 0; i0 < 250; i0++) {
      ch4f[i0] = ch1->data[(i0 + k) - 250];
    }

    // EOGCFILT EOG Filter for conversion to C. All inputs must be constant.
    //  Vectorize:
    //  Sampling Frequency = 250;
    // BW for 10Hz upper bound, Order of 3.
    // BW filt for 2Hz lower bound, Order of 3:
    filtfilt(ch4f, dv0);
    b_filtfilt(dv0, ch1f);
    k = ch2->size[0];
    for (i0 = 0; i0 < 250; i0++) {
      ch4f[i0] = ch2->data[(i0 + k) - 250];
    }

    // EOGCFILT EOG Filter for conversion to C. All inputs must be constant.
    //  Vectorize:
    //  Sampling Frequency = 250;
    // BW for 10Hz upper bound, Order of 3.
    // BW filt for 2Hz lower bound, Order of 3:
    filtfilt(ch4f, dv1);
    b_filtfilt(dv1, ch2f);
    k = ch3->size[0];
    for (i0 = 0; i0 < 250; i0++) {
      ch4f[i0] = ch3->data[(i0 + k) - 250];
    }

    // EOGCFILT EOG Filter for conversion to C. All inputs must be constant.
    //  Vectorize:
    //  Sampling Frequency = 250;
    // BW for 10Hz upper bound, Order of 3.
    // BW filt for 2Hz lower bound, Order of 3:
    filtfilt(ch4f, dv2);
    b_filtfilt(dv2, ch3f);
    k = ch4->size[0];
    for (i0 = 0; i0 < 250; i0++) {
      ch4f[i0] = ch4->data[(i0 + k) - 250];
    }

    emxInit_real_T(&r0, 2);
    emxInit_real_T(&r1, 2);
    emxInit_real_T(&r2, 2);
    emxInit_real_T(&r3, 2);

    // EOGCFILT EOG Filter for conversion to C. All inputs must be constant.
    //  Vectorize:
    //  Sampling Frequency = 250;
    // BW for 10Hz upper bound, Order of 3.
    // BW filt for 2Hz lower bound, Order of 3:
    // Extract EOG features: (1s window)
    featureExtractionEOG(ch1f, r0);
    featureExtractionEOG(ch2f, r1);
    featureExtractionEOG(ch3f, r2);
    filtfilt(ch4f, dv3);
    b_filtfilt(dv3, dv4);
    featureExtractionEOG(dv4, r3);

    // Combine features:
    // Boolean DB: represents presence of a double blink.
    for (i0 = 0; i0 < 10; i0++) {
      dv5[i0] = r0->data[i0];
    }

    emxFree_real_T(&r0);
    for (i0 = 0; i0 < 10; i0++) {
      dv5[i0 + 10] = r1->data[i0];
    }

    emxFree_real_T(&r1);
    for (i0 = 0; i0 < 10; i0++) {
      dv5[i0 + 20] = r2->data[i0];
    }

    emxFree_real_T(&r2);
    for (i0 = 0; i0 < 10; i0++) {
      dv5[i0 + 30] = r3->data[i0];
    }

    emxFree_real_T(&r3);
    Y[0] = knn(dv5);
    DB = false;
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k + 1 < 8)) {
      if (fabs(Y[k] - 1.0) < eps(Y[k] / 2.0)) {
        DB = true;
        exitg1 = true;
      } else {
        k++;
      }
    }

    // IF 1 is a member of B
    //  SSVEP CLASSIFICATION:
    //  PRECONDITIONS: EOG must not have been triggered.
    //  starts with 1/2 second analysis and moves up.
    //  Output can be one of the four SSVEP classes [10 12 15 16]
    emxInit_real_T(&sch1f, 2);
    emxInit_real_T(&sch2f, 2);
    emxInit_real_T(&sch3f, 2);
    emxInit_real_T(&r4, 2);
    emxInit_real_T1(&b_ch3, 1);
    emxInit_real_T1(&b_ch2, 1);
    emxInit_real_T1(&b_ch1, 1);
    if ((!DB) && (!EOGOnly)) {
      // If no double blink has been detected in final second of data.
      //  Use a decision tree.
      // Y = knn(tsX, tX, tY, 1); %Fine KNN
      // limit size:
      // ln = min([size(ch1,1) size(ch2,1) size(ch3,1)]);
      //  Make sure len is even:
      // if mod(ln,2)~=0
      //     ln=ln-1;
      // end
      // Filter:
      k = ch1->size[0];
      i0 = b_ch1->size[0];
      b_ch1->size[0] = k;
      emxEnsureCapacity((emxArray__common *)b_ch1, i0, (int)sizeof(double));
      for (i0 = 0; i0 < k; i0++) {
        b_ch1->data[i0] = ch1->data[i0];
      }

      eegcfilt(b_ch1, sch1f);
      k = ch2->size[0];
      i0 = b_ch2->size[0];
      b_ch2->size[0] = k;
      emxEnsureCapacity((emxArray__common *)b_ch2, i0, (int)sizeof(double));
      for (i0 = 0; i0 < k; i0++) {
        b_ch2->data[i0] = ch2->data[i0];
      }

      eegcfilt(b_ch2, sch2f);
      k = ch3->size[0];
      i0 = b_ch3->size[0];
      b_ch3->size[0] = k;
      emxEnsureCapacity((emxArray__common *)b_ch3, i0, (int)sizeof(double));
      for (i0 = 0; i0 < k; i0++) {
        b_ch3->data[i0] = ch3->data[i0];
      }

      eegcfilt(b_ch3, sch3f);
      featureExtractionSSVEP(sch1f, sch2f, sch3f, Fs, F);
      fSSVEPnew(sch1f, sch2f, sch3f, Fs, r4);
      k = r4->size[1];
      for (i0 = 0; i0 < k; i0++) {
        F2->data[F2->size[0] * i0] = r4->data[r4->size[0] * i0];
      }

      //  Extract SSVEP Features (Part 1 from individual channels):
      // Different window lengths correspond to different classifiers.
      k = ch1->size[0];
      if (k < 500) {
        // numFeats <= 53
        treeClassifier(F, F2, Y);
      } else {
        k = ch1->size[0];
        if (k >= 500) {
          //  Window length is 500 or more samples.
          //  use separate classifier for this.
          treeClassifier(F, F2, Y);
        }
      }
    }

    emxFree_real_T(&b_ch1);
    emxFree_real_T(&b_ch2);
    emxFree_real_T(&b_ch3);
    emxFree_real_T(&r4);
    emxFree_real_T(&sch3f);
    emxFree_real_T(&sch2f);
    emxFree_real_T(&sch1f);
  }
  emxFree_real_T(&F2);
}
//
// Arguments    : void
// Return Type  : void
//
void fullHybridClassifier_initialize()
{
  rt_InitInfAndNaN(8U);
//  omp_init_nest_lock(&emlrtNestLockGlobal);
}

//
// Arguments    : void
// Return Type  : void
//
void fullHybridClassifier_terminate()
{
//  omp_destroy_nest_lock(&emlrtNestLockGlobal);
}



//
// File trailer for fullHybridClassifier.cpp
//
// [EOF]
//
