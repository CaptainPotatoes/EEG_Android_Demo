//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fullHybridClassifier.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Mar-2017 00:59:20
//

// Include Files
#include "rt_nonfinite.h"
#include "fullHybridClassifier.h"
#include <cmath>
#include <cstdlib>
/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define LOG_TAG "fullHybridClassifier-cpp"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)
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
omp_nest_lock_t emlrtNestLockGlobal;

// Function Declarations
static double WCountMax(const double X[250]);
static double WCountMin(const double X[250]);
static double Wmax(const double X[250]);
static double Wmean(const double X[250]);
static double Wmin(const double X[250]);
static double Wstd(const double X[250]);
static void assignOutputs(const double y[250], const double x[250], const double
  iPk_data[], const int iPk_size[1], double YpkOut_data[], int YpkOut_size[2],
  double XpkOut_data[], int XpkOut_size[2]);
static void b_abs(const creal_T x[2048], double y[2048]);
static void b_assignOutputs(const double y[1025], const double x[1025], const
  double iPk_data[], const int iPk_size[1], double YpkOut_data[], int
  YpkOut_size[2], double XpkOut_data[], int XpkOut_size[2]);
static void b_fft(const emxArray_real_T *x, emxArray_creal_T *y);
static void b_filter(const double x_data[], const int x_size[1], const double
                     zi[10], double y_data[], int y_size[1]);
static void b_filtfilt(const double x_in[250], double y_out[250]);
static void b_findLocalMaxima(const double yTemp[1025], double iPk_data[], int
  iPk_size[1], double iInflect_data[], int iInflect_size[1]);
static void b_findpeaks(const double Yin[1025], double Ypk_data[], int Ypk_size
  [2], double Xpk_data[], int Xpk_size[2]);
static void b_getAllPeaks(const double y[1025], double iPk_data[], int iPk_size
  [1], double iInf_data[], int iInf_size[1], double iInflect_data[], int
  iInflect_size[1]);
static void b_keepAtMostNpPeaks(double idx_data[], int idx_size[1]);
static double b_knn(const double tsX[42], const emxArray_real_T *tX, const
                    emxArray_real_T *tY);
static double b_mean(const double x[4]);
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int
                    np, int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static double b_mod(double x);
static void b_r2br_r2dit_trig(const emxArray_real_T *x, int n1_unsigned, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T *y);
static void b_removePeaksBelowMinPeakHeight(const double Y[1025], double
  iPk_data[], int iPk_size[1]);
static void b_removePeaksBelowThreshold(const double Y[1025], double iPk_data[],
  int iPk_size[1]);
static void b_repmat(double a, double varargin_1, emxArray_real_T *b);
static void b_sign(emxArray_real_T *x);
static void b_sort(emxArray_real_T *x, emxArray_int32_T *idx);
static void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
static double b_sum(const emxArray_real_T *x);
static void bluestein_setup(int nRows, emxArray_creal_T *wwc);
static void c_assignOutputs(const emxArray_real_T *y, const emxArray_real_T *x,
  const emxArray_real_T *iPk, emxArray_real_T *YpkOut, emxArray_real_T *XpkOut);
static void c_findLocalMaxima(emxArray_real_T *yTemp, emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void c_findPeaksSeparatedByMoreThanM(const int iPk_size[1], double
  idx_data[], int idx_size[1]);
static void c_findpeaks(const emxArray_real_T *Yin, emxArray_real_T *Ypk,
  emxArray_real_T *Xpk);
static void c_getAllPeaks(const emxArray_real_T *y, emxArray_real_T *iPk,
  emxArray_real_T *iInf, emxArray_real_T *iInflect);
static void c_keepAtMostNpPeaks(emxArray_real_T *idx, double Np);
static void c_r2br_r2dit_trig(const emxArray_creal_T *x, int n1_unsigned, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T *y);
static void c_removePeaksBelowMinPeakHeight(const emxArray_real_T *Y,
  emxArray_real_T *iPk);
static void c_removePeaksBelowThreshold(const emxArray_real_T *Y,
  emxArray_real_T *iPk);
static void c_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx);
static void d_findPeaksSeparatedByMoreThanM(const int iPk_size[1], double
  idx_data[], int idx_size[1]);
static void d_r2br_r2dit_trig(const emxArray_creal_T *x, int n1_unsigned, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T *y);
static void diff(const emxArray_real_T *x, emxArray_real_T *y);
static void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib);
static void dobluesteinfft(const emxArray_real_T *x, int N2, int n1, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, const emxArray_real_T *
  sintabinv, emxArray_creal_T *y);
static void e_findPeaksSeparatedByMoreThanM(const emxArray_real_T *iPk,
  emxArray_real_T *idx);
static void eegcfilt(double X_data[], int X_size[1], emxArray_real_T *Y);
static void eml_fft(const emxArray_real_T *x, creal_T y[2048]);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
static void emxFree_creal_T(emxArray_creal_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
static void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static void eogcfilt(const double X[250], double Y[250]);
static void featureExtractionEOG(const double samplesX[250], double F_data[],
  int F_size[2]);
static void featureExtractionSSVEP(const emxArray_real_T *fch1, const
  emxArray_real_T *fch2, const emxArray_real_T *fch3, double Fs, double F[42]);
static void fft(const emxArray_real_T *x, creal_T y[2048]);
static void filter(const double b[4], const double a[4], const double x[268],
                   const double zi[3], double y[268]);
static void filtfilt(const double x_in[250], double y_out[250]);
static void findLocalMaxima(const double yTemp[250], double iPk_data[], int
  iPk_size[1], double iInflect_data[], int iInflect_size[1]);
static int findbin(double x, const emxArray_real_T *bin_edges);
static void findpeaks(const double Yin[250], double Ypk_data[], int Ypk_size[2],
                      double Xpk_data[], int Xpk_size[2]);
static void flipud(double x_data[], int x_size[1]);
static void generate_twiddle_tables(int nRows, boolean_T useRadix2,
  emxArray_real_T *costab, emxArray_real_T *sintab, emxArray_real_T *sintabinv);
static void getAllPeaks(const double y[250], double iPk_data[], int iPk_size[1],
  double iInf_data[], int iInf_size[1], double iInflect_data[], int
  iInflect_size[1]);
static void get_algo_sizes(int n1, boolean_T useRadix2, int *N2blue, int *nRows);
static void get_nfft_data(const emxArray_real_T *X, double Fs, double f[1025],
  double C[1025]);
static void hannWin(double x, emxArray_real_T *w);
static boolean_T isequal(const double varargin_1[4], const double varargin_2[4]);
static void keepAtMostNpPeaks(double idx_data[], int idx_size[1]);
static double knn(const double tsX[40], const emxArray_real_T *tX, const
                  emxArray_real_T *tY);
static double mean(const emxArray_real_T *x);
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_pow2_block(emxArray_int32_T *idx, emxArray_real_T *x, int
  offset);
static int nonSingletonDim(const emxArray_real_T *x);
static void orderPeaks(const double Y[1025], const double iPk_data[], double
  idx_data[], int idx_size[1]);
static void parse_inputs(const emxArray_real_T *Yin, emxArray_real_T *y,
  emxArray_real_T *x, double *NpOut);
static void power(const emxArray_real_T *a, emxArray_real_T *y);
static void r2br_r2dit_trig(const emxArray_real_T *x, const double costab[1025],
  const double sintab[1025], creal_T y[2048]);
static void r2br_r2dit_trig_impl(const emxArray_creal_T *x, int unsigned_nRows,
  const emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T
  *y);
static void removePeaksBelowMinPeakHeight(const double Y[250], double iPk_data[],
  int iPk_size[1]);
static void removePeaksBelowThreshold(const double Y[250], double iPk_data[],
  int iPk_size[1]);
static void repmat(const emxArray_real_T *a, emxArray_real_T *b);
static double rt_hypotd_snf(double u0, double u1);
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);
static void sort(emxArray_real_T *x, emxArray_int32_T *idx);
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
static double sum(const boolean_T x[250]);
static double trapz(const double x[250]);
static void welch_psd(const emxArray_real_T *signals, double fs, emxArray_real_T
                      *window, emxArray_real_T *CSM, emxArray_real_T
                      *frequencies);

// Function Definitions

//
// Arguments    : const double X[250]
// Return Type  : double
//
static double WCountMax(const double X[250])
{
  double Y;
  boolean_T x[250];
  int k;
  for (k = 0; k < 250; k++) {
    x[k] = (X[k] > 8.4999999999999993E-5);
  }

  Y = x[0];
  for (k = 0; k < 249; k++) {
    Y += (double)x[k + 1];
  }

  return Y;
}

//
// Arguments    : const double X[250]
// Return Type  : double
//
static double WCountMin(const double X[250])
{
  double Y;
  boolean_T x[250];
  int k;
  for (k = 0; k < 250; k++) {
    x[k] = (X[k] < -0.0001);
  }

  Y = x[0];
  for (k = 0; k < 249; k++) {
    Y += (double)x[k + 1];
  }

  return Y;
}

//
// Arguments    : const double X[250]
// Return Type  : double
//
static double Wmax(const double X[250])
{
  int ixstart;
  double mtmp;
  int ix;
  boolean_T exitg1;
  ixstart = 1;
  mtmp = X[0];
  if (rtIsNaN(X[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 251)) {
      ixstart = ix;
      if (!rtIsNaN(X[ix - 1])) {
        mtmp = X[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 250) {
    while (ixstart + 1 < 251) {
      if (X[ixstart] > mtmp) {
        mtmp = X[ixstart];
      }

      ixstart++;
    }
  }

  return mtmp;
}

//
// Arguments    : const double X[250]
// Return Type  : double
//
static double Wmean(const double X[250])
{
  double y;
  int k;
  y = X[0];
  for (k = 0; k < 249; k++) {
    y += X[k + 1];
  }

  return y / 250.0;
}

//
// Arguments    : const double X[250]
// Return Type  : double
//
static double Wmin(const double X[250])
{
  int ixstart;
  double mtmp;
  int ix;
  boolean_T exitg1;
  ixstart = 1;
  mtmp = X[0];
  if (rtIsNaN(X[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 251)) {
      ixstart = ix;
      if (!rtIsNaN(X[ix - 1])) {
        mtmp = X[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 250) {
    while (ixstart + 1 < 251) {
      if (X[ixstart] < mtmp) {
        mtmp = X[ixstart];
      }

      ixstart++;
    }
  }

  return mtmp;
}

//
// Arguments    : const double X[250]
// Return Type  : double
//
static double Wstd(const double X[250])
{
  int ix;
  double xbar;
  int k;
  double r;
  double y;
  ix = 0;
  xbar = X[0];
  for (k = 0; k < 249; k++) {
    ix++;
    xbar += X[ix];
  }

  xbar /= 250.0;
  ix = 0;
  r = X[0] - xbar;
  y = r * r;
  for (k = 0; k < 249; k++) {
    ix++;
    r = X[ix] - xbar;
    y += r * r;
  }

  y /= 249.0;
  return std::sqrt(y);
}

//
// Arguments    : const double y[250]
//                const double x[250]
//                const double iPk_data[]
//                const int iPk_size[1]
//                double YpkOut_data[]
//                int YpkOut_size[2]
//                double XpkOut_data[]
//                int XpkOut_size[2]
// Return Type  : void
//
static void assignOutputs(const double y[250], const double x[250], const double
  iPk_data[], const int iPk_size[1], double YpkOut_data[], int YpkOut_size[2],
  double XpkOut_data[], int XpkOut_size[2])
{
  int loop_ub;
  int i4;
  YpkOut_size[0] = 1;
  YpkOut_size[1] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i4 = 0; i4 < loop_ub; i4++) {
    YpkOut_data[i4] = y[(int)iPk_data[i4] - 1];
  }

  XpkOut_size[0] = 1;
  XpkOut_size[1] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i4 = 0; i4 < loop_ub; i4++) {
    XpkOut_data[i4] = x[(int)iPk_data[i4] - 1];
  }
}

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
// Arguments    : const double y[1025]
//                const double x[1025]
//                const double iPk_data[]
//                const int iPk_size[1]
//                double YpkOut_data[]
//                int YpkOut_size[2]
//                double XpkOut_data[]
//                int XpkOut_size[2]
// Return Type  : void
//
static void b_assignOutputs(const double y[1025], const double x[1025], const
  double iPk_data[], const int iPk_size[1], double YpkOut_data[], int
  YpkOut_size[2], double XpkOut_data[], int XpkOut_size[2])
{
  int loop_ub;
  int i12;
  YpkOut_size[0] = 1;
  YpkOut_size[1] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i12 = 0; i12 < loop_ub; i12++) {
    YpkOut_data[i12] = y[(int)iPk_data[i12] - 1];
  }

  XpkOut_size[0] = 1;
  XpkOut_size[1] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i12 = 0; i12 < loop_ub; i12++) {
    XpkOut_data[i12] = x[(int)iPk_data[i12] - 1];
  }
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_creal_T *y
// Return Type  : void
//
static void b_fft(const emxArray_real_T *x, emxArray_creal_T *y)
{
  emxArray_real_T *costab;
  emxArray_real_T *sintab;
  int N2blue;
  emxArray_real_T *sintabinv;
  boolean_T useRadix2;
  int nRows;
  if (x->size[0] == 0) {
    N2blue = y->size[0];
    y->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)y, N2blue, (int)sizeof(creal_T));
  } else {
    emxInit_real_T(&costab, 2);
    emxInit_real_T(&sintab, 2);
    emxInit_real_T(&sintabinv, 2);
    useRadix2 = ((x->size[0] & (x->size[0] - 1)) == 0);
    get_algo_sizes(x->size[0], useRadix2, &N2blue, &nRows);
    generate_twiddle_tables(nRows, useRadix2, costab, sintab, sintabinv);
    if (useRadix2) {
      b_r2br_r2dit_trig(x, x->size[0], costab, sintab, y);
    } else {
      dobluesteinfft(x, N2blue, x->size[0], costab, sintab, sintabinv, y);
    }

    emxFree_real_T(&sintabinv);
    emxFree_real_T(&sintab);
    emxFree_real_T(&costab);
  }
}

//
// Arguments    : const double x_data[]
//                const int x_size[1]
//                const double zi[10]
//                double y_data[]
//                int y_size[1]
// Return Type  : void
//
static void b_filter(const double x_data[], const int x_size[1], const double
                     zi[10], double y_data[], int y_size[1])
{
  double dbuffer[11];
  int j;
  int k;
  static const double dv6[11] = { 2.13961520749732E-5, 0.0,
    -0.000106980760374866, 0.0, 0.000213961520749732, 0.0, -0.000213961520749732,
    0.0, 0.000106980760374866, 0.0, -2.13961520749732E-5 };

  static const double dv7[11] = { 1.0, -8.77043379286888, 35.0068378010024,
    -83.7229808056309, 132.845833785487, -146.117834417428, 112.823239428442,
    -60.389449129414, 21.4471017127118, -4.56451967201817, 0.442209182399621 };

  y_size[0] = (short)x_size[0];
  memcpy(&dbuffer[1], &zi[0], 10U * sizeof(double));
  for (j = 0; j + 1 <= x_size[0]; j++) {
    for (k = 0; k < 10; k++) {
      dbuffer[k] = dbuffer[k + 1];
    }

    dbuffer[10] = 0.0;
    for (k = 0; k < 11; k++) {
      dbuffer[k] += x_data[j] * dv6[k];
    }

    for (k = 0; k < 10; k++) {
      dbuffer[k + 1] -= dbuffer[0] * dv7[k + 1];
    }

    y_data[j] = dbuffer[0];
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
  double a[3];
  double b_y[268];
  static const double b_a[3] = { -0.95097188792826548, 1.9019437758560462,
    -0.95097188792780118 };

  static const double dv4[4] = { 0.950971887923409, -2.85291566377023,
    2.85291566377023, -0.950971887923409 };

  static const double dv5[4] = { 1.0, -2.89947959461186, 2.803947977383,
    -0.904347531392409 };

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

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 268U * sizeof(double));
  filter(dv4, dv5, b_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&c_y[0], &y[0], 268U * sizeof(double));
  filter(dv4, dv5, c_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  memcpy(&y_out[0], &y[9], 250U * sizeof(double));
}

//
// Arguments    : const double yTemp[1025]
//                double iPk_data[]
//                int iPk_size[1]
//                double iInflect_data[]
//                int iInflect_size[1]
// Return Type  : void
//
static void b_findLocalMaxima(const double yTemp[1025], double iPk_data[], int
  iPk_size[1], double iInflect_data[], int iInflect_size[1])
{
  double b_yTemp[1027];
  boolean_T yFinite[1027];
  int ii;
  boolean_T x[1026];
  int idx;
  short ii_data[1026];
  boolean_T exitg3;
  boolean_T guard3 = false;
  int nx;
  short tmp_data[1027];
  int i10;
  int i11;
  double yTemp_data[1027];
  short iTemp_data[1027];
  int yTemp_size[1];
  emxArray_real_T *r3;
  emxArray_real_T b_yTemp_data;
  int s_size[1];
  emxArray_boolean_T *b_x;
  double s_data[1026];
  emxArray_real_T b_s_data;
  emxArray_int32_T *b_ii;
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

  idx = 0;
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 1027)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      ii_data[idx - 1] = (short)ii;
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

  if (1 > idx) {
    nx = 0;
  } else {
    nx = idx;
  }

  tmp_data[0] = 1;
  for (i10 = 0; i10 < nx; i10++) {
    tmp_data[i10 + 1] = (short)(ii_data[i10] + 1);
  }

  if (1 > idx) {
    i11 = 0;
  } else {
    i11 = idx;
  }

  ii = 1 + i11;
  for (i10 = 0; i10 < ii; i10++) {
    iTemp_data[i10] = tmp_data[i10];
  }

  yTemp_size[0] = 1 + nx;
  nx++;
  for (i10 = 0; i10 < nx; i10++) {
    yTemp_data[i10] = b_yTemp[iTemp_data[i10] - 1];
  }

  emxInit_real_T1(&r3, 1);
  b_yTemp_data.data = (double *)&yTemp_data;
  b_yTemp_data.size = (int *)&yTemp_size;
  b_yTemp_data.allocatedSize = 1027;
  b_yTemp_data.numDimensions = 1;
  b_yTemp_data.canFreeData = false;
  diff(&b_yTemp_data, r3);
  b_sign(r3);
  s_size[0] = r3->size[0];
  nx = r3->size[0];
  for (i10 = 0; i10 < nx; i10++) {
    s_data[i10] = r3->data[i10];
  }

  emxInit_boolean_T(&b_x, 1);
  b_s_data.data = (double *)&s_data;
  b_s_data.size = (int *)&s_size;
  b_s_data.allocatedSize = 1026;
  b_s_data.numDimensions = 1;
  b_s_data.canFreeData = false;
  diff(&b_s_data, r3);
  i10 = b_x->size[0];
  b_x->size[0] = r3->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i10, (int)sizeof(boolean_T));
  nx = r3->size[0];
  for (i10 = 0; i10 < nx; i10++) {
    b_x->data[i10] = (r3->data[i10] < 0.0);
  }

  emxFree_real_T(&r3);
  emxInit_int32_T(&b_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i10 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i10, (int)sizeof(int));
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
      i10 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i10, (int)sizeof(int));
    }
  } else {
    i10 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i10, (int)sizeof(int));
  }

  if (1 > s_size[0] - 1) {
    nx = 0;
  } else {
    nx = s_size[0] - 1;
  }

  if (2 > s_size[0]) {
    i10 = 0;
  } else {
    i10 = 1;
  }

  ii = b_x->size[0];
  b_x->size[0] = nx;
  emxEnsureCapacity((emxArray__common *)b_x, ii, (int)sizeof(boolean_T));
  for (ii = 0; ii < nx; ii++) {
    b_x->data[ii] = (s_data[ii] != s_data[i10 + ii]);
  }

  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i10 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i10, (int)sizeof(int));
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
      i10 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i10, (int)sizeof(int));
    }
  } else {
    i10 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i10, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  iInflect_size[0] = c_ii->size[0];
  nx = c_ii->size[0];
  for (i10 = 0; i10 < nx; i10++) {
    iInflect_data[i10] = (double)iTemp_data[c_ii->data[i10]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  iPk_size[0] = b_ii->size[0];
  nx = b_ii->size[0];
  for (i10 = 0; i10 < nx; i10++) {
    iPk_data[i10] = (double)iTemp_data[b_ii->data[i10]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
}

//
// Arguments    : const double Yin[1025]
//                double Ypk_data[]
//                int Ypk_size[2]
//                double Xpk_data[]
//                int Xpk_size[2]
// Return Type  : void
//
static void b_findpeaks(const double Yin[1025], double Ypk_data[], int Ypk_size
  [2], double Xpk_data[], int Xpk_size[2])
{
  double iFinite_data[1025];
  int iFinite_size[1];
  double iInfite_data[1025];
  int iInfite_size[1];
  double iInflect_data[1025];
  int iInflect_size[1];
  int iPk_size[1];
  int loop_ub;
  int i9;
  double iPk_data[2050];
  emxArray_real_T *iPkOut;
  emxArray_int32_T *ia;
  emxArray_int32_T *ib;
  emxArray_real_T b_iPk_data;
  emxArray_real_T b_iInfite_data;
  double iPkOut_data[2050];
  int iPkOut_size[1];
  b_getAllPeaks(Yin, iFinite_data, iFinite_size, iInfite_data, iInfite_size,
                iInflect_data, iInflect_size);
  b_removePeaksBelowMinPeakHeight(Yin, iFinite_data, iFinite_size);
  iPk_size[0] = iFinite_size[0];
  loop_ub = iFinite_size[0];
  for (i9 = 0; i9 < loop_ub; i9++) {
    iPk_data[i9] = iFinite_data[i9];
  }

  iFinite_size[0] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i9 = 0; i9 < loop_ub; i9++) {
    iFinite_data[i9] = iPk_data[i9];
  }

  b_removePeaksBelowThreshold(Yin, iFinite_data, iFinite_size);
  iPk_size[0] = iFinite_size[0];
  loop_ub = iFinite_size[0];
  for (i9 = 0; i9 < loop_ub; i9++) {
    iPk_data[i9] = iFinite_data[i9];
  }

  emxInit_real_T1(&iPkOut, 1);
  emxInit_int32_T(&ia, 1);
  emxInit_int32_T(&ib, 1);
  b_iPk_data.data = (double *)&iPk_data;
  b_iPk_data.size = (int *)&iPk_size;
  b_iPk_data.allocatedSize = 2050;
  b_iPk_data.numDimensions = 1;
  b_iPk_data.canFreeData = false;
  b_iInfite_data.data = (double *)&iInfite_data;
  b_iInfite_data.size = (int *)&iInfite_size;
  b_iInfite_data.allocatedSize = 1025;
  b_iInfite_data.numDimensions = 1;
  b_iInfite_data.canFreeData = false;
  do_vectors(&b_iPk_data, &b_iInfite_data, iPkOut, ia, ib);
  d_findPeaksSeparatedByMoreThanM(iPkOut->size, iPk_data, iPk_size);
  orderPeaks(Yin, iPkOut->data, iPk_data, iPk_size);
  b_keepAtMostNpPeaks(iPk_data, iPk_size);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&ia);
  for (i9 = 0; i9 < 1025; i9++) {
    iFinite_data[i9] = 1.0 + (double)i9;
  }

  iPkOut_size[0] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i9 = 0; i9 < loop_ub; i9++) {
    iPkOut_data[i9] = iPkOut->data[(int)iPk_data[i9] - 1];
  }

  emxFree_real_T(&iPkOut);
  b_assignOutputs(Yin, iFinite_data, iPkOut_data, iPkOut_size, Ypk_data,
                  Ypk_size, Xpk_data, Xpk_size);
}

//
// Arguments    : const double y[1025]
//                double iPk_data[]
//                int iPk_size[1]
//                double iInf_data[]
//                int iInf_size[1]
//                double iInflect_data[]
//                int iInflect_size[1]
// Return Type  : void
//
static void b_getAllPeaks(const double y[1025], double iPk_data[], int iPk_size
  [1], double iInf_data[], int iInf_size[1], double iInflect_data[], int
  iInflect_size[1])
{
  boolean_T x[1025];
  int idx;
  short ii_data[1025];
  int ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  double yTemp[1025];
  for (idx = 0; idx < 1025; idx++) {
    x[idx] = (rtIsInf(y[idx]) && (y[idx] > 0.0));
  }

  idx = 0;
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii < 1026)) {
    guard1 = false;
    if (x[ii - 1]) {
      idx++;
      ii_data[idx - 1] = (short)ii;
      if (idx >= 1025) {
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

  if (1 > idx) {
    idx = 0;
  }

  iInf_size[0] = idx;
  for (ii = 0; ii < idx; ii++) {
    iInf_data[ii] = ii_data[ii];
  }

  memcpy(&yTemp[0], &y[0], 1025U * sizeof(double));
  for (ii = 0; ii < idx; ii++) {
    ii_data[ii] = (short)iInf_data[ii];
  }

  for (ii = 0; ii < idx; ii++) {
    yTemp[ii_data[ii] - 1] = rtNaN;
  }

  b_findLocalMaxima(yTemp, iPk_data, iPk_size, iInflect_data, iInflect_size);
}

//
// Arguments    : double idx_data[]
//                int idx_size[1]
// Return Type  : void
//
static void b_keepAtMostNpPeaks(double idx_data[], int idx_size[1])
{
  double b_idx_data[2050];
  int i27;
  if (idx_size[0] > 1025) {
    memcpy(&b_idx_data[0], &idx_data[0], 1025U * sizeof(double));
    idx_size[0] = 1025;
    for (i27 = 0; i27 < 1025; i27++) {
      idx_data[i27] = b_idx_data[i27];
    }
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
// Arguments    : const double tsX[42]
//                const emxArray_real_T *tX
//                const emxArray_real_T *tY
// Return Type  : double
//
static double b_knn(const double tsX[42], const emxArray_real_T *tX, const
                    emxArray_real_T *tY)
{
  double yfit;
  emxArray_int32_T *idx;
  int na;
  unsigned int outsize_idx_0;
  int pEnd;
  int i2;
  emxArray_int32_T *iwork;
  int n;
  emxArray_real_T *Uc;
  int k;
  boolean_T p;
  int i;
  int b_p;
  int j;
  int q;
  int qEnd;
  int kEnd;
  double x;
  int exitg3;
  double absxk;
  emxArray_real_T *a;
  int exponent;
  int i21;
  emxArray_real_T *y;
  int b_k;
  int c_k;
  emxArray_real_T *dist;
  emxArray_int32_T *pos;
  emxArray_real_T *edges;
  int exitg2;
  int exitg1;
  emxArray_real_T *nn;
  int b_exponent;
  emxInit_int32_T(&idx, 1);
  na = tY->size[0];
  outsize_idx_0 = (unsigned int)tY->size[0];
  pEnd = idx->size[0];
  idx->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, pEnd, (int)sizeof(int));
  i2 = (int)outsize_idx_0;
  for (pEnd = 0; pEnd < i2; pEnd++) {
    idx->data[pEnd] = 0;
  }

  if (tY->size[0] == 0) {
  } else {
    emxInit_int32_T(&iwork, 1);
    n = tY->size[0] + 1;
    pEnd = iwork->size[0];
    iwork->size[0] = (int)outsize_idx_0;
    emxEnsureCapacity((emxArray__common *)iwork, pEnd, (int)sizeof(int));
    for (k = 1; k <= n - 2; k += 2) {
      if ((tY->data[k - 1] <= tY->data[k]) || rtIsNaN(tY->data[k])) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((tY->size[0] & 1) != 0) {
      idx->data[tY->size[0] - 1] = tY->size[0];
    }

    i = 2;
    while (i < n - 1) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < n; pEnd = qEnd + i) {
        b_p = j;
        q = pEnd - 1;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if ((tY->data[idx->data[b_p - 1] - 1] <= tY->data[idx->data[q] - 1]) ||
              rtIsNaN(tY->data[idx->data[q] - 1])) {
            p = true;
          } else {
            p = false;
          }

          if (p) {
            iwork->data[k] = idx->data[b_p - 1];
            b_p++;
            if (b_p == pEnd) {
              while (q + 1 < qEnd) {
                k++;
                iwork->data[k] = idx->data[q];
                q++;
              }
            }
          } else {
            iwork->data[k] = idx->data[q];
            q++;
            if (q + 1 == qEnd) {
              while (b_p < pEnd) {
                k++;
                iwork->data[k] = idx->data[b_p - 1];
                b_p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&iwork);
  }

  emxInit_real_T1(&Uc, 1);
  outsize_idx_0 = (unsigned int)tY->size[0];
  pEnd = Uc->size[0];
  Uc->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)Uc, pEnd, (int)sizeof(double));
  for (k = 0; k + 1 <= na; k++) {
    Uc->data[k] = tY->data[idx->data[k] - 1];
  }

  emxFree_int32_T(&idx);
  k = 0;
  while ((k + 1 <= na) && rtIsInf(Uc->data[k]) && (Uc->data[k] < 0.0)) {
    k++;
  }

  b_p = k;
  k = tY->size[0];
  while ((k >= 1) && rtIsNaN(Uc->data[k - 1])) {
    k--;
  }

  pEnd = tY->size[0] - k;
  while ((k >= 1) && rtIsInf(Uc->data[k - 1]) && (Uc->data[k - 1] > 0.0)) {
    k--;
  }

  i = (tY->size[0] - k) - pEnd;
  q = -1;
  if (b_p > 0) {
    q = 0;
  }

  i2 = (b_p + k) - b_p;
  while (b_p + 1 <= i2) {
    x = Uc->data[b_p];
    do {
      exitg3 = 0;
      b_p++;
      if (b_p + 1 > i2) {
        exitg3 = 1;
      } else {
        absxk = std::abs(x / 2.0);
        if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
          if (absxk <= 2.2250738585072014E-308) {
            absxk = 4.94065645841247E-324;
          } else {
            frexp(absxk, &exponent);
            absxk = std::ldexp(1.0, exponent - 53);
          }
        } else {
          absxk = rtNaN;
        }

        if ((std::abs(x - Uc->data[b_p]) < absxk) || (rtIsInf(Uc->data[b_p]) &&
             rtIsInf(x) && ((Uc->data[b_p] > 0.0) == (x > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (!p) {
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);

    q++;
    Uc->data[q] = x;
  }

  if (i > 0) {
    q++;
    Uc->data[q] = Uc->data[i2];
  }

  b_p = i2 + i;
  for (j = 1; j <= pEnd; j++) {
    q++;
    Uc->data[q] = Uc->data[(b_p + j) - 1];
  }

  emxInit_real_T(&a, 2);
  pEnd = Uc->size[0];
  if (1 > q + 1) {
    i21 = -1;
  } else {
    i21 = q;
  }

  Uc->size[0] = i21 + 1;
  emxEnsureCapacity((emxArray__common *)Uc, pEnd, (int)sizeof(double));
  i = tY->size[0];
  pEnd = a->size[0] * a->size[1];
  a->size[0] = i;
  a->size[1] = 42;
  emxEnsureCapacity((emxArray__common *)a, pEnd, (int)sizeof(double));
  for (pEnd = 0; pEnd < i; pEnd++) {
    for (i2 = 0; i2 < 42; i2++) {
      a->data[pEnd + a->size[0] * i2] = tX->data[pEnd + tX->size[0] * i2] -
        tsX[i2];
    }
  }

  emxInit_real_T(&y, 2);
  pEnd = y->size[0] * y->size[1];
  y->size[0] = a->size[0];
  y->size[1] = 42;
  emxEnsureCapacity((emxArray__common *)y, pEnd, (int)sizeof(double));
  i = a->size[0] * 42;

//#pragma omp parallel for \
// num_threads(omp_get_max_threads()) \
// private(c_k)

  for (b_k = 1; b_k <= i; b_k++) {
    c_k = b_k;
    y->data[c_k - 1] = a->data[c_k - 1] * a->data[c_k - 1];
  }

  emxFree_real_T(&a);
  emxInit_real_T1(&dist, 1);
  pEnd = dist->size[0];
  dist->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)dist, pEnd, (int)sizeof(double));
  i = y->size[0];
  for (j = 0; j + 1 <= i; j++) {
    absxk = y->data[j];
    for (k = 0; k < 41; k++) {
      absxk += y->data[j + (k + 1) * i];
    }

    dist->data[j] = absxk;
  }

  emxFree_real_T(&y);
  emxInit_int32_T(&pos, 1);
  sort(dist, pos);
  pEnd = dist->size[0];
  dist->size[0] = pos->size[0];
  emxEnsureCapacity((emxArray__common *)dist, pEnd, (int)sizeof(double));
  i2 = pos->size[0];
  for (pEnd = 0; pEnd < i2; pEnd++) {
    dist->data[pEnd] = pos->data[pEnd];
  }

  emxFree_int32_T(&pos);
  emxInit_real_T(&edges, 2);
  i2 = Uc->size[0];
  pEnd = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = (int)(i2 + 1U);
  emxEnsureCapacity((emxArray__common *)edges, pEnd, (int)sizeof(double));
  k = 0;
  do {
    exitg2 = 0;
    i2 = Uc->size[0];
    if (k <= i2 - 2) {
      edges->data[1 + k] = Uc->data[k] + (Uc->data[1 + k] - Uc->data[k]) / 2.0;
      k++;
    } else {
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  edges->data[0] = rtMinusInf;
  edges->data[edges->size[1] - 1] = rtInf;
  k = 1;
  do {
    exitg1 = 0;
    i2 = Uc->size[0];
    if (k - 1 <= i2 - 2) {
      absxk = std::abs(edges->data[k]);
      if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
        if (absxk <= 2.2250738585072014E-308) {
          absxk = 4.94065645841247E-324;
        } else {
          frexp(absxk, &b_exponent);
          absxk = std::ldexp(1.0, b_exponent - 53);
        }
      } else {
        absxk = rtNaN;
      }

      edges->data[k] += absxk;
      k++;
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxInit_real_T1(&nn, 1);
  outsize_idx_0 = (unsigned int)edges->size[1];
  pEnd = nn->size[0];
  nn->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)nn, pEnd, (int)sizeof(double));
  i2 = (int)outsize_idx_0;
  for (pEnd = 0; pEnd < i2; pEnd++) {
    nn->data[pEnd] = 0.0;
  }

  i = 0;
  for (k = 0; k < 5; k++) {
    i2 = findbin(tY->data[(int)dist->data[i] - 1], edges);
    if (i2 > 0) {
      nn->data[i2 - 1]++;
    }

    i++;
  }

  emxFree_real_T(&dist);
  pEnd = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = nn->size[0] - 1;
  emxEnsureCapacity((emxArray__common *)edges, pEnd, (int)sizeof(double));
  for (k = 0; k <= nn->size[0] - 2; k++) {
    edges->data[k] = nn->data[k];
  }

  if (nn->size[0] - 1 > 0) {
    edges->data[edges->size[1] - 1] += nn->data[nn->size[0] - 1];
  }

  emxFree_real_T(&nn);
  n = edges->size[1];
  absxk = edges->data[0];
  i = 0;
  if (edges->size[1] > 1) {
    for (i2 = 1; i2 + 1 <= n; i2++) {
      if (edges->data[i2] > absxk) {
        absxk = edges->data[i2];
        i = i2;
      }
    }
  }

  emxFree_real_T(&edges);
  yfit = Uc->data[i];
  emxFree_real_T(&Uc);
  return yfit;
}

//
// Arguments    : const double x[4]
// Return Type  : double
//
static double b_mean(const double x[4])
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
  if (nq == 0) {
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
          n = (iout - p) + 1;
          while (p + 1 <= np) {
            idx->data[n + p] = iwork->data[p];
            x->data[n + p] = xwork->data[p];
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
// Return Type  : double
//
static double b_mod(double x)
{
  return x - std::floor(x / 2.0) * 2.0;
}

//
// Arguments    : const emxArray_real_T *x
//                int n1_unsigned
//                const emxArray_real_T *costab
//                const emxArray_real_T *sintab
//                emxArray_creal_T *y
// Return Type  : void
//
static void b_r2br_r2dit_trig(const emxArray_real_T *x, int n1_unsigned, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T *y)
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
  if (x->size[0] <= n1_unsigned) {
    j = x->size[0];
  } else {
    j = n1_unsigned;
  }

  nRowsD2 = n1_unsigned / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = n1_unsigned;
  emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
  if (n1_unsigned > x->size[0]) {
    iDelta = y->size[0];
    iy = y->size[0];
    y->size[0] = iDelta;
    emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
    for (iy = 0; iy < iDelta; iy++) {
      y->data[iy].re = 0.0;
      y->data[iy].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy].re = x->data[ix];
    y->data[iy].im = 0.0;
    iDelta = n1_unsigned;
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy].re = x->data[ix];
  y->data[iy].im = 0.0;
  if (n1_unsigned > 1) {
    for (i = 0; i <= n1_unsigned - 2; i += 2) {
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
// Arguments    : const double Y[1025]
//                double iPk_data[]
//                int iPk_size[1]
// Return Type  : void
//
static void b_removePeaksBelowMinPeakHeight(const double Y[1025], double
  iPk_data[], int iPk_size[1])
{
  int end;
  int trueCount;
  int i;
  int partialTrueCount;
  if (!(iPk_size[0] == 0)) {
    end = iPk_size[0] - 1;
    trueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y[(int)iPk_data[i] - 1] > rtMinusInf) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y[(int)iPk_data[i] - 1] > rtMinusInf) {
        iPk_data[partialTrueCount] = iPk_data[i];
        partialTrueCount++;
      }
    }

    iPk_size[0] = trueCount;
  }
}

//
// Arguments    : const double Y[1025]
//                double iPk_data[]
//                int iPk_size[1]
// Return Type  : void
//
static void b_removePeaksBelowThreshold(const double Y[1025], double iPk_data[],
  int iPk_size[1])
{
  emxArray_real_T *maxval;
  short csz_idx_0;
  int k;
  int trueCount;
  int i;
  int partialTrueCount;
  emxInit_real_T1(&maxval, 1);
  csz_idx_0 = (short)iPk_size[0];
  k = maxval->size[0];
  maxval->size[0] = (short)iPk_size[0];
  emxEnsureCapacity((emxArray__common *)maxval, k, (int)sizeof(double));
  for (k = 0; k + 1 <= csz_idx_0; k++) {
    if ((Y[(int)(iPk_data[k] - 1.0) - 1] >= Y[(int)(iPk_data[k] + 1.0) - 1]) ||
        rtIsNaN(Y[(int)(iPk_data[k] + 1.0) - 1])) {
      maxval->data[k] = Y[(int)(iPk_data[k] - 1.0) - 1];
    } else {
      maxval->data[k] = Y[(int)(iPk_data[k] + 1.0) - 1];
    }
  }

  k = iPk_size[0] - 1;
  trueCount = 0;
  for (i = 0; i <= k; i++) {
    if (Y[(int)iPk_data[i] - 1] - maxval->data[i] >= 0.0) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (i = 0; i <= k; i++) {
    if (Y[(int)iPk_data[i] - 1] - maxval->data[i] >= 0.0) {
      iPk_data[partialTrueCount] = iPk_data[i];
      partialTrueCount++;
    }
  }

  emxFree_real_T(&maxval);
  iPk_size[0] = trueCount;
}

//
// Arguments    : double a
//                double varargin_1
//                emxArray_real_T *b
// Return Type  : void
//
static void b_repmat(double a, double varargin_1, emxArray_real_T *b)
{
  int i15;
  int loop_ub;
  i15 = b->size[0];
  b->size[0] = (int)varargin_1;
  emxEnsureCapacity((emxArray__common *)b, i15, (int)sizeof(double));
  loop_ub = (int)varargin_1;
  for (i15 = 0; i15 < loop_ub; i15++) {
    b->data[i15] = a;
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
//                emxArray_int32_T *idx
// Return Type  : void
//
static void b_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int i25;
  i25 = nonSingletonDim(x);
  c_sort(x, i25, idx);
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
static void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_real_T *b_x;
  int ib;
  int nNonNaN;
  int m;
  int n;
  double x4[4];
  int idx4[4];
  emxArray_real_T *xwork;
  emxArray_int32_T *iwork;
  emxArray_real_T *b_xwork;
  int nNaNs;
  int k;
  int wOffset;
  signed char perm[4];
  emxArray_int32_T *b_iwork;
  int i4;
  emxInit_real_T1(&b_x, 1);
  ib = x->size[0];
  nNonNaN = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, nNonNaN, (int)sizeof(double));
  m = x->size[0];
  for (nNonNaN = 0; nNonNaN < m; nNonNaN++) {
    b_x->data[nNonNaN] = x->data[nNonNaN];
  }

  nNonNaN = idx->size[0];
  idx->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)idx, nNonNaN, (int)sizeof(int));
  for (nNonNaN = 0; nNonNaN < ib; nNonNaN++) {
    idx->data[nNonNaN] = 0;
  }

  n = x->size[0];
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  emxInit_real_T1(&xwork, 1);
  emxInit_int32_T(&iwork, 1);
  emxInit_real_T1(&b_xwork, 1);
  nNonNaN = iwork->size[0];
  iwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)iwork, nNonNaN, (int)sizeof(int));
  ib = x->size[0];
  nNonNaN = b_xwork->size[0];
  b_xwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)b_xwork, nNonNaN, (int)sizeof(double));
  nNonNaN = xwork->size[0];
  xwork->size[0] = b_xwork->size[0];
  emxEnsureCapacity((emxArray__common *)xwork, nNonNaN, (int)sizeof(double));
  m = b_xwork->size[0];
  emxFree_real_T(&b_xwork);
  for (nNonNaN = 0; nNonNaN < m; nNonNaN++) {
    xwork->data[nNonNaN] = 0.0;
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
          nNonNaN = 3;
          i4 = 4;
        } else {
          nNonNaN = 4;
          i4 = 3;
        }

        if (x4[m - 1] >= x4[nNonNaN - 1]) {
          if (x4[wOffset - 1] >= x4[nNonNaN - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)nNonNaN;
            perm[3] = (signed char)i4;
          } else if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)nNonNaN;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)m;
            perm[1] = (signed char)nNonNaN;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else if (x4[m - 1] >= x4[i4 - 1]) {
          if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)nNonNaN;
            perm[1] = (signed char)m;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)nNonNaN;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else {
          perm[0] = (signed char)nNonNaN;
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
    nNonNaN = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[n - k];
    idx->data[n - k] = nNonNaN;
    b_x->data[wOffset + k] = xwork->data[n - k];
    b_x->data[n - k] = xwork->data[wOffset + k];
  }

  if ((nNaNs & 1) != 0) {
    b_x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
  }

  emxInit_int32_T(&b_iwork, 1);
  nNonNaN = b_iwork->size[0];
  b_iwork->size[0] = iwork->size[0];
  emxEnsureCapacity((emxArray__common *)b_iwork, nNonNaN, (int)sizeof(int));
  m = iwork->size[0];
  emxFree_int32_T(&iwork);
  for (nNonNaN = 0; nNonNaN < m; nNonNaN++) {
    b_iwork->data[nNonNaN] = 0;
  }

  nNonNaN = x->size[0] - nNaNs;
  ib = 2;
  if (nNonNaN > 1) {
    if (x->size[0] >= 256) {
      m = nNonNaN >> 8;
      if (m > 0) {
        for (ib = 1; ib <= m; ib++) {
          merge_pow2_block(idx, b_x, (ib - 1) << 8);
        }

        ib = m << 8;
        m = nNonNaN - ib;
        if (m > 0) {
          b_merge_block(idx, b_x, ib, m, 2, b_iwork, xwork);
        }

        ib = 8;
      }
    }

    b_merge_block(idx, b_x, 0, nNonNaN, ib, b_iwork, xwork);
  }

  if ((nNaNs > 0) && (nNonNaN > 0)) {
    for (k = 0; k + 1 <= nNaNs; k++) {
      xwork->data[k] = b_x->data[nNonNaN + k];
      b_iwork->data[k] = idx->data[nNonNaN + k];
    }

    for (k = nNonNaN - 1; k + 1 > 0; k--) {
      b_x->data[nNaNs + k] = b_x->data[k];
      idx->data[nNaNs + k] = idx->data[k];
    }

    for (k = 0; k + 1 <= nNaNs; k++) {
      b_x->data[k] = xwork->data[k];
      idx->data[k] = b_iwork->data[k];
    }
  }

  emxFree_real_T(&xwork);
  emxFree_int32_T(&b_iwork);
  nNonNaN = x->size[0];
  x->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)x, nNonNaN, (int)sizeof(double));
  m = b_x->size[0];
  for (nNonNaN = 0; nNonNaN < m; nNonNaN++) {
    x->data[nNonNaN] = b_x->data[nNonNaN];
  }

  emxFree_real_T(&b_x);
}

//
// Arguments    : const emxArray_real_T *x
// Return Type  : double
//
static double b_sum(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[0]; k++) {
      y += x->data[k - 1];
    }
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
// Arguments    : const emxArray_real_T *y
//                const emxArray_real_T *x
//                const emxArray_real_T *iPk
//                emxArray_real_T *YpkOut
//                emxArray_real_T *XpkOut
// Return Type  : void
//
static void c_assignOutputs(const emxArray_real_T *y, const emxArray_real_T *x,
  const emxArray_real_T *iPk, emxArray_real_T *YpkOut, emxArray_real_T *XpkOut)
{
  int i20;
  int loop_ub;
  i20 = YpkOut->size[0] * YpkOut->size[1];
  YpkOut->size[0] = 1;
  YpkOut->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)YpkOut, i20, (int)sizeof(double));
  loop_ub = iPk->size[0];
  for (i20 = 0; i20 < loop_ub; i20++) {
    YpkOut->data[YpkOut->size[0] * i20] = y->data[(int)iPk->data[i20] - 1];
  }

  i20 = XpkOut->size[0] * XpkOut->size[1];
  XpkOut->size[0] = 1;
  XpkOut->size[1] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)XpkOut, i20, (int)sizeof(double));
  loop_ub = iPk->size[0];
  for (i20 = 0; i20 < loop_ub; i20++) {
    XpkOut->data[XpkOut->size[0] * i20] = x->data[(int)iPk->data[i20] - 1];
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
  emxArray_real_T *r8;
  int i19;
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
  emxArray_int32_T *r9;
  emxArray_real_T *b_iTemp;
  emxArray_real_T *b_yTemp;
  emxArray_real_T *s;
  emxArray_real_T *r10;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *b_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxInit_real_T1(&r8, 1);
  i19 = r8->size[0];
  r8->size[0] = 2 + yTemp->size[0];
  emxEnsureCapacity((emxArray__common *)r8, i19, (int)sizeof(double));
  r8->data[0] = rtNaN;
  cdiff = yTemp->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    r8->data[i19 + 1] = yTemp->data[i19];
  }

  r8->data[1 + yTemp->size[0]] = rtNaN;
  i19 = yTemp->size[0];
  yTemp->size[0] = r8->size[0];
  emxEnsureCapacity((emxArray__common *)yTemp, i19, (int)sizeof(double));
  cdiff = r8->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    yTemp->data[i19] = r8->data[i19];
  }

  emxFree_real_T(&r8);
  ndbl = (int)std::floor(((double)yTemp->size[0] - 1.0) + 0.5);
  apnd = ndbl + 1;
  cdiff = (ndbl - yTemp->size[0]) + 1;
  absb = yTemp->size[0];
  if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
    ndbl++;
    apnd = yTemp->size[0];
  } else if (cdiff > 0) {
    apnd = ndbl;
  } else {
    ndbl++;
  }

  emxInit_real_T(&y, 2);
  i19 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, i19, (int)sizeof(double));
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
  i19 = iTemp->size[0];
  iTemp->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)iTemp, i19, (int)sizeof(double));
  cdiff = y->size[1];
  for (i19 = 0; i19 < cdiff; i19++) {
    iTemp->data[i19] = y->data[y->size[0] * i19];
  }

  emxFree_real_T(&y);
  emxInit_boolean_T(&yFinite, 1);
  i19 = yFinite->size[0];
  yFinite->size[0] = yTemp->size[0];
  emxEnsureCapacity((emxArray__common *)yFinite, i19, (int)sizeof(boolean_T));
  cdiff = yTemp->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    yFinite->data[i19] = rtIsNaN(yTemp->data[i19]);
  }

  i19 = yFinite->size[0];
  emxEnsureCapacity((emxArray__common *)yFinite, i19, (int)sizeof(boolean_T));
  cdiff = yFinite->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    yFinite->data[i19] = !yFinite->data[i19];
  }

  emxInit_boolean_T(&x, 1);
  cdiff = yTemp->size[0] - 2;
  i19 = x->size[0];
  x->size[0] = cdiff + 1;
  emxEnsureCapacity((emxArray__common *)x, i19, (int)sizeof(boolean_T));
  for (i19 = 0; i19 <= cdiff; i19++) {
    x->data[i19] = ((yTemp->data[i19] != yTemp->data[1 + i19]) && (yFinite->
      data[i19] || yFinite->data[1 + i19]));
  }

  emxFree_boolean_T(&yFinite);
  emxInit_int32_T(&ii, 1);
  absb = x->size[0];
  ndbl = 0;
  i19 = ii->size[0];
  ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)ii, i19, (int)sizeof(int));
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
      i19 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i19, (int)sizeof(int));
    }
  } else {
    i19 = ii->size[0];
    if (1 > ndbl) {
      ii->size[0] = 0;
    } else {
      ii->size[0] = ndbl;
    }

    emxEnsureCapacity((emxArray__common *)ii, i19, (int)sizeof(int));
  }

  emxInit_int32_T(&r9, 1);
  i19 = r9->size[0];
  r9->size[0] = 1 + ii->size[0];
  emxEnsureCapacity((emxArray__common *)r9, i19, (int)sizeof(int));
  r9->data[0] = 1;
  cdiff = ii->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    r9->data[i19 + 1] = ii->data[i19] + 1;
  }

  emxInit_real_T1(&b_iTemp, 1);
  i19 = b_iTemp->size[0];
  b_iTemp->size[0] = r9->size[0];
  emxEnsureCapacity((emxArray__common *)b_iTemp, i19, (int)sizeof(double));
  cdiff = r9->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    b_iTemp->data[i19] = iTemp->data[r9->data[i19] - 1];
  }

  emxFree_int32_T(&r9);
  i19 = iTemp->size[0];
  iTemp->size[0] = b_iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i19, (int)sizeof(double));
  cdiff = b_iTemp->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    iTemp->data[i19] = b_iTemp->data[i19];
  }

  emxFree_real_T(&b_iTemp);
  emxInit_real_T1(&b_yTemp, 1);
  i19 = b_yTemp->size[0];
  b_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)b_yTemp, i19, (int)sizeof(double));
  cdiff = iTemp->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    b_yTemp->data[i19] = yTemp->data[(int)iTemp->data[i19] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_real_T1(&r10, 1);
  diff(b_yTemp, s);
  b_sign(s);
  diff(s, r10);
  i19 = x->size[0];
  x->size[0] = r10->size[0];
  emxEnsureCapacity((emxArray__common *)x, i19, (int)sizeof(boolean_T));
  cdiff = r10->size[0];
  emxFree_real_T(&b_yTemp);
  for (i19 = 0; i19 < cdiff; i19++) {
    x->data[i19] = (r10->data[i19] < 0.0);
  }

  emxFree_real_T(&r10);
  absb = x->size[0];
  ndbl = 0;
  i19 = ii->size[0];
  ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)ii, i19, (int)sizeof(int));
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
      i19 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i19, (int)sizeof(int));
    }
  } else {
    i19 = ii->size[0];
    if (1 > ndbl) {
      ii->size[0] = 0;
    } else {
      ii->size[0] = ndbl;
    }

    emxEnsureCapacity((emxArray__common *)ii, i19, (int)sizeof(int));
  }

  if (1 > s->size[0] - 1) {
    cdiff = 0;
  } else {
    cdiff = s->size[0] - 1;
  }

  if (2 > s->size[0]) {
    i19 = 0;
  } else {
    i19 = 1;
  }

  absb = x->size[0];
  x->size[0] = cdiff;
  emxEnsureCapacity((emxArray__common *)x, absb, (int)sizeof(boolean_T));
  for (absb = 0; absb < cdiff; absb++) {
    x->data[absb] = (s->data[absb] != s->data[i19 + absb]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&b_ii, 1);
  absb = x->size[0];
  ndbl = 0;
  i19 = b_ii->size[0];
  b_ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i19, (int)sizeof(int));
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
      i19 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i19, (int)sizeof(int));
    }
  } else {
    i19 = b_ii->size[0];
    if (1 > ndbl) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = ndbl;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i19, (int)sizeof(int));
  }

  emxFree_boolean_T(&x);
  i19 = iInflect->size[0];
  iInflect->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i19, (int)sizeof(double));
  cdiff = b_ii->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    iInflect->data[i19] = iTemp->data[b_ii->data[i19]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  i19 = iPk->size[0];
  iPk->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i19, (int)sizeof(double));
  cdiff = ii->size[0];
  for (i19 = 0; i19 < cdiff; i19++) {
    iPk->data[i19] = iTemp->data[ii->data[i19]] - 1.0;
  }

  emxFree_int32_T(&ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : const int iPk_size[1]
//                double idx_data[]
//                int idx_size[1]
// Return Type  : void
//
static void c_findPeaksSeparatedByMoreThanM(const int iPk_size[1], double
  idx_data[], int idx_size[1])
{
  int ndbl;
  int apnd;
  int cdiff;
  emxArray_real_T *y;
  int k;
  if (iPk_size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)iPk_size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - iPk_size[0]) + 1;
    if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)iPk_size[0])
    {
      ndbl++;
      apnd = iPk_size[0];
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
      cdiff = (ndbl - 1) / 2;
      for (k = 1; k < cdiff; k++) {
        y->data[k] = 1.0 + (double)k;
        y->data[(ndbl - k) - 1] = apnd - k;
      }

      if (cdiff << 1 == ndbl - 1) {
        y->data[cdiff] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[cdiff] = 1.0 + (double)cdiff;
        y->data[cdiff + 1] = apnd - cdiff;
      }
    }
  }

  idx_size[0] = y->size[1];
  cdiff = y->size[1];
  for (k = 0; k < cdiff; k++) {
    idx_data[k] = y->data[y->size[0] * k];
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
  emxArray_real_T *iPk;
  emxArray_real_T *y;
  emxArray_real_T *x;
  emxArray_real_T *iFinite;
  emxArray_real_T *iInflect;
  emxArray_real_T *iInf;
  double maxN;
  int i17;
  int loop_ub;
  emxArray_real_T *b_iPk;
  emxArray_int32_T *ia;
  emxArray_int32_T *ib;
  emxArray_real_T *b_x;
  emxArray_real_T *c_iPk;
  emxArray_int32_T *iidx;
  emxArray_real_T *d_iPk;
  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&y, 1);
  emxInit_real_T1(&x, 1);
  emxInit_real_T1(&iFinite, 1);
  emxInit_real_T1(&iInflect, 1);
  emxInit_real_T1(&iInf, 1);
  parse_inputs(Yin, y, x, &maxN);
  c_getAllPeaks(y, iFinite, iInf, iInflect);
  c_removePeaksBelowMinPeakHeight(y, iFinite);
  i17 = iPk->size[0];
  iPk->size[0] = iFinite->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i17, (int)sizeof(double));
  loop_ub = iFinite->size[0];
  emxFree_real_T(&iInflect);
  for (i17 = 0; i17 < loop_ub; i17++) {
    iPk->data[i17] = iFinite->data[i17];
  }

  i17 = iFinite->size[0];
  iFinite->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)iFinite, i17, (int)sizeof(double));
  loop_ub = iPk->size[0];
  for (i17 = 0; i17 < loop_ub; i17++) {
    iFinite->data[i17] = iPk->data[i17];
  }

  c_removePeaksBelowThreshold(y, iFinite);
  i17 = iPk->size[0];
  iPk->size[0] = iFinite->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i17, (int)sizeof(double));
  loop_ub = iFinite->size[0];
  for (i17 = 0; i17 < loop_ub; i17++) {
    iPk->data[i17] = iFinite->data[i17];
  }

  emxFree_real_T(&iFinite);
  emxInit_real_T1(&b_iPk, 1);
  emxInit_int32_T(&ia, 1);
  emxInit_int32_T(&ib, 1);
  do_vectors(iPk, iInf, b_iPk, ia, ib);
  e_findPeaksSeparatedByMoreThanM(b_iPk, iPk);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&ia);
  emxFree_real_T(&iInf);
  if (iPk->size[0] == 0) {
  } else {
    emxInit_real_T1(&b_x, 1);
    i17 = b_x->size[0];
    b_x->size[0] = iPk->size[0];
    emxEnsureCapacity((emxArray__common *)b_x, i17, (int)sizeof(double));
    loop_ub = iPk->size[0];
    for (i17 = 0; i17 < loop_ub; i17++) {
      b_x->data[i17] = y->data[(int)b_iPk->data[(int)iPk->data[i17] - 1] - 1];
    }

    emxInit_int32_T(&iidx, 1);
    emxInit_real_T1(&d_iPk, 1);
    b_sort(b_x, iidx);
    i17 = d_iPk->size[0];
    d_iPk->size[0] = iidx->size[0];
    emxEnsureCapacity((emxArray__common *)d_iPk, i17, (int)sizeof(double));
    loop_ub = iidx->size[0];
    emxFree_real_T(&b_x);
    for (i17 = 0; i17 < loop_ub; i17++) {
      d_iPk->data[i17] = iPk->data[iidx->data[i17] - 1];
    }

    emxFree_int32_T(&iidx);
    i17 = iPk->size[0];
    iPk->size[0] = d_iPk->size[0];
    emxEnsureCapacity((emxArray__common *)iPk, i17, (int)sizeof(double));
    loop_ub = d_iPk->size[0];
    for (i17 = 0; i17 < loop_ub; i17++) {
      iPk->data[i17] = d_iPk->data[i17];
    }

    emxFree_real_T(&d_iPk);
  }

  emxInit_real_T1(&c_iPk, 1);
  c_keepAtMostNpPeaks(iPk, maxN);
  i17 = c_iPk->size[0];
  c_iPk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)c_iPk, i17, (int)sizeof(double));
  loop_ub = iPk->size[0];
  for (i17 = 0; i17 < loop_ub; i17++) {
    c_iPk->data[i17] = b_iPk->data[(int)iPk->data[i17] - 1];
  }

  emxFree_real_T(&b_iPk);
  emxFree_real_T(&iPk);
  c_assignOutputs(y, x, c_iPk, Ypk, Xpk);
  emxFree_real_T(&c_iPk);
  emxFree_real_T(&x);
  emxFree_real_T(&y);
}

//
// Arguments    : const emxArray_real_T *y
//                emxArray_real_T *iPk
//                emxArray_real_T *iInf
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void c_getAllPeaks(const emxArray_real_T *y, emxArray_real_T *iPk,
  emxArray_real_T *iInf, emxArray_real_T *iInflect)
{
  emxArray_boolean_T *r6;
  int i18;
  int ii;
  emxArray_boolean_T *x;
  emxArray_int32_T *b_ii;
  int nx;
  int idx;
  boolean_T exitg1;
  boolean_T guard1 = false;
  emxArray_real_T *yTemp;
  emxArray_int32_T *r7;
  emxArray_real_T *b_yTemp;
  emxInit_boolean_T(&r6, 1);
  i18 = r6->size[0];
  r6->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)r6, i18, (int)sizeof(boolean_T));
  ii = y->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    r6->data[i18] = rtIsInf(y->data[i18]);
  }

  emxInit_boolean_T(&x, 1);
  i18 = x->size[0];
  x->size[0] = r6->size[0];
  emxEnsureCapacity((emxArray__common *)x, i18, (int)sizeof(boolean_T));
  ii = r6->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    x->data[i18] = (r6->data[i18] && (y->data[i18] > 0.0));
  }

  emxFree_boolean_T(&r6);
  emxInit_int32_T(&b_ii, 1);
  nx = x->size[0];
  idx = 0;
  i18 = b_ii->size[0];
  b_ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i18, (int)sizeof(int));
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
      i18 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i18, (int)sizeof(int));
    }
  } else {
    i18 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i18, (int)sizeof(int));
  }

  emxFree_boolean_T(&x);
  i18 = iInf->size[0];
  iInf->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInf, i18, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    iInf->data[i18] = b_ii->data[i18];
  }

  emxInit_real_T1(&yTemp, 1);
  i18 = yTemp->size[0];
  yTemp->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)yTemp, i18, (int)sizeof(double));
  ii = y->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    yTemp->data[i18] = y->data[i18];
  }

  emxInit_int32_T(&r7, 1);
  i18 = r7->size[0];
  r7->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r7, i18, (int)sizeof(int));
  ii = b_ii->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    r7->data[i18] = b_ii->data[i18];
  }

  emxFree_int32_T(&b_ii);
  ii = r7->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    yTemp->data[r7->data[i18] - 1] = rtNaN;
  }

  emxFree_int32_T(&r7);
  emxInit_real_T1(&b_yTemp, 1);
  i18 = b_yTemp->size[0];
  b_yTemp->size[0] = yTemp->size[0];
  emxEnsureCapacity((emxArray__common *)b_yTemp, i18, (int)sizeof(double));
  ii = yTemp->size[0];
  for (i18 = 0; i18 < ii; i18++) {
    b_yTemp->data[i18] = yTemp->data[i18];
  }

  emxFree_real_T(&yTemp);
  c_findLocalMaxima(b_yTemp, iPk, iInflect);
  emxFree_real_T(&b_yTemp);
}

//
// Arguments    : emxArray_real_T *idx
//                double Np
// Return Type  : void
//
static void c_keepAtMostNpPeaks(emxArray_real_T *idx, double Np)
{
  int loop_ub;
  emxArray_real_T *b_idx;
  int i28;
  if (idx->size[0] > Np) {
    if (1.0 > Np) {
      loop_ub = 0;
    } else {
      loop_ub = (int)Np;
    }

    emxInit_real_T1(&b_idx, 1);
    i28 = b_idx->size[0];
    b_idx->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_idx, i28, (int)sizeof(double));
    for (i28 = 0; i28 < loop_ub; i28++) {
      b_idx->data[i28] = idx->data[i28];
    }

    i28 = idx->size[0];
    idx->size[0] = b_idx->size[0];
    emxEnsureCapacity((emxArray__common *)idx, i28, (int)sizeof(double));
    loop_ub = b_idx->size[0];
    for (i28 = 0; i28 < loop_ub; i28++) {
      idx->data[i28] = b_idx->data[i28];
    }

    emxFree_real_T(&b_idx);
  }
}

//
// Arguments    : const emxArray_creal_T *x
//                int n1_unsigned
//                const emxArray_real_T *costab
//                const emxArray_real_T *sintab
//                emxArray_creal_T *y
// Return Type  : void
//
static void c_r2br_r2dit_trig(const emxArray_creal_T *x, int n1_unsigned, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T *y)
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
  if (x->size[0] <= n1_unsigned) {
    j = x->size[0];
  } else {
    j = n1_unsigned;
  }

  nRowsD2 = n1_unsigned / 2;
  nRowsD4 = nRowsD2 / 2;
  iy = y->size[0];
  y->size[0] = n1_unsigned;
  emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
  if (n1_unsigned > x->size[0]) {
    iDelta = y->size[0];
    iy = y->size[0];
    y->size[0] = iDelta;
    emxEnsureCapacity((emxArray__common *)y, iy, (int)sizeof(creal_T));
    for (iy = 0; iy < iDelta; iy++) {
      y->data[iy].re = 0.0;
      y->data[iy].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy] = x->data[ix];
    iDelta = n1_unsigned;
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
  if (n1_unsigned > 1) {
    for (i = 0; i <= n1_unsigned - 2; i += 2) {
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
// Return Type  : void
//
static void c_removePeaksBelowMinPeakHeight(const emxArray_real_T *Y,
  emxArray_real_T *iPk)
{
  int end;
  int trueCount;
  int i;
  int partialTrueCount;
  if (!(iPk->size[0] == 0)) {
    end = iPk->size[0] - 1;
    trueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y->data[(int)iPk->data[i] - 1] > rtMinusInf) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y->data[(int)iPk->data[i] - 1] > rtMinusInf) {
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
// Return Type  : void
//
static void c_removePeaksBelowThreshold(const emxArray_real_T *Y,
  emxArray_real_T *iPk)
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
    if (Y->data[(int)iPk->data[c] - 1] - base->data[c] >= 0.0) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (c = 0; c <= k; c++) {
    if (Y->data[(int)iPk->data[c] - 1] - base->data[c] >= 0.0) {
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
// Arguments    : emxArray_real_T *x
//                int dim
//                emxArray_int32_T *idx
// Return Type  : void
//
static void c_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
  int i26;
  emxArray_real_T *vwork;
  int vstride;
  int x_idx_0;
  int j;
  emxArray_int32_T *iidx;
  if (dim <= 1) {
    i26 = x->size[0];
  } else {
    i26 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  vstride = vwork->size[0];
  vwork->size[0] = i26;
  emxEnsureCapacity((emxArray__common *)vwork, vstride, (int)sizeof(double));
  x_idx_0 = x->size[0];
  vstride = idx->size[0];
  idx->size[0] = x_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, vstride, (int)sizeof(int));
  vstride = 1;
  x_idx_0 = 1;
  while (x_idx_0 <= dim - 1) {
    vstride *= x->size[0];
    x_idx_0 = 2;
  }

  j = 0;
  emxInit_int32_T(&iidx, 1);
  while (j + 1 <= vstride) {
    for (x_idx_0 = 0; x_idx_0 + 1 <= i26; x_idx_0++) {
      vwork->data[x_idx_0] = x->data[j + x_idx_0 * vstride];
    }

    b_sortIdx(vwork, iidx);
    for (x_idx_0 = 0; x_idx_0 + 1 <= i26; x_idx_0++) {
      x->data[j + x_idx_0 * vstride] = vwork->data[x_idx_0];
      idx->data[j + x_idx_0 * vstride] = iidx->data[x_idx_0];
    }

    j++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

//
// Arguments    : const int iPk_size[1]
//                double idx_data[]
//                int idx_size[1]
// Return Type  : void
//
static void d_findPeaksSeparatedByMoreThanM(const int iPk_size[1], double
  idx_data[], int idx_size[1])
{
  int ndbl;
  int apnd;
  int cdiff;
  emxArray_real_T *y;
  int k;
  if (iPk_size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)iPk_size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - iPk_size[0]) + 1;
    if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)iPk_size[0])
    {
      ndbl++;
      apnd = iPk_size[0];
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
      cdiff = (ndbl - 1) / 2;
      for (k = 1; k < cdiff; k++) {
        y->data[k] = 1.0 + (double)k;
        y->data[(ndbl - k) - 1] = apnd - k;
      }

      if (cdiff << 1 == ndbl - 1) {
        y->data[cdiff] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[cdiff] = 1.0 + (double)cdiff;
        y->data[cdiff + 1] = apnd - cdiff;
      }
    }
  }

  idx_size[0] = y->size[1];
  cdiff = y->size[1];
  for (k = 0; k < cdiff; k++) {
    idx_data[k] = y->data[y->size[0] * k];
  }

  emxFree_real_T(&y);
}

//
// Arguments    : const emxArray_creal_T *x
//                int n1_unsigned
//                const emxArray_real_T *costab
//                const emxArray_real_T *sintab
//                emxArray_creal_T *y
// Return Type  : void
//
static void d_r2br_r2dit_trig(const emxArray_creal_T *x, int n1_unsigned, const
  emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T *y)
{
  int j;
  int nRowsD2;
  int nRowsD4;
  int iDelta2;
  int iy;
  int ix;
  int ju;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  double r;
  double twid_im;
  int ihi;
  if (x->size[0] <= n1_unsigned) {
    j = x->size[0];
  } else {
    j = n1_unsigned;
  }

  nRowsD2 = n1_unsigned / 2;
  nRowsD4 = nRowsD2 / 2;
  iDelta2 = y->size[0];
  y->size[0] = n1_unsigned;
  emxEnsureCapacity((emxArray__common *)y, iDelta2, (int)sizeof(creal_T));
  if (n1_unsigned > x->size[0]) {
    iy = y->size[0];
    iDelta2 = y->size[0];
    y->size[0] = iy;
    emxEnsureCapacity((emxArray__common *)y, iDelta2, (int)sizeof(creal_T));
    for (iDelta2 = 0; iDelta2 < iy; iDelta2++) {
      y->data[iDelta2].re = 0.0;
      y->data[iDelta2].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < j; i++) {
    y->data[iy] = x->data[ix];
    iDelta2 = n1_unsigned;
    tst = true;
    while (tst) {
      iDelta2 >>= 1;
      ju ^= iDelta2;
      tst = ((ju & iDelta2) == 0);
    }

    iy = ju;
    ix++;
  }

  y->data[iy] = x->data[ix];
  if (n1_unsigned > 1) {
    for (i = 0; i <= n1_unsigned - 2; i += 2) {
      temp_re = y->data[i + 1].re;
      temp_im = y->data[i + 1].im;
      y->data[i + 1].re = y->data[i].re - y->data[i + 1].re;
      y->data[i + 1].im = y->data[i].im - y->data[i + 1].im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }
  }

  iy = 2;
  iDelta2 = 4;
  ix = 1 + ((nRowsD4 - 1) << 2);
  while (nRowsD4 > 0) {
    for (i = 0; i < ix; i += iDelta2) {
      temp_re = y->data[i + iy].re;
      temp_im = y->data[i + iy].im;
      y->data[i + iy].re = y->data[i].re - temp_re;
      y->data[i + iy].im = y->data[i].im - temp_im;
      y->data[i].re += temp_re;
      y->data[i].im += temp_im;
    }

    ju = 1;
    for (j = nRowsD4; j < nRowsD2; j += nRowsD4) {
      r = costab->data[j];
      twid_im = sintab->data[j];
      i = ju;
      ihi = ju + ix;
      while (i < ihi) {
        temp_re = r * y->data[i + iy].re - twid_im * y->data[i + iy].im;
        temp_im = r * y->data[i + iy].im + twid_im * y->data[i + iy].re;
        y->data[i + iy].re = y->data[i].re - temp_re;
        y->data[i + iy].im = y->data[i].im - temp_im;
        y->data[i].re += temp_re;
        y->data[i].im += temp_im;
        i += iDelta2;
      }

      ju++;
    }

    nRowsD4 /= 2;
    iy = iDelta2;
    iDelta2 <<= 1;
    ix -= iy;
  }

  if (y->size[0] > 1) {
    r = 1.0 / (double)y->size[0];
    iDelta2 = y->size[0];
    emxEnsureCapacity((emxArray__common *)y, iDelta2, (int)sizeof(creal_T));
    iy = y->size[0];
    for (iDelta2 = 0; iDelta2 < iy; iDelta2++) {
      y->data[iDelta2].re *= r;
      y->data[iDelta2].im *= r;
    }
  }
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
  double work_data_idx_0;
  int m;
  double tmp1;
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
      orderForDim = x->size[0] - 1;
      iyLead = y->size[0];
      y->size[0] = orderForDim;
      emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
      if (!(y->size[0] == 0)) {
        orderForDim = 1;
        iyLead = 0;
        work_data_idx_0 = x->data[0];
        for (m = 2; m <= x->size[0]; m++) {
          tmp1 = x->data[orderForDim];
          tmp2 = work_data_idx_0;
          work_data_idx_0 = tmp1;
          tmp1 -= tmp2;
          orderForDim++;
          y->data[iyLead] = tmp1;
          iyLead++;
        }
      }
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
  double absxk;
  emxArray_int32_T *b_ia;
  int exponent;
  emxArray_int32_T *b_ib;
  boolean_T p;
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
    absxk = std::abs(bk / 2.0);
    if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
      if (absxk <= 2.2250738585072014E-308) {
        absxk = 4.94065645841247E-324;
      } else {
        frexp(absxk, &exponent);
        absxk = std::ldexp(1.0, exponent - 53);
      }
    } else {
      absxk = rtNaN;
    }

    if ((std::abs(bk - ak) < absxk) || (rtIsInf(ak) && rtIsInf(bk) && ((ak > 0.0)
          == (bk > 0.0)))) {
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
  int k;
  int xidx;
  double wwc_re;
  double wwc_im;
  emxArray_creal_T *fy;
  emxArray_creal_T *fv;
  double fv_re;
  double fv_im;
  double b_wwc_re;
  double b_fv_im;
  double b_wwc_im;
  double b_fv_re;
  emxInit_creal_T(&wwc, 1);
  bluestein_setup(n1, wwc);
  if (n1 <= x->size[0]) {
    minNrowsNx = n1;
  } else {
    minNrowsNx = x->size[0];
  }

  k = y->size[0];
  y->size[0] = n1;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(creal_T));
  if (n1 > x->size[0]) {
    xidx = y->size[0];
    k = y->size[0];
    y->size[0] = xidx;
    emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(creal_T));
    for (k = 0; k < xidx; k++) {
      y->data[k].re = 0.0;
      y->data[k].im = 0.0;
    }
  }

  xidx = 0;
  for (k = 0; k + 1 <= minNrowsNx; k++) {
    wwc_re = wwc->data[(n1 + k) - 1].re;
    wwc_im = wwc->data[(n1 + k) - 1].im;
    y->data[k].re = wwc_re * x->data[xidx];
    y->data[k].im = wwc_im * -x->data[xidx];
    xidx++;
  }

  while (minNrowsNx + 1 <= n1) {
    y->data[minNrowsNx].re = 0.0;
    y->data[minNrowsNx].im = 0.0;
    minNrowsNx++;
  }

  emxInit_creal_T(&fy, 1);
  emxInit_creal_T(&fv, 1);
  r2br_r2dit_trig_impl(y, N2, costab, sintab, fy);
  c_r2br_r2dit_trig(wwc, N2, costab, sintab, fv);
  k = fy->size[0];
  emxEnsureCapacity((emxArray__common *)fy, k, (int)sizeof(creal_T));
  xidx = fy->size[0];
  for (k = 0; k < xidx; k++) {
    wwc_re = fy->data[k].re;
    wwc_im = fy->data[k].im;
    fv_re = fv->data[k].re;
    fv_im = fv->data[k].im;
    fy->data[k].re = wwc_re * fv_re - wwc_im * fv_im;
    fy->data[k].im = wwc_re * fv_im + wwc_im * fv_re;
  }

  d_r2br_r2dit_trig(fy, N2, costab, sintabinv, fv);
  xidx = 0;
  k = n1 - 1;
  emxFree_creal_T(&fy);
  while (k + 1 <= wwc->size[0]) {
    wwc_re = wwc->data[k].re;
    fv_re = fv->data[k].re;
    wwc_im = wwc->data[k].im;
    fv_im = fv->data[k].im;
    b_wwc_re = wwc->data[k].re;
    b_fv_im = fv->data[k].im;
    b_wwc_im = wwc->data[k].im;
    b_fv_re = fv->data[k].re;
    y->data[xidx].re = wwc_re * fv_re + wwc_im * fv_im;
    y->data[xidx].im = b_wwc_re * b_fv_im - b_wwc_im * b_fv_re;
    xidx++;
    k++;
  }

  emxFree_creal_T(&fv);
  emxFree_creal_T(&wwc);
}

//
// Arguments    : const emxArray_real_T *iPk
//                emxArray_real_T *idx
// Return Type  : void
//
static void e_findPeaksSeparatedByMoreThanM(const emxArray_real_T *iPk,
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
    if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
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
// EOGCFILT EEG filter for conversion to C.
//  Vectorize:
// Arguments    : double X_data[]
//                int X_size[1]
//                emxArray_real_T *Y
// Return Type  : void
//
static void eegcfilt(double X_data[], int X_size[1], emxArray_real_T *Y)
{
  int x_size_idx_0;
  int loop_ub;
  double x_data[1250];
  int i6;
  double d2;
  double d3;
  double y_data[1310];
  double a[10];
  double b_y_data[1310];
  static const double b_a[10] = { -2.1396152021655335E-5, -2.1396152489276133E-5,
    8.558460975207999E-5, 8.5584605288149449E-5, -0.00012837690837852629,
    -0.00012837691616921775, 8.5584610596008311E-5, 8.5584607376171939E-5,
    -2.1396151855180404E-5, -2.1396152098550849E-5 };

  int y_size[1];
  int b_y_size[1];
  int c_y_size[1];
  int tmp_data[1250];

  //  Fs = 250, N = 5
  //  flim = [8 18], bandpass
  if (X_size[0] == 1) {
    x_size_idx_0 = 1;
    x_data[0] = X_data[0];
  } else {
    x_size_idx_0 = X_size[0];
    loop_ub = X_size[0];
    for (i6 = 0; i6 < loop_ub; i6++) {
      x_data[i6] = X_data[i6];
    }
  }

  if (x_size_idx_0 == 0) {
    i6 = Y->size[0] * Y->size[1];
    Y->size[0] = 0;
    Y->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)Y, i6, (int)sizeof(double));
  } else {
    d2 = 2.0 * x_data[0];
    d3 = 2.0 * x_data[x_size_idx_0 - 1];
    for (i6 = 0; i6 < 30; i6++) {
      y_data[i6] = d2 - x_data[30 - i6];
    }

    for (i6 = 0; i6 < x_size_idx_0; i6++) {
      y_data[i6 + 30] = x_data[i6];
    }

    for (i6 = 0; i6 < 30; i6++) {
      y_data[(i6 + x_size_idx_0) + 30] = d3 - x_data[(x_size_idx_0 - i6) - 2];
    }

    for (i6 = 0; i6 < 10; i6++) {
      a[i6] = b_a[i6] * y_data[0];
    }

    y_size[0] = 60 + x_size_idx_0;
    loop_ub = 60 + x_size_idx_0;
    for (i6 = 0; i6 < loop_ub; i6++) {
      b_y_data[i6] = y_data[i6];
    }

    b_filter(b_y_data, y_size, a, y_data, b_y_size);
    flipud(y_data, b_y_size);
    for (i6 = 0; i6 < 10; i6++) {
      a[i6] = b_a[i6] * y_data[0];
    }

    c_y_size[0] = b_y_size[0];
    loop_ub = b_y_size[0];
    for (i6 = 0; i6 < loop_ub; i6++) {
      b_y_data[i6] = y_data[i6];
    }

    b_filter(b_y_data, c_y_size, a, y_data, b_y_size);
    flipud(y_data, b_y_size);
    if (X_size[0] == 1) {
      i6 = Y->size[0] * Y->size[1];
      Y->size[0] = 1;
      Y->size[1] = x_size_idx_0;
      emxEnsureCapacity((emxArray__common *)Y, i6, (int)sizeof(double));
      for (i6 = 0; i6 < x_size_idx_0; i6++) {
        Y->data[Y->size[0] * i6] = y_data[30 + i6];
      }
    } else {
      i6 = Y->size[0] * Y->size[1];
      Y->size[0] = x_size_idx_0;
      Y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)Y, i6, (int)sizeof(double));
      for (i6 = 0; i6 < x_size_idx_0; i6++) {
        tmp_data[i6] = 31 + i6;
      }

      for (i6 = 0; i6 < x_size_idx_0; i6++) {
        Y->data[i6] = y_data[tmp_data[i6] - 1];
      }
    }
  }
}

//
// Arguments    : const emxArray_real_T *x
//                creal_T y[2048]
// Return Type  : void
//
static void eml_fft(const emxArray_real_T *x, creal_T y[2048])
{
  int i;
  int b_x[1];
  emxArray_real_T c_x;
  static const double dv10[1025] = { 1.0, 0.99999529380957619,
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

  static const double dv11[1025] = { 0.0, -0.0030679567629659761,
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

  i = x->size[0];
  if (i == 0) {
    i = x->size[0];
    if (2048 > i) {
      for (i = 0; i < 2048; i++) {
        y[i].re = 0.0;
        y[i].im = 0.0;
      }
    }
  } else {
    b_x[0] = x->size[0];
    c_x = *x;
    c_x.size = (int *)&b_x;
    c_x.numDimensions = 1;
    r2br_r2dit_trig(&c_x, dv10, dv11, y);
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
// EOGCFILT EOG Filter for conversion to C. All inputs must be constant.
//  Vectorize:
// Arguments    : const double X[250]
//                double Y[250]
// Return Type  : void
//
static void eogcfilt(const double X[250], double Y[250])
{
  double dv1[250];

  //  Sampling Frequency = 250;
  // BW for 10Hz upper bound, Order of 3.
  // BW filt for 2Hz lower bound, Order of 3:
  filtfilt(X, dv1);
  b_filtfilt(dv1, Y);
}

//
// featureExtraction Summary of this function goes here if I ever feel like
// writing one up.
//  samplesX = samplesX(:);
// Arguments    : const double samplesX[250]
//                double F_data[]
//                int F_size[2]
// Return Type  : void
//
static void featureExtractionEOG(const double samplesX[250], double F_data[],
  int F_size[2])
{
  double peaks_data[500];
  int peaks_size[2];
  double loc_data[500];
  int loc_size[2];
  int T_count_findpeaks;
  double T_findpeaks_distX_data[1];
  boolean_T b_samplesX[250];
  int i0;
  findpeaks(samplesX, peaks_data, peaks_size, loc_data, loc_size);
  if (peaks_size[1] == 0) {
    T_count_findpeaks = 0;
    T_findpeaks_distX_data[0] = 0.0;
  } else {
    T_count_findpeaks = peaks_size[1];
    if (peaks_size[1] > 1) {
      T_findpeaks_distX_data[0] = loc_data[loc_size[1] - 1] - loc_data[0];

      // TODO: TAKE AVG, NOT MAX-MIN
    } else {
      T_findpeaks_distX_data[0] = 0.0;
    }
  }

  //  F = horzcat(T_mean, T_stdv, T_max, T_min, T_countmin_1, T_countmin_2, T_countmax, T_Integrate); 
  for (i0 = 0; i0 < 250; i0++) {
    b_samplesX[i0] = ((samplesX[i0] < -9.9999999999999991E-6) && (samplesX[i0] >
      -0.0001));
  }

  F_size[0] = 1;
  F_size[1] = 10;
  F_data[0] = Wmean(samplesX);
  F_data[1] = Wstd(samplesX);
  F_data[2] = Wmax(samplesX);
  F_data[3] = Wmin(samplesX);
  F_data[4] = sum(b_samplesX);
  F_data[5] = WCountMin(samplesX);
  F_data[6] = WCountMax(samplesX);
  F_data[7] = trapz(samplesX);
  F_data[8] = T_count_findpeaks;
  for (i0 = 0; i0 < 1; i0++) {
    F_data[9] = T_findpeaks_distX_data[0];
  }
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
//                double F[42]
// Return Type  : void
//
static void featureExtractionSSVEP(const emxArray_real_T *fch1, const
  emxArray_real_T *fch2, const emxArray_real_T *fch3, double Fs, double F[42])
{
  double threshFFT[8];
  int i7;
  signed char wLFFT[4];
  int i;
  double threshPSD[8];
  signed char wLPSD[4];
  emxArray_real_T *fchw;
  int fch1_idx_0;
  emxArray_int32_T *r1;
  emxArray_real_T *b_fch1;
  emxArray_real_T *b_fch2;
  emxArray_real_T *b_fch3;
  double FFT_Ltop[8];
  double FFT_MMM[4];
  double FFT_PkRatio[4];
  emxArray_real_T *PSD;
  double y;
  double PSD_Ltop[8];
  double PSD_MMM[4];
  double PSD_PkRatio[4];
  emxArray_real_T *hW;
  emxArray_real_T *fPSD;
  double f[1025];
  emxArray_real_T *PSD_PKS;
  emxArray_real_T *PSD_L;
  emxArray_real_T *r2;
  emxArray_real_T *b_PSD;
  emxArray_real_T *b_fchw;
  emxArray_real_T *c_fchw;
  static double FFT[4100];
  int ch;
  double b_FFT[1025];
  double dv8[1025];
  static double FFT_PKS_data[2050];
  int FFT_PKS_size[2];
  static double FFT_L_data[2050];
  int FFT_L_size[2];
  emxArray_real_T *c_PSD;
  int w;
  boolean_T exitg2;
  boolean_T exitg4;
  emxArray_real_T *d_PSD;
  boolean_T exitg1;
  boolean_T exitg3;
  double FFTPeaks1[4];
  int chn;
  boolean_T b0;
  boolean_T b3;
  double b_wLFFT[4];
  boolean_T b1;
  double c_wLFFT[4];
  boolean_T b4;
  double b_wLPSD[4];
  double c_wLPSD[4];
  double F0[32];
  double b_FFT_Ltop[4];
  double b_PSD_Ltop[4];

  // -% windows around certain target frequencies
  for (i7 = 0; i7 < 2; i7++) {
    threshFFT[i7 << 2] = 9.5 + 1.1300000000000008 * (double)i7;
    threshFFT[1 + (i7 << 2)] = 11.9 + 0.79999999999999893 * (double)i7;
    threshFFT[2 + (i7 << 2)] = 14.6 + 0.90000000000000036 * (double)i7;
    threshFFT[3 + (i7 << 2)] = 16.2 + 0.53999999999999915 * (double)i7;
  }

  for (i = 0; i < 4; i++) {
    wLFFT[i] = 0;
  }

  // ----PSD----%
  for (i7 = 0; i7 < 2; i7++) {
    threshPSD[i7 << 2] = 9.5 + (double)i7;
    threshPSD[1 + (i7 << 2)] = 12.0 + (double)i7;
    threshPSD[2 + (i7 << 2)] = 14.9 + 0.19999999999999929 * (double)i7;
    threshPSD[3 + (i7 << 2)] = 16.0 + (double)i7;
  }

  for (i = 0; i < 4; i++) {
    wLPSD[i] = 0;
  }

  emxInit_real_T(&fchw, 2);

  // ----PREALLOCATE----%
  i7 = fchw->size[0] * fchw->size[1];
  fchw->size[0] = 3;
  emxEnsureCapacity((emxArray__common *)fchw, i7, (int)sizeof(double));
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  i7 = fchw->size[0] * fchw->size[1];
  fchw->size[1] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)fchw, i7, (int)sizeof(double));
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  i = 3 * fch1_idx_0;
  for (i7 = 0; i7 < i; i7++) {
    fchw->data[i7] = 0.0;
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

  emxInit_int32_T(&r1, 1);
  i7 = r1->size[0];
  r1->size[0] = i;
  emxEnsureCapacity((emxArray__common *)r1, i7, (int)sizeof(int));
  for (i7 = 0; i7 < i; i7++) {
    r1->data[i7] = i7;
  }

  emxInit_real_T1(&b_fch1, 1);
  i7 = b_fch1->size[0];
  b_fch1->size[0] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch1, i7, (int)sizeof(double));
  for (i7 = 0; i7 < fch1_idx_0; i7++) {
    b_fch1->data[i7] = fch1->data[i7];
  }

  i = r1->size[0];
  for (i7 = 0; i7 < i; i7++) {
    fchw->data[fchw->size[0] * r1->data[i7]] = b_fch1->data[i7];
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

  i7 = r1->size[0];
  r1->size[0] = i;
  emxEnsureCapacity((emxArray__common *)r1, i7, (int)sizeof(int));
  for (i7 = 0; i7 < i; i7++) {
    r1->data[i7] = i7;
  }

  emxInit_real_T1(&b_fch2, 1);
  i7 = b_fch2->size[0];
  b_fch2->size[0] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch2, i7, (int)sizeof(double));
  for (i7 = 0; i7 < fch1_idx_0; i7++) {
    b_fch2->data[i7] = fch2->data[i7];
  }

  i = r1->size[0];
  for (i7 = 0; i7 < i; i7++) {
    fchw->data[1 + fchw->size[0] * r1->data[i7]] = b_fch2->data[i7];
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

  i7 = r1->size[0];
  r1->size[0] = i;
  emxEnsureCapacity((emxArray__common *)r1, i7, (int)sizeof(int));
  for (i7 = 0; i7 < i; i7++) {
    r1->data[i7] = i7;
  }

  emxInit_real_T1(&b_fch3, 1);
  i7 = b_fch3->size[0];
  b_fch3->size[0] = fch1_idx_0;
  emxEnsureCapacity((emxArray__common *)b_fch3, i7, (int)sizeof(double));
  for (i7 = 0; i7 < fch1_idx_0; i7++) {
    b_fch3->data[i7] = fch3->data[i7];
  }

  i = r1->size[0];
  for (i7 = 0; i7 < i; i7++) {
    fchw->data[2 + fchw->size[0] * r1->data[i7]] = b_fch3->data[i7];
  }

  emxFree_real_T(&b_fch3);
  emxFree_int32_T(&r1);

  //  all windows should be the same size:
  memset(&FFT_Ltop[0], 0, sizeof(double) << 3);
  for (i = 0; i < 4; i++) {
    FFT_MMM[i] = 0.0;
    FFT_PkRatio[i] = 0.0;
  }

  emxInit_real_T(&PSD, 2);
  fch1_idx_0 = fch1->size[0] * fch1->size[1];
  y = (double)fch1_idx_0 / 2.0;
  i7 = PSD->size[0] * PSD->size[1];
  PSD->size[0] = 4;
  PSD->size[1] = (int)y;
  emxEnsureCapacity((emxArray__common *)PSD, i7, (int)sizeof(double));
  i = (int)y << 2;
  for (i7 = 0; i7 < i; i7++) {
    PSD->data[i7] = 0.0;
  }

  memset(&PSD_Ltop[0], 0, sizeof(double) << 3);
  for (i = 0; i < 4; i++) {
    PSD_MMM[i] = 0.0;
    PSD_PkRatio[i] = 0.0;
  }

  emxInit_real_T1(&hW, 1);
  emxInit_real_T(&fPSD, 2);

  // 0?
  //  Data is already filtered:
  // between 250?500dp
  //  Preallocate for spd:
  emxInit_real_T(&PSD_PKS, 2);
  emxInit_real_T(&PSD_L, 2);
  emxInit_real_T(&r2, 2);
  emxInit_real_T(&b_PSD, 2);
  emxInit_real_T(&b_fchw, 2);
  emxInit_real_T(&c_fchw, 2);
  for (ch = 0; ch < 3; ch++) {
    //  #1 Take FFT:
    i = fchw->size[1];
    i7 = c_fchw->size[0] * c_fchw->size[1];
    c_fchw->size[0] = 1;
    c_fchw->size[1] = i;
    emxEnsureCapacity((emxArray__common *)c_fchw, i7, (int)sizeof(double));
    for (i7 = 0; i7 < i; i7++) {
      c_fchw->data[c_fchw->size[0] * i7] = fchw->data[ch + fchw->size[0] * i7];
    }

    get_nfft_data(c_fchw, Fs, f, dv8);

    //  #1.1 Find Peaks and M/I
    for (i7 = 0; i7 < 1025; i7++) {
      FFT[ch + (i7 << 2)] = dv8[i7];
      b_FFT[i7] = FFT[ch + (i7 << 2)];
    }

    b_findpeaks(b_FFT, FFT_PKS_data, FFT_PKS_size, FFT_L_data, FFT_L_size);
    if (FFT_PKS_size[1] > 1) {
      // Peak max minus min
      for (i7 = 0; i7 < 2; i7++) {
        FFT_Ltop[ch + (i7 << 2)] = f[(int)FFT_L_data[i7] - 1];
      }

      w = 0;
      exitg4 = false;
      while ((!exitg4) && (w < 4)) {
        if ((FFT_Ltop[ch] > threshFFT[w]) && (FFT_Ltop[ch] < threshFFT[4 + w]))
        {
          FFT_MMM[ch] = FFT_PKS_data[0] - FFT_PKS_data[1];
          FFT_PkRatio[ch] = FFT_PKS_data[0] / FFT_PKS_data[1];
          wLFFT[ch] = (signed char)(1 + w);
          exitg4 = true;
        } else {
          FFT_MMM[ch] = 0.0;
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
    i7 = b_fchw->size[0] * b_fchw->size[1];
    b_fchw->size[0] = 1;
    b_fchw->size[1] = i;
    emxEnsureCapacity((emxArray__common *)b_fchw, i7, (int)sizeof(double));
    for (i7 = 0; i7 < i; i7++) {
      b_fchw->data[b_fchw->size[0] * i7] = fchw->data[ch + fchw->size[0] * i7];
    }

    welch_psd(b_fchw, Fs, hW, r2, fPSD);
    i = r2->size[1];
    for (i7 = 0; i7 < i; i7++) {
      PSD->data[ch + PSD->size[0] * i7] = r2->data[r2->size[0] * i7];
    }

    // fin-start
    //  #2.2 Find Peaks and Max
    i = PSD->size[1];
    i7 = b_PSD->size[0] * b_PSD->size[1];
    b_PSD->size[0] = 1;
    b_PSD->size[1] = i;
    emxEnsureCapacity((emxArray__common *)b_PSD, i7, (int)sizeof(double));
    for (i7 = 0; i7 < i; i7++) {
      b_PSD->data[b_PSD->size[0] * i7] = PSD->data[ch + PSD->size[0] * i7];
    }

    c_findpeaks(b_PSD, PSD_PKS, PSD_L);
    if (PSD_PKS->size[1] > 1) {
      for (i7 = 0; i7 < 2; i7++) {
        PSD_Ltop[ch + (i7 << 2)] = fPSD->data[(int)PSD_L->data[i7] - 1];
      }

      w = 0;
      exitg3 = false;
      while ((!exitg3) && (w < 4)) {
        if ((PSD_Ltop[ch] >= threshPSD[w]) && (PSD_Ltop[ch] <= threshPSD[4 + w]))
        {
          PSD_MMM[ch] = PSD_PKS->data[0] - PSD_PKS->data[1];
          PSD_PkRatio[ch] = PSD_PKS->data[0] / PSD_PKS->data[1];
          wLPSD[ch] = (signed char)(1 + w);
          exitg3 = true;
        } else {
          PSD_MMM[ch] = 0.0;
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
  emxFree_real_T(&r2);
  emxFree_real_T(&hW);
  emxFree_real_T(&fchw);

  // Combine data into 'fourth' channel:
  for (i7 = 0; i7 < 1025; i7++) {
    FFT[3 + (i7 << 2)] = (FFT[i7 << 2] + FFT[1 + (i7 << 2)]) + FFT[2 + (i7 << 2)];
    b_FFT[i7] = FFT[3 + (i7 << 2)];
  }

  b_findpeaks(b_FFT, FFT_PKS_data, FFT_PKS_size, FFT_L_data, FFT_L_size);
  if (FFT_PKS_size[1] > 1) {
    for (i7 = 0; i7 < 2; i7++) {
      FFT_Ltop[3 + (i7 << 2)] = f[(int)FFT_L_data[i7] - 1];
    }

    w = 0;
    exitg2 = false;
    while ((!exitg2) && (w < 4)) {
      if ((FFT_Ltop[3] > threshFFT[w]) && (FFT_Ltop[3] < threshFFT[4 + w])) {
        FFT_MMM[3] = FFT_PKS_data[0] - FFT_PKS_data[1];
        FFT_PkRatio[3] = FFT_PKS_data[0] / FFT_PKS_data[1];
        wLFFT[3] = (signed char)(1 + w);
        exitg2 = true;
      } else {
        FFT_MMM[3] = 0.0;
        FFT_PkRatio[3] = 0.0;
        wLFFT[3] = 0;
        w++;
      }
    }
  }

  emxInit_real_T(&c_PSD, 2);
  i = PSD->size[1];
  i7 = c_PSD->size[0] * c_PSD->size[1];
  c_PSD->size[0] = 1;
  c_PSD->size[1] = i;
  emxEnsureCapacity((emxArray__common *)c_PSD, i7, (int)sizeof(double));
  for (i7 = 0; i7 < i; i7++) {
    c_PSD->data[c_PSD->size[0] * i7] = (PSD->data[PSD->size[0] * i7] + PSD->
      data[1 + PSD->size[0] * i7]) + PSD->data[2 + PSD->size[0] * i7];
  }

  i = c_PSD->size[1];
  for (i7 = 0; i7 < i; i7++) {
    PSD->data[3 + PSD->size[0] * i7] = c_PSD->data[c_PSD->size[0] * i7];
  }

  emxFree_real_T(&c_PSD);
  emxInit_real_T(&d_PSD, 2);
  i = PSD->size[1];
  i7 = d_PSD->size[0] * d_PSD->size[1];
  d_PSD->size[0] = 1;
  d_PSD->size[1] = i;
  emxEnsureCapacity((emxArray__common *)d_PSD, i7, (int)sizeof(double));
  for (i7 = 0; i7 < i; i7++) {
    d_PSD->data[d_PSD->size[0] * i7] = PSD->data[3 + PSD->size[0] * i7];
  }

  emxFree_real_T(&PSD);
  c_findpeaks(d_PSD, PSD_PKS, PSD_L);
  emxFree_real_T(&d_PSD);
  if (PSD_PKS->size[1] > 1) {
    //          PSD_L(4,:) = PSD_L(4,:)(:);
    for (i7 = 0; i7 < 2; i7++) {
      PSD_Ltop[3 + (i7 << 2)] = fPSD->data[(int)PSD_L->data[i7] - 1];
    }

    w = 0;
    exitg1 = false;
    while ((!exitg1) && (w < 4)) {
      if ((PSD_Ltop[3] >= threshPSD[w]) && (PSD_Ltop[3] <= threshPSD[4 + w])) {
        PSD_MMM[3] = PSD_PKS->data[0] - PSD_PKS->data[1];
        PSD_PkRatio[3] = PSD_PKS->data[0] / PSD_PKS->data[1];
        wLPSD[3] = (signed char)(1 + w);
        exitg1 = true;
      } else {
        PSD_MMM[3] = 0.0;
        PSD_PkRatio[3] = 0.0;
        wLPSD[3] = 0;
        w++;
      }
    }
  }

  emxFree_real_T(&PSD_L);
  emxFree_real_T(&PSD_PKS);
  emxFree_real_T(&fPSD);

  //  fprintf('FFT Matching Class %d %d %d %d \n',wLFFT(1),wLFFT(2)...
  //      ,wLFFT(3),wLFFT(4));
  //  fprintf('PSD Matching Class %d %d %d %d \n',wLPSD(1),wLPSD(2)...
  //      ,wLPSD(3),wLPSD(4));
  //  fprintf('FFT Peak Ratio: %1.3f\n',FFT_PkRatio);
  //  fprintf('PSD Peak Ratio: %1.3f\n',PSD_PkRatio);
  for (chn = 0; chn < 4; chn++) {
    FFTPeaks1[chn] = FFT_Ltop[chn];
  }

  //  fprintf('Avg FFTL: %1.3f \n',averageFFTPeak);
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
    b3 = isequal(b_wLFFT, c_wLFFT);
  } else {
    b3 = false;
  }

  //  fprintf('Avg PSDL: %1.3f \n',averagePSDPeak);
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

  //  fprintf('Booleans: [%d %d %d %d] \n',b1,b2,b3,b4);
  // /windowLength>=500
  //     %% Collect Feature data into 'F'
  // First separate features by channel: (row vects)
  //  first to remove: *FFT_Ltop(2) ... not sure how I will use this
  //  Also remove FFTPeaks2 and averageFFTPeak2
  //      F_1 = [Ch{1}.FFT_Ltop(1) Ch{1}.FFT_Ltop(2) Ch{1}.FFT_MMM Ch{1}.FFT_PkRatio Ch{1}.wLFFT ... 
  //          Ch{1}.PSD_Ltop(1) Ch{1}.PSD_Ltop(2) Ch{1}.PSD_MMM Ch{1}.PSD_PkRatio Ch{1}.wLPSD]; %10 features 
  for (ch = 0; ch < 4; ch++) {
    F0[ch] = FFT_Ltop[ch];
    F0[4 + ch] = FFT_MMM[ch];
    F0[8 + ch] = FFT_PkRatio[ch];
    F0[12 + ch] = wLFFT[ch];
    F0[16 + ch] = PSD_Ltop[ch];
    F0[20 + ch] = PSD_MMM[ch];
    F0[24 + ch] = PSD_PkRatio[ch];
    F0[28 + ch] = wLPSD[ch];
  }

  //      Extras = [FFTPeaks1 FFTPeaks2 averageFFTPeak averageFFTPeak2 averagePSDPeak b1 b2]; 
  b_FFT_Ltop[0] = FFT_Ltop[0];
  b_FFT_Ltop[1] = FFT_Ltop[1];
  b_FFT_Ltop[2] = FFT_Ltop[2];
  b_FFT_Ltop[3] = FFT_Ltop[3];
  b_PSD_Ltop[0] = PSD_Ltop[0];
  b_PSD_Ltop[1] = PSD_Ltop[1];
  b_PSD_Ltop[2] = PSD_Ltop[2];
  b_PSD_Ltop[3] = PSD_Ltop[3];
  for (i7 = 0; i7 < 8; i7++) {
    F[i7] = F0[i7 << 2];
    F[i7 + 8] = F0[1 + (i7 << 2)];
    F[i7 + 16] = F0[2 + (i7 << 2)];
    F[i7 + 24] = F0[3 + (i7 << 2)];
  }

  for (i7 = 0; i7 < 4; i7++) {
    F[i7 + 32] = FFTPeaks1[i7];
  }

  F[36] = b_mean(b_FFT_Ltop);
  F[37] = b_mean(b_PSD_Ltop);
  F[38] = b0;
  F[39] = b1;
  F[40] = b3;
  F[41] = b4;

  // END FUNCTION
}

//
// Arguments    : const emxArray_real_T *x
//                creal_T y[2048]
// Return Type  : void
//
static void fft(const emxArray_real_T *x, creal_T y[2048])
{
  int b_x[2];
  emxArray_real_T c_x;
  creal_T b_y1[2048];
  b_x[0] = x->size[1];
  b_x[1] = 1;
  c_x = *x;
  c_x.size = (int *)&b_x;
  c_x.numDimensions = 1;
  eml_fft(&c_x, b_y1);
  memcpy(&y[0], &b_y1[0], sizeof(creal_T) << 11);
}

//
// Arguments    : const double b[4]
//                const double a[4]
//                const double x[268]
//                const double zi[3]
//                double y[268]
// Return Type  : void
//
static void filter(const double b[4], const double a[4], const double x[268],
                   const double zi[3], double y[268])
{
  double dbuffer[4];
  int k;
  int j;
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
  double a[3];
  double b_y[268];
  static const double b_a[3] = { 0.99843298964950811, -1.5048763860936776,
    0.60567670985792155 };

  static const double dv2[4] = { 0.00156701035058832, 0.00470103105176495,
    0.00470103105176495, 0.00156701035058832 };

  static const double dv3[4] = { 1.0, -2.49860834469118, 2.11525412700316,
    -0.604109699507275 };

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

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 268U * sizeof(double));
  filter(dv2, dv3, b_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&c_y[0], &y[0], 268U * sizeof(double));
  filter(dv2, dv3, c_y, a, y);
  for (i = 0; i < 134; i++) {
    xtmp = y[i];
    y[i] = y[267 - i];
    y[267 - i] = xtmp;
  }

  memcpy(&y_out[0], &y[9], 250U * sizeof(double));
}

//
// Arguments    : const double yTemp[250]
//                double iPk_data[]
//                int iPk_size[1]
//                double iInflect_data[]
//                int iInflect_size[1]
// Return Type  : void
//
static void findLocalMaxima(const double yTemp[250], double iPk_data[], int
  iPk_size[1], double iInflect_data[], int iInflect_size[1])
{
  double b_yTemp[252];
  boolean_T yFinite[252];
  int ii;
  boolean_T x[251];
  int idx;
  unsigned char ii_data[251];
  boolean_T exitg3;
  boolean_T guard3 = false;
  int nx;
  unsigned char tmp_data[252];
  int i2;
  int i3;
  double yTemp_data[252];
  unsigned char iTemp_data[252];
  int yTemp_size[1];
  emxArray_real_T *r0;
  emxArray_real_T b_yTemp_data;
  int s_size[1];
  emxArray_boolean_T *b_x;
  double s_data[251];
  emxArray_real_T b_s_data;
  emxArray_int32_T *b_ii;
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

  idx = 0;
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 252)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      ii_data[idx - 1] = (unsigned char)ii;
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

  if (1 > idx) {
    nx = 0;
  } else {
    nx = idx;
  }

  tmp_data[0] = 1;
  for (i2 = 0; i2 < nx; i2++) {
    tmp_data[i2 + 1] = (unsigned char)(ii_data[i2] + 1);
  }

  if (1 > idx) {
    i3 = 0;
  } else {
    i3 = idx;
  }

  ii = 1 + i3;
  for (i2 = 0; i2 < ii; i2++) {
    iTemp_data[i2] = tmp_data[i2];
  }

  yTemp_size[0] = 1 + nx;
  nx++;
  for (i2 = 0; i2 < nx; i2++) {
    yTemp_data[i2] = b_yTemp[iTemp_data[i2] - 1];
  }

  emxInit_real_T1(&r0, 1);
  b_yTemp_data.data = (double *)&yTemp_data;
  b_yTemp_data.size = (int *)&yTemp_size;
  b_yTemp_data.allocatedSize = 252;
  b_yTemp_data.numDimensions = 1;
  b_yTemp_data.canFreeData = false;
  diff(&b_yTemp_data, r0);
  b_sign(r0);
  s_size[0] = r0->size[0];
  nx = r0->size[0];
  for (i2 = 0; i2 < nx; i2++) {
    s_data[i2] = r0->data[i2];
  }

  emxInit_boolean_T(&b_x, 1);
  b_s_data.data = (double *)&s_data;
  b_s_data.size = (int *)&s_size;
  b_s_data.allocatedSize = 251;
  b_s_data.numDimensions = 1;
  b_s_data.canFreeData = false;
  diff(&b_s_data, r0);
  i2 = b_x->size[0];
  b_x->size[0] = r0->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i2, (int)sizeof(boolean_T));
  nx = r0->size[0];
  for (i2 = 0; i2 < nx; i2++) {
    b_x->data[i2] = (r0->data[i2] < 0.0);
  }

  emxFree_real_T(&r0);
  emxInit_int32_T(&b_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i2 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i2, (int)sizeof(int));
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
      i2 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i2, (int)sizeof(int));
    }
  } else {
    i2 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i2, (int)sizeof(int));
  }

  if (1 > s_size[0] - 1) {
    nx = 0;
  } else {
    nx = s_size[0] - 1;
  }

  if (2 > s_size[0]) {
    i2 = 0;
  } else {
    i2 = 1;
  }

  ii = b_x->size[0];
  b_x->size[0] = nx;
  emxEnsureCapacity((emxArray__common *)b_x, ii, (int)sizeof(boolean_T));
  for (ii = 0; ii < nx; ii++) {
    b_x->data[ii] = (s_data[ii] != s_data[i2 + ii]);
  }

  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i2 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i2, (int)sizeof(int));
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
      i2 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i2, (int)sizeof(int));
    }
  } else {
    i2 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i2, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  iInflect_size[0] = c_ii->size[0];
  nx = c_ii->size[0];
  for (i2 = 0; i2 < nx; i2++) {
    iInflect_data[i2] = (double)iTemp_data[c_ii->data[i2]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  iPk_size[0] = b_ii->size[0];
  nx = b_ii->size[0];
  for (i2 = 0; i2 < nx; i2++) {
    iPk_data[i2] = (double)iTemp_data[b_ii->data[i2]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
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
  if (!rtIsNaN(x)) {
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
//                double Ypk_data[]
//                int Ypk_size[2]
//                double Xpk_data[]
//                int Xpk_size[2]
// Return Type  : void
//
static void findpeaks(const double Yin[250], double Ypk_data[], int Ypk_size[2],
                      double Xpk_data[], int Xpk_size[2])
{
  double iFinite_data[250];
  int iFinite_size[1];
  double iInfite_data[250];
  int iInfite_size[1];
  double iInflect_data[250];
  int iInflect_size[1];
  int iPk_size[1];
  int loop_ub;
  int i1;
  double iPk_data[500];
  emxArray_real_T *iPkOut;
  emxArray_int32_T *ia;
  emxArray_int32_T *ib;
  emxArray_real_T b_iPk_data;
  emxArray_real_T b_iInfite_data;
  double iPkOut_data[500];
  int iPkOut_size[1];
  getAllPeaks(Yin, iFinite_data, iFinite_size, iInfite_data, iInfite_size,
              iInflect_data, iInflect_size);
  removePeaksBelowMinPeakHeight(Yin, iFinite_data, iFinite_size);
  iPk_size[0] = iFinite_size[0];
  loop_ub = iFinite_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    iPk_data[i1] = iFinite_data[i1];
  }

  iFinite_size[0] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    iFinite_data[i1] = iPk_data[i1];
  }

  removePeaksBelowThreshold(Yin, iFinite_data, iFinite_size);
  iPk_size[0] = iFinite_size[0];
  loop_ub = iFinite_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    iPk_data[i1] = iFinite_data[i1];
  }

  emxInit_real_T1(&iPkOut, 1);
  emxInit_int32_T(&ia, 1);
  emxInit_int32_T(&ib, 1);
  b_iPk_data.data = (double *)&iPk_data;
  b_iPk_data.size = (int *)&iPk_size;
  b_iPk_data.allocatedSize = 500;
  b_iPk_data.numDimensions = 1;
  b_iPk_data.canFreeData = false;
  b_iInfite_data.data = (double *)&iInfite_data;
  b_iInfite_data.size = (int *)&iInfite_size;
  b_iInfite_data.allocatedSize = 250;
  b_iInfite_data.numDimensions = 1;
  b_iInfite_data.canFreeData = false;
  do_vectors(&b_iPk_data, &b_iInfite_data, iPkOut, ia, ib);
  c_findPeaksSeparatedByMoreThanM(iPkOut->size, iPk_data, iPk_size);
  keepAtMostNpPeaks(iPk_data, iPk_size);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&ia);
  for (i1 = 0; i1 < 250; i1++) {
    iFinite_data[i1] = 1.0 + (double)i1;
  }

  iPkOut_size[0] = iPk_size[0];
  loop_ub = iPk_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    iPkOut_data[i1] = iPkOut->data[(int)iPk_data[i1] - 1];
  }

  emxFree_real_T(&iPkOut);
  assignOutputs(Yin, iFinite_data, iPkOut_data, iPkOut_size, Ypk_data, Ypk_size,
                Xpk_data, Xpk_size);
}

//
// Arguments    : double x_data[]
//                int x_size[1]
// Return Type  : void
//
static void flipud(double x_data[], int x_size[1])
{
  int m;
  int md2;
  int i;
  double xtmp;
  m = x_size[0];
  md2 = x_size[0] >> 1;
  for (i = 1; i <= md2; i++) {
    xtmp = x_data[i - 1];
    x_data[i - 1] = x_data[m - i];
    x_data[m - i] = xtmp;
  }
}

//
// Arguments    : int nRows
//                boolean_T useRadix2
//                emxArray_real_T *costab
//                emxArray_real_T *sintab
//                emxArray_real_T *sintabinv
// Return Type  : void
//
static void generate_twiddle_tables(int nRows, boolean_T useRadix2,
  emxArray_real_T *costab, emxArray_real_T *sintab, emxArray_real_T *sintabinv)
{
  emxArray_real_T *costab1q;
  double e;
  int nRowsD4;
  int nd2;
  int k;
  int n2;
  emxInit_real_T(&costab1q, 2);
  e = 6.2831853071795862 / (double)nRows;
  nRowsD4 = nRows / 2 / 2;
  nd2 = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = nRowsD4 + 1;
  emxEnsureCapacity((emxArray__common *)costab1q, nd2, (int)sizeof(double));
  costab1q->data[0] = 1.0;
  nd2 = nRowsD4 / 2;
  for (k = 1; k <= nd2; k++) {
    costab1q->data[k] = std::cos(e * (double)k);
  }

  for (k = nd2 + 1; k < nRowsD4; k++) {
    costab1q->data[k] = std::sin(e * (double)(nRowsD4 - k));
  }

  costab1q->data[nRowsD4] = 0.0;
  if (!useRadix2) {
    nRowsD4 = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)costab, nd2, (int)sizeof(double));
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)sintab, nd2, (int)sizeof(double));
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    nd2 = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)sintabinv, nd2, (int)sizeof(double));
    for (k = 1; k <= nRowsD4; k++) {
      sintabinv->data[k] = costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      sintabinv->data[k] = costab1q->data[k - nRowsD4];
    }

    for (k = 1; k <= nRowsD4; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - nRowsD4];
    }
  } else {
    nRowsD4 = costab1q->size[1] - 1;
    n2 = (costab1q->size[1] - 1) << 1;
    nd2 = costab->size[0] * costab->size[1];
    costab->size[0] = 1;
    costab->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)costab, nd2, (int)sizeof(double));
    nd2 = sintab->size[0] * sintab->size[1];
    sintab->size[0] = 1;
    sintab->size[1] = n2 + 1;
    emxEnsureCapacity((emxArray__common *)sintab, nd2, (int)sizeof(double));
    costab->data[0] = 1.0;
    sintab->data[0] = 0.0;
    for (k = 1; k <= nRowsD4; k++) {
      costab->data[k] = costab1q->data[k];
      sintab->data[k] = -costab1q->data[nRowsD4 - k];
    }

    for (k = costab1q->size[1]; k <= n2; k++) {
      costab->data[k] = -costab1q->data[n2 - k];
      sintab->data[k] = -costab1q->data[k - nRowsD4];
    }

    nd2 = sintabinv->size[0] * sintabinv->size[1];
    sintabinv->size[0] = 1;
    sintabinv->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)sintabinv, nd2, (int)sizeof(double));
  }

  emxFree_real_T(&costab1q);
}

//
// Arguments    : const double y[250]
//                double iPk_data[]
//                int iPk_size[1]
//                double iInf_data[]
//                int iInf_size[1]
//                double iInflect_data[]
//                int iInflect_size[1]
// Return Type  : void
//
static void getAllPeaks(const double y[250], double iPk_data[], int iPk_size[1],
  double iInf_data[], int iInf_size[1], double iInflect_data[], int
  iInflect_size[1])
{
  boolean_T x[250];
  int idx;
  unsigned char ii_data[250];
  int ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  double yTemp[250];
  for (idx = 0; idx < 250; idx++) {
    x[idx] = (rtIsInf(y[idx]) && (y[idx] > 0.0));
  }

  idx = 0;
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii < 251)) {
    guard1 = false;
    if (x[ii - 1]) {
      idx++;
      ii_data[idx - 1] = (unsigned char)ii;
      if (idx >= 250) {
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

  if (1 > idx) {
    idx = 0;
  }

  iInf_size[0] = idx;
  for (ii = 0; ii < idx; ii++) {
    iInf_data[ii] = ii_data[ii];
  }

  memcpy(&yTemp[0], &y[0], 250U * sizeof(double));
  for (ii = 0; ii < idx; ii++) {
    ii_data[ii] = (unsigned char)iInf_data[ii];
  }

  for (ii = 0; ii < idx; ii++) {
    yTemp[ii_data[ii] - 1] = rtNaN;
  }

  findLocalMaxima(yTemp, iPk_data, iPk_size, iInflect_data, iInflect_size);
}

//
// Arguments    : int n1
//                boolean_T useRadix2
//                int *N2blue
//                int *nRows
// Return Type  : void
//
static void get_algo_sizes(int n1, boolean_T useRadix2, int *N2blue, int *nRows)
{
  int nn1m1;
  int pmax;
  int pmin;
  boolean_T exitg1;
  int p;
  int pow2p;
  *N2blue = 1;
  if (useRadix2) {
    *nRows = n1;
  } else {
    nn1m1 = (n1 + n1) - 1;
    pmax = 31;
    if (nn1m1 <= 1) {
      pmax = 0;
    } else {
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        p = (pmin + pmax) >> 1;
        pow2p = 1 << p;
        if (pow2p == nn1m1) {
          pmax = p;
          exitg1 = true;
        } else if (pow2p > nn1m1) {
          pmax = p;
        } else {
          pmin = p;
        }
      }
    }

    *N2blue = 1 << pmax;
    *nRows = *N2blue;
  }
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
  static creal_T A[2048];
  creal_T b_A[2048];
  int i8;
  double dv9[2048];
  fft(X, A);
  for (i8 = 0; i8 < 2048; i8++) {
    if (A[i8].im == 0.0) {
      b_A[i8].re = A[i8].re / 2048.0;
      b_A[i8].im = 0.0;
    } else if (A[i8].re == 0.0) {
      b_A[i8].re = 0.0;
      b_A[i8].im = A[i8].im / 2048.0;
    } else {
      b_A[i8].re = A[i8].re / 2048.0;
      b_A[i8].im = A[i8].im / 2048.0;
    }
  }

  b_abs(b_A, dv9);
  memcpy(&C[0], &dv9[0], 1025U * sizeof(double));
  for (i8 = 0; i8 < 2048; i8++) {
    if (A[i8].im == 0.0) {
      b_A[i8].re = A[i8].re / 2048.0;
      b_A[i8].im = 0.0;
    } else if (A[i8].re == 0.0) {
      b_A[i8].re = 0.0;
      b_A[i8].im = A[i8].im / 2048.0;
    } else {
      b_A[i8].re = A[i8].re / 2048.0;
      b_A[i8].im = A[i8].im / 2048.0;
    }
  }

  b_abs(b_A, dv9);
  for (i8 = 0; i8 < 1023; i8++) {
    C[1 + i8] = 2.0 * dv9[1 + i8];
  }

  for (i8 = 0; i8 < 1025; i8++) {
    f[i8] = Fs * (double)i8 / 2048.0;
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
    if (std::abs(cdiff) < 4.4408920985006262E-16 * std::abs(x - 1.0)) {
      ndbl++;
      apnd = x - 1.0;
    } else if (cdiff > 0.0) {
      apnd = ndbl - 1.0;
    } else {
      ndbl++;
    }

    n = (int)ndbl;
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
// Arguments    : double idx_data[]
//                int idx_size[1]
// Return Type  : void
//
static void keepAtMostNpPeaks(double idx_data[], int idx_size[1])
{
  double b_idx_data[500];
  int i22;
  if (idx_size[0] > 250) {
    memcpy(&b_idx_data[0], &idx_data[0], 250U * sizeof(double));
    idx_size[0] = 250;
    for (i22 = 0; i22 < 250; i22++) {
      idx_data[i22] = b_idx_data[i22];
    }
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
//                const emxArray_real_T *tX
//                const emxArray_real_T *tY
// Return Type  : double
//
static double knn(const double tsX[40], const emxArray_real_T *tX, const
                  emxArray_real_T *tY)
{
  double yfit;
  emxArray_int32_T *idx;
  int na;
  unsigned int outsize_idx_0;
  int pEnd;
  int i2;
  emxArray_int32_T *iwork;
  int n;
  emxArray_real_T *Uc;
  int k;
  boolean_T p;
  int i;
  int b_p;
  int j;
  int q;
  int qEnd;
  int kEnd;
  double x;
  int exitg3;
  double absxk;
  emxArray_real_T *a;
  int exponent;
  int i5;
  emxArray_real_T *y;
  int b_k;
  int c_k;
  emxArray_real_T *dist;
  emxArray_int32_T *pos;
  emxArray_real_T *edges;
  int exitg2;
  int exitg1;
  emxArray_real_T *nn;
  int b_exponent;
  emxInit_int32_T(&idx, 1);
  na = tY->size[0];
  outsize_idx_0 = (unsigned int)tY->size[0];
  pEnd = idx->size[0];
  idx->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, pEnd, (int)sizeof(int));
  i2 = (int)outsize_idx_0;
  for (pEnd = 0; pEnd < i2; pEnd++) {
    idx->data[pEnd] = 0;
  }

  if (tY->size[0] == 0) {
  } else {
    emxInit_int32_T(&iwork, 1);
    n = tY->size[0] + 1;
    pEnd = iwork->size[0];
    iwork->size[0] = (int)outsize_idx_0;
    emxEnsureCapacity((emxArray__common *)iwork, pEnd, (int)sizeof(int));
    for (k = 1; k <= n - 2; k += 2) {
      if ((tY->data[k - 1] <= tY->data[k]) || rtIsNaN(tY->data[k])) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((tY->size[0] & 1) != 0) {
      idx->data[tY->size[0] - 1] = tY->size[0];
    }

    i = 2;
    while (i < n - 1) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < n; pEnd = qEnd + i) {
        b_p = j;
        q = pEnd - 1;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if ((tY->data[idx->data[b_p - 1] - 1] <= tY->data[idx->data[q] - 1]) ||
              rtIsNaN(tY->data[idx->data[q] - 1])) {
            p = true;
          } else {
            p = false;
          }

          if (p) {
            iwork->data[k] = idx->data[b_p - 1];
            b_p++;
            if (b_p == pEnd) {
              while (q + 1 < qEnd) {
                k++;
                iwork->data[k] = idx->data[q];
                q++;
              }
            }
          } else {
            iwork->data[k] = idx->data[q];
            q++;
            if (q + 1 == qEnd) {
              while (b_p < pEnd) {
                k++;
                iwork->data[k] = idx->data[b_p - 1];
                b_p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&iwork);
  }

  emxInit_real_T1(&Uc, 1);
  outsize_idx_0 = (unsigned int)tY->size[0];
  pEnd = Uc->size[0];
  Uc->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)Uc, pEnd, (int)sizeof(double));
  for (k = 0; k + 1 <= na; k++) {
    Uc->data[k] = tY->data[idx->data[k] - 1];
  }

  emxFree_int32_T(&idx);
  k = 0;
  while ((k + 1 <= na) && rtIsInf(Uc->data[k]) && (Uc->data[k] < 0.0)) {
    k++;
  }

  b_p = k;
  k = tY->size[0];
  while ((k >= 1) && rtIsNaN(Uc->data[k - 1])) {
    k--;
  }

  pEnd = tY->size[0] - k;
  while ((k >= 1) && rtIsInf(Uc->data[k - 1]) && (Uc->data[k - 1] > 0.0)) {
    k--;
  }

  i = (tY->size[0] - k) - pEnd;
  q = -1;
  if (b_p > 0) {
    q = 0;
  }

  i2 = (b_p + k) - b_p;
  while (b_p + 1 <= i2) {
    x = Uc->data[b_p];
    do {
      exitg3 = 0;
      b_p++;
      if (b_p + 1 > i2) {
        exitg3 = 1;
      } else {
        absxk = std::abs(x / 2.0);
        if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
          if (absxk <= 2.2250738585072014E-308) {
            absxk = 4.94065645841247E-324;
          } else {
            frexp(absxk, &exponent);
            absxk = std::ldexp(1.0, exponent - 53);
          }
        } else {
          absxk = rtNaN;
        }

        if ((std::abs(x - Uc->data[b_p]) < absxk) || (rtIsInf(Uc->data[b_p]) &&
             rtIsInf(x) && ((Uc->data[b_p] > 0.0) == (x > 0.0)))) {
          p = true;
        } else {
          p = false;
        }

        if (!p) {
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);

    q++;
    Uc->data[q] = x;
  }

  if (i > 0) {
    q++;
    Uc->data[q] = Uc->data[i2];
  }

  b_p = i2 + i;
  for (j = 1; j <= pEnd; j++) {
    q++;
    Uc->data[q] = Uc->data[(b_p + j) - 1];
  }

  emxInit_real_T(&a, 2);
  pEnd = Uc->size[0];
  if (1 > q + 1) {
    i5 = -1;
  } else {
    i5 = q;
  }

  Uc->size[0] = i5 + 1;
  emxEnsureCapacity((emxArray__common *)Uc, pEnd, (int)sizeof(double));
  i = tY->size[0];
  pEnd = a->size[0] * a->size[1];
  a->size[0] = i;
  a->size[1] = 40;
  emxEnsureCapacity((emxArray__common *)a, pEnd, (int)sizeof(double));
  for (pEnd = 0; pEnd < i; pEnd++) {
    for (i2 = 0; i2 < 40; i2++) {
      a->data[pEnd + a->size[0] * i2] = tX->data[pEnd + tX->size[0] * i2] -
        tsX[i2];
    }
  }

  emxInit_real_T(&y, 2);
  pEnd = y->size[0] * y->size[1];
  y->size[0] = a->size[0];
  y->size[1] = 40;
  emxEnsureCapacity((emxArray__common *)y, pEnd, (int)sizeof(double));
  i = a->size[0] * 40;

//#pragma omp parallel for \
// num_threads(omp_get_max_threads()) \
// private(c_k)

  for (b_k = 1; b_k <= i; b_k++) {
    c_k = b_k;
    y->data[c_k - 1] = a->data[c_k - 1] * a->data[c_k - 1];
  }

  emxFree_real_T(&a);
  emxInit_real_T1(&dist, 1);
  pEnd = dist->size[0];
  dist->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)dist, pEnd, (int)sizeof(double));
  i = y->size[0];
  for (j = 0; j + 1 <= i; j++) {
    absxk = y->data[j];
    for (k = 0; k < 39; k++) {
      absxk += y->data[j + (k + 1) * i];
    }

    dist->data[j] = absxk;
  }

  emxFree_real_T(&y);
  emxInit_int32_T(&pos, 1);
  sort(dist, pos);
  pEnd = dist->size[0];
  dist->size[0] = pos->size[0];
  emxEnsureCapacity((emxArray__common *)dist, pEnd, (int)sizeof(double));
  i2 = pos->size[0];
  for (pEnd = 0; pEnd < i2; pEnd++) {
    dist->data[pEnd] = pos->data[pEnd];
  }

  emxFree_int32_T(&pos);
  emxInit_real_T(&edges, 2);
  i2 = Uc->size[0];
  pEnd = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = (int)(i2 + 1U);
  emxEnsureCapacity((emxArray__common *)edges, pEnd, (int)sizeof(double));
  k = 0;
  do {
    exitg2 = 0;
    i2 = Uc->size[0];
    if (k <= i2 - 2) {
      edges->data[1 + k] = Uc->data[k] + (Uc->data[1 + k] - Uc->data[k]) / 2.0;
      k++;
    } else {
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  edges->data[0] = rtMinusInf;
  edges->data[edges->size[1] - 1] = rtInf;
  k = 1;
  do {
    exitg1 = 0;
    i2 = Uc->size[0];
    if (k - 1 <= i2 - 2) {
      absxk = std::abs(edges->data[k]);
      if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
        if (absxk <= 2.2250738585072014E-308) {
          absxk = 4.94065645841247E-324;
        } else {
          frexp(absxk, &b_exponent);
          absxk = std::ldexp(1.0, b_exponent - 53);
        }
      } else {
        absxk = rtNaN;
      }

      edges->data[k] += absxk;
      k++;
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxInit_real_T1(&nn, 1);
  outsize_idx_0 = (unsigned int)edges->size[1];
  pEnd = nn->size[0];
  nn->size[0] = (int)outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)nn, pEnd, (int)sizeof(double));
  i2 = (int)outsize_idx_0;
  for (pEnd = 0; pEnd < i2; pEnd++) {
    nn->data[pEnd] = 0.0;
  }

  i = 0;
  for (k = 0; k < 3; k++) {
    i2 = findbin(tY->data[(int)dist->data[i] - 1], edges);
    if (i2 > 0) {
      nn->data[i2 - 1]++;
    }

    i++;
  }

  emxFree_real_T(&dist);
  pEnd = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = nn->size[0] - 1;
  emxEnsureCapacity((emxArray__common *)edges, pEnd, (int)sizeof(double));
  for (k = 0; k <= nn->size[0] - 2; k++) {
    edges->data[k] = nn->data[k];
  }

  if (nn->size[0] - 1 > 0) {
    edges->data[edges->size[1] - 1] += nn->data[nn->size[0] - 1];
  }

  emxFree_real_T(&nn);
  n = edges->size[1];
  absxk = edges->data[0];
  i = 0;
  if (edges->size[1] > 1) {
    for (i2 = 1; i2 + 1 <= n; i2++) {
      if (edges->data[i2] > absxk) {
        absxk = edges->data[i2];
        i = i2;
      }
    }
  }

  emxFree_real_T(&edges);
  yfit = Uc->data[i];
  emxFree_real_T(&Uc);
  return yfit;
}

//
// Arguments    : const emxArray_real_T *x
// Return Type  : double
//
static double mean(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[0]; k++) {
      y += x->data[k - 1];
    }
  }

  y /= (double)x->size[0];
  return y;
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
  if (nq == 0) {
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
      if (xwork->data[p] <= xwork->data[n]) {
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
          n = (iout - p) + 1;
          while (p + 1 <= np) {
            idx->data[n + p] = iwork->data[p];
            x->data[n + p] = xwork->data[p];
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
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
// Return Type  : void
//
static void merge_pow2_block(emxArray_int32_T *idx, emxArray_real_T *x, int
  offset)
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
        iwork[q - 1] = idx->data[blockOffset + q];
        xwork[q - 1] = x->data[blockOffset + q];
      }

      p = 0;
      q = bLen;
      do {
        exitg1 = 0;
        blockOffset++;
        if (xwork[p] >= xwork[q]) {
          idx->data[blockOffset] = iwork[p];
          x->data[blockOffset] = xwork[p];
          if (p + 1 < bLen) {
            p++;
          } else {
            exitg1 = 1;
          }
        } else {
          idx->data[blockOffset] = iwork[q];
          x->data[blockOffset] = xwork[q];
          if (q + 1 < bLen2) {
            q++;
          } else {
            q = blockOffset - p;
            while (p + 1 <= bLen) {
              idx->data[(q + p) + 1] = iwork[p];
              x->data[(q + p) + 1] = xwork[p];
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
// Arguments    : const emxArray_real_T *x
// Return Type  : int
//
static int nonSingletonDim(const emxArray_real_T *x)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  return dim;
}

//
// Arguments    : const double Y[1025]
//                const double iPk_data[]
//                double idx_data[]
//                int idx_size[1]
// Return Type  : void
//
static void orderPeaks(const double Y[1025], const double iPk_data[], double
  idx_data[], int idx_size[1])
{
  emxArray_real_T *x;
  int i24;
  int loop_ub;
  emxArray_int32_T *iidx;
  double b_idx_data[2050];
  int idx_size_idx_0;
  if (idx_size[0] == 0) {
  } else {
    emxInit_real_T1(&x, 1);
    i24 = x->size[0];
    x->size[0] = idx_size[0];
    emxEnsureCapacity((emxArray__common *)x, i24, (int)sizeof(double));
    loop_ub = idx_size[0];
    for (i24 = 0; i24 < loop_ub; i24++) {
      x->data[i24] = Y[(int)iPk_data[(int)idx_data[i24] - 1] - 1];
    }

    emxInit_int32_T(&iidx, 1);
    b_sort(x, iidx);
    idx_size_idx_0 = iidx->size[0];
    loop_ub = iidx->size[0];
    emxFree_real_T(&x);
    for (i24 = 0; i24 < loop_ub; i24++) {
      b_idx_data[i24] = idx_data[iidx->data[i24] - 1];
    }

    emxFree_int32_T(&iidx);
    idx_size[0] = idx_size_idx_0;
    for (i24 = 0; i24 < idx_size_idx_0; i24++) {
      idx_data[i24] = b_idx_data[i24];
    }
  }
}

//
// Arguments    : const emxArray_real_T *Yin
//                emxArray_real_T *y
//                emxArray_real_T *x
//                double *NpOut
// Return Type  : void
//
static void parse_inputs(const emxArray_real_T *Yin, emxArray_real_T *y,
  emxArray_real_T *x, double *NpOut)
{
  int cdiff;
  int nm1d2;
  int ndbl;
  int apnd;
  emxArray_real_T *b_y;
  cdiff = y->size[0];
  y->size[0] = Yin->size[1];
  emxEnsureCapacity((emxArray__common *)y, cdiff, (int)sizeof(double));
  nm1d2 = Yin->size[1];
  for (cdiff = 0; cdiff < nm1d2; cdiff++) {
    y->data[cdiff] = Yin->data[cdiff];
  }

  nm1d2 = Yin->size[1];
  if (nm1d2 < 1) {
    ndbl = 0;
    apnd = 0;
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

    if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)nm1d2) {
      ndbl++;
      apnd = Yin->size[1];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }
  }

  emxInit_real_T(&b_y, 2);
  cdiff = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)b_y, cdiff, (int)sizeof(double));
  if (ndbl > 0) {
    b_y->data[0] = 1.0;
    if (ndbl > 1) {
      b_y->data[ndbl - 1] = apnd;
      nm1d2 = (ndbl - 1) / 2;
      for (cdiff = 1; cdiff < nm1d2; cdiff++) {
        b_y->data[cdiff] = 1.0 + (double)cdiff;
        b_y->data[(ndbl - cdiff) - 1] = apnd - cdiff;
      }

      if (nm1d2 << 1 == ndbl - 1) {
        b_y->data[nm1d2] = (1.0 + (double)apnd) / 2.0;
      } else {
        b_y->data[nm1d2] = 1.0 + (double)nm1d2;
        b_y->data[nm1d2 + 1] = apnd - nm1d2;
      }
    }
  }

  cdiff = x->size[0];
  x->size[0] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)x, cdiff, (int)sizeof(double));
  nm1d2 = b_y->size[1];
  for (cdiff = 0; cdiff < nm1d2; cdiff++) {
    x->data[cdiff] = b_y->data[b_y->size[0] * cdiff];
  }

  emxFree_real_T(&b_y);
  nm1d2 = Yin->size[1];
  *NpOut = nm1d2;
}

//
// Arguments    : const emxArray_real_T *a
//                emxArray_real_T *y
// Return Type  : void
//
static void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  emxArray_real_T *x;
  int i16;
  int loop_ub;
  int k;
  int b_k;
  emxInit_real_T1(&x, 1);
  i16 = x->size[0];
  x->size[0] = a->size[0];
  emxEnsureCapacity((emxArray__common *)x, i16, (int)sizeof(double));
  loop_ub = a->size[0];
  for (i16 = 0; i16 < loop_ub; i16++) {
    x->data[i16] = a->data[i16];
  }

  loop_ub = a->size[0];
  i16 = y->size[0];
  y->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)y, i16, (int)sizeof(double));
  loop_ub = a->size[0];

//#pragma omp parallel for \
// num_threads(omp_get_max_threads()) \
// private(b_k)

  for (k = 1; k <= loop_ub; k++) {
    b_k = k;
    y->data[b_k - 1] = x->data[b_k - 1] * x->data[b_k - 1];
  }

  emxFree_real_T(&x);
}

//
// Arguments    : const emxArray_real_T *x
//                const double costab[1025]
//                const double sintab[1025]
//                creal_T y[2048]
// Return Type  : void
//
static void r2br_r2dit_trig(const emxArray_real_T *x, const double costab[1025],
  const double sintab[1025], creal_T y[2048])
{
  int iDelta2;
  int ix;
  int i;
  int ju;
  int iy;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  if (x->size[0] <= 2048) {
    iDelta2 = x->size[0];
  } else {
    iDelta2 = 2048;
  }

  if (2048 > x->size[0]) {
    for (i = 0; i < 2048; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 1; i < iDelta2; i++) {
    y[iy].re = x->data[ix];
    y[iy].im = 0.0;
    iy = 2048;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy].re = x->data[ix];
  y[iy].im = 0.0;
  for (i = 0; i <= 2047; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  iDelta2 = 4;
  ix = 512;
  ju = 2045;
  while (ix > 0) {
    for (i = 0; i < ju; i += iDelta2) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ix; j < 1024; j += ix) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += iDelta2;
      }

      istart++;
    }

    ix /= 2;
    iy = iDelta2;
    iDelta2 <<= 1;
    ju -= iy;
  }
}

//
// Arguments    : const emxArray_creal_T *x
//                int unsigned_nRows
//                const emxArray_real_T *costab
//                const emxArray_real_T *sintab
//                emxArray_creal_T *y
// Return Type  : void
//
static void r2br_r2dit_trig_impl(const emxArray_creal_T *x, int unsigned_nRows,
  const emxArray_real_T *costab, const emxArray_real_T *sintab, emxArray_creal_T
  *y)
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

  ix = 0;
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
// Arguments    : const double Y[250]
//                double iPk_data[]
//                int iPk_size[1]
// Return Type  : void
//
static void removePeaksBelowMinPeakHeight(const double Y[250], double iPk_data[],
  int iPk_size[1])
{
  int end;
  int trueCount;
  int i;
  int partialTrueCount;
  if (!(iPk_size[0] == 0)) {
    end = iPk_size[0] - 1;
    trueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y[(int)iPk_data[i] - 1] > 6.5E-5) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y[(int)iPk_data[i] - 1] > 6.5E-5) {
        iPk_data[partialTrueCount] = iPk_data[i];
        partialTrueCount++;
      }
    }

    iPk_size[0] = trueCount;
  }
}

//
// Arguments    : const double Y[250]
//                double iPk_data[]
//                int iPk_size[1]
// Return Type  : void
//
static void removePeaksBelowThreshold(const double Y[250], double iPk_data[],
  int iPk_size[1])
{
  emxArray_real_T *maxval;
  unsigned char csz_idx_0;
  int k;
  int trueCount;
  int i;
  int partialTrueCount;
  emxInit_real_T1(&maxval, 1);
  csz_idx_0 = (unsigned char)iPk_size[0];
  k = maxval->size[0];
  maxval->size[0] = (unsigned char)iPk_size[0];
  emxEnsureCapacity((emxArray__common *)maxval, k, (int)sizeof(double));
  for (k = 0; k + 1 <= csz_idx_0; k++) {
    if ((Y[(int)(iPk_data[k] - 1.0) - 1] >= Y[(int)(iPk_data[k] + 1.0) - 1]) ||
        rtIsNaN(Y[(int)(iPk_data[k] + 1.0) - 1])) {
      maxval->data[k] = Y[(int)(iPk_data[k] - 1.0) - 1];
    } else {
      maxval->data[k] = Y[(int)(iPk_data[k] + 1.0) - 1];
    }
  }

  k = iPk_size[0] - 1;
  trueCount = 0;
  for (i = 0; i <= k; i++) {
    if (Y[(int)iPk_data[i] - 1] - maxval->data[i] >= 0.0) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (i = 0; i <= k; i++) {
    if (Y[(int)iPk_data[i] - 1] - maxval->data[i] >= 0.0) {
      iPk_data[partialTrueCount] = iPk_data[i];
      partialTrueCount++;
    }
  }

  emxFree_real_T(&maxval);
  iPk_size[0] = trueCount;
}

//
// Arguments    : const emxArray_real_T *a
//                emxArray_real_T *b
// Return Type  : void
//
static void repmat(const emxArray_real_T *a, emxArray_real_T *b)
{
  int outsize_idx_0;
  int i14;
  outsize_idx_0 = a->size[0];
  i14 = b->size[0];
  b->size[0] = outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)b, i14, (int)sizeof(double));
  if ((!(a->size[0] == 0)) && (!(outsize_idx_0 == 0))) {
    for (outsize_idx_0 = 0; outsize_idx_0 + 1 <= a->size[0]; outsize_idx_0++) {
      b->data[outsize_idx_0] = a->data[outsize_idx_0];
    }
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
  a = std::abs(u0);
  b = std::abs(u1);
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
  double absxk;
  int exponent;
  boolean_T p;
  xk = x->data[*k - 1];
  exitg1 = false;
  while ((!exitg1) && (*k < x->size[0])) {
    absxk = std::abs(xk / 2.0);
    if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
      if (absxk <= 2.2250738585072014E-308) {
        absxk = 4.94065645841247E-324;
      } else {
        frexp(absxk, &exponent);
        absxk = std::ldexp(1.0, exponent - 53);
      }
    } else {
      absxk = rtNaN;
    }

    if ((std::abs(xk - x->data[*k]) < absxk) || (rtIsInf(x->data[*k]) && rtIsInf
         (xk) && ((x->data[*k] > 0.0) == (xk > 0.0)))) {
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
  int i23;
  emxArray_real_T *vwork;
  int j;
  int vstride;
  int k;
  emxArray_int32_T *iidx;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  if (dim <= 1) {
    i23 = x->size[0];
  } else {
    i23 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  j = vwork->size[0];
  vwork->size[0] = i23;
  emxEnsureCapacity((emxArray__common *)vwork, j, (int)sizeof(double));
  vstride = x->size[0];
  j = idx->size[0];
  idx->size[0] = vstride;
  emxEnsureCapacity((emxArray__common *)idx, j, (int)sizeof(int));
  vstride = 1;
  k = 1;
  while (k <= dim - 1) {
    vstride *= x->size[0];
    k = 2;
  }

  j = 0;
  emxInit_int32_T(&iidx, 1);
  while (j + 1 <= vstride) {
    for (k = 0; k + 1 <= i23; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, iidx);
    for (k = 0; k + 1 <= i23; k++) {
      x->data[j + k * vstride] = vwork->data[k];
      idx->data[j + k * vstride] = iidx->data[k];
    }

    j++;
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_real_T *b_x;
  int ib;
  int wOffset;
  int m;
  int n;
  double x4[4];
  int idx4[4];
  emxArray_int32_T *iwork;
  emxArray_real_T *xwork;
  int nNaNs;
  int k;
  signed char perm[4];
  int nNonNaN;
  int i3;
  int i4;
  int nBlocks;
  int b_iwork[256];
  double b_xwork[256];
  int bLen2;
  int nPairs;
  int exitg1;
  emxInit_real_T1(&b_x, 1);
  ib = x->size[0];
  wOffset = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, wOffset, (int)sizeof(double));
  m = x->size[0];
  for (wOffset = 0; wOffset < m; wOffset++) {
    b_x->data[wOffset] = x->data[wOffset];
  }

  wOffset = idx->size[0];
  idx->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)idx, wOffset, (int)sizeof(int));
  for (wOffset = 0; wOffset < ib; wOffset++) {
    idx->data[wOffset] = 0;
  }

  n = x->size[0];
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  emxInit_int32_T(&iwork, 1);
  wOffset = iwork->size[0];
  iwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)iwork, wOffset, (int)sizeof(int));
  m = iwork->size[0];
  wOffset = iwork->size[0];
  iwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)iwork, wOffset, (int)sizeof(int));
  for (wOffset = 0; wOffset < m; wOffset++) {
    iwork->data[wOffset] = 0;
  }

  emxInit_real_T1(&xwork, 1);
  m = x->size[0];
  wOffset = xwork->size[0];
  xwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)xwork, wOffset, (int)sizeof(double));
  m = xwork->size[0];
  wOffset = xwork->size[0];
  xwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)xwork, wOffset, (int)sizeof(double));
  for (wOffset = 0; wOffset < m; wOffset++) {
    xwork->data[wOffset] = 0.0;
  }

  nNaNs = 1;
  ib = 0;
  for (k = 0; k + 1 <= n; k++) {
    if (rtIsNaN(b_x->data[k])) {
      idx->data[n - nNaNs] = k + 1;
      xwork->data[n - nNaNs] = b_x->data[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = b_x->data[k];
      if (ib == 4) {
        ib = k - nNaNs;
        if (x4[0] <= x4[1]) {
          m = 1;
          wOffset = 2;
        } else {
          m = 2;
          wOffset = 1;
        }

        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }

        if (x4[m - 1] <= x4[i3 - 1]) {
          if (x4[wOffset - 1] <= x4[i3 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)i3;
            perm[3] = (signed char)i4;
          } else if (x4[wOffset - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else if (x4[m - 1] <= x4[i4 - 1]) {
          if (x4[wOffset - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)m;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else {
          perm[0] = (signed char)i3;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)m;
          perm[3] = (signed char)wOffset;
        }

        idx->data[ib - 2] = idx4[perm[0] - 1];
        idx->data[ib - 1] = idx4[perm[1] - 1];
        idx->data[ib] = idx4[perm[2] - 1];
        idx->data[ib + 1] = idx4[perm[3] - 1];
        b_x->data[ib - 2] = x4[perm[0] - 1];
        b_x->data[ib - 1] = x4[perm[1] - 1];
        b_x->data[ib] = x4[perm[2] - 1];
        b_x->data[ib + 1] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  wOffset = x->size[0] - nNaNs;
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
      idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
      b_x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
    }
  }

  m = (nNaNs - 1) >> 1;
  for (k = 1; k <= m; k++) {
    ib = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[n - k];
    idx->data[n - k] = ib;
    b_x->data[wOffset + k] = xwork->data[n - k];
    b_x->data[n - k] = xwork->data[wOffset + k];
  }

  if (((nNaNs - 1) & 1) != 0) {
    b_x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
  }

  nNonNaN = (x->size[0] - nNaNs) + 1;
  m = 2;
  if (nNonNaN > 1) {
    if (x->size[0] >= 256) {
      nBlocks = nNonNaN >> 8;
      if (nBlocks > 0) {
        for (i3 = 1; i3 <= nBlocks; i3++) {
          i4 = ((i3 - 1) << 8) - 1;
          for (nNaNs = 0; nNaNs < 6; nNaNs++) {
            n = 1 << (nNaNs + 2);
            bLen2 = n << 1;
            nPairs = 256 >> (nNaNs + 3);
            for (k = 1; k <= nPairs; k++) {
              m = i4 + (k - 1) * bLen2;
              for (ib = 1; ib <= bLen2; ib++) {
                b_iwork[ib - 1] = idx->data[m + ib];
                b_xwork[ib - 1] = b_x->data[m + ib];
              }

              wOffset = 0;
              ib = n;
              do {
                exitg1 = 0;
                m++;
                if (b_xwork[wOffset] <= b_xwork[ib]) {
                  idx->data[m] = b_iwork[wOffset];
                  b_x->data[m] = b_xwork[wOffset];
                  if (wOffset + 1 < n) {
                    wOffset++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx->data[m] = b_iwork[ib];
                  b_x->data[m] = b_xwork[ib];
                  if (ib + 1 < bLen2) {
                    ib++;
                  } else {
                    ib = m - wOffset;
                    while (wOffset + 1 <= n) {
                      idx->data[(ib + wOffset) + 1] = b_iwork[wOffset];
                      b_x->data[(ib + wOffset) + 1] = b_xwork[wOffset];
                      wOffset++;
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

  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
  wOffset = x->size[0];
  x->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)x, wOffset, (int)sizeof(double));
  m = b_x->size[0];
  for (wOffset = 0; wOffset < m; wOffset++) {
    x->data[wOffset] = b_x->data[wOffset];
  }

  emxFree_real_T(&b_x);
}

//
// Arguments    : const boolean_T x[250]
// Return Type  : double
//
static double sum(const boolean_T x[250])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 249; k++) {
    y += (double)x[k + 1];
  }

  return y;
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
  int apnd;
  int i13;
  int varargin_1;
  double y;
  double bnew;
  int ndbl;
  double cdiff;
  emxArray_real_T *b_y;
  int k;
  emxArray_real_T *data_taper;
  emxArray_real_T *S;
  double back_shift;
  double z;
  double number_of_blocks;
  int a;
  emxArray_real_T *Data_Block;
  emxArray_real_T *P;
  emxArray_creal_T *b_Data_Block;
  emxArray_real_T *r4;
  emxArray_creal_T *c_Data_Block;
  emxArray_int32_T *r5;
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
  apnd = window->size[0];
  i13 = window->size[0];
  window->size[0] = apnd;
  emxEnsureCapacity((emxArray__common *)window, i13, (int)sizeof(double));
  varargin_1 = window->size[0];
  y = (double)window->size[0] / 2.0;
  bnew = (double)window->size[0] / 2.0 - 1.0;
  if (y - 1.0 < 0.0) {
    ndbl = 0;
    bnew = y - 1.0;
  } else {
    ndbl = (int)std::floor(bnew + 0.5);
    apnd = ndbl;
    cdiff = (double)ndbl - (y - 1.0);
    if (std::abs(cdiff) < 4.4408920985006262E-16 * std::abs(bnew)) {
      ndbl++;
      bnew = y - 1.0;
    } else if (cdiff > 0.0) {
      bnew = (double)ndbl - 1.0;
    } else {
      ndbl++;
      bnew = apnd;
    }
  }

  emxInit_real_T(&b_y, 2);
  i13 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)b_y, i13, (int)sizeof(double));
  if (ndbl > 0) {
    b_y->data[0] = 0.0;
    if (ndbl > 1) {
      b_y->data[ndbl - 1] = bnew;
      apnd = (ndbl - 1) / 2;
      for (k = 1; k < apnd; k++) {
        b_y->data[k] = k;
        b_y->data[(ndbl - k) - 1] = bnew - (double)k;
      }

      if (apnd << 1 == ndbl - 1) {
        b_y->data[apnd] = bnew / 2.0;
      } else {
        b_y->data[apnd] = apnd;
        b_y->data[apnd + 1] = bnew - (double)apnd;
      }
    }
  }

  i13 = frequencies->size[0] * frequencies->size[1];
  frequencies->size[0] = 1;
  frequencies->size[1] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)frequencies, i13, (int)sizeof(double));
  apnd = window->size[0];
  k = b_y->size[0] * b_y->size[1];
  for (i13 = 0; i13 < k; i13++) {
    frequencies->data[i13] = b_y->data[i13] * fs / (double)apnd;
  }

  emxFree_real_T(&b_y);
  emxInit_real_T1(&data_taper, 1);
  emxInit_real_T1(&S, 1);

  // must be even, best if 2^n
  back_shift = (double)window->size[0] / 2.0;

  // ORIGINAL;
  apnd = signals->size[1];
  z = 2.0 * (double)apnd / (double)window->size[0];
  number_of_blocks = std::floor(z) - 1.0;
  repmat(window, data_taper);

  //  Data segmentation into blocks of size block_samples:
  y = (double)window->size[0] / 2.0;
  i13 = S->size[0];
  S->size[0] = (int)y;
  emxEnsureCapacity((emxArray__common *)S, i13, (int)sizeof(double));
  k = (int)y;
  for (i13 = 0; i13 < k; i13++) {
    S->data[i13] = 0.0;
  }

  // ORIGINAL
  //  S = zeros(ceil(block_samples/2),number_of_signals.^2);
  a = 0;
  emxInit_real_T1(&Data_Block, 1);
  emxInit_real_T1(&P, 1);
  emxInit_creal_T(&b_Data_Block, 1);
  emxInit_real_T1(&r4, 1);
  emxInit_creal_T(&c_Data_Block, 1);
  while (a <= (int)number_of_blocks - 1) {
    //  Retrieve current data block
    bnew = ((1.0 + (double)a) - 1.0) * back_shift + 1.0;
    cdiff = (double)varargin_1 + ((1.0 + (double)a) - 1.0) * back_shift;
    if (bnew > cdiff) {
      i13 = 0;
      ndbl = 0;
    } else {
      i13 = (int)bnew - 1;
      ndbl = (int)cdiff;
    }

    apnd = Data_Block->size[0];
    Data_Block->size[0] = ndbl - i13;
    emxEnsureCapacity((emxArray__common *)Data_Block, apnd, (int)sizeof(double));
    k = ndbl - i13;
    for (ndbl = 0; ndbl < k; ndbl++) {
      Data_Block->data[ndbl] = signals->data[i13 + ndbl];
    }

    b_repmat(mean(Data_Block), (double)varargin_1, r4);
    i13 = Data_Block->size[0];
    emxEnsureCapacity((emxArray__common *)Data_Block, i13, (int)sizeof(double));
    k = Data_Block->size[0];
    for (i13 = 0; i13 < k; i13++) {
      Data_Block->data[i13] -= r4->data[i13];
    }

    i13 = Data_Block->size[0];
    emxEnsureCapacity((emxArray__common *)Data_Block, i13, (int)sizeof(double));
    k = Data_Block->size[0];
    for (i13 = 0; i13 < k; i13++) {
      Data_Block->data[i13] *= data_taper->data[i13];
    }

    // Taper it
    b_fft(Data_Block, b_Data_Block);

    // FFT it,
    //  bilateral DFT
    //  viii
    bnew = (double)varargin_1 / 2.0;
    if (1.0 > bnew) {
      k = 0;
    } else {
      k = (int)bnew;
    }

    i13 = c_Data_Block->size[0];
    c_Data_Block->size[0] = k;
    emxEnsureCapacity((emxArray__common *)c_Data_Block, i13, (int)sizeof(creal_T));
    for (i13 = 0; i13 < k; i13++) {
      c_Data_Block->data[i13] = b_Data_Block->data[i13];
    }

    i13 = b_Data_Block->size[0];
    b_Data_Block->size[0] = c_Data_Block->size[0];
    emxEnsureCapacity((emxArray__common *)b_Data_Block, i13, (int)sizeof(creal_T));
    k = c_Data_Block->size[0];
    for (i13 = 0; i13 < k; i13++) {
      b_Data_Block->data[i13] = c_Data_Block->data[i13];
    }

    // ORIGINAL
    //  Data_Block = Data_Block(1:ceil(block_samples/2),:);
    // All spectral combinations:
    y = (double)varargin_1 / 2.0;
    i13 = P->size[0];
    P->size[0] = (int)y;
    emxEnsureCapacity((emxArray__common *)P, i13, (int)sizeof(double));
    k = (int)y;
    for (i13 = 0; i13 < k; i13++) {
      P->data[i13] = 0.0;
    }

    // ORIGINAL
    //  P = zeros(ceil(block_samples/2)/2,number_of_signals.^2);
    //  P(:,c) = Data_Block(:,b).*conj(Data_Block(:,aa)); % THIS
    //  IS FOR WIND TUNNEL EESC-USP BEAMFORMING CODE
    //  P(:,c) = Data_Block(:,aa).*conj(Data_Block(:,b)); % THIS IS THE ORIGINAL 
    //  LINE
    k = b_Data_Block->size[0] - 1;
    for (i13 = 0; i13 <= k; i13++) {
      bnew = b_Data_Block->data[i13].re;
      cdiff = -b_Data_Block->data[i13].im;
      bnew = b_Data_Block->data[i13].re * bnew - b_Data_Block->data[i13].im *
        cdiff;
      P->data[i13] = bnew;
    }

    //  IS FOR FAN RIG BEAMFORMING CODE
    //  Sum the spectrums up ...
    i13 = S->size[0];
    emxEnsureCapacity((emxArray__common *)S, i13, (int)sizeof(double));
    k = S->size[0];
    for (i13 = 0; i13 < k; i13++) {
      S->data[i13] += P->data[i13];
    }

    a++;
  }

  emxFree_creal_T(&c_Data_Block);
  emxFree_creal_T(&b_Data_Block);
  emxFree_real_T(&P);
  emxFree_real_T(&Data_Block);
  power(window, r4);
  bnew = b_sum(r4) * fs * (std::floor(z) - 1.0);
  i13 = data_taper->size[0];
  data_taper->size[0] = S->size[0];
  emxEnsureCapacity((emxArray__common *)data_taper, i13, (int)sizeof(double));
  k = S->size[0];
  emxFree_real_T(&r4);
  for (i13 = 0; i13 < k; i13++) {
    data_taper->data[i13] = S->data[i13] * 2.0 / bnew;
  }

  i13 = S->size[0];
  emxEnsureCapacity((emxArray__common *)S, i13, (int)sizeof(double));
  k = S->size[0];
  for (i13 = 0; i13 < k; i13++) {
    S->data[i13] = S->data[i13] * 2.0 / bnew;
  }

  //  Average them out
  i13 = CSM->size[0] * CSM->size[1];
  CSM->size[0] = 1;
  CSM->size[1] = data_taper->size[0];
  emxEnsureCapacity((emxArray__common *)CSM, i13, (int)sizeof(double));
  k = data_taper->size[0];
  emxFree_real_T(&data_taper);
  for (i13 = 0; i13 < k; i13++) {
    CSM->data[i13] = 0.0;
  }

  emxInit_int32_T(&r5, 1);

  //  for a = 1:sensors
  apnd = S->size[0];
  i13 = r5->size[0];
  r5->size[0] = apnd;
  emxEnsureCapacity((emxArray__common *)r5, i13, (int)sizeof(int));
  for (i13 = 0; i13 < apnd; i13++) {
    r5->data[i13] = i13;
  }

  emxInit_real_T1(&b_S, 1);
  k = S->size[0];
  i13 = b_S->size[0];
  b_S->size[0] = k;
  emxEnsureCapacity((emxArray__common *)b_S, i13, (int)sizeof(double));
  for (i13 = 0; i13 < k; i13++) {
    b_S->data[i13] = S->data[i13];
  }

  emxFree_real_T(&S);
  apnd = r5->size[0];
  for (i13 = 0; i13 < apnd; i13++) {
    CSM->data[CSM->size[0] * r5->data[i13]] = b_S->data[i13];
  }

  emxFree_real_T(&b_S);
  emxFree_int32_T(&r5);

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
  emxInit_real_T(&emx, numDimensions);
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
  emxInit_real_T(&emx, numDimensions);
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
  emxInit_real_T(&emx, 2);
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
  emxInit_real_T(&emx, 2);
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
  emxInit_real_T(pEmxArray, numDimensions);
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
//  tXSSVEP    = SSVEP Training Data (not sure if I will include just yet).
//  tYSSVEP    = SSVEP Training Classes
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
// Arguments    : const double ch1_data[]
//                const int ch1_size[1]
//                const double ch2_data[]
//                const int ch2_size[1]
//                const double ch3_data[]
//                const int ch3_size[1]
//                const double ch4_data[]
//                const int ch4_size[1]
//                const emxArray_real_T *tXEOG
//                const emxArray_real_T *tYEOG
//                const emxArray_real_T *tXSSVEP
//                const emxArray_real_T *tYSSVEP
//                double Fs
// Return Type  : double
//
double fullHybridClassifier(const double ch1_data[], const int ch1_size[1],
  const double ch2_data[], const int ch2_size[1], const double ch3_data[], const
  int ch3_size[1], const double ch4_data[], const int ch4_size[1], const
  emxArray_real_T *tXEOG, const emxArray_real_T *tYEOG, const emxArray_real_T
  *tXSSVEP, const emxArray_real_T *tYSSVEP, double Fs)
{
  double Y;
  boolean_T DB;
  double ch1f[250];
  double ch2f[250];
  double ch3f[250];
  double ch4f[250];
  double tmp_data[10];
  int tmp_size[2];
  double b_tmp_data[10];
  double c_tmp_data[10];
  double d_tmp_data[10];
  double dv0[40];
  int mtmp;
  short varargin_1[3];
  int ix;
  int ln;
  int loop_ub;
  int b_loop_ub;
  double b_ch1_data[1250];
  int b_ch1_size[1];
  emxArray_real_T *b_ch1f;
  int b_ch2_size[1];
  emxArray_real_T *b_ch2f;
  int b_ch3_size[1];
  emxArray_real_T *b_ch3f;
  double F[42];

  //  Window length??
  Y = 0.0;

  // Default Value.
  DB = false;

  // Temporary:
  //  F = zeros(1,100);
  if (ch1_size[0] >= 250) {
    //  Filter using optimized EOG filter:
    //  (take last second, regardless of actual window length):
    eogcfilt(*(double (*)[250])&ch1_data[ch1_size[0] - 250], ch1f);
    eogcfilt(*(double (*)[250])&ch2_data[ch2_size[0] - 250], ch2f);
    eogcfilt(*(double (*)[250])&ch3_data[ch3_size[0] - 250], ch3f);
    eogcfilt(*(double (*)[250])&ch4_data[ch4_size[0] - 250], ch4f);

    // Extract EOG features: (1s window)
    featureExtractionEOG(ch1f, tmp_data, tmp_size);
    featureExtractionEOG(ch2f, b_tmp_data, tmp_size);
    featureExtractionEOG(ch3f, c_tmp_data, tmp_size);
    featureExtractionEOG(ch4f, d_tmp_data, tmp_size);

    // Combine features:
    // Boolean DB: represents presence of a double blink.
    for (mtmp = 0; mtmp < 10; mtmp++) {
      dv0[mtmp] = tmp_data[mtmp];
      dv0[mtmp + 10] = b_tmp_data[mtmp];
      dv0[mtmp + 20] = c_tmp_data[mtmp];
      dv0[mtmp + 30] = d_tmp_data[mtmp];
    }

    Y = knn(dv0, tXEOG, tYEOG);
    if (Y == 1.0) {
      DB = true;
    }

    //  SSVEP CLASSIFICATION:
    //  PRECONDITIONS: EOG must not have been triggered.
    //  starts with 1/2 second analysis and moves up.
    //  Output can be one of the four SSVEP classes [10 12 15 16]
    if (!DB) {
      // If no double blink has been detected in final second of data.
      //  Use a decision tree.
      // Y = knn(tsX, tX, tY, 1); %Fine KNN
      // limit size:
      varargin_1[0] = (short)ch1_size[0];
      varargin_1[1] = (short)ch2_size[0];
      varargin_1[2] = (short)ch3_size[0];
      mtmp = (short)ch1_size[0];
      for (ix = 1; ix + 1 < 4; ix++) {
        if (varargin_1[ix] < mtmp) {
          mtmp = varargin_1[ix];
        }
      }

      ln = mtmp;

      //  Make sure len is even:
      if (b_mod((double)mtmp) != 0.0) {
        ln = mtmp - 1;
      }

      if (1 > ln) {
        loop_ub = 0;
      } else {
        loop_ub = ln;
      }

      if (1 > ln) {
        b_loop_ub = 0;
      } else {
        b_loop_ub = ln;
      }

      if (1 > ln) {
        ix = 0;
      } else {
        ix = ln;
      }

      // Filter:
      b_ch1_size[0] = loop_ub;
      for (mtmp = 0; mtmp < loop_ub; mtmp++) {
        b_ch1_data[mtmp] = ch1_data[mtmp];
      }

      emxInit_real_T(&b_ch1f, 2);
      eegcfilt(b_ch1_data, b_ch1_size, b_ch1f);
      b_ch2_size[0] = b_loop_ub;
      for (mtmp = 0; mtmp < b_loop_ub; mtmp++) {
        b_ch1_data[mtmp] = ch2_data[mtmp];
      }

      emxInit_real_T(&b_ch2f, 2);
      eegcfilt(b_ch1_data, b_ch2_size, b_ch2f);
      b_ch3_size[0] = ix;
      for (mtmp = 0; mtmp < ix; mtmp++) {
        b_ch1_data[mtmp] = ch3_data[mtmp];
      }

      emxInit_real_T(&b_ch3f, 2);
      eegcfilt(b_ch1_data, b_ch3_size, b_ch3f);

      //  Extract SSVEP Features (Part 1 from individual channels):
      featureExtractionSSVEP(b_ch1f, b_ch2f, b_ch3f, Fs, F);

      // Different window lengths correspond to different classifiers.
      emxFree_real_T(&b_ch3f);
      emxFree_real_T(&b_ch2f);
      emxFree_real_T(&b_ch1f);
      if (ch1_size[0] < 500) {
        // numFeats <= 53
        //  window length is less than 500 samples
        //  Verify # of features:
        //              if size(F,1) == 42
        //              end
        Y = b_knn(F, tXSSVEP, tYSSVEP);
      } else {
        // numFeats>53
        //  Window length is 500 or more samples.
        //              Y = knn(F,tXSSVEPlong,tYSSVEPlong,5);
      }

      //  Analysis: Use Tree-based classification or CCA-FKNN:
      // %% ^ TODO: use self-programmed decision tree for l = 250;
      // %% ^ Use CCA-FKNN for larger datasets (for l>=500);
      //          numFeatures = length(F);
      //  Decisions
      //          if ~SSVEP_PRESENT %TEMP (obviously)
      //              Y = 0;
      //          end
    }
  }

  return Y;
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
