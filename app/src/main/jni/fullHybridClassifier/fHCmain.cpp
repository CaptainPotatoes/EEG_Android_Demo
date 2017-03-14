//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: main.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 14-Mar-2017 00:59:20
//

//***********************************************************************
// This automatically generated example C main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************
// Include Files
#include "rt_nonfinite.h"
#include "fullHybridClassifier.h"
#include "fHCmain.h"
/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define LOG_TAG "fHCmain-cpp"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)
// Function Declarations
static emxArray_real_T *argInit_Unboundedx1_real_T();
static emxArray_real_T *argInit_Unboundedx40_real_T();
static emxArray_real_T *argInit_Unboundedx42_real_T();
static void argInit_d1250x1_real_T(double result_data[], int result_size[1]);
static double argInit_real_T();
static void main_fullHybridClassifier();

// Function Definitions

//
// Arguments    : void
// Return Type  : emxArray_real_T *
//
static emxArray_real_T *argInit_Unboundedx1_real_T()
{
  emxArray_real_T *result;
  static int iv1[1] = { 2 };

  int idx0;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreateND_real_T(1, *(int (*)[1])&iv1[0]);

  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result->data[idx0] = argInit_real_T();
  }

  return result;
}

//
// Arguments    : void
// Return Type  : emxArray_real_T *
//
static emxArray_real_T *argInit_Unboundedx40_real_T()
{
  emxArray_real_T *result;
  static int iv0[2] = { 2, 40 };

  int idx0;
  int idx1;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreateND_real_T(2, *(int (*)[2])&iv0[0]);

  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < 40; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

//
// Arguments    : void
// Return Type  : emxArray_real_T *
//
static emxArray_real_T *argInit_Unboundedx42_real_T()
{
  emxArray_real_T *result;
    //TODO: CHANGE: '2' in following func to initialize correct training data size.
  static int iv2[2] = { 2, 42 };

  int idx0;
  int idx1;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreateND_real_T(2, *(int (*)[2])&iv2[0]);

  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < 42; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

//
// Arguments    : double result_data[]
//                int result_size[1]
// Return Type  : void
//
static void argInit_d1250x1_real_T(double result_data[], int result_size[1])
{
  int idx0;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result_size[0] = 2;

  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < 2; idx0++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result_data[idx0] = argInit_real_T();
  }
}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_fullHybridClassifier()
{
  double ch1_data[1250];
  int ch1_size[1];
  double ch2_data[1250];
  int ch2_size[1];
  double ch3_data[1250];
  int ch3_size[1];
  double ch4_data[1250];
  int ch4_size[1];
  emxArray_real_T *tXEOG;
  emxArray_real_T *tYEOG;
  emxArray_real_T *tXSSVEP;
  emxArray_real_T *tYSSVEP;
  double Y;

  // Initialize function 'fullHybridClassifier' input arguments.
  // Initialize function input argument 'ch1'.
  argInit_d1250x1_real_T(ch1_data, ch1_size);

  // Initialize function input argument 'ch2'.
  argInit_d1250x1_real_T(ch2_data, ch2_size);

  // Initialize function input argument 'ch3'.
  argInit_d1250x1_real_T(ch3_data, ch3_size);

  // Initialize function input argument 'ch4'.
  argInit_d1250x1_real_T(ch4_data, ch4_size);

  // Initialize function input argument 'tXEOG'.
  tXEOG = argInit_Unboundedx40_real_T();

  // Initialize function input argument 'tYEOG'.
  tYEOG = argInit_Unboundedx1_real_T();

  // Initialize function input argument 'tXSSVEP'.
  tXSSVEP = argInit_Unboundedx42_real_T();

  // Initialize function input argument 'tYSSVEP'.
  tYSSVEP = argInit_Unboundedx1_real_T();

  // Call the entry-point 'fullHybridClassifier'.
  Y = fullHybridClassifier(ch1_data, ch1_size, ch2_data, ch2_size, ch3_data,
    ch3_size, ch4_data, ch4_size, tXEOG, tYEOG, tXSSVEP, tYSSVEP, argInit_real_T
    ());
  emxDestroyArray_real_T(tYSSVEP);
  emxDestroyArray_real_T(tXSSVEP);
  emxDestroyArray_real_T(tYEOG);
  emxDestroyArray_real_T(tXEOG);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // Initialize the application.
  // You do not need to do this more than one time.
  fullHybridClassifier_initialize();

  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_fullHybridClassifier();

  // Terminate the application.
  // You do not need to do this more than one time.
  fullHybridClassifier_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
