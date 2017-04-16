//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: main.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 22-Mar-2017 11:10:49
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
#include "fhcmain.h"
/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define  LOG_TAG "fullHybridClassifierMain-cpp"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

// Function Declarations
static boolean_T argInit_boolean_T();
static emxArray_real_T *argInit_d4996x1_real_T();
static double argInit_real_T();
static void main_fullHybridClassifier();

// Function Definitions

//
// Arguments    : void
// Return Type  : boolean_T
//
static boolean_T argInit_boolean_T()
{
  return false;
}

//
// Arguments    : void
// Return Type  : emxArray_real_T *
//
static emxArray_real_T *argInit_d4996x1_real_T()
{
  emxArray_real_T *result;
  static int iv6[1] = { 2 };

  int idx0;

  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreateND_real_T(1, *(int (*)[1])&iv6[0]);

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
  emxArray_real_T *ch1;
  emxArray_real_T *ch2;
  emxArray_real_T *ch3;
  emxArray_real_T *ch4;
  double Y[7];

  // Initialize function 'fullHybridClassifier' input arguments.
  // Initialize function input argument 'ch1'.
  ch1 = argInit_d4996x1_real_T();

  // Initialize function input argument 'ch2'.
  ch2 = argInit_d4996x1_real_T();

  // Initialize function input argument 'ch3'.
  ch3 = argInit_d4996x1_real_T();

  // Initialize function input argument 'ch4'.
  ch4 = argInit_d4996x1_real_T();

  // Call the entry-point 'fullHybridClassifier'.
   LOGE("RUNNING FHC STARTUP TEST!");
  fullHybridClassifier(ch1, ch2, ch3, ch4, argInit_real_T(), argInit_boolean_T(),
                       Y);
    LOGE("Y:[ %f %f %f %f %f %f %f ]",Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6]);
  emxDestroyArray_real_T(ch4);
  emxDestroyArray_real_T(ch3);
  emxDestroyArray_real_T(ch2);
  emxDestroyArray_real_T(ch1);
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

extern "C" {
JNIEXPORT jint JNICALL
Java_com_mahmoodms_bluetooth_eegssvepdemo_DeviceControlActivity_jmainFHC(
        JNIEnv *env, jobject obj, jboolean terminate) {
  if(!(bool)terminate) {
    fullHybridClassifier_initialize();
    main_fullHybridClassifier();
//    fir_combined_initialize();
//    main_fir_combined();
    return 0;
  } else {
      fullHybridClassifier_terminate();
//    fir_combined_terminate();
    return -1;
  }
}
}


//
// File trailer for main.cpp
//
// [EOF]
//
