//JNI Code:

//TODO: NEED TO INITIALIZE IN fHCmain.cpp

#include <cstdlib>
/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define LOG_TAG "fullHybridClassifier-cpp"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

//TODO: Will need to pass all data (4 arrays), and length of array.
    // How do I get array length in C++?

extern "C" {
JNIEXPORT jdoubleArray JNICALL
//Call this function with (data, data, data, data, datalen, Fs);
Java_com_mahmoodms_bluetooth_ecgfallsensordemo_DeviceControlActivity_jfullHybridClassifier(
        JNIEnv *env, jobject jobject1, jdoubleArray array1, jdoubleArray array2, jdoubleArray array3, jdoubleArray array4, jdouble dataLength) {
    jdouble  *c_array;
    c_array = env->GetDoubleArrayElements(array1, NULL);
    if (c_array==NULL) LOGE("ERROR - C_ARRAY IS NULL");
    int len = (int) dataLength;
    double r_array[len];
//    fir_combined(c_array, r_array);
    double *result = r_array;
    jdoubleArray m_result;
    m_result = env->NewDoubleArray(len);
    env->SetDoubleArrayRegion(m_result, 0, len, result);
    return m_result;
}
}