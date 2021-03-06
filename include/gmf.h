/***********************************************************************************************************************************
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME : gmf_extension.h

PURPOSE : Declares the functions for the GMF filter method implemented in c++ for its use within python.

FILE REFERENCES :
Name        I / O        Description
None----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES :
None
************************************************************************************************************************************/
#ifndef GMF_DLL_H_INCLUDED
#define GMF_DLL_H_INCLUDED

#ifdef BUILDING_PYTHON_MODULE
    #include <Python.h>
    #include <numpy/ndarraytypes.h>
    #include <numpy/ufuncobject.h>
    #include <numpy/npy_3kcompat.h>
    #define GMF_DLL_PUBLIC
    #define GMF_DLL_LOCAL
#else
    #if defined(_WIN32) || defined(_WIN64)
        #ifdef BUILDING_GMF_DLL
            #ifdef __GNUC__
                #define GMF_DLL_PUBLIC __attribute__ ((dllexport))
            #else
                #define GMF_DLL_PUBLIC __declspec(dllexport)
            #endif
        #else
            #ifdef __GNUC__
                #define GMF_DLL_PUBLIC __attribute__ ((dllimport))
            #else
                #define GMF_DLL_PUBLIC __declspec(dllimport)
            #endif
        #endif
        #define GMF_DLL_LOCAL
    #else
        #if __GNUC__ >= 4
            #define GMF_DLL_PUBLIC __attribute__ ((visibility ("default")))
            #define GMF_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
        #else
            #define GMF_DLL_PUBLIC
            #define GMF_DLL_LOCAL
        #endif
    #endif
#endif
#ifdef COMPILE_ACML
#include <acml.h>
#define ft_complex doublecomplex
#define ft_real(POINTER) (POINTER)->real
#define ft_imag(POINTER) (POINTER)->imag
#define ft_real_assign(POINTER) (POINTER)->real
#define ft_imag_assign(POINTER) (POINTER)->imag
#define allocate_ft_complex(ARRAY_SIZE) (doublecomplex*)malloc(ARRAY_SIZE*sizeof(doublecomplex))
#define allocate_ft_complex_pointers(ARRAY_SIZE) (doublecomplex**)malloc(ARRAY_SIZE*sizeof(doublecomplex*))
#define allocate_ft_complex_pointers_of_pointers(ARRAY_SIZE) (doublecomplex***)malloc(ARRAY_SIZE*sizeof(doublecomplex**))
#define deallocate_ft_complex(COMPLEX_ARRAY) free(COMPLEX_ARRAY)
#define ft_variables(HEIGHT, WIDTH) double* communication_work_array = (double*)malloc((4*WIDTH + 6*HEIGHT + 300) * sizeof(double)); int information_integer
#define ft_forward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) dzfft2d(0, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_backward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) zdfft2d(0, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_forward(HEIGHT, WIDTH, INPUT, OUTPUT) dzfft2d(1, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_backward(HEIGHT, WIDTH, INPUT, OUTPUT) zdfft2d(1, HEIGHT, WIDTH, INPUT, OUTPUT, communication_work_array, &information_integer)
#define ft_release_forward 
#define ft_release_backward 
#define ft_close free(communication_work_array)
#else
#include <fftw3.h>
#define ft_complex fftw_complex
#define ft_real(POINTER) *(*(POINTER))
#define ft_imag(POINTER) *(*(POINTER) + 1)
#define ft_real_assign(POINTER) *(*(POINTER))
#define ft_imag_assign(POINTER) *(*(POINTER) + 1)
#define allocate_ft_complex(ARRAY_SIZE) (fftw_complex*)fftw_malloc(ARRAY_SIZE*sizeof(fftw_complex))
#define allocate_ft_complex_pointers(ARRAY_SIZE) (fftw_complex**)fftw_malloc(ARRAY_SIZE*sizeof(fftw_complex*))
#define allocate_ft_complex_pointers_of_pointers(ARRAY_SIZE) (fftw_complex***)fftw_malloc(ARRAY_SIZE*sizeof(fftw_complex**))
#define deallocate_ft_complex(COMPLEX_ARRAY) fftw_free(COMPLEX_ARRAY)
#define ft_variables(HEIGHT, WIDTH) fftw_plan forward_plan, backward_plan; char forward_plan_active = 0, backward_plan_active = 0;
#define ft_forward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) forward_plan = fftw_plan_dft_r2c_2d(HEIGHT, WIDTH, INPUT, OUTPUT,  FFTW_ESTIMATE); forward_plan_active = 1
#define ft_backward_setup(HEIGHT, WIDTH, INPUT, OUTPUT) backward_plan = fftw_plan_dft_c2r_2d(HEIGHT, WIDTH, INPUT, OUTPUT, FFTW_ESTIMATE); backward_plan_active = 1
#define ft_forward(HEIGHT, WIDTH, INPUT, OUTPUT) fftw_execute(forward_plan)
#define ft_backward(HEIGHT, WIDTH, INPUT, OUTPUT) fftw_execute(backward_plan)
#define ft_release_forward fftw_destroy_plan(forward_plan); forward_plan_active = 0
#define ft_release_backward fftw_destroy_plan(backward_plan); backward_plan_active = 0
#define ft_close if (forward_plan_active){fftw_destroy_plan(forward_plan);} if (backward_plan_active){fftw_destroy_plan(backward_plan);}
#endif


#ifndef NDEBUG
#define DEBMSG(MESSAGE) printf(MESSAGE)
#define DEBNUMMSG(MESSAGE, NUM) printf(MESSAGE, NUM);
#else
#define DEBMSG(MESSAGE) 
#define DEBNUMMSG(MESSAGE, NUM) 
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MY_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164

#define MY_INF 179769313486231570814527423731704356798070567525844996598917476803157260780028538760589558632766878171540458953514382464234321326889464182768467546703537516986049910576551282076245490090389328944075868508455133942304583236903222948165808559332123348274797826204144723168738177180919299881250404026184124858368.0


/* C implementations: */
static inline double bicubicInterpolation_impl(double* src, const double x, const double y, const unsigned int src_height, const unsigned int src_width);

void GMF_DLL_LOCAL rotateBicubic_impl(double* src, double* dst, const double theta, const int src_height, const int src_width);

double * GMF_DLL_LOCAL generateTemplate(double * template_src, const int par_T, const int par_L, unsigned int *kernel_height, unsigned int *kernel_width, const unsigned char untrimmed_kernels);
double * GMF_DLL_LOCAL generateGMFTemplate(const int par_T, const int par_L, const double par_sigma, unsigned int *kernel_height, unsigned int *kernel_width, const unsigned char untrimmed_kernels);
double **  GMF_DLL_LOCAL generateFilterbank_space(double * base_kernel, const int par_K, const int nearest_2p_dim, const unsigned int kernel_height, const unsigned int kernel_width);

ft_complex ** GMF_DLL_LOCAL generateFilterbank(double * base_kernel, const int par_K, const int nearest_2p_dim, const unsigned int kernel_height, const unsigned int kernel_width);

void GMF_DLL_LOCAL applyFilterAngles(ft_complex * fft_img_src, double* img_dst, double* ang_dst, const int nearest_2p_dim, const int height, const int width, const unsigned int s_scales, const int kernel_height, const int kernel_width, const int par_K, ft_complex ** fft_filter_bank);
void GMF_DLL_LOCAL multiscaleFilterAngles(double * raw_input, double * output, double * angles, const unsigned int n_inputs, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels, double ** template_src);

void GMF_DLL_LOCAL applyFilter(ft_complex * fft_img_src, double* img_dst, const int nearest_2p_dim, const int height, const int width, const int kernel_height, const int kernel_width, const int par_K, ft_complex ** fft_filter_bank);
void GMF_DLL_LOCAL multiscaleFilter(double * raw_input, double * output, const unsigned int n_inputs, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels, double ** template_src);

#ifdef BUILDING_PYTHON_MODULE
static PyObject* gmfFilter(PyObject *self, PyObject *args);
#endif

#ifdef BUILDING_PYTHON_MODULE
static PyObject* gmfFilterWithAngles(PyObject *self, PyObject *args);
#endif

#ifdef BUILDING_PYTHON_MODULE
static PyObject* gmfFilterBank(PyObject *self, PyObject *args);
#endif

#endif //GMF_DLL_H_INCLUDED
