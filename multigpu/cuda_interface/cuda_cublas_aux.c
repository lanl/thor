/*--------------------------------------------------------------------------~*
 * Copyright (c) 2025 Triad National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~*/

/*
 * @file cuda_cublas_aux.c
 * @author Oleg Korobkin
 * @date  May 2025
 * @brief Convenience C wrappers for cuBLAS library functions
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <cublas_v2.h>

#include <mpi.h>
/*
#ifdef USE_CAL_MPI
    #include <cal_mpi.h>
#endif
*/
#include <cal.h>

int set_cublasHandle(cublasHandle_t *handle, int VRBZ){
    if(VRBZ>0) printf("[+][IN][set_cublasHandle]\n");
    cublasStatus_t cublasStat;
    cublasStat = cublasCreate(handle);
    if(VRBZ>0) printf("[+][set_cublasHandle]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


int destroy_cublasHandle(cublasHandle_t *handle, int VRBZ){
    if(VRBZ>0) printf("[+][IN][destroy_cublasHandle]\n");
    cublasStatus_t cublasStat;
    cublasStat = cublasDestroy(*handle);
    if(VRBZ>0) printf("[+][destroy_cublasHandle]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


// Returns index of the maximum element in an array (base 1)
int compute_cublasIdamax(cublasHandle_t *handle, int n,
                         const double *x, int incx, int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasIdamax]\n");
    cublasStatus_t cublasStat;
    int result;
    cublasStat = cublasIdamax(*handle, n, x, incx, &result);
    if(VRBZ>0) printf("[+][compute_cublasIdamax]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return result;
}


// Returns index of the minimum element in an array (base 1)
int compute_cublasIdamin(cublasHandle_t *handle, int n,
                         const double *x, int incx, int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasIdamin]\n");
    cublasStatus_t cublasStat;
    int result;
    cublasStat = cublasIdamin(*handle, n, x, incx, &result);
    if(VRBZ>0) printf("[+][compute_cublasIdamin]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return result;
}


// Computes the sum of the absolute values of the elements of vector
double compute_cublasDasum(cublasHandle_t *handle, int n,
                           const double *x, int incx, int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDasum]\n");
    cublasStatus_t cublasStat;
    double result;
    cublasStat = cublasDasum(*handle, n, x, incx, &result);
    if(VRBZ>0) printf("[+][compute_cublasDasum]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return result;
}


// Computes a*x + y -> y
int compute_cublasDaxpy(cublasHandle_t *handle, int n,
                        const double alpha,
                        const double *x, int incx,
                        double       *y, int incy,
                        int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDaxpy]\n");
    cublasStatus_t cublasStat;
    cublasStat = cublasDaxpy(*handle, n, &alpha, x, incx, y, incy);
    if(VRBZ>0) printf("[+][compute_cublasDaxpy]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


// Copies x -> y
int compute_cublasDcopy(cublasHandle_t *handle, int n,
                        const double *x, int incx,
                        double       *y, int incy,
                        int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDcopy]\n");
    cublasStatus_t cublasStat;
    cublasStat = cublasDcopy(*handle, n, x, incx, y, incy);
    if(VRBZ>0) printf("[+][compute_cublasDcopy]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


// Computes a dot product (x, y)
double compute_cublasDdot(cublasHandle_t *handle, int n,
                          const double *x, int incx,
                          double       *y, int incy,
                          int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDdot]\n");
    cublasStatus_t cublasStat;
    double result;
    cublasStat = cublasDdot(*handle, n, x, incx, y, incy, &result);
    if(VRBZ>0) printf("[+][compute_cublasDdot]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return result;
}


// Computes the Euclidean norm of the vector x (= sqrt(sum[|x_i|^2]))
double compute_cublasDnrm2(cublasHandle_t *handle, int n,
                           const double *x, int incx, int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDnrm2]\n");
    cublasStatus_t cublasStat;
    double result;
    cublasStat = cublasDnrm2(*handle, n, x, incx, &result);
    if(VRBZ>0) printf("[+][compute_cublasDnrm2]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return result;
}


// Applies Givens rotation matrix to a pair of vectors {x, y}
int compute_cublasDrot(cublasHandle_t *handle, int n,
                       double *x, int incx,
                       double *y, int incy,
                       const double c, const double s,
                       int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDrot]\n");
    cublasStatus_t cublasStat;
    cublasStat = cublasDrot(*handle, n, x, incx, y, incy, &c, &s);
    if(VRBZ>0) printf("[+][compute_cublasDrot]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


// Scales a vector by the specified factor
int compute_cublasDscal(cublasHandle_t *handle, int n,
                        const double alpha,
                        double *x, int incx,
                        int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDscal]\n");
    cublasStatus_t cublasStat;
    cublasStat = cublasDscal(*handle, n, &alpha, x, incx);
    if(VRBZ>0) printf("[+][compute_cublasDscal]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


// Swaps two vectors x <-> y
int compute_cublasDswap(cublasHandle_t *handle, int n,
                        double *x, int incx,
                        double *y, int incy,
                        int VRBZ) {
    if(VRBZ>0) printf("[+][IN][compute_cublasDswap");
    cublasStatus_t cublasStat;
    cublasStat = cublasDswap(*handle, n, x, incx, y, incy);
    if(VRBZ>0) printf("[+][compute_cublasDswap]: stat=%d\n", (int)cublasStat);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return (int)cublasStat;
}


