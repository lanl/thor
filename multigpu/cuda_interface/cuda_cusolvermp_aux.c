/*--------------------------------------------------------------------------~*
 * Copyright (c) 2025 Triad National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~*/

/*
 * @file cuda_cusolvermp_aux.c
 * @author Ismael Djibrilla Boureima, Oleg Korobkin
 * @date  December 2023
 * @brief Convenience C wrappers for cuSolverMp library functions
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <cusolverMp.h>
#include <mpi.h>
/*
#ifdef USE_CAL_MPI
    #include <cal_mpi.h>
#endif
*/
#include <nccl.h>

// MODULE PROTOTYPES
static void print_host_matrix(int64_t M, int64_t N, double* A, int64_t lda, const char* msg);
static inline int getLocalRank();


int show_cusolverMpHandle(cusolverMpHandle_t *cusolverMpHandle){
  printf("[+][show_cusolverMpHandle]: cusolverMpHandle=%x \n", cusolverMpHandle);
  return 0;
}


int show_cusolverMpGrid(cusolverMpGrid_t *cusolverMpGrid){
  printf("[+][show_cusolverMpGrid]: cusolverMpGrid=%x \n", cusolverMpGrid);
  return 0;
}

int show_cusolverMpMatrixDescriptor(cusolverMpMatrixDescriptor_t *cusolverMpMatrixDescriptor){
  printf("[+][show_cusolverMpMatrixDescriptor]: cusolverMpMatrixDescriptor=%x \n", cusolverMpMatrixDescriptor);
  return 0;
}

int show_device_int_array(int *d_array, int size_array){
    int *h_array;
    int DSIZE, dt;
    dt =1;
    DSIZE = size_array * sizeof(dt);
    h_array = (int *)malloc(DSIZE);
    cudaMemcpy(h_array, d_array, DSIZE, cudaMemcpyDeviceToHost);
    printf("[+][show_device_int_array]: device_int_array={");
    for (int i=0; i< size_array; i++){
       printf("%d ", h_array[i]);
    }
    printf("} \n");
    return 0;
}

int show_device_real_array(float *d_array, int size_array){
    float *h_array;
    float dt;
    dt =1.0;
    int DSIZE = size_array * sizeof(dt);
    h_array = (float *)malloc(DSIZE);
    cudaMemcpy(h_array, d_array, DSIZE, cudaMemcpyDeviceToHost);
    printf("[+][show_device_real_array]: device_real_array={");
    for (int i=0; i< size_array; i++){
       printf("%12.5f ", h_array[i]);
    }
    printf("} \n");
    return 0;
}

int show_device_double_array(double *d_array, int size_array){
    double *h_array;
    double dt =1.0;
    int DSIZE = size_array * sizeof(dt);
    h_array = (double *)malloc(DSIZE);
    cudaMemcpy(h_array, d_array, DSIZE, cudaMemcpyDeviceToHost);
    printf("[+][show_device_double_array]: device_double_array={");
    for (int i=0; i< size_array; i++){
       printf("%12.5f ", h_array[i]);
    }
    printf("} \n");
    return 0;
}


// Print matrix
static void print_host_matrix(int64_t M, int64_t N, double* A, int64_t lda, const char* msg)
{
    if (M * N > 2000) return;
    printf("print_host_matrix : %s\n", msg);

    for (int64_t i = 0; i < M; i++)
    {
        for (int64_t j = 0; j < N; j++)
        {
            printf("%.2lf  ", A[i + j * lda]);
        }
        printf("\n");
    }
}

static inline int getLocalRank()
{
    int      localRank;
    MPI_Comm localComm;

    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &localComm);
    MPI_Comm_rank(localComm, &localRank);
    MPI_Comm_free(&localComm);

    return localRank;
}


int set_cusolverMpHandle(cusolverMpHandle_t *cusolverMpHandle, int localDeviceId, cudaStream_t *localStream, int VRBZ){
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    /* Initialize cusolverMp library handle */
    if(VRBZ>0) printf("[+][set_cusolverMpHandle()][INPUTs]: localDeviceId=%d | localStream=%d  \n", localDeviceId, *localStream);
    cusolverStat = cusolverMpCreate(cusolverMpHandle, localDeviceId, *localStream);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if(VRBZ>0) {
        if(cusolverStat==CUSOLVER_STATUS_SUCCESS)
            printf("[+][set_cusolverMpHandle]: cusolverMpHandle INITIALIZED OK: cusolverStat=%x \n", cusolverStat);
        else
            printf("[!][set_cusolverMpHandle][ERR]: UNABLE to INITIALIZE: cusolverStat=%x \n", cusolverStat);
    }
    return 0;
}

int destroy_cusolverMpHandle(cusolverMpHandle_t *cusolverMpHandle, int VRBZ){
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    /* Destroy cusolverMp handle */
    if(VRBZ>0) printf("[+][IN][destroy_cusolverMpHandle()]\n");
    cusolverStat = cusolverMpDestroy(*cusolverMpHandle);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][destroy_cusolverMpHandle()]: cusolverMpHandle DESTROYED OK: cusolverStat=%x \n", cusolverStat);
    return 0;
}

// cusolverMp grids
int set_cusolverMpGrid(cusolverMpHandle_t *cusolverMpHandle, cusolverMpGrid_t *gridA, ncclComm_t *comm,
                       int64_t numRowDevices, int64_t numColDevices, int VRBZ){
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    cusolverStat = cusolverMpCreateDeviceGrid(*cusolverMpHandle, gridA, *comm, (int32_t)numRowDevices, (int32_t)numColDevices,
                                               CUSOLVERMP_GRID_MAPPING_COL_MAJOR);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][set_cusolverMpGrid()]: cusolverMpGrid CREATED OK: cusolverStat=%x \n", cusolverStat);
    if(VRBZ>0) printf("[+][OUT][set_cusolverMpGrid()] numRowDevices=%d | numColDevices=%d| comm= %x  \n",
                                                     numRowDevices, numColDevices, comm);
    return 0;

}

int destroy_cusolverMpGrid(cusolverMpGrid_t *gridA, int VRBZ){
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    /* Destroy cusolverMpGrid */
    if(VRBZ>0) printf("[+][IN][destroy_cusolverMpGrid()]\n");
    cusolverStat = cusolverMpDestroyGrid(*gridA);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][destroy_cusolverMpGrid()]: cusolverMpGrid DESTROYED OK: cusolverStat=%x \n", cusolverStat);
    return 0;
}

// cusolverMp matrix descriptor
int set_cusolverMpMatrixDesc(cusolverMpMatrixDescriptor_t *descrA,
          cusolverMpGrid_t *gridA, int64_t M, int64_t N,
          int64_t MA, int64_t NA,  int64_t numRowDevices, int dev_row_idx, int VRBZ){
    if(VRBZ>0) printf("[+][IN][set_cusolverMpMatrixDesc()] M=%d| N=%d| MA=%d| NA=%d| numRowDevices=%d| dev_row_idx=%d \n",
                               M,N,MA,NA,numRowDevices, dev_row_idx);
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    const uint32_t RSRCA = 0;
    const uint32_t CSRCA = 0;
    const int64_t LLDA   = cusolverMpNUMROC(M, MA, dev_row_idx, RSRCA, numRowDevices);
    if(VRBZ>0) printf("[+][IN][set_cusolverMpMatrixDesc()]: llda=%d|\n", LLDA);

    // descr    Host  Out Matrix descriptor object initialized by this function.
    // dataType Host  In  Data type of the matrix A.
    // M_A      Host  In  Number of rows in the global matrix A.
    // N_A      Host  In  Number of columns in the global matrix A.
    // MB_A     Host  In  Blocking factor used to distribute the rows of the global matrix A.
    // NB_A     Host  In  Blocking factor used to distribute the columns of the global matrix A.
    // RSRC_A   Host  In  Process row over which the first row of the matrix A is distributed. Only the value of 0 is currently supported.
    // CSRC_A   Host  In  Process column over which the first row of the matrix A is distributed. Only the value of 0 is currently supported.
    // LLD_A    Host  In  Leading dimension of the local matrix.
    cusolverStat = cusolverMpCreateMatrixDesc(descrA, *gridA, CUDA_R_64F, M, N, MA, NA, RSRCA, CSRCA, LLDA);
    if(VRBZ>0) printf("[+][set_cusolverMpMatrixDesc()]: cusolverMpMatrixDescriptor  CREATED OK: cusolverStat=%x \n", cusolverStat);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    return 0;
}

int destroy_cusolverMpMatrixDesc(cusolverMpMatrixDescriptor_t *descrA, int VRBZ){
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    /* Destroy cusolverMpMatrixDescriptor */
    if(VRBZ>0) printf("[+][IN][destroy_cusolverMpMatrixDesc()]\n");
    cusolverStat = cusolverMpDestroyMatrixDesc(*descrA);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][destroy_cusolverMpMatrixDesc()]: cusolverMpMatrixDescriptor DESTROYED OK: cusolverStat=%x \n", cusolverStat);
    return 0;
}



// Compute local matrix dimensions: ML and NL
int set_dist_mat_row_col_params(int64_t M, int64_t N, int64_t MA, int64_t NA,
                    int rankId, int64_t numRowDevices, int64_t numColDevices,
                    int64_t *ML, int64_t *NL, int VRBZ) {
    const uint32_t RSRCA = 0;
    if(VRBZ>0) printf("[+][rank%d][set_dist_mat_row_col_params()][INPUTs]:"
                      " M=%d| N=%d| MA=%d| NA=%d |"
                      " numRowDevices=%d| numColDevices=%d|\n",
                      rankId,M,N,MA,NA,numRowDevices,numColDevices);
    *ML = M;
    *NL = N;

    if (numColDevices > 1) {
       // Current implementation only allows RSRC,CSRC = 0
       const uint32_t CSRCA = 0;
       *NL = (int)(cusolverMpNUMROC(N, NA, rankId, CSRCA, numColDevices));
    }
    if (numRowDevices > 1) {
       const uint32_t RSRCA = 0;
       *ML = (int)(cusolverMpNUMROC(M, MA, rankId, RSRCA, numRowDevices));
    }
    if(VRBZ>0) printf("[+][rank%d][set_dist_mat_row_col_params()][OUTPUT]:"
                      " ML=%d| NL=%d\n", rankId, *ML, *NL);
    return 0;
}


int check_qr_buffsize(cusolverMpHandle_t *cusolverMpHandle,
                      const int64_t M, const int64_t N,
                      const int64_t IA, const int64_t JA,
                      cusolverMpMatrixDescriptor_t *descrA, double* d_A){
    printf("[+][check_qr_buffsize][INPUTs] cusolverMpHandle=%x| M=%d,N=%d,IA=%d,JA=%d | descrA=%x| d_A=%x\n", cusolverMpHandle, M,N,IA,JA, descrA, d_A);

    /* Error codes */
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    /* size of workspace on device */
    size_t workspaceInBytesOnDevice_geqrf = 0;
    /* size of workspace on host */
    size_t workspaceInBytesOnHost_geqrf = 0;
    cusolverStat = cusolverMpGeqrf_bufferSize(*cusolverMpHandle,
                                              M,
                                              N,
                                              d_A,
                                              IA,
                                              JA,
                                              *descrA,
                                              CUDA_R_64F,
                                              &workspaceInBytesOnDevice_geqrf,
                                              &workspaceInBytesOnHost_geqrf);
    printf("[+][check_qr_buffsize][STAT][cusolverMpGeqrf_bufferSize]: cusolverStat = %d \n", cusolverStat);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    return 0;
}

int dist_geqrf(const cudaStream_t * const localStream,
               const cusolverMpHandle_t * const cusolverMpHandle,
               const cusolverMpMatrixDescriptor_t * const descrA,
               double* d_A,
               const int64_t M,
               const int64_t N,
               double* d_tau,
               int VRBZ) {

    // Error codes
    cusolverStatus_t cusolverStat;
    calError_t       calStat;
    cudaError_t      cudaStat;

    cusolverStat = CUSOLVER_STATUS_SUCCESS;
    calStat      = CAL_OK;
    cudaStat     = cudaSuccess;

    // localDeviceId is the deviceId from rank's point of view. This is
    // system-dependent. For example, setting one device per process,
    // Summit always sees the local device as device 0.
    int localDeviceId;
    localDeviceId = getLocalRank();
    if(VRBZ>1)
        printf("[i][%d][dist_geqrf]: localStream=%x | cusolverMpHandle=%x |"
               " descrA=%x | d_A=%x | M = %d | N = %d | d_tau=%x \n", localDeviceId,
               localStream, cusolverMpHandle, descrA, d_A, M, N, d_tau);

    cudaStat = cudaSetDevice(localDeviceId);
    assert(cudaStat == cudaSuccess);
    cudaStat = cudaFree(0);
    assert(cudaStat == cudaSuccess);

    // Distributed device workspace
    void* d_work_geqrf;
    d_work_geqrf = NULL;

    // Distributed host workspace
    void* h_work_geqrf;
    h_work_geqrf = NULL;

    // size of workspace on device
    size_t workspaceInBytesOnDevice_geqrf;
    workspaceInBytesOnDevice_geqrf = 0;

    // size of workspace on host
    size_t workspaceInBytesOnHost_geqrf;
    workspaceInBytesOnHost_geqrf = 0;

    // error codes from cusolverMp (device)
    int* d_info_geqrf;
    d_info_geqrf = NULL;

    // Allocate d_info
    cudaStat = cudaMalloc((void**)&d_info_geqrf, sizeof(int));
    assert(cudaStat == cudaSuccess);

    // Reset d_info
    cudaStat = cudaMemset(d_info_geqrf, 0, sizeof(int));
    assert(cudaStat == cudaSuccess);

    // Query workspace size for mp routines
    if(VRBZ>1)
        printf("[i][%d][dist_geqrf]: calling cusolverMpGeqrf_bufferSize\n", localDeviceId);
    cusolverStat = cusolverMpGeqrf_bufferSize(*cusolverMpHandle,
                                              M,
                                              N,
                                              d_A,
                                              1,
                                              1,
                                              *descrA,
                                              CUDA_R_64F,
                                              &workspaceInBytesOnDevice_geqrf,
                                              &workspaceInBytesOnHost_geqrf);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if(VRBZ>1)
        printf("[i][%d][cusolverMpGeqrf_bufferSize]: workspaceInBytesOnDevice_geqrf=%d | "
               "workspaceInBytesOnHost_geqrf=%d\n", localDeviceId,
               workspaceInBytesOnDevice_geqrf, workspaceInBytesOnHost_geqrf);

    // allocate geqrf workspace
    cudaStat = cudaMalloc((void**)&d_work_geqrf,
                          workspaceInBytesOnDevice_geqrf);
    assert(cudaStat == cudaSuccess);

    h_work_geqrf = (void*)malloc(workspaceInBytesOnHost_geqrf);
    assert(h_work_geqrf != NULL);

    // call geqrf
    if(VRBZ>1)
        printf("[i][%d][dist_geqrf]: calling cusolverMpGeqrf\n", localDeviceId);
    cusolverStat = cusolverMpGeqrf(*cusolverMpHandle,
                                   M,
                                   N,
                                   d_A,
                                   1,
                                   1,
                                   *descrA,
                                   d_tau,
                                   CUDA_R_64F,
                                   d_work_geqrf,
                                   workspaceInBytesOnDevice_geqrf,
                                   h_work_geqrf,
                                   workspaceInBytesOnHost_geqrf,
                                   d_info_geqrf);
    // sync after cusolverMpgeqrf
    cudaStat = cudaStreamSynchronize(*localStream);
    assert(cudaStat == cudaSuccess);

    // deallocate device workspace
    if (d_work_geqrf != NULL) {
        cudaStat = cudaFree(d_work_geqrf);
        assert(cudaStat == cudaSuccess);
        d_work_geqrf = NULL;
    }

    if (d_info_geqrf != NULL) {
        cudaStat = cudaFree(d_info_geqrf);
        assert(cudaStat == cudaSuccess);
        d_info_geqrf = NULL;
    }

    // deallocate host workspace
    if (h_work_geqrf) {
        free(h_work_geqrf);
        h_work_geqrf = NULL;
    }

    // barrier
    MPI_Barrier(MPI_COMM_WORLD);
    if(VRBZ>1) printf("[i][%d][dist_geqrf]: exit\n", localDeviceId);

    return 0;

}

//
// Wrapper for cuSOLVERMp multiply of a matrix C by a matrix Q represented by
// Huseholder reflectors, as returned by cusolverMpGeqrf in matrix A
//
int dist_ormqr(cudaStream_t *localStream,
        cusolverMpHandle_t *cusolverMpHandle,
        cusolverMpMatrixDescriptor_t *descA, // matrix A descriptor
        double* d_A,                         // matrix data
        const int iside,   // which side to multiply:
                           // - 0: multiply on the left
                           // - else: multiply on the right
        const int itrans,  // is A transposed?
                           // - 0: no operation
                           // - else: transpose
        const int64_t M,   // number of rows in global matrix A
        const int64_t N,   // number of columns in global matrix A
        cusolverMpMatrixDescriptor_t *descC, // matrix C descriptor
        double* d_C,                         // matrix data
        const int64_t K,   // Number of Householder reflectors (in A) defining Q
        const double *d_tau,   // Pointer into the local memory to an array of
                               // dimension LOCc(JA+N-1). This array contains
                               // scalar factors of the Householder reflectors
        const int verbose) {   // verbosity flag: 0: silent, 1-..:more verbose

    const cublasSideMode_t side = iside ? CUBLAS_SIDE_RIGHT : CUBLAS_SIDE_LEFT;
    const cublasOperation_t trans = itrans ? CUBLAS_OP_T : CUBLAS_OP_N;

    // Error codes
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    calError_t       calStat      = CAL_OK;
    cudaError_t      cudaStat     = cudaSuccess;

    assert(cudaStat == cudaSuccess);
    cudaStat = cudaFree(0);
    assert(cudaStat == cudaSuccess);

    // Distributed device workspace
    void* d_work_ormqr = NULL;

    // Distributed host workspace
    void* h_work_ormqr = NULL;

    // size of workspace on device
    size_t workspaceInBytesOnDevice_ormqr = 0;

    // size of workspace on host
    size_t workspaceInBytesOnHost_ormqr = 0;

    // error codes from cusolverMp (device)
    int* d_info_ormqr = NULL;

    // Allocate d_info
    cudaStat = cudaMalloc((void**)&d_info_ormqr, sizeof(int));
    assert(cudaStat == cudaSuccess);

    // Reset d_info
    cudaStat = cudaMemset(d_info_ormqr, 0, sizeof(int));
    assert(cudaStat == cudaSuccess);

    if (verbose>1) {
        printf("[i][dist_ormqr] before calling cusolverMpOrmqr_bufferSize\n");
    }
    // Query workspace size for mp routines
    cusolverStat = cusolverMpOrmqr_bufferSize(*cusolverMpHandle,
                                              side,
                                              trans,
                                              M,
                                              N,
                                              K,
                                              d_A,
                                              1,
                                              1,
                                              *descA,
                                              d_tau,
                                              d_C,
                                              1,
                                              1,
                                              *descC,
                                              CUDA_R_64F,
                                              &workspaceInBytesOnDevice_ormqr,
                                              &workspaceInBytesOnHost_ormqr);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);

    if (verbose>1) {
        printf("[i][dist_ormqr] after calling cusolverMpOrmqr_bufferSize:"
                               " workspace on host, device = %d, %d\n",
                               workspaceInBytesOnHost_ormqr, workspaceInBytesOnDevice_ormqr);
    }

    // allocate geqrf workspace
    cudaStat = cudaMalloc((void**)&d_work_ormqr,
                          workspaceInBytesOnDevice_ormqr);
    assert(cudaStat == cudaSuccess);

    h_work_ormqr = (void*)malloc(workspaceInBytesOnHost_ormqr);
    assert(h_work_ormqr != NULL);

    // call ormqr
    if (verbose>1) printf("[i][dist_ormqr] before calling cusolverMpOrmqr\n");
    if (verbose>2) {
        printf("[d][dist_ormqr] input parameters for cusolverMpOrmqr:\n"
               "[d][dist_ormqr] side = %d, trans = %d,\n"
               "[d][dist_ormqr] M=%d, N=%d, K=%d, d_A=%x, d_tau=%x, d_C=%x\n"
               "[d][dist_ormqr] wDM=%d, wH=%d, d_work=%x, h_work=%x\n",
               side, trans, M, N, K, d_A, d_tau, d_C,
               workspaceInBytesOnDevice_ormqr,
               workspaceInBytesOnHost_ormqr,
               d_work_ormqr, h_work_ormqr);
    }
    cusolverStat = cusolverMpOrmqr(*cusolverMpHandle,
                                   side,
                                   trans,
                                   M,
                                   N,
                                   K,
                                   d_A,
                                   1,
                                   1,
                                   *descA,
                                   d_tau,
                                   d_C,
                                   1,
                                   1,
                                   *descC,
                                   CUDA_R_64F,
                                   d_work_ormqr,
                                   workspaceInBytesOnDevice_ormqr,
                                   h_work_ormqr,
                                   workspaceInBytesOnHost_ormqr,
                                   d_info_ormqr);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    if (verbose>1) printf("[i][dist_ormqr] after calling cusolverMpOrmqr\n");

    // sync after cusolverMpOrmqr
    cudaStat = cudaStreamSynchronize(*localStream);
    assert(cudaStat == cudaSuccess);

    // deallocate device workspace
    if (d_work_ormqr != NULL) {
        cudaStat = cudaFree(d_work_ormqr);
        assert(cudaStat == cudaSuccess);
        d_work_ormqr = NULL;
    }

    if (d_info_ormqr != NULL) {
        cudaStat = cudaFree(d_info_ormqr);
        assert(cudaStat == cudaSuccess);
        d_info_ormqr = NULL;
    }

    // deallocate host workspace
    if (h_work_ormqr) {
        free(h_work_ormqr);
        h_work_ormqr = NULL;
    }

    // barrier
    MPI_Barrier(MPI_COMM_WORLD);

    return 0;

}


//
// Wrapper for cusolverMpSytrd function that computes tridiagonal
// decomposition of a square symmetric matrix: A -> Q @ T @ Q'
//
int dist_mpsytrd(cudaStream_t *localStream,
        cusolverMpHandle_t *cusolverMpHandle,
        const int64_t N,   // number of columns in global matrix A
        cusolverMpMatrixDescriptor_t *descA, // matrix A descriptor
        double* d_A,       // matrix data
        double* d_d,       // output: diagonal
        double* d_e,       // output: subdiagonal, e(i) = A(i+1,i)
        double* d_tau) {   // output: scalar factors of the Householder reflectors

    // Error codes
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    calError_t       calStat      = CAL_OK;
    cudaError_t      cudaStat     = cudaSuccess;
    const int localDeviceId = getLocalRank();

    cudaStat = cudaSetDevice(localDeviceId);
    assert(cudaStat == cudaSuccess);
    cudaStat = cudaFree(0);
    assert(cudaStat == cudaSuccess);

    // size of workspace on device
    size_t workspaceInBytesOnDevice_sytrd = 0;

    // size of workspace on host
    size_t workspaceInBytesOnHost_sytrd = 0;

    // Query workspace size for mp routines
    cusolverStat = cusolverMpSytrd_bufferSize(*cusolverMpHandle,
                                              CUBLAS_FILL_MODE_LOWER,
                                              N,
                                              d_A,
                                              1,
                                              1,
                                              *descA,
                                              d_d,
                                              d_e,
                                              d_tau,
                                              CUDA_R_64F,
                                              &workspaceInBytesOnDevice_sytrd,
                                              &workspaceInBytesOnHost_sytrd);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);

    // error codes from cusolverMp (device)
    int* d_info_sytrd = NULL;

    // Allocate d_info
    cudaStat = cudaMalloc((void**)&d_info_sytrd, sizeof(int));
    assert(cudaStat == cudaSuccess);

    // Reset d_info
    cudaStat = cudaMemset(d_info_sytrd, 0, sizeof(int));
    assert(cudaStat == cudaSuccess);

    // Distributed device workspace
    void* d_work_sytrd = NULL;

    // allocate geqrf workspace
    cudaStat = cudaMalloc((void**)&d_work_sytrd,
                          workspaceInBytesOnDevice_sytrd);
    assert(cudaStat == cudaSuccess);

    // Distributed host workspace
    void* h_work_sytrd = NULL;

    h_work_sytrd = (void*)malloc(workspaceInBytesOnHost_sytrd);
    assert(h_work_sytrd != NULL);

    // call cusolverMpSytrd
    cusolverStat = cusolverMpSytrd(*cusolverMpHandle,
                                   CUBLAS_FILL_MODE_LOWER,
                                   N,
                                   d_A,
                                   1,
                                   1,
                                   *descA,
                                   d_d,
                                   d_e,
                                   d_tau,
                                   CUDA_R_64F,
                                   d_work_sytrd,
                                   workspaceInBytesOnDevice_sytrd,
                                   h_work_sytrd,
                                   workspaceInBytesOnHost_sytrd,
                                   d_info_sytrd);

    // sync after cusolverMpOrmqr
    cudaStat = cudaStreamSynchronize(*localStream);
    assert(cudaStat == cudaSuccess);

    // deallocate device workspace
    if (d_work_sytrd != NULL) {
        cudaStat = cudaFree(d_work_sytrd);
        assert(cudaStat == cudaSuccess);
        d_work_sytrd = NULL;
    }

    if (d_info_sytrd != NULL) {
        cudaStat = cudaFree(d_info_sytrd);
        assert(cudaStat == cudaSuccess);
        d_info_sytrd = NULL;
    }

    // deallocate host workspace
    if (h_work_sytrd) {
        free(h_work_sytrd);
        h_work_sytrd = NULL;
    }

    // barrier
    MPI_Barrier(MPI_COMM_WORLD);

    return 0;

}


//
// Wrapper for cusolverMpStedc function that computes eigenvalues
// of a symmetric tridiagonal matrix
//
int dist_mpstedc(cusolverMpHandle_t *cusolverMpHandle,
        const int64_t N,   // number of columns in global matrix A
        double* d_d,       // diagonal elements
        double* d_e,       // subdiagonal, e(i) = A(i+1,i)
        cusolverMpMatrixDescriptor_t *descQ, // descriptor of the rotation matrix Q
        double* d_Q,       // matrix Q
        int verbose) {     // verbose (1=yes, 0=no)

    // Error codes
    cusolverStatus_t cusolverStat = CUSOLVER_STATUS_SUCCESS;
    calError_t       calStat      = CAL_OK;
    cudaError_t      cudaStat     = cudaSuccess;
    const int localDeviceId = getLocalRank();
    if (verbose) {
        printf("[i][dist_mpstedc] localDeviceId = %d\n", localDeviceId);
    }

    cudaStat = cudaSetDevice(localDeviceId);
    assert(cudaStat == cudaSuccess);
    cudaStat = cudaFree(0);
    assert(cudaStat == cudaSuccess);

    // size of workspace on device
    size_t workspaceInBytesOnDevice_stedc = 0;

    // size of workspace on host
    size_t workspaceInBytesOnHost_stedc = 0;

    // error codes from cusolverMp (device)
    int* d_info_stedc = NULL;
    int h_info;
    cudaStat = cudaMalloc((void**)&d_info_stedc, sizeof(int));

    // A single-character flat to run in eigenvalue-only regime
    char compz = 'I'; // 'N' is not implemented (had 'Z' working here before - O.K.)

    // Query workspace size for mp routines
    cusolverStat = cusolverMpStedc_bufferSize(*cusolverMpHandle,
                                              &compz,
                                              N,
                                              d_d,
                                              d_e,
                                              d_Q,
                                              1,
                                              1,
                                              *descQ,
                                              CUDA_R_64F,
                                              &workspaceInBytesOnDevice_stedc,
                                              &workspaceInBytesOnHost_stedc,
                                              &h_info);
    assert(cusolverStat == CUSOLVER_STATUS_SUCCESS);
    //workspaceInBytesOnDevice_stedc *= 2;
    //workspaceInBytesOnHost_stedc   *= 2;
    if (verbose) {
        printf("[i][dist_mpstedc] workspaceInBytesOnDevice,OnHost = %d, %d\n",
            workspaceInBytesOnDevice_stedc,
            workspaceInBytesOnHost_stedc
        );
    }

    // Distributed device workspace
    void* d_work_stedc = NULL;

    //// Allocate d_info
    cudaStat = cudaMalloc((void**)&d_info_stedc, sizeof(int));
    assert(cudaStat == cudaSuccess);

    // allocate geqrf workspace
    cudaStat = cudaMalloc((void**)&d_work_stedc,
                          workspaceInBytesOnDevice_stedc);
    assert(cudaStat == cudaSuccess);

    // Distributed host workspace
    void* h_work_stedc = NULL;

    h_work_stedc = (void*)malloc(workspaceInBytesOnHost_stedc);
    assert(h_work_stedc != NULL);

    if (verbose) {
        printf("[i][dist_mpstedc] pointers before calling:\n");
        printf("[i][dist_mpstedc]  - d_info_stedc = %x\n", d_info_stedc);
        printf("[i][dist_mpstedc]  - d_work_stedc = %x\n", d_work_stedc);
        printf("[i][dist_mpstedc]  - h_work_stedc = %x\n", h_work_stedc);
        printf("[i][dist_mpstedc]  - d_d = %x\n", d_d);
        printf("[i][dist_mpstedc]  - d_e = %x\n", d_e);
        printf("[i][dist_mpstedc]  - d_Q = %x\n", d_Q);
    }
    // call cusolverMpSytrd
    cusolverStat = cusolverMpStedc(*cusolverMpHandle,
                                   &compz,
                                   N,
                                   d_d,
                                   d_e,
                                   d_Q,
                                   1,
                                   1,
                                   *descQ,
                                   CUDA_R_64F,
                                   d_work_stedc,
                                   workspaceInBytesOnDevice_stedc,
                                   h_work_stedc,
                                   workspaceInBytesOnHost_stedc,
                                   d_info_stedc);
                                   //&h_info);
    assert(cudaStat == cudaSuccess);
    if (verbose) {
        printf("[i][dist_stedc] h_info = %d\n", h_info);
        printf("[+][dist_stedc] AFTER call to stedc: \n");
        show_device_double_array(d_d, N);
    }

    // deallocate device workspace
    if (d_work_stedc != NULL) {
        cudaStat = cudaFree(d_work_stedc);
        assert(cudaStat == cudaSuccess);
        d_work_stedc = NULL;
    }

    if (d_info_stedc != NULL) {
        cudaStat = cudaFree(d_info_stedc);
        assert(cudaStat == cudaSuccess);
        d_info_stedc = NULL;
    }

    // deallocate host workspace
    if (h_work_stedc) {
        free(h_work_stedc);
        h_work_stedc = NULL;
    }

    // barrier
    MPI_Barrier(MPI_COMM_WORLD);

    return 0;

}

