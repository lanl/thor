/*--------------------------------------------------------------------------~*
 * Copyright (c) 2025 Triad National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~*/

/*
 * @file cuda_cublasmp_aux.c
 * @author Ismael Djibrilla Boureima, Oleg Korobkin
 * @date  January 2024
 * @brief Convenience C wrappers for cuBLASMp library
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <cublasmp.h>
#include <mpi.h>
/*
#ifdef USE_CAL_MPI
    #include <cal_mpi.h>
#endif
*/
#include <cal.h>

/////////////////////////////////////////////////////////////////////////////////////
//
//                CUBLAS TYPES and ENUMS
//
/////////////////////////////////////////////////////////////////////////////////////
// Reference: https://docs.nvidia.com/cuda/archive/8.0/cublas/index.html
//
// Type                Constant                    Description
//
// cublasOperation_t   CUBLAS_OP_N                 the non-transpose operation is selected
//                     CUBLAS_OP_T                 the transpose operation is selected
//                     CUBLAS_OP_C                 the conjugate transpose operation is selected
//
// cublasFillMode_t    CUBLAS_FILL_MODE_LOWER      lower part filled ('L')
//                     CUBLAS_FILL_MODE_UPPER      upper part filled ('U')
//
// cublasDiagType_t    CUBLAS_DIAG_NON_UNIT        the matrix diagonal has non-unit elements
//                     CUBLAS_DIAG_UNIT            the matrix diagonal has unit elements
//
// cublasSideMode_t    CUBLAS_SIDE_LEFT            the matrix is on the left side in the equation
//                     CUBLAS_SIDE_RIGHT           the matrix is on the right side in the equation
//
int show_cublasMpHandle(cublasMpHandle_t *handle){
  printf("[+][show_cublasMpHandle]: handle=>%x \n", handle);
  return 0;
}

int show_cublasMpGrid(cublasMpGrid_t *grid){
  printf("[+][show_cublasMpGrid]: grid=> %x \n", grid);
  return 0;
}

int show_cublasMpMatrixDescriptor(cublasMpMatrixDescriptor_t *descriptor){
  printf("[+][show_cublasMpMatrixDescriptor]: cublasMpMatrixDescriptor=%x \n", descriptor);
  return 0;
}
// ////////////////////////////////////////////////////////////////////////////////////
int set_cublasMpHandle(cublasMpHandle_t *handle, cudaStream_t *localStream, int VRBZ){
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    /* Initialize cublasMp library handle */
    if(VRBZ>0) printf("[+][set_cublasMpHandle()][INPUTs]: localStream=%d  \n", *localStream);
    //cublasStatus_t cublasMpCreate(cublasMpHandle_t *handle, cudaStream_t stream);
    cublasStat =     cublasMpCreate(handle, *localStream);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][set_cublasMpHandle]: cublasMpHandle INITIALIZED OK\n");
    return 0;
};

int destroy_cublasMpHandle(cublasMpHandle_t *handle, int VRBZ){
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    /* Destroy cublasMp handle */
    if(VRBZ>0) printf("[+][IN][destroy_cublasMpHandle()]\n");
    //cublasStatus_t cublasMpDestroy( cublasMpHandle_t handle);
    cublasStat =     cublasMpDestroy(*handle);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][destroy_cublasMpHandle()]: cublasMpHandle DESTROYED OK\n");
    return 0;
};

// cublasMp device grid
int set_cublasMpGrid(cublasMpGrid_t *grid, cal_comm_t *cal_comm,
                     int64_t numRowDevices, int64_t numColDevices, int VRBZ){
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    if(VRBZ>0) printf("[+][IN][set_cublasMpGrid()] numRowDevices=%d | numColDevices=%d| cal_comm=%x\n",
                                                   numRowDevices, numColDevices, cal_comm);
    // In CUDA 11.8 version:
    //cublasStatus_t cublasMpGridCreate(cublasMpHandle_t handle,
    //                                  int64_t nprow, int64_t npcol, int64_t myprow,
    //                                  cublasMpGridLayout_t layout,
    //                                  cal_comm_t comm, cublasMpGrid_t* grid);
    //
    // In CUDA 12.5 version:
    // cublasMpStatus_t cublasMpGridCreate(int64_t nprow, int64_t npcol,
	// 									cublasMpGridLayout_t layout,
	// 									ncclComm_t comm,
	// 									cublasMpGrid_t* grid);
    cublasStat= cublasMpGridCreate(numRowDevices, numColDevices,
                                   CUBLASMP_GRID_LAYOUT_COL_MAJOR, *cal_comm, grid);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][OUT][set_cublasMpGrid()]: grid=> %x CREATED OK\n", grid);
    return 0;
}

int destroy_cublasMpGrid(cublasMpGrid_t *grid,  int VRBZ){
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    if(VRBZ>0) printf("[+][IN][destroy_cublasMpGrid()]\n");
    //cublasStat = cublasMpGridDestroy(*handle, *grid); // CUDA 11.8 version
    cublasStat = cublasMpGridDestroy(*grid); // CUDA 12.6 version
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][destroy_cublasMpGrid()]: cublasMpGrid DESTROYED OK\n");
    return 0;
}


// cuBLASMp matrix descriptor
int set_cublasMpMatrixDesc(cublasMpHandle_t *handle, cublasMpMatrixDescriptor_t *descrA,
                           cublasMpGrid_t *gridA, int64_t M, int64_t N,
                           int64_t MB_A, int64_t NB_A,
                           int64_t numRowDevices,
                           int64_t dev_row_idx, int VRBZ) {
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    const uint64_t RSRCA = 0;
    const uint64_t CSRCA = 0;
    const int64_t llda   = cublasMpNumroc(M, MB_A, dev_row_idx, RSRCA, numRowDevices);
    if(VRBZ>0) {
        printf("[+][IN][set_cublasMpMatrixDesc()]: M=%d| N=%d| MB_A=%d| NB_A=%d| llda=%d\n"
               "                                   numRowDevices=%d | dev_row_idx=%d\n",
               M, N, MB_A, NB_A, llda, numRowDevices, dev_row_idx);
    }
    // hndl  Host  In  cuBLASMp handle
    // M     Host  In  Number of rows in the global matrix
    // N     Host  In  Number of columns in the global matrix
    // MB    Host  In  Blocking factor used to distribute the rows of the global matrix
    // NB    Host  In  Blocking factor used to distribute the columns of the global matrix
    // RSRC  Host  In  Row rank of the process who owns the first row block of the global matrix
    // CSRC  Host  In  Column rank of the process who owns the 1st column block of the global matrix
    // LLD   Host  In  Leading dimension of the local matrix
    // TYPE  Host  In  Data type of the matrix
    // GRID  Host  In  Grid object associated with the matrix descriptor
    // DESC  Host  Out Matrix descriptor object initialized by this function
    // // CUDA 11.8 syntax
    // cublasStat = cublasMpMatrixDescriptorCreate(*handle, M, N, MB_A, NB_A, RSRCA, CSRCA,
    //                                             llda, CUDA_R_64F, *gridA, descrA);
    cublasStat = cublasMpMatrixDescriptorCreate(M, N, MB_A, NB_A, RSRCA, CSRCA,
                                                llda, CUDA_R_64F, *gridA, descrA);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][OUT][set_cublasMpMatrixDesc()]: descriptor=> %x CREATED OK\n", descrA);
    return 0;
}

int destroy_cublasMpMatrixDesc(cublasMpHandle_t *handle, cublasMpMatrixDescriptor_t *descrA, int VRBZ){
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    if(VRBZ>0) printf("[+][IN][destroy_cublasMpMatrixDesc()]\n");
    //cublasStatus_t cublasMpMatrixDescriptorDestroy(cublasMpHandle_t handle,
    //                                               cublasMpMatrixDescriptor_t desc);
    // // CUDA 11.8 version:
    // cublasStat =     cublasMpMatrixDescriptorDestroy(*handle, *descrA);
    cublasStat =     cublasMpMatrixDescriptorDestroy(*descrA);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][destroy_cublasMpMatrixDesc()]: cublasMpMatrixDescriptor DESTROYED OK\n");
    return 0;
}

//////////////////////////////////////////////////////////////////////////////
//
//                            MATRIX SCATTER
//
/////////////////////////////////////////////////////////////////////////////
//
// Redistributes general rectangular matrix A
// according to the distribution properties of matrix B.
//
int redistribute_rectMatrix(cublasMpHandle_t *handle,
        cal_comm_t *cal_comm, cudaStream_t *localStream,
        int64_t M, int64_t N,
        cublasMpMatrixDescriptor_t *descrA, double* d_A,
        cublasMpMatrixDescriptor_t *descrB, double* d_B, int VRBZ) {

    if(VRBZ>0) printf("[+][IN][redstribute_rectMatrix]: M=%d| N=%d \n", M, N);
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    calError_t       calStat  = CAL_OK;
    cudaError_t      cudaStat = cudaSuccess;
    const uint64_t RSRCA = 0;
    const uint64_t CSRCA = 0;

    size_t workspaceSizeInBytesOnDevice = 0;
    size_t workspaceSizeInBytesOnHost   = 0;

    void* d_work = NULL;
    void* h_work = NULL;

    cublasStat = cublasMpGemr2D_bufferSize(*handle, M, N,
                     d_A, 1, 1, *descrA, d_B, 1, 1, *descrB,
                     &workspaceSizeInBytesOnDevice,
                     &workspaceSizeInBytesOnHost, *cal_comm);

    calStat = cal_stream_sync(*cal_comm, *localStream);
    if (VRBZ>0) printf("[%c][redisribute_rectMatrix]:"
                       " cublasMpGemr2D_bufferSize exit status: %d = %s\n",
                       (cublasStat == CUBLAS_STATUS_SUCCESS ? '+' : '-'),
                       cublasStat, cublasGetStatusString(cublasStat));
    assert(calStat == CAL_OK);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);

    cudaStat = cudaMalloc((void**)&d_work, workspaceSizeInBytesOnDevice);
    assert(cudaStat == cudaSuccess);
    h_work = (void*)malloc(workspaceSizeInBytesOnHost);
    if (VRBZ>0)
        printf("[+][redisribute_rectMatrix]: d_work and h_work alloc OK \n");

    cublasStat = cublasMpGemr2D(*handle, M, N,
                     d_A, 1, 1, *descrA, d_B, 1, 1, *descrB,
                     d_work, workspaceSizeInBytesOnDevice,
                     h_work, workspaceSizeInBytesOnHost, *cal_comm);

    calStat = cal_stream_sync(*cal_comm, *localStream);
    if (VRBZ>0) printf("[%c][redisribute_rectMatrix]:"
                       " cublasMpGemr2D exit status: %d = %s\n",
                       (cublasStat == CUBLAS_STATUS_SUCCESS ? '+' : '-'),
                       cublasStat, cublasGetStatusString(cublasStat));
    assert(calStat == CAL_OK);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return 0;
};

///////////////////////////////////////////////////////////////////////////////
//
//  MATRIX Addition/Subtraction
//
////////////////////////////////////////////////////////////////////////////
// ===============================================   C += A : C = A + C
// /////////////////////////////////////////////////////////////////////////
int cublasMpGeadd_wrapper(cublasMpHandle_t *handle,
          cal_comm_t *cal_comm,
          cudaStream_t *localStream,
          int64_t M, int64_t N, int trA, double alpha, double beta,
          cublasMpMatrixDescriptor_t *descrA, double* d_A,
          cublasMpMatrixDescriptor_t *descrC, double* d_C, int VRBZ) {
    //
    // This function performs the general matrix-matrix addition
    //
    //   C = alpha*OP(A) + beta*C
    //
    // where A and C are matrices with dimensions M x N.
    // trA = 0: no transpose
    // trA = 1: yes transpose (matrix A is stored in transposed state)
    //
    if(VRBZ>0) printf("[+][IN][cublasMpGeadd_wrapper]: M=%d| N=%d |\n", M,N);
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    calError_t       calStat  = CAL_OK;
    cudaError_t      cudaStat = cudaSuccess;
    const uint64_t RSRCA  = 0;
    const uint64_t CSRCA  = 0;
    if(VRBZ>0) printf("[+][cublasMpGeadd_wrapper]: d_A => %x | d_C => %x\n", d_A, d_C);
    size_t workspaceSizeInBytesOnDevice = 0;
    size_t workspaceSizeInBytesOnHost   = 0;
    void* d_work = NULL;
    void* h_work = NULL;
    cublasComputeType_t computeType = CUBLAS_COMPUTE_64F;

    // Compute workspace buffer sizes
    // cublasMpStatus_t cublasMpGeadd_bufferSize(
    //       cublasMpHandle_t handle,
    //       cublasOperation_t trans, int64_t m, int64_t n,
    //       const void* alpha, const void* a, int64_t ia, int64_t ja,
    //       cublasMpMatrixDescriptor_t descA,
    //       const void* beta,  void* c,       int64_t ic, int64_t jc,
    //       cublasMpMatrixDescriptor_t descC,
    //       size_t* workspaceSizeInBytesOnDevice,
    //       size_t* workspaceSizeInBytesOnHost);
    cublasStat = cublasMpGeadd_bufferSize(*handle,
                     CUBLAS_OP_N, M, N,
                     &alpha, d_A, 1, 1, *descrA,
                     &beta,  d_C, 1, 1, *descrC,
                     &workspaceSizeInBytesOnDevice,
                     &workspaceSizeInBytesOnHost);

    // sync wait for data to arrive
    calStat = cal_stream_sync(*cal_comm, *localStream);
    assert(calStat == CAL_OK);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][cublasMpGeadd_wrapper]: "
                          "cublasMpGeadd_bufferSize executed OK\n");

    // Allocate work buffers
    cudaStat = cudaMalloc((void**)&d_work, workspaceSizeInBytesOnDevice);
    assert(cudaStat == cudaSuccess);
    h_work = (void*)malloc(workspaceSizeInBytesOnHost);
    assert(h_work != NULL);
    if(VRBZ>0) printf("[+][cublasMpGeadd_wrapper]: "
                          "d_work and h_work allocated OK\n");
    // Calculate C = alpha*OP(A) + beta*C:/
    // cublasMpStatus_t cublasMpGeadd(
    //       cublasMpHandle_t handle,
    //       cublasOperation_t trans, int64_t m, int64_t n,
    //       const void* alpha, const void* a, int64_t ia, int64_t ja,
    //       cublasMpMatrixDescriptor_t descA,
    //       const void* beta,        void* c, int64_t ic, int64_t jc,
    //       cublasMpMatrixDescriptor_t descC,
    //       void* d_work, size_t workspaceSizeInBytesOnDevice,
    //       void* h_work, size_t workspaceSizeInBytesOnHost);
    cublasStat = cublasMpGeadd(*handle,
                     (trA == 0) ? CUBLAS_OP_N : CUBLAS_OP_T, M, N,
                     &alpha, d_A, 1, 1, *descrA,
                     &beta,  d_C, 1, 1, *descrC,
                     d_work, workspaceSizeInBytesOnDevice,
                     h_work, workspaceSizeInBytesOnHost);

    calStat = cal_stream_sync(*cal_comm, *localStream);
    assert(calStat == CAL_OK);
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    if(VRBZ>0) printf("[+][cublasMpGeadd_wrapper]: cublasMpGeadd executed OK\n");
    return 0;
};

/////////////////////////////////////////////////////////////////////////////
//
//           MATRIX Multiplication: C = A.OP @ B.OP
//
/////////////////////////////////////////////////////////////////////////////
int cublasMpGemm_wrapper(cublasMpHandle_t *handle,
        cal_comm_t *cal_comm, cudaStream_t *localStream,
        int64_t M, int64_t N, int64_t K, int trA, int trB,
        cublasMpMatrixDescriptor_t *descrA, double* d_A,
        cublasMpMatrixDescriptor_t *descrB, double* d_B,
        cublasMpMatrixDescriptor_t *descrC, double* d_C, int VRBZ){
    if(VRBZ>0) {
        printf("[+][IN][A.%cxB.%c]: M=%d |N=%d |K=%d |trA=%d |trB=%d\n"
               "               cublasmp_handle=%x| cal_comm=%x| stream=%x\n",
               (trA ? 'T' : 'N'), (trB ? 'T' : 'N'),
               M, N, K, trA, trB, handle, cal_comm, localStream);
    }
    // TODO: check MB_A, NB_A, MB_B, NB_B, MB_C, NB_C
    cublasStatus_t cublasStat = CUBLAS_STATUS_SUCCESS;
    calError_t       calStat  = CAL_OK;
    cudaError_t      cudaStat = cudaSuccess;
    const uint64_t RSRCA = 0;
    const uint64_t CSRCA = 0;
    if(VRBZ>0) printf("[+][A.%cxB.%c]: d_A=%x |d_B=%x |d_C=%x "
                      "|descrA=%x |descrB=%x |descrC=%x\n",
                      (trA ? 'T' : 'N'), (trB ? 'T' : 'N'),
                      d_A, d_B, d_C, descrA, descrB, descrC);

    size_t workspaceSizeInBytesOnDevice = 0;
    size_t workspaceSizeInBytesOnHost   = 0;
    double* d_work = NULL;
    double* h_work = NULL;
    double alpha   = 1.0;
    double beta    = 0.0;

    cublasComputeType_t computeType = CUBLAS_COMPUTE_64F;
    cublasStat = cublasMpGemm_bufferSize(*handle,
                       (trA ? CUBLAS_OP_T : CUBLAS_OP_N),
                       (trB ? CUBLAS_OP_T : CUBLAS_OP_N),
                       M,  N,  K,  &alpha,
                       d_A, 1, 1, *descrA,
                       d_B, 1, 1, *descrB, &beta,
                       d_C, 1, 1, *descrC,
                       computeType, //CUDA_R_64F,
                       &workspaceSizeInBytesOnDevice,
                       &workspaceSizeInBytesOnHost);
    if (VRBZ>0) printf("[%c][A.%cxB.%c]: cublasMpGemm_bufferSize exit status: %d = %s\n",
                       (cublasStat == CUBLAS_STATUS_SUCCESS ? '+' : '-'),
                       (trA ? 'T' : 'N'), (trB ? 'T' : 'N'),
                       cublasStat, cublasGetStatusString(cublasStat));
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    workspaceSizeInBytesOnDevice *= 2;

    // allocate   workspaces:
    cudaStat = cudaMallocAsync((void**)&d_work, workspaceSizeInBytesOnDevice, *localStream);
    assert(cudaStat == cudaSuccess);
    h_work = (double*)malloc(workspaceSizeInBytesOnHost);
    assert(h_work != NULL);
    calStat = cal_stream_sync(*cal_comm, *localStream);
    cal_comm_barrier(*cal_comm, *localStream);
    assert(calStat == CAL_OK);
    if(VRBZ>0) printf("[+][A.%cxB.%c]: cublasMpGemm_bufferSize executed OK:\n"
                      "           workspaceInBytesOnDevice  = %d\n"
                      "           workspaceInBytesOnHost    = %d\n",
                      (trA ? 'T' : 'N'), (trB ? 'T' : 'N'),
                      workspaceSizeInBytesOnDevice, workspaceSizeInBytesOnHost);

    // calculate C = 1.0*A@B.T + 0*C
    cublasStat = cublasMpGemm(*handle,
                       (trA ? CUBLAS_OP_T : CUBLAS_OP_N),
                       (trB ? CUBLAS_OP_T : CUBLAS_OP_N),
                        M,  N,  K,  &alpha,
                        d_A, 1, 1, *descrA,
                        d_B, 1, 1, *descrB, &beta,
                        d_C, 1, 1, *descrC,
                        computeType, //CUDA_R_64F,
                        d_work, workspaceSizeInBytesOnDevice,
                        h_work, workspaceSizeInBytesOnHost);
    calStat = cal_stream_sync(*cal_comm, *localStream);
    assert(calStat == CAL_OK);
    if (VRBZ>0) printf("[%c][A.%cxB.%c]: cublasMpGemm exit status: %d = %s\n",
                       (cublasStat == CUBLAS_STATUS_SUCCESS ? '+' : '-'),
                       (trA ? 'T' : 'N'), (trB ? 'T' : 'N'),
                       cublasStat, cublasGetStatusString(cublasStat));
    assert(cublasStat == CUBLAS_STATUS_SUCCESS);
    return 0;
}

