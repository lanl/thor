/*--------------------------------------------------------------------------~*
 * Copyright (c) 2025 Triad National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~*/

/*
 * @file cuda_cal_aux.c
 * @author Ismael Djibrilla Boureima, Oleg Korobkin
 * @date  October 2023
 * @brief Convenience C wrappers for CAL (CUDA Communication Abstraction Library)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>  // NEEDED FOR BOOL type
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <cal.h>

calError_t allgather(void* src_buf, void* recv_buf, size_t size, void* data, void** request);
calError_t request_test(void* request);
calError_t request_free(void* request);
calError_t cal_comm_create_mpi(MPI_Comm mpi_comm, int rank, int nranks, int local_device, cal_comm_t* comm);
calError_t cal_comm_create_mpi_fortran(MPI_Fint mpi_comm, int rank, int nranks, int local_device, cal_comm_t* comm);
int mpi_is_initialized();
int CUDA_CHECK(cudaError_t cudaStat);
int set_cal_comm(cal_comm_t *cal_comm, int *rank, int *commSize, int *localRank, int VRBZ);
cal_comm_t get_cal_comm(int rank, int commSize, int localRank, int VRBZ);


int set_localStream(cudaStream_t *localStream, int VRBZ){
    cudaError_t  cudaStat = cudaSuccess;
    /* Create local stream */
      //printf("[+][set_localStream()] VRBZ=%d \n", VRBZ);
      if(VRBZ>0) printf("[+][IN set_localStream()]\n");
      cudaStat = cudaStreamCreate(localStream);
      assert(cudaStat == cudaSuccess);
      if(VRBZ>0) printf("[+][set_localStream()]: cudaStat = %x | localStream=%d \n", cudaStat, *localStream);
    return 0;
}

int show_localStream(cudaStream_t *localStream){
  printf("[+][show_localStream]: localStream=> %x \n", localStream);
  return 0;
}

calError_t allgather(void* src_buf, void* recv_buf, size_t size, void* data, void** request)
{
    MPI_Request req;
    int err = MPI_Iallgather(src_buf, size, MPI_BYTE, recv_buf, size, MPI_BYTE, (MPI_Comm)data, &req);
    if (err != MPI_SUCCESS)
    {
        return CAL_ERROR;
    }
    *request = (void*)req;
    return CAL_OK;
}

calError_t request_test(void* request)
{
    MPI_Request req = (MPI_Request)request;
    int         completed;
    int         err = MPI_Test(&req, &completed, MPI_STATUS_IGNORE);
    if (err != MPI_SUCCESS)
    {
        return CAL_ERROR;
    }
    return completed ? CAL_OK : CAL_ERROR_INPROGRESS;
}

calError_t request_free(void* request)
{
    return CAL_OK;
}

calError_t cal_comm_create_mpi(MPI_Comm mpi_comm, int rank, int nranks, int local_device, cal_comm_t* comm)
{
    cal_comm_create_params_t params;
    params.allgather = allgather;
    params.req_test = request_test;
    params.req_free = request_free;
    params.data = (void*)mpi_comm;
    params.rank = rank;
    params.nranks = nranks;
    params.local_device = local_device;
    return cal_comm_create(params, comm);
}

calError_t cal_comm_create_mpi_fortran(MPI_Fint mpi_comm, int rank, int nranks, int local_device, cal_comm_t* comm)
{
  return cal_comm_create_mpi(MPI_Comm_f2c(mpi_comm), rank, nranks, local_device, comm);
}

int mpi_is_initialized(){
    int mpi_stat;
    MPI_Initialized(&mpi_stat);
    return mpi_stat;
}

int CUDA_CHECK(cudaError_t cudaStat){
    assert(cudaStat == cudaSuccess);
    return 0;
}

int get_hw_topo(int *rank, int *commSize, MPI_Comm *localComm, int *localRank, int *localCommSize,
                MPI_Comm *ctxComm, int *ctxRank, int *ctxCommSize, int *nDevLocal, int ctx_color, int *color, int VRBZ){

    int stat;
    if(! mpi_is_initialized()){
        MPI_Init(NULL, NULL);
        if(VRBZ>0) printf("[+][get_hw_topo]: MPI Initialized OK\n");
    }
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, commSize);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, localComm);
    MPI_Comm_rank(*localComm, localRank);
    MPI_Comm_size(*localComm, localCommSize);

    /* Check local device inventory  & defines cuda ctx*/
    CUDA_CHECK( cudaGetDeviceCount( nDevLocal) );
    if( *localRank < *nDevLocal){
        CUDA_CHECK(cudaSetDevice(*localRank));
        CUDA_CHECK(cudaFree(0));
        *color    =  ctx_color;
    }
    else{
        *color    =  ctx_color + 1;
    }
    if(VRBZ>0) printf("[L%d/G%d][get_hw_topo]: nDevLocal=%d | ctx_color = %d | color=%d \n",
                      *localRank, *rank, *nDevLocal, ctx_color, *color);
    MPI_Barrier(*localComm);
    stat = MPI_Comm_split(*localComm, *color, *localRank, ctxComm);
    if(VRBZ>0 && stat==0) printf("[L%d/G%d][get_hw_topo]: CUDA CTX DEFINED OK\n",
                      *localRank, *rank);
    return stat;
}

int init_cal_comm(cal_comm_t *cal_comm, int *rank, int *commSize, int *localRank, int VRBZ){
    // Create CAL communicator
    calError_t calStat = CAL_OK;
    calStat = cal_comm_create_mpi(MPI_COMM_WORLD, *rank, *commSize, *localRank, cal_comm);
    if (VRBZ>0) {
        printf("[+][init_cal_comm] cal_comm_create_mpi called OK \n");
        printf("[L%d/G%d][init_cal_comm] calStat= %x | cal_comm= %x \n",
                *localRank, *rank, calStat, cal_comm);
    }
    assert(calStat == CAL_OK);
    if(VRBZ>0) printf("[L%d/G%d][init_cal_comm]: CAL comm initialized OK\n",
                       *localRank, *rank);
    return 0;
}

int show_cal_comm(cal_comm_t *cal_comm){
  printf("[+][show_cal_comm]: cal_comm=> %x \n", cal_comm);
  return 0;
}


int destroy_cal_comm(cal_comm_t *cal_comm,  int VRBZ){
        calError_t calStat = CAL_OK;
        calStat            =  cal_comm_destroy( *cal_comm);
        assert(calStat == CAL_OK);
        if(VRBZ>0) printf("[+][destroy_cal_comm()]:  Cal COMM DESTROYED OK| calStat  = %x \n", calStat);
        if(VRBZ>0) {
        if(calStat==CAL_OK)
            printf("[+][destroy_cal_comm()]:  Cal COMM DESTROYED OK| calStat = %x \n", calStat);
        else
            printf("[!][destroy_cal_comm()][ERR] UNABLE to DESTROYED Cal COMM: calStat = %x \n", calStat);
    }
    return 0;
}

