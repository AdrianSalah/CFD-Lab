#ifndef __PARALLEL_HPP__
#define __PARALLEL_HPP__
#include <mpi.h>
#include "datastructures.hpp"
#include "enums.hpp"

enum CommunicationType {
    TO_RIGHT    = 0,
    TO_LEFT     = 1,
    TO_TOP      = 2,
    TO_BOTTOM   = 3
};


void init_parallel(
    int iproc,
    int jproc,
    int imax,
    int jmax,
    int* myrank,
    int* il,
    int* ir,
    int* jb,
    int* jt,
    int* rank_l,
    int* rank_r,
    int* rankn_b,
    int* rank_t,
    int* omg_i,
    int* omg_j,
    int num_proc);


void pressure_comm(
        matrix<double>& P,
    int il,
    int ir,
    int jb,
    int jt,
    int rank_l,
    int rank_r,
    int rank_b,
    int rank_t,
    double* bufSend,
    double* bufRecv,
    MPI_Status* status,
    int chunk);

void u_velocity_comm(
        matrix<double> U,
        int il,
        int ir,
        int jb,
        int jt,
        int rank_l,
        int rank_r,
        int rank_b,
        int rank_t,
        double* bufSend,
        double* bufRecv,
        MPI_Status* status,
        int chunk);

void v_velocity_comm(
        matrix<double> V,
        int il,
        int ir,
        int jb,
        int jt,
        int rank_l,
        int rank_r,
        int rank_b,
        int rank_t,
        double* bufSend,
        double* bufRecv,
        MPI_Status* status,
        int chunk);

void f_comm(
        matrix <double> F,
        int il,
        int ir,
        int jb,
        int jt,
        int rank_l,
        int rank_r,
        int rank_b,
        int rank_t,
        double* bufSend,
        double* bufRecv,
        MPI_Status* status,
        int chunk);

void g_comm(
        matrix <double> G,
        int il,
        int ir,
        int jb,
        int jt,
        int rank_l,
        int rank_r,
        int rank_b,
        int rank_t,
        double* bufSend,
        double* bufRecv,
        MPI_Status* status,
        int chunk);


#endif // __PARALLEL_HPP__