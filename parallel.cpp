#include "parallel.hpp"
#include <math.h>
#include <mpi.h>

void init_parallel(
        int iproc,
        int jproc,
        int imax,
        int jmax,
        int *myrank,
        int *il,
        int *ir,
        int *jb,
        int *jt,
        int *rank_l,
        int *rank_r,
        int *rank_b,
        int *rank_t,
        int *omg_i,
        int *omg_j,
        int num_proc) {


    // compute indices for subdomains
    float i_delta_f = static_cast<float>(imax + 2) / static_cast<float>(iproc);
    float j_delta_f = static_cast<float>(jmax + 2) / static_cast<float>(jproc);
    int i_delta = ceil(i_delta_f);
    int j_delta = ceil(j_delta_f);
    int i_delta_last_col = (imax + 2) - (iproc - 1) * i_delta;
    int j_delta_last_row = (jmax + 2) - (jproc - 1) * j_delta;

    int ir_arr[iproc], il_arr[iproc], jb_arr[jproc], jt_arr[jproc];

    /* omg_i = {1, 2,..., iproc}, omg_j = {1, 2,..., jproc} */
    int omg_i_arr[iproc], omg_j_arr[jproc];

    // store indices for columns
    for (int i = 0; i < iproc - 1; i++) {
        il_arr[i] = i * i_delta;
        ir_arr[i] = i * i_delta + i_delta;

        omg_i_arr[i] = i + 1;
    }
    // store indices for last column
    il_arr[iproc - 1] = (iproc - 1) * i_delta;
    ir_arr[iproc - 1] = (iproc - 1) * i_delta + i_delta_last_col;

    omg_i_arr[iproc - 1] = iproc;

    // store indices for rows
    for (int j = 0; j < jproc - 1; j++) {
        jb_arr[j] = j * j_delta;
        jt_arr[j] = j * j_delta + j_delta;

        omg_j_arr[j] = j + 1;
    }

    //store indices for last row
    jb_arr[jproc - 1] = (jproc - 1) * j_delta;
    jt_arr[jproc - 1] = (jproc - 1) * j_delta + j_delta_last_row;

    omg_j_arr[jproc - 1] = jproc;

    // send the indices from rank 0 (master process) to the assigned processes
    if (*myrank == 0) {
        for (int j = 0; j < jproc; j++) {
            for (int i = 0; i < iproc; i++) {
                int proc = (i + j * iproc); //% num_proc; //start sending to process 0
                //int tag = (i + j * iproc);

                *omg_i = omg_i_arr[i];
                *omg_j = omg_j_arr[j];

                *il = il_arr[i];
                *ir = ir_arr[i];
                *jb = jb_arr[j];
                *jt = jt_arr[j];


                // check weather neighbour process is outside of boundary or not and specify ranks of neighbour processes
                if (*il == 0){*rank_l = MPI_PROC_NULL;}
                else{ *rank_l = proc - 1;}
                if (*jb == 0){*rank_b = MPI_PROC_NULL;}
                else{*rank_b = proc - iproc;}
                if (*ir == imax+2){*rank_r = MPI_PROC_NULL;}
                else{*rank_r = proc + 1;}
                if (*jt == jmax+2){*rank_t = MPI_PROC_NULL;}
                else{*rank_t = proc + iproc;}



                MPI_Send(omg_i, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(omg_j, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);

                MPI_Send(il, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(ir, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(jb, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(jt, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);

                MPI_Send(rank_l, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(rank_r, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(rank_b, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                MPI_Send(rank_t, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);

            }
        }
    }
    // receive indices from rank 0
    for (int j = 0; j < jproc; j++) {
        for (int i = 0; i < iproc; i++) {
            //int tag = (i + j * iproc);
            int proc = (i + j * iproc); //% num_proc;
            if (*myrank == proc) {
                MPI_Recv(omg_i, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(omg_j, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                MPI_Recv(il, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(ir, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(jb, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(jt, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                MPI_Recv(rank_l, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rank_r, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rank_b, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(rank_t, 1, MPI_INT, 0, proc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
}









void pressure_comm(
        double **P,
        int il,
        int ir,
        int jb,
        int jt,
        int rank_l,
        int rank_r,
        int rank_b,
        int rank_t,
        double *bufSend,
        double *bufRecv,
        MPI_Status *status,
        int chunk){

    int num_ghostcells_j = jt - jb + 3; //inner cells + 2 ghost layers in y-direction
    int num_ghostcells_i = ir - il + 3; //inner cells + 2 ghost layers in y-direction

    bufSend[num_ghostcells_j] = {0};
    //fill send buffer with left ghost cells
    for(int j = 0; j < num_ghostcells_j; j++){
        bufSend[j] = P[0][j];
    }
    /* ---- send to the left - receive from the right ---- */
    MPI_Sendrecv(bufSend, num_ghostcells_j, MPI_DOUBLE, rank_l, 1, bufRecv, num_ghostcells_j,
                 MPI_DOUBLE, chunk, 2, MPI_COMM_WORLD, status );


    bufSend[num_ghostcells_j] = {0};
    //fill send buffer with right ghost cells
    for(int j = 0; j < num_ghostcells_j; j++){
        bufSend[j] = P[ir+1][j];
    }
    /* ---- send to right - receive from left ---- */
    MPI_Sendrecv(bufSend, num_ghostcells_j, MPI_DOUBLE, rank_r, 2, bufRecv, num_ghostcells_j,
            MPI_DOUBLE, chunk, 1, MPI_COMM_WORLD, status );


    bufSend[num_ghostcells_i] = {0};
    //fill send buffer with bottom ghost cells
    for(int i = 0; i < num_ghostcells_i; i++){
        bufSend[i] = P[i][0];
    }
    /* ---- send to bottom - receive from top ---- */
    MPI_Sendrecv(bufSend, num_ghostcells_i, MPI_DOUBLE, rank_b, 3, bufRecv, num_ghostcells_i,
                 MPI_DOUBLE, chunk, 4, MPI_COMM_WORLD, status );


    bufSend[num_ghostcells_i] = {0};
    //fill send buffer with bottom ghost cells
    for(int i = 0; i < num_ghostcells_i; i++){
        bufSend[i] = P[i][jt + 1];
    }
    /* ---- send to top - receive from bottom ---- */
    MPI_Sendrecv(bufSend, num_ghostcells_i, MPI_DOUBLE, rank_t, 4, bufRecv, num_ghostcells_i,
                 MPI_DOUBLE, chunk, 3, MPI_COMM_WORLD, status );





}
