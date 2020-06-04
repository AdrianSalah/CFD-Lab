#include "parallel.hpp"
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
        int num_proc){

    //*il = 0, *jb = 0;
    int i_delta = (jmax + 2) / iproc, j_delta = (imax + 2) / jproc;
    *omg_i = 0, *omg_j = 0;

    if (num_proc > 1) {

        //il, ir, jb, jt, rank_r, rank_b, ...are all send buffers
        if (*myrank == 0) {
            for (int j = 0; j < jproc; j++) {
                for (int i = 0; i < iproc; i++) {
                    int proc = i + j * iproc + 1 ; //start sending to proc = 1

                    *il = i * i_delta;
                    *ir = i * i_delta + i_delta;

                    *jb = j * j_delta;
                    *jt = j * j_delta + j_delta;

                    *omg_i = i;
                    *omg_j = j;


                    if (*il == 0){*rank_l = MPI_PROC_NULL;}
                    if (*jb == 0){*rank_b = MPI_PROC_NULL;}
                    if (*ir == imax + 1){*rank_r = MPI_PROC_NULL;}
                    if (*jt == jmax + 1){*rank_t = MPI_PROC_NULL;}
                    else{
                        //...
                    }
                    if(proc < num_proc) {
                        MPI_Send(il, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                        MPI_Send(ir, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                        MPI_Send(jb, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                        MPI_Send(jt, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                        MPI_Send(omg_i, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                        MPI_Send(omg_j, 1, MPI_INT, proc, proc, MPI_COMM_WORLD);
                    }
                }
            }
        }

        int x_cells_tot, y_cells_tot;
        int x_layer_full_width;
        int y_layer_full_width;
        int x_layer_res, y_layer_res;
        int x_uniform_add;
        int x_uniform_res;
        int x_layer_width_final;
        int x_layer_width_final_last_col;

        x_layer_full_width = x_cells_tot/ iproc;
        x_layer_res = x_cells_tot % iproc;
        x_uniform_add = x_layer_res / (iproc - 1);

        x_layer_width_final = x_layer_full_width + x_uniform_add;
        x_layer_width_final_last_col = x_layer_width_final + x_layer_res % (iproc - 1);



        for (int p = 1; p < num_proc; p++){
            if (*myrank == p) {
                MPI_Recv(il, 1, MPI_INT, 0, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(ir, 1, MPI_INT, 0, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(jb, 1, MPI_INT, 0, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(jt, 1, MPI_INT, 0, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(omg_i, 1, MPI_INT, 0, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(omg_j, 1, MPI_INT, 0, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

    int num_ghostcells_j = (jt +1) - (jb -1);
    int num_ghostcells_i = (ir +1) - (il -1);
    //fill send buffer for right ghost cell
    for(int k; k < num_ghostcells_j; k++){
        bufSend[k] = P[ir+1][jb-1 + k];
    }
    //send ghost cells to right, receive ghost cells from left
    MPI_Sendrecv(bufSend, num_ghostcells_j, MPI_DOUBLE, rank_r, 123, bufRecv, num_ghostcells_j,
            MPI_DOUBLE, chunk, 123, MPI_COMM_WORLD, status );

}