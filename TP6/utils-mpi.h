#ifndef TDP_UTILS_MPI_H
#define TDP_UTILS_MPI_H

#include <mpi.h>

struct group {
    MPI_Comm comm;
    int size;
    int rank;
};

double mytimer(void);

void output_board(int N, int *board, int ldboard, int loop);

int generate_initial_board(int N, int *board, int ldboard);

void init_group(struct group* group);

#endif //TDP_UTILS_MPI_H
