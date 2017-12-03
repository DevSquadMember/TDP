#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "cblas.h"

void init_group(struct group* group) {
    MPI_Comm_rank(group->comm, &group->rank);
    MPI_Comm_size(group->comm, &group->size);
}

long get_flops(int size, int nb_threads) {
    return (long) (1.0 * (size * size * size * 2) / nb_threads);
}

// MATRICES

void print_matrix(double* matrix, int size) {
    for (int i = 0 ; i < size ; i++) {
        for (int j = 0 ; j < size ; j++) {
            printf("%lf ", matrix[j*size + i]);
        }
        printf("\n");
    }
}

void free_matrix(double* matrix) {
    if (matrix != NULL) {
        free(matrix);
    }
}

double calcul_norm2(const double* Cpara, const double* Cseq, int size) {
    double* matrix = malloc(size * size * sizeof(double));
    int i;
    for(i = 0 ; i < size*size ; i++)
        matrix[i] = Cpara[i] - Cseq[i];
    double norme2 = cblas_dnrm2(size*size, matrix, 1);
    free(matrix);
    return norme2;
}
