#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "cblas.h"

#define MAX 100000
#define MIN (-100000)

void init_group(struct group* group) {
    MPI_Comm_rank(group->comm, &group->rank);
    MPI_Comm_size(group->comm, &group->size);
}

long long get_flops(int size, int nb_threads) {
    return (long long int) (1.0 * (size * size * size * 2) / nb_threads);
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

double* generate_matrix(int size) {
    double* matrix = malloc(size * size * sizeof(double));

    for (int i = 0 ; i < size ; i++) {
        for (int j = 0 ; j < size ; j++) {
            matrix[i*size + j] = ((rand()/(double)RAND_MAX) * (2*MAX) + (rand()/(double)RAND_MAX)) + MIN;
        }
    }

    return matrix;
}