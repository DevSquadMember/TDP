#ifndef TDP_UTILS_H
#define TDP_UTILS_H

#include <mpi.h>

struct group {
    MPI_Comm comm;
    int size;
    int rank;
};

// Initialiser un groupe avec la taille et le rang
void init_group(struct group* group);

// Nombre de flops en fonction de la taille des matrices
long long get_flops(int size, int nb_threads);

// Afficher une matrice
void print_matrix(double* matrix, int size);

void free_matrix(double* matrix);

double calcul_norm2(const double* Cpara, const double* Cseq, int size);

#endif //TDP_UTILS_H
