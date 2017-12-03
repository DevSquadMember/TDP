#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "parser.h"
#include "cblas.h"
#include "perf.h"
#include "utils.h"

#define TRUE 1
#define FALSE 0

#define TAG 99

#define NB_DIMENSIONS 2

#define A_FILE "a.txt"
#define B_FILE "b.txt"

int size;

/* Matrices A, B et C complètes sur le process 0 */
double* matrix_a = NULL;
double* matrix_b = NULL;
double* matrix_c = NULL;
double* matrix_seq_c = NULL;

/* Matrices A, B et C locales (morceaux) */
double* matrix_local_a = NULL;
double* matrix_local_b = NULL;
double* matrix_local_c = NULL;

// Récupération des matrices
int init_matrices(int nb_blocs) {
    int sizeB;

    // Chargement de la matrice A
    ///printf("Matrice A\n");
    size = parse_size(A_FILE);
    if (size == -1) {
        return -1;
    }
    matrix_a = malloc(sizeof(double) * size * size);
    parse_matrix(matrix_a, size);

    // Chargement de la matrice B
    ///printf("Matrice B\n");
    sizeB = parse_size(B_FILE);
    if (sizeB == -1) {
        return -1;
    } else if (sizeB != size) {
        printf("Les deux matrices n'ont pas la même taille\n");
        return -1;
    }
    matrix_b = malloc(sizeof(double) * size * size);
    parse_matrix(matrix_b, size);

    matrix_c = malloc(sizeof(double) * size * size);
    matrix_seq_c = malloc(sizeof(double) * size * size);

    ///printf("Matrices de taille : %d\n", size);

    return size/nb_blocs;
}

void end() {
    free_matrix(matrix_a);
    free_matrix(matrix_b);
    free_matrix(matrix_c);

    free_matrix(matrix_seq_c);

    free_matrix(matrix_local_a);
    free_matrix(matrix_local_b);
    free_matrix(matrix_local_c);

    MPI_Finalize();
}

int main(int argc, char **argv) {
    struct group world_group, grid_group, col_group, row_group;
    int nb_blocs;
    int dimensions[NB_DIMENSIONS], periods[NB_DIMENSIONS], remain_dims[NB_DIMENSIONS], grid_coords[NB_DIMENSIONS];
    int reorder;

    perf_t start, stop, start_scatter, stop_scatter, start_calcul, stop_calcul, start_gather, stop_gather;

    MPI_Status status;
    MPI_Init(NULL, NULL);

    world_group.comm = MPI_COMM_WORLD;
    init_group(&world_group);

    nb_blocs = (int) sqrt(world_group.size);

    if (world_group.rank == 0) {
        printf("> Lancement du calcul en parallèle sur %d processeurs\n", world_group.size);
        printf("-- Début des calculs\n");
    }

    if (nb_blocs*nb_blocs != world_group.size) {
        printf("> Impossible de découper la matrice avec %d processeurs, %d n'est pas un carré\n", world_group.size, world_group.size);
        end();
        return EXIT_SUCCESS;
    }

    // La grille/matrice est de taille nb_procs*nb_procs
    dimensions[0] = nb_blocs;
    dimensions[1] = nb_blocs;
    periods[0] = FALSE;
    periods[1] = TRUE;
    reorder = FALSE;

    // Création de la topologie cartésienne
    if (MPI_Cart_create(MPI_COMM_WORLD, NB_DIMENSIONS, dimensions, periods, reorder, &(grid_group.comm))  != MPI_SUCCESS) {
        printf("Erreur pour grid_group\n");
    }
    init_group(&grid_group);

    // Groupe pour la colonne
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    if (MPI_Cart_sub(grid_group.comm, remain_dims, &(col_group.comm)) != MPI_SUCCESS) {
        printf("Erreur pour col_group\n");
    }
    init_group(&col_group);

    // Groupe pour la ligne
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    if (MPI_Cart_sub(grid_group.comm, remain_dims, &(row_group.comm)) != MPI_SUCCESS) {
        printf("Erreur pour row_group\n");
    }
    init_group(&row_group);

    // Récupération des coordonnées
    MPI_Cart_coords(grid_group.comm, grid_group.rank, NB_DIMENSIONS, grid_coords);

    int local_size = 0;
    int sendcounts[grid_group.size];
    int displs[grid_group.size];

    //printf("Proc %d (%d, %d) is %d in row_group and %d in col_groud\n", grid_group.rank, grid_coords[0], grid_coords[1], row_group.rank, col_group.rank);

    // Le processeur 0 charge les matrices et calcule la taille locale des blocs pour chaque processeur
    if (grid_group.rank == 0) {
        local_size = init_matrices(nb_blocs);
    }

    MPI_Bcast(&local_size, 1, MPI_INT, 0, grid_group.comm);

    if (local_size == -1) {
        end();
        return EXIT_SUCCESS;
    }

    for (int i = 0 ; i < grid_group.size ; i++) {
        sendcounts[i] = 1;
        //displs[i] = (int) (floor(i / nb_blocs) * nb_blocs * local_size + i%nb_blocs);
        displs[i] = i % nb_blocs * nb_blocs * local_size + i / nb_blocs;
        //printf("displs[%d] = %d\n", i, displs[i]);
    }

    double* matrix_buffer_a = malloc(local_size*local_size* sizeof(double));
    matrix_local_a = malloc(local_size*local_size* sizeof(double));
    matrix_local_b = malloc(local_size*local_size* sizeof(double));
    matrix_local_c = malloc(local_size*local_size* sizeof(double));

    // Initialisation de C
    /*for (int i = 0 ; i < local_size ; i++) {
        for (int j = 0 ; j < local_size ; j++) {
            matrix_local_c[j*local_size + i] = 0;
        }
    }*/

    // Définition du bloc pour chaque processeur
    MPI_Datatype bloc;
    MPI_Type_vector(local_size, local_size, local_size*nb_blocs, MPI_DOUBLE, &bloc); // nb_blocs
    MPI_Type_create_resized(bloc, 0, local_size * sizeof(double), &bloc);
    MPI_Type_commit(&bloc);

    perf(&start);

    perf(&start_scatter);

    // Découpage des matrices A et B et envoi d'un bloc sur chaque processus
    MPI_Scatterv(matrix_a, sendcounts, displs, bloc, matrix_local_a, local_size*local_size, MPI_DOUBLE, 0, grid_group.comm);
    MPI_Scatterv(matrix_b, sendcounts, displs, bloc, matrix_local_b, local_size*local_size, MPI_DOUBLE, 0, grid_group.comm);

    perf(&stop_scatter);

    /** Multiplication de Fox **/

    perf(&start_calcul);

    //printf("Proc %d is at %d/%d on grid - %d/%d\n", grid_group.rank, grid_coords[0], grid_coords[1], row_group.rank, col_group.rank);

    int dest = (col_group.rank - 1 + nb_blocs) % nb_blocs;
    int src = (col_group.rank + 1) % nb_blocs;

    //printf("Proc at %d/%d (%d) is sending B to %d and receiving from %d in col group\n", grid_coords[0], grid_coords[1], col_group.rank, dest, src);

    for (int i = 0 ; i < nb_blocs ; i++ ) {

        // Broadcast de A sur la ligne
        if (row_group.rank == (col_group.rank + i)%nb_blocs) {
            printf("Iteration %d - Le process %d va envoyer son A=%lf\n", i, grid_group.rank, matrix_local_a[0]);
	    // Le processeur qui va partager sa partie locale de A la charge dans le buffer
            memcpy(matrix_buffer_a, matrix_local_a, local_size * local_size * sizeof(double));
        }
        MPI_Bcast(matrix_buffer_a, local_size * local_size, MPI_DOUBLE, (col_group.rank + i )%nb_blocs, row_group.comm);

        // Au premier tour, le bloc B nécessaire vient d'être récupéré avec le scatterv
        if (i == 0) {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, local_size, local_size, 1., matrix_buffer_a, local_size, matrix_local_b, local_size, 0., matrix_local_c, local_size);
        } else {
            // Echange de B --> décalage vers le haut
            MPI_Sendrecv_replace(matrix_local_b, local_size*local_size, MPI_DOUBLE, dest, TAG, src, TAG, col_group.comm, &status);

            // Multiplication des blocs locaux A et B
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, local_size, local_size, 1., matrix_buffer_a, local_size, matrix_local_b, local_size, 1., matrix_local_c, local_size);
        }
	if (grid_group.rank == 0) {
            printf("Iteration %d - A is %lf and B is %lf for rank %d\n", i, matrix_buffer_a[0], matrix_local_b[0], grid_group.rank); 
        }
    }

    perf(&stop_calcul);

    /** Fin de la Multiplication de Fox **/

    perf(&start_gather);

    // Rassemblement des blocs de C sur le processeur 0
    MPI_Gatherv(matrix_local_c, local_size*local_size, MPI_DOUBLE, matrix_c, sendcounts, displs, bloc, 0, grid_group.comm);

    perf(&stop_gather);

    perf(&stop);

    if (grid_group.rank == 0) {
        printf("Matrice C\n");
        print_matrix(matrix_c, size);

        perf_diff(&start_scatter, &stop_scatter);
        printf("Temps pour le scatter : ");
        perf_printmicro(&stop_scatter);

        perf_diff(&start_calcul, &stop_calcul);
        printf("Temps pour le calcul : ");
        perf_printmicro(&stop_calcul);

        perf_diff(&start_gather, &stop_gather);
        printf("Temps pour le gather : ");
        perf_printmicro(&stop_gather);

        perf_diff(&start, &stop);
        double perf_par = perf_mflops(&stop, get_flops(size, grid_group.size));
        printf("Performance parallèle : %lf\n", perf_par);

        printf("Temps pour le parallèle : ");
        perf_printmicro(&stop);

        // Séquentiel

        perf(&start);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1., matrix_a, size, matrix_b, size, 0., matrix_seq_c, size);
        perf(&stop);
        printf("Vraie Matrice C\n");
        print_matrix(matrix_seq_c, size);

        perf_diff(&start, &stop);
        double perf_seq = perf_mflops(&stop, get_flops(size, 1));
        printf("Performance séquentielle : %lf\n", perf_seq);

        printf("Temps pour le séquentiel : ");
        perf_printmicro(&stop);

        double norm = calcul_norm2(matrix_c, matrix_seq_c, size);
        printf("Norme : %lf\n", norm);
    }

    MPI_Type_free(&bloc);
    MPI_Comm_free(&(row_group.comm));
    MPI_Comm_free(&(col_group.comm));
    MPI_Comm_free(&(grid_group.comm));
    end();

    return EXIT_SUCCESS;
}
