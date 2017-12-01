#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "parser.h"
#include "cblas.h"
#include <string.h>

#define TRUE 1
#define FALSE 0

#define TAG 99

#define NB_DIMENSIONS 2

#define A_FILE "a.txt"
#define B_FILE "b.txt"

struct group {
    MPI_Comm comm;
    int size;
    int rank;
};

void init_group(struct group* group) {
    MPI_Comm_rank(group->comm, &group->rank);
    MPI_Comm_size(group->comm, &group->size);
}

void print_matrix(double* matrix, int size) {
    for (int i = 0 ; i < size ; i++) {
        for (int j = 0 ; j < size ; j++) {
            printf("%lf ", matrix[j*size + i]);
        }
        printf("\n");
    }
}

int size;

/* Matrices A, B et C complètes sur le process 0 */
double* matrix_a = NULL;
double* matrix_b = NULL;
double* matrix_c = NULL;

/* Matrices A, B et C locales (morceaux) */
double* matrix_local_a = NULL;
double* matrix_local_b = NULL;
double* matrix_local_c = NULL;

// Récupération des matrices
int init_matrices(int nb_blocs) {
    int sizeB;

    // Chargement de la matrice A
    printf("Matrice A\n");
    size = parse_size(A_FILE);
    if (size == -1) {
        return -1;
    }
    matrix_a = malloc(sizeof(double) * size * size);
    parse_matrix(matrix_a, size);

    // Chargement de la matrice B
    printf("Matrice B\n");
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

    return size/nb_blocs;
}

void free_matrix(double* matrix) {
    if (matrix != NULL) {
        free(matrix);
    }
}

void end() {
    free_matrix(matrix_a);
    free_matrix(matrix_b);
    free_matrix(matrix_c);

    free_matrix(matrix_local_a);
    free_matrix(matrix_local_b);
    free_matrix(matrix_local_c);

    MPI_Finalize();
}

int main(int argc, char **argv) {
    struct group world_group, grid_group, col_group, row_group;
    int nb_blocs;
    int dimensions[NB_DIMENSIONS], periods[NB_DIMENSIONS], remain_dims[NB_DIMENSIONS], grid_coords[NB_DIMENSIONS];
    //int col_coords[1], row_coords[1];
    int reorder;

    MPI_Status status;
    MPI_Request request_recv;
    MPI_Request request_send;
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
    reorder = TRUE;

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
        displs[i] = (int) (floor(i / nb_blocs) * nb_blocs * local_size + i%nb_blocs);
        //displs[i] = (int) (floor(i / local_size) * nb_blocs + i);
        ///printf("displs %d is %d\n", i, displs[i]);
    }

    double* matrix_buffer_a = malloc(local_size*local_size* sizeof(double));
    matrix_local_a = malloc(local_size*local_size* sizeof(double));
    matrix_local_b = malloc(local_size*local_size* sizeof(double));
    matrix_local_c = malloc(local_size*local_size* sizeof(double));

    // Définition du bloc pour chaque processeur
    MPI_Datatype bloc;
    MPI_Type_vector(local_size, local_size, local_size*nb_blocs, MPI_DOUBLE, &bloc); // nb_blocs
    MPI_Type_create_resized(bloc, 0, local_size * sizeof(double), &bloc);
    MPI_Type_commit(&bloc);

    // Découpage des matrices A et B et envoi d'un bloc sur un processus
    MPI_Scatterv(matrix_a, sendcounts, displs, bloc, matrix_local_a, local_size*local_size, MPI_DOUBLE, 0, grid_group.comm);
    MPI_Scatterv(matrix_b, sendcounts, displs, bloc, matrix_local_b, local_size*local_size, MPI_DOUBLE, 0, grid_group.comm);

    /** Calcul **/
    for (int i = 0 ; i < local_size ; i++) {
        for (int j = 0 ; j < local_size ; j++) {
            printf("Process %d has for A(%d, %d) : %lf\n", grid_group.rank, i, j, matrix_local_a[j*local_size + i]);
            printf("Process %d has for B(%d, %d) : %lf\n", grid_group.rank, i, j, matrix_local_b[j*local_size + i]);
            matrix_local_c[j*local_size + i] = grid_group.rank;
        }
    }
    /** Fin des calculs **/

    /** multiplication de fox **/
    for (int i = 0 ; i < nb_blocs ; i++ ) {
        if(row_group.rank == (col_group.rank + i)%nb_blocs){
	        memcpy(matrix_buffer_a, matrix_local_a, local_size*local_size*sizeof(double));
        }
        // Broadcast de A sur la ligne
        MPI_Bcast(matrix_buffer_a, local_size * local_size, MPI_DOUBLE, (col_group.rank + i )%nb_blocs, row_group.comm);
        // récupération de B
        if (i != 0) {
            MPI_Isend(matrix_local_b, local_size * local_size, MPI_DOUBLE, ((col_group.rank - 1) + nb_blocs) % nb_blocs, TAG, col_group.comm, &request_send);
            MPI_Irecv(matrix_local_b, local_size * local_size, MPI_DOUBLE, (col_group.rank + 1) % nb_blocs, TAG, col_group.comm, &request_recv);

            // multiplication des blocs
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, local_size, local_size, local_size, 1., matrix_buffer_a, local_size, matrix_local_b, local_size, 1., matrix_local_c, local_size);
        }
    }


    /**if (col_rank == 1) {
        W = rank;
        V[0] = rank;
        V[1] = rank;
        V[2] = rank;
    }

    // Envoi de V aux voisins avec le rang de la deuxième cellule
    //MPI_Scatter(V, 1, MPI_INT, &W, 1, MPI_INT, 1, column);

    // Communication avec les voisins (dessus et dessous)
    if (col_rank == 1) {
        int rank_top, rank_bottom;
        MPI_Cart_shift(column, 0, 1, &rank_bottom, &rank_top);
        MPI_Send(&W, 1, MPI_INT, rank_bottom, 99, column);
        MPI_Send(&W, 1, MPI_INT, rank_top, 99, column);
    } else {
        MPI_Recv(&W, 1, MPI_INT, 1, 99, column, &status);
    }

    printf("Proc %d is %d/%d (%d;%d) in grid and %d/%d (%d) in col - W : %d\n", rank, grid_rank, grid_size-1, grid_coords[0], grid_coords[1], col_rank, col_size-1, col_coords[0], W);
**/

    MPI_Gatherv(matrix_local_c, local_size*local_size, MPI_DOUBLE, matrix_c, sendcounts, displs, bloc, 0, grid_group.comm);

    if (grid_group.rank == 0) {
        printf("Matrice C\n");
        print_matrix(matrix_c, size);
    }

    MPI_Type_free(&bloc);
    MPI_Comm_free(&(row_group.comm));
    MPI_Comm_free(&(col_group.comm));
    MPI_Comm_free(&(grid_group.comm));
    end();

    return EXIT_SUCCESS;
}