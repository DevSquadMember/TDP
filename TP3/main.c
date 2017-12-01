#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "parser.h"
#include "cblas.h"
#include <string.h>

#define TRUE 1
#define FALSE 0

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

int main(int argc, char **argv) {
    struct group world_group, grid_group, col_group, row_group;
    int nb_blocs;
    int dimensions[NB_DIMENSIONS], periods[NB_DIMENSIONS], remain_dims[NB_DIMENSIONS], grid_coords[NB_DIMENSIONS];
    //int col_coords[1], row_coords[1];
    int reorder;

    MPI_Status status;
    MPI_Init(NULL, NULL);

    world_group.comm = MPI_COMM_WORLD;
    init_group(&world_group);

    nb_blocs = (int) sqrt(world_group.size);

    if (world_group.rank == 0) {
        printf("Lancement du calcul en parallèle sur %d processeurs\n", world_group.size);
        printf("-- Début des calculs\n");
    }

    if (nb_blocs*nb_blocs != world_group.size) {
        printf("Impossible de découper la matrice avec %d processeurs\n", world_group.size);
        return EXIT_FAILURE;
    }

    // La grille/matrice est de taille nb_procs*nb_procs
    dimensions[0] = nb_blocs;
    dimensions[1] = nb_blocs;
    periods[0] = FALSE;
    periods[1] = TRUE;
    reorder = TRUE;

    // Création de la topologie cartésienne
    MPI_Cart_create(MPI_COMM_WORLD, NB_DIMENSIONS, dimensions, periods, reorder, &grid_group.comm);
    init_group(&grid_group);

    // Groupe pour la colonne
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(grid_group.comm, remain_dims, &col_group.comm);
    init_group(&col_group);

    // Groupe pour la ligne
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(grid_group.comm, remain_dims, &row_group.comm);
    init_group(&row_group);

    // Récupération des coordonnées
    MPI_Cart_coords(grid_group.comm, grid_group.rank, NB_DIMENSIONS, grid_coords);

    int size, local_size;
    int sendcounts[world_group.size];
    int displs[world_group.size];
    double* matrix_a = malloc(sizeof(double));
    double* matrix_b = malloc(sizeof(double));
    double* matrix_c = malloc(sizeof(double));

    if (world_group.rank == 0) {
        // Récupération des matrices
        printf("Matrice A\n");
        parse(A_FILE, &matrix_a, &size);
        printf("Matrice B\n");
        parse(B_FILE, &matrix_b, &size);
        local_size = size/nb_blocs;
    }

    MPI_Bcast(&local_size, 1, MPI_INT, 0, grid_group.comm);
    for (int i = 0 ; i < world_group.size ; i++) {
        sendcounts[i] = 1;
        displs[i] = (int) (floor(i / nb_blocs) * nb_blocs * local_size + i%nb_blocs);
        //displs[i] = (int) (floor(i / local_size) * nb_blocs + i);
        ///printf("displs %d is %d\n", i, displs[i]);
    }

    double* matrix_buffer_a = malloc(local_size*local_size* sizeof(double));
    double* matrix_local_a = malloc(local_size*local_size* sizeof(double));
    double* matrix_local_b = malloc(local_size*local_size* sizeof(double));
    double* matrix_local_c = malloc(local_size*local_size* sizeof(double));

    // Définition du bloc pour chaque processeur
    MPI_Datatype bloc;
    MPI_Type_vector(local_size, local_size, local_size*nb_blocs, MPI_DOUBLE, &bloc); // nb_blocs
    MPI_Type_create_resized(bloc, 0, local_size * sizeof(double), &bloc);
    MPI_Type_commit(&bloc);

    // Découpage des matrices A et B et envoi d'un bloc sur un processus
    MPI_Scatterv(matrix_a, sendcounts, displs, bloc, matrix_local_a, local_size*local_size, MPI_DOUBLE, 0, world_group.comm);
    MPI_Scatterv(matrix_b, sendcounts, displs, bloc, matrix_local_b, local_size*local_size, MPI_DOUBLE, 0, world_group.comm);

    /** Calcul **/
    for (int i = 0 ; i < local_size ; i++) {
        for (int j = 0 ; j < local_size ; j++) {
            printf("Process %d has for A(%d, %d) : %lf\n", world_group.rank, i, j, matrix_local_a[j*local_size + i]);
            printf("Process %d has for B(%d, %d) : %lf\n", world_group.rank, i, j, matrix_local_b[j*local_size + i]);
            matrix_local_c[j*local_size + i] = world_group.rank;
        }
    }
    /** Fin des calculs **/

    /** multiplication de fox **/
    for (int i =0 ; i < nb_blocs : i++ ){
      if(row_group.rank == (col_group.rank + i)%nb_blocs){
	memcpy(matrix_buffer_a,matrix_local_a,local_size*local_size*sizeof(double));
      }
      // Broadcast de A sur la ligne
      MPI_Bcast(&matrix_buffer_a,local_size*local_size,MPI_DOUBLE,(col_group.rank + i )%nb_blocs,row_group.comm);
      // récupération de B
      if(i != 0){
	MPI_Isend(matrix_local_b,local_size*local_size,MPI_DOUBLE,((col_group.rank - 1) + nb_blocs)%nb_blocs ,99,col_group.comm);
	MPI_Irecv(matrix_local_b,local_size*local_size,MPI_DOUBLE,(col_group.rank + 1)%nb_blocs,99,col_group.comm);

      // multiplication des blocs
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,local_size,local_size,local_size,1.,matrix_buffer_a,local_size,matrix_local_b,local_size,1.,matric_local_c,local_size);
      
      


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

    MPI_Gatherv(matrix_local_c, local_size*local_size, MPI_DOUBLE, matrix_c, sendcounts, displs, bloc, 0, MPI_COMM_WORLD);

    if (world_group.rank == 0) {
        print_matrix(matrix_c, size);
    }

    free(matrix_a);
    free(matrix_b);
    free(matrix_c);
    free(matrix_local_a);
    free(matrix_local_b);
    free(matrix_local_c);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
