#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "utils-mpi.h"

#define RENDERING 1

#define TRUE 1
#define FALSE 0

#define TAG 0

int BS;

#define cell( _i_, _j_ ) local_board[ local_ldboard * (_j_) + (_i_) ]
#define ngb( _i_, _j_ )  local_nbngb[ local_ldnbngb * ((_j_) - 1) + ((_i_) - 1 ) ]

int rendering;

void end() {
    MPI_Finalize();
}

int main(int argc, char* argv[]) {
    int i, j, loop, num_alive, maxloop;
    int ldboard;
    int local_ldboard, local_ldnbngb;
	double t1, t2;
	double temps;

    int *board;
    int *nbngb;

    int *local_board;
    int *local_nbngb;

    struct group world_group, grid_group, col_group, row_group;
    int dimensions[2], periods[2], remain_dims[2], grid_coords[2];
    int reorder;

    MPI_Status status;
    MPI_Init(NULL, NULL);

    world_group.comm = MPI_COMM_WORLD;
    init_group(&world_group);

    int nb_blocs = (int) sqrt(world_group.size);

    if (nb_blocs*nb_blocs != world_group.size) {
        printf("> Impossible de découper le plateau avec %d processeurs, %d n'est pas un carré\n", world_group.size, world_group.size);
        end();
        return EXIT_SUCCESS;
    }

    // La grille des processeurs est de taille nb_procs*nb_procs
    dimensions[0] = nb_blocs;
    dimensions[1] = nb_blocs;
    periods[0] = FALSE;
    periods[1] = TRUE;
    reorder = FALSE;

	rendering = RENDERING;

	if (argc < 3) {
		printf("Usage: %s nb_iterations size [rendering=0/1]\n", argv[0]);
		return EXIT_SUCCESS;
	} else {
		maxloop = atoi(argv[1]);
		BS = atoi(argv[2]);
		if (argc > 3) {
			rendering = atoi(argv[3]);
		}
		//printf("Running MPI sync version, grid of size %d, %d iterations\n", BS, maxloop);
	}

    if (BS % nb_blocs != 0) {
        printf("> Impossible de découper le plateau de %d de côté avec %d blocs-processeurs\n", BS, nb_blocs);
        end();
        return EXIT_SUCCESS;
    }

    /** TOPOLOGIE CARTÉSIENNE **/

    // Création de la topologie cartésienne
    if (MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, reorder, &(grid_group.comm)) != MPI_SUCCESS) {
        printf("Erreur pour grid_group\n");
    }
    init_group(&grid_group);

    // Groupe pour la ligne
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    if (MPI_Cart_sub(grid_group.comm, remain_dims, &(row_group.comm)) != MPI_SUCCESS) {
        printf("Erreur pour row_group\n");
    }
    init_group(&row_group);

    // Groupe pour la colonne
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    if (MPI_Cart_sub(grid_group.comm, remain_dims, &(col_group.comm)) != MPI_SUCCESS) {
        printf("Erreur pour col_group\n");
    }
    init_group(&col_group);

    int sendcounts[grid_group.size];
    int displs[grid_group.size];
    int displs2[grid_group.size];

    // Récupération des coordonnées
    MPI_Cart_coords(grid_group.comm, grid_group.rank, 2, grid_coords);

    /** FIN TOPOLOGIE CARTÉSIENNE **/

    /* Leading dimension of the board array */
	ldboard = BS + 2;

    local_ldboard = BS/nb_blocs + 2;
    local_ldnbngb = BS/nb_blocs;

    // Le processeur 0 charge le plateau de jeu
    if (grid_group.rank == 0) {
        board = calloc(ldboard * ldboard, sizeof(int));
        nbngb = calloc(BS * BS, sizeof(int));

        num_alive = generate_initial_board(BS, &board[ldboard + 1], ldboard);
        if (rendering) {
            output_board(BS + 2, &board[0], ldboard, 0);
        }
        printf("%d cells are alive\n", num_alive);
        printf("Starting number of living cells = %d\n", num_alive);
    }

    local_board = calloc(local_ldboard * local_ldboard, sizeof(int));
    local_nbngb = calloc(local_ldnbngb * local_ldnbngb, sizeof(int));

    for (i = 0 ; i < grid_group.size ; i++) {
        sendcounts[i] = 1;
        displs[i] = i % nb_blocs * nb_blocs * (local_ldboard - 2) + i / nb_blocs;
        displs2[i] = i % nb_blocs * nb_blocs * (local_ldboard - 2) + i / nb_blocs;

        //printf("displs[%d] = %d\n", i, displs[i]);
    }

    // Définition du bloc pour chaque processeur
    MPI_Datatype bloc;
    MPI_Type_vector(local_ldboard - 1, local_ldboard, ldboard, MPI_INT, &bloc);
    //MPI_Type_create_resized(bloc, 0, (local_ldboard - 2) * sizeof(int), &bloc);
    MPI_Type_create_resized(bloc, 0, (local_ldboard - 2) * sizeof(int), &bloc);
    MPI_Type_commit(&bloc);

    MPI_Datatype row;
    MPI_Type_contiguous(local_ldboard, MPI_INT, &row);
    MPI_Type_commit(&row);

    MPI_Datatype col;
    MPI_Type_vector(local_ldboard - 2, 1, local_ldboard, MPI_INT, &col);
    MPI_Type_commit(&col);

    MPI_Datatype rec_bloc;
    MPI_Type_vector(local_ldboard - 2, local_ldboard - 2, BS, MPI_INT, &rec_bloc);
    MPI_Type_create_resized(rec_bloc, 0, (local_ldboard - 2) * sizeof(int), &rec_bloc);
    MPI_Type_commit(&rec_bloc);

    // Distribution des blocs sur les processeurs
    MPI_Scatterv(board, sendcounts, displs, bloc, &(local_board[0]), local_ldboard * (local_ldboard - 1), MPI_INT, 0, grid_group.comm);

    int top_rank = (col_group.rank - 1 + nb_blocs) % nb_blocs;
    int bottom_rank = (col_group.rank + 1) % nb_blocs;
    int left_rank = (row_group.rank - 1 + nb_blocs) % nb_blocs;
    int right_rank = (row_group.rank + 1) % nb_blocs;

    for (i = 0 ; i < grid_group.size ; i++) {
        if (grid_group.rank == i) {
            printf("PROC %d - COORDS : %d, %d\n", i, grid_coords[0], grid_coords[1]);
            output_board(local_ldboard, &(cell(0, 0)), local_ldboard, 0);
        }
        MPI_Barrier(grid_group.comm);
    }

    t1 = mytimer();

    for (loop = 1 ; loop <= maxloop ; loop++) {

        /// Echange des bordures

        // Envoi de la colonne de gauche
        MPI_Sendrecv(&(cell(1, 1)), 1, col, left_rank, TAG, &(cell(local_ldboard - 1, 1)), 1, col, right_rank, TAG, row_group.comm, &status);
        // Envoi de la colonne de droite
        MPI_Sendrecv(&(cell(local_ldboard - 2, 1)), 1, col, right_rank, TAG, &(cell(0, 1)), 1, col, left_rank, TAG, row_group.comm, &status);

        // Envoi de la ligne du dessus
        MPI_Sendrecv(&(cell(0, 1)), 1, row, top_rank, TAG, &(cell(0, local_ldboard - 1)), 1, row, bottom_rank, TAG, col_group.comm, &status);
        // Envoi de la ligne du dessous
        MPI_Sendrecv(&(cell(0, local_ldboard - 2)), 1, row, bottom_rank, TAG, &(cell(0, 0)), 1, row, top_rank, TAG, col_group.comm, &status);

        /// Calcul du nombre de voisins

        for (j = 1; j <= local_ldnbngb; j++) {
            for (i = 1; i <= local_ldnbngb; i++) {
                ngb( i, j ) =
                        cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
                        cell( i-1, j   ) +                  cell( i+1, j   ) +
                        cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
            }
        }

        for (i = 0 ; i < grid_group.size ; i++) {
            if (grid_group.rank == i) {
                printf("PROC %d - COORDS : %d, %d\n", i, grid_coords[0], grid_coords[1]);
                output_board(local_ldboard, &(cell(0, 0)), local_ldboard, 0);
            }
            MPI_Barrier(grid_group.comm);
        }

        /// Mise à jour des cellules

        num_alive = 0;
        for (j = 1; j <= local_ldnbngb; j++) {
            for (i = 1; i <= local_ldnbngb; i++) {
                if ( (ngb( i, j ) < 2) ||
                     (ngb( i, j ) > 3) ) {
                    cell(i, j) = 0;
                }
                else {
                    if ((ngb( i, j )) == 3)
                        cell(i, j) = 1;
                }
                if (cell(i, j) == 1) {
                    num_alive ++;
                }
            }
        }
    }

    for (i = 0 ; i < grid_group.size ; i++) {
        if (grid_group.rank == i) {
            printf("PROC %d - COORDS : %d, %d\n", i, grid_coords[0], grid_coords[1]);
            output_board(local_ldboard, &(cell(0, 0)), local_ldboard, 0);
        }
        MPI_Barrier(grid_group.comm);
    }

    t2 = mytimer();
	temps = t2 - t1;

    for (j = 1; j <= local_ldnbngb; j++) {
        for (i = 1; i <= local_ldnbngb; i++) {
            ngb(i, j) = cell(i, j);
        }
    }

    MPI_Gatherv(local_nbngb, local_ldnbngb * local_ldnbngb, MPI_INT, nbngb, sendcounts, displs2, rec_bloc, 0, grid_group.comm);

    double max_time;
    MPI_Reduce(&temps, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, grid_group.comm);

    int total_alive;
    MPI_Reduce(&num_alive, &total_alive, 1, MPI_INT, MPI_SUM, 0, grid_group.comm);

    if (grid_group.rank == 0) {
        if (rendering) {
            output_board(BS, &nbngb[0], ldboard - 2, loop - 1);
            printf("%d cells are alive\n", total_alive);
        }
        printf("Final number of living cells = %d\n", total_alive);
        printf("%.2lf\n",(double)temps * 1.e3);
        free(board);
        free(nbngb);
    }

    free(local_board);
    free(local_nbngb);

    MPI_Type_free(&bloc);
    MPI_Type_free(&row);
    MPI_Type_free(&col);
    MPI_Type_free(&rec_bloc);

    MPI_Comm_free(&(row_group.comm));
    MPI_Comm_free(&(col_group.comm));
    MPI_Comm_free(&(grid_group.comm));
    end();

	return EXIT_SUCCESS;
}

