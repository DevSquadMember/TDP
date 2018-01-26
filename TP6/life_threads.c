#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <semaphore.h>
#include "barrier.h"
#include "utils.h"

#define RENDERING 1

int BS;

#define cell( _i_, _j_ ) board[ ldboard * (_j_) + (_i_) ]
#define ngb( _i_, _j_ )  nbngb[ ldnbngb * ((_j_) - 1) + ((_i_) - 1 ) ]

int rendering;

int nb_threads;

int* nums_alive; /* Tableau du nombre de cellules vivantes pour chaque thread */
int ldboard, ldnbngb;
int *board; /* Tableau du Plateau du jeu */
int *nbngb; /* Tableau des compteurs de voisins array */

int maxloop;

struct barrier barrier_loop, barrier_ngb, barrier_sync, barrier_render;
sem_t* left_sem;
sem_t* right_sem;

int compute_num_alive(int nb_threads) {
    int num_alive = 0;
    for (int i = 0 ; i < nb_threads ; i++) {
        num_alive += nums_alive[i];
    }
    return num_alive;
}

void *thread_run(void* arg) {
    int rank = *((int *)arg);

    int nb_cols = BS/nb_threads;
    int start_col = 1 + nb_cols * rank;
    int last_col = start_col + nb_cols - 1;

    for (int loop = 1; loop <= maxloop; loop++) {
        printf("LOOP %d for rank %d\n", loop, rank);

        /** RECOPIE DES BORDURES POUR CHAQUE THREAD **/

        // Chaque thread recopie les bordures pour ses colonnes
        for (int i = start_col ; i <= last_col ; i++) {
            // recopie du bas en haut
            cell(i, 0) = cell(i, BS);
            // recopie du haut en bas
            cell(i, BS + 1) = cell(i, 1);
        }

        if (rank == 0) {
            cell(0, 0) = cell(BS, BS);
            cell(0, BS + 1) = cell(BS, 1);
            for (int i = 1 ; i <= BS ; i++) {
                cell(0, i) = cell(BS, i);
            }
        } else if (rank == nb_threads - 1) {
            for (int i = 1 ; i <= BS ; i++) {
                cell(BS + 1, i) = cell(1, i);
            }
            cell(BS + 1, 0) = cell(1, BS);
            cell(BS + 1, BS + 1) = cell(1, 1);
        }

        /** FIN DE LA RECOPIE DES BORDURES POUR CHAQUE THREAD **/

        printf("BARRIER SYNC - %d\n", rank);
        // Barrière de synchronisation
        barrier_stop(&barrier_sync);
        printf("OK for %d\n", rank);

        /** CALCUL DES CELLULES EN FONCTION DE SES VOISINS **/

        // Calcul de la première colonne
        for (int j = 1; j <= BS; j++) {
            ngb( start_col, j ) =
                    cell( start_col-1, j-1 ) + cell( start_col, j-1 ) + cell( start_col+1, j-1 ) +
                    cell( start_col-1, j   ) +                          cell( start_col+1, j   ) +
                    cell( start_col-1, j+1 ) + cell( start_col, j+1 ) + cell( start_col+1, j+1 );
        }
        sem_post(&(left_sem[rank]));

        // Calcul de la dernière colonne
        for (int j = 1; j <= BS; j++) {
                ngb( last_col, j ) =
                        cell( last_col-1, j-1 ) + cell( last_col, j-1 ) + cell( last_col+1, j-1 ) +
                        cell( last_col-1, j   ) +                         cell( last_col+1, j   ) +
                        cell( last_col-1, j+1 ) + cell( last_col, j+1 ) + cell( last_col+1, j+1 );
        }
        sem_post(&(right_sem[rank]));

        // Barrière de synchronisation pour le calcul des colonnes du bord
        printf("BARRIER NGB - %d\n", rank);
        barrier_stop(&barrier_ngb);
        printf("OK for %d\n", rank);

        for (int j = 1; j <= BS; j++) {
            for (int i = start_col + 1 ; i < last_col ; i++) {
                ngb( i, j ) =
                        cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
                        cell( i-1, j   ) +                  cell( i+1, j   ) +
                        cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
            }
        }

        /** FIN DU CALCUL DES CELLULES EN FONCTION DE SES VOISINS **/

        nums_alive[rank] = 0;

        for (int j = 1 ; j <= BS ; j++) {
            for (int i = start_col + 1 ; i < last_col ; i++) {
                if ( (ngb( i, j ) < 2) ||
                     (ngb( i, j ) > 3) ) {
                    cell(i, j) = 0;
                }
                else {
                    if ((ngb( i, j )) == 3)
                        cell(i, j) = 1;
                }
                if (cell(i, j) == 1) {
                    nums_alive[rank]++;
                }
            }
        }

        // première colonne à mettre à jour
        sem_wait(&(right_sem[rank - 1]));
        for (int j = 1 ; j <= BS ; j++) {
            if ( (ngb( start_col, j ) < 2) ||
                 (ngb( start_col, j ) > 3) ) {
                cell(start_col, j) = 0;
            }
            else {
                if ((ngb( start_col, j )) == 3)
                    cell(start_col, j) = 1;
            }
            if (cell(start_col, j) == 1) {
                nums_alive[rank]++;
            }
        }
        sem_post(&(right_sem[rank - 1]));

        // dernière colonne à mettre à jour
        sem_wait(&(left_sem[rank + 1]));
        for (int j = 1 ; j <= BS ; j++) {
            if ( (ngb( last_col, j ) < 2) ||
                 (ngb( last_col, j ) > 3) ) {
                cell(last_col, j) = 0;
            }
            else {
                if ((ngb( last_col, j )) == 3)
                    cell(last_col, j) = 1;
            }
            if (cell(last_col, j) == 1) {
                nums_alive[rank]++;
            }
        }
        sem_post(&(left_sem[rank + 1]));

        if (rendering) {
            barrier_stop(&barrier_render);
            if (rank == 0) {
                // Avec juste les "vraies" cellules: on commence à l'élément (1,1)
                output_board(BS+2, &(cell(0, 0)), ldboard, loop);

                int num_alive = compute_num_alive(nb_threads);
                printf("%d cells are alive\n", num_alive);
            }
        }
        sem_wait(&(right_sem[rank]));
        sem_wait(&(left_sem[rank]));

        printf("BARRIER STOP - %d\n", rank);
        barrier_stop(&barrier_loop);
    }

    printf("EXIT %d\n", rank);
    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
	int i, num_alive;
	double t1, t2;
	double temps;

	rendering = RENDERING;

	if (argc < 4) {
		printf("Usage: %s nb_iterations size num_threads [rendering=0/1]\n", argv[0]);
		return EXIT_SUCCESS;
	} else {
		maxloop = atoi(argv[1]);
		BS = atoi(argv[2]);
        nb_threads = atoi(argv[3]);
		if (argc > 4) {
			rendering = atoi(argv[4]);
		}
		//printf("Running sequential version, grid of size %d, %d iterations\n", BS, maxloop);
	}

    /** Création des threads **/
    pthread_t threads[nb_threads];
    int ranks[nb_threads];
    nums_alive = malloc(sizeof(int) * nb_threads);
    left_sem = malloc(sizeof(sem_t) * nb_threads);
    right_sem = malloc(sizeof(sem_t) * nb_threads);

    // Initialisation des barrières de synchronisation
    barrier_init(&barrier_sync, nb_threads);
    barrier_init(&barrier_ngb, nb_threads);
    barrier_init(&barrier_loop, nb_threads);
    barrier_init(&barrier_render, nb_threads);

	/* Leading dimension of the board array */
	ldboard = BS + 2;
	/* Leading dimension of the neigbour counters array */
	ldnbngb = BS;

	board = malloc( ldboard * ldboard * sizeof(int) );
	nbngb = malloc( ldnbngb * ldnbngb * sizeof(int) );

	num_alive = generate_initial_board( BS, &(cell(1, 1)), ldboard );

    // Avec juste les "vraies" cellules: on commence à l'élément (1,1)
    if (rendering) {
        output_board(BS + 2, &(cell(0, 0)), ldboard, 0);
    }
    printf("%d cells are alive\n", num_alive);

	printf("Starting number of living cells = %d\n", num_alive);
    t1 = mytimer();

    for (i = 0 ; i < nb_threads ; i++) {
        sem_init(&(left_sem[i]), 0, 0);
        sem_init(&(right_sem[i]), 0, 0);
    }

    /** TRAVAIL DES THREADS **/
    for (i = 0 ; i < nb_threads ; i++) {
        ranks[i] = i;
        pthread_create(&(threads[i]), NULL, thread_run, &ranks[i]);
    }

    // Attente des threads
    for (i = 0 ; i < nb_threads ; i++) {
        pthread_join(threads[i], NULL);
    }

	t2 = mytimer();
	temps = t2 - t1;
    num_alive = compute_num_alive(nb_threads);
	printf("Final number of living cells = %d\n", num_alive);
	printf("%.2lf\n",(double)temps * 1.e3);

    barrier_free(&barrier_render);
    barrier_free(&barrier_sync);
    barrier_free(&barrier_ngb);
    barrier_free(&barrier_loop);
	free(board);
	free(nbngb);
    free(nums_alive);
    free(left_sem);
    free(right_sem);

	return EXIT_SUCCESS;
}

