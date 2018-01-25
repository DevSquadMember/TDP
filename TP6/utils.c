#include <stdio.h>
#include <sys/time.h>
#include "utils.h"

#define cell( _i_, _j_ ) board[ ldboard * (_j_) + (_i_) ]

double mytimer(void) {
    struct timeval tp;
    gettimeofday( &tp, NULL );
    return tp.tv_sec + 1e-6 * tp.tv_usec;
}

void output_board(int N, int *board, int ldboard, int loop) {
    int i,j;
    printf("loop %d\n", loop);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if ( cell( j, i ) == 1)
                printf("X");
            else
                printf(".");
        }
        printf("\n");
    }
}

/**
 * This function generates the initial board with one row and one
 * column of living cells in the middle of the board
 */
int generate_initial_board(int N, int *board, int ldboard) {
    int i, j, num_alive = 0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == N/2 || j == N/2) {
                cell(i, j) = 1;
                num_alive ++;
            }
            else {
                cell(i, j) = 0;
            }
        }
    }

    return num_alive;
}