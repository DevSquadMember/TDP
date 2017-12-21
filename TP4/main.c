#include <stdlib.h>
#include <stdio.h>
#include "lib_utils.h"
#include "lib_matrix.h"
#include "mpi.h"
#include "perf.h"
#include "utils.h"

/*
 * dgtrf_nopiv
 *   for j = 1 to n ; += b
 *     dgetf2 - OK
 *     dtrsm(L)
 *     dgemm
 */

int rank, size;

void solve_seq(struct matrix* A, struct vector* X, struct vector* B) {
    if (A->nb_rows < MAX_SIZE_TO_PRINT) {
        printf("\nMatrix A before\n");
        matrix_show(A);

        printf("\nVector B before\n");
        vector_show(B);
    }

    solve_sequential(A, X, B);

    if (A->nb_rows < MAX_SIZE_TO_PRINT) {
        printf("\nVector X after\n");
        vector_show(X);
    }
}

void solve_par(struct matrix* A, struct vector* X, struct vector* B) {
    solve_parallel(A, X, B);

    if (rank == 0 && A->nb_rows < MAX_SIZE_TO_PRINT) {
        printf("\nSolution X is \n");
        vector_show(X);
    }

    matrix_free(A);
    vector_free(X);
    vector_free(B);
}

int main(int argc, char** argv) {
    int matrix_size = MATRIX_SIZE;

    if (argc > 1) {
        matrix_size = atoi(argv[1]);
    }

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    struct matrix A;
    struct vector X, B;
    matrix_init(&A, matrix_size, matrix_size);
    vector_init(&X, matrix_size);
    vector_init(&B, matrix_size);

    if (rank == 0) {
        srand((unsigned int) time(NULL));

        matrix_generate(&A);
        vector_generate(&B);

        if (matrix_size < MAX_SIZE_TO_PRINT) {
            printf("\nVECTOR B\n");
            vector_show(&B);

            printf("\nMatrix A\n");
            matrix_show(&A);
        }

        printf("\t\n---------- SEQUENTIAL ----------\n");

        struct matrix a;
        struct vector x, b;

        matrix_copy(&A, &a);
        vector_copy(&B, &b);
        solve_seq(&a, &x, &b);

        matrix_free(&a);
        vector_free(&x);
        vector_free(&b);

        printf("\t\n---------- PARALLEL (%d) ----------\n", size);
    }

    solve_par(&A, &X, &B);

    matrix_free(&A);
    vector_free(&X);
    vector_free(&B);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
