#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "lib_matrix.h"
#include "timer.h"
#include "mpi.h"

/*
 * dgtrf_nopiv
 *   for j = 1 to n ; += b
 *     dgetf2 - OK
 *     dtrsm(L)
 *     dgemm
 */

void matrix_load(struct matrix* m) {
    double value;
    for (int i = 0 ; i < m->nb_rows ; i++) {
        for (int j = 0 ; j < m->nb_cols ; j++) {
            if (j == 0) {
                value = 1;
            } else if (j == i+1) {
                value = i+1;
            } else {
                value = 0;
            }
            matrix_set(m, i, j, value);
        }
    }
}

void vector_load(struct vector* v) {
    for (int i = 0 ; i < v->nb_values ; i++) {
        vector_set(v, i, (i+1) % (v->nb_values) + 1);
    }
}

void check_trsm_example() {
    struct matrix m;
    struct matrix t;

    matrix_init(&m,3,6);
    m.values[0] = 1;
    m.values[1] = 1;
    m.values[2] = 2;

    m.values[3] = 3;
    m.values[4] = 1;
    m.values[5] = 5;

    m.values[6] = 6;
    m.values[7] = 7;
    m.values[8] = 1;

    m.values[9]  = 1;
    m.values[10] = 2;
    m.values[11] = 3;

    m.values[12] = 4;
    m.values[13] = 5;
    m.values[14] = 6;

    m.values[15] = 7;
    m.values[16] = 8;
    m.values[17] = 9;

    matrix_sub(&t, 3, 3, 0, 3, &m);

    printf("\nCHECK TRSM\n");
    printf("\nMatrix M before\n");
    matrix_show(&m);
    printf("\nMatrix T before\n");
    matrix_show(&t);

    dtrsm('D', 'L', 'N', 'U', 3, 6, 1, &m, 3, &t, 3);

    printf("\nMatrix M after\n");
    matrix_show(&m);
    printf("\nMatrix T after\n");
    matrix_show(&t);
}

void check_lu_example() {
    struct matrix m;

    matrix_init(&m, 3, 3);
    m.values[0] = 2;
    m.values[1] = 1;
    m.values[2] = 1;

    m.values[3] = 4;
    m.values[4] = 3;
    m.values[5] = 5;

    m.values[6] = 4;
    m.values[7] = 1;
    m.values[8] = 6;

    printf("\nCHECK LU\n");
    printf("\nMatrix before\n");

    matrix_show(&m);

    // DGETF2
    dgetf2(&m);

    printf("\nMatrix after\n");
    matrix_show(&m);

    matrix_free(&m);
}

void check_solving_example(int n) {
    struct matrix m;
    struct vector x, b;

    matrix_init(&m, n, n);
    vector_init(&x, n);
    vector_init(&b, n);

    matrix_load(&m);
    vector_load(&b);

    printf("\nSolving problem\n");
    printf("\nMatrix A before\n");
    matrix_show(&m);

    printf("\nVector X before\n");
    vector_show(&x);
    printf("\nVector B before\n");
    vector_show(&b);

    // DGETF2
    dgetf2(&m);
    // Descente-remontée
    matrix_solve(&m, &x, &b);

    printf("\nVector X after\n");
    vector_show(&x);
}

void check_extract_matrix() {
    struct matrix m;
    struct matrix sub;

    matrix_init(&m, 16, 16);

    for (int i = 0 ; i < m.nb_rows ; i++) {
        for (int j = 0 ; j < m.nb_cols ; j++) {
            matrix_set(&m, i, j, 1);
        }
    }

    printf("BEFORE M\n");
    matrix_show(&m);
    matrix_sub(&sub, 4, 4, 6, 8, &m);

    for (int i = 0 ; i < sub.nb_rows ; i++) {
        for (int j = 0 ; j < sub.nb_cols ; j++) {
            matrix_set(&m, i, j, 2);
        }
    }
    printf("AFTER SUB\n");
    matrix_show(&sub);

    printf("AFTER M\n");
    matrix_show(&m);

    matrix_free(&m);
    matrix_free(&sub);
}

int main(int argc, char** argv) {
    /*struct matrix m;
    struct vector v;*/

    /*matrix_init(&m, 4, 4);

    printf("\nMatrice\n");
    matrix_load(&m);

    matrix_show(&m);

    printf("\nAfter FACTO LU\n");
    facto_lu(&m);

    matrix_show(&m);

    vector_init(&v, 4);

    printf("\nVector\n");
    vector_load(&v);

    vector_show(&v);

    matrix_free(&m);
    vector_free(&v);*/

    check_lu_example();
    check_trsm_example();
    check_solving_example(5);
    /*timer_begin();
    check_lu_example();
    double t = timer_end();
    printf("Time : %lf microseconds\n", t);*/

    /**int rank, size;

    MPI_Status status;
    MPI_Request request_recv;
    MPI_Request request_send;
    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int BLOC_LENGTH = 5;
    int M_SIZE = 20;

    MPI_Datatype bloc;
    MPI_Type_vector(M_SIZE, BLOC_LENGTH, BLOC_LENGTH, MPI_DOUBLE, &bloc); // nb_blocs
    //MPI_Type_create_resized(bloc, 0, local_size * sizeof(double), &bloc);
    MPI_Type_commit(&bloc);**/

    // Découpage des matrices A et B et envoi d'un bloc sur un processus
    //MPI_Scatterv(matrix_a, sendcounts, displs, bloc, matrix_local_a, local_size*local_size, MPI_DOUBLE, 0, grid_group.comm);


    //check_extract_matrix();

    return EXIT_SUCCESS;
}
