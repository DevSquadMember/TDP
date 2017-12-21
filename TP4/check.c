#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "lib_utils.h"
#include "lib_matrix.h"
#include "mpi.h"
#include "perf.h"
#include "utils.h"

int rank, size;
struct timeval perf_begin, perf_end;

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

    perf(&perf_begin);
    dtrsm('D', 'L', 'N', 'U', 3, 6, 1, &m, &t);
    perf(&perf_end);

    perf_diff(&perf_begin, &perf_end);
    printf("DTRSM séquentiel sur matrice de taille 3 6 : ");
    perf_printmicro(&perf_end);

    printf("\nMatrix M after\n");
    matrix_show(&m);
    printf("\nMatrix T after\n");
    matrix_show(&t);
}

void check_lu_example() {
    struct matrix m;
    int m_size = 3;

    matrix_init(&m, m_size, m_size);
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

    perf(&perf_begin);
    // DGETF2
    dgetf2(&m);
    perf(&perf_end);
    perf_diff(&perf_begin, &perf_end);
    printf("DGETF2 séquentiel sur matrice de taille %d %d : ", m_size, m_size);
    perf_printmicro(&perf_end);

    printf("\nMatrix after\n");
    matrix_show(&m);

    matrix_free(&m);
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

void check_solve_seq(int n) {
    struct matrix m;
    struct vector x, b;

    matrix_init(&m, n, n);
    vector_init(&x, n);
    vector_init(&b, n);

    matrix_load(&m);
    vector_load(&b);

    if (n < MAX_SIZE_TO_PRINT) {
        printf("\nMatrix A before\n");
        matrix_show(&m);

        printf("\nVector B before\n");
        vector_show(&b);
    }

    solve_sequential_dgetf2(&m, &x, &b);
    //solve_sequential_dgetrf(&m, &x, &b);

    // Vérification du résultat
    check_correctness(&x);

    if (n < MAX_SIZE_TO_PRINT) {
        printf("\nVector X after\n");
        vector_show(&x);
    }
}

void check_solve_par(int n) {
    struct matrix a;
    struct vector x, b;
    matrix_init(&a, n, n);
    vector_init(&x, n);
    vector_init(&b, n);

    if (rank == 0) {
        matrix_load(&a);

        /*1.000000 1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
        1.000000 -1.000000 2.000000 0.000000 0.000000 0.000000 0.000000 0.000000
        1.000000 1.000000 -2.000000 3.000000 0.000000 0.000000 0.000000 0.000000
        1.000000 1.000000 1.000000 -3.000000 4.000000 0.000000 0.000000 0.000000
        1.000000 1.000000 1.000000 1.000000 -4.000000 5.000000 0.000000 0.000000
        1.000000 1.000000 1.000000 1.000000 1.000000 -5.000000 6.000000 0.000000
        1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 -6.000000 7.000000
        1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 -7.000000*/
        /*for (int i = 0 ; i < m.nb_cols ; i++) {
            for (int j = 0 ; j < m.nb_rows ; j++) {
                double value;
                if (i < j) {
                    value = 1;
                } else if (i > j+1) {
                    value = 0.;
                } else if (i == j+1) {
                    value = i;
                } else if (i == 0 && j == 0) {
                    value = 1.;
                } else {
                    value = -1*i;
                }
                matrix_set(&m, j, i, value);
            }
        }*/

        vector_load(&b);
        if (n < MAX_SIZE_TO_PRINT) {
            printf("\nVECTOR B\n");
            vector_show(&b);

            printf("\nMatrix A\n");
            matrix_show(&a);
        }
    }

    solve_parallel(&a, &x, &b);

    if (rank == 0) {
        // Vérification du résultat
        check_correctness(&x);

        if (n < MAX_SIZE_TO_PRINT) {
            printf("\nSolution X is \n");
            vector_show(&x);
        }
    }

    matrix_free(&a);
    vector_free(&x);
    vector_free(&b);
}

int main(int argc, char** argv) {
    int matrix_size = MATRIX_SIZE;

    if (argc > 1) {
        matrix_size = atoi(argv[1]);
    }

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        srand((unsigned int) time(NULL));
        printf("\t\n---------- SEQUENTIAL ----------\n");
        ///check_lu_example();
        ///check_trsm_example();
        //check_extract_matrix();
        check_solve_seq(matrix_size);
        /*timer_begin();
        check_lu_example();
        double t = timer_end();
        printf("Time : %lf microseconds\n", t);*/
        printf("\t\n---------- PARALLEL (%d) ----------\n", size);
    }

    check_solve_par(matrix_size);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
