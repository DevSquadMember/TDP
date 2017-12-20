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

int rank, size;

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

    printf("\nVector B before\n");
    vector_show(&b);

    // DGETF2
    dgetf2(&m);
    /*printf("\nMatrix A after\n");
    matrix_show(&m);*/
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

void solve_parallel(struct matrix* A, struct vector* X, struct vector* B) {
    struct matrix local_a;
    struct vector local_x, local_b;
    int m_size = A->nb_rows;

    int local_size = m_size/size;
    matrix_init(&local_a, local_size, m_size);
    vector_init(&local_x, local_size);
    vector_init(&local_b, local_size);

    // Type col : colonne de la matrice à envoyer une par une sur chaque processeur pour la matrice A
    MPI_Datatype col;
    MPI_Type_vector(local_size, m_size, m_size*size, MPI_DOUBLE, &col);
    MPI_Type_create_resized(col, 0, m_size * sizeof(double), &col);
    MPI_Type_commit(&col);

    // Type cell : case du vecteur à envoyer une par une sur chaque processeur pour les vecteurs
    MPI_Datatype cell;
    MPI_Type_vector(local_size, 1, size, MPI_DOUBLE, &cell);
    MPI_Type_create_resized(cell, 0, 1 * sizeof(double), &cell);
    MPI_Type_commit(&cell);

    // Envoi des colonnes de la matrice A et des cases du vecteur B
    MPI_Scatter(A->values, 1, col, local_a.values, m_size*local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B->values, 1, cell, local_b.values, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Décomposition LU
    // TODO

    // Descente-remontée LU
    MPI_matrix_solve(&local_a, &local_x, &local_b);

    // Récupération du résultat dans le vecteur X
    MPI_Gather(local_x.values, local_size, MPI_DOUBLE, X->values, 1, cell, 0, MPI_COMM_WORLD);

    matrix_free(&local_a);
    vector_free(&local_x);
    vector_free(&local_b);
}

void check_solve_parallel() {
    struct matrix m;
    struct vector x, b;
    int m_size = 8;
    matrix_init(&m, m_size, m_size);
    vector_init(&x, m_size);
    vector_init(&b, m_size);

    if (rank == 0) {
        /*for (int i = 0 ; i < m.nb_cols ; i++) {
            for (int j = 0 ; j < m.nb_rows ; j++) {
                matrix_set(&m, j, i, i%size + 1);
            }
        }*/
        for (int i = 0 ; i < m.nb_cols ; i++) {
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
        }

        vector_load(&b);
        printf("\nVECTOR B\n");
        vector_show(&b);

        /*1.000000 1.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
        1.000000 -1.000000 2.000000 0.000000 0.000000 0.000000 0.000000 0.000000
        1.000000 1.000000 -2.000000 3.000000 0.000000 0.000000 0.000000 0.000000
        1.000000 1.000000 1.000000 -3.000000 4.000000 0.000000 0.000000 0.000000
        1.000000 1.000000 1.000000 1.000000 -4.000000 5.000000 0.000000 0.000000
        1.000000 1.000000 1.000000 1.000000 1.000000 -5.000000 6.000000 0.000000
        1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 -6.000000 7.000000
        1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 -7.000000*/

        printf("\nMatrix A\n");
        matrix_show(&m);
    }

    solve_parallel(&m, &x, &b);

    if (rank == 0) {
        printf("\nSolution X is \n");
        vector_show(&x);
    }

    matrix_free(&m);
    vector_free(&x);
    vector_free(&b);
}

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        /**check_lu_example();
        check_trsm_example();*/
        //check_extract_matrix();
        //check_solving_example(8);
        /*timer_begin();
        check_lu_example();
        double t = timer_end();
        printf("Time : %lf microseconds\n", t);*/
    }

    check_solve_parallel();

    MPI_Finalize();

    return EXIT_SUCCESS;
}
