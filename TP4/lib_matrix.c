#include <math.h>
#include <printf.h>
#include "lib_matrix.h"
#include "lib_utils.h"
#include "mpi.h"
#include "perf.h"

void dgetf2(struct matrix* A) {
    for (int k = 0 ; k < A->nb_rows ; k++) {
        for (int i = k + 1 ; i < A->nb_rows ; i++) {

            double value = matrix_setdiv(A, i, k, matrix_get(A, k, k));

            for (int j = k + 1 ; j < A->nb_cols ; j++) {
                matrix_setsub(A, i, j, value * matrix_get(A, k, j));
            }
        }
    }
}

int dgetf2_complexity(int size) {
    return size * size * (size - 1)/2;
}

void MPI_dgetf2(struct matrix* A) {
    int rank, size, current;
    double diag, value;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int k = 0 ; k < A->nb_rows ; k++) {
        int col = (k-rank)/size;

        current = k%size;
        if (rank == current) {
            diag = matrix_get(A, k, col);
        }
        MPI_Bcast(&diag, 1, MPI_DOUBLE, current, MPI_COMM_WORLD);

        for (int i = k + 1 ; i < A->nb_rows ; i++) {
            if (rank == current) {
                value = matrix_setdiv(A, i, col, diag);
            }
            MPI_Bcast(&value, 1, MPI_DOUBLE, current, MPI_COMM_WORLD);

            int begin = (int) floor(1.0 * (k + size - rank) / size);
            for (int j = begin ; j < A->nb_cols ; j++) {
                matrix_setsub(A, i, j, value * matrix_get(A, k, j));
            }
        }
    }
}

void dtrsm(char side, char uplo, char trans, char unit, int m, int n, int alpha, struct matrix* A, struct matrix* B) {
    int i, j, k;

    if (uplo == 'L') {
        for (k = 0 ; k < n ; k++) {
            for (i = 0 ; i < m ; i++) {
                for (j = 0 ; j < i ; j++) {
                    matrix_setsub(B, i, k, matrix_get(A, i, j) * matrix_get(B, j, k));
                }
                if (unit != 'U') {
                    matrix_setdiv(B, i, k, matrix_get(A, i, i));
                }
            }
        }
    } else {
        for (k = 0 ; k < n ; k++) {
            for (i = m-1 ; i > 0 ; i--) {
                for (j = i ; j < m ; j++) {
                    matrix_setsub(B, k, i, matrix_get(A, i, j) * matrix_get(B, k, j));
                }
                if (unit != 'U') {
                    matrix_setdiv(B, k, i, matrix_get(A, i, i));
                }
            }
        }
    }
}

void dgemm(char transa, char transb, int M, int N, int K, double alpha, struct matrix* A, struct matrix* B, double beta, struct matrix* C) {
    for (int i = 0 ; i < M ; i++) {
        for (int j = 0 ; j < N ; j++) {
            matrix_setmul(C, i, j, beta);
            for (int k = 0 ; k < K ; k++) {
                matrix_set(C, i, j, alpha * matrix_get(A, i, k) * matrix_get(B, k, j));
            }
        }
    }
}

void dgetrf(struct matrix* A) {
    int bloc_length = 4;//10;
    int nb_blocs = (int) ceil(1. * A->nb_cols / bloc_length);
    if (nb_blocs == 1) {
        dgetf2(A);
    } else {
        printf("NB BLOCS : %d\n", nb_blocs);
        struct matrix LuBloc, TrsmVt, TrsmHz, GemmBloc;
        for (int i = 0 ; i < nb_blocs ; i++) {
            printf("ROUND %d\n", i);
            matrix_sub(&LuBloc, bloc_length, bloc_length, i*bloc_length, i*bloc_length, A);
            matrix_sub(&TrsmVt, bloc_length, A->nb_rows - bloc_length*(i+1), i*bloc_length, (i+1)*bloc_length, A);
            matrix_sub(&TrsmHz, A->nb_cols - bloc_length*(i+1), bloc_length, (i+1)*bloc_length, i*bloc_length, A);
            matrix_sub(&GemmBloc, A->nb_cols - bloc_length*(i+1), A->nb_cols - bloc_length*(i+1), (i+1)*bloc_length, (i+1)*bloc_length, A);

            printf("\nLuBloc \n");
            matrix_show(&LuBloc);

            printf("\nTrsmVt \n");
            matrix_show(&TrsmVt);

            printf("\nTrsmHz \n");
            matrix_show(&TrsmHz);

            printf("\nGemmbloc \n");
            matrix_show(&GemmBloc);

            printf("DGETF2\n");
            dgetf2(&LuBloc);

            printf("DTRSM\n");
            dtrsm('a', 'L', 'a', 'U', bloc_length, A->nb_cols - bloc_length*(i+1), 1, &LuBloc, &TrsmHz);
            printf("\nAFFFFTTTER --- TrsmHz \n");
            matrix_show(&TrsmHz);
            printf("TWO\n");
            dtrsm('a', 'U', 'a', 'U', bloc_length, A->nb_cols - bloc_length*(i+1), 1, &LuBloc, &TrsmVt);
            printf("\nAFFFFTTTER --- TrsmVt \n");
            matrix_show(&TrsmVt);

            printf("DGEMM\n");
            //dgemm('a', 'a', bloc_length, bloc_length, A->nb_rows - bloc_length*(i+1), 1, &TrsmVt, &TrsmHz, 1, &GemmBloc);

            printf("\nGemmbloc \n");
            matrix_show(&GemmBloc);
        }
    }
}

void MPI_dgetrf(struct matrix* A) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int bloc_length = 1;
    int nb_blocs = (int) ceil(1. * A->nb_cols / bloc_length);
    if (nb_blocs == 1) {
        dgetf2(A);
    } else {
        printf("NB BLOCS : %d\n", nb_blocs);
        struct matrix LuBloc, TrsmVt, TrsmHz, GemmBloc;
        for (int i = 0 ; i < nb_blocs ; i++) {
            printf("ROUND %d\n", i);
            matrix_sub(&LuBloc, bloc_length, bloc_length, i*bloc_length, i*bloc_length, A);
            matrix_sub(&TrsmVt, bloc_length, A->nb_rows - bloc_length*(i+1), i*bloc_length, (i+1)*bloc_length, A);
            matrix_sub(&TrsmHz, A->nb_cols - bloc_length*(i+1), bloc_length, (i+1)*bloc_length, i*bloc_length, A);
            matrix_sub(&GemmBloc, A->nb_cols - bloc_length*(i+1), A->nb_cols - bloc_length*(i+1), (i+1)*bloc_length, (i+1)*bloc_length, A);

            printf("\nLuBloc \n");
            matrix_show(&LuBloc);

            printf("\nTrsmVt \n");
            matrix_show(&TrsmVt);

            printf("\nTrsmHz \n");
            matrix_show(&TrsmHz);

            printf("\nGemmbloc \n");
            matrix_show(&GemmBloc);

            printf("DGETF2\n");
            dgetf2(&LuBloc);

            printf("DTRSM\n");
            dtrsm('a', 'L', 'a', 'a', bloc_length, A->nb_cols - bloc_length*(i+1), 1, &LuBloc, &TrsmHz);
            printf("\nAFFFFTTTER --- TrsmHz \n");
            matrix_show(&TrsmHz);
            printf("TWO\n");
            dtrsm('a', 'U', 'a', 'U', bloc_length, A->nb_cols - bloc_length*(i+1), 1, &LuBloc, &TrsmVt);
            printf("\nAFFFFTTTER --- TrsmVt \n");
            matrix_show(&TrsmVt);

            printf("DGEMM\n");
            //dgemm('a', 'a', bloc_length, bloc_length, A->nb_rows - bloc_length*(i+1), 1, &TrsmVt, &TrsmHz, 1, &GemmBloc);

            printf("\nGemmbloc \n");
            matrix_show(&GemmBloc);
        }
    }
}

/**
 * Descente-remontée de LU pour trouver le vecteur X solution du système
 *
 * @param A matrice A
 * @param X vecteur solution du système
 * @param B vecteur résultat du système
 */
void matrix_solve(struct matrix* A, struct vector* X, struct vector* B) {
    // descente
    for (int i = 0 ; i < B->nb_values ; i++) {
        double sum = 0.;
        for (int j = 0 ; j < i ; j++) {
            sum += matrix_get(A, i, j) * vector_get(X, j);
        }
        vector_set(X, i, vector_get(B, i) - sum);
    }

    // remontée
    for (int i = B->nb_values - 1 ; i >= 0 ; i--) {
        double sum = 0.;
        for (int j = i + 1 ; j < A->nb_cols ; j++) {
            sum += matrix_get(A, i, j) * vector_get(X, j);
        }
        vector_set(X, i, (vector_get(X, i) - sum)/matrix_get(A, i, i));
    }
}

int matrix_solve_complexity(int size) {
    return 2*size*(size-1) + 3*size;
}

int get_index(int i, int rank, int size) {
    return (int) floor(1.0 * (i - rank) / size);
}

void MPI_matrix_solve(struct matrix* A, struct vector* X, struct vector* B) {
    int rank, size, current, index;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // descente
    for (int i = 0; i < A->nb_rows; i++) {
        double local_sum = 0.;
        double sum = 0.;
        current = i % size;
        index = rank;
        for (int j = 0 ; j < i; j++) {
            if (index >= i) {
                break;
            }
            local_sum += matrix_get(A, i, j) * vector_get(X, j);
            index += size;
        }
        MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, current, MPI_COMM_WORLD);
        if (rank == current) {
            index = get_index(i, rank, size);
            vector_set(X, index, vector_get(B, index) - sum);
        }
    }

    // remontée
    for (int i = A->nb_rows - 1; i >= 0; i--) {
        double local_sum = 0.;
        double sum = 0.;
        current = i % size;
        int begin = (int) floor(1.0 * (i + size - rank) / size);
        for (int j = begin ; j < A->nb_cols; j++) {
            local_sum += matrix_get(A, i, j) * vector_get(X, j);
        }
        MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, current, MPI_COMM_WORLD);
        if (rank == current) {
            index = (int) floor((i - rank) / size);
            vector_set(X, index, (vector_get(X, index) - sum) / matrix_get(A, i, index));
        }
    }
}

void solve_parallel(struct matrix* A, struct vector* X, struct vector* B) {
    int rank, size;
    struct timeval p_begin, p_end;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    perf(&p_begin);

    // Envoi des colonnes de la matrice A et des cases du vecteur B
    MPI_Scatter(A->values, 1, col, local_a.values, m_size*local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B->values, 1, cell, local_b.values, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    perf(&p_end);
    if (rank == 0) {
        perf_diff(&p_begin, &p_end);
        printf("Scatter parallèle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
        perf_printmicro(&p_end);
    }

    perf(&p_begin);

    // Décomposition LU
    MPI_dgetf2(&local_a);

    perf(&p_end);

    if (rank == 0) {
        perf_diff(&p_begin, &p_end);
        printf("Décomposition LU parallèle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
        perf_printmicro(&p_end);
        printf("Performance : %lf flop/s\n", perf_mflops(&p_end, dgetf2_complexity(A->nb_rows)));
    }

    perf(&p_begin);

    // Descente-remontée LU
    MPI_matrix_solve(&local_a, &local_x, &local_b);

    perf(&p_end);
    if (rank == 0) {
        perf_diff(&p_begin, &p_end);
        printf("Descente-remontée parallèle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
        perf_printmicro(&p_end);
        printf("Performance : %lf flop/s\n", perf_mflops(&p_end, matrix_solve_complexity(A->nb_rows)));
    }

    perf(&p_begin);

    // Récupération du résultat dans le vecteur X
    MPI_Gather(local_x.values, local_size, MPI_DOUBLE, X->values, 1, cell, 0, MPI_COMM_WORLD);

    perf(&p_end);
    if (rank == 0) {
        perf_diff(&p_begin, &p_end);
        printf("Gather parallèle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
        perf_printmicro(&p_end);
    }

    matrix_free(&local_a);
    vector_free(&local_x);
    vector_free(&local_b);
}

void solve_sequential_dgetf2(struct matrix* A, struct vector* X, struct vector* B) {
    struct timeval p_begin, p_end;
    perf(&p_begin);

    // DGETF2
    dgetf2(A);

    perf(&p_end);
    perf_diff(&p_begin, &p_end);
    printf("DGETF2 séquentiel sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
    perf_printmicro(&p_end);
    printf("Performance : %lf flop/s\n", perf_mflops(&p_end, dgetf2_complexity(A->nb_rows)));

    ///printf("Matrix A is now :\n");
    ///matrix_show(A);

    perf(&p_begin);

    // Descente-remontée
    matrix_solve(A, X, B);

    perf(&p_end);
    perf_diff(&p_begin, &p_end);
    printf("Descente-remontée séquentielle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
    perf_printmicro(&p_end);
    printf("Performance : %lf flop/s\n", perf_mflops(&p_end, matrix_solve_complexity(A->nb_rows)));
}

void solve_sequential_dgetrf(struct matrix* A, struct vector* X, struct vector* B) {
    struct timeval p_begin, p_end;
    perf(&p_begin);

    // DGETRF
    dgetrf(A);

    perf(&p_end);
    perf_diff(&p_begin, &p_end);
    printf("DGETF2 séquentiel sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
    perf_printmicro(&p_end);

    ///printf("Matrix A is now :\n");
    ///matrix_show(A);

    perf(&p_begin);

    // Descente-remontée
    matrix_solve(A, X, B);

    perf(&p_end);
    perf_diff(&p_begin, &p_end);
    printf("Descente-remontée séquentielle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
    perf_printmicro(&p_end);
}
