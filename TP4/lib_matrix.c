#include <math.h>
#include <printf.h>
#include "lib_matrix.h"
#include "lib_utils.h"
#include "mpi.h"
#include "perf.h"

/**
 * Factorisation LU in-place :
 * la matrice A est modifiée sur place pour accueillir les matrices :
 * - L : triangulaire inférieure avec diagonale identité (la diagonale n'est pas stockée)
 * - U : triangulaire supérieure
 *
 * @param A la matrice A à factoriser en LU
 */
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

void dtrsm(char side,char uplo, char trans,char unit,int m,int n,int alpha, struct matrix* A,int lda,struct matrix* B, int ldb){
    int i,j,k;

    if(uplo == 'L'){
        for(k=0;k<n;k++){
            for(i=0;i<m;i++){
                for(j=0;j<i;j++){
                    matrix_set(B,i,k, matrix_get(B,i,k)-matrix_get(A,i,j)*matrix_get(B,j,k));
                }
                if(unit != 'U'){
                    matrix_setdiv(B,i,k,matrix_get(A,i,i));
                }
            }
        }
    } else {
        for(k=0;k<n;k++){
            for(i=m-1;i>0;i--){
                for(j=i;j<m;j--){
                    matrix_set(B,k,i,matrix_get(B,k,i)-matrix_get(A,i,j)*matrix_get(B,k,j));
                }
                if(unit != 'U'){
                    matrix_setdiv(B,k,i,matrix_get(A,i,i));
                }
            }
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
    // TODO

    perf(&p_end);
    if (rank == 0) {
        perf_diff(&p_begin, &p_end);
        printf("Décomposition LU parallèle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
        perf_printmicro(&p_end);
    }

    perf(&p_begin);

    // Descente-remontée LU
    MPI_matrix_solve(&local_a, &local_x, &local_b);

    perf(&p_end);
    if (rank == 0) {
        perf_diff(&p_begin, &p_end);
        printf("Descente-remontée parallèle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
        perf_printmicro(&p_end);
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

void solve_sequential(struct matrix* A, struct vector* X, struct vector* B) {
    struct timeval p_begin, p_end;
    perf(&p_begin);

    // DGETF2
    dgetf2(A);

    perf(&p_end);
    perf_diff(&p_begin, &p_end);
    printf("DGETF2 séquentiel sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
    perf_printmicro(&p_end);

    perf(&p_begin);

    // Descente-remontée
    matrix_solve(A, X, B);

    perf(&p_end);
    perf_diff(&p_begin, &p_end);
    printf("Descente-remontée séquentielle sur matrice de taille %d %d : ", A->nb_rows, A->nb_cols);
    perf_printmicro(&p_end);
}
