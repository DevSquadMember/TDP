#ifndef TDP_LIB_MATRIX_H
#define TDP_LIB_MATRIX_H

struct matrix;
struct vector;

void dgetrf(struct matrix* A);
void dgetf2(struct matrix* A);
void dgemm(char transa, char transb, int M, int N, int K,double alpha,struct matrix* A, struct matrix* B,  double beta, struct matrix* C);
void dtrsm(char side,char uplo, char trans,char unit,int m,int n,int alpha, struct matrix* A,struct matrix* B);
void matrix_solve(struct matrix* A, struct vector* X, struct vector* B);
void MPI_matrix_solve(struct matrix* A, struct vector* X, struct vector* B);

/**
 * Résolution du système AX = B en parallèle
 * @param A matrice A
 * @param X vecteur solution X
 * @param B vecteur résultat B
 */
void solve_parallel(struct matrix* A, struct vector* X, struct vector* B);

/**
 * Résolution du système AX = B en séquentiel
 * @param A matrice A
 * @param X vecteur solution X
 * @param B vecteur résultat B
 */
void solve_sequential(struct matrix* A, struct vector* X, struct vector* B);

#endif //TDP_LIB_MATRIX_H
