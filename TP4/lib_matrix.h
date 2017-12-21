#ifndef TDP_LIB_MATRIX_H
#define TDP_LIB_MATRIX_H

struct matrix;
struct vector;

/**
 * Factorisation LU in-place :
 * la matrice A est modifiée sur place pour accueillir les matrices :
 * - L : triangulaire inférieure avec diagonale identité (la diagonale n'est pas stockée)
 * - U : triangulaire supérieure
 *
 * @param A la matrice A à factoriser en LU
 */
void dgetf2(struct matrix* A);

long long dgetf2_complexity(int size);

/**
 * Factorisation LU in-place en parallèle:
 * la matrice A est modifiée sur place pour accueillir les matrices :
 * - L : triangulaire inférieure avec diagonale identité (la diagonale n'est pas stockée)
 * - U : triangulaire supérieure
 *
 * @param A la sous-matrice A du processeur à factoriser en LU
 */
void MPI_dgetf2(struct matrix* A);

void dtrsm(char side, char uplo, char trans, char unit, int m, int n, int alpha, struct matrix* A, struct matrix* B);
void dgemm(char transa, char transb, int M, int N, int K, double alpha, struct matrix* A, struct matrix* B, double beta, struct matrix* C);
void dgetrf(struct matrix* A);

/**
 * Descente-remontée pour calculer le vecteur solution X
 *
 * @param A matrice stockant les matrices L et U (donc matrice A factorisée)
 * @param X vecteur solution
 * @param B vecteur résultat du système
 */
void matrix_solve(struct matrix* A, struct vector* X, struct vector* B);

long long matrix_solve_complexity(int size);

/**
 * Descente-remontée parallèle pour calculer le vecteur solution X
 *
 * @param A sous-matrice du processeur stockant les matrices L et U (donc matrice A factorisée)
 * @param X sous-vecteur solution du processeur
 * @param B sous-vecteur résultat du système du processeur
 */
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
void solve_sequential_dgetf2(struct matrix* A, struct vector* X, struct vector* B);

void solve_sequential_dgetrf(struct matrix* A, struct vector* X, struct vector* B);

#endif //TDP_LIB_MATRIX_H
