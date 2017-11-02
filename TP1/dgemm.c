#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"

#define BLOC_SIZE 100

/*
 *  * CBLAS_ORDER Order : ordre par ligne ou colonne
 *   * CBLAS_TRANSPOSE TransA (char) : N/n pour normal - T/t - C/c pour la transposée
 *    * CBLAS_TRANSPOSE TransB (char) : N/n pour norlal - T/t - C/c pour la transposée
 *     * M : nombre de lignes de A et C
 *      * N : nombre de colonnes de B et C
 *       * K : nombre de colonnes de A et nombre de lignes de B
 *        * ALPHA : scalaire alpha (1)
 *         * A : matrice A
 *          * lda : décalage entre chaque ligne de A
 *           * B : matrice B
 *            * ldb : décalage entre chaque ligne de B
 *             * BETA : scalaire beta (0)
 *              * C : matrice C
 *               * ldc : décalage entre chaque ligne de C
 *                *
 *                 * TransA :
 *                  * - CblasNoTrans = (i, j, k) --> size * size * size * 2
 *                   * - CblasTrans = (k, i, j) --> size * size * size * 2
 *                    * - CblasConjTrans = (j, i, k) --> size * size * size * 2
 *                     */
void cblas_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
	int n, m, k, val, b_pos, a_pos, c_pos;
	if (TransA == CblasNoTrans) {
		if (lda == 1) { 
			for (m = 0 ; m < M ; m++) {
				for (n = 0 ; n < N ; n++) {
					c_pos = n*ldc + m;
					val = C[c_pos];
					b_pos = n*ldb;
					for (k = 0 ; k < K ; k++) {
						val += A[k + m] * B[b_pos + k];
					} 
					if (alpha == 1) {
						C[c_pos] = val;
					} else {
						C[c_pos] = alpha*val;
					}
				}
			}
		} else {
			for (m = 0 ; m < M ; m++) {
				for (n = 0 ; n < N ; n++) {
					c_pos = n*ldc + m;
					val = C[c_pos];
					b_pos = n*ldb;
					for (k = 0 ; k < K ; k++) {
						val += A[k*lda + m] * B[b_pos + k];
					} 
					if (alpha == 1) {
						C[c_pos] = val;
					} else {
						C[c_pos] = alpha*val;
					}
				}
			}
		}
	} else if (TransA == CblasTrans) {
		if (alpha == 1) {
			for (k = 0 ; k < K ; k++) {
				for (m = 0 ; m < M ; m++) {
					a_pos = k*lda + m;
					for (n = 0 ; n < N ; n++) {
						C[n*ldc + m] += A[a_pos] * B[n*ldb + k];
					} 
				}
			}
		} else {
			for (k = 0 ; k < K ; k++) {
				for (m = 0 ; m < M ; m++) {
					a_pos = k*lda + m;
					for (n = 0 ; n < N ; n++) {
						C[n*ldc + m] += alpha*(A[a_pos] * B[n*ldb + k]);
					} 
				}
			}
		}
	} else {
		if (lda == 1) { 
			for (n = 0 ; n < N ; n++) {
				b_pos = n*ldb;
				a_pos = n*ldc;
				for (m = 0 ; m < M ; m++) {
					c_pos = a_pos + m;
					val = C[c_pos];
					if (alpha == 1) {
						for (k = 0 ; k < K ; k++) {
							val += A[k + m] * B[b_pos + k];
						} 
					} else {
						for (k = 0 ; k < K ; k++) {
							val += alpha * A[k + m] * B[b_pos + k];
						} 
					}
					C[c_pos] = val;
				}
			}
		} else {
			for (n = 0 ; n < N ; n++) {
				b_pos = n*ldb;
				a_pos = n*ldc;
				for (m = 0 ; m < M ; m++) {
					c_pos = a_pos + m;
					val = C[c_pos];
					if (alpha == 1) {
						for (k = 0 ; k < K ; k++) {
							val += A[k*lda + m] * B[b_pos + k];
						} 
					} else {
						for (k = 0 ; k < K ; k++) {
							val += alpha * A[k*lda + m] * B[b_pos + k];
						} 
					}
					C[c_pos] = val;
				}
			}
		}
	}
}

void cblas_dgemm_scalar_opti(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
	int n, m, k, val, b_pos, a_pos, c_pos;
	if (alpha == 1) {
		for (k = 0 ; k < K ; k++) {
			for (m = 0 ; m < M ; m++) {
				a_pos = k*lda + m;
				for (n = 0 ; n < N ; n++) {
					C[n*ldc + m] += A[a_pos] * B[n*ldb + k];
				} 
			}
		}
	} else {
		for (k = 0 ; k < K ; k++) {
			for (m = 0 ; m < M ; m++) {
				a_pos = k*lda + m;
				for (n = 0 ; n < N ; n++) {
					C[n*ldc + m] += alpha*(A[a_pos] * B[n*ldb + k]);
				} 
			}
		}
	}
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
	if (M <= BLOC_SIZE && N <= BLOC_SIZE && K <= BLOC_SIZE) {
		cblas_dgemm_scalar_opti(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
		return;
	}

	int i, j, k;
	int i_pos, j_pos, k_pos;
	int nb_m_blocks, nb_k_blocks, nb_n_blocks;
	int m_size, k_size, n_size;
	int last_m_size, last_k_size, last_n_size;
	int last_m_block, last_k_block, last_n_block;
	int a_index, b_index, c_index;
	nb_m_blocks = ceil(M/BLOC_SIZE);
	nb_k_blocks = ceil(K/BLOC_SIZE);
	nb_n_blocks = ceil(N/BLOC_SIZE);
	last_m_size = M - BLOC_SIZE*(nb_m_blocks-1);
	last_k_size = K - BLOC_SIZE*(nb_k_blocks-1);
	last_n_size = N - BLOC_SIZE*(nb_n_blocks-1);

	last_m_block = nb_m_blocks - 1;
	last_n_block = nb_n_blocks - 1;
	last_k_block = nb_k_blocks - 1;

	for (j = 0 ; j < nb_n_blocks ; j++) {
		if (j < last_n_block) {
			n_size = BLOC_SIZE;
		} else {
			n_size = last_n_size;
		}
		j_pos = j*BLOC_SIZE;
		for (i = 0 ; i < nb_m_blocks ; i++) {
			if (i < last_m_block) {
				m_size = BLOC_SIZE;
			} else {
				m_size = last_m_size;
			}
			i_pos = i*BLOC_SIZE;
			c_index = i_pos + j_pos*M;
			for (k = 0 ; k < nb_k_blocks ; k++) {
				k_pos = k*BLOC_SIZE;
				a_index = i_pos + k_pos*M;
				b_index = k_pos + j_pos*K;
				if (k < last_k_block) {	
					k_size = BLOC_SIZE;
				} else {
					k_size = last_k_size;
				}
				cblas_dgemm_scalar_opti(Order, TransA, TransB, m_size, n_size, k_size, alpha, &(A[a_index]), lda, &(B[b_index]), ldb, beta, &(C[c_index]), ldc);
			}
		} 
	}
}

