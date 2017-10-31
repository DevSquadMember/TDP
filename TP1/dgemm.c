#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"

/**
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
 *                  * - CblasNoTrans = (i, j, k)
 *                   * - CblasTrans = (k, i, j)
 *                    * - CblasConjTrans = (j, i, k)
 *                     */
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
	int n, m, k, val, b_pos, a_pos, c_pos;
    if (TransA == CblasNoTrans) {
    	if (lda == 1) { 
			for (m = 0 ; m < M ; m++) {
				for (n = 0 ; n < N ; n++) {
					val = 0;
					b_pos = n*ldb;
					for (k = 0 ; k < K ; k++) {
						val += A[k + m] * B[b_pos + k];
					} 
					C[n*ldc + m] = val;
				}
			}
		} else {
			for (m = 0 ; m < M ; m++) {
				for (n = 0 ; n < N ; n++) {
					val = 0;
					b_pos = n*ldb;
					for (k = 0 ; k < K ; k++) {
						val += A[k*lda + m] * B[b_pos + k];
					} 
					C[n*ldc + m] = val;
				}
			}
		}
	} else if (TransA == CblasTrans) {
		for (k = 0 ; k < K ; k++) {
			if (k == 0) {
				for (m = 0 ; m < M ; m++) {
					a_pos = k*lda + m;
					for (n = 0 ; n < N ; n++) {
						C[n*ldc + m] = A[a_pos] * B[n*ldb + k];
					} 
				}
			} else {
				for (m = 0 ; m < M ; m++) {
					a_pos = k*lda + m;
					for (n = 0 ; n < N ; n++) {
						C[n*ldc + m] += A[a_pos] * B[n*ldb + k];
					} 
				}
			}
		}
	} else {
		if (lda == 1) { 
			for (n = 0 ; n < N ; n++) {
				b_pos = n*ldb;
				c_pos = n*ldc;
				for (m = 0 ; m < M ; m++) {
					val = 0;
					for (k = 0 ; k < K ; k++) {
						val += A[k + m] * B[b_pos + k];
					} 
					C[c_pos + m] = val;
				}
			}
		} else {
			for (n = 0 ; n < N ; n++) {
				b_pos = n*ldb;
				c_pos = n*ldc;
				for (m = 0 ; m < M ; m++) {
					val = 0;
					for (k = 0 ; k < K ; k++) {
						val += A[k*lda + m] * B[b_pos + k];
					} 
					C[c_pos + m] = val;
				}
			}
		}
	}
}

