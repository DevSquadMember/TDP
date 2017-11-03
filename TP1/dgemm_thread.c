#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "dgemm_thread.h"

#define DEFAULT 1

struct thread_params {
	int M, N, K;
	double alpha;
	double* A;
	int lda;
	double* B;
	int ldb;
	double beta;
	double* C;
	int ldc;
	int pos;
	int nb_threads;
};

/*
Exemple de répartition pour trois threads :

     0 0 0
     1 1 1
     2 2 2
0 1 2
0 1 2
0 1 2
Complexité : M*N*K*2 (ou 3 si alpha et beta)
*/
void gemm_thread_block(struct thread_params* params) {
	int b, i, j, k;
	int k_size;
	int k_begin, k_end;
	int a_pos, b_pos, c_pos;

	k_size = floor(params->K/params->nb_threads);

	k_begin = k_size * params->pos;
	if (params->pos == params->nb_threads -1)
		k_size = params->K - (params->nb_threads-1)*k_size;
	k_end = k_begin + k_size;
	
	if (params->alpha == 1) {
		for (k = k_begin ; k < k_end ; k++) {
			a_pos = k*params->lda;
			c_pos = 0;
			for (j = 0 ; j < params->N ; j++) {
				b_pos = 0;
				for (i = 0 ; i < params->M ; i++) {
					params->C[c_pos + i] += params->A[a_pos + i] * params->B[b_pos + k];
				}
				c_pos += params->ldc;
				b_pos += params->ldb;
			}
			a_pos += params->lda;
		}
	} else {
		for (k = k_begin ; k < k_end ; k++) {
			a_pos = k*params->lda;
			c_pos = 0;
			for (j = 0 ; j < params->N ; j++) {
				b_pos = 0;
				for (i = 0 ; i < params->M ; i++) {
					params->C[c_pos + i] += params->alpha * params->A[a_pos + i] * params->B[b_pos + k];
				}
				c_pos += params->ldc;
				b_pos += params->ldb;
			}
			a_pos += params->lda;
		}
	}

	pthread_exit(NULL);
}


/*
Exemple de répartition pour trois threads :

     0 0 0   1 1 1   2 2 2
     0 0 0   1 1 1   2 2 2
     0 0 0   1 1 1   2 2 2
0 0 0
1 1 1
2 2 2
Complexité : M*N*K*2 (ou 3 si alpha et beta)
*/
void gemm_thread_line(struct thread_params* params) {
	int b, i, j, k;
	int m_size;
	int m_begin, m_end, n_begin, n_end;
	int a_pos, b_pos, c_pos;

	m_size = floor(params->M/params->nb_threads);

	m_begin = m_size * params->pos;
	if (params->pos == params->nb_threads -1) {
		m_size = params->M - (params->nb_threads-1)*m_size;
	}
	m_end = m_begin + m_size;
	
	if (params->alpha == 1) {
		for (k = 0 ; k < params->K ; k++) {
			a_pos = 0;
			c_pos = 0;
			for (j = 0 ; j < params->N ; j++) {
				b_pos = k;
				for (i = m_begin ; i < m_end ; i++) {
					params->C[c_pos + i] += params->A[a_pos + i] * params->B[b_pos];
				}
				c_pos += params->ldc;
				b_pos += params->ldb;
			}
			a_pos += params->lda;
		}
	} else {
		for (k = 0 ; k < params->K ; k++) {
			a_pos = 0;
			c_pos = 0;
			for (j = 0 ; j < params->N ; j++) {
				b_pos = k;
				for (i = m_begin ; i < m_end ; i++) {
					params->C[c_pos + i] += params->alpha * params->A[a_pos + i] * params->B[b_pos];
				}
				c_pos += params->ldc;
				b_pos += params->ldb;
			}
			a_pos += params->lda;
		}
	}

	pthread_exit(NULL);
}

void cblas_dgemm_thread(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc) {
	int i, nb_threads;

  /* NOMBRE DE THREADS DANS LA VARIABLE 'MYLIB_NUM_THREADS' */
	char* num = getenv("MYLIB_NUM_THREADS");
	if (num != NULL) {
		nb_threads = atoi(num);
	} else {
		nb_threads = DEFAULT;
	}

	pthread_t threads[nb_threads];
	struct thread_params params[nb_threads];

	for (i = 0 ; i < nb_threads ; i++) {
		params[i].M = M;
		params[i].N = N;
		params[i].K = K;
		params[i].alpha = alpha;
		params[i].A = A;
		params[i].lda = lda;
		params[i].B = B;
		params[i].ldb = ldb;
		params[i].beta = beta;
		params[i].C = C;
		params[i].ldc = ldc;
		params[i].pos = i;
		params[i].nb_threads = nb_threads;
		pthread_create(&(threads[i]), NULL, (void *(*) (void *))gemm_thread_block, &(params[i]));
		/*pthread_create(&(threads[i]), NULL, (void *(*) (void *))gemm_thread_line, &(params[i]));*/
	}

	for (i = 0 ; i < nb_threads ; i++) {
		pthread_join(threads[i], NULL);
	}
}