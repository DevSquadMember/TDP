#ifndef DGEMM_THREAD_H
#define DGEMM_THREAD_H

#include "cblas.h"

void cblas_dgemm_thread(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);

#endif // "DGEMM_THREAD_H"