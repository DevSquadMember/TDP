#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"
#include "util.h"
#include "perf.h"
#include "dgemm_thread.h"

#define MIN 50
#define MAX 1400
#define INNER 100
#define COEF 1.25
#define ADD 0

void cblas_dgemm_scalar(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
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

long compute_ddot(long size, perf_t* start, perf_t* stop) {
  int i;
  double *dx = d_vector_create(size);
  double *dy = d_vector_create(size);

  d_vector_fill(dx, size, 2);
  d_vector_fill(dy, size, 5);

  perf(start);
  for (i = 0 ; i < INNER ; ++i) {
    cblas_ddot(size, dx, 1, dy, 1); 
  }
  perf(stop);

  d_vector_free(dx);
  d_vector_free(dy);  

  return (size + size - 1) * INNER;
}

long compute_gemm(long size, perf_t* start, perf_t* stop, const enum CBLAS_TRANSPOSE trans) {
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);

  perf(start);
  cblas_dgemm_scalar(CblasColMajor, trans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  return size * size * size * 2;
}

long compute_gemm_block(long size, perf_t* start, perf_t* stop) {
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);
  d_matrix_fill(C, size, size, 0);

  perf(start);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  int nb = ceil(1.0*size/100);
  if (size <= 100) {
    return size * size * size * 2;
  } else {
    return (nb*nb*nb + size * size * size * 2);
  }
}

long compute_gemm_thread(long size, perf_t* start, perf_t* stop) {
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);
  d_matrix_fill(C, size, size, 0);

  char* num = getenv("MYLIB_NUM_THREADS");
  int nb_threads;
  if (num != NULL) {
    nb_threads = atoi(num);
  } else {
    nb_threads = 1;
  }

  perf(start);
  cblas_dgemm_thread(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  return 1.0*(size * size * size * 2)/nb_threads;
}

long compute_gemm_mkl(long size, perf_t* start, perf_t* stop) {
  double *A = d_matrix_create(size, size);
  double *B = d_matrix_create(size, size);
  double *C = d_matrix_create(size, size);

  d_matrix_fill(A, size, size, 2);
  d_matrix_fill(B, size, size, 5);
  d_matrix_fill(C, size, size, 0);

  char* num = getenv("MKL_NUM_THREADS");
  int nb_threads;
  if (num != NULL) {
    nb_threads = atoi(num);
  } else {
    nb_threads = 1;
  }

  perf(start);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1, A, size, B, size, 0, C, size);
  perf(stop);

  d_matrix_free(A);
  d_matrix_free(B);
  d_matrix_free(C);

  return 1.0*(size * size * size * 2)/nb_threads;
}

int main(int argc, char** argv) {
  perf_t start, stop;
  double performance, perf1, perf2, perf3, perf4;
  long flop;
  long size = MIN;
  FILE* file = fopen("res.dot", "w");

  while (size < MAX) { 

    /* DDOT */

    /*flop = compute_ddot(size, &start, &stop);
 *     perf_diff(&start, &stop);
 *         performance = perf_mflops(&stop, flop);
 *
 *             printf("Size : %ld - Mflops/s : %lf\n", size, performance);
 *                 fprintf(file, "%ld %lf\n", size, performance);
 *                     */
    /* GEMM */

    /*flop = compute_gemm(size, &start, &stop, CblasNoTrans);
 *     perf_diff(&start, &stop);
 *         perf1 = perf_mflops(&stop, flop);
 *
 *             flop = compute_gemm(size, &start, &stop, CblasTrans);
 *                 perf_diff(&start, &stop);
 *                     perf2 = perf_mflops(&stop, flop);
 *
 *                         flop = compute_gemm(size, &start, &stop, CblasConjTrans);
 *                             perf_diff(&start, &stop);
 *                                 perf3 = perf_mflops(&stop, flop);
 *
 *                                     printf("Size : %ld - Mflop/s : %lf - %lf - %lf\n", size, perf1, perf2, perf3);
 *                                         fprintf(file, "%ld %lf %lf %lf\n", size, perf1, perf2, perf3);
 *                                             */

    /* GEMM BLOCK */
    
    /*flop = compute_gemm_block(size, &start, &stop);
 *      perf_diff(&start, &stop);
 *           performance = perf_mflops(&stop, flop);
 *
 *                printf("Size : %ld - Mflop/s : %lf\n", size, performance);
 *                     fprintf(file, "%ld %lf\n", size, performance);*/

    /* GEMM VS GEMM BLOCK VS GEMM THREAD */

    /*flop = compute_gemm(size, &start, &stop, CblasNoTrans);
 *     perf_diff(&start, &stop);
 *         perf1 = perf_mflops(&stop, flop);
 *
 *             flop = compute_gemm(size, &start, &stop, CblasTrans);
 *                 perf_diff(&start, &stop);
 *                     perf2 = perf_mflops(&stop, flop);
 *
 *                         flop = compute_gemm(size, &start, &stop, CblasConjTrans);
 *                             perf_diff(&start, &stop);
 *                                 perf3 = perf_mflops(&stop, flop);
 *
 *                                     flop = compute_gemm_block(size, &start, &stop);
 *                                         perf_diff(&start, &stop);
 *                                             performance = perf_mflops(&stop, flop);
 *
 *                                                 flop = compute_gemm_thread(size, &start, &stop);
 *                                                     perf_diff(&start, &stop);
 *                                                         perf4 = perf_mflops(&stop, flop);
 *
 *                                                             printf("Size : %ld - Mflop/s : %lf - %lf - %lf - %lf - %lf\n", size, perf1, perf2, perf3, performance, perf4);
 *                                                                 fprintf(file, "%ld %lf %lf %lf %lf %lf\n", size, perf1, perf2, perf3, performance, perf4);*/

    /* MKL */

    /*flop = compute_gemm_mkl(size, &start, &stop);
 *     perf_diff(&start, &stop);
 *         performance = perf_mflops(&stop, flop);
 *
 *             printf("Size : %ld - Mflops/s : %lf\n", size, performance);
 *                 fprintf(file, "%ld %lf\n", size, performance);*/

    /* MKL VS THREADS */
    flop = compute_gemm_mkl(size, &start, &stop);
    perf_diff(&start, &stop);
    perf1 = perf_mflops(&stop, flop);

    flop = compute_gemm_thread(size, &start, &stop);
    perf_diff(&start, &stop);
    perf2 = perf_mflops(&stop, flop);

    printf("Size : %ld - Mflops/s : %lf %lf\n", size, perf1, perf2);
    fprintf(file, "%ld %lf %lf\n", size, perf1, perf2);
    
    size = size * COEF + ADD;
  }
  fclose(file);
  return 0;
}
