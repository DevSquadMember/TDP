#ifndef TDP_LIB_MATRIX_H
#define TDP_LIB_MATRIX_H

struct matrix;

void dgetf2(struct matrix* A);
void dtrsm(char side,char uplo, char trans,char unit,int m,int n,int alpha, struct matrix* A,int lda,struct matrix* B, int ldb);

#endif //TDP_LIB_MATRIX_H
