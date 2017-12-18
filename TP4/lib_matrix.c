#include "lib_matrix.h"
#include "utils.h"

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
