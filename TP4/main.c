#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

/*
 * dgtrf_nopiv
 *   for j = 1 to n ; += b
 *     dgetf2 - OK
 *     dtrsm(L)
 *     dgemm
 */

void dtrsm(char side,char uplo, char trans,char unit,int m,int n,int alpha, struct matrix* A,int lda,struct matrix* B, int ldb){
  int i,j,k;

  if(uplo == "L"){
    for(k=0;k<n;k++){
      for(i=0;i<m;i++){
	for(j=0;j<i;j++){
	  matrix_set(B,i,k, matrix_get(B,i,k)-matrix_get(A,i,j)*matrix_get(B,j,k));
	}
	if(unit != "U"){
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
	if(unit != "U"){
	  matrix_setdiv(B,k,i,matrix_get(A,i,i));
	}
      }
    }
  }
}

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

void matrix_load(struct matrix* m) {
    double value;
    for (int i = 0 ; i < m->nb_rows ; i++) {
        for (int j = 0 ; j < m->nb_cols ; j++) {
            if (j == 0) {
                value = 1;
            } else if (j == i+1) {
                value = i+1;
            } else {
                value = 0;
            }
            matrix_set(m, i, j, value);
        }
    }
}

void vector_load(struct vector* v) {
    for (int i = 0 ; i < v->nb_values ; i++) {
        vector_set(v, i, (i+1) % (v->nb_values) + 1);
    }
}

void check_trsm_example() {
  struct matrix m;
  struct matrix t;

  matrix_init(&m,3,6);
  m.values[0] = 1;
  m.values[1] = 1;
  m.values[2] = 2;
  
  m.values[3] = 3;
  m.values[4] = 1;
  m.values[5] = 5;

  m.values[6] = 6;
  m.values[7] = 7;
  m.values[8] = 1;
  
  m.values[9] = 1;
  m.values[10] = 2;
  m.values[11] = 3;
  
  m.values[12] = 4;
  m.values[13] = 5;
  m.values[14] = 6;
  
  m.values[15] = 7;
  m.values[16] = 8;
  m.values[17] = 9;

  matrix_sub(&t,3,3,0,3,&m);

  dtrsm("D","L","N","U",3,6,1,&m,3,&t,3);
}

void check_lu_example() {
    struct matrix m;

    matrix_init(&m, 3, 3);
    m.values[0] = 2;
    m.values[1] = 1;
    m.values[2] = 1;

    m.values[3] = 4;
    m.values[4] = 3;
    m.values[5] = 5;

    m.values[6] = 4;
    m.values[7] = 1;
    m.values[8] = 6;

    printf("\nMatrix before\n");

    matrix_show(&m);

    // DGETF2
    dgetf2(&m);

    printf("\nMatrix after\n");
    matrix_show(&m);

    matrix_free(&m);
}

void check_extract_matrix() {
    struct matrix m;
    struct matrix sub;

    matrix_init(&m, 16, 16);

    for (int i = 0 ; i < m.nb_rows ; i++) {
        for (int j = 0 ; j < m.nb_cols ; j++) {
            matrix_set(&m, i, j, 1);
        }
    }

    printf("BEFORE M\n");
    matrix_show(&m);
    matrix_sub(&sub, 4, 4, 6, 8, &m);

    for (int i = 0 ; i < sub.nb_rows ; i++) {
        for (int j = 0 ; j < sub.nb_cols ; j++) {
            matrix_set(&m, i, j, 2);
        }
    }
    printf("AFTER SUB\n");
    matrix_show(&sub);

    printf("AFTER M\n");
    matrix_show(&m);

    matrix_free(&m);
    matrix_free(&sub);
}

int main(int argc, char** argv) {
    struct matrix m;
    struct vector v;

    /*matrix_init(&m, 4, 4);

    printf("\nMatrice\n");
    matrix_load(&m);

    matrix_show(&m);

    printf("\nAfter FACTO LU\n");
    facto_lu(&m);

    matrix_show(&m);

    vector_init(&v, 4);

    printf("\nVector\n");
    vector_load(&v);

    vector_show(&v);

    matrix_free(&m);
    vector_free(&v);*/

    check_lu_example();

    //check_extract_matrix();

    return EXIT_SUCCESS;
}
