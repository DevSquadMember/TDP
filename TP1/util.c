#include <stdio.h>
#include <stdlib.h>
#include "util.h"

/* Allouer un vecteur d'entiers de taille n */
int* s_vector_create(int n) {
  return (int*) malloc(n * sizeof(int));
}

double* d_vector_create(int n) {
  return (double*) malloc(n * sizeof(double));
}

void s_vector_fill(int* vector, int n, int value) {
  int i;
  for (i = 0 ; i < n ; i++) {
    vector[i] = value;
  }
}

void d_vector_fill(double* vector, int n, double value) {
  int i;
  for (i = 0 ; i < n ; i++) {
    vector[i] = value;
  }
}

/* Libérer le vecteur */
void s_vector_free(int* vector) {
  free(vector);
}

void d_vector_free(double* vector) {
  free(vector);
}

/* Allouer une matrice d'entiers de taille n * m */
int* s_matrix_create(int n, int m) {
  return (int*) calloc(n * m, sizeof(int));
}

double* d_matrix_create(int n, int m) {
  return (double*) calloc(n * m, sizeof(double));
}

void d_matrix_fill(double* matrix, int n, int m, double value) {
  int i, j, col;
  for (j = 0 ; j < n ; j++) {
    col = j*m;
    for (i = 0 ; i < m ; i++) {
      matrix[col + i] = value;
    }
  }
}

/* Libérer la matrice */
void s_matrix_free(int* matrix) {
  free(matrix);
}

void d_matrix_free(double* matrix) {
  free(matrix);
}

void affiche(int m, int n, double* a, int lda, FILE* flux) {
  int i, j;
  for (i = 0 ; i < m ; i++) {
    for (j = 0 ; j < n ; j++) {
      fprintf(flux, "%f ", a[j * lda + i]);
    }
    fprintf(flux, "\n");
  }
}
