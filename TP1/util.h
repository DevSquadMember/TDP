#ifndef UTIL_H
#define UTIL_H

int* s_vector_create(int n);

double* d_vector_create(int n);

void s_vector_fill(int* vector, int n, int value);

void d_vector_fill(double* vector, int n, double value);

void s_vector_free(int* vector);

void d_vector_free(double* vector);

int* s_matrix_create(int n, int m);

double* d_matrix_create(int n, int m);

void s_matrix_free(int* matrix);

void d_matrix_free(double* matrix);

#endif

