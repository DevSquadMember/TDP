#include <stdlib.h>
#include <stdio.h>
#include "lib_utils.h"

#define MAX 100000
#define MIN (-100000)

/** MATRIX **/

void matrix_init(struct matrix* matrix, int nb_cols, int nb_rows) {
    matrix->nb_cols = nb_cols;
    matrix->nb_rows = nb_rows;
    matrix->ld = nb_rows;
    matrix->values = malloc(sizeof(double) * nb_cols * nb_rows);
}

void matrix_generate(struct matrix* matrix) {
    for (int i = 0 ; i < matrix->nb_rows ; i++) {
        for (int j = 0; j < matrix->nb_cols; j++) {
            matrix->values[j * matrix->ld + i] =
                    ((rand() / (double) RAND_MAX) * (2 * MAX) + (rand() / (double) RAND_MAX)) + MIN;
        }
    }
}

void matrix_copy(struct matrix* src, struct matrix* dest) {
    matrix_init(dest, src->nb_cols, src->nb_rows);
    dest->ld = src->ld;

    for (int i = 0 ; i < src->nb_rows ; i++) {
        for (int j = 0; j < src->nb_cols; j++) {
            dest->values[j * dest->ld + i] = src->values[j * src->ld + i];
        }
    }
}

void matrix_sub(struct matrix* matrix, int nb_cols, int nb_rows, int i, int j, struct matrix* initial) {
    matrix->nb_cols = nb_cols;
    matrix->nb_rows = nb_rows;
    matrix->ld = initial->ld;
    matrix->values = &(initial->values[i * initial->ld + j]);
}

double matrix_get(struct matrix* matrix, int i, int j) {
    return matrix->values[j * matrix->ld + i];
}

void matrix_set(struct matrix* matrix, int i, int j, double value) {
    matrix->values[j * matrix->ld + i] = value;
}

double matrix_setadd(struct matrix* matrix, int i, int j, double value) {
    double* res = &(matrix->values[j * matrix->ld + i]);
    (*res) += value;
    return *res;
}

double matrix_setsub(struct matrix* matrix, int i, int j, double value) {
    double* res = &(matrix->values[j * matrix->ld + i]);
    (*res) -= value;
    return *res;
}

double matrix_setmul(struct matrix* matrix, int i, int j, double value) {
    double* res = &(matrix->values[j * matrix->ld + i]);
    (*res) *= value;
    return *res;
}

double matrix_setdiv(struct matrix* matrix, int i, int j, double value) {
    double* res = &(matrix->values[j * matrix->ld + i]);
    (*res) /= value;
    return *res;
}

void matrix_show(struct matrix* matrix) {
    for (int i = 0 ; i < matrix->nb_rows ; i++) {
        for (int j = 0 ; j < matrix->nb_cols ; j++) {
            printf("%lf ", matrix->values[j*matrix->ld + i]);
        }
        printf("\n");
    }
}

void matrix_free(struct matrix* matrix) {
    if (matrix->ld == matrix->nb_rows) {
        free(matrix->values);
    }
}

/** VECTOR **/

void vector_init(struct vector* vector, int nb_values) {
    vector->nb_values = nb_values;
    vector->ld = 1;
    vector->values = malloc(sizeof(double) * nb_values);
}

void vector_generate(struct vector* vector) {
    for (int i = 0 ; i < vector->nb_values ; i++) {
        vector->values[vector->ld * i] = ((rand() / (double) RAND_MAX) * (2 * MAX) + (rand() / (double) RAND_MAX)) + MIN;
    }
}

void vector_copy(struct vector* src, struct vector* dest) {
    vector_init(dest, src->nb_values);
    dest->ld = src->ld;

    for (int i = 0 ; i < src->nb_values ; i++) {
        dest->values[i * dest->ld] = src->values[i * src->ld];
    }
}

void vector_set(struct vector* vector, int i, double value) {
    vector->values[vector->ld * i] = value;
}

double vector_setsub(struct vector* vector, int i, double value) {
    double* res = &(vector->values[vector->ld * i]);
    (*res) -= value;
    return *res;
}

double vector_get(struct vector* vector, int i) {
    return vector->values[vector->ld * i];
}

void vector_show(struct vector* vector) {
    for (int i = 0 ; i < vector->nb_values ; i++) {
        printf("%lf ", vector->values[vector->ld * i]);
    }
    printf("\n");
}

void vector_free(struct vector* vector) {
    free(vector->values);
}
