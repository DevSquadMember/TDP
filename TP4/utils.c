#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

void matrix_init(struct matrix* matrix, int nb_cols, int nb_rows) {
    matrix->nb_cols = nb_cols;
    matrix->nb_rows = nb_rows;
    matrix->values = malloc(sizeof(double) * nb_cols * nb_rows);
}

void matrix_show(struct matrix* matrix) {
    for (int i = 0 ; i < matrix->nb_rows ; i++) {
        for (int j = 0 ; j < matrix->nb_cols ; j++) {
            printf("%lf ", matrix->values[j*matrix->nb_cols + i]);
        }
        printf("\n");
    }
}

void matrix_free(struct matrix* matrix) {
    free(matrix->values);
}

void vector_init(struct vector* vector, int nb_values) {
    vector->nb_values = nb_values;
    vector->values = malloc(sizeof(double) * nb_values);
}

void vector_show(struct vector* vector) {
    for (int i = 0 ; i < vector->nb_values ; i++) {
        printf("%lf ", vector->values[i]);
    }
    printf("\n");
}

void vector_free(struct vector* vector) {
    free(vector->values);
}
