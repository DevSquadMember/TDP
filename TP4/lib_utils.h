#ifndef TDP_LIBUTILS_H
#define TDP_LIBUTILS_H

#include <stdbool.h>

struct matrix {
    int nb_cols;
    int nb_rows;

    int ld;

    double* values;
};

struct vector {
    int nb_values;

    int ld;

    double* values;
};

void matrix_init(struct matrix* matrix, int nb_cols, int nb_rows);
void matrix_generate(struct matrix* matrix);
void matrix_copy(struct matrix* src, struct matrix* dest);
void matrix_sub(struct matrix* matrix, int nb_cols, int nb_rows, int i, int j, struct matrix* initial);
double matrix_get(struct matrix* matrix, int i, int j);
void matrix_set(struct matrix* matrix, int i, int j, double value);
double matrix_setadd(struct matrix* matrix, int i, int j, double value);
double matrix_setsub(struct matrix* matrix, int i, int j, double value);
double matrix_setmul(struct matrix* matrix, int i, int j, double value);
double matrix_setdiv(struct matrix* matrix, int i, int j, double value);
void matrix_show(struct matrix* matrix);
void matrix_free(struct matrix* matrix);

void vector_init(struct vector* vector, int nb_values);
void vector_generate(struct vector* vector);
void vector_copy(struct vector* src, struct vector* dest);
void vector_set(struct vector* vector, int i, double value);
double vector_setsub(struct vector* vector, int i, double value);
double vector_get(struct vector* vector, int i);
void vector_show(struct vector* vector);
void vector_free(struct vector* vector);

#endif //TDP_LIBUTILS_H
