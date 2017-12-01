#ifndef TDP_UTILS_H
#define TDP_UTILS_H

struct matrix {
    int nb_cols;
    int nb_rows;
    double* values;
};

struct vector {
    int nb_values;
    double* values;
};

void matrix_init(struct matrix* matrix, int nb_cols, int nb_rows);
void matrix_show(struct matrix* matrix);
void matrix_free(struct matrix* matrix);

void vector_init(struct vector* vector, int nb_values);
void vector_show(struct vector* vector);
void vector_free(struct vector* vector);

#endif //TDP_UTILS_H
