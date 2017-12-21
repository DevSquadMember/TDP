#ifndef TDP_UTILS_H
#define TDP_UTILS_H

#define MAX_SIZE_TO_PRINT 20
#define MATRIX_SIZE 8

struct matrix;
struct vector;

void matrix_load(struct matrix* m);
void vector_load(struct vector* v);

void check_correctness(struct vector* X);
void check_correctness_2(struct vector* X, struct vector* ref);

#endif //TDP_UTILS_H
