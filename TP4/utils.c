#include "utils.h"
#include "lib_utils.h"
#include <math.h>
#include <printf.h>

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

void check_correctness(struct vector* X) {
    double precision, max = 0.;
    // Vérification du résultat
    for (int i = 0 ; i < X->nb_values ; i++) {
        precision = fabs(vector_get(X, i) - 1.0);
        if (precision > max) {
            max = precision;
        }
    }
    printf("Erreur : %E\n", max);
}

void check_correctness_2(struct vector* X, struct vector* ref) {
    double precision, max = 0.;
    // Vérification du résultat
    for (int i = 0 ; i < X->nb_values ; i++) {
        precision = fabs(vector_get(X, i) - vector_get(ref, i));
        if (precision > max) {
            max = precision;
        }
    }
    printf("Erreur : %E\n", max);
}
