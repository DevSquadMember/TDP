#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

int main(int argc, char** argv) {
    struct matrix m;
    struct vector v;

    matrix_init(&m, 4, 4);

    printf("\nMatrice\n");
    for (int i = 0 ; i < m.nb_cols * m.nb_rows ; i++) {
        m.values[i] = i;
    }

    matrix_show(&m);

    vector_init(&v, 4);

    printf("\nVector\n");
    for (int i = 0 ; i < v.nb_values ; i++) {
        v.values[i] = i;
    }

    vector_show(&v);

    matrix_free(&m);
    vector_free(&v);

    return EXIT_SUCCESS;
}
