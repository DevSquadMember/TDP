#ifndef TDP_LIBUTILS_H
#define TDP_LIBUTILS_H

#include <stdbool.h>

// Structure matrice
struct matrix {
    int nb_cols;
    int nb_rows;

    int ld;

    double* values;
};

// Structure vecteur
struct vector {
    int nb_values;

    int ld;

    double* values;
};

// Initialisation de la matrice
void matrix_init(struct matrix* matrix, int nb_cols, int nb_rows);
// Génération aléatoire d'une matrice
void matrix_generate(struct matrix* matrix);
// Copie de la matrice (pas besoin de faire matrix_init)
void matrix_copy(struct matrix* src, struct matrix* dest);
// Création d'une sous-matrice (éléments partagés)
void matrix_sub(struct matrix* matrix, int nb_cols, int nb_rows, int i, int j, struct matrix* initial);
// Récupérer l'élément de la matrice
double matrix_get(struct matrix* matrix, int i, int j);
// Affecter l'élément de la matrice
void matrix_set(struct matrix* matrix, int i, int j, double value);
// Additionner une valeur à l'élément de la matrice
double matrix_setadd(struct matrix* matrix, int i, int j, double value);
// Soustraire une valeur à l'élément de la matrice
double matrix_setsub(struct matrix* matrix, int i, int j, double value);
// Multiplier une valeur à l'élément de la matrice
double matrix_setmul(struct matrix* matrix, int i, int j, double value);
// Diviser une valeur à l'élément de la matrice
double matrix_setdiv(struct matrix* matrix, int i, int j, double value);
// Afficher la matrice
void matrix_show(struct matrix* matrix);
// Libérer la mémoire allouée par la matrice
void matrix_free(struct matrix* matrix);

// Initialisation du vecteur
void vector_init(struct vector* vector, int nb_values);
// Génération aléatoire d'un vecteur
void vector_generate(struct vector* vector);
// Copie du vecteur (pas besoin de faire vector_init)
void vector_copy(struct vector* src, struct vector* dest);
// Affecter l'élément du vecteur
void vector_set(struct vector* vector, int i, double value);
// Soustraire une valeur à l'élément du vecteur
double vector_setsub(struct vector* vector, int i, double value);
// Récupérer l'élément du vecteur
double vector_get(struct vector* vector, int i);
// Afficher le vecteur
void vector_show(struct vector* vector);
// Libérer la mémoire allouée par le vecteur
void vector_free(struct vector* vector);

#endif //TDP_LIBUTILS_H
