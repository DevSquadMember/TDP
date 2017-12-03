#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define A_FILENAME "a.txt"
#define B_FILENAME "b.txt"
#define SIZE 100

#define MAX 100000
#define MIN (-100000)

void generate_matrix(int size, char* filename) {
    FILE* file;

    file = fopen(filename, "w");
    fprintf(file, "%d\n", size);
    double value;

    for (int i = 0 ; i < size ; i++) {
        for (int j = 0 ; j < size ; j++) {
            value = ((rand()/(double)RAND_MAX) * (2*MAX) + (rand()/(double)RAND_MAX)) + MIN;
            fprintf(file, "%lf ", value);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

int main(int argc, char** argv) {
    int size = SIZE;

    if (argc > 1) {
        size = atoi(argv[1]);
    }

    srand((unsigned int) time(NULL));

    printf("Génération des fichiers de matrices A et B \n");
    printf("-- Taille des matrices : %d\n", size);

    generate_matrix(size, A_FILENAME);
    generate_matrix(size, B_FILENAME);

    return EXIT_SUCCESS;
}
