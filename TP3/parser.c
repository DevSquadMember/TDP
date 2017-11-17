#include <stdio.h>
#include <stdlib.h>
#include "parser.h"

void getLine(FILE *file, char *buffer, int buffer_size) {
    while (fgets(buffer, buffer_size, file) != NULL && buffer[0] == '%');
}

void parse(char* filename, double** matrix, int* matrix_size) {
    int buffer_size = 255;
    char buffer[buffer_size];
    FILE* file = fopen(filename, "r");

    // Check if the file can be read
    if (file == NULL) {
        printf("Le fichier %s n'a pas pu être ouvert en lecture\n", filename);
        return;
    }

    getLine(file, buffer, buffer_size);
    sscanf(buffer, "%d", matrix_size);

    *matrix = realloc(*matrix, sizeof(double) * (*matrix_size) * (*matrix_size));

    double value;
    for (int i = 0 ; i < *matrix_size ; i++) {
        for (int j = 0 ; j < *matrix_size ; j++) {
            if (!fscanf(file, "%lf", &value))
                break;
            (*matrix)[j*(*matrix_size)+i] = value;
            printf("%lf ",(*matrix)[j*(*matrix_size)+i]);
        }
        printf("\n");
    }
    fclose(file);
}