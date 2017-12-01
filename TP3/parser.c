#include <stdio.h>
#include <stdlib.h>
#include "parser.h"

FILE* file;
const int buffer_size = 255;
char buffer[buffer_size];

void getLine(FILE *file, char *buffer, int buffer_size) {
    while (fgets(buffer, buffer_size, file) != NULL && buffer[0] == '%');
}

int parse_size(char* filename) {
    file = fopen(filename, "r");
    int value = -1;

    // Check if the file can be read
    if (file == NULL) {
        printf("Le fichier %s n'a pas pu être ouvert en lecture\n", filename);
        return value;
    }

    getLine(file, buffer, buffer_size);
    sscanf(buffer, "%d", &value);

    return value;
}

void parse_matrix(double* matrix, int matrix_size) {
    double value;
    for (int i = 0 ; i < matrix_size ; i++) {
        for (int j = 0 ; j < matrix_size ; j++) {
            if (!fscanf(file, "%lf", &value))
                break;
            matrix[j*(matrix_size)+i] = value;
            printf("%lf ", matrix[j * matrix_size + i]);
        }
        printf("\n");
    }
    fclose(file);
}