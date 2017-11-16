#include <stdio.h>
#include <stdlib.h>
#include "parser.h"

int nb_planets_in_file = 0;
FILE* file;

void getLine(FILE *file, char *buffer, int buffer_size) {
    while (fgets(buffer, buffer_size, file) != NULL && buffer[0] == '%');
}

int parser_nb_planets(char* filename) {
    int buffer_size = 255;
    char buffer[buffer_size];
    file = fopen(filename, "r");

    // Check if the file can be read
    if (file == NULL) {
        printf("Le fichier %s n'a pas pu être ouvert en lecture\n", filename);
        return -1;
    }

    getLine(file, buffer, buffer_size);
    sscanf(buffer, "%d", &nb_planets_in_file);
    return nb_planets_in_file;
}

int parser_load(planet* planets, int start, int nb_planets, int verbose) {
    int buffer_size = 255;
    char buffer[buffer_size];

    double m, px, py, vx, vy;

    for (int i = 0 ; i < start ; i++) {
        getLine(file, buffer, buffer_size);
    }

    for (int i = 0 ; i < nb_planets ; i++) {
        getLine(file, buffer, buffer_size);

        sscanf(buffer, "%lf %lf %lf %lf %lf", &m, &px, &py, &vx, &vy);
        if (nb_planets < 20 && verbose)
            printf("Planète n°%d chargée - Masse : %lf, Position : (%lf, %lf), Vitesse : (%lf, %lf)\n", start + i + 1, m, px, py, vx, vy);
        planets[i].mass = m;
        planets[i].pos.x = px;
        planets[i].pos.y = py;
        planets[i].speed.x = vx;
        planets[i].speed.y = vy;
        planets[i].acc.x = 0;
        planets[i].acc.y = 0;
    }

    fclose(file);
    return EXIT_SUCCESS;
}
