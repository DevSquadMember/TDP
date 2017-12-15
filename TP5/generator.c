#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"

#define FILENAME "particles.txt"
#define NB_PARTICULES 12

int main(int argc, char** argv) {
    int nb_particules = NB_PARTICULES;
    char* filename;

    if (argc > 1) {
        nb_particules = atoi(argv[1]);
        if (argc > 2) {
            filename = argv[2];
        } else {
            filename = FILENAME;
        }
    } else {
        filename = FILENAME;
    }

    printf("Génération d'un fichier de particules\n");
    printf("-- Nombre de particules : %d\n", nb_particules);
    printf("-- Nom du fichier généré : %s\n", filename);

    generate_particules(nb_particules, filename);

    return EXIT_SUCCESS;
}

