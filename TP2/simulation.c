#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "simulation.h"
#include "parser.h"
#include "saver.h"

#define MAX_SIZE_FOR_RENDERING 12
#define MIN_MASS 100

void generate_particules(int nb_particules, char* filename) {
    FILE* file;

    srand((unsigned int) time(NULL));

    file = fopen(filename, "w");
    fprintf(file, "%d\n", nb_particules);
    double mass, posx, posy, vitx, vity;

    for (int i = 0 ; i < nb_particules ; i++) {
        mass = (rand()/(double)RAND_MAX) * (nb_particules * 100 - MIN_MASS) + MIN_MASS;
        posx = (rand()/(double)RAND_MAX) * (nb_particules * 10000);
        posy = (rand()/(double)RAND_MAX) * (nb_particules * 10000);
        vitx = (rand()/(double)RAND_MAX) * (1000) - 500;
        vity = (rand()/(double)RAND_MAX) * (1000) - 500;

        fprintf(file, "%lf %lf %lf %lf %lf\n", mass, posx, posy, vitx, vity);
    }

    fclose(file);
}

void launch_sequential_simulation(int nb_iterations, int rendering, char* filename) {
    int i, j;
    double dtmin;
    int generate_graph;

    int nb_planets = parser_nb_planets(filename);

    if (rendering) {
        printf("Lancement du calcul en séquentiel\n");
        printf("-- Nombre d'itérations : %d\n", nb_iterations);
        printf("-- Fichier des planètes : %s\n", filename);
        printf("-- Nombre de planètes trouvées : %d\n", nb_planets);
        printf("-- Début des calculs\n");
    }

    generate_graph = nb_planets <= MAX_SIZE_FOR_RENDERING && rendering;

    planet myplanets[nb_planets];
    parser_load(myplanets, 0, nb_planets, rendering != 0);

    point *forcebuf = malloc(nb_planets * sizeof(point));
    double *dmin = malloc(nb_planets * sizeof(double));
    for (i = 0; i < nb_planets; i++) {
        dmin[i] = MAX_DOUBLE;
    }
    if (generate_graph)
        save_seq(myplanets, nb_planets);
    for (i = 0; i < nb_iterations; i++) {
        for (j = 0; j < nb_planets; j++) {
            forcebuf[j].x = 0;
            forcebuf[j].y = 0;
            dmin[j] = MAX_DOUBLE;
        }
        calcul_force_seq(myplanets, nb_planets, forcebuf, dmin);
        dtmin = calcul_dtmin(myplanets, forcebuf, dmin, nb_planets);
        calcul_newpos(myplanets, forcebuf, nb_planets, dtmin);
        if (generate_graph)
            save_seq(myplanets, nb_planets);
    }
    if (generate_graph) {
        save_close();
        render_seq(nb_planets, "Calcul séquentiel");
    }

    free(forcebuf);
    free(dmin);
}
