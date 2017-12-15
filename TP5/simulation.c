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

void fill_planets(planet* planets, int start, int end, planet *bloc) {
    int j = 0;
    for (int i = start ; i < end ; i++) {
        copy_planet(&(planets[i]), &(bloc[j]));
        j++;
    }
}

void launch_sequential_simulation_blocs(int nb_iterations, int rendering, char* filename) {
    int i, j;
    double dtmin, dtmin1, dtmin2;
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

    // Chargement des planètes
    planet myplanets[nb_planets];
    parser_load(myplanets, 0, nb_planets, rendering != 0);

    // Chargement des deux blocs de planètes
    int nb_planets_bloc = nb_planets/2;
    planet bloc1[nb_planets_bloc];
    planet bloc2[nb_planets_bloc];
    fill_planets(myplanets, 0, nb_planets_bloc, bloc1);
    fill_planets(myplanets, nb_planets_bloc, nb_planets, bloc2);

    /** CALCUL DU BLOC TOTAL **/
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
        printf("dtmin is %lf\n", dtmin);
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

    /** CALCUL DES DEUX BLOCS **/
    /*point *forcebuf1 = malloc(nb_planets_bloc * sizeof(point));
    point *forcebuf2 = malloc(nb_planets_bloc * sizeof(point));
    double* dmin1 = malloc(nb_planets_bloc * sizeof(double));
    double* dmin2 = malloc(nb_planets_bloc * sizeof(double));
    for (i = 0; i < nb_planets_bloc; i++) {
        dmin1[i] = MAX_DOUBLE;
        dmin2[i] = MAX_DOUBLE;
    }
    if (generate_graph)
        save_seq(bloc1, nb_planets_bloc);
    for (i = 0; i < nb_iterations; i++) {
        for (j = 0; j < nb_planets_bloc; j++) {
            forcebuf1[j].x = 0;
            forcebuf1[j].y = 0;
            forcebuf2[j].x = 0;
            forcebuf2[j].y = 0;
            dmin1[j] = MAX_DOUBLE;
            dmin2[j] = MAX_DOUBLE;
        }
        // Calcul des forces du bloc 1
        calcul_force_own(bloc1, nb_planets_bloc, forcebuf1, dmin1);

        // Calcul des forces du bloc 2
        calcul_force_own(bloc2, nb_planets_bloc, forcebuf2, dmin2);

        for (j = 0; j < nb_planets_bloc; j++) {
            printf("Forcebuf 1 [%d] = (%lf, %lf)\n", j, forcebuf1[j].x, forcebuf1[j].y);
            printf("Forcebuf 2 [%d] = (%lf, %lf)\n", j, forcebuf2[j].x, forcebuf2[j].y);
        }

        // Calcul des intéractions
        calcul_force_both(bloc1, bloc2, nb_planets_bloc, forcebuf1, forcebuf2, dmin1, dmin2);

        dtmin1 = calcul_dtmin(bloc1, forcebuf1, dmin1, nb_planets_bloc);
        dtmin2 = calcul_dtmin(bloc2, forcebuf2, dmin2, nb_planets_bloc);
        printf("dtmin1 is %lf and dtmin2 is %lf\n", dtmin1, dtmin2);
        dtmin = dtmin1 < dtmin2 ? dtmin1 : dtmin2;
        calcul_newpos(bloc1, forcebuf1, nb_planets_bloc, dtmin);
        calcul_newpos(bloc2, forcebuf2, nb_planets_bloc, dtmin);

        if (generate_graph)
            save_seq(bloc1, nb_planets_bloc);
    }
    if (generate_graph) {
        save_close();
        render_seq(nb_planets_bloc, "Calcul séquentiel");
    }

    free(forcebuf1);
    free(forcebuf2);
    free(dmin1);
    free(dmin2);*/
}
