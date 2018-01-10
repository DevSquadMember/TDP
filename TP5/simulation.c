#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "simulation.h"
#include "parser.h"
#include "saver.h"

#define MAX_SIZE_FOR_RENDERING 12
#define MAX_MASS 1000000000
#define MIN_MASS 100

#define THRESHOLD 10

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

void box_generate_particules(struct box* box, struct box* ref_box, int* current, int nb_particules) {
    double mass, posx, posy, vitx, vity;
    double min_x, min_y;

    min_x = box->pos.x * box->size;
    min_y = box->pos.y * box->size;

    for (int i = 0 ; i < nb_particules ; i++) {
        mass = (rand()/(double)RAND_MAX) * (MAX_MASS - MIN_MASS) + MIN_MASS;
        posx = (rand()/(double)RAND_MAX) * (box->size) + min_x;
        posy = (rand()/(double)RAND_MAX) * (box->size) + min_y;
        vitx = (rand()/(double)RAND_MAX) * (1000) - 500;
        vity = (rand()/(double)RAND_MAX) * (1000) - 500;

        box->planets[i].mass = mass;
        box->planets[i].pos.x = posx;
        box->planets[i].pos.y = posy;
        box->planets[i].speed.x = vitx;
        box->planets[i].speed.y = vity;

        copy_planet(&(box->planets[i]), &(ref_box->planets[*current]));

        *current = *current + 1;
    }
}

void generate_boxes(box* ref, box* boxes, int nb_boxes, int box_size, int nb_particules) {
    int total = 0;
    int row, current;
    // Box référente
    ref->size = box_size*nb_boxes;
    ref->pos.x = 0;
    ref->pos.y = 0;
    printf("Generating boxes\n");

    // Initialisation des Box
    for (int i = 0 ; i < nb_boxes ; i++) {
        row = i*nb_boxes;
        for (int j = 0 ; j < nb_boxes ; j++) {
            current = row + j;
            printf("Generating box n°%d - %d particules\n", current, nb_particules);
            // Initialisation de la Box
            box_init(&(boxes[current]), nb_particules);
            printf("OK\n");
            boxes[current].size = box_size;
            boxes[current].pos.x = row;
            boxes[current].pos.y = j;

            // Génération des particules dans la Box
            box_generate_particules(&(boxes[current]), ref, &total, nb_particules);
        }
    }
}

void fill_planets(planet* planets, int start, int end, planet *bloc) {
    int j = 0;
    for (int i = start ; i < end ; i++) {
        copy_planet(&(planets[i]), &(bloc[j]));
        j++;
    }
}

double compute_error(double ref, double val) {
    double error = val - ref;
    double percent = val == 0 ? 0 : 100.0*error/val;

    if (percent < 0)
        percent *= -1;

    if (error != 0.) {
        printf("Valeur attendue : %lf, valeur calculée : %lf\n", ref, val);
        printf("Le pourcentage d'erreur est de : %e%%\n", percent);
    }

    return percent;
}

/**
 * Vérifier que le calcul sur les boîtes donne bien des valeurs correspondant à celles obtenues dans la boîte de référence
 * @param ref la boîte de référence où toutes les planètes sont présentes
 * @param nb_total_planets le nombre total de planètes
 * @param boxes les boîtes de calcul
 */
void check_boxes(box* ref, int nb_total_planets, box* boxes) {
    printf("\n\tVérification de la cohérence des vecteurs de force\n");
    int current_box = -1;
    int nb_errors = 0;
    double percent;
    double max_percent = 0.;
    int j = 0;
    int nb_particules = boxes[0].nb_planets;

    for (int i = 0 ; i < nb_total_planets ; i++) {
        if (i % nb_particules == 0) {
            current_box++;
            j = 0;
        }

        percent = compute_error(ref->force[i].x, boxes[current_box].force[j].x);
        if (percent != 0.) {
            nb_errors++;
            if (percent > max_percent)
                max_percent = percent;
        }

        percent = compute_error(ref->force[i].y, boxes[current_box].force[j].y);
        if (percent != 0.) {
            nb_errors++;
            if (percent > max_percent)
                max_percent = percent;
        }

        j++;
    }
    printf("\n\tRésultat : %d erreurs trouvée(s) sur %d valeurs, avec un pourcentage d'erreur maximal de %e%%\n", nb_errors, nb_total_planets*2, max_percent);
}

/**
 * Lancement de la simulation sur plusieurs boîtes dont comparaison avec boîte de référence
 * @param nb_boxes nombre de boîtes à utiliser
 * @param nb_particules nombre de particules dans chaque boîte
 * @param size taille en km du monde (donc de la boîte de référence)
 * @param rendering 1 = affichage des infos, 0 = mode silencieux
 */
void launch_sequential_simulation_box(int nb_boxes, int nb_particules, int size, int rendering) {
    int nb_total_boxes = nb_boxes * nb_boxes;
    int nb_planets = nb_boxes * nb_boxes * nb_particules;

    if (rendering) {
        printf("Lancement de la simulation par boîtes\n");
        printf("-- Nombre de boîtes : %d\n", nb_total_boxes);
        printf("-- Nombre de planètes par boîte : %d\n", nb_particules);
        printf("-- Nombre de planètes total : %d\n", nb_planets);
        printf("-- Taille du monde (en km) : %d\n", size);
        printf("-- Début des calculs\n");
    }

    // Chargement des boîtes
    box ref;
    box* boxes = malloc(sizeof(box) * nb_total_boxes);

    box_init(&ref, nb_planets);
    generate_boxes(&ref, boxes, nb_boxes, size/nb_boxes, nb_particules);

    /** CALCUL DU BLOC TOTAL **/
    for (int i = 0; i < nb_planets; i++) {
        ref.force[i].x = 0;
        ref.force[i].y = 0;
    }
    calcul_force_own(&ref);

    /** CALCUL DES DEUX BLOCS **/
    for (int i = 0; i < nb_total_boxes; i++) {
        for (int j = 0 ; j < nb_particules ; j++) {
            boxes[i].force[j].x = 0;
            boxes[i].force[j].y = 0;
        }
    }

    // Calcul des forces du bloc 1
    calcul_force_own(&(boxes[0]));

    // Calcul des forces du bloc 2
    calcul_force_own(&(boxes[1]));

    // Calcul des intéractions
    calcul_force_two_boxes(&(boxes[0]), &(boxes[1]), THRESHOLD);

    /** VÉRIFICATION DES DONNÉES **/
    check_boxes(&ref, nb_planets, boxes);

    for (int i = 0 ; i < nb_total_boxes ; i++) {
        box_free(&boxes[i]);
    }
    free(boxes);
    box_free(&ref);
}

/** VIEILLE VERSION **/

void check_blocs(point* force_ref, int size, point* force1, point* force2) {
    printf("\n\tVérification de la cohérence des vecteurs de force\n");
    point* force_val = force1;
    int nb_errors = 0;
    double percent;
    double max_percent = 0.;
    int j = 0;
    for (int i = 0 ; i < size ; i++) {
        if (i == size/2) {
            force_val = force2;
            j = 0;
        }

        percent = compute_error(force_ref[i].x, force_val[j].x);
        if (percent != 0.) {
            nb_errors++;
            if (percent > max_percent)
                max_percent = percent;
        }

        percent = compute_error(force_ref[i].y, force_val[j].y);
        if (percent != 0.) {
            nb_errors++;
            if (percent > max_percent)
                max_percent = percent;
        }

        j++;
    }
    printf("\n\tRésultat : %d erreurs trouvée(s) sur %d valeurs, avec un pourcentage d'erreur maximal de %e%%\n", nb_errors, size*2, max_percent);
}

void launch_sequential_simulation_blocs(int rendering, char* filename) {
    int i, j;
    int generate_graph;

    int nb_planets = parser_nb_planets(filename);

    if (rendering) {
        printf("Lancement du calcul en séquentiel\n");
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
    box box1;
    box box2;
    box_init(&box1, nb_planets_bloc);
    box_init(&box2, nb_planets_bloc);

    fill_planets(myplanets, 0, nb_planets_bloc, box1.planets);
    fill_planets(myplanets, nb_planets_bloc, nb_planets, box2.planets);

    /** CALCUL DU BLOC TOTAL **/
    point *forcebuf = malloc(nb_planets * sizeof(point));
    double *dmin = malloc(nb_planets * sizeof(double));
    for (i = 0; i < nb_planets; i++) {
        dmin[i] = MAX_DOUBLE;
    }
    if (generate_graph)
        save_seq(myplanets, nb_planets);

    // Calcul
    for (j = 0; j < nb_planets; j++) {
        forcebuf[j].x = 0;
        forcebuf[j].y = 0;
        dmin[j] = MAX_DOUBLE;
    }
    calcul_force_seq(myplanets, nb_planets, forcebuf, dmin);

    if (generate_graph) {
        save_seq(myplanets, nb_planets);
        save_close();
        render_seq(nb_planets, "Calcul séquentiel");
    }

    /** CALCUL DES DEUX BLOCS **/
    point *forcebuf1 = malloc(nb_planets_bloc * sizeof(point));
    point *forcebuf2 = malloc(nb_planets_bloc * sizeof(point));
    double* dmin1 = malloc(nb_planets_bloc * sizeof(double));
    double* dmin2 = malloc(nb_planets_bloc * sizeof(double));
    for (i = 0; i < nb_planets_bloc; i++) {
        dmin1[i] = MAX_DOUBLE;
        dmin2[i] = MAX_DOUBLE;
    }
    if (generate_graph) {
        save(0, box1.planets, nb_planets_bloc);
        save(1, box2.planets, nb_planets_bloc);
    }

    // Calcul
    for (j = 0; j < nb_planets_bloc; j++) {
        forcebuf1[j].x = 0;
        forcebuf1[j].y = 0;
        forcebuf2[j].x = 0;
        forcebuf2[j].y = 0;
        dmin1[j] = MAX_DOUBLE;
        dmin2[j] = MAX_DOUBLE;
    }
    // Calcul des forces du bloc 1
    /**calcul_force_own(box1.planets, nb_planets_bloc, forcebuf1, dmin1);

    // Calcul des forces du bloc 2
    calcul_force_own(box2.planets, nb_planets_bloc, forcebuf2, dmin2);

    // Calcul des intéractions
    //calcul_force_two_boxes(&box1, &box2, THRESHOLD);

    calcul_force_other(box1.planets, box2.planets, nb_planets_bloc, forcebuf1, dmin1);
    calcul_force_other(box2.planets, box1.planets, nb_planets_bloc, forcebuf2, dmin2);

    calcul_newpos(box1.planets, forcebuf1, nb_planets_bloc, 1);
    calcul_newpos(box2.planets, forcebuf2, nb_planets_bloc, 1);**/

    if (generate_graph) {
        save(0, box1.planets, nb_planets_bloc);
        save(1, box2.planets, nb_planets_bloc);
        save_close();
        render(2, nb_planets, "Calcul 2 Blocs");
    }

    /** VÉRIFICATION DES DONNÉES **/
    check_blocs(forcebuf, nb_planets, forcebuf1, forcebuf2);

    free(forcebuf);
    free(forcebuf1);
    free(forcebuf2);
    free(dmin);
    free(dmin1);
    free(dmin2);
}
