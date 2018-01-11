#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "simulation.h"
#include "parser.h"
#include "saver.h"
#include "perf.h"
#include "utils.h"

#define MAX_SIZE_FOR_RENDERING 12
#define MAX_MASS 1000000000
#define MIN_MASS 100

#define THRESHOLD 1000000

/**
 * Generate particles and storing them in a file
 * @param nb_particules the number of particles to generate
 * @param filename the name of the file where the particles whill be saved
 */
void generate_particules(int nb_particles, char* filename) {
    FILE* file;

    srand((unsigned int) time(NULL));

    file = fopen(filename, "w");
    fprintf(file, "%d\n", nb_particles);
    double mass, posx, posy, vitx, vity;

    for (int i = 0 ; i < nb_particles ; i++) {
        mass = (rand()/(double)RAND_MAX) * (nb_particles * 100 - MIN_MASS) + MIN_MASS;
        posx = (rand()/(double)RAND_MAX) * (nb_particles * 10000);
        posy = (rand()/(double)RAND_MAX) * (nb_particles * 10000);
        vitx = (rand()/(double)RAND_MAX) * (1000) - 500;
        vity = (rand()/(double)RAND_MAX) * (1000) - 500;

        fprintf(file, "%lf %lf %lf %lf %lf\n", mass, posx, posy, vitx, vity);
    }

    fclose(file);
}

/**
 * Generate particles for a box
 * @param box the box where the particles should be stored
 * @param ref_box the reference box where the particles should be copied (it represents the world)
 * @param current the current amount of particles in the world
 * @param nb_particules the amount of particles to create in the box
 */
void box_generate_particules(struct box* box, struct box* ref_box, int* current, int nb_particules, node* leaf) {
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
        if (leaf != NULL) {
            copy_planet(&(box->planets[i]), &(leaf->box.planets[i]));
        }

        *current = *current + 1;
    }

    calcul_center_mass(box);
    if (leaf != NULL) {
        calcul_center_mass(&(leaf->box));
    }
}

/**
 * Generate particles in the boxes
 * @param ref The reference box containing the world
 * @param boxes the array of boxes where to create the particles
 * @param nb_boxes the number of boxes (width or height of the grid)
 * @param box_size the size of the box
 * @param nb_particules the number of particles to create in each box
 * @param tree The tree to create for Barnet Huts
 */
void generate_boxes(box* ref, box* boxes, int nb_boxes, int box_size, int nb_particules, node* tree) {
    int total = 0;
    int row, current;
    // Box référente
    ref->size = box_size*nb_boxes;
    ref->pos.x = 0;
    ref->pos.y = 0;
    node* leaf = NULL;

    int nb_levels = quad_tree_get_size(nb_boxes*nb_boxes);
    if (nb_levels > 0) {
        tree->box.size = ref->size;
        tree->parent = NULL;
        tree->box.pos.x = 0;
        tree->box.pos.y = 0;

        quad_tree_create(tree, nb_levels, nb_boxes * nb_boxes * nb_particules);
    } else {
        tree->nb_children = 0;
    }

    // Initialisation des Box
    for (int i = 0 ; i < nb_boxes ; i++) {
        row = i*nb_boxes;
        for (int j = 0 ; j < nb_boxes ; j++) {
            current = row + j;
            // Initialisation de la Box
            box_init(&(boxes[current]), nb_particules);
            boxes[current].size = box_size;
            boxes[current].pos.x = row;
            boxes[current].pos.y = j;

            if (tree->nb_children > 0) {
                leaf = quad_tree_find_leaf(tree, row, j);
            }

            // Génération des particules dans la Box
            box_generate_particules(&(boxes[current]), ref, &total, nb_particules, leaf);
        }
    }

    // Les feuilles ont été créées, il faut faire remonter les particules dans les boîtes parentes
    if (tree->nb_children > 0) {
        quad_tree_fill_boxes(tree, nb_levels);
    }
}

double compute_error(double ref, double val) {
    double error = val - ref;
    double percent = val == 0 ? 0 : 100.0*error/val;

    if (percent < 0)
        percent *= -1;

    /*if (error != 0.) {
        printf("Valeur attendue : %lf (%e), valeur calculée : %lf (%e)\n", ref, ref, val, val);
        printf("Le pourcentage d'erreur est de : %e%%\n", percent);
    }*/

    return percent;
}

void handle_error(double ref, double val, int* nb_errors, double* max_percent, double* r_ref, double* r_val) {
    double percent = compute_error(ref, val);
    if (percent != 0.) {
        (*nb_errors)++;
        if (percent > *max_percent) {
            (*max_percent) = percent;
            *r_ref = ref;
            *r_val = val;
        }
    }
}

/**
 * Vérifier que le calcul sur les boîtes donne bien des valeurs correspondant à celles obtenues dans la boîte de référence
 * @param ref la boîte de référence où toutes les planètes sont présentes
 * @param nb_total_planets le nombre total de planètes
 * @param boxes les boîtes de calcul
 */
void check_leafs(box* ref, int nb_total_planets, node** leafs) {
    printf("\n\tVérification de la cohérence des vecteurs de force - Barnes Hut\n");
    int current_box = -1;
    int nb_errors = 0;
    double max_percent = 0.;
    double r_ref = 0., r_val = 0.;
    int j = 0;
    int nb_particules = leafs[0]->box.nb_planets;

    for (int i = 0 ; i < nb_total_planets ; i++) {
        if (i % nb_particules == 0) {
            current_box++;
            j = 0;
        }

        handle_error(ref->force[i].x, leafs[current_box]->box.force[j].x, &nb_errors, &max_percent, &r_ref, &r_val);
        handle_error(ref->force[i].y, leafs[current_box]->box.force[j].y, &nb_errors, &max_percent, &r_ref, &r_val);

        j++;
    }
    printf("-- Résultat : %d erreurs trouvée(s) sur %d valeurs, avec un pourcentage d'erreur maximal de %e%% (attendu : %e, calculé : %e)\n", nb_errors, nb_total_planets*2, max_percent, r_ref, r_val);
}

/**
 * Vérifier que le calcul sur les boîtes donne bien des valeurs correspondant à celles obtenues dans la boîte de référence
 * @param ref la boîte de référence où toutes les planètes sont présentes
 * @param nb_total_planets le nombre total de planètes
 * @param boxes les boîtes de calcul
 */
void check_boxes(box* ref, int nb_total_planets, box* boxes) {
    printf("\n\tVérification de la cohérence des vecteurs de force - Séquentiel par boîtes\n");
    int current_box = -1;
    int nb_errors = 0;
    double max_percent = 0.;
    double r_ref = 0., r_val = 0.;
    int j = 0;
    int nb_particules = boxes[0].nb_planets;

    for (int i = 0 ; i < nb_total_planets ; i++) {
        if (i % nb_particules == 0) {
            current_box++;
            j = 0;
        }

        handle_error(ref->force[i].x, boxes[current_box].force[j].x, &nb_errors, &max_percent, &r_ref, &r_val);
        handle_error(ref->force[i].y, boxes[current_box].force[j].y, &nb_errors, &max_percent, &r_ref, &r_val);

        j++;
    }
    printf("-- Résultat : %d erreurs trouvée(s) sur %d valeurs, avec un pourcentage d'erreur maximal de %e%% (attendu : %e, calculé : %e)\n", nb_errors, nb_total_planets*2, max_percent, r_ref, r_val);
}

/**
 * Vérifier que les valeurs obtenues pour les forces correspond bien à celles obtenues dans la boîte de référence
 * @param ref le tableau des forces de référence
 * @param val le tableau des forces obtenues
 * @param nb_particles nombre de particules
 */
void check_forces(point* ref, point* val, int nb_particles) {
    printf("\n\tVérification de la cohérence des vecteurs de force - Parallèle\n");
    int nb_errors = 0;
    double max_percent = 0.;
    double r_ref = 0., r_val = 0.;

    for (int i = 0 ; i < nb_particles ; i++) {
        handle_error(ref[i].x, val[i].x, &nb_errors, &max_percent, &r_ref, &r_val);
        handle_error(ref[i].y, val[i].y, &nb_errors, &max_percent, &r_ref, &r_val);
    }
    printf("-- Résultat : %d erreurs trouvée(s) sur %d valeurs, avec un pourcentage d'erreur maximal de %e%% (attendu : %e, calculé : %e)\n", nb_errors, nb_particles*2, max_percent, r_ref, r_val);
}

void load_boxes(box* ref, box* boxes, node* tree, int nb_boxes, int nb_total_boxes, int nb_particules, int nb_planets, int world_size, int rendering) {
    if (rendering) {
        printf("\tLancement de la simulation par boîtes\n");
        printf("-- Nombre de boîtes : %d\n", nb_total_boxes);
        printf("-- Nombre de planètes par boîte : %d\n", nb_particules);
        printf("-- Nombre de planètes total : %d\n", nb_planets);
        printf("-- Taille du monde (en km) : %d\n", world_size);
        printf("-- Génération des boîtes\n");
    }

    box_init(ref, nb_planets);
    generate_boxes(ref, boxes, nb_boxes, world_size/nb_boxes, nb_particules, tree);
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

    // Chargement des boîtes
    box ref;
    box* boxes = malloc(sizeof(box) * nb_total_boxes);
    node tree;

    load_boxes(&ref, boxes, &tree, nb_boxes, nb_total_boxes, nb_particules, nb_planets, size, rendering);

    launch_sequential_simulation_box_on(&ref, boxes, &tree, nb_boxes, nb_particules, rendering);

    for (int i = 0 ; i < nb_total_boxes ; i++) {
        box_free(&boxes[i]);
    }
    free(boxes);
    box_free(&ref);
}

void get_leafs_and_compute(node** leafs, node* current, int *nb_leafs) {
    for (int i = 0 ; i < current->box.nb_planets ; i++) {
        current->box.force[i].x = 0;
        current->box.force[i].y = 0;
    }
    if (current->nb_children != 0) {
        for (int i = 0; i < current->nb_children ; i++) {
            get_leafs_and_compute(leafs, current->nodes[i], nb_leafs);
        }
    } else {
        leafs[*nb_leafs] = current;
        (*nb_leafs)++;
        calcul_force_own(&(current->box));
    }
}

/**
 * Lancement de la simulation sur plusieurs boîtes dont comparaison avec boîte de référence
 * @param ref Boîte de référence
 * @param boxes Boîtes de travail
 * @param tree Arbre pour Barnes-Hut
 * @param nb_boxes nombre de boîtes
 * @param nb_particules nombre de particules par boîte
 * @param rendering 1 = affichage des infos, 0 = mode silencieux
 */
void launch_sequential_simulation_box_on(box* ref, box* boxes, node* tree, int nb_boxes, int nb_particules, int rendering) {
    perf_t gen_start, gen_end, box_start, box_end, barnes_start, barnes_end;

    int nb_total_boxes = nb_boxes * nb_boxes;
    int nb_planets = nb_boxes * nb_boxes * nb_particules;

    node* leafs[nb_total_boxes];
    int nb_leafs = 0;

    if (rendering) {
        printf("\n\tLancement des calculs en séquentiel\n");
    }

    /** CALCUL DU BLOC TOTAL OU MODE STANDARD **/
    for (int i = 0; i < nb_planets; i++) {
        ref->force[i].x = 0;
        ref->force[i].y = 0;
    }

    perf(&gen_start);
    calcul_force_own(ref);
    perf(&gen_end);

    /** CALCUL DES DEUX BLOCS **/
    perf(&box_start);

    for (int i = 0; i < nb_total_boxes; i++) {
        for (int j = 0 ; j < nb_particules ; j++) {
            boxes[i].force[j].x = 0;
            boxes[i].force[j].y = 0;
        }

        // Calcul des forces du bloc
        calcul_force_own(&(boxes[i]));
    }

    for (int i = 0 ; i < nb_boxes ; i++) {
        for (int j = 0 ; j < nb_boxes ; j++) {
            if (i != j) {
                calcul_force_two_boxes(&(boxes[i]), &(boxes[j]), THRESHOLD);
            }
        }
    }

    perf(&box_end);

    /** CALCUL AVEC BARNES-HUT **/

    if (tree->nb_children != 0) {

        perf(&barnes_start);

        get_leafs_and_compute(leafs, tree, &nb_leafs);

        for (int i = 0; i < nb_total_boxes; i++) {
            calcul_force_Barnes_Hut(leafs[i], tree, THRESHOLD);
        }

        perf(&barnes_end);

        /// Vérification
        check_leafs(ref, nb_planets, leafs);
    }

    /** VÉRIFICATION DES DONNÉES **/
    check_boxes(ref, nb_planets, boxes);


    perf_diff(&gen_start, &gen_end);
    perf_diff(&box_start, &box_end);
    if (tree->nb_children != 0) {
        perf_diff(&barnes_start, &barnes_end);
    }

    if (rendering) {
        printf("Temps séquentiel - version basique : ");
        perf_printmicro(&gen_end);

        printf("Temps séquentiel - version box (sur %d boîtes) : ", nb_total_boxes);
        perf_printmicro(&box_end);

        if (tree->nb_children != 0) {
            printf("Temps séquentiel - version Barnes-Hut : ");
            perf_printmicro(&barnes_end);
        }
    }
}