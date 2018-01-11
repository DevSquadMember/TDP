#ifndef TDP_SIMULATION_H
#define TDP_SIMULATION_H

#include "physics.h"

void generate_particules(int nb_particules, char* filename);

void generate_boxes(box* ref, box* boxes, int nb_boxes, int box_size, int nb_particules, node* tree);

void check_boxes(box* ref, int nb_total_planets, box* boxes);

void check_forces(point* ref, point* val, int nb_particles);

void load_boxes(box* ref, box* boxes, node* tree, int nb_boxes, int nb_total_boxes, int nb_particules, int nb_planets, int world_size, int rendering);

void launch_sequential_simulation_box(int nb_boxes, int nb_particules, int size, int rendering);

void launch_sequential_simulation_box_on(box* ref, box* boxes, node* tree, int nb_boxes, int nb_particules, int rendering);

#endif //TDP_SIMULATION_H
