#ifndef TDP_SIMULATION_H
#define TDP_SIMULATION_H

#include "physics.h"

void generate_particules(int nb_particules, char* filename);

void generate_boxes(box* ref, box* boxes, int nb_boxes, int box_size, int nb_particules);

void check_boxes(box* ref, int nb_total_planets, box* boxes);

void load_boxes(box* ref, box* boxes, int nb_boxes, int nb_total_boxes, int nb_particules, int nb_planets, int world_size, int rendering);

void launch_sequential_simulation_box(int nb_boxes, int nb_particules, int size, int rendering);

void launch_sequential_simulation_box_on(box* ref, box* boxes, int nb_boxes, int nb_particules, int rendering);

void launch_sequential_simulation_blocs(int rendering, char* filename);

#endif //TDP_SIMULATION_H
