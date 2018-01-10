#ifndef TDP_SIMULATION_H
#define TDP_SIMULATION_H

void generate_particules(int nb_particules, char* filename);

void launch_sequential_simulation_box(int nb_boxes, int nb_particules, int size, int rendering);

void launch_sequential_simulation_blocs(int rendering, char* filename);

#endif //TDP_SIMULATION_H
