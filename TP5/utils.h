#ifndef TDP_UTILS_H
#define TDP_UTILS_H

#include "physics.h"

void box_init(box* box, int nb_planets);
void box_copy(box* src, box* dest);
void box_free(box* box);

void copy_planet(planet* src, planet* dest);

int quad_tree_get_size(int nb_boxes);
void quad_tree_create(node* tree, int levels, int nb_particules);
void quad_tree_free(node* tree);
struct node* quad_tree_find_leaf(node* tree, int x, int y);
void quad_tree_fill_boxes(node* tree, int nb_levels);

#endif //TDP_UTILS_H
