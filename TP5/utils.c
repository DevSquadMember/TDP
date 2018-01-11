#include <stdlib.h>
#include "utils.h"

void box_init(box* box, int nb_planets) {
    box->pos.x = 0.;
    box->pos.y = 0.;
    box->size = 0;
    box->force = malloc(sizeof(point) * nb_planets);
    box->nb_planets = nb_planets;
    box->planets = malloc(sizeof(planet) * nb_planets);
}

void box_copy(box* src, box* dest) {
    box_init(dest, src->nb_planets);
    for (int i = 0 ; i < src->nb_planets ; i++) {
        copy_planet(&(src->planets[i]), &(dest->planets[i]));
    }
}

void box_free(box* box) {
    free(box->force);
    free(box->planets);
}

void copy_planet(planet* src, planet* dest) {
    dest->acc.x = src->acc.x;
    dest->acc.y = src->acc.y;
    dest->mass = src->mass;
    dest->pos.x = src->pos.x;
    dest->pos.y = src->pos.y;
    dest->speed.x = src->speed.x;
    dest->speed.y = src->speed.y;
}

int quad_tree_get_size(int nb_boxes) {
    int size = 0;
    if (nb_boxes == 0)
        return size;
    while (nb_boxes != 1) {
        if(nb_boxes%4 != 0)
            return 0;
        nb_boxes = nb_boxes/4;
        size++;
    }
    return size;
}

void quad_tree_create(node* tree, int levels, int nb_particules) {
    box_init(&(tree->box), nb_particules);
    if (levels > 0) {
        tree->nb_children = 4;
        tree->nodes = malloc(sizeof(node*) * 4);
        int half = tree->box.size/2;
        for (int i = 0 ; i < 4 ; i++) {
            node* n = malloc(sizeof(node));
            n->parent = tree;
            n->box.pos.x = tree->box.pos.x + (i/2) * half;
            n->box.pos.y = tree->box.pos.y + (i%2) * half;
            n->box.size = half;
            tree->nodes[i] = n;
            quad_tree_create(n, levels - 1, nb_particules/4);
        }
    } else {
        tree->nb_children = 0;
    }
}

void quad_tree_free(node* tree) {
    for (int i = 0 ; i < tree->nb_children ; i++) {
        quad_tree_free(tree->nodes[i]);
    }
    box_free(&(tree->box));
    free(tree->nodes);
}

struct node* quad_tree_find_leaf(node* tree, int x, int y) {
    struct node* current = tree;
    while (current->nb_children != 0) {
        int i = 0;
        if (x > current->box.pos.x) {
            i += 2;
        }
        if (y > current->box.pos.y) {
            i += 1;
        }
        current = current->nodes[i];
    }
    return current;
}

void copy_leafs(node* tree) {
    int current = 0;
    for (int i = 0 ; i < tree->nb_children ; i++) {
        for (int j = 0 ; j < tree->nodes[i]->nb_children ; j++) {
            copy_planet(&(tree->nodes[i]->box.planets[j]), &(tree->box.planets[current]));
            current++;
        }
    }
    calcul_center_mass(&(tree->box));
    if (tree->parent != NULL) {
        copy_leafs(tree->parent);
    }
}

void quad_tree_fill_boxes(node* tree, int nb_levels) {
    if (nb_levels > 1) {
        for (int j = 0; j < tree->nb_children; j++) {
            quad_tree_fill_boxes(tree->nodes[j], nb_levels - 1);
        }
    } else if (nb_levels == 1) {
        copy_leafs(tree);
    }
}
