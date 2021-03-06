#ifndef PHYSICS
#define PHYSICS

#define MAX_DOUBLE 10000000000000000000000000.0

typedef struct point {
    double x;
    double y;
} point;

typedef struct planet{
    double mass;
    point acc;
    point speed;
    point pos;
} planet;

typedef struct box {
    planet* planets;
    int nb_planets;
    point* force;
    int size;
    point pos;
    planet center;
} box;

typedef struct node {
    struct box box;
    struct node** nodes;
    int nb_children;
    struct node* parent;
} node;

struct planet_handle {
    double m, px, py;
};

void calcul_force_seq(planet*, int, point*,double*);
void calcul_force_first_loop(planet*, struct planet_handle*,int, point*,double*);
void calcul_force_own(box* box);
void calcul_force(planet*, struct planet_handle*, int, point*,double*);
void calcul_force_other(planet* myplanets, planet* planets, int size, point* forcebuf, double* dmin);
void calcul_force_Barnes_Hut(struct node* node, struct node* remote_node, double threshold);
void calcul_force_two_boxes(box* my_box, box* remote_box,double threshold);
void calcul_force_center(planet* center,planet* myplanets,int size, point* forcebuf);
void calcul_force_complete(planet* myplanets, int sizea, planet* planets, int sizeb, point* forcebuf1);
void calcul_center_mass(box* box);

#endif
