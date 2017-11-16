#ifndef PHYSICS
#define PHYSICS

#define MAX_DOUBLE 100000

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

struct planet_handle {
    double m, px, py;
};

void calcul_force_seq(planet*, planet*,int, point*,double*);
void calcul_force_first_loop(planet*, struct planet_handle*,int, point*,double*);
void calcul_force(planet*, struct planet_handle*, int, point*,double*);
double calcul_dtmin(planet*,point*,double*,int);
void calcul_newpos(planet*,point*,int,double);

#endif
