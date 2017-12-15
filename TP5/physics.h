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

typedef struct box{
  planet* planet_list;
  point* forcebuf;
  int size;
  point box_pos;
  planet center;
} box;

struct planet_handle {
    double m, px, py;
};

void copy_planet(planet* src, planet* dest);

void calcul_force_seq(planet*, int, point*,double*);
void calcul_force_first_loop(planet*, struct planet_handle*,int, point*,double*);
void calcul_force_own(planet* myplanets, int size, point* forcebuf, double* dmin);
void calcul_force(planet*, struct planet_handle*, int, point*,double*);
void calcul_force_other(planet* myplanets, planet* planets, int size, point* forcebuf, double* dmin);
void calcul_force_both(planet* myplanets, planet* planets, int size, point* forcebuf1, point* forcebuf2, double* dmin1, double* dmin2);
double calcul_dtmin(planet*,point*,double*,int);
void calcul_newpos(planet*,point*,int,double);

#endif
