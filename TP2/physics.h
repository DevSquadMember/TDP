#ifndef PHYSICS
#define PHYSICS
#define MAX_DOUBLE 100000

typedef struct point{
  double x;
  double y;
} point;

typedef struct planet{
  double mass;
  point acc;
  point speed;
  point pos;
} planet;

void calcul_force_first_loop(planet*, planet*,int, point*,double*);
void calcul_force(planet*, planet*, int, point*,double*);
double calcul_dtmin(planet*,point*,double*,int);
void calcul_newpos(planet*,point*,int,double);

#endif
