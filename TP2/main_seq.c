#include "physics.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>


int main(int argc,char ** argv){
  int nb_it = atoi(argv[1]);
  int i;
  double dtmin;
  planet myplanets[3];

  myplanets[0].mass = 9000;
  myplanets[0].pos.x = 0;
  myplanets[0].pos.y = 0;
  myplanets[0].speed.x = 0;
  myplanets[0].speed.y = 0;
  
  myplanets[1].mass = 100;
  myplanets[1].pos.x = 100;
  myplanets[1].pos.y = 0;
  myplanets[1].speed.x = 0;
  myplanets[1].speed.y = 100;
  
  myplanets[2].mass = 100;
  myplanets[2].pos.x = 0;
  myplanets[2].pos.y = 100;
  myplanets[2].speed.x = 0;
  myplanets[2].speed.y = 100;
  
  point* forcebuf = malloc(3*sizeof(point));
  double* dmin = malloc(3*sizeof(double));
  for(i=0;i<3;i++){
    dmin[i] = MAX_DOUBLE;
  }
  for(i=0;i<nb_it;i++){
    for(i=0;i<3;i++){
      forcebuf[i].x = 0;
      forcebuf[i].y = 0;
    }
    calcul_force_first_loop(myplanets,myplanets,3,forcebuf,dmin);
    dtmin = calcul_dt_min(myplanets, dmin,3);
    calcul_newpos(myplanets,forcebuf, 3, dtmin);
  }


  free(forcebuf);
  free(dmin);
  return 0;
}
