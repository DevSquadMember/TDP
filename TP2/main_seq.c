#include "physics.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>


int main(int argc,char ** argv){
  int nb_it = atoi(argv[1]);
  int i,j;
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
  myplanets[1].speed.y = 1;
  
  myplanets[2].mass = 100;
  myplanets[2].pos.x = 0;
  myplanets[2].pos.y = 100;
  myplanets[2].speed.x = 0;
  myplanets[2].speed.y = 3;
  
  printf("Bonjour\n");

  point* forcebuf = malloc(3*sizeof(point));
  double* dmin = malloc(3*sizeof(double));
  for(i=0;i<3;i++){
    dmin[i] = MAX_DOUBLE;
  }
  for(i=0;i<nb_it;i++){
    for(j=0;j<3;j++){
      forcebuf[j].x = 0;
      forcebuf[j].y = 0;
    }
    calcul_force_first_loop(myplanets,myplanets,3,forcebuf,dmin);
    dtmin = calcul_dtmin(myplanets,forcebuf,dmin,3);
    //printf("dtmin it %d : %f\n",i,dtmin);
    calcul_newpos(myplanets,forcebuf, 3, dtmin);
    printf("fin it %d : 0 : posx : %f, posy : %f, 1 : posx : %f ,posy : %f, 2 : posx : %f, posy : %f\n", i,myplanets[0].pos.x,myplanets[0].pos.y,myplanets[1].pos.x,myplanets[1].pos.y,myplanets[2].pos.x,myplanets[2].pos.y);
  }


  free(forcebuf);
  free(dmin);
  return 0;
}
