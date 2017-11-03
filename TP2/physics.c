#include <math.h>
#include <limits.h>
#include "physics.h"
#define G 6.67e-11

void calcul_force_first_loop(planet* myplanets, planet* bufplanets,int size, point* forcebuf,double* dmin){
  int i,j;
  double sqdisty,sqdistx,norm_dist;
  for (i=0; i<size;i++){
    for(j=0;j<size;j++){ 
      if(i!=j){
	sqdistx = pow(abs(myplanets[i].pos.x-bufplanets[j].pos.x),2);
	sqdisty = pow(abs(myplanets[i].pos.y-bufplanets[j].pos.y),2);
	norm_dist = sqrt(sqdistx+sqdisty);
	forcebuf[i].x += (G*myplanets[i].mass*bufplanets[j].mass)/sqdistx;
	forcebuf[i].y += (G*myplanets[i].mass*bufplanets[j].mass)/sqdisty;
	if(dmin[i] > norm_dist){
	  dmin[i] = norm_dist;
	}
      }
    }
  }
}



void calcul_force(planet* myplanets, planet* bufplanets,int size, point* forcebuf,double* dmin){
  int i,j;
  double sqdisty,sqdistx,norm_dist;
  for (i=0; i<size;i++){
    for(j=0;j<size;j++){ 
      sqdistx = pow(abs(myplanets[i].pos.x-bufplanets[j].pos.x),2);
      sqdisty = pow(abs(myplanets[i].pos.y-bufplanets[j].pos.y),2);
      norm_dist = sqrt(sqdistx+sqdisty);
      forcebuf[i].x = (G*myplanets[i].mass*bufplanets[j].mass)/sqdistx;
      forcebuf[i].y = (G*myplanets[i].mass*bufplanets[j].mass)/sqdisty;
      if(dmin[i] > norm_dist){
	dmin[i] = norm_dist;
      }
    }
  }
}
      

double calcul_dtmin(planet* myplanets,double* dmin, int size){
  int i;
  double dtmin = MAX_DOUBLE;
  double sqspeed,norm_acc,delta;
  for(i=0;i<size;i++){
    sqspeed = pow(myplanets[i].speed.x,2)+pow(myplanets[i].speed.y,2);
    norm_acc = sqrt(pow(myplanets[i].acc.x,2)+pow(myplanets[i].acc.y,2));
    delta = sqspeed - 2*norm_acc*(-0.1*dmin[i]);
    dmin[i] = (-sqrt(sqspeed)+sqrt(delta))/norm_acc;
  }
  for(i=0;i<size;i++){
    if(dtmin > dmin[i]){
      dtmin = dmin[i];
    }
  }
  return dtmin;
}


void calcul_newpos(planet* myplanets,point* forcebuf,int size,double dt){
  int i,ax,ay;
  for (i=0;i<size;i++){
    ax = forcebuf[i].x/myplanets[i].mass;
    ay = forcebuf[i].y/myplanets[i].mass;
    myplanets[i].speed.x += ax*dt;
    myplanets[i].speed.y += ay*dt;
    myplanets[i].pos.x += myplanets[i].speed.x*dt + (ax*pow(dt,2))/2;
    myplanets[i].pos.y += myplanets[i].speed.y*dt + (ay*pow(dt,2))/2;
  }
}
