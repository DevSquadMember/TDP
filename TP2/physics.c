#include <math.h>
#include <limits.h>
#include <stdio.h>
#include "physics.h"
#define G 1

void calcul_force_first_loop(planet* myplanets, planet* bufplanets,int size, point* forcebuf,double* dmin){
  int i,j,sidex,sidey;
  double distx,disty,angle,force,sqdist,dist;
  for(i=0; i<size;i++){
    for(j=0;j<size;j++){ 
      if(i!=j){
	
	sidex = sidey = 1;
	


	distx = bufplanets[j].pos.x-myplanets[i].pos.x;
	disty = bufplanets[j].pos.y-myplanets[i].pos.y;
        
	if(distx < 0){
	  distx = -distx;
	  sidex = -1;
	}
	if(disty < 0){
	  disty = -disty;
	  sidey = -1;
	}

	angle = atan(disty/distx);

	sqdist = pow(distx,2)+pow(disty,2);
	dist = sqrt(sqdist);
	force = (G*myplanets[i].mass*bufplanets[j].mass)/sqdist;


	
	forcebuf[i].x += sidex*force*cos(angle);
	forcebuf[i].y += sidey*force*sin(angle);

	if(dmin[i] > dist){
	  dmin[i] = dist;
	}
       
      }
    }
  }
}



void calcul_force(planet* myplanets, planet* bufplanets,int size, point* forcebuf,double* dmin){
  int i,j,sidex,sidey;
  double distx,disty,angle,force,sqdist,dist;
  for (i=0; i<size;i++){
    for(j=0;j<size;j++){ 
      	
      sidex = sidey = 1;
	

      distx = bufplanets[j].pos.x-myplanets[i].pos.x;
      disty = bufplanets[j].pos.y-myplanets[i].pos.y;
        
      if(distx < 0){
	distx = -distx;
	sidex = -1;
      }
      if(disty < 0){
	disty = -disty;
	sidey = -1;
      }

      angle = atan(disty/distx);

      sqdist = pow(distx,2)+pow(disty,2);
      dist = sqrt(sqdist);
      force = (G*myplanets[i].mass*bufplanets[j].mass)/sqdist;    
	
      forcebuf[i].x += sidex*force*cos(angle);
      forcebuf[i].y += sidey*force*sin(angle);

      if(dmin[i] > dist){
	dmin[i] = dist;
      }

    }
  }
}
      

double calcul_dtmin(planet* myplanets,point* forcebuf,double* dmin, int size){
  int i;
  double dtmin = MAX_DOUBLE;
  double sqspeed,norm_acc,delta;

  for(i=0;i<size;i++){
    sqspeed = pow(myplanets[i].speed.x,2)+pow(myplanets[i].speed.y,2);
    norm_acc = sqrt(pow(forcebuf[i].x/myplanets[i].mass,2)+pow(forcebuf[i].y/myplanets[i].mass,2));
    delta = sqspeed - 2*norm_acc*(-0.1*dmin[i]);
    //printf("planete %d, sqspeed :  %f, acc : %f, delta: %f\n",i,sqspeed,norm_acc,delta);
    dmin[i] = (-sqrt(sqspeed)+sqrt(delta))/norm_acc;
  }
  for(i=0;i<size;i++){
    if(dtmin > dmin[i]){
      dtmin = dmin[i];
    }
  }
  return 1;
}


void calcul_newpos(planet* myplanets,point* forcebuf,int size,double dt){
  int i;
  double ax,ay;
  for (i=0;i<size;i++){
    ax = forcebuf[i].x/myplanets[i].mass;
    ay = forcebuf[i].y/myplanets[i].mass;
    myplanets[i].speed.x += ax*dt;
    myplanets[i].speed.y += ay*dt;
    myplanets[i].pos.x += myplanets[i].speed.x*dt + (ax*pow(dt,2))/2;
    myplanets[i].pos.y += myplanets[i].speed.y*dt + (ay*pow(dt,2))/2;
  }
}
