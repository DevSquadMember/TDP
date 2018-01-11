#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "physics.h"
#define G 6.67e-11

void box_init(box* box, int nb_planets) {
    box->pos.x = 0.;
    box->pos.y = 0.;
    box->size = 0;
    box->force = malloc(sizeof(point) * nb_planets);
    box->nb_planets = nb_planets;
    box->planets = malloc(sizeof(planet) * nb_planets);
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

/** Calcul des forces en séquentiel sur ses propres planètes **/
void calcul_force_seq(planet* myplanets, int size, point* forcebuf, double* dmin) {
    int i,j,sidex,sidey;
    double distx,disty,angle,force,sqdist,dist;
    for(i=0; i<size;i++) {
        for(j=0;j<size;j++) {
            if (i != j) {

                sidex = sidey = 1;

                distx = myplanets[j].pos.x - myplanets[i].pos.x;
                disty = myplanets[j].pos.y - myplanets[i].pos.y;

                if(distx < 0){
                    distx = -distx;
                    sidex = -1;
                }
                if(disty < 0){
                    disty = -disty;
                    sidey = -1;
                }

                angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

                sqdist = pow(distx,2)+pow(disty,2);
                dist = sqrt(sqdist);
                force = (G*myplanets[i].mass * myplanets[j].mass)/sqdist;

                forcebuf[i].x += sidex*force*cos(angle);
                forcebuf[i].y += sidey*force*sin(angle);

                if(dmin[i] > dist){
                    dmin[i] = dist;
                }
            }
        }
    }
}

/** Calcul des forces lors du démarrage du parcours de l'anneau : on calcule les forces sur le propre groupe de planètes **/
void calcul_force_first_loop(planet* myplanets, struct planet_handle* bufplanets,int size, point* forcebuf,double* dmin){
    int i,j,sidex,sidey;
    double distx,disty,angle,force,sqdist,dist;

    for(i=0; i<size;i++){
        for(j=0;j<size;j++){
            if(i!=j){

                sidex = sidey = 1;

                distx = bufplanets[j].px - myplanets[i].pos.x;
                disty = bufplanets[j].py - myplanets[i].pos.y;

                if(distx < 0){
                    distx = -distx;
                    sidex = -1;
                }
                if(disty < 0){
                    disty = -disty;
                    sidey = -1;
                }

                angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

                sqdist = pow(distx,2)+pow(disty,2);
                dist = sqrt(sqdist);
                force = (G*myplanets[i].mass*bufplanets[j].m)/sqdist;

                forcebuf[i].x += sidex*force*cos(angle);
                forcebuf[i].y += sidey*force*sin(angle);

                if(dmin[i] > dist){
                    dmin[i] = dist;
                }
            }
        }
    }
}

void calcul_force_own(box* box) {
    int i, j, sidex, sidey;
    double distx, disty, angle, force, sqdist;

    for(i=0; i<box->nb_planets;i++){
        for(j=0;j<box->nb_planets;j++){
            if(i!=j){

                sidex = sidey = 1;

                distx = box->planets[j].pos.x - box->planets[i].pos.x;
                disty = box->planets[j].pos.y - box->planets[i].pos.y;

                if(distx < 0){
                    distx = -distx;
                    sidex = -1;
                }
                if(disty < 0){
                    disty = -disty;
                    sidey = -1;
                }

                angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

                sqdist = pow(distx,2)+pow(disty,2);
                force = (G*box->planets[i].mass*box->planets[j].mass)/sqdist;

                box->force[i].x += sidex*force*cos(angle);
                box->force[i].y += sidey*force*sin(angle);
            }
        }
    }
}


/** Calcul des forces entre deux groupes différents **/
void calcul_force(planet* myplanets, struct planet_handle* bufplanets,int size, point* forcebuf,double* dmin){
    int i,j,sidex,sidey;
    double distx,disty,angle,force,sqdist,dist;
    for (i=0; i<size;i++){
        for(j=0;j<size;j++){

            sidex = sidey = 1;

            distx = bufplanets[j].px-myplanets[i].pos.x;
            disty = bufplanets[j].py-myplanets[i].pos.y;

            if(distx < 0){
                distx = -distx;
                sidex = -1;
            }
            if(disty < 0){
                disty = -disty;
                sidey = -1;
            }

            angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

            sqdist = pow(distx, 2) + pow(disty, 2);
            dist = sqrt(sqdist);
            force = (G * myplanets[i].mass * bufplanets[j].m) / sqdist;

            forcebuf[i].x += sidex*force*cos(angle);
            forcebuf[i].y += sidey*force*sin(angle);

            if(dmin[i] > dist){
                dmin[i] = dist;
            }
        }
    }
}

void calcul_force_other(planet* myplanets, planet* planets, int size, point* forcebuf, double* dmin){
    int i,j,sidex,sidey;
    double distx,disty,angle,force,sqdist,dist;

    for (i=0; i<size;i++){
        for(j=0;j<size;j++){

            sidex = sidey = 1;

            distx = planets[j].pos.x-myplanets[i].pos.x;
            disty = planets[j].pos.y-myplanets[i].pos.y;

            if(distx < 0){
                distx = -distx;
                sidex = -1;
            }
            if(disty < 0){
                disty = -disty;
                sidey = -1;
            }

            angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

            sqdist = pow(distx, 2) + pow(disty, 2);
            dist = sqrt(sqdist);
            force = (G * myplanets[i].mass * planets[j].mass) / sqdist;

            forcebuf[i].x += sidex*force*cos(angle);
            forcebuf[i].y += sidey*force*sin(angle);

            if(dmin[i] > dist){
                dmin[i] = dist;
            }
        }
    }
}

void calcul_force_Barnes_Hut(box* my_box, box* remote_box, double threshold){
  int i;
  if (remote_box->isLeaf /* booleen qui dit si on est une feuille */ == 1){
    calcul_force_complete(my_box->planets, my_box->nb_planets, remote_box->planets, remote_box->nb_planets, my_box->force);
  } else {
    double dist;
    for(i = 0 ; i < 4 ; i++){
      dist = sqrt(pow(my_box->center.pos.x - /* pos x du centre de enfant i */, 2) + pow(my_box->center.pos.y - /* pos y du centre de enfant i */, 2));
      if (remote_box->size/(2*dist) > threshold){
	calcul_force_center(&(/* centre de l enfant */), my_box->planets, my_box->nb_planets, my_box->force);
      } else {
	calcul_force_Barnes_Hut(my_box,/* enfant */, threshold);
      }
    }
  }
}

void calcul_force_two_boxes(box* my_box, box* remote_box,double threshold){
    double dist = sqrt(pow(my_box->center.pos.x - remote_box->center.pos.x, 2) + pow(my_box->center.pos.y - remote_box->center.pos.y, 2));
    if (remote_box->size/dist > threshold){
        calcul_force_center(&(remote_box->center), my_box->planets, my_box->nb_planets, my_box->force);
    } else {
        calcul_force_complete(my_box->planets, my_box->nb_planets, remote_box->planets, remote_box->nb_planets, my_box->force);
    }
}

void calcul_force_center(planet* center, planet* myplanets, int size, point* forcebuf){
    int i, sidex, sidey;
    double distx, disty, angle, force, dist, sqdist;
    for( i = 0 ; i < size ; i++ ){
        sidex = sidey = 1;
        distx = center->pos.x-myplanets[i].pos.x;
        disty = center->pos.y-myplanets[i].pos.y;

        if(distx < 0){
            distx = -distx;
            sidex = -1;
        }
        if(disty < 0){
            disty = -disty;
            sidey = -1;
        }
        angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

        sqdist = pow(distx, 2) + pow(disty, 2);
        dist = sqrt(sqdist);
        force = (G * myplanets[i].mass * center->mass) / sqdist;

        forcebuf[i].x = sidex*force*cos(angle);
        forcebuf[i].y = sidey*force*sin(angle);
    }
}

void calcul_force_complete(planet* myplanets, int sizea, planet* planets,int sizeb, point* forcebuf1) {
    int i, j, sidex, sidey;
    double distx, disty, angle, force, sqdist, dist;
    double valx, valy;
    for (i=0; i<sizea;i++){
        for(j=0;j<sizeb;j++){

            sidex = sidey = 1;

            distx = planets[j].pos.x-myplanets[i].pos.x;
            disty = planets[j].pos.y-myplanets[i].pos.y;

            if(distx < 0){
                distx = -distx;
                sidex = -1;
            }
            if(disty < 0){
                disty = -disty;
                sidey = -1;
            }

            angle = (distx == 0) ? M_PI/2 : atan(disty/distx);

            sqdist = pow(distx, 2) + pow(disty, 2);
            dist = sqrt(sqdist);
            force = (G * myplanets[i].mass * planets[j].mass) / sqdist;

            valx = sidex*force*cos(angle);
            valy = sidey*force*sin(angle);

            forcebuf1[i].x += valx;
            forcebuf1[i].y += valy;
        }
    }
}

/** Calcul du centre de masse **/
void calcul_center_mass(box* box){
    box->center.mass = 0;
    box->center.pos.x = 0;
    box->center.pos.y = 0;
    int i;
    for (i = 0 ; i < box->nb_planets ; i++){
        box->center.mass += box->planets[i].mass;
        box->center.pos.x += box->planets[i].pos.x*box->planets[i].mass;
        box->center.pos.y += box->planets[i].pos.y*box->planets[i].mass;
    }
    box->center.pos.x /= box->center.mass;
    box->center.pos.y /= box->center.mass;
}
