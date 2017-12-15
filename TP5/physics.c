#include <math.h>
#include <limits.h>
#include <stdio.h>
#include "physics.h"
//#define G 1
#define G 6.67e-11

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

void calcul_force_own(planet* myplanets, int size, point* forcebuf, double* dmin){
    int i,j,sidex,sidey;
    double distx,disty,angle,force,sqdist,dist;

    for(i=0; i<size;i++){
        for(j=0;j<size;j++){
            if(i!=j){

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
                force = (G*myplanets[i].mass*myplanets[j].mass)/sqdist;

                forcebuf[i].x += sidex*force*cos(angle);
                forcebuf[i].y += sidey*force*sin(angle);

                if(dmin[i] > dist){
                    dmin[i] = dist;
                }
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

void calcul_force_both(planet* myplanets, planet* planets, int size, point* forcebuf1, point* forcebuf2, double* dmin1, double* dmin2) {
    int i, j, sidex, sidey;
    double distx, disty, angle, force, sqdist, dist;
    double valx, valy;
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

            valx = sidex*force*cos(angle);
            valy = sidey*force*sin(angle);

            forcebuf1[i].x += valx;
            forcebuf1[i].y += valy;

            forcebuf2[j].x -= valx;
            forcebuf2[j].y -= valy;

            if(dmin1[i] > dist){
                dmin1[i] = dist;
            }
            if(dmin2[j] > dist){
                dmin2[j] = dist;
            }
        }
    }
}

/** Calcul de dt min **/
double calcul_dtmin(planet* myplanets, point* forcebuf, double* dmin, int size){
    int i;
    double dtmin = MAX_DOUBLE;
    double sqspeed,norm_acc,delta;

    for(i = 0 ; i < size ; i++) {
        sqspeed = pow(myplanets[i].speed.x, 2) + pow(myplanets[i].speed.y, 2);
        norm_acc = sqrt(pow(forcebuf[i].x / myplanets[i].mass, 2) +  pow(forcebuf[i].y / myplanets[i].mass, 2));
        delta = sqspeed - 2 * norm_acc * (-0.1 * dmin[i]);
        dmin[i] = (-sqrt(sqspeed) + sqrt(delta))/norm_acc;

        if(dtmin > dmin[i]){
            dtmin = dmin[i];
        }
    }

    return dtmin;
}

/** Calcul des nouvelles positions pour les planètes **/
void calcul_newpos(planet* myplanets, point* forcebuf, int size, double dt) {
    int i;
    double ax, ay;
    for (i = 0 ; i < size ; i++) {
        ax = 1.0*forcebuf[i].x / myplanets[i].mass;
        ay = 1.0*forcebuf[i].y / myplanets[i].mass;

        myplanets[i].pos.x += myplanets[i].speed.x * dt + (ax*pow(dt,2))/2;
        myplanets[i].pos.y += myplanets[i].speed.y * dt + (ay*pow(dt,2))/2;

        myplanets[i].speed.x += ax*dt;
        myplanets[i].speed.y += ay*dt;
    }
}


/** Calcul du centre de masse **/
void calcul_center_mass(planet * myplanets, int size, planet* center){
  center->mass = 0;
  center->pos.x = 0;
  center->pos.y = 0;
  int i;
  for (i = 0 ; i < size ; i++){
    center->mass += myplanets[i].mass;
    center->pos.x += myplanets[i].pos.x*myplanets[i].mass;
    center->pos.y += myplanets[i].pos.y*myplanets[i].mass;
  }
  center->pos.x /= center->mass;
  center->pos.y /= center->mass;
}