#include "physics.h"
#include "parser.h"
#include "saver.h"
#include <stdio.h>
#include <stdlib.h>

#define TEST_FILE "test_planets.txt"
#define NB_ITERATIONS 1

int main(int argc,char ** argv){
    int nb_it = NB_ITERATIONS;
    int i,j;
    double dtmin;

    char* filename = TEST_FILE;

    if (argc > 1) {
        nb_it = atoi(argv[1]);
        if (argc > 2) {
            filename = argv[2];
        }
    }
    int nb_planets = parser_nb_planets(filename);

    planet myplanets[nb_planets];
    parser_load(myplanets, 0, nb_planets);

    point* forcebuf = malloc(nb_planets * sizeof(point));
    double* dmin = malloc(nb_planets * sizeof(double));
    for(i=0;i<nb_planets;i++){
        dmin[i] = MAX_DOUBLE;
    }
    save_seq(myplanets, nb_planets);
    for(i=0;i<nb_it;i++){
        for(j=0;j<nb_planets;j++){
            forcebuf[j].x = 0;
            forcebuf[j].y = 0;
            dmin[j] = MAX_DOUBLE;
        }
        calcul_force_seq(myplanets, myplanets, nb_planets, forcebuf, dmin);
        dtmin = calcul_dtmin(myplanets, forcebuf, dmin, nb_planets);
        //printf("dtmin it %d : %f\n",i,dtmin);
        calcul_newpos(myplanets, forcebuf, nb_planets, dtmin);
        ///printf("fin it %d : 0 : posx : %f, posy : %f, 1 : posx : %f ,posy : %f, 2 : posx : %f, posy : %f\n", i,myplanets[0].pos.x,myplanets[0].pos.y,myplanets[1].pos.x,myplanets[1].pos.y,myplanets[2].pos.x,myplanets[2].pos.y);
        save_seq(myplanets, nb_planets);
    }
    save_close();
    render_seq(nb_planets, "Calcul séquentiel");

    free(forcebuf);
    free(dmin);
    return 0;
}
