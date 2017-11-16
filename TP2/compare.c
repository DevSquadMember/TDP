#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "simulation.h"

#define NB_PROC_MAX 4
#define NB_PARTICULES 1200
#define NB_MIN_ITERATIONS 5
#define NB_MAX_ITERATIONS 10
#define COEF 1.25

#define PARTICULE_FILENAME "particles.txt"
#define FILENAME "res.dot"

void measure_parallel_time(FILE *file, struct timeval* t_start, struct timeval* t_end, int nb_iterations, int nb_procs) {
    char command[100];
    sprintf(command, "make NP=%d run_par %d %s 0", nb_procs, nb_iterations, PARTICULE_FILENAME);
    char answer[10];
    FILE* f;

    gettimeofday(t_start, NULL);
    f = popen(command, "w");
    // exécution de la commande
    fgets(answer, 10, f);
    pclose(f);
    //launch_parrallel_simulation(nb_iterations, 0, PARTICULE_FILENAME);

    gettimeofday(t_end, NULL);
    fprintf(file, "%lf ", ((t_end->tv_usec + t_end->tv_sec*1e6) - (t_start->tv_usec+t_start->tv_sec*1e6)));
}

void measure_sequential_time(FILE *file, struct timeval* t_start, struct timeval* t_end, int nb_iterations) {
    gettimeofday(t_start, NULL);

    launch_sequential_simulation(nb_iterations, 0, PARTICULE_FILENAME);

    gettimeofday(t_end, NULL);
    fprintf(file, "%lf ", ((t_end->tv_usec + t_end->tv_sec*1e6) - (t_start->tv_usec+t_start->tv_sec*1e6)));
}

int main(int argc, char** argv) {
    FILE *file = fopen(FILENAME, "w");
    int *array;
    struct timeval t_start, t_end;

    generate_particules(NB_PARTICULES, PARTICULE_FILENAME);

    printf("Comparaison version séquentielle / parallèle jusqu'à %d processeurs\n", NB_PROC_MAX);
    printf("-- Nombre de particules : %d\n", NB_PARTICULES);
    printf("-- Nombre d'itérations min : %d\n", NB_MIN_ITERATIONS);
    printf("-- Nombre d'itérations max : %d\n", NB_MAX_ITERATIONS);

    int size = NB_MIN_ITERATIONS;
    while (size < NB_MAX_ITERATIONS) {
        printf("-- %d itérations en cours\n", size);
        fprintf(file, "%d ", size);
        measure_sequential_time(file, &t_start, &t_end, size);
        for (int i = 2 ; i <= NB_PROC_MAX ; i++) {
            measure_parallel_time(file, &t_start, &t_end, size, i);
        }
        fprintf(file, "\n");

        size *= COEF;
    }

    FILE *gnuplot = popen("gnuplot -persistent", "w");

    fprintf(gnuplot, "set title \"Comparaison version sequentielle/parallele pour %d particules\"\n", NB_PARTICULES);
    fprintf(gnuplot, "set xlabel \"%s\"\n", "Nombre d'iterations");
    fprintf(gnuplot, "set ylabel \"%s\"\n", "Temps en microsecondes");
    fprintf(gnuplot, "plot \"%s\" using 1:2 title \"Version sequentielle\" with lines", FILENAME);
    for (int i = 2 ; i <= NB_PROC_MAX ; i++)
        fprintf(gnuplot, ", \"%s\" using 1:%d title \"Version parallele (nb procs=%d)\" with lines", FILENAME, i+1, i);
    fprintf(gnuplot, "\n");

    fclose(file);
    fclose(gnuplot);
}

