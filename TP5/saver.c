#include <stdio.h>
#include <unistd.h>
#include "saver.h"

#define FINAL_DATA_FILENAME "res.dat"
#define FINAL_BLOCS_FILENAME "res-bloc.dat"
#define DATA_FILENAME "res.dat%d"

FILE* res_file = NULL;
FILE* res_file1 = NULL;
FILE* res_file2 = NULL;
int counter = 0;
int counter1 = 0;
int counter2 = 0;

void render_plot(int nb_total_planets, char*title, char* filename, char* resname) {
    char args[20] = "";
    char final_title[100];
    sprintf(final_title, "%s en %d itérations", title, counter-1);

    // Ouverture du shell et lancement de gnuplot
    FILE* f = popen("gnuplot", "w");

    // exécution de la commande gnuplot
    fprintf(f, "set term png\n");
    fprintf(f, "set title \"%s\"\n", final_title);
    fprintf(f, "set output \"%s\"\n", resname);
    fprintf(f, "plot ");
    for (int i = 0 ; i < nb_total_planets ; i++) {
        sprintf(args, "%d:%d", 2*(i+1), 2*(i+1) + 1);
        fprintf(f, "\"%s\" using %s title \"Planète %d\"with lines", filename, args, (i+1));
        if (i < nb_total_planets - 1) {
            fprintf(f, ", ");
        }
    }
    fprintf(f, "\n");
    fflush(f);
    // terminer l'envoi de commandes et fermer gnuplot
    sleep(1);
    pclose(f);

    printf("Graphe généré et sauvegardé dans le fichier %s\n", resname);
}

void render_seq(int nb_total_planets, char* title) {
    render_plot(nb_total_planets, title, FINAL_DATA_FILENAME, "res.png");
}

// join file1 file2 | join - file3 > output
void render(int nb_procs, int nb_total_planets, char* title, char* res_filename) {
    if (nb_procs > 1) {
        char command[200] = "join";
        char filename[10];
        for (int i = 0; i < nb_procs; i++) {
            sprintf(filename, DATA_FILENAME, i);
            printf("Fichier intermédiaire %s généré\n", filename);
            if (i < 2) {
                sprintf(command, "%s %s", command, filename);
            } else {
                sprintf(command, "%s | join - %s", command, filename);
            }
        }
        sprintf(command, "%s > %s\n", command, res_filename);

        FILE *f = popen(command, "r");
        sleep(1);
        pclose(f);
    }

    render_plot(nb_total_planets, title, res_filename, "res-bloc.png");
}

void save_seq(planet* planets, int nb_planets) {
    if (res_file == NULL) {
        res_file = fopen(FINAL_DATA_FILENAME, "w");
    }
    fprintf(res_file, "%d ", counter);
    for (int i = 0 ; i < nb_planets ; i++) {
        fprintf(res_file, "%lf %lf ", planets[i].pos.x, planets[i].pos.y);
    }
    fprintf(res_file, "\n");
    counter++;
}

void save(int proc_id, planet* planets, int nb_planets) {
    FILE** file = proc_id == 1 ? &res_file1 : &res_file2;
    int* c = proc_id == 1 ? &counter1 : &counter2;
    if (*file == NULL) {
        char filename[10];
        sprintf(filename, DATA_FILENAME, proc_id);
        *file = fopen(filename, "w");
    }
    fprintf(*file, "%d ", *c);
    for (int i = 0 ; i < nb_planets ; i++) {
        fprintf(*file, "%lf %lf ", planets[i].pos.x, planets[i].pos.y);
    }
    fprintf(*file, "\n");
    *c = *c+1;
}

void save_close() {
    if (res_file != NULL) {
        fclose(res_file);
    }
    if (res_file1 != NULL) {
        fclose(res_file1);
    }
    if (res_file2 != NULL) {
        fclose(res_file2);
    }
}
