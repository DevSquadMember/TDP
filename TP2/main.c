#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "parser.h"
#include "saver.h"

#define FILENAME "planets.txt"
#define NB_ITERATIONS 1

#define TAG_BUFFER_SEND 1
#define TAG_BUFFER_RECV 2

#define TRUE 1
#define FALSE 0

struct planet_handle {
    double m, px, py;
};

void declare_types(MPI_Datatype* planet_type, MPI_Datatype* planets_type, struct planet_handle* handle, int nb_planets) {
    // Déclaration du type PLANET_TYPE selon la structure planet_handle
    MPI_Aint i1, i2;
    MPI_Datatype type[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
    int block_len[3] = { 1, 1, 1};
    MPI_Aint disp[3];

    MPI_Get_address(handle, &i1);
    MPI_Get_address(&(handle->m), &i2);
    disp[0] = i2 - i1;
    MPI_Get_address(&(handle->px), &i2);
    disp[1] = i2 - i1;
    MPI_Get_address(&(handle->py), &i2);
    disp[2] = i2 - i1;

    MPI_Type_create_struct(3, block_len, disp, type, planet_type);
    MPI_Type_commit(planet_type);

    // Déclaration du type PLANETS_TYPE
    MPI_Type_vector(nb_planets, 1, 1, *planet_type, planets_type);
    MPI_Type_commit(planets_type);
}

void ring_iteration() {
    
}

int main( int argc, char **argv ) {
    int i, j;
    int ring_id, ring_size;
    int iteration, nb_iterations = NB_ITERATIONS;

    // Données pour MPI
    MPI_Datatype planet_type, planets_type;
    MPI_Status status_send, status_recv;
    MPI_Request request_send, request_recv;

    // Initialisation de MPI
    MPI_Init( NULL, NULL );
    MPI_Comm_rank(MPI_COMM_WORLD, &ring_id);
    MPI_Comm_size(MPI_COMM_WORLD, &ring_size);

    // Chargement du fichier de données
    char* filename = FILENAME;
    if (argc > 1) {
        nb_iterations = atoi(argv[1]);
        if (argc > 2) {
            filename = argv[2];
        }
    }
    int nb_total_planets = parser_nb_planets(filename);
    int nb_planets = nb_total_planets/ring_size;

    if (nb_total_planets > 0) {
        // Déclaration du groupe de planètes locales
        planet planets[nb_planets];

        // Structures buffer pour l'envoi et la réception des groupes de planètes MPI
        struct planet_handle handle1[nb_planets], handle2[nb_planets];
        struct planet_handle *sender, *receiver;
        
        // Chargement des données des planètes locales
        parser_load(planets, ring_id, nb_planets);

        // Déclaration des types pour MPI
        declare_types(&planet_type, &planets_type, &(handle1[0]), nb_planets);

        // Copie des données des planètes locales dans le buffer d'envoi
        for (i = 0; i < nb_planets; i++) {
            handle1[i].m = planets[i].mass;
            handle1[i].px = planets[i].pos.x;
            handle1[i].py = planets[i].pos.y;
        }

        //printf("My own - Rank %d\n", ring_id);
        //printf("Planet %lf (%lf, %lf - %lf, %lf)\n", planets[ring_id].mass, planets[ring_id].position.x, planets[ring_id].position.y, planets[ring_id].vitesse.x, planets[ring_id].vitesse.y);

        // Initialisation des pointeurs de buffer envoi/réception
        sender = handle1;
        receiver = handle2;

        if (ring_size > 1) {
            save(ring_id, planets, nb_planets);
        } else {
            save_seq(planets, nb_planets);
        }

        for (iteration = 0; iteration < nb_iterations; iteration++) {
            // Parcours de l'anneau pour une itération
            for (i = 0; i < ring_size; i++) {
                // Lancement des communications
                MPI_Isend(sender, 1, planets_type, (ring_id + 1) % ring_size, TAG_BUFFER_SEND, MPI_COMM_WORLD,
                          &request_send);
                MPI_Irecv(receiver, 1, planets_type, (ring_id - 1 + ring_size) % ring_size, TAG_BUFFER_SEND,
                          MPI_COMM_WORLD, &request_recv);
                //printf("Process %d sent to %d and received from %d\n", ring_id, (ring_id + 1)%ring_size, (ring_id - 1 + ring_size)%ring_size);

                // Partie calcul
                /// TODO

                if (ring_id == 0) {
                    printf("-- Run %d for process %d\n", i, ring_id);
                    for (j = 0; j < nb_planets; j++) {
                        printf("Process %d sent Planet %d (%lf (%lf, %lf))\n", ring_id, j, sender[j].m, sender[j].px,
                               sender[j].py);
                    }
                }

                // Attente des communications
                MPI_Wait(&request_send, &status_send);
                MPI_Wait(&request_recv, &status_recv);

                // Echange des buffers des anneaux
                sender = receiver;
                if (receiver == handle1) {
                    receiver = handle2;
                } else {
                    receiver = handle1;
                }
            }

            //// TEST
            for (i = 0 ; i < nb_planets ; i++) {
                planets[i].pos.x += planets[i].pos.x;
                planets[i].pos.y += planets[i].pos.y;
            }

            if (ring_size > 1) {
                save(ring_id, planets, nb_planets);
            } else {
                save_seq(planets, nb_planets);
            }
        }
        save_close();
        MPI_Bcast(&i, 1, MPI_INT, 1, MPI_COMM_WORLD);

        if (ring_id == 0) {
            char title[100];
            sprintf(title, "Calcul parallèle sur %d noeuds\n", ring_size);
            render(ring_size, nb_total_planets, title);
        }

        MPI_Type_free(&planets_type);
        MPI_Type_free(&planet_type);
    }

    MPI_Finalize();
}