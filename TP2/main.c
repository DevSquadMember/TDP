#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include "simulation.h"
#include "parser.h"
#include "saver.h"

#define FILENAME "planets.txt"
#define NB_ITERATIONS 1
#define MAX_SIZE_FOR_RENDERING 12

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

void launch_parrallel_simulation(int nb_iterations, int rendering, char* filename) {
    int i, j;
    int ring_id, ring_size;
    int iteration;
    double dtmin;
    int generate_graph;

    // Données pour MPI
    MPI_Datatype planet_type, planets_type;
    MPI_Status status_send, status_recv;
    MPI_Request request_send, request_recv;

    // Initialisation de MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &ring_id);
    MPI_Comm_size(MPI_COMM_WORLD, &ring_size);

    int nb_total_planets = parser_nb_planets(filename);
    int nb_planets = nb_total_planets/ring_size;

    if (ring_id == 0 && rendering) {
        printf("Lancement du calcul en parallèle sur %d processeurs\n", ring_size);
        printf("-- Nombre d'itérations : %d\n", nb_iterations);
        printf("-- Fichier des planètes : %s\n", filename);
        printf("-- Nombre de planètes trouvées : %d\n", nb_total_planets);
        printf("-- Dispersion des planètes par groupe de : %d\n", nb_planets);
        printf("-- Début des calculs\n");
    }

    generate_graph = nb_planets <= MAX_SIZE_FOR_RENDERING && rendering;

    if (nb_total_planets > 0) {
        // Déclaration du groupe de planètes locales
        planet planets[nb_planets];

        // Structures buffer pour l'envoi et la réception des groupes de planètes MPI
        struct planet_handle handle1[nb_planets], handle2[nb_planets];
        struct planet_handle *sender, *receiver;

        // Chargement des données des planètes locales
        parser_load(planets, ring_id*nb_planets, nb_planets, rendering != 0);

        // Déclaration des types pour MPI
        declare_types(&planet_type, &planets_type, &(handle1[0]), nb_planets);

        // Initialisation des pointeurs de buffer envoi/réception
        sender = handle1;
        receiver = handle2;

        // Sauvegarde dans le fichier de résultats de la position initiale des planètes
        if (generate_graph) {
            if (ring_size > 1) {
                save(ring_id, planets, nb_planets);
            } else {
                save_seq(planets, nb_planets);
            }
        }

        point forcebuf[nb_planets];
        double dmin[nb_planets];
        for (i = 0 ; i < nb_planets ; i++) {
            dmin[i] = MAX_DOUBLE;
        }

        // Position des voisins dans l'anneau
        int ring_left = (ring_id - 1 + ring_size) % ring_size;
        int ring_right = (ring_id + 1) % ring_size;

        for (iteration = 0; iteration < nb_iterations; iteration++) {
            // Copie des données des planètes locales dans le buffer d'envoi
            for (i = 0; i < nb_planets; i++) {
                sender[i].m = planets[i].mass;
                sender[i].px = planets[i].pos.x;
                sender[i].py = planets[i].pos.y;
            }
            for (j = 0 ; j < nb_planets ; j++) {
                forcebuf[j].x = 0;
                forcebuf[j].y = 0;
                dmin[j] = MAX_DOUBLE;
            }

            // Parcours de l'anneau pour une itération
            for (i = 0; i < ring_size; i++) {
                // Lancement des communications
                MPI_Isend(sender, 1, planets_type, ring_right, ring_right, MPI_COMM_WORLD, &request_send);
                MPI_Irecv(receiver, 1, planets_type, ring_left, ring_id, MPI_COMM_WORLD, &request_recv);

                /** Partie calcul des forces **/

                if (i == 0) { // on applique le calcul sur le groupe des planètes locales
                    calcul_force_first_loop(planets, sender, nb_planets, forcebuf, dmin);
                } else { // on applique le calcul sur un autre groupe
                    calcul_force(planets, sender, nb_planets, forcebuf, dmin);
                }

                /** Fin de la partie de calcul des forces **/

                // Attente des communications
                MPI_Wait(&request_send, &status_send);
                MPI_Wait(&request_recv, &status_recv);

                // Echange des buffers des anneaux
                if (receiver == handle1) {
                    sender = handle1;
                    receiver = handle2;
                } else {
                    sender = handle2;
                    receiver = handle1;
                }
            }
            // Calcul du dtmin local
            double local_dtmin = calcul_dtmin(planets, forcebuf, dmin, nb_planets);

            // Récupérer la plus petite valeur de dtmin
            MPI_Allreduce(&local_dtmin, &dtmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

            // Calcul de la nouvelle position
            calcul_newpos(planets, forcebuf, nb_planets, dtmin);

            // Ecriture dans le fichier de la position des planètes
            if (generate_graph) {
                if (ring_size > 1) {
                    save(ring_id, planets, nb_planets);
                } else {
                    save_seq(planets, nb_planets);
                }
            }
        }
        // Fermeture du fichier
        if (generate_graph)
            save_close();

        // On attend que tout le monde ait fermé le fichier des résultats
        MPI_Barrier(MPI_COMM_WORLD);

        // Le processus 0 va générer le graphe de résultat
        if (ring_id == 0 && generate_graph) {
            char title[100];
            sprintf(title, "Calcul parallèle sur %d noeuds", ring_size);
            render(ring_size, nb_total_planets, title);
        }

        MPI_Type_free(&planets_type);
        MPI_Type_free(&planet_type);
    }

    MPI_Finalize();
}


int main( int argc, char **argv ) {
    int nb_iterations = NB_ITERATIONS;
    int rendering = 1;

    // Chargement du fichier de données
    char* filename = FILENAME;
    if (argc > 1) {
        nb_iterations = atoi(argv[1]);
        if (argc > 2) {
            filename = argv[2];
            if (argc > 3) {
                rendering = atoi(argv[3]);
            }
        }
    }

    launch_parrallel_simulation(nb_iterations, rendering, filename);

    return EXIT_SUCCESS;
}