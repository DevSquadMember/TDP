#include "simulation.h"
#include "perf.h"
#include <stdlib.h>
#include <mpi.h>
#include <printf.h>

#define TEST_FILE "test_planets.txt"
#define NB_PARTICLES 5
#define WORLD_SIZE 1000000000

#define THRESHOLD 1000000

MPI_Datatype point_type;
MPI_Datatype points_type;
MPI_Datatype planet_type;
MPI_Datatype planets_type;
MPI_Datatype box_type;

void declare_types(int nb_planets) {
    struct point p;
    struct planet pl;
    struct box b;
    MPI_Aint i1, i2;

    /// Déclaration du type POINT

    MPI_Datatype type[2] = { MPI_DOUBLE, MPI_DOUBLE };
    int block_len[2] = { 1, 1};
    MPI_Aint disp[2];

    MPI_Get_address(&p, &i1);
    MPI_Get_address(&(p.x), &i2);
    disp[0] = i2 - i1;
    MPI_Get_address(&(p.y), &i2);
    disp[1] = i2 - i1;

    MPI_Type_create_struct(2, block_len, disp, type, &point_type);
    MPI_Type_commit(&point_type);

    /// Déclaration du type PLANET

    MPI_Datatype type2[4] = { MPI_DOUBLE, point_type, point_type, point_type };
    int block_len2[4] = { 1, 1, 1, 1};
    MPI_Aint disp2[4];

    MPI_Get_address(&pl, &i1);
    MPI_Get_address(&(pl.mass), &i2);
    disp2[0] = i2 - i1;
    MPI_Get_address(&(pl.acc), &i2);
    disp2[1] = i2 - i1;
    MPI_Get_address(&(pl.speed), &i2);
    disp2[2] = i2 - i1;
    MPI_Get_address(&(pl.pos), &i2);
    disp2[3] = i2 - i1;

    MPI_Type_create_struct(4, block_len2, disp2, type2, &planet_type);
    MPI_Type_commit(&planet_type);

    /// Déclaration du type PLANETS_TYPE

    MPI_Type_vector(nb_planets, 1, 1, planet_type, &planets_type);
    MPI_Type_commit(&planets_type);

    /// Déclaration du type POINTS_TYPE

    MPI_Type_vector(nb_planets, 1, 1, point_type, &points_type);
    MPI_Type_commit(&points_type);

    /// Déclaration du type BOX

    MPI_Datatype type3[3] = { MPI_INT, point_type, planet_type };
    int block_len3[3] = { 1, 1, 1 };
    MPI_Aint disp3[3];

    MPI_Get_address(&b, &i1);
    MPI_Get_address(&(b.size), &i2);
    disp3[0] = i2 - i1;
    MPI_Get_address(&(b.pos), &i2);
    disp3[1] = i2 - i1;
    MPI_Get_address(&(b.center), &i2);
    disp3[2] = i2 - i1;

    MPI_Type_create_struct(3, block_len3, disp3, type3, &box_type);
    MPI_Type_create_resized(box_type, 0, sizeof(b), &box_type);
    MPI_Type_commit(&box_type);
}

void print_boxes(int size, struct box* boxes) {
    for (int i = 0 ; i < size ; i++) {
        printf("BOX - %d\n", i);
        printf("NB_PLANETS - %d\n", boxes[i].nb_planets);
        printf("SIZE - %d\n", boxes[i].size);
        printf("POS - %lf, %lf\n", boxes[i].pos.x, boxes[i].pos.y);

        for (int j = 0 ; j < boxes[i].nb_planets ; j++) {
            printf("Planet %d - Mass (%lf), Pos(%lf, %lf)\n", j, boxes[i].planets[j].mass, boxes[i].planets[j].pos.x, boxes[i].planets[j].pos.y);
        }
    }
}

void print_local_box(int rank, int size, struct box local_box) {
    for (int i = 0 ; i < size ; i++) {
        if (i == rank) {
            printf("BOX - %d\n", i);
            printf("NB_PLANETS - %d\n", local_box.nb_planets);
            printf("SIZE - %d\n", local_box.size);
            printf("POS - %lf, %lf\n", local_box.pos.x, local_box.pos.y);

            for (int j = 0 ; j < local_box.nb_planets ; j++) {
                printf("Planet %d - Mass (%lf), Pos(%lf, %lf)\n", j, local_box.planets[j].mass, local_box.planets[j].pos.x, local_box.planets[j].pos.y);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc,char ** argv) {
    int rank, size;
    int nb_particles = NB_PARTICLES;
    int world_size = WORLD_SIZE;

    int rendering = 1;

    //char* filename = TEST_FILE;

    if (argc > 1) {
        nb_particles = atoi(argv[1]);
        if (argc > 2) {
            rendering = atoi(argv[2]);
        }
    }

    perf_t scatter_start, scatter_end, par_start, par_end;

    MPI_Status status_p_send, status_p_recv, status_b_send, status_b_recv;
    MPI_Request request_p_send, request_p_recv, request_b_send, request_b_recv;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    declare_types(nb_particles);

    struct box ref_box; // boîte de référence
    struct box* boxes;

    struct box local_box; // boîte locale de chaque processeur
    box_init(&local_box, nb_particles);

    struct box buffer1_box; // boîte de travail 1
    box_init(&buffer1_box, nb_particles);
    struct box buffer2_box; // boîte de travail 2
    box_init(&buffer2_box, nb_particles);

    struct box *sender, *receiver;

    // Chargement de la boîte de référence et des boîtes pour les processeurs
    if (rank == 0) {
        int nb_total_boxes = size; // size * size;
        int nb_planets = size * nb_particles; // size * size * nb_particles;
        boxes = malloc(sizeof(struct box) * nb_total_boxes);
        load_boxes(&ref_box, boxes, size, nb_total_boxes, nb_particles, nb_planets, world_size, rendering);
        /*box_init(&ref_box, nb_planets);
        generate_boxes(&ref_box, boxes, size, world_size/size, nb_particles);*/
    }

    perf(&scatter_start);

    // Envoi des boîtes sur chaque processeur
    MPI_Scatter(boxes, 1, box_type, &local_box, 1, box_type, 0, MPI_COMM_WORLD);

    // Envoi des planètes dans les boîtes sur chaque processeur
    MPI_Scatter(ref_box.planets, nb_particles, planet_type, local_box.planets, nb_particles, planet_type, 0, MPI_COMM_WORLD);

    perf(&scatter_end);

    /** CALCUL EN SÉQUENTIEL **/
    if (rank == 0) {
        launch_sequential_simulation_box_on(&ref_box, boxes, size, nb_particles, rendering);
    }
    /** FIN DU CALCUL EN SÉQUENTIEL **/

    /** CALCUL EN PARALLÈLE **/
    // Position des voisins dans l'anneau
    int ring_left = (rank - 1 + size) % size;
    int ring_right = (rank + 1) % size;

    // Initialisation des pointeurs de buffer envoi/réception
    sender = &buffer1_box;
    receiver = &buffer2_box;

    // Copie des données des planètes locales dans le buffer d'envoi
    for (int i = 0; i < nb_particles; i++) {
        sender->planets[i].mass = local_box.planets[i].mass;
        sender->planets[i].acc = local_box.planets[i].acc;
        sender->planets[i].speed = local_box.planets[i].speed;
        sender->planets[i].pos = local_box.planets[i].pos;

        local_box.force[i].x = 0;
        local_box.force[i].y = 0;
    }

    perf(&par_start);

    // Parcours de l'anneau pour une itération
    for (int i = 0; i < size; i++) {
        // Lancement des communications - boîtes et particules
        MPI_Isend(sender, 1, box_type, ring_right, ring_right, MPI_COMM_WORLD, &request_b_send);
        MPI_Isend(sender->planets, nb_particles, planet_type, ring_right, ring_right, MPI_COMM_WORLD, &request_p_send);
        MPI_Irecv(receiver, 1, box_type, ring_left, rank, MPI_COMM_WORLD, &request_b_recv);
        MPI_Irecv(receiver->planets, nb_particles, planet_type, ring_left, rank, MPI_COMM_WORLD, &request_p_recv);

        if (i == 0) {
            calcul_force_own(&local_box);
        } else {
            calcul_force_two_boxes(&local_box, sender, THRESHOLD);
        }

        // Attente des communications
        MPI_Wait(&request_b_send, &status_b_send);
        MPI_Wait(&request_b_recv, &status_b_recv);
        MPI_Wait(&request_p_send, &status_p_send);
        MPI_Wait(&request_p_recv, &status_p_recv);

        // Echange des buffers des anneaux
        if (receiver == &buffer1_box) {
            sender = &buffer1_box;
            receiver = &buffer2_box;
        } else {
            sender = &buffer2_box;
            receiver = &buffer1_box;
        }
    }

    perf(&par_end);

    /** FIN DU CALCUL EN PARALLÈLE **/

    if (rank == 0) {
        perf_diff(&scatter_start, &scatter_end);
        perf_diff(&par_start, &par_end);

        if (rendering) {
            printf("Temps parallèle - scatter : ");
            perf_printmicro(&scatter_end);

            printf("Temps parallèle - calcul : ");
            perf_printmicro(&par_end);
        }

        check_boxes(&ref_box, size*nb_particles, boxes);
    }

    MPI_Type_free(&point_type);
    MPI_Type_free(&points_type);
    MPI_Type_free(&planet_type);
    MPI_Type_free(&planets_type);
    MPI_Type_free(&box_type);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
