#include "simulation.h"
#include "perf.h"
#include <stdlib.h>
#include <mpi.h>
#include <printf.h>

#define TEST_FILE "test_planets.txt"
#define NB_PARTICLES 5
#define WORLD_SIZE 1000000000

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

    perf_t scatter_start, scatter_end;

    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    declare_types(nb_particles);

    struct box ref_box; // boîte de référence
    struct box* boxes;

    struct box local_box; // boîte locale de chaque processeur
    box_init(&local_box, nb_particles);

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

    // nb_boxes, nb_particles_per_box, world_size, rendering (1 = true)
    ///launch_sequential_simulation_box(2, 5, 1000000000, 1);

    // Lancement du calcul en séquentiel
    if (rank == 0) {
        launch_sequential_simulation_box_on(&ref_box, boxes, size, nb_particles, rendering);
    }

    if (rank == 0) {
        perf_diff(&scatter_start, &scatter_end);
        printf("Temps de scatter : ");
        perf_printmicro(&scatter_end);
    }

    MPI_Type_free(&point_type);
    MPI_Type_free(&points_type);
    MPI_Type_free(&planet_type);
    MPI_Type_free(&planets_type);
    MPI_Type_free(&box_type);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
