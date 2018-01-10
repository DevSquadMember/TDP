#include "simulation.h"
#include <stdlib.h>

#define TEST_FILE "test_planets.txt"

int main(int argc,char ** argv){
    int rendering = 1;

    char* filename = TEST_FILE;

    if (argc > 1) {
        filename = argv[1];
        if (argc > 2) {
            rendering = atoi(argv[2]);
        }
    }

    //launch_sequential_simulation_blocs(rendering, filename);

    // nb_boxes, nb_particles_per_box, world_size, rendering (1 = true)
    launch_sequential_simulation_box(2, 5, 1000000000, 1);

    return EXIT_SUCCESS;
}
