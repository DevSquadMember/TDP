#include "simulation.h"
#include <stdlib.h>

#define TEST_FILE "test_planets.txt"
#define NB_ITERATIONS 1

int main(int argc,char ** argv){
    int nb_iterations = NB_ITERATIONS;
    int rendering = 1;

    char* filename = TEST_FILE;

    if (argc > 1) {
        nb_iterations = atoi(argv[1]);
        if (argc > 2) {
            filename = argv[2];
            if (argc > 3) {
                rendering = atoi(argv[3]);
            }
        }
    }

    launch_sequential_simulation_blocs(nb_iterations, rendering, filename);

    return EXIT_SUCCESS;
}
