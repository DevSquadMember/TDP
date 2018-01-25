#include <stdio.h>
#include "utils-mpi.h"

void init_group(struct group* group) {
    MPI_Comm_rank(group->comm, &group->rank);
    MPI_Comm_size(group->comm, &group->size);
}