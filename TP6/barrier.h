#ifndef TDP_BARRIER_H
#define TDP_BARRIER_H

#include <pthread.h>

typedef struct barrier {
    int nb_stopped;
    int nb_threads;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} barrier;

void barrier_init(struct barrier* barrier, int nb_threads);
void barrier_free(struct barrier* barrier);
void barrier_STOP(int nb_threads);
void barrier_stop(struct barrier* barrier);

#endif //TDP_BARRIER_H
