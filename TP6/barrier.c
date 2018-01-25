#include "barrier.h"

int BS;
int nb_stopped;
int barrier_id;

pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
pthread_mutex_t mut_barrier = PTHREAD_MUTEX_INITIALIZER;

void barrier_init(struct barrier* barrier, int nb_threads) {
    barrier->nb_threads = nb_threads;
    barrier->nb_stopped = 0;
    pthread_cond_init(&barrier->cond, NULL);
    pthread_mutex_init(&barrier->mutex, NULL);
}

void barrier_free(struct barrier* barrier) {
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
}

void barrier_STOP(int nb_threads) {
    pthread_mutex_lock(&mut_barrier);
    nb_stopped++;

    if (nb_stopped == nb_threads) {
        nb_stopped = 0;
        pthread_cond_broadcast(&cond);
    }
    else{
        pthread_cond_wait(&cond, &mut_barrier);
    }
    pthread_mutex_unlock(&mut_barrier);
}

void barrier_stop(struct barrier* barrier) {
    pthread_mutex_lock(&barrier->mutex);
    barrier->nb_stopped++;

    if (barrier->nb_stopped == barrier->nb_threads) {
        barrier->nb_stopped = 0;
        pthread_cond_broadcast(&barrier->cond);
    }
    else{
        pthread_cond_wait(&barrier->cond, &barrier->mutex);
    }
    pthread_mutex_unlock(&barrier->mutex);
}
