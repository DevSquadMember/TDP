int BS;
int barrier;
int barrier_id;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
pthread_mutex_t mut_barrier = PTHREAD_MUTEX_INITIALIZER;

void barrier_STOP(int nb_threads){
    pthread_mutex_lock(&mut_barrier);
    barrier++;

    if (barrier == nb_threads) {
	barrier = 0;
	pthread_cond_broadcast(&cond);
    }
    else{
	pthread_cond_wait(&cond, &mut_barrier);
    }
    pthread_mutex_unlock(&mut_barrier);
}
