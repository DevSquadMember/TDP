#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include "timer.h"

struct timer {
    struct timeval t_start;
    struct timeval t_end;
    double time;
};

/*void timer_begin() {
    gettimeofday(&t_start, NULL);
}

double timer_end() {
    gettimeofday(&t_end, NULL);
    return ((t_end.tv_usec + t_end.tv_sec*1e6) - (t_start.tv_usec+t_start.tv_sec*1e6));
}*/

struct timer* timer_create() {
    struct timer* t = malloc(sizeof(struct timer));
    t->time = 0;
    return t;
}

void timer_start(struct timer* timer) {
    gettimeofday(&(timer->t_start), NULL);
}

void timer_stop(struct timer* timer) {
    gettimeofday(&(timer->t_end), NULL);
    timer->time += ((timer->t_end.tv_usec + timer->t_end.tv_sec*1e6) - (timer->t_start.tv_usec+timer->t_start.tv_sec*1e6));
}

double timer_get_time(struct timer* timer) {
    return timer->time;
}

void timer_free(struct timer* timer) {
    free(timer);
}
