#ifndef TDP_TIMER_H
#define TDP_TIMER_H

struct timer;

struct timer* timer_create();
void timer_start(struct timer* timer);
void timer_stop(struct timer* timer);
double timer_get_time(struct timer* timer);
void timer_free(struct timer* timer);

#endif //TDP_TIMER_H
