#ifndef TDP_SAVER_H
#define TDP_SAVER_H

#include "physics.h"

void render_seq(int nb_total_planets, char* title);
void render(int nb_procs, int nb_total_planets, char* title);

void save_seq(planet* planets, int nb_planets);
void save(int proc_id, planet* planets, int nb_planets);

void save_close();

#endif //TDP_SAVER_H
