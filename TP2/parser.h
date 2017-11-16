#ifndef TDP_PARSER_H
#define TDP_PARSER_H

#include "physics.h"

int parser_nb_planets(char* filename);
int parser_load(planet* planets, int start, int nb_planets);

#endif //TDP_PARSER_H
