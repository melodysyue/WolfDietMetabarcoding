#ifndef THERMOSTATS_H_
#define THERMOSTATS_H_

#include "../libecoprimer/ecoprimer.h"

void getThermoProperties (ppair_t* pairs, size_t count, poptions_t options);
word_t extractSite(char* sequence, size_t begin, size_t length, bool_t strand);

#endif