/*
 * sortmatch.c
 *
 *  Created on: 15 déc. 2008
 *      Author: coissac
 */

/*
 * sortword.c
 *
 *
 *  Created on: 6 nov. 2008
 *      Author: coissac
 */

#include "ecoprimer.h"
#include <math.h>

void su_smoothsort(void *base, uint32_t r, uint32_t N,
		   int (*less)(void *m, uint32_t a, uint32_t b),
		   void (*swap)(void *m, uint32_t a, uint32_t b));

static int less(void *m, uint32_t a, uint32_t b);
static void swap(void *m, uint32_t a, uint32_t b);


void sortmatch(pprimermatch_t table,uint32_t N)
{
	su_smoothsort((void*)table,0,N,less,swap);
}

int less(void *m, uint32_t a, uint32_t b)
{
	pprimermatch_t t;

	t = (pprimermatch_t)m;

	return t[a].position <= t[b].position;
}

void swap(void *m, uint32_t a, uint32_t b)
{
  primermatch_t  tmp;
  pprimermatch_t t;

  t = (pprimermatch_t)m;
  tmp = t[a];
  t[a]= t[b];
  t[b]= tmp;
}

