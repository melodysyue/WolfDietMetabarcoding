/*
 * ahocorasick.h
 *
 *  Created on: 26 march 2011
 *      Author: tiayyba
 */

#ifndef H_ahocorasick
#define H_ahocorasick

#include "ecoprimer.h"

typedef struct aho_output_t{
	uint32_t wordidx; //index of strict word (dont save the word of 64B)
	bool_t isdirect; //we need to find both direct and reverse words so we must know which one is it
}aho_output;

typedef struct aho_output_count_t{
	uint32_t count;
	aho_output *out_set;
}aho_output_count;

typedef struct aho_state_t{
	int32_t id;
	struct aho_state_t *next[4]; //for labels A=0,C=1,G=2 and T=3
	struct aho_state_t *fail;
	aho_output_count output;
}aho_state;

typedef struct queue_node_t {
	aho_state *state_node;
	struct queue_node_t *next;
}queue_node;

typedef struct{
	queue_node *first;
	queue_node *last;
}aho_queue;

pprimercount_t ahoc_lookforStrictPrimers (pecodnadb_t database, uint32_t seqdbsize,uint32_t exampleCount,
		                     pwordcount_t words,poptions_t options);
#endif /* H_ahocorasick */

