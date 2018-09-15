#ifndef PRIMERSETS_H_
#define PRIMERSETS_H_

#include "ecoprimer.h"

#define PRIMERS_IN_SET_COUNT 10

typedef struct {
	int		*set_wellIdentifiedTaxa;
	int32_t set_pairs[PRIMERS_IN_SET_COUNT];
	float set_specificity;
	float set_coverage;
	float set_lmean;
	float set_lcov;
	float set_score;
	int32_t set_intaxa;
	int32_t set_wi_cnt;
}pairset;

typedef struct{
	ppair_t* sortedpairs;
	int32_t sorted_count;
	pecodnadb_t seqdb;
	poptions_t options;
}SetParams;

typedef struct{
	float t_spc; //specificity contribution
	float t_cov; //coverage contribution
	float t_lmd; //link spread difference
	float len; //length
	float score; //score
}primerscore;

void add_pair_in_set (pairset *pair_set, int32_t pset_idx, int32_t prb_idx, SetParams *pparams);
void get_next_pair_options (int *pair_wi_count_sorted_ids, pairset *pair_set, SetParams *pparams);
float get_links_distribution (int prb_idx, pairset *prob_set, SetParams *pparams);
pairset build_primers_set_greedy_spc (SetParams *pparams);
void get_set_mean_cov_stats (pairset *prob_set, SetParams *pparams);
void some_other_set_possibilities (pairset *pair_set,
		ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options);
void sets_by_SimulatedAnealing (pairset *pair_set,
		ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options);
void sets_by_TabuSearch (pairset *pair_set,
		ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options);
pairset * sets_by_BruteForce (ppair_t * sortedpairs, 
		int32_t sorted_count, pecodnadb_t seqdb, poptions_t options);
pairset * extend_set_randomly (pairset *pair_set, SetParams *params, int extend_to_cnt);
void build_and_print_sets (ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options);
int32_t get_next_option_increasing_cov (pairset *pair_set, SetParams *pparams);
void reset_set_props (pairset *pair_set, SetParams *pparams);
void primers_graph_graphviz (ppair_t * sortedpairs,
		int32_t sorted_count, poptions_t options);
size_t primers_changeSortedArray (ppair_t ** pairs,
		size_t sorted_count, poptions_t options);
size_t primers_filterWithGivenLinks (ppair_t ** pairs,
		size_t sorted_count, poptions_t options);
#endif 
