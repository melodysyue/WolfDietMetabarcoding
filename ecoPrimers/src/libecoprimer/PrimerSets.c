#include <stdint.h>
#include <stdlib.h>
#include <time.h>

//#include "ecoprimer.h"
#include "PrimerSets.h"

int TabuList[PRIMERS_IN_SET_COUNT];
int next_tabu_slot = -1;
int total_pairs = -1;
//int32_t total_wi = -1;

int32_t counttaxon(int32_t taxid);
int find_in_tabu (int index);

int32_t count_taxons (int32_t taxid)
{
	static int32_t count = 0;
	static int32_t slots = 0;
	static int32_t *taxon_array = NULL;
	int32_t i;

	if (taxid == -1)
	{
		if (taxon_array)
			ECOFREE (taxon_array, "Could not free memory for taxon array");
		taxon_array = NULL;
		slots = 0;
		count = 0;
	}
	else
	{
		for (i = 0; i < count; i++)
		{
			if (taxid == taxon_array[i]) return count;
		}
		
		if (count == slots)
		{
			slots += 500;
			
			if (taxon_array == NULL)
			{
				taxon_array = (int32_t *) ECOMALLOC(slots*sizeof (int32_t),
				            "Could not allocate memory for taxon array");
			}
			else
				taxon_array = (int32_t *) ECOREALLOC(taxon_array, slots*sizeof (int32_t),
								            "Could not reallocate memory for taxon array");
		}
		taxon_array[count] = taxid;
		count++;
	}
	return count;
}

float get_set_coverage (pairset *p_set, SetParams *pparams, int32_t pidx_toexclude)
{
	int32_t i, j;
	float cov;
	int32_t s_intaxa = 0;
	int32_t seqcount;
	
	//counttaxon(-1);
	count_taxons (-1);
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (p_set->set_pairs[i] == -1 || pidx_toexclude == i) continue;
		
		seqcount=pparams->sortedpairs[p_set->set_pairs[i]]->pcr.ampcount;
		for (j=0; j < seqcount; j++)
			if (pparams->sortedpairs[p_set->set_pairs[i]]->pcr.amplifias[j].sequence->isexample
					&& pparams->sortedpairs[p_set->set_pairs[i]]->pcr.amplifias[j].sequence->ranktaxonid > 0 )
				s_intaxa = count_taxons(pparams->sortedpairs[p_set->set_pairs[i]]->pcr.amplifias[j].sequence->ranktaxonid);

	}
	//fprintf(stderr, "%d/%d\n", s_intaxa, pparams->options->intaxa);
	p_set->set_intaxa = s_intaxa;
	cov = s_intaxa*1.0/(pparams->options->intaxa*1.0);
	count_taxons (-1);
	return cov;
}

void set_cov_spc (pairset *p_set, SetParams *pparams)
{
	int32_t i;
	int32_t ssp = 0;
	
	count_taxons (-1);
	for (i = 0; i < pparams->options->dbsize; i++)
		if (p_set->set_wellIdentifiedTaxa[i] == 1)
			ssp = count_taxons (pparams->seqdb[i]->ranktaxonid);
		
	//set coverage
	p_set->set_coverage = get_set_coverage (p_set, pparams, -1);
	
	//set specificity
	p_set->set_specificity = ((float)ssp)/((float)p_set->set_intaxa);
	p_set->set_wi_cnt = ssp;
	count_taxons (-1);
}

//currently used only to open dead lock
void tabu_aspiration (pairset *pair_set)
{
	int i;
	int count = 0;
	
	for (i=0; i<PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] != -1)
			count++;
		if (TabuList[i] != -1)
			count++;
	}
	
	if (count == total_pairs)
		for (i=0; i<PRIMERS_IN_SET_COUNT; i++)
			TabuList[i] = -1;
}

int ok_to_add (int id, pairset *pair_set, SetParams *pparams)
{
	int i;
	int secnt = 0;
	float lnk_prcnt = 0;
	
	static int effency_switch = 0;
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (pair_set->set_pairs[i] == id) return 0;
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (pair_set->set_pairs[i] != -1) secnt++;
	
	//if (secnt == PRIMERS_IN_SET_COUNT) return 0;
	if (secnt == 0) return 1;
	
	lnk_prcnt = 1.0;
	if (secnt > pparams->options->links_cnt)
		lnk_prcnt = (pparams->options->links_cnt*1.0)/(secnt*1.0);
	
	//TR 6/2/11: new elements must have some links with atleast one elem of set
	if (get_links_distribution (id, pair_set, pparams) < lnk_prcnt) return 0;
	
	//if in tabu search search tabu list as well
	if (next_tabu_slot != -1)
	{
		//effency_switch is only there to avoid tabu_aspiration execution
		//each time
		effency_switch++;
		if ((effency_switch%5) == 0)
		{
			effency_switch=1;
			tabu_aspiration(pair_set);
		}
		
		i = find_in_tabu (id);
		if (i != -1) return 0;
	}
	return 1;
}

pairset build_primers_set_greedy_cov (SetParams *params)
{
	pairset pair_set;
	int32_t i;
	int32_t pset_idx;
	int32_t prb_idx;
	
	memset (&pair_set, 0, sizeof(pairset));
	pair_set.set_wellIdentifiedTaxa = (int *) ECOMALLOC(params->options->dbsize*sizeof (int),
            "Could not allocate memory for pair set");

	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		pair_set.set_pairs[i] = -1;
	
	pset_idx = 0;
	prb_idx = 0;
	//add first pair by default, this should be the one having highiest specificty
	add_pair_in_set (&pair_set, pset_idx, prb_idx, params);
	pset_idx++;

	while (pset_idx < PRIMERS_IN_SET_COUNT)
	{
		if (pair_set.set_coverage == 1.0) break;
		prb_idx = get_next_option_increasing_cov (&pair_set, params);
		if (prb_idx == 0) break;
		add_pair_in_set (&pair_set, pset_idx, prb_idx, params);
		pset_idx++;
	}
	//get_set_mean_cov_stats (&pair_set, &params);
	reset_set_props (&pair_set, params);

	return pair_set;
}

int32_t get_next_option_increasing_cov (pairset *pair_set, SetParams *pparams)
{
	float ini_set_cov;
	float set_cov;
	float cov_diff = 0.;
	int32_t i, id_slot = -1, max_id = 0;
	
	//coverage already 1, dont proceed.
	if (pair_set->set_coverage == 1.0) return 0;
	
	//find next available slot in the set
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (pair_set->set_pairs[i] == -1)
		{
			id_slot = i;
			break;
		}
	//set already full
	if (id_slot == -1) return 0;
	//save original set coverage
	ini_set_cov = pair_set->set_coverage;
	for (i = 1; i < pparams->sorted_count; i++)
	{
		if (ok_to_add (i, pair_set, pparams) == 0) continue;
		pair_set->set_pairs[id_slot] = i;
		set_cov = get_set_coverage (pair_set, pparams, -1);
		if ((set_cov - ini_set_cov) > cov_diff)
		{
			cov_diff = set_cov - ini_set_cov;
			max_id = i;
		}
	}
	pair_set->set_pairs[id_slot] = -1;
	return max_id;
}


//1. Add in set the first pair having highiest specificity
//2. For the sequences not WI by primers in set, calculate count for each pair equal to number of seqs WI by it
//3. Take the pair with highiest such count and see its links with primers already in set, 
//4. If no/insufficient links, take next pair else add pair in set
//5. repeate 3,4 untill pair gets added in set
//6. Repeate 2 to 5 until set completes

pairset build_primers_set_greedy_spc (SetParams *params)
{
	pairset pair_set;
	int32_t i;
	int32_t pset_idx;
	int32_t prb_idx;
	int *pair_wi_count_sorted_ids;
	
	memset (&pair_set, 0, sizeof(pairset));
	pair_set.set_wellIdentifiedTaxa = (int *) ECOMALLOC(params->options->dbsize*sizeof (int),
            "Could not allocate memory for pair set");
	
	pair_wi_count_sorted_ids = (int *) ECOMALLOC(params->sorted_count*sizeof (int),
            "Could not allocate memory for pair_wi_count_sorted");
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		pair_set.set_pairs[i] = -1;
	
	pset_idx = 0;
	prb_idx = 0;
	//add first pair by default, this should be the one having highiest specificty
	add_pair_in_set (&pair_set, pset_idx, prb_idx, params);
	pset_idx++;
	
	while (pset_idx < PRIMERS_IN_SET_COUNT)
	{
		//get a sorted list of pair ids with the pair well identifying most of the remaining seqs at top
		get_next_pair_options (pair_wi_count_sorted_ids, &pair_set, params);
		
		if (pair_wi_count_sorted_ids[0] == 0)
		{
			fprintf (stderr, "No further pair found, total so far %d\n", pset_idx);
			break;
		}		
		
		for (prb_idx = 0; prb_idx < params->sorted_count; prb_idx++)
		{
			if (pair_wi_count_sorted_ids[prb_idx])
			if (ok_to_add (pair_wi_count_sorted_ids[prb_idx], &pair_set, params))
			{
				//fprintf (stderr, "Oktoadd\n");
				add_pair_in_set (&pair_set, pset_idx, pair_wi_count_sorted_ids[prb_idx], params);
				pset_idx++;
			}
			if (pset_idx == PRIMERS_IN_SET_COUNT) break;
		}
		
		if (prb_idx == params->sorted_count)
		{
			fprintf (stderr, "No further pair found, total so far %d\n", pset_idx);
			break;
		}
		
		if (pair_set.set_specificity == 1.0)
		{
			fprintf (stderr, "Set Complete with total primers: %d\n", pset_idx);
			break;
		}

	}
	//get_set_mean_cov_stats (&pair_set, &params);
	reset_set_props (&pair_set, params);
	return pair_set;
}

float get_links_distribution (int prb_idx, pairset *pair_set, SetParams *pparams)
{
	int i, j;
	int *pair_link_count;
	int *pwi;
	int *pswi;
	float pcnt = 0.0;
	float pscnt = 0.0;
	
	pair_link_count = (int *) ECOMALLOC(PRIMERS_IN_SET_COUNT*sizeof (int),
	            "Could not allocate memory for pair_link_count");
	
	pwi = pparams->sortedpairs[prb_idx]->wellIdentifiedSeqs;
	//fprintf(stderr, "a,");
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		pair_link_count[i] = 0;
		if (pair_set->set_pairs[i] != -1)
		{
			//fprintf(stderr, "b,");
			pswi = pparams->sortedpairs[pair_set->set_pairs[i]]->wellIdentifiedSeqs;
			for (j = 0; j < pparams->options->dbsize; j++)
				if (pwi[j] == 1 && pwi[j] == pswi[j])
					pair_link_count[i] += 1;
		}
	}
	//fprintf(stderr, "c,");
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] != -1)
			pscnt++;
		if (pair_link_count[i] > 0)
			pcnt++;
	}
	ECOFREE (pair_link_count, "Could not free memory for pair_link_count");
	//fprintf(stderr, "d,");
	return (pcnt/pscnt);
}


void get_next_pair_options (int *pair_wi_count_sorted_ids, pairset *pair_set, SetParams *pparams)
{
	int *pair_count;
	int32_t i, j;
	int max;
	int tmp;
	
	pair_count = (int *) ECOMALLOC(pparams->sorted_count*sizeof (int),
            "Could not allocate memory for pair_count");
	
	memset (pair_wi_count_sorted_ids, 0, pparams->sorted_count*sizeof(int));
	
	for (i = 0; i < pparams->options->dbsize; i++)
	{
		if (pair_set->set_wellIdentifiedTaxa[i] == 1) continue;
		
		for (j = 0; j < pparams->sorted_count; j++)
		{
			if (pparams->sortedpairs[j]->wellIdentifiedSeqs[i] == 1)
				pair_count[j] += 1;
		}
	}
	
	//set pair ids
	for (j = 0; j < pparams->sorted_count; j++)
		pair_wi_count_sorted_ids[j] = j;
	
	//set count of primers already in set to zero (it should already be zero) just for testing
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (pair_set->set_pairs[i] != -1)
			pair_count[pair_set->set_pairs[i]] = 0;
	
	//sort two arrays in descending wrt count
	for (i = 0; i < pparams->sorted_count - 1; i++)
	{
		max = i;
		for (j = i + 1; j < pparams->sorted_count; j++)
			if (pair_count[max] < pair_count[j])
				max = j;
		
		if (max > i)
		{
			tmp = pair_count[i];
			pair_count[i] = pair_count[max];
			pair_count[max] = tmp;
			
			tmp = pair_wi_count_sorted_ids[i];
			pair_wi_count_sorted_ids[i] = pair_wi_count_sorted_ids[max];
			pair_wi_count_sorted_ids[max] = tmp;
		}
	}
	
	for (i = 0; i < pparams->sorted_count - 1; i++)
		if (pair_count[i] == 0)
			pair_wi_count_sorted_ids[i] = 0;	
		//else
		//	fprintf (stderr, "%d:%d, ", i, pair_count[i]);
	
	ECOFREE (pair_count, "Could not free memory for pair_count");
}

void add_pair_in_set (pairset *pair_set, int32_t pset_idx, int32_t prb_idx, SetParams *pparams)
{
	int *pwi;
	int32_t i;
	
	if (prb_idx < 0 || prb_idx >= pparams->sorted_count) return;
	pair_set->set_pairs[pset_idx] = prb_idx;
	
//	fprintf (stderr, "%d:", prb_idx);
	//fprintf (stderr, "%d:", pparams->sortedpairs[prb_idx]);
	//fprintf (stderr, "%d:", pparams->sortedpairs[prb_idx]->wellIdentifiedSeqs);
	//fprintf (stderr, "%d:%d, ", i, pair_count[i]);
	
	pwi = pparams->sortedpairs[prb_idx]->wellIdentifiedSeqs;
	for (i = 0; i < pparams->options->dbsize; i++)
		if (pwi[i] == 1)
			pair_set->set_wellIdentifiedTaxa[i] = 1;
	
	set_cov_spc (pair_set, pparams);
}

int isthisset (pairset *pair_set)
{
	int set1[PRIMERS_IN_SET_COUNT];
	int set2[PRIMERS_IN_SET_COUNT];
	int i,j=0,k;
	for (i=0; i<PRIMERS_IN_SET_COUNT; i++)
	{
		set1[i]=-1;set2[i]=-1;
		if (pair_set->set_pairs[i] != -1) set2[j++]=pair_set->set_pairs[i];
	}
	set1[0]=0;set1[1]=2;set1[2]=7;
	if (j==3)
	{
		for(i=0;i<3;i++)
		{
			for(k=0;k<3;k++)
				if (set2[k]==set1[i])break;
			if (k==3)break;
		}
		if(i==3) return 1;//found
	}
	
	set1[0]=0;set1[1]=2;set1[2]=37;
	if (j==3)
	{
		for(i=0;i<3;i++)
		{
			for(k=0;k<3;k++)
				if (set2[k]==set1[i])break;
			if (k==3)break;
		}
		if(i==3) return 1;//found
	}
	return 0;
}

void get_set_mean_cov_stats (pairset *pair_set, SetParams *pparams)
{
	int32_t i, j, k;
	int interseq_vals[PRIMERS_IN_SET_COUNT*PRIMERS_IN_SET_COUNT];
	int interseq_cnt = 0;
	double msum;
	int *p1wi;
	int *p2wi;
	int dbg=0;
	
	dbg=isthisset(pair_set);
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] == -1) continue;
		p1wi = pparams->sortedpairs[pair_set->set_pairs[i]]->wellIdentifiedSeqs;
		
		if(dbg)
		{
			printf ("\n\nWellIdentified for primer pair:%d\n",pair_set->set_pairs[i]);
			for (k = 0; k < pparams->options->dbsize; k++)
				if(p1wi[k]==1)
				printf("%d,",k);
			printf("\n");
		}
		
		for (j = i+1; j < PRIMERS_IN_SET_COUNT; j++)
		{
			if (pair_set->set_pairs[j] == -1) continue;
			p2wi = pparams->sortedpairs[pair_set->set_pairs[j]]->wellIdentifiedSeqs;
			interseq_vals[interseq_cnt] = 0;
			
			if (dbg)
			{
				printf ("Intersection for %d and %d:\n", pair_set->set_pairs[i], pair_set->set_pairs[j]);
			}
			
			for (k = 0; k < pparams->options->dbsize; k++)
				if (p1wi[k] == 1 && p2wi[k] == 1)
				{
					interseq_vals[interseq_cnt]++;
					if(dbg)
						printf("%d,",k);
				}
			if(dbg)
				printf("\n");
			interseq_cnt++;
		}
	}
	
	//calculate mean
	msum = 0;
	pair_set->set_score = 0;
	pair_set->set_lmean = 0;
	pair_set->set_lcov = -1;
	if (interseq_cnt == 0) return;
	
	for (i = 0; i < interseq_cnt; i++)
		msum += interseq_vals[i];
	pair_set->set_lmean = msum/interseq_cnt;
	
	msum = 0;
	for (i = 0; i < interseq_cnt; i++)
		msum += (interseq_vals[i] - pair_set->set_lmean)*(interseq_vals[i] - pair_set->set_lmean);
	pair_set->set_lcov = msum/interseq_cnt;
	
	if (pair_set->set_lcov != 0)
		//pair_set->set_score = (pair_set->set_coverage*pair_set->set_specificity*pair_set->set_specificity*(sqrt(pair_set->set_lmean)))/sqrt(sqrt(pair_set->set_lcov));
		pair_set->set_score = (pair_set->set_coverage*pair_set->set_specificity*(pair_set->set_lmean))/sqrt(pair_set->set_lcov);
}

void get_set_mean_cov_normalised_stats (pairset *pair_set, SetParams *pparams)
{
	int32_t i, j, k;
	int interseq_vals[PRIMERS_IN_SET_COUNT*PRIMERS_IN_SET_COUNT];
	int interseq_cnt = 0;
	double msum;
	int *p1wi;
	int *p2wi;
	int dbg=0;

	dbg=isthisset(pair_set);

	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] == -1) continue;
		p1wi = pparams->sortedpairs[pair_set->set_pairs[i]]->wellIdentifiedSeqs;

		if(dbg)
		{
			printf ("\n\nWellIdentified for primer pair:%d\n",pair_set->set_pairs[i]);
			for (k = 0; k < pparams->options->dbsize; k++)
				if(p1wi[k]==1)
				printf("%d,",k);
			printf("\n");
		}

		for (j = i+1; j < PRIMERS_IN_SET_COUNT; j++)
		{
			if (pair_set->set_pairs[j] == -1) continue;
			p2wi = pparams->sortedpairs[pair_set->set_pairs[j]]->wellIdentifiedSeqs;
			interseq_vals[interseq_cnt] = 0;

			if (dbg)
			{
				printf ("Intersection for %d and %d:\n", pair_set->set_pairs[i], pair_set->set_pairs[j]);
			}

			for (k = 0; k < pparams->options->dbsize; k++)
				if (p1wi[k] == 1 && p2wi[k] == 1)
				{
					interseq_vals[interseq_cnt]++;
					if(dbg)
						printf("%d,",k);
				}
			if(dbg)
				printf("\n");
			interseq_cnt++;
		}
	}

	//calculate mean
	msum = 0;
	pair_set->set_score = 0;
	pair_set->set_lmean = 0;
	pair_set->set_lcov = -1;
	if (interseq_cnt == 0) return;

	for (i = 0; i < interseq_cnt; i++)
		msum += interseq_vals[i];
	pair_set->set_lmean = msum/interseq_cnt;

	msum = 0;
	for (i = 0; i < interseq_cnt; i++)
		msum += (interseq_vals[i] - pair_set->set_lmean)*(interseq_vals[i] - pair_set->set_lmean);
	pair_set->set_lcov = msum/interseq_cnt;

	if (pair_set->set_lcov != 0)
	{
		//normalised links
		double nl = pair_set->set_lmean/sqrt (pair_set->set_lcov);
		nl = nl/pparams->options->insamples; //max links cannot be more than insample value
		pair_set->set_score = pair_set->set_coverage*pair_set->set_specificity*nl;
	}
}

void reset_set_props (pairset *pair_set, SetParams *pparams)
{
	int *pwi;
	int i, j;
	
	memset (pair_set->set_wellIdentifiedTaxa, 0, pparams->options->dbsize*sizeof(int));
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] == -1) continue;
		pwi = pparams->sortedpairs[pair_set->set_pairs[i]]->wellIdentifiedSeqs;
		for (j = 0; j < pparams->options->dbsize; j++)
			if (pwi[j] == 1)
				pair_set->set_wellIdentifiedTaxa[j] = 1;
	}
	
	set_cov_spc (pair_set, pparams);
	
	//TR 6/2/11: commented following, now score is just product of spc and cov
	//get_set_mean_cov_stats (pair_set, pparams);
	//get_set_mean_cov_normalised_stats (pair_set, pparams);
	//pair_set->set_score = pair_set->set_coverage*pair_set->set_specificity;
	pair_set->set_score = pair_set->set_coverage;
	//pair_set->set_score = pair_set->set_specificity;
}

void print_set_info (pairset *pair_set, SetParams *pparams)
{
	int i;
	int printed1st = 0;
	
	//TR 6/2/11: commented following, now score is just product of spc and cov
	//printf ("%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t", pair_set->set_specificity,
	//		pair_set->set_coverage,pair_set->set_lmean,
	//		sqrt(pair_set->set_lcov),pair_set->set_score);
	
	printf ("%4.3f\t%4.3f\t%4.3f\t", pair_set->set_coverage,
			pair_set->set_specificity,pair_set->set_score);
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] == -1) continue;
		
		if (printed1st)
			printf (":%d", pair_set->set_pairs[i]);
		else
			printf ("%d", pair_set->set_pairs[i]);
		printed1st = 1;
	}
	printf ("\t%d\t%d", pair_set->set_intaxa, pair_set->set_wi_cnt);
	printf ("\n");
}

void some_other_set_possibilities (pairset *pair_set,
		ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options)
{
	SetParams params;
	pairset tmp_pair_set;
	int i, j;
	int *wi;
	
	params.sortedpairs = sortedpairs;
	params.sorted_count = sorted_count;
	params.seqdb = seqdb;
	params.options = options;
	wi = (int *) ECOMALLOC(options->dbsize*sizeof (int),
					            "Could not allocate memory for pair set wi");
	//print stats for first original set
	printf ("\nspecificity\tcoverage\tmean\tcovariance\tscore\tprimers\n");
	print_set_info (pair_set, &params);
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (pair_set->set_pairs[i] == -1) continue;
			
		for (j = 0; j < sorted_count; j++)
		{
			if (ok_to_add (j, pair_set, &params))
			{
				tmp_pair_set = *pair_set;
				tmp_pair_set.set_pairs[i] = j;
				memset (wi, 0, options->dbsize*sizeof (int));
				tmp_pair_set.set_wellIdentifiedTaxa = wi;
				reset_set_props (&tmp_pair_set, &params);
				print_set_info (&tmp_pair_set, &params);
			}
		}
	}
	ECOFREE (wi, "Could not free memory for pair set wi");
}

pairset clone_set (pairset *s, int32_t dbsize)
{
	pairset clone;
	
	clone = *s;
	clone.set_wellIdentifiedTaxa = (int *) ECOMALLOC(dbsize*sizeof (int),
            "Could not allocate memory for pair set");
	memcpy (clone.set_wellIdentifiedTaxa, s->set_wellIdentifiedTaxa, dbsize*sizeof (int));
	return clone;
}

void add_in_tabu (int index)
{
	if (next_tabu_slot == -1) return;
	
	TabuList[next_tabu_slot] = index;
	next_tabu_slot++;
	if (next_tabu_slot >= PRIMERS_IN_SET_COUNT) next_tabu_slot = 0;
}

int find_in_tabu (int index)
{
	int i;
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (TabuList[i] == index) return i;
	
	return -1;
}

//do random changes in the seed set to generate new set
pairset get_neighbor (pairset *set, SetParams *params)
{
	int pinset = 0;
	int i, j, id, cnt;
	int how_many_to_replace;
	pairset nset;
	int replace_idx;
	
	//take the seed set as next neighbour sets
	nset = clone_set (set, params->options->dbsize);
	//see how many elements are in this set
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (nset.set_pairs[i] != -1) pinset++;
	}
	
	if (pinset == params->sorted_count) return nset;
	
	//Randomly get number of elements to be replaced
	//with new unused elements
	how_many_to_replace = rand ()%pinset + 1;
	//replace these many elements in the seed set
	for (i = 0; i < how_many_to_replace; i++)
	{
		do 
		{
			//we wil choose a random unused element as new element
			id = rand ()%params->sorted_count;
		}while (ok_to_add (id, &nset, params) == 0);
		//again choose a random element in the set to replace
		replace_idx = rand ()%pinset+1;
		cnt = 0;
		for (j = 0; j < PRIMERS_IN_SET_COUNT; j++)
		{
			if (nset.set_pairs[j] != -1) cnt++;
			if (cnt == replace_idx)
			{
				if (next_tabu_slot != -1)
					add_in_tabu (nset.set_pairs[j]);
					
				nset.set_pairs[j] = id;
				break;
			}
		}
	}
	reset_set_props (&nset, params);
	return nset;
}

//remove a random number of least contributing elements
//from the seed set with random elements from the remaining
pairset get_neighbor4 (pairset *set, SetParams *params)
{
	int pinset = 0;
	int i, j, id, id2;
	int how_many_to_replace;
	pairset nset;
	//int replace_idx;
	float contribution[PRIMERS_IN_SET_COUNT];
	int usedids[PRIMERS_IN_SET_COUNT];
	float sscore;
	float leastContri;
	//float lastLeastcontri;
	int k, l;

	//take the seed set as next neighbour sets
	nset = clone_set (set, params->options->dbsize);
	//see how many elements are in this set
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (nset.set_pairs[i] != -1) pinset++;
	}
	
	if (pinset == params->sorted_count) return nset;
	
	//Randomly get number of elements to be replaced
	//with new unused elements
	how_many_to_replace = rand ()%pinset + 1;
	
	//calculate contribution of each element in the set
	sscore = nset.set_score;
	//fprintf (stderr, "{%f-", sscore);
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		contribution[i] = 10000.0;
		if (nset.set_pairs[i] != -1)
		{
			id = nset.set_pairs[i];
			nset.set_pairs[i] = -1;
			reset_set_props (&nset, params);
			contribution[i] = sscore - nset.set_score;
			//fprintf (stderr, "[%f;%f]", nset.set_score, contribution[i]);
			nset.set_pairs[i] = id;
			reset_set_props (&nset, params);
			//fprintf (stderr, "%f:, ", nset.set_score);
		}
		usedids[i] = -1;
	}
	
	//lastLeastcontri = 10000.0;
	k=0;
	//replace these many elements in the seed set
	//fprintf (stderr, "} (%d) ", how_many_to_replace);
	for (i = 0; i < how_many_to_replace; i++)
	{
		do 
		{
			//we wil choose a random unused element as new element
			id = rand ()%params->sorted_count;
		}while (ok_to_add (id, &nset, params) == 0);
		
		leastContri = 10000.0;
		id2 = -1;
		for (j = 0; j < PRIMERS_IN_SET_COUNT; j++)
		{
			if (nset.set_pairs[j] == -1) continue;
			for (l = 0; l < k; l++)
				if (usedids[l] == j) break;
			
			if (leastContri > contribution[j] /*&& leastContri >= lastLeastcontri*/ && l == k)
			{
				leastContri = contribution[j];
				id2 = j;
				//fprintf (stderr, "%f:%d, ", leastContri, id2);
			}
		}
		
		if (id2 != -1)
		{
			usedids[k++] = id2;
			//lastLeastcontri = leastContri;
			if (next_tabu_slot != -1)
				add_in_tabu (nset.set_pairs[id2]);

			nset.set_pairs[id2] = id;
		}
	}
	//fprintf (stderr, "\n");
	reset_set_props (&nset, params);
	return nset;
}

//by replacing one element of the set with next from unused
int which_one_to_replace;
int last_replaced_with;
pairset get_neighbor2 (pairset *set, SetParams *params)
{
	int pinset = 0;
	int i, j;
	pairset nset;
	
	//take the seed set as next neighbour sets
	nset = clone_set (set, params->options->dbsize);
	//see how many elements are in this set
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (nset.set_pairs[i] != -1) pinset++;
	}
	
	if (pinset == params->sorted_count) return nset;
	
	if (which_one_to_replace == pinset) which_one_to_replace = 0;
	if (last_replaced_with == params->sorted_count) last_replaced_with = 0;
	
	j = -1;
	for (i = 0; i < pinset; i++)
	{
		if (nset.set_pairs[i] != -1) j++;
		if (j == which_one_to_replace) break;
	}
	
	for (j = last_replaced_with; j < params->sorted_count; j++)
		if (ok_to_add (j, &nset, params) == 1) break;
	
	if (j < params->sorted_count)
	{
		if (next_tabu_slot != -1)
			add_in_tabu (nset.set_pairs[i]);

		nset.set_pairs[i] = j;
		reset_set_props (&nset, params);
	}
	last_replaced_with++;
	which_one_to_replace++;
	
	return nset;
}

//replace element having least contribution with the next from the list
pairset get_neighbor3 (pairset *set, SetParams *params)
{
	int pinset = 0;
	int i, j, id;
	pairset nset;
	float least_contribution;
	float sscore;
	
	//take the seed set as next neighbour sets
	nset = clone_set (set, params->options->dbsize);
	//see how many elements are in this set
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (nset.set_pairs[i] != -1) pinset++;
	}
	
	if (pinset == params->sorted_count) return nset;
	
	if (last_replaced_with == params->sorted_count) last_replaced_with = 0;
	
	sscore = nset.set_score;
	which_one_to_replace = -1;
	least_contribution = 1000.0; //impossible
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		if (nset.set_pairs[i] != -1)
		{
			id = nset.set_pairs[i];
			nset.set_pairs[i] = -1;
			reset_set_props (&nset, params);
			if ((sscore - nset.set_score) < least_contribution) 
			{
				which_one_to_replace = i;
				least_contribution = sscore - nset.set_score;
			}
			nset.set_pairs[i] = id;
		}
	}
	
	for (j = last_replaced_with; j < params->sorted_count; j++)
		if (ok_to_add (j, &nset, params) == 1) break;
	
	if (j < params->sorted_count && which_one_to_replace != -1)
	{
		if (next_tabu_slot != -1)
			add_in_tabu (nset.set_pairs[which_one_to_replace]);

		nset.set_pairs[which_one_to_replace] = j;
		reset_set_props (&nset, params);
	}
	last_replaced_with++;
	
	return nset;
}

float sa_P (float score, float nscore, float tem)
{
	if (nscore >= score) return 0.95;
	if (nscore < score && tem > 0.3) return 0.1;
	if (nscore < score && tem > 0.1) return 0.0001;
	return 0.0;
}

void extend_set_to (pairset *pair_set, SetParams *params, int extend_to_cnt)
{
	int pinset = 0;
	int i, j;
	if (extend_to_cnt > PRIMERS_IN_SET_COUNT) return;
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (pair_set->set_pairs[i] != -1) pinset++;
	
	if (pinset >= extend_to_cnt) return;
	
	for (i = 0; i < params->sorted_count; i++)
	{
		if (ok_to_add (i, pair_set, params) == 1)
		{
			for (j = 0; j < PRIMERS_IN_SET_COUNT; j++)
				if (pair_set->set_pairs[j] == -1) break;
			add_pair_in_set (pair_set, j, i, params);
			reset_set_props (pair_set, params);
			pinset++;
			if (pinset >= extend_to_cnt) break;
		}
	}
}

pairset * extend_set_randomly (pairset *pset, SetParams *params, int extend_to_cnt)
{
	int pinset = 0;
	int i = 0;
	int id;
	pairset *pair_set;
	
	if (extend_to_cnt > PRIMERS_IN_SET_COUNT) return pset;
	
	if (pset == NULL)
	{
		pair_set =  (pairset *) ECOMALLOC(sizeof (pairset),
			            "Could not allocate memory for pair");
		pair_set->set_wellIdentifiedTaxa = (int *) ECOMALLOC(params->options->dbsize*sizeof (int),
			            "Could not allocate memory for pair set WI");
		
		for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
			pair_set->set_pairs[i] = -1;
	}
	else
		pair_set = pset;
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		if (pair_set->set_pairs[i] != -1) pinset++;
	
	if (pinset >= extend_to_cnt) return pair_set;
	
	i = 0;
	if (pinset == 0)
	{
		id = rand ()%params->sorted_count;
		add_pair_in_set (pair_set, i, id, params);
		i++;
		pinset++;
	}
	
	while (pinset < extend_to_cnt)
	{
		do 
		{
			//we wil choose a random unused element as new element
			id = rand ()%params->sorted_count;
		}while (ok_to_add (id, pair_set, params) == 0);
		add_pair_in_set (pair_set, i, id, params);
		i++;
		pinset++;			
	}
	return pair_set;
}

void set_reduce_to_best (pairset *pair_set, SetParams *params)
{
	int original_members[PRIMERS_IN_SET_COUNT];
	int i;
	int mcnt = 0;
	float max_score;
	int m_to_remove = -1;
	int tmem;
	
	for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
	{
		original_members[i] = pair_set->set_pairs[i];
		if (original_members[i] != -1) mcnt++;
	}
	
	while (mcnt > 3)
	{
		max_score = pair_set->set_score;
		m_to_remove = -1;
		
		for (i = 0; i < PRIMERS_IN_SET_COUNT; i++)
		{
			if (pair_set->set_pairs[i] == -1) continue;
			tmem = pair_set->set_pairs[i];
			pair_set->set_pairs[i] = -1; //remove current element temporarily
			reset_set_props (pair_set, params);
			pair_set->set_pairs[i] = tmem; //restore
			
			if (max_score <= pair_set->set_score)
			{
				max_score = pair_set->set_score;
				m_to_remove = i;
			}
		}
		
		if (m_to_remove != -1)
		{
			pair_set->set_pairs[m_to_remove] = -1; //remove element
			reset_set_props (pair_set, params);
			mcnt--;
		}
		else 
		{
			reset_set_props (pair_set, params);
			break;
		}
	}
}

float sa_temp (int k, int kmax)
{
	return ((kmax - k)*1.0)/(kmax*1.0);
}

void sets_by_SimulatedAnealing (pairset *pair_set,
		ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options)
{
	SetParams params;
	pairset set, bset, nset;
	float score, bscore, nscore;
	int k = 0, mcnt=0;
	int kmax = sorted_count*10;
	float max_score = 10000.0;
	float min_spc, min_cov;
	float max_spc, max_cov;
	double randval = 0;
	
	int32_t count_per_size;
	int32_t size_counter = 1;
	
	srand ( time(NULL) );
	params.sortedpairs = sortedpairs;
	params.sorted_count = sorted_count;
	params.seqdb = seqdb;
	params.options = options;
	which_one_to_replace = 0;
	last_replaced_with = 0;
	next_tabu_slot = -1; //dont use tabu search props
	total_pairs = sorted_count;
	
	if (pair_set == NULL)
	{
		pair_set = extend_set_randomly (NULL, &params, 3);
		printf("\nStart Random seed set for Simulated :\n");
		print_set_info (&pair_set, &params);
	}
	min_spc = max_spc = pair_set->set_specificity;
	min_cov = max_cov = pair_set->set_coverage;
	
	for (k = 0; k < PRIMERS_IN_SET_COUNT; k++)
	{
		if (pair_set->set_pairs[k] != -1) mcnt++;
	}
	count_per_size = kmax/(PRIMERS_IN_SET_COUNT-mcnt);	
	k = 1;
	
	set = clone_set (pair_set, options->dbsize);
	
	/*if (mcnt < 5)
	{
		printf ("Set before extension:\n");
		print_set_info (&set, &params);
		extend_set_to (&set, &params, 5);
		printf ("Set after extension:\n");
		print_set_info (&set, &params);
		set_reduce_to_best (&set, &params);
		printf ("Set after reduction to best size:\n");
		print_set_info (&set, &params);
	}*/
	mcnt++;
	extend_set_to (&set, &params, mcnt);
	
	bset = clone_set (&set, options->dbsize);
	score = bset.set_score;
	bscore = score;
	nset.set_wellIdentifiedTaxa = NULL;

	//srand ( time(NULL) );
	
	while (k <= kmax && score < max_score)
	{
		if (k == (size_counter*count_per_size))
		{
			size_counter++;
			mcnt++;
			extend_set_to (&set, &params, mcnt);
		}
		//nset = get_neighbor (&set, &params); //all random
		//nset = get_neighbor2 (&set, &params); //replace next with next available
		nset = get_neighbor3 (&set, &params); //replace the one with least contribution with next
		//nset = get_neighbor4 (&set, &params); //replace randome no of least contributing elements with random elements in the remaining set

		if (nset.set_specificity < min_spc)
			min_spc = nset.set_specificity;
		
		if (nset.set_specificity > max_spc)
			max_spc = nset.set_specificity;
		
		if (nset.set_coverage < min_cov)
			min_cov = nset.set_coverage;
		
		if (nset.set_coverage > max_cov)
			max_cov = nset.set_coverage;
		
		nscore = nset.set_score;
		printf ("Neighbor: ");
		print_set_info (&nset, &params);

		if (nscore > bscore)
		{
			ECOFREE (bset.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
			bset = clone_set (&nset, options->dbsize);
			bscore = nscore;
			printf ("best: ");
			print_set_info (&nset, &params);
		}
		
		randval = (double)rand()/(double)RAND_MAX;
		if (sa_P (score, nscore, sa_temp (k,kmax)) > randval)
		{
			ECOFREE (set.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
			set = clone_set (&nset, options->dbsize);
			score = nscore;
			//which_one_to_replace = 0;
			//last_replaced_with = 0;
			//printf ("Seed Set: ");
			//print_set_info (&set, &params);
		}
		k++;
		ECOFREE (nset.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
	}
	printf ("Minimum specificity: %0.3f, Maximum specificity: %0.3f, range: %0.3f\n",min_spc, max_spc, max_spc-min_spc);
	printf ("Minimum coverage: %0.3f, Maximum coverage: %0.3f, range: %0.3f\n",min_cov, max_cov, max_cov-min_cov);
}


void sets_by_TabuSearch (pairset *pair_set,
		ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options)
{
	SetParams params;
	pairset set, bset, nset;
	float bscore, nscore;
	int k = 0, mcnt=0;
	int kmax = sorted_count*10;
	float max_score = 1000.0;
	float min_spc, min_cov;
	float max_spc, max_cov;

	int32_t count_per_size;
	int32_t size_counter = 1;

	srand ( time(NULL) );
	params.sortedpairs = sortedpairs;
	params.sorted_count = sorted_count;
	params.seqdb = seqdb;
	params.options = options;
	which_one_to_replace = 0;
	last_replaced_with = 0;
	next_tabu_slot = 0; //use tabu search
	total_pairs = sorted_count;
	
	if (pair_set == NULL)
		pair_set = extend_set_randomly (NULL, &params, 3);
	min_spc = max_spc = pair_set->set_specificity;
	min_cov = max_cov = pair_set->set_coverage;

	for (k = 0; k < PRIMERS_IN_SET_COUNT; k++)
	{
		if (pair_set->set_pairs[k] != -1) mcnt++;
	}
	count_per_size = kmax/(PRIMERS_IN_SET_COUNT-mcnt);	

	set = clone_set (pair_set, options->dbsize);
	/*if (mcnt < 5)
	{
		printf ("Set before extension:\n");
		print_set_info (&set, &params);
		extend_set_to (&set, &params, 5);
		printf ("Set after extension:\n");
		print_set_info (&set, &params);
		set_reduce_to_best (&set, &params);
		printf ("Set after reduction to best size:\n");
		print_set_info (&set, &params);
	}*/
	mcnt++;
	extend_set_to (&set, &params, mcnt);

	bset = clone_set (&set, options->dbsize);
	bscore = bset.set_score;
	nset.set_wellIdentifiedTaxa = NULL;
	
	for (k = 0; k < PRIMERS_IN_SET_COUNT; k++)
		TabuList[k] = -1;

	k = 1;
	while (k < kmax && bscore < max_score)
	{
		if (k == (size_counter*count_per_size))
		{
			size_counter++;
			mcnt++;
			extend_set_to (&set, &params, mcnt);
		}

		//nset = get_neighbor (&set, &params); //all random
		//nset = get_neighbor2 (&set, &params); //replace next with next available
		nset = get_neighbor3 (&set, &params); //replace the one with least contribution with next
		//nset = get_neighbor4 (&set, &params); //replace randome no of least contributing elements with random elements in the remaining set

		if (nset.set_specificity < min_spc)
			min_spc = nset.set_specificity;
		
		if (nset.set_specificity > max_spc)
			max_spc = nset.set_specificity;
		
		if (nset.set_coverage < min_cov)
			min_cov = nset.set_coverage;
		
		if (nset.set_coverage > max_cov)
			max_cov = nset.set_coverage;

		nscore = nset.set_score;
		printf ("Neighbor: ");
		print_set_info (&nset, &params);

		if (nscore > bscore)
		{
			ECOFREE (bset.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
			bset = clone_set (&nset, options->dbsize);
			bscore = nscore;
			printf ("best: ");
			print_set_info (&nset, &params);
		}
		ECOFREE (set.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
		set = nset;
		k++;
	}
	printf ("Minimum specificity: %0.3f, Maximum specificity: %0.3f, range: %0.3f\n",min_spc, max_spc, max_spc-min_spc);
	printf ("Minimum coverage: %0.3f, Maximum coverage: %0.3f, range: %0.3f\n",min_cov, max_cov, max_cov-min_cov);
}

pairset *sets_by_BruteForce (ppair_t * sortedpairs,
//void sets_by_BruteForce (ppair_t * sortedpairs, 
		int32_t sorted_count, pecodnadb_t seqdb, poptions_t options)
{
	SetParams params;	
	params.sortedpairs = sortedpairs;
	params.sorted_count = sorted_count;
	params.seqdb = seqdb;
	params.options = options;
	
	pairset set;
	pairset *pset = NULL;
	if (sorted_count < 3)
	{
		printf ("Too few primer pairs to find a pair set.\n");
		return NULL;
	}
	
	set.set_wellIdentifiedTaxa = (int *) ECOMALLOC(options->dbsize*sizeof (int),
	            "Could not allocate memory for pair set WI");
	int32_t set_indeces_array[PRIMERS_IN_SET_COUNT];
	//int start_elements = 3;
	int end_elements = 3;
	int current_elements = 3;
	int32_t i, j;
	float maxscore = -1000.0;
	int maxcount = 2000;
	int counter = 0;
	
	if (sorted_count <= PRIMERS_IN_SET_COUNT)
		end_elements = sorted_count;
	
	if (end_elements < sorted_count)
	{
		pset =  (pairset *) ECOMALLOC(sizeof (pairset),
	            "Could not allocate memory for pair");
		pset->set_wellIdentifiedTaxa = (int *) ECOMALLOC(options->dbsize*sizeof (int),
	            "Could not allocate memory for pair set WI");
	}
	
	while (current_elements <= end_elements)
	{
		for (i=0; i<PRIMERS_IN_SET_COUNT; i++)
		{
			if (i < current_elements)
				set_indeces_array[i] = i;
			else
				set_indeces_array[i] = -1;
		}
		
		while (TRUE)
		{
			memset (set.set_wellIdentifiedTaxa, 0, options->dbsize*sizeof (int));
			for (i=0; i<PRIMERS_IN_SET_COUNT; i++)
			{
				if (i < current_elements)
					set.set_pairs[i] = set_indeces_array[i];
				else
					set.set_pairs[i] = -1;
			}
			reset_set_props (&set, &params);
			
			if (set.set_score > maxscore)
			{
				maxscore = set.set_score;
				printf ("best: ");
				print_set_info (&set, &params);
				
				ECOFREE (pset->set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
				*pset = clone_set (&set, options->dbsize);
			}
			else
			{
				printf ("set:  ");
				print_set_info (&set, &params);
			}
			
			for (i=current_elements-1; i>0; i--)
			{
				set_indeces_array[i]++;
				if (set_indeces_array[i] == (sorted_count+(i-current_elements+1)))
				{
					set_indeces_array[i] = set_indeces_array[i-1]+2;
					for (j=i+1; j<current_elements; j++)
						set_indeces_array[j]=set_indeces_array[j-1]+1;
				}
				else break;
			}
			if (i == 0)
				set_indeces_array[i]++;
			
			for (i=0; i<current_elements; i++)
				if (set_indeces_array[i] >= sorted_count)
				{
					break;
				}
			if (i < current_elements) //above loop broken?
			{
				current_elements++;
				break;
			}
			counter++;
			if (counter > maxcount)
				break;
		}
		if (counter > maxcount)
			break;
	}
	return pset;
}

void build_and_print_sets (ppair_t * sortedpairs, int32_t sorted_count, pecodnadb_t seqdb, poptions_t options)
{
	SetParams params;	
	pairset pair_set;
	pairset *pset;
	
	params.sortedpairs = sortedpairs;
	params.sorted_count = sorted_count;
	params.seqdb = seqdb;
	params.options = options;
	
	pair_set = build_primers_set_greedy_spc (&params);
	printf("Greedy algorithm results based on specificity:\n");
	print_set_info (&pair_set, &params);
	if (pair_set.set_wellIdentifiedTaxa)
		ECOFREE (pair_set.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
	
	pair_set = build_primers_set_greedy_cov (&params);
	printf("\nGreedy algorithm results based on coverage:\n");
	print_set_info (&pair_set, &params);
	if (pair_set.set_wellIdentifiedTaxa)
		ECOFREE (pair_set.set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
	
	pset = extend_set_randomly (NULL, &params, 3);
	printf("\nStart Random seed set:\n");
	print_set_info (pset, &params);

	printf("\nResults from simulated Anealing:\n");
	sets_by_SimulatedAnealing (pset, sortedpairs, sorted_count, seqdb, options);
	printf("\nResults from Tabu Search:\n");
	sets_by_TabuSearch (pset, sortedpairs, sorted_count, seqdb, options);
}

void primers_graph_graphviz (ppair_t * sortedpairs,
		int32_t sorted_count, poptions_t options)
{
	int32_t i, j, k, total_links;
	char fileName[100] = "PrimerLinks";
	int *owi;
	int *iwi;
	int allowedtaxa;
	FILE *of;

	srand ( time(NULL) );
	sprintf (fileName, "PrimerLinks_%d.gv", rand ());
	of = fopen (fileName, "w");

	fprintf (of, "graph primerlinks {\n");
	for (i=0; i<sorted_count; i++)
	{
		owi = sortedpairs[i]->wellIdentifiedSeqs;

		for (j=i+1; j<sorted_count; j++)
		{
			iwi = sortedpairs[j]->wellIdentifiedSeqs;
			total_links = 0;

			for (k=0; k<options->dbsize; k++)
				if (owi[k] == 1 && iwi[k] == 1)
					total_links++;

			//if (total_links > 0 && ((total_links*1.0) <= (options->max_links_percent*options->dbsize)))
			//if (total_links > 0 && total_links <= options->max_links_percent)

			if((sortedpairs[i]->intaxa - sortedpairs[i]->notwellidentifiedtaxa) < (sortedpairs[j]->intaxa - sortedpairs[j]->notwellidentifiedtaxa))
				allowedtaxa = (sortedpairs[i]->intaxa - sortedpairs[i]->notwellidentifiedtaxa)/2;
			else
					allowedtaxa = (sortedpairs[j]->intaxa - sortedpairs[j]->notwellidentifiedtaxa)/2;


			//if (total_links > 5 && (total_links <= options->max_links_percent || options->max_links_percent==-1))
			//	fprintf (of, "\t%d -- %d [label=\"%d: %0.2f: %0.2f\"];\n", i, j, total_links,sortedpairs[i]->bc,sortedpairs[j]->bc );

			allowedtaxa = options->max_links_percent;
			if (total_links > 5 && total_links < allowedtaxa)
				fprintf (of, "\t%d -- %d [label=\"%d: %0.2f: %0.2f\"];\n", i, j, total_links,sortedpairs[i]->bc,sortedpairs[j]->bc );
				//fprintf (of, "\t%d\t%d\t%d\n", i, j, total_links);

			//fprintf (of, "\t%d\t%d\t%d\n", i, j, total_links);
				//fprintf (of, "\t%d -- %d;\n", i, j, total_links);
		}
	}
	fprintf (of, "}\n");
	fclose (of);
}

int32_t *addinset (int32_t *set, int32_t i, int32_t j, int32_t* slots, int32_t *index)
{
	int32_t k;

	if (*index == *slots)
	{
		*slots += 50;
		set = ECOREALLOC(set, (*slots)*sizeof (int32_t),
				            "Could not allocate memory for index set.");
	}

	if (i > -1)
	{
		for (k=0; k<*index; k++)
			if (set[k] == i) break;
		if (k== *index)
			set[(*index)++] = i;
	}

	if (j > -1)
	{
		for (k=0; k<*index; k++)
			if (set[k] == j) break;
		if (k== *index)
			set[(*index)++] = j;
	}

	return set;
}

size_t primers_changeSortedArray (ppair_t ** pairs,
		size_t sorted_count, poptions_t options)
{
	int32_t i, j, k, l, total_links;
	int *owi;
	int *iwi;
	int allowedtaxa;
	ppair_t *sortedpairs = *pairs;
	bool_t passed;

	int32_t *idx_set = NULL;
	int32_t slots=50, index=0;

	idx_set = ECOMALLOC(slots*sizeof (int32_t),
		            "Could not allocate memory for index set.");

	for (i=0; i<sorted_count; i++)
	{
		owi = sortedpairs[i]->wellIdentifiedSeqs;
		passed = FALSE;

		for (j=0; j<sorted_count; j++)
		{
			if (i == j) continue;

			iwi = sortedpairs[j]->wellIdentifiedSeqs;
			total_links = 0;

			for (k=0; k<options->dbsize; k++)
				if (owi[k] == 1 && iwi[k] == 1)
					total_links++;

			//if (total_links > 0 && ((total_links*1.0) <= (options->max_links_percent*options->dbsize)))
			//if (total_links > 0 && total_links <= options->max_links_percent)

			if((sortedpairs[i]->intaxa - sortedpairs[i]->notwellidentifiedtaxa) < (sortedpairs[j]->intaxa - sortedpairs[j]->notwellidentifiedtaxa))
				allowedtaxa = (sortedpairs[i]->intaxa - sortedpairs[i]->notwellidentifiedtaxa)/2;
			else
				allowedtaxa = (sortedpairs[j]->intaxa - sortedpairs[j]->notwellidentifiedtaxa)/2;


			//if (total_links > 5 && (total_links <= options->max_links_percent || options->max_links_percent==-1))
			//	fprintf (of, "\t%d -- %d [label=\"%d: %0.2f: %0.2f\"];\n", i, j, total_links,sortedpairs[i]->bc,sortedpairs[j]->bc );

			if (options->max_links_percent > 0)
			{
				allowedtaxa = options->max_links_percent;
				if (total_links > allowedtaxa)
					passed = TRUE;
					break;
			}
			else
			if (!(total_links > 5 && total_links <= allowedtaxa))
			{
				//idx_set = addinset (idx_set, i, j, &slots, &index);
				passed = TRUE;
				break;
			}
		}
		if (passed == TRUE)
			idx_set = addinset (idx_set, i, -1, &slots, &index);
	}

	i=-1;
	for (j=0; j<sorted_count; j++)
	{
		for (k=0; k<index; k++)
			if (j == idx_set[k]) break;
		//need to remove this element
		if (k == index)
		{
			ECOFREE (sortedpairs[j]->wellIdentifiedSeqs, "Cannot free wi for changing sorted array");
			ECOFREE (sortedpairs[j]->pcr.amplifias, "Cannot free wi for changing sorted array");
			if (i == -1) i = j;
		}
		else
		{
			if (i != -1)
				sortedpairs[i++] = sortedpairs[j];
		}
	}
	ECOFREE (idx_set, "Cannot free index set.");
	if (i != -1)
	{
		*pairs = ECOREALLOC (*pairs, i*sizeof(pair_t), "Cannot free wi for changing sorted array");
	}
	else i=sorted_count;
	return i;
}


int32_t *addinset_withLinks (int32_t *set, int32_t i, int32_t* slots, int32_t *index, ppair_t *pairs, poptions_t options)
{
	int32_t j, k, total_links;
	int *owi;
	int *iwi;
	bool_t passed = TRUE;
	int allowedtaxa;

	//see if we need to extend the set array
	if (*index == *slots)
	{
		*slots += 50;
		set = ECOREALLOC(set, (*slots)*sizeof (int32_t),
				            "Could not allocate memory for index set.");
	}

	//find no of links of current element i with all the elements
	//in the set to see that they are within limit
	owi = pairs[i]->coveredSeqs;
	for (j=0; j<*index; j++)
	{
		iwi = pairs[set[j]]->coveredSeqs;
		total_links = 0;

		for (k=0; k<options->dbsize; k++)
			if (owi[k] == 1 && iwi[k] == 1)
				total_links++;

		//if((pairs[i]->intaxa - pairs[i]->notwellidentifiedtaxa) < (pairs[set[j]]->intaxa - pairs[set[j]]->notwellidentifiedtaxa))
		//	allowedtaxa = (pairs[i]->intaxa - pairs[i]->notwellidentifiedtaxa)/2;
		//else
		//		allowedtaxa = (pairs[set[j]]->intaxa - pairs[set[j]]->notwellidentifiedtaxa)/2;

		if(pairs[i]->intaxa < pairs[set[j]]->intaxa)
			allowedtaxa = pairs[i]->intaxa/2;
		else
				allowedtaxa = pairs[set[j]]->intaxa/2;

		if (!(total_links > 5 && total_links <= allowedtaxa))
			passed = FALSE;
	}

	//links respect the limits with all set elements
	if (passed)
	{
		for (k=0; k<*index; k++)
			if (set[k] == i) break;
		if (k== *index)
			set[(*index)++] = i;
	}

	return set;
}

size_t primers_filterWithGivenLinks (ppair_t ** pairs,
		size_t sorted_count, poptions_t options)
{
	int32_t i, j, k;
	ppair_t *sortedpairs = *pairs;
	bool_t passed;

	int32_t *idx_set = NULL;
	int32_t slots=50, index=0;
	int *cov = ECOMALLOC(options->dbsize*sizeof (int),
            "Could not allocate memory for index set.");

	idx_set = ECOMALLOC(slots*sizeof (int32_t),
		            "Could not allocate memory for index set.");

	for (i=sorted_count-1; i>=0; i--)
	{
		idx_set = addinset_withLinks (idx_set, i, &slots, &index, sortedpairs, options);
	}

	i=-1;
	for (j=0; j<sorted_count; j++)
	{
		for (k=0; k<index; k++)
			if (j == idx_set[k]) break;
		//need to remove this element
		if (k == index)
		{
			ECOFREE (sortedpairs[j]->coveredSeqs, "Cannot free wi for changing sorted array");
			ECOFREE (sortedpairs[j]->wellIdentifiedSeqs, "Cannot free wi for changing sorted array");
			ECOFREE (sortedpairs[j]->pcr.amplifias, "Cannot free wi for changing sorted array");
			if (i == -1) i = j;
		}
		else
		{
			if (i != -1)
				sortedpairs[i++] = sortedpairs[j];
		}
	}
	ECOFREE (idx_set, "Cannot free index set.");
	if (i != -1)
	{
		*pairs = ECOREALLOC (*pairs, i*sizeof(pair_t), "Cannot free wi for changing sorted array");
	}
	else i=sorted_count;


	for (j=0; j<i; j++)
		for (k=0; k<options->dbsize; k++)
			if ((*pairs)[j]->coveredSeqs[k] == 1)
				cov[k] = 1;
	j=0;
	for (k=0; k<options->dbsize; k++)
		if (cov[k] == 1)
			j++;
	fprintf (stderr, "\nALL ELEMENTS COVERAGE: (%d/%d) %0.2f\n", j, options->intaxa, j*1.0/options->intaxa);
	ECOFREE (cov, "Cannot free cov");

	return i;
}
