/*
 * ecoprimer.c
 *
 *  Created on: 7 nov. 2008
 *      Author: coissac
 */

#include "libecoprimer/ecoprimer.h"
#include "libecoprimer/PrimerSets.h"
#include "libecoprimer/ahocorasick.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <dlfcn.h>
#include"libthermo/nnparams.h"
#include"libthermo/thermostats.h"


#define VERSION "0.4"
  /* TR: by default, statistics are made on species level*/
#define DEFAULTTAXONRANK "species"

static int cmpprintedpairs(const void* p1,const void* p2);
//float _Z27calculateMeltingTemperature_ (char * seq1, char * seq2);
pwordcount_t reduce_words_to_debug (pwordcount_t words, poptions_t options);
void print_wordwith_positions (primer_t prm, uint32_t seqdbsize, poptions_t options);

void* lib_handle = NULL;
float (*calcMelTemp)(char*, char*);

/* ----------------------------------------------- */
/* printout help                                   */
/* ----------------------------------------------- */
#define PP fprintf(stdout,

static void PrintHelp()
{
			  PP      "------------------------------------------\n");
			  PP      " ecoPrimer Version %s\n", VERSION);
			  PP      "------------------------------------------\n");
              PP      "synopsis : finding primers and measureing the quality of primers and barcode region\n");
              PP      "usage: ./ecoPrimer [options] \n");
              PP      "------------------------------------------\n");
              PP      "options:\n");
              PP      "-d    : [D]atabase : to match the expected format, the database\n");
              PP      "        has to be formated first by the ecoPCRFormat.py program located.\n");
              PP      "        in the ecoPCR/tools directory.\n");
              PP      "        ecoPCRFormat.py creates three file types :\n");
              PP      "            .sdx : contains the sequences\n");
              PP      "            .tdx : contains information concerning the taxonomy\n");
              PP      "            .rdx : contains the taxonomy rank\n\n");
              PP      "        ecoPrimer needs all the file type. As a result, you have to write the\n");
              PP      "        database radical without any extension. For example /ecoPrimerDB/fstvert\n\n");
              PP      "-e    : [E]rror : max error allowed by oligonucleotide (0 by default)\n\n");
              PP      "-h    : [H]elp - print <this> help\n\n");
              PP      "-i    : [I]gnore the given taxonomy id (define the counterexample taxon set).\n\n");
              PP      "-l    : minimum [L]ength : define the minimum amplication length. \n\n");
              PP      "-L    : maximum [L]ength : define the maximum amplicationlength. \n\n");
              PP      "-r    : [R]estricts the search to the given taxonomic id (restrict the example taxon set).\n\n");
              PP      "-E    : [E]xception taxid allows to indicate than some subclade of example sequences are conterexamples.\n\n");
              PP      "-c    : Consider that the database sequences are [c]ircular\n\n");
              PP      "-3 	 : Three prime strict match\n\n");
              PP      "-q    : Strict matching [q]uorum, percentage of the sequences in which strict primers are found. By default it is 70\n\n");
              PP      "-s    : [S]ensitivity quorum\n\n");
              PP      "-t    : required [t]axon level for results, by default the results are computed at species level\n\n");
              PP      "-x    : false positive quorum\n\n");
              PP      "-D    : set in [d]ouble strand mode\n\n");
              PP      "-O    : set the primer length (default 18) \n\n");
              PP      "-S    : Set in [s]ingle strand mode\n\n");
              PP      "-m    : Salt correction method for Tm computation (SANTALUCIA : 1 or OWCZARZY:2, default=1)\n\n");
              PP      "-a    : Salt contentration in M for Tm computation (default 0.05 M)\n\n");
              PP      "-U    : No multi match\n\n");
              PP      "-R    : Define the [R]eference sequence identifier (must be part of example set)\n\n");
              PP      "-A    : Print the list of all identifier of sequences present in the database\n\n");
              PP      "-f    : Remove data mining step during  strict primer identification\n\n");
              PP      "-v    : Store statistic file about memory usage during strict primer identification\n\n");
              PP      "-p    : Print sets of primers (may take several minutes after primers have been designed!)\n\n");
              PP      "-T    : Ignore pairs having specificity below this Threshold\n\n");
              PP      "\n");
              PP      "------------------------------------------\n");
              PP      "Table result description : \n");
              PP      "column  1 : serial number\n");
              PP      "column  2 : primer1\n");
              PP      "column  3 : primer2\n");
              PP      "column  4 : primer1 Tm without mismatch\n");
              PP      "column  5 : primer1 lowest Tm against exemple sequences\n");
              PP      "column  6 : primer2 Tm without mismatch\n");
              PP      "column  7 : primer2 lowest Tm against exemple sequences\n");
              PP      "column  8 : primer1 G+C count\n");
              PP      "column  9 : primer2 G+C count\n");
              PP      "column 10 : good/bad\n");
              PP      "column 11 : amplified example sequence count\n");
              PP      "column 12 : amplified counterexample sequence count\n");
              PP      "column 13 : yule\n");
              PP      "column 14 : amplified example taxa count\n");
              PP      "column 15 : amplified counterexample taxa count\n");
              PP      "column 16 : ratio of amplified example taxa versus all example taxa (Bc index)\n");
              PP      "column 17 : unambiguously identified example taxa count\n");
              PP      "column 18 : ratio of specificity unambiguously identified example taxa versus all example taxa (Bs index)\n");
              PP      "column 19 : minimum amplified length\n");
              PP      "column 20 : maximum amplified length\n");
              PP      "column 21 : average amplified length\n");
              PP      "------------------------------------------\n");
              PP		" http://www.grenoble.prabi.fr/trac/ecoPrimer/\n");
              PP      "------------------------------------------\n\n");
              PP      "\n");

}

static void ExitUsage(int stat)
{
        PP      "usage: ecoprimer [-d database] [-l value] [-L value] [-e value] [-r taxid] [-i taxid] [-R rank] [-t taxon level]\n");
        PP      "type \"ecoprimer -h\" for help\n");

        if (stat)
            exit(stat);
}

#undef  PP

void initoptions(poptions_t options)
{
	options->statistics=FALSE;
	options->filtering=TRUE;
	options->lmin=0;                   //< Amplifia minimal length
	options->lmax=1000;                   //< Amplifia maximal length
	options->error_max=3;              //**< maximum error count in fuzzy search
	options->primer_length=18;          //**< minimal length of the primers
	options->restricted_taxid=NULL;      //**< limit amplification below these taxid
	options->ignored_taxid=NULL;         //**< no amplification below these taxid
	options->exception_taxid=NULL;         //**< no amplification below these taxid
	options->prefix=NULL;
	options->reference=NULL;
	options->refseq=NULL;
	options->circular=0;
	options->doublestrand=1;
	options->strict_quorum=0.7;
	options->strict_exclude_quorum=0.1;
	options->sensitivity_quorum=0.9;
	options->false_positive_quorum=0.1;
	options->strict_three_prime=0;
	options->r=0;
	options->g=0;
	options->e=0;
	options->no_multi_match=FALSE;
	options->pnparm = NULL;
	strcpy(options->taxonrank, DEFAULTTAXONRANK);			/*taxon level for results, species by default*/
	options->saltmethod = SALT_METHOD_SANTALUCIA;
	options->salt = DEF_SALT;
	options->printAC=FALSE;
	options->print_sets_of_primers = FALSE;
	options->specificity_threshold = 0.6;
	options->links_cnt = 1;
	options->max_links_percent = -1;   /*graph only those primers having maximum 15% links*/
	options->filter_on_links = FALSE;
}

void printapair(int32_t index,ppair_t pair, poptions_t options)
{
	bool_t asdirect1=pair->asdirect1;
	bool_t asdirect2=pair->asdirect2;
	bool_t asdirecttmp;
	word_t w1=pair->p1->word;
	word_t w2=pair->p2->word;
	word_t wtmp;
	bool_t good1=pair->p1->good;
	bool_t good2=pair->p2->good;
	bool_t goodtmp;
	bool_t strand;
	uint32_t i, j;
	float temp;
	CNNParams nnparams;

	//nparam_InitParams(&nnparams, DEF_CONC_PRIMERS,DEF_CONC_SEQUENCES,DEF_SALT,SALT_METHOD_SANTALUCIA);

	char *c;
	char p1[32];
	char p2[32];

	if (!asdirect1)
		w1=ecoComplementWord(w1,options->primer_length);

	if (!asdirect2)
		w2=ecoComplementWord(w2,options->primer_length);


	if (w2 < w1)
	{
		wtmp=w1;
		w1=w2;
		w2=wtmp;

		asdirecttmp=asdirect1;
		asdirect1=asdirect2;
		asdirect2=asdirecttmp;

		goodtmp=good1;
		good1=good2;
		good2=goodtmp;
	}

	//print serial number
	printf("%6d\t",index);

	c = ecoUnhashWord(w1,options->primer_length);
	strcpy (p1, c);
	c = ecoUnhashWord(w2,options->primer_length);
	strcpy (p2, c);

	//print primer1
	printf("%s\t", p1);

	//print primer2
	printf("%s", p2);

	//print primer1 melting temperature
	printf ("\t%3.1f", pair->p1temp);

	//print minimum melting temperature of approximate versions of primer1
	printf ("\t%3.1f", pair->p1mintemp);

	//print primer2 melting temperature
	printf ("\t%3.1f", pair->p2temp);

	//print minimum melting temperature of approximate versions of primer2
	printf ("\t%3.1f", pair->p2mintemp);

	//print gc contents of primer1
	printf ("\t%d",nparam_CountGCContent(p1));

	//print gc contents of primer2
	printf ("\t%d",nparam_CountGCContent(p2));

	//print good/bad pair indicator
	printf("\t%c%c", "bG"[(int)good1],"bG"[(int)good2]);

	//print inexample count
	printf("\t%d", pair->inexample);

	//print out example count
	printf("\t%d", pair->outexample);

	//print yule
	printf("\t%4.3f", pair->yule);

	//print in taxa count
	printf("\t%d", pair->intaxa);

	//print out taxa count
	printf("\t%d", pair->outtaxa);

	//print coverage
	printf("\t%4.3f", (float)pair->bc);

	//print well identified taxa count
	printf("\t%d", pair->intaxa - pair->notwellidentifiedtaxa);

	//print specificity
	printf("\t%4.3f", pair->bs);

	//print min amplifia lenght
	printf("\t%d", pair->mind);

	//print max amplifia lenght
	printf("\t%d", pair->maxd);

	//print average amplifia lenght
	printf("\t%3.2f", (float)pair->sumd/pair->amplifiacount);

	//print amplifia information about reference sequence if specified
	if (options->refseq && pair->refsequence >=0)
	{
		printf("\t%s:",options->reference);
		strand = pair->pcr.amplifias[pair->refsequence].strand;

		if (strand)
			printf("join(");
		else
			printf("complement(");

		printf("%d..%d,%d..%d",pair->pcr.amplifias[pair->refsequence].begin - options->primer_length + 1,
				                 pair->pcr.amplifias[pair->refsequence].begin,
							  	 pair->pcr.amplifias[pair->refsequence].end + 2,
							  	 pair->pcr.amplifias[pair->refsequence].end + options->primer_length + 1
							  	 );
		printf(")");
		printf("\t");

		for (c=pair->pcr.amplifias[pair->refsequence].amplifia,
		     i=pair->pcr.amplifias[pair->refsequence].begin;
			 i<=pair->pcr.amplifias[pair->refsequence].end;
			 i++,
			 c+=(strand)? 1:-1)
			printf("%c","acgt"[(strand)? (*c):(~*c)&3]);


	}
	else
		printf("\t\t");

/*	j=0;
	for (i=0; i<options->dbsize; i++)
		if (pair->wellIdentifiedSeqs[i] == 1)
			j++;
	printf("%d", j);*/

	printf("\n");

}

static int cmpprintedpairs(const void* p1,const void* p2)
{
	float s1,s2;
	ppair_t pair1,pair2;

	pair1=*((ppair_t*)p1);
	pair2=*((ppair_t*)p2);

	s1 = pair1->yule * pair1->bs;
	s2 = pair2->yule * pair2->bs;

//	fprintf(stderr,"s1 : %4.3f  %4.3f   %4.3f\n",pair1->yule , pair1->bs,s1);
//	fprintf(stderr,"s2 : %4.3f  %4.3f   %4.3f\n\n",pair2->yule , pair2->bs,s2);

	if (s1 > s2) return -1;
	if (s1 < s2) return 1;
	return 0;
}

uint32_t filterandsortpairs(ppair_t* sortedpairs,uint32_t count, poptions_t options, pecodnadb_t seqdb)
{
	uint32_t i,j;
	float q,qfp;

	for (i=0,j=0;i < count;i++)
	{
		if (options->insamples)
			q = (float)sortedpairs[i]->inexample/options->insamples;
		else q=1.0;

		if (options->outsamples)
			qfp = (float)sortedpairs[i]->outexample/options->outsamples;
		else qfp=0.0;

		sortedpairs[i]->wellIdentifiedSeqs = NULL; //TR 05/09/10 - wellIdentified needed for primer sets
		sortedpairs[i]->coveredSeqs = NULL; //TR 05/09/10 - wellIdentified needed for primer sets
		sortedpairs[i]->quorumin = q;
		sortedpairs[i]->quorumout = qfp;
		sortedpairs[i]->yule = q - qfp;
		sortedpairs[j]=sortedpairs[i];

		if (q > options->sensitivity_quorum &&
			qfp < options->false_positive_quorum)
		{
			//TR 05/09/10 - wellIdentified needed for primer sets
			sortedpairs[j]->wellIdentifiedSeqs = ECOMALLOC(options->dbsize * sizeof(int),"Cannot allocate well_identified_array");
			sortedpairs[j]->coveredSeqs = ECOMALLOC(options->dbsize * sizeof(int),"Cannot allocate well_identified_array");
			(void)taxonomycoverage(sortedpairs[j],options, seqdb, options->dbsize);
			taxonomyspecificity(sortedpairs[j], seqdb, options->dbsize);
			//j++;
			//if specificity less than user provieded threshold (default 60%) then ignore this pair
			if (sortedpairs[j]->bs >= options->specificity_threshold)
				j++;
		}

	}
	qsort(sortedpairs,j,sizeof(ppair_t),cmpprintedpairs);
	return j;
}


void printpairs (ppairtree_t pairs, poptions_t options,ecotaxonomy_t *taxonomy, pecodnadb_t seqdb)
{
   ppair_t* sortedpairs;
   ppair_t* index;
   ppairlist_t pl;
   size_t i,j;
   size_t count;
   char *taxon[]={"taxon","taxa"};
   ecotx_t    *current_taxon;
   //pairset pair_sets;
   pairset *pset = NULL;

   //printf("Index\tPrimer1\tPrimer2\tGB\tInexampleCount\tOutexampleCount\tYule\tIntaxaCount\tOuttaxaCount\tCoverage\tSpecificity\tMinAmplifiedLength\tMaxAmplifiedLength\tAvgAmplifiedLength\n");

   fprintf(stderr,"Total pair count : %d\n",pairs->count);

   sortedpairs = ECOMALLOC(pairs->count*sizeof(ppair_t),"Cannot Allocate ordered pairs");
   index=sortedpairs;
   pl=pairs->first;
   j=0;
   while(pl->next)
   {
	   for (i=0;i<pl->paircount;i++,j++)
		   sortedpairs[j]=pl->pairs+i;
	   pl=pl->next;
   }

   for (i=0;i<pl->paircount;i++,j++)
	   sortedpairs[j]=pl->pairs+i;
   
   count=filterandsortpairs(sortedpairs,pairs->count,options, seqdb);
   getThermoProperties(sortedpairs, count, options);

   fprintf(stderr,"Total good pair count : %u\n",(uint32_t)count);

	printf("#\n");
	printf("# ecoPrimer version %s\n",VERSION);
	printf("# Rank level optimisation : %s\n", options->taxonrank);
	printf("# max error count by oligonucleotide : %d\n",options->error_max);
	printf("#\n");

	if (options->r)
	{
		printf("# Restricted to %s:\n",taxon[(options->r>1) ? 1:0]);
		for(i=0;i<(uint32_t)options->r;i++)
		{
			current_taxon=eco_findtaxonbytaxid(taxonomy,options->restricted_taxid[i]);
			printf("#     %d : %s (%s)\n", current_taxon->taxid,
										   current_taxon->name,
										   taxonomy->ranks->label[current_taxon->rank]
										 );
		}
		printf("#\n");
	}
	if (options->g)
	{
		printf("# Ignore  %s:\n",taxon[(options->g>1) ? 1:0]);
		for(i=0;i<(uint32_t)options->g;i++)
		{
			current_taxon=eco_findtaxonbytaxid(taxonomy,options->ignored_taxid[i]);
			printf("#     %d : %s (%s)\n", current_taxon->taxid,
										   current_taxon->name,
										   taxonomy->ranks->label[current_taxon->rank]
										 );

		}
		printf("#\n");
	}

	printf("# strict primer quorum  : %3.2f\n",options->strict_quorum);
	printf("# example quorum        : %3.2f\n",options->sensitivity_quorum);
	if (options->g + options->r)
		printf("# counterexample quorum : %3.2f\n",options->false_positive_quorum);


	printf("#\n");
	printf("# database : %s\n",options->prefix);
    printf("# Database is constituted of %5d examples        corresponding to %5d %s\n",options->insamples,
    		options->intaxa,options->taxonrank);
    printf("#                        and %5d counterexamples corresponding to %5d %s\n",options->outsamples,
    		options->outtaxa,options->taxonrank);
	printf("#\n");

	if (options->lmin && options->lmax)
		printf("# amplifiat length between [%d,%d] bp\n",options->lmin,options->lmax);
	else if (options->lmin)
		printf("# amplifiat length larger than %d bp\n",options->lmin);
	else if (options->lmax)
		printf("# amplifiat length smaller than %d bp\n",options->lmax);
	if (options->circular)
		printf("# DB sequences are considered as circular\n");
	else
		printf("# DB sequences are considered as linear\n");
	printf("# Pairs having specificity less than %0.2f will be ignored\n", options->specificity_threshold);
	printf("#\n");
	

   for (i=0;i < count;i++)
	   printapair(i,sortedpairs[i],options);

   if (options->filter_on_links)
   {
		fprintf (stderr, "Old size: %d, ", count);
		count = primers_changeSortedArray (&sortedpairs, count, options);
		//count = primers_filterWithGivenLinks (&sortedpairs, count, options);
		fprintf (stderr, "New size: %d\n", count);

		if (count == 0)
		{
			fprintf (stderr, "No pairs passed the links constraints.\n");
			printf ("No pairs passed the links constraints.\n");
			return;
		}

		for (i=0;i < count;i++)
			printapair(i,sortedpairs[i],options);
   }
   
   if (options->print_sets_of_primers == TRUE)
   {
	   /*pair_sets = build_primers_set (sortedpairs, count, seqdb, options);
	   printf("Results from Greedy Algorithm and some other possibilities:\n");
	   some_other_set_possibilities (&pair_sets, sortedpairs, count, seqdb, options);
	   printf("Results from simulated Anealing:\n");
	   sets_by_SimulatedAnealing (&pair_sets, sortedpairs, count, seqdb, options);
	   printf("Results from Tabu Search:\n");
	   sets_by_TabuSearch (&pair_sets, sortedpairs, count, seqdb, options);*/
	   //pset = sets_by_BruteForce (sortedpairs, count, seqdb, options);
	   //if (pset)
	   /*/{
		   printf("Results from simulated Anealing:\n");
		   sets_by_SimulatedAnealing (pset, sortedpairs, count, seqdb, options);
		   printf("Results from Tabu Search:\n");
		   sets_by_TabuSearch (pset, sortedpairs, count, seqdb, options);
		   
		   if (pset)
		   {
			   ECOFREE (pset->set_wellIdentifiedTaxa, "Could not free memory for pair set wi");
			   ECOFREE (pset, "Could not free memory for pair");
		   }
	   }*/
	   build_and_print_sets (sortedpairs, count, seqdb, options);
   }
   //primers_graph_graphviz (sortedpairs, count, options);
}


/*updateseqparams: This function counts the insample and outsample sequences
 *  and with each sequences adds a tag of the taxon to which the sequence beongs*/

void updateseqparams (pecodnadb_t seqdb, uint32_t seqdbsize, ecotaxonomy_t *taxonomy,
		poptions_t options, int32_t *insamples, int32_t *outsamples)
{
	uint32_t i;
	int32_t taxid;
	ecotx_t  *tmptaxon;

    for (i=0;i<seqdbsize;i++)
    {
    	seqdb[i]->isexample=isExampleTaxon(taxonomy,seqdb[i]->taxid,options);
    	if (seqdb[i]->isexample)
    		(*insamples)++;
    	else
    		(*outsamples)++;

    	taxid = taxonomy->taxons->taxon[seqdb[i]->taxid].taxid;
		tmptaxon = eco_findtaxonbytaxid(taxonomy, taxid);
		if (tmptaxon)
			tmptaxon = eco_findtaxonatrank(tmptaxon, options->taxonrankidx);
		if (tmptaxon)
			seqdb[i]->ranktaxonid = tmptaxon->taxid;
    }
}

void setresulttaxonrank (ecotaxonomy_t *taxonomy, poptions_t options)
{
	int32_t i;

    /*set taxon rank for which result is to be given*/
    for (i = 0; i < taxonomy->ranks->count; i++)
    {
    	if (strcmp(taxonomy->ranks->label[i], options->taxonrank) == 0)
    	{
    		options->taxonrankidx = i;
    		break;
    	}
    }

    if (i == taxonomy->ranks->count)
    {
    	fprintf(stderr,"\nUnknown taxon level: '%s'\n", options->taxonrank);
    	exit(0);
    }
}
/* to get db stats, totals of species, genus etc....*/

static void printAC(pecodnadb_t seqdb,uint32_t seqdbsize)
{
	uint32_t i;

	for (i=0; i< seqdbsize;i++)
		printf("%15s (%8d bp ): %s\n",seqdb[i]->AC,seqdb[i]->SQ_length,seqdb[i]->DE);
}

int main(int argc, char **argv)
{
	pecodnadb_t   seqdb; /* of type ecoseq_t */
	uint32_t        seqdbsize=0;
	ecotaxonomy_t *taxonomy;

	options_t     options;
	int           carg;
	int32_t       errflag=0;

	int32_t       insamples=0;
	int32_t       outsamples=0;
	uint32_t        i;

	pwordcount_t    words;
//	pwordcount_t    words2;
	pprimercount_t  primers;
	ppairtree_t		pairs;

	int32_t		  rankdbstats = 0;

	CNNParams nnparams;

	initoptions(&options);

    while ((carg = getopt(argc, argv, "hAfvcUDSpbE:d:l:L:e:i:r:R:q:3:s:x:t:O:m:a:T:k:M:")) != -1) {

     switch (carg) {
								 /* ---------------------------- */
		 case 'v':               /* set in single strand mode    */
								 /* ---------------------------- */
			options.statistics=TRUE;
			break;

			 /* ---------------------------- */
		 case 'f':               /* set in single strand mode    */
			 /* ---------------------------- */
			 options.filtering=FALSE;
			 break;

			 /* ---------------------------- */
		 case 'A':               /* set in single strand mode    */
			 /* ---------------------------- */
			 options.printAC=TRUE;
			 break;

                                /* -------------------- */
        case 'd':               /* database name        */
                                /* -------------------- */
        	options.prefix = ECOMALLOC(strlen(optarg)+1,
                              "Error on prefix allocation");
           strcpy(options.prefix,optarg);
           break;

                                /* -------------------- */
        case 'h':               /* help                 */
                                /* -------------------- */
        	PrintHelp();
        	 exit(0);
        	 break;

                                /* ------------------------- */
        case 'l':               /* min amplification lenght  */
                                /* ------------------------- */
           sscanf(optarg,"%d",&(options.lmin));
           break;

                                /* -------------------------- */
        case 'L':               /* max amplification lenght   */
                                /* -------------------------- */
          sscanf(optarg,"%d",&(options.lmax));
          break;

                                /* -------------------- */
        case 'e':               /* error max            */
                                /* -------------------- */
          sscanf(optarg,"%d",&(options.error_max));
          break;


								/* ------------------------ */
		case '3':               /* three prime strict match */
                                /* ------------------------ */
			sscanf(optarg,"%d",&(options.strict_three_prime));
			break;

                                /* -------------------- */
        case 'q':               /* strict matching quorum           */
                                /* -------------------- */
          sscanf(optarg,"%f",&(options.strict_quorum));
          break;

								/* -------------------- */
        case 's':               /* strict matching quorum      */
								/* -------------------- */
        	sscanf(optarg,"%f",&(options.sensitivity_quorum));
        	break;

        						/* -------------------- */
        case 't':               /* required taxon level for results           */
        						/* -------------------- */
        	strncpy(options.taxonrank, optarg, 19);
        	options.taxonrank[19] = 0;
            break;

        						/* -------------------- */
        case 'x':               /* strict matching quorum           */
								/* -------------------- */
        	sscanf(optarg,"%f",&(options.false_positive_quorum));
        	break;

							    /* ---------------------------- */
		case 'D':               /* set in double strand mode    */
							    /* ---------------------------- */
			options.doublestrand=1;
			break;

								/* ---------------------------- */
		case 'S':               /* set in single strand mode    */
								/* ---------------------------- */
			options.doublestrand=0;
			break;

			                    /* ---------------------------- */
        case 'U':               /* set in single strand mode    */
			                    /* ---------------------------- */
            options.no_multi_match=TRUE;
            break;

         						/* ------------------------------------------ */
        case 'r':               /* stores the restricting search taxonomic id */
        						/* ------------------------------------------ */
          options.restricted_taxid = ECOREALLOC(options.restricted_taxid,sizeof(int32_t)*(options.r+1),
          									"Error on restricted_taxid reallocation");
          sscanf(optarg,"%d",&(options.restricted_taxid[options.r]));
          options.r++;
          break;

								/* ------------------------------------------ */
        case 'E':               /* stores the restricting search taxonomic id */
								/* ------------------------------------------ */
        	options.exception_taxid = ECOREALLOC(options.exception_taxid,sizeof(int32_t)*(options.e+1),
						            "Error on exception_taxid reallocation");
			sscanf(optarg,"%d",&(options.exception_taxid[options.e]));
			options.e++;
			break;

								/* -------------------- */
        case 'R':               /* reference sequence   */
                                /* -------------------- */
          options.reference = ECOMALLOC(strlen(optarg)+1,
                                     "Error on prefix allocation");
          strcpy(options.reference,optarg);
          break;

         						/* --------------------------------- */
        case 'i':               /* stores the taxonomic id to ignore */
        						/* --------------------------------- */
          options.ignored_taxid = ECOREALLOC(options.ignored_taxid,sizeof(int32_t)*(options.g+1),
          									"Error on excluded_taxid reallocation");
          sscanf(optarg,"%d",&(options.ignored_taxid[options.g]));
          options.g++;
          break;

			                    /* --------------------------------- */
        case 'O':               /* set primer size                   */
        						/* --------------------------------- */
			sscanf(optarg,"%d",&(options.primer_length));
			break;

								/* --------------------------------- */
        case 'm':               /* set salt method                   */
								/* --------------------------------- */
        	sscanf(optarg,"%d",&(options.saltmethod));
        	break;

        						/* --------------------------------- */
        case 'a':               /* set salt 	                     */
								/* --------------------------------- */
        	sscanf(optarg,"%f",&(options.salt));
        	break;

        						/* -------------------- */
        case 'c':               /* sequences are circular */
        						/* --------------------------------- */
		  options.circular = 1;
          break;

								/* -------------------- */
        case 'p':               /* print sets of primers */
								/* --------------------------------- */
        	//options.print_sets_of_primers = TRUE;
        break;
		
        						/* --------------------------------- */
        case 'T':               /* Ignore pairs having specificity below this Threshold     */
								/* --------------------------------- */
        	sscanf(optarg,"%f",&(options.specificity_threshold));
        	break;

								/* --------------------------------- */
		case 'M':               /* Max link percentage for graph     */
								/* --------------------------------- */
			sscanf(optarg,"%f",&(options.max_links_percent));
		break;

								/* --------------------------------- */
		case 'k':               /* links count                   */
								/* --------------------------------- */
			sscanf(optarg,"%d",&(options.links_cnt));
		break;

		case 'b':
			options.filter_on_links = TRUE;
		break;

        case '?':               /* bad option           */
                                /* -------------------- */
            errflag++;
   		}

	}
    options.pnparm = &nnparams;
    if (options.saltmethod != 2) //if not SALT_METHOD_OWCZARZY
    	options.saltmethod = SALT_METHOD_SANTALUCIA; //then force SALT_METHOD_SANTALUCIA

    if (options.salt < 0.01 || options.salt > 0.3) //if salt value out of literature values
    	options.salt = DEF_SALT; //set to default

	nparam_InitParams(&nnparams, DEF_CONC_PRIMERS,DEF_CONC_SEQUENCES,options.salt,options.saltmethod);

    fprintf(stderr,"Reading taxonomy database ...");
    taxonomy = read_taxonomy(options.prefix,0);
    fprintf(stderr,"Ok\n");

    setresulttaxonrank(taxonomy, &options); /*TR: set rank level for statistics*/

    fprintf(stderr,"Reading sequence database ...\n");

    seqdb = readdnadb(options.prefix,taxonomy,&seqdbsize, &options);

	if (options.printAC)
	{
		printAC(seqdb,seqdbsize);
		exit(0);
	}
    if (options.reference)
    	for (i=0; i < seqdbsize;i++)
    		if (strcmp(seqdb[i]->AC,options.reference)==0)
    		{
    			options.refseq=seqdb[i];
    			options.refseqid=i;
    			fprintf(stderr,"Reference sequence %s identified\n",options.reference);
    		}

    fprintf(stderr,"Ok\n");
    fprintf(stderr,"Sequence read : %d\n",(int32_t)seqdbsize);

    updateseqparams(seqdb, seqdbsize, taxonomy, &options, &insamples , &outsamples);
    options.dbsize=seqdbsize;
    options.insamples=insamples;
    options.outsamples=outsamples;

    rankdbstats = getrankdbstats(seqdb, seqdbsize, taxonomy, &options);

    fprintf(stderr,"Database is constituted of %5d examples        corresponding to %5d %s\n",insamples,
    		options.intaxa,options.taxonrank);
    fprintf(stderr,"                       and %5d counterexamples corresponding to %5d %s\n",outsamples,
    		options.outtaxa,options.taxonrank);
    fprintf(stderr,"Total distinct %s count %d\n",options.taxonrank, rankdbstats);

    fprintf(stderr,"\nIndexing words in sequences\n");

    words = lookforStrictPrimer(seqdb,seqdbsize,insamples,&options);
    fprintf(stderr,"\n  Strict primer count : %d\n",words->size);
    
    /*/TR Testing
    fprintf(stderr,"\nReducing for debugging\n");
    words = reduce_words_to_debug (words, &options);
    ///*/
//    options.filtering=FALSE;
//    words2= lookforStrictPrimer(seqdb,seqdbsize,insamples,&options);
//    fprintf(stderr,"\n  Strict primer count : %d\n",words2->size);
//
//    fprintf(stderr,"\n\n  Primer sample : \n");
//    for (i=0; i<words->size; i++)
//    	fprintf(stderr,"  + Primer : %s   sequence count : %d\n",ecoUnhashWord(words->words[i],options.primer_length),words->strictcount[i]);
//    fprintf(stderr,"\n\n  Primer sample : \n");
//    for (i=0; i<words2->size; i++)
//    	fprintf(stderr,"  + Primer : %s   sequence count : %d\n",ecoUnhashWord(words2->words[i],options.primer_length),words2->strictcount[i]);

    if (options.no_multi_match)
    {
    	(void)filterMultiStrictPrimer(words);
    	fprintf(stderr,"\n  Strict primer with single match count : %d\n",words->size);
    }


    fprintf(stderr,"\n\n  Primer sample : \n");
    for (i=0; i<MINI(10,words->size); i++)
    	fprintf(stderr,"  + Primer : %s   sequence count : %d\n",ecoUnhashWord(words->words[i],options.primer_length),words->strictcount[i]);

    fprintf(stderr,"\nEncoding sequences for fuzzy pattern matching...\n");
    for (i=0;i<seqdbsize;i++)
    {
    	encodeSequence(seqdb[i]);
    	fprintf(stderr,"  Encoded sequences %5d/%5d          \r",(int32_t)i+1,(int32_t)seqdbsize);
    }

    ECOFREE(words->strictcount,"Free strict primer count table");

    if (options.error_max == 0)//aho, if(options.error_max == 0 && 0) old
    	primers = ahoc_lookforStrictPrimers (seqdb,seqdbsize,insamples,words,&options);
    else
    	primers = lookforAproxPrimer(seqdb,seqdbsize,insamples,words,&options);
    
    //for (i=0; i<primers->size; i++)
    //	print_wordwith_positions (primers->primers[i], seqdbsize, &options);

    ECOFREE(words->words,"Free strict primer table");
    ECOFREE(words,"Free strict primer structure");
    fprintf(stderr,"\n\n  Approximate repeats :%d \n", primers->size);

    fprintf(stderr,"\n\n  Primer sample : \n");
    for (i=0; i<MINI(10,primers->size); i++)
    	fprintf(stderr,"  + Primer : %s   example sequence count : %5d counterexample sequence count : %5d  status : %s\n",ecoUnhashWord(primers->primers[i].word,options.primer_length),
                       primers->primers[i].inexample,
                       primers->primers[i].outexample,
                       primers->primers[i].good ? "good":"bad");

    fprintf(stderr,"\n");

   
    	pairs = buildPrimerPairs(seqdb, seqdbsize, primers, &options);
    	printpairs (pairs, &options,taxonomy, seqdb);
    
    return 0;
}

#define DEBUG_WORDS_CNT 14
pwordcount_t reduce_words_to_debug (pwordcount_t words, poptions_t options)
{
	uint32_t i, k;
	pwordcount_t new_words;
	char *rwrd;
	char dwrd[20];
	/*char *strict_words[DEBUG_WORDS_CNT] = {"GAGTCTCTGCACCTATCC", "GCAATCCTGAGCCAAATC", "ACCCCTAACCACAACTCA", 
			"TCCGAACCGACTGATGTT", "GAAGCTTGGGTGAAACTA", "GGAGAACCAGCTAGCTCT", "GCTGGTTCTCCCCGAAAT",
			"TCGATTTGGTACCGCTCT", "AAAGGAGAGAGAGGGATT", "GGATTGCTAATCCGTTGT", "CCCCCATCGTCTCACTGG",
			"TGAGGCGCAGCAGTTGAC", "GCGCTACGGCGCTGAAGT", "TTTCCTGGGAGTATGGCA"};*/
	char *strict_words[DEBUG_WORDS_CNT] = {"CTCCGGTCTGAACTCAGA", "TGTTGGATCAGGACATCC", "TAGATAGAAACCGACCTG", 
			"TGGTGCAGCCGCTATTAA", "AGATAGAAACTGACCTGG", "TGGTGCAGCCGCTATTAA", "CTAATGGTGCAGCCGCTA",
			"TAGAAACTGACCTGGATT", "AGATAGAAACCGACCTGG", "ATGGTGCAGCCGCTATTA", "ATAGATAGAAACCGACCT",
			"GCCGCTATTAAGGGTTCG", "GGTGCAGCCGCTATTAAG", "TAGAAACTGACCTGGATT"};
	int word_seen[DEBUG_WORDS_CNT];
	                   
	
	new_words = ECOMALLOC(sizeof(wordcount_t),"Cannot allocate memory for word count structure");
	new_words->inseqcount = words->inseqcount;
	new_words->outseqcount = words->outseqcount;
	new_words->size = DEBUG_WORDS_CNT;
	new_words->strictcount = ECOMALLOC((new_words->size*sizeof(uint32_t)), "Cannot allocate memory for word count table");
	new_words->words = ECOMALLOC(new_words->size*sizeof(word_t), "I cannot allocate memory for debug words");
			
	for (k = 0; k < DEBUG_WORDS_CNT; k++)
		word_seen[k] = 0;
	
	for (i=0; i < words->size; i++)
	{
		rwrd = ecoUnhashWord(words->words[i],options->primer_length);
		strcpy (dwrd, rwrd);
		rwrd = ecoUnhashWord(ecoComplementWord(words->words[i],options->primer_length),options->primer_length);
		for (k = 0; k < DEBUG_WORDS_CNT; k++)
		{
			if (strcmp (dwrd, strict_words[k]) == 0) break;
			if (strcmp (rwrd, strict_words[k]) == 0) break;
		}
		
		if (k < DEBUG_WORDS_CNT)
		{
			if (word_seen[k] == 0)
			{
				new_words->words[k] = words->words[i];
				new_words->strictcount[k] = words->strictcount[i];
			}
			word_seen[k]++;
		}
	}
	
	fprintf (stderr, "Debug Words Info:\n");
	for (k = 0; k < DEBUG_WORDS_CNT; k++)
		fprintf (stderr, "%s:%d\n", strict_words[k], word_seen[k]);
	
	
	//clean input wods;
	ECOFREE(words->words,"Clean  word table");
	ECOFREE(words->strictcount,"Clean  word count table");
	ECOFREE(words,"Clean word structure");
	
	return new_words;
}

void print_wordwith_positions (primer_t prm, uint32_t seqdbsize, poptions_t options)
{
	char *wrd;
	uint32_t i, j;
	char *twrd = "GCCTGTTTACCAAAAACA";

	wrd = ecoUnhashWord(prm.word,options->primer_length);

	if (strcmp (twrd, wrd) == 0)
	{
		printf ("Positions for Word: %s\n", wrd);
		for (i=0; i<seqdbsize; i++)
		{
			if (prm.directCount[i] > 0)
			{
				printf ("%d:", i);
				if (prm.directCount[i] == 1)
					printf ("%d", prm.directPos[i].value);
				else
					for (j=0; j<prm.directCount[i]; j++)
						printf ("%d,", prm.directPos[i].pointer[j]);
				printf (" ");
			}
		}
		printf ("\n");
		for (i=0; i<seqdbsize; i++)
		{
			if (prm.reverseCount[i] > 0)
			{
				printf ("%d:", i);
				if (prm.reverseCount[i] == 1)
					printf ("%d", prm.reversePos[i].value);
				else
					for (j=0; j<prm.reverseCount[i]; j++)
						printf ("%d,", prm.reversePos[i].pointer[j]);
				printf (" ");
			}
		}
		printf ("\n");
	}
}
