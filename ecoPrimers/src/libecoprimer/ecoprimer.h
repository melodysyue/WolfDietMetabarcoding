/*
 * epsort.h
 *
 *  Created on: 6 nov. 2008
 *      Author: coissac
 */

#ifndef EPSORT_H_
#define EPSORT_H_

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include "ecotype.h"
#include "../libecoPCR/ecoPCR.h"
#include "../libthermo/nnparams.h"
#include "apat.h"

#define DEBUG
#include "debug.h"

/****
 * Word format used :
 *
 *     bit 63    : bad word -> this word should not be used
 *     bit 62    : multi word -> this word is not uniq in at least one seq
 *     bits 0-61 : hashed dna word of max size 31 pb
 *                 code used for a : 00
 *                 code used for c : 01
 *                 code used for g : 10
 *                 code used for t : 11
 */

typedef uint64_t word_t, *pword_t;

#define WORD(x)           ((x) & 0x3FFFFFFFFFFFFFFFLLU)
#define WORD(x)           ((x) & 0x3FFFFFFFFFFFFFFFLLU)

#define ISBADWORD(x)      (((x) & 0x8000000000000000LLU) >> 63)
#define SETBADWORD(x)     ((x)  | 0x8000000000000000LLU)
#define RESETBADWORD(x)   ((x)  & 0x7FFFFFFFFFFFFFFFLLU)

#define ISMULTIWORD(x)    (((x) & 0x4000000000000000LLU) >> 62)
#define SETMULTIWORD(x)   ((x)  | 0x4000000000000000LLU)
#define RESETMULTIWORD(x) ((x)  & 0xBFFFFFFFFFFFFFFFLLU)


#define WORDMASK(s)       ((1LLU << ((s) * 2)) -1)
#define LSHIFTWORD(x,s)    (((x) << 2) & WORDMASK(s))
#define RSHIFTWORD(x,s)    (((x) & WORDMASK(s))>> 2)
#define ERRORMASK(s)       ((int32_t)((1LLU << (s)) -1))

#define RAPPENDBASE(x,s,c) (LSHIFTWORD((x),(s)) | (word_t)(c))
#define LAPPENDBASE(x,s,c) (RSHIFTWORD((x),(s)) | ((word_t)((~(c)) & 3) << (((s)-1) *2)))


#define ECO_ASSERT(x,message) if (!(x)) \
	                           { \
	                              fprintf(stderr,"Assertion Error in %s (line %d): %s\n", \
	                                          __FILE__,\
	                                          __LINE__,\
	                                          message\
	                                     ); \
	                             exit(ECO_ASSERT_ERROR); \
	                           }

#define MINI(x,y) (((x) < (y)) ? (x):(y))
#define MAXI(x,y) (((x) < (y)) ? (y):(x))

#define FWORDSIZE (13)
#define FWORDMASK WORDMASK(FWORDSIZE)
#define FILTERWORD(x) ((uint32_t)((x) & FWORDMASK))
#define CFILTERWORD(x,s) ((uint32_t)(((x) >> (((s)-FWORDSIZE)*2)) & FWORDMASK))



typedef struct {
	pword_t    words;
	uint32_t    *strictcount;
	uint32_t     inseqcount;
	uint32_t     outseqcount;
	uint64_t     size;
} wordcount_t, *pwordcount_t;


typedef union {
	uint32_t *pointer;
	uint32_t value;
} poslist_t, *pposlist_t;


/**
 * primer_t structure store fuzzy match positions for a primer
 * on all sequences
 */

typedef struct {
    word_t     word;              //< code for the primer
	uint32_t   *directCount;      //< Occurrence count on direct strand
	pposlist_t directPos;        //< list of position list on direct strand

	uint32_t   *reverseCount;     //< Occurrence count on reverse strand
	pposlist_t reversePos;       //< list of position list on reverse strand

	bool_t     good;              //< primer match more than quorum example and no
	                              //  more counterexample quorum.

	uint32_t   inexample;         //< count of example sequences matching primer
	uint32_t   outexample;        //< count of counterexample sequences matching primer
} primer_t, *pprimer_t;

/**
 * primercount_t structure store fuzzy match positions for all primers
 * on all sequences as a list of primer_t
 */
typedef struct {
	pprimer_t   primers;
	uint32_t    size;
} primercount_t, *pprimercount_t;

typedef struct {
	pprimer_t   primer;
	uint32_t    position;
	bool_t      strand;
} primermatch_t, *pprimermatch_t;

/*TR: Added*/
typedef struct {
	pprimermatch_t matches;
	uint32_t	matchcount;
} primermatchcount_t, *pprimermatchcount_t;

typedef struct {
	pecoseq_t  sequence;
	bool_t     strand;
	const char *amplifia;
	int32_t	   length;
	uint32_t   begin;
	uint32_t   end;
} amplifia_t, *pamplifia_t;

typedef struct {
	pamplifia_t amplifias;
	uint32_t    ampcount;
	uint32_t	ampslot;
} amplifiacount_t, *pamplifiacount_t;

typedef struct {
	char *amplifia;
	int32_t *taxonids;
	uint32_t seqidcount;
	uint32_t seqidindex;
} ampseqset_t, *pampseqset_t;

typedef struct {
	int32_t taxonid;
	char **amplifia;
	uint32_t amplifiacount;
	uint32_t amplifiaindex;
} taxampset_t, *ptaxampset_t;

typedef struct {
	pprimer_t 		p1;
	bool_t			asdirect1;
	pprimer_t 		p2;
	bool_t			asdirect2;

	amplifiacount_t pcr;

	uint32_t   		inexample;         //< example sequence count
	uint32_t   		outexample;        //< counterexample sequence count
	uint32_t   		intaxa;            //< example taxa count
	uint32_t   		outtaxa;            //< counterexample taxa count
	uint32_t		notwellidentifiedtaxa;

	int				*wellIdentifiedSeqs; //< an array having elements equla to total seqs
									 // values are either 0 or 1, if seq is well identified
									 // its 1 else 0
	int				*coveredSeqs; //< an array having elements equal to total seqs, 1 if seq is covered else 0

					// these statistics are relative to inexample sequences

	uint32_t 		mind;				//< minimum distance between primers
	uint32_t		maxd;				//< maximum distance between primers
	uint32_t        sumd;				//< distance sum
	uint32_t        amplifiacount;
	float			yule;
	float			quorumin;
	float           quorumout;
	float           bs;
	float           bc;
	int32_t        refsequence;
//
//	uint32_t		taxsetcount;
//	uint32_t		taxsetindex;
//	ptaxampset_t	taxset;
//
//	uint32_t		oktaxoncount;
	uint32_t		curseqid;
	float			p1temp;			//strict primer1 melting temperature
	float			p1mintemp;		//approx primer1 minimum melting temperature
	float			p2temp;			//strict primer2 melting temperature
	float			p2mintemp;		//approx primer2 minimum melting temperature
} pair_t, *ppair_t;

/*TR: Added*/

typedef struct {
	size_t  	paircount;
	size_t      pairslots;
    void*       next;
	pair_t 	    pairs[1];
} pairlist_t, *ppairlist_t;

typedef struct {
	ppairlist_t first;
	ppairlist_t last;
	void       *tree;
	int32_t     count;
} pairtree_t, *ppairtree_t;

typedef struct {
	pword_t     words;
	uint32_t   *count;
	uint32_t    push;
	uint32_t    pop;
	uint32_t    size;
	bool_t      empty;
	bool_t      full;
} queue_t, *pqueue_t;

typedef struct {
	pword_t     words;
	uint32_t   *count;
	uint32_t    write;
	uint32_t    read1;
	uint32_t    read2;
	uint32_t    size;
} merge_t, *pmerge_t;

typedef struct {
	const char 		*amplifia;
	bool_t     	strand;
	int32_t	   	length;
	int32_t		taxoncount;
	void		*taxontree;
}amptotaxon_t, *pamptotaxon_t;

typedef struct {
	int32_t 	taxid;
	void		*amptree;
}taxontoamp_t, *ptaxontoamp_t;

typedef struct {
	bool_t         printAC;
	bool_t         statistics;
	bool_t         filtering;
	uint32_t        lmin;                   //**< Amplifia minimal length
	uint32_t        lmax;                   //**< Amplifia maximal length
	uint32_t        error_max;              //**< maximum error count in fuzzy search
	uint32_t        primer_length;          //**< minimal length of the primers
	int32_t		  *restricted_taxid;      //**< limit amplification below these taxid
	int32_t       *ignored_taxid;         //**< no amplification below these taxid
	int32_t		  *exception_taxid;
	char          *prefix;
	char          *reference;
	pecoseq_t      refseq;
	uint32_t       refseqid;
	uint32_t		  circular;
	uint32_t        doublestrand;
	float         strict_quorum;
	float         strict_exclude_quorum;
	float         sensitivity_quorum;
	float         false_positive_quorum;
	uint32_t      strict_three_prime;
	int32_t		  r;                      //**< count of restrited taxa (restricted_taxid array size)
	int32_t		  g;                      //**< count of ignored taxa (ignored_taxid array size)
	int32_t		  e;                      //**< count of ignored taxa (ignored_taxid array size)
	bool_t        no_multi_match;
	char		  taxonrank[20];  //TR to count ranks against a pair
	int32_t		  taxonrankidx;   //TR to count ranks against a pair

	// Some statistics useful for options filters

	int32_t       dbsize;
	int32_t		  insamples;
	int32_t		  outsamples;
	int32_t       intaxa;
	int32_t       outtaxa;
	int			  saltmethod;
	float		  salt;
	PNNParams	  pnparm;
	bool_t		  print_sets_of_primers;
	float		  specificity_threshold;
	int			  links_cnt;
	float		  max_links_percent;
	bool_t		  filter_on_links;
} options_t, *poptions_t;

typedef ecoseq_t  **pecodnadb_t;

void sortword(pword_t table,uint32_t N);


pecodnadb_t readdnadb(const char *name, ecotaxonomy_t *taxonomy, uint32_t *size,poptions_t options);

int isGoodTaxon(ecotaxonomy_t *taxonomy,int32_t taxon,poptions_t options);
int isExampleTaxon(ecotaxonomy_t *taxonomy,int32_t taxon,poptions_t options);
int isCounterExampleTaxon(ecotaxonomy_t *taxonomy,int32_t taxon,poptions_t options);

uint32_t ecoWordCount(uint32_t wordsize, uint32_t circular, ecoseq_t *seq);
pword_t ecoHashSequence(pword_t dest, uint32_t wordsize, uint32_t circular, uint32_t doublestrand, ecoseq_t *seq,uint32_t *size,int32_t *neededWords,uint32_t neededWordCount,
	    int32_t quorum);
uint32_t ecoCompactHashSequence(pword_t dest,uint32_t size);
const char* ecoUnhashWord(word_t word,uint32_t size);
word_t ecoComplementWord(word_t word,uint32_t size);
uint32_t ecoFindWord(pwordcount_t table,word_t word);


void ecomerge(pwordcount_t data,uint32_t s1,uint32_t s2,uint32_t remainingSeq,uint32_t seqQuorum);
pwordcount_t initCountTable(pwordcount_t table, uint32_t wordsize, uint32_t circular, uint32_t doublestrand,uint32_t seqQuorum,ecoseq_t *seq,int32_t *neededWords,uint32_t neededWordCount);
void addSeqToWordCountTable(pwordcount_t table, uint32_t wordsize, uint32_t circular, uint32_t doublestrand,uint32_t exampleCount,uint32_t seqQuorum,ecoseq_t *seq,int32_t *neededWords,uint32_t neededWordCount);

pqueue_t newQueue(pqueue_t queue, uint32_t size);
pqueue_t resizeQueue(pqueue_t queue, uint32_t size);

void pop(pqueue_t queue);
void push(pqueue_t queue, word_t word, uint32_t count);

pqueue_t cleanQueue(pqueue_t queue);

pwordcount_t lookforStrictPrimer(pecodnadb_t database, uint32_t seqdbsize,
		                         uint32_t exampleCount, poptions_t options);
uint32_t filterMultiStrictPrimer(pwordcount_t strictprimers);

void encodeSequence(ecoseq_t *seq);
ppattern_t buildPatternFromWord(word_t word, uint32_t patlen);

pprimercount_t lookforAproxPrimer(pecodnadb_t database, uint32_t seqdbsize,uint32_t exampleCount,
		                     pwordcount_t words,poptions_t options);

void sortmatch(pprimermatch_t table,uint32_t N);

ppairtree_t initpairtree(ppairtree_t tree);
ppair_t pairintree (pair_t key,ppairtree_t pairlist);
ppair_t insertpair(pair_t key,ppairtree_t list);


/*TR: Added*/
ppairtree_t buildPrimerPairs(pecodnadb_t seqdb,uint32_t seqdbsize,pprimercount_t primers,poptions_t options);

int32_t counttaxon(int32_t taxid);
int32_t getrankdbstats(pecodnadb_t seqdb,
					   uint32_t seqdbsize,
					   ecotaxonomy_t *taxonomy,
					   poptions_t options);
float taxonomycoverage(ppair_t pair, poptions_t options, pecodnadb_t seqdb,uint32_t seqdbsize);
char ecoComplementChar(char base);
void taxonomyspecificity (ppair_t pair, pecodnadb_t seqdb,uint32_t seqdbsize);

int32_t *filteringSeq(pecodnadb_t database, uint32_t seqdbsize,
					 uint32_t exampleCount,poptions_t options,uint32_t *size,int32_t  sequenceQuorum);

void printSeqTest(pecodnadb_t seqdb,uint32_t seqdbsize);

#endif /* EPSORT_H_ */
