/*
 * aproxpattern.c
 *
 *  Created on: 20 nov. 2008
 *      Author: coissac
 */


#include  "ecoprimer.h"
#include  "apat.h"
#include  <math.h>

static uint8_t encoder[] = {0,                                            // A
		                   4,                                           // b
		                   1,                                            // C
	       	               4,4,4,                                     // d, e, f
		                   2,                                            // G
		                   4,4,4,4,4,4,4,4,4,4,4,4,          // h,i,j,k,l,m,n,o,p,q,r,s
		                   3,3,                                           // T,U
		                   4,4,4,4,4};                              // v,w,x,y,z


ppattern_t buildPatternFromWord(word_t word, uint32_t patlen)
{
	static pattern_t       pattern;
	uint32_t i;

    for (i = 0 ; i < ALPHA_LEN ; i++)
       pattern[i] = 0x0;

    for (i=0;i < patlen; i++)
	{
		pattern[word & 3LLU] |= 1 << i;
		word>>=2;
	}

	return pattern;

}


#ifdef IS_UPPER
#undef IS_UPPER
#endif

/* -------------------------------------------- */
/* encode sequence                              */
/* IS_UPPER is slightly faster than isupper     */
/* -------------------------------------------- */

#define IS_UPPER(c) (((c) >= 'A') && ((c) <= 'Z'))

void encodeSequence(ecoseq_t *seq)
{
     int   i;
     uint8_t *data;
     char  *cseq;

     data = (uint8_t*)(seq->SQ);
     cseq = seq->SQ;

     for (i=0;i<seq->SQ_length;i++,data++,cseq++)
     {
         *data = encoder[(IS_UPPER(*cseq) ? *cseq : 'Z') - 'A'];
     }
}

pprimercount_t lookforAproxPrimer(pecodnadb_t database, uint32_t seqdbsize,uint32_t exampleCount,
		                     pwordcount_t words,poptions_t options)
{
	pprimer_t  data;
	pprimercount_t primers;
	ppattern_t pattern;
	patternParam_t params;
	uint32_t   i;
	uint32_t   w;
	uint32_t   j;
	Stacki     positions;
	uint32_t   count=1;
	uint32_t   goodPrimers=0;

	uint32_t  inSequenceQuorum;
	uint32_t  outSequenceQuorum;
	bool_t    conserved = TRUE;

	//poslist_t ttt;


	inSequenceQuorum = (uint32_t)floor((float)exampleCount * options->sensitivity_quorum);
	outSequenceQuorum = (uint32_t)floor((float)(seqdbsize-exampleCount) * options->false_positive_quorum);

	fprintf(stderr,"  Primers should be at least present in %d/%d example sequences\n",inSequenceQuorum,exampleCount);
	fprintf(stderr,"  Primers should not be present in more than %d/%d counterexample sequences\n",outSequenceQuorum,(seqdbsize-exampleCount));

	data = ECOMALLOC(words->size * sizeof(primer_t),
			         "Cannot allocate memory for fuzzy matching results");

	params.circular = options->circular;
	params.maxerr   = options->error_max;
//	params.omask    = (1 << options->strict_three_prime) -1;
	params.omask    = 0;
	params.patlen   = options->primer_length;

	positions.val=NULL;

	for (i=0,w=0; i < words->size; i++)
	{
		data[w].word=WORD(words->words[i]);
		data[w].inexample = 0;
		data[w].outexample= 0;
		count = 1;

		if (conserved)
		{
			data[w].directCount=ECOMALLOC(seqdbsize * sizeof(uint32_t),
                                      "Cannot allocate memory for primer position");
			data[w].directPos = ECOMALLOC(seqdbsize * sizeof(poslist_t),
				                      "Cannot allocate memory for primer position");
			data[w].reverseCount=ECOMALLOC(seqdbsize * sizeof(uint32_t),
	                                      "Cannot allocate memory for primer position");
			data[w].reversePos = ECOMALLOC(seqdbsize * sizeof(poslist_t),
					                       "Cannot allocate memory for primer position");
		}

		pattern = buildPatternFromWord(data[w].word,options->primer_length);
		positions.val=NULL;

		for (j=0; j < seqdbsize && (count < 2 || !options->no_multi_match); j++)
		{
			positions.cursor=0;
			positions.top =0;
			if (!positions.val)
			{
				positions.size=1;
				positions.val = ECOMALLOC(sizeof(uint32_t),
                                      "Cannot allocate memory for primer position");
			}


			count = ManberAll(database[j],pattern,&params,&positions);
			data[w].directCount[j]=count;


			if (count>1)
			{
				data[w].directPos[j].pointer = (uint32_t*)positions.val;
				positions.val=NULL;
			}
			else
			{
				data[w].directPos[j].pointer=NULL;
				if (count==1)
					data[w].directPos[j].value = (uint32_t)*(positions.val);
			}


		}

		pattern = buildPatternFromWord(ecoComplementWord(data[w].word,options->primer_length),
				                       options->primer_length);

		for (j=0; j < seqdbsize && (count < 2 || !options->no_multi_match); j++)
		{
			positions.cursor=0;
			positions.top =0;
			if (!positions.val)
			{
				positions.size=1;
				positions.val = ECOMALLOC(sizeof(uint32_t),
                                      "Cannot allocate memory for primer position");
			}

			count = ManberAll(database[j],pattern,&params,&positions);
			data[w].reverseCount[j]=count;

			if (count>1)
			{
				data[w].reversePos[j].pointer = (uint32_t*)positions.val;
				positions.val=NULL;
			}
			else
			{
				data[w].reversePos[j].pointer=NULL;
				if (count==1)
					data[w].reversePos[j].value = (uint32_t)*(positions.val);
			}

			if (database[j]->isexample)
			{
				data[w].inexample+=(data[w].directCount[j] || data[w].reverseCount[j])? 1:0;
			}
			else
			{
				data[w].outexample+=(data[w].directCount[j] || data[w].reverseCount[j])? 1:0;

			}

            count+=data[w].directCount[j];
		}

		data[w].good = data[w].inexample >= inSequenceQuorum && data[w].outexample <= outSequenceQuorum;
		goodPrimers+=data[w].good? 1:0;

		fprintf(stderr,"Primers %5d/%d analyzed  => sequence : %s in %d example and %d counterexample sequences       \r",
				i+1,words->size,ecoUnhashWord(data[w].word,options->primer_length),
				data[w].inexample,data[w].outexample);


		conserved=data[w].inexample >= inSequenceQuorum;
		conserved=conserved && (count < 2 || !options->no_multi_match);

		if (conserved)
			w++;
	}

	if (positions.val)
		ECOFREE(positions.val,"Free stack position pointer");

	if (!conserved)
	{
		ECOFREE(data[w].directCount,"Free direct count table");
		ECOFREE(data[w].directPos,"Free direct count table");
		ECOFREE(data[w].reverseCount,"Free direct count table");
		ECOFREE(data[w].reversePos,"Free direct count table");
	}

	fprintf(stderr,"\n\nOn %d analyzed primers %d respect quorum conditions\n",words->size,goodPrimers);
	fprintf(stderr,"Conserved primers for further analysis : %d/%d\n",w,words->size);

	primers = ECOMALLOC(sizeof(primercount_t),"Cannot allocate memory for primer table");
	primers->primers=ECOREALLOC(data,
			                    w * sizeof(primer_t),
								"Cannot reallocate memory for fuzzy matching results");
	primers->size=w;

	return primers;
}
