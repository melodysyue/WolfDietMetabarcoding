/*
 * strictprimers.c
 *
 *  Created on: 7 nov. 2008
 *      Author: coissac
 */

#define _GNU_SOURCE
#include "ecoprimer.h"
#include <string.h>
#include <math.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

#ifndef RUSAGE_SELF
#define   RUSAGE_SELF     0
#define   RUSAGE_CHILDREN     -1
#endif

static double timeval_subtract (struct timeval *x, struct timeval *y);


 /* Subtract the `struct timeval' values X and Y,
	Return elapsed secondes as a double.  */

double timeval_subtract (struct timeval *x, struct timeval *y)
{
   struct timeval result;

   /* Perform the carry for the later subtraction by updating y. */
   if (x->tv_usec < y->tv_usec) {
	 int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
	 y->tv_usec -= 1000000 * nsec;
	 y->tv_sec += nsec;
   }
   if (x->tv_usec - y->tv_usec > 1000000) {
	 int nsec = (x->tv_usec - y->tv_usec) / 1000000;
	 y->tv_usec += 1000000 * nsec;
	 y->tv_sec -= nsec;
   }

   /* Compute the time remaining to wait.
	  tv_usec is certainly positive. */
   result.tv_sec = x->tv_sec - y->tv_sec;
   result.tv_usec = x->tv_usec - y->tv_usec;

   return (double)result.tv_sec + (double)result.tv_usec/1e6;
 }

pwordcount_t initCountTable(pwordcount_t table, uint32_t wordsize, uint32_t circular, uint32_t doublestrand,uint32_t seqQuorum,ecoseq_t *seq,int32_t *neededWords,uint32_t neededWordCount)
{
	uint32_t i;
	uint32_t buffsize;
	//wordcount_t t;

	if (!table)
		table = ECOMALLOC(sizeof(wordcount_t),"Cannot allocate memory for word count structure");

	table->words=NULL;
	table->size =0;
	table->outseqcount=0;
	table->inseqcount=0;
	table->strictcount =0;

	if (seq)
	{
		table->words = ecoHashSequence(NULL,wordsize,circular,doublestrand,seq,&buffsize,neededWords,neededWordCount,seqQuorum);
		table->size  = ecoCompactHashSequence(table->words,buffsize);

		table->inseqcount=1;
		table->strictcount =ECOMALLOC((table->size*sizeof(uint32_t)),
							  "Cannot allocate memory for word count table"
							 );

		for (i=0; i < table->size; i++) table->strictcount[i]=1;
	}

	return table;
}

void addSeqToWordCountTable(pwordcount_t table, uint32_t wordsize, uint32_t circular, uint32_t doublestrand,uint32_t exampleCount,uint32_t seqQuorum,ecoseq_t *seq,int32_t *neededWords,uint32_t neededWordCount)
{
	uint32_t buffersize;
	pword_t newtable;
	uint32_t  newsize;
	uint32_t  i;

	buffersize = table->size + ecoWordCount(wordsize,circular,seq);

	table->words = ECOREALLOC(table->words,buffersize*sizeof(word_t),
							  "\n\nCannot allocate memory to extend word table" );

	/*
	 * newtable is a pointer on the memory planed to be used for the new sequence (ecoWordCount new hash codes max)
	 */
	newtable = table->words + table->size;

//	DEBUG_LOG("Words = %x (%u) new = %x", table->words,table->size,newtable);

	(void)ecoHashSequence(newtable,wordsize,circular,doublestrand,seq,&newsize,neededWords,neededWordCount,seqQuorum);
//	DEBUG_LOG("new seq wordCount : %d",newsize);

	/*
	 * at this stage, new hash codes have been added in the table but the table is not sorted
	 */

	newsize = ecoCompactHashSequence(newtable,newsize);

	/*
	 * new hash codes have now been sorted BUT the whole table is not.
	 * MULTIWORDS have been tagged (and compacted)
	 */

//	DEBUG_LOG("compacted wordCount : %d",newsize);
	buffersize = table->size + newsize;

	/*
	 * buffersize is now set to the REAL size used by the table (but the memory chunck may be larger)
	 */

	// resize the count buffer

	table->inseqcount++;

	//fprintf (stderr, "\nOldAddress: %x", table->strictcount);
	table->strictcount = ECOREALLOC(table->strictcount,(buffersize+5000)*sizeof(uint32_t),
								"Cannot allocate memory to extend example word count table");
	//fprintf (stderr, " NewAddress: %x\n", table->strictcount);
	

	for (i=table->size; i < buffersize; i++)
		table->strictcount[i]=1;

	/*
	 * new words in the table are set to a count of ONE
	 */

	// Now we have to merge in situ the two tables

	ecomerge(table,table->size,newsize,exampleCount - table->inseqcount,seqQuorum);
//	DEBUG_LOG("Dictionnary size : %d",table->size);

}

pwordcount_t lookforStrictPrimer(pecodnadb_t database, uint32_t seqdbsize,
								 uint32_t exampleCount,poptions_t options)
{
	struct rusage start;
	struct rusage usage;
	double seconde;
	char   *logfilename;
	FILE   *logfile;
	uint32_t i, j;
	bool_t first=TRUE;
	pwordcount_t strictprimers=NULL;
	uint64_t totallength=0;
	uint32_t  sequenceQuorum = (uint32_t)floor((float)exampleCount * options->strict_quorum);
	int32_t *neededWords;
	uint32_t neededWordCount;

	fprintf(stderr,"Filtering... ");

	if (options->filtering)
		neededWords = filteringSeq(database,seqdbsize,exampleCount,options,&neededWordCount,(int32_t)sequenceQuorum);
	else
	{
		neededWordCount=0;
		neededWords=NULL;
	}

	if (options->statistics)
	{
		asprintf(&logfilename,"ecoprimer_%d.log",getpid());
		logfile = fopen(logfilename,"w");
		fprintf(logfile,"# seq\tlength\tsize\ttime\tspeed\n");
		fclose(logfile);
	}


	fprintf(stderr,"  Primers should be at least present in %d/%d example sequences\n",sequenceQuorum,exampleCount);

	strictprimers = initCountTable(NULL,options->primer_length,
                                                 options->circular,
                                                 options->doublestrand,
                                                 0,
                                   NULL,NULL,0);


	getrusage(RUSAGE_SELF,&start);

    for (i=0;i<seqdbsize;i++)
    {
    	if (database[i]->isexample && database[i]->SQ_length > options->primer_length)
    	{

    		if (first)
    		{
    			strictprimers = initCountTable(strictprimers,options->primer_length,
                                                             options->circular,
                                                             options->doublestrand,
                                                             sequenceQuorum,
                                               database[i],neededWords,neededWordCount);
    			first=FALSE;
    		}
    		else
    		{
    			uint32_t s;
    			s = strictprimers->size;
//    			DEBUG_LOG("stack size : %u",s);
    			addSeqToWordCountTable(strictprimers,options->primer_length,
    					                             options->circular,
    					                             options->doublestrand,
    					                             exampleCount,
    					                             sequenceQuorum,
    					               database[i],neededWords,neededWordCount);
    		};
    		totallength+=database[i]->SQ_length;
    		getrusage(RUSAGE_SELF,&usage);
    		if (options->statistics)
    		{
    			asprintf(&logfilename,"ecoprimer_%d.log",getpid());
    			logfile = fopen(logfilename,"a");
    			seconde = timeval_subtract(&(usage.ru_utime),&(start.ru_utime)) +
    			          timeval_subtract(&(usage.ru_stime),&(start.ru_stime));
    			fprintf(logfile,"%d\t%llu\t%lu\t%8.3f\t%8.3e\n",i,
    					(long long unsigned)totallength,
					strictprimers->size*(sizeof(int64_t)+sizeof(int32_t)),
    					seconde,seconde/(double)totallength);
    			fclose(logfile);
    		}
    	}
    	else
    		strictprimers->outseqcount++;

    	fprintf(stderr,"  Indexed sequences %5d/%5d  : considered words %-10llu  \r",
				(int32_t)i+1,(int32_t)seqdbsize,
				(long long unsigned)strictprimers->size);

//    	DEBUG_LOG("First  word : %s ==> %d",ecoUnhashWord(strictprimers->words[0],18),strictprimers->incount[0])
//    	DEBUG_LOG("Second word : %s ==> %d",ecoUnhashWord(strictprimers->words[1],18),strictprimers->incount[1])
    }

    strictprimers->strictcount = ECOREALLOC(strictprimers->strictcount,
    		                            sizeof(uint32_t)*strictprimers->size,
    		                            "Cannot reallocate strict primer count table");
    strictprimers->words   = ECOREALLOC(strictprimers->words,
										sizeof(word_t)*strictprimers->size,
										"Cannot reallocate strict primer table");

    if (neededWords)
    	ECOFREE(neededWords,"Clean needed word table");

	//TR: Somehow for some primers strictcount value is extremely large hence invalid
	//we need to remove these primers from the list
	j = strictprimers->size+1;
	for (i=0; i<strictprimers->size; i++)
	{
		if (strictprimers->strictcount[i] > seqdbsize)
		{
			if (j == (strictprimers->size+1))
				j = i;
		}
		
		if (j < i && strictprimers->strictcount[i] <= seqdbsize)
		{
			strictprimers->words[j] = strictprimers->words[i];
			strictprimers->strictcount[j] = strictprimers->strictcount[i];
			j++;
		}
	}
	
	if (j < strictprimers->size)
	{
		strictprimers->size = j;
		strictprimers->strictcount = ECOREALLOC(strictprimers->strictcount,
												sizeof(uint32_t)*strictprimers->size,
												"Cannot reallocate strict primer count table");
		strictprimers->words   = ECOREALLOC(strictprimers->words,
											sizeof(word_t)*strictprimers->size,
											"Cannot reallocate strict primer table");		
	}
	
	return strictprimers;
}

uint32_t filterMultiStrictPrimer(pwordcount_t strictprimers)
{
	uint32_t i;
	uint32_t w;

	for (i=0,w=0;i < strictprimers->size;i++)
	{
		if (w < i)
		{
			strictprimers->words[w]=strictprimers->words[i];
			strictprimers->strictcount[w]=strictprimers->strictcount[i];
		}
		if (! ISMULTIWORD(strictprimers->words[w]))
			w++;
	}

	strictprimers->size=w;
    strictprimers->strictcount = ECOREALLOC(strictprimers->strictcount,
    		                            sizeof(uint32_t)*strictprimers->size,
    		                            "Cannot reallocate strict primer count table");
    strictprimers->words   = ECOREALLOC(strictprimers->words,
										sizeof(word_t)*strictprimers->size,
										"Cannot reallocate strict primer table");

    return w;
}
