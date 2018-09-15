/*
 * filtering.c
 *
 *  Created on: 12 mai 2009
 *      Author: coissac
 */

#include "ecoprimer.h"
#include <string.h>
#include <math.h>

#include "hashencoder.h"

static int32_t *ecoFilteringHashSequence(int32_t *dest,
										 uint32_t circular,
										 uint32_t doublestrand,
										 ecoseq_t *seq,
										 uint32_t *size);





static	int32_t *ecoFilteringHashSequence(int32_t *dest,
										  uint32_t circular,
										  uint32_t doublestrand,
										  ecoseq_t *seq,
										  uint32_t *size)
{

	/*
	 * This function aims at building a vector of count for every possible hash codes
	 *
	 * The function must be applied once on each sequence
	 *
	 * The function allocates memory on the first call for the dest table
	 * The function also allocates memory for the static temporary table in_last_seq and
	 * the function must be called with *dest == -1 in order to free this temporary table
	 *
	 */
	static char    *in_last_seq=NULL;
	uint32_t i=0;
	uint32_t j;
	char *base;
	int8_t code;
	int32_t error=0;
	word_t word=0;
	word_t antiword=0;
	uint32_t goodword;
	uint32_t lmax=0;

	// run on the first call;

	if (dest==(void*)-1)
	{
		if (in_last_seq) ECOFREE(in_last_seq,"Free in last seq table");
		return NULL;
	}


	/* FWORDSIZE = 13 => *size =    67 108 864
	 * FWORDSIZE = 14 => *size =   268 435 456
	 * FWORDSIZE = 15 => *size = 1 073 741 824
	 */

	*size = pow(4,FWORDSIZE);

	/*
	 * in_last_seq is a vector of char as it is just to avoid counting twice (or more) a hash code (a DNA word)
	 * it is set to zero (with memset) and then filled with ones for each word belonging to the sequence
	 */

	if (!in_last_seq)
		in_last_seq = ECOMALLOC(*size*sizeof(char),
								"Cannot allocate filtering hash table");

	memset(in_last_seq,0,*size*sizeof(char));

	/*
	 * Allocate (on first call) the memory for the table of counts
	 */

	if (!dest)
	{
		dest = ECOMALLOC(*size*sizeof(int32_t),
				               "Cannot allocate filtering hash table");
		memset(dest,0,*size*sizeof(int32_t));
	}

	lmax = seq->SQ_length;
	if (!circular)
       lmax-= FWORDSIZE-1;



//	DEBUG_LOG("Sequence %s @ %d : %18.18s",seq->AC,i,(seq->SQ+i));

	/*
	 * Compute first word of seq
	 */

	for (i=0, base = seq->SQ; i < FWORDSIZE && i < lmax; i++,base++)
	{
		error<<= 1;
		error&=ERRORMASK(FWORDSIZE);

		code = encoder[(*base) - 'A'];
		if (code <0)
		{
			code = 0;
			error|= 1;
		}


		word=RAPPENDBASE(word,FWORDSIZE,code);
		if (doublestrand)
			antiword=LAPPENDBASE(antiword,FWORDSIZE,code);
	}

	if (!error && i==FWORDSIZE)
	{

		goodword=(uint32_t)((doublestrand) ? MINI(word,antiword):word);

		/*
		 * FB: I don't think the test is necessary as the table has just been initialized
		 */
		if (!in_last_seq[goodword])
		{
			in_last_seq[goodword]=1;
			dest[goodword]++;
		}
	}

	/*
	 * compute and store counts (avoid counting twice a word) for the other words of the seq
	 */

	for (j=1; j < lmax; j++,i++,base++)
	{

//		DEBUG_LOG("Sequence %s @ %d : %18.18s",seq->AC,j,(seq->SQ+j));

							/* roll over the sequence for circular ones */
		if (i==(uint32_t)seq->SQ_length) base=seq->SQ;

		error<<= 1;
		error&=ERRORMASK(FWORDSIZE);

		//code = -1;
		//if((*base) >= 'A' && (*base) <= 'Z')
		code = encoder[(*base) - 'A'];
		if (code <0)
		{
			code = 0;
			error|= 1;
		}

		word=RAPPENDBASE(word,FWORDSIZE,code);
		if (doublestrand)
			antiword=LAPPENDBASE(antiword,FWORDSIZE,code);

		if (!error)
		{
			if (doublestrand)
				goodword=(uint32_t)MINI(word,antiword);
			else
				goodword=word;
			if (!in_last_seq[goodword])
			{
				in_last_seq[goodword]=1;
				dest[goodword]++;
			}
		}

	}


	return dest;

}


int32_t *filteringSeq(pecodnadb_t database, uint32_t seqdbsize,
					 uint32_t exampleCount,poptions_t options,uint32_t *size,int32_t  sequenceQuorum)
{
	int32_t *wordscount=NULL;
	int32_t keep=0;
	uint32_t i,j=0;

	for (i=0;i<seqdbsize;i++)
    {
    	if (database[i]->isexample && database[i]->SQ_length > options->primer_length)
    	{
    		j++;
    		wordscount=ecoFilteringHashSequence(wordscount,
									 options->circular,
                                     options->doublestrand,
                                     database[i],
                                     size);
    	}
    	fprintf(stderr,"  Filtered sequences %5u/%5u          \r",j,exampleCount);

    }

	fprintf(stderr,"\n");

	for (i=0;i<*size;i++)
		if (wordscount[i] >= sequenceQuorum)
			keep++;


	(void)ecoFilteringHashSequence((int32_t*)-1,
									options->circular,
	                                options->doublestrand,
	                                NULL,
	                                NULL);

	fprintf(stderr,"ok\n Considered word of size %d for filtering : %d\n",FWORDSIZE,keep);
	return wordscount;

}
