/*
 * hashsequence.c
 *
 *  Created on: 7 nov. 2008
 *      Author: coissac
 */


#include "ecoprimer.h"

static int cmpword(const void *x,const void *y);

#include "hashencoder.h"

uint32_t ecoWordCount(uint32_t wordsize, uint32_t circular, ecoseq_t *seq)
{
	uint32_t wordcount;

	wordcount = seq->SQ_length;

	if (!circular) wordcount-=wordsize-1;

	return wordcount;
}

pword_t ecoHashSequence(pword_t dest,
					    uint32_t wordsize,
					    uint32_t circular,
					    uint32_t doublestrand,
					    ecoseq_t *seq,
					    uint32_t *size,
					    int32_t  *neededWords,
					    uint32_t neededWordCount,
					    int32_t quorum)
{

	/*
	 * dest            / out    : words of the hashed sequence position per position
	 * wordsize        / in     : size of the word to be hashed (record error for that size) BUT not equal to FWORDSIZE ...
	 *                            ... the size of the word REALLY returned as a result
	 * circular        / in     : is the sequence circular
	 * doublestrand    / in     : if we have to hash on both strands of the sequence
	 * seq             / in     : the sequence in ecoseq format
	 * size            / out    : number of hashed words (size of the dest vector)
	 * neededWordCount / in     : table hash codes of word counts in the full DB (used to filter words)
	 * quorum          / in     : minimum quorum used to filter words based on the neededWordCount table
	 */

	uint32_t i=0;
	uint32_t j;
	char *base;
	int8_t code;
	int32_t error=0;
	word_t word=0;
	word_t antiword=0;
	word_t goodword;
	uint32_t lmax=0;

	(*size)=0;

	lmax = seq->SQ_length;
	if (!circular)
       lmax-= wordsize-1;

	if (!dest)
		dest = ECOMALLOC(lmax*sizeof(word_t),
				         "I cannot allocate memory for sequence hashing"
				        );

	//DEBUG_LOG("Sequence %s @ %d : %18.18s",seq->AC,i,(seq->SQ+i));

	for (i=0, base = seq->SQ; i < wordsize && i < lmax; i++,base++)
	{

		error<<= 1;
		error&=ERRORMASK(wordsize);

		code = encoder[(*base) - 'A'];
		if (code <0)
		{
			code = 0;
			error|= 1;
		}


		word=RAPPENDBASE(word,wordsize,code);

		if (doublestrand)
			antiword=LAPPENDBASE(antiword,wordsize,code);

		if (neededWordCount && i>=(FWORDSIZE-1))
		{

			goodword = (doublestrand) ? MINI(FILTERWORD(word),CFILTERWORD(antiword,wordsize)):FILTERWORD(word);
			if (neededWords[(uint32_t)goodword]<quorum)
				error|= (1 << (FWORDSIZE-1));

		}

	}


	if (!error && i==wordsize)
	{
		dest[*size]=(doublestrand) ? MINI(word,antiword):word;
		(*size)++;
	}


	for (j=1; j < lmax; j++,i++,base++)
	{

		//DEBUG_LOG("Sequence %s @ %d : %18.18s",seq->AC,j,(seq->SQ+j));

							/* roll over the sequence for circular ones */

		if (i==(uint32_t)seq->SQ_length) base=seq->SQ;

		error<<= 1;
		error&=ERRORMASK(wordsize);

		code = encoder[(*base) - 'A'];
		if (code <0)
		{
			code = 0;
			error|= 1;
		}

		word=RAPPENDBASE(word,wordsize,code);
		if (doublestrand)
			antiword=LAPPENDBASE(antiword,wordsize,code);

		if (neededWordCount)
		{
			goodword = (doublestrand) ? MINI(FILTERWORD(word),CFILTERWORD(antiword,wordsize)):FILTERWORD(word);
			if (neededWords[(uint32_t)goodword]<quorum)
				error|= (1 << (FWORDSIZE-1));
//			else
//				DEBUG_LOG("%s goodword = %p %d/%d (pos:%d error:%d)",seq->AC,goodword,neededWords[(uint32_t)goodword],quorum,i,error);

		}


		if (!error)
		{
			dest[*size]=(doublestrand) ? MINI(word,antiword):word;
			(*size)++;
		}

	}
	//DEBUG_LOG("%s goodword = %d",seq->AC,*size);
	return dest;

}

uint32_t ecoCompactHashSequence(pword_t table,uint32_t size)
{
	/*
	 *
	 * MULTIWORD is a word occurring more than once in a sequence
	 *
	 */

	uint32_t i,j;
	word_t  current;
//	bool_t here=FALSE;

	sortword(table,size);

	current = 0;
	current=SETMULTIWORD(current);   /* build impossible word for the first loop cycle */

//	if (strcmp(ecoUnhashWord(table[size-1],18),"GTTTGTTCAACGATTAAA")==0)
//		here=TRUE;

	for (i=0,j=0; j < size;j++)
	{
		if (WORD(table[j])!=current)
		{
			current =table[j];
			table[i]=current;
			i++;
		}
		else
			table[i]=SETMULTIWORD(table[i]);
	}

//		if (strcmp(ecoUnhashWord(WORD(table[i-1]),18),"TACGACCTCGATGTTGGA")==0)
//			DEBUG_LOG("winner %d",i)

	return i;
}

const char* ecoUnhashWord(word_t word,uint32_t size)
{
	static char buffer[32];
	static char decode[]="ACGT";

	uint32_t i;

	for (i=0; i < size; i++)
	{
		buffer[i]=decode[(word >> (2 * (size - 1 -i))) & 3];
	}

	buffer[size]=0;

	return buffer;
}

word_t ecoComplementWord(word_t word,uint32_t size)
{
	word_t rep=0;
	uint32_t i;

//	DEBUG_LOG("%llx  %llx",word,~word);
	word=(~word) & WORDMASK(size);
	for (i=0;i < size; i++)
	{

		rep = RAPPENDBASE(rep,size,word & 3LLU);
//		DEBUG_LOG("%016llx  %016llx  %016llx",word,word & 3LLU,rep);
		word>>=2;
	}
//	DEBUG_LOG("Complemented = %s",ecoUnhashWord(rep,18));
	return rep;

}

static int cmpword(const void *x,const void *y)
{
	word_t w1 = *(pword_t)x;
	word_t w2 = *(pword_t)y;

	w1 = WORD(w1);
	w2 = WORD(w2);

	if (w1 < w2)
		return -1;
	if (w1 > w2)
		return +1;

	return 0;
}

uint32_t ecoFindWord(pwordcount_t table,word_t word)
{
	pword_t dest;

	dest = (pword_t)bsearch((const void*)&word,(const void*)table->words,table->size,sizeof(word_t),cmpword);

	if (dest)
		return dest - table->words;
	else
		return ~0;
}

char ecoComplementChar(char base)
{
	return (base < 4)? !base & 3: 4;
}

