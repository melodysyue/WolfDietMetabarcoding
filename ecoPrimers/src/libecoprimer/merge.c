/*
 * merge.c
 *
 *  Created on: 11 nov. 2008
 *      Author: coissac
 */

#include "ecoprimer.h"

static pmerge_t mergeInit(pmerge_t merge,pwordcount_t data,uint32_t s1,uint32_t s2);


static pmerge_t mergeInit(pmerge_t merge, pwordcount_t data, uint32_t s1, uint32_t s2)
{
	merge->words = data->words;
	merge->count = data->strictcount;
    merge->write = 0;
    merge->read1 = 0;
    merge->read2 = s1;
    merge->size  = s1+s2;
    return merge;
}


typedef enum {S1=1,S2=2,STACK=3} source_t;

void ecomerge(pwordcount_t data,uint32_t s1,uint32_t s2,uint32_t remainingSeq,uint32_t seqQuorum)
{

	/*
	 * data          / in out : the table that contains the two parts to be merged
	 * s1            / in     : end of the first part of the table
	 * s2            / in     : end of the second part of the table
	 * remainingSeq  / in     : the number of remaining seqs to be added to the table
	 * seqQuorum     / in     : the minimum number of sequences in which a pattern must appear
	 */

	merge_t   merged;
	source_t  source;
	word_t    currentword,tmpword;
	uint32_t  currentcount,tmpcount;
	int       same;
	queue_t   queue;
	int       nsame=0;
	uint32_t  maxcount=0;
	bool_t    writed=TRUE;

//	DEBUG_LOG("Coucou %p  s1= %d s2= %d",data,s1,s2)

	/*
	 * init the merged structure (used only for coding convenience, never returned, allocated on the C-stack)
	 * note that :
	 *     merged.words : hashcodes           (initialized to data->words)
	 *     merged.count : counts of each word (initialized to data->strictcount)
	 *     merged.read1 : index of the first word of the first subtable   (initialized to 0)
	 *     merged.read1 : index of the first word of the first subtable   (initialized to 0)
	 *     merged.read2 : index of the first word of the second subtable  (initialized to s1)
	 *     merged.size  : total size of the table (initialized to s1+s2)
	 *
	 * allocate a new stack of size min(s1, s2)
	 */

	(void)mergeInit(&merged,data,s1,s2);
	(void)newQueue(&queue,MINI(s1,s2));


	/* true until
	 * merged.read1 == s1 AND merged.read2 == merged.size, i.e. ALL words have been processed
	 */
	while (merged.read1 < s1 || merged.read2 < merged.size)
	{
		/*
		 * initialize current{word,count} from either STACK (if not empty) or first table (S1)
		 */
		if (! queue.empty)
		{
			currentword  = queue.words[queue.pop];
			currentcount = queue.count[queue.pop];
			source=STACK;
		}
		else
		{
			currentword  = merged.words[merged.read1];
			currentcount = merged.count[merged.read1];
			source=S1;
		}

		/*
		 * IF there are some words in the second subtable remaining to be processed AND
		 *    its first word is lower than current word
		 * THEN initialize current{word,count} from the second table (S2)
		 *
		 */
		if (merged.read2 < merged.size &&
				WORD(currentword) > WORD(merged.words[merged.read2]))
		{
			currentword  = merged.words[merged.read2];
			currentcount = merged.count[merged.read2];
			source  = S2;
		}

		/*
		 * record if the two words in the both subtable are the same
		 */
		same = (source != S2) && (WORD(currentword) == WORD(merged.words[merged.read2]));
		nsame+=same;

//		DEBUG_LOG("Merging : r1 = %d s1 = %d r2 = %d size = %d word = %s source=%u same=%u",merged.read1,s1,merged.read2-s1,merged.size,ecoUnhashWord(currentword,18),source,same)


		/*
		 * merge step (AND apply the quorum property)
		 * update merged.read1 AND/OR merged.read2
		 * record the word and its count in the table
		 */
		tmpword = merged.words[merged.write];
		tmpcount= merged.count[merged.write];

		merged.words[merged.write] = currentword;
		merged.count[merged.write] = currentcount;

		if (source != S2)
		{
			if (same)
			{
				merged.count[merged.write]+=merged.count[merged.read2];

				if (ISMULTIWORD(currentword) || ISMULTIWORD(merged.words[merged.read2]))
					merged.words[merged.write]=SETMULTIWORD(currentword);

				merged.read2++;
			}

			if (source==STACK)
				pop(&queue);
			merged.read1++;
		}
		else
			merged.read2++;

		if (writed && merged.read1 <= merged.write && merged.write < s1)
			push(&queue,tmpword,tmpcount);

		if (merged.count[merged.write] > maxcount)
			maxcount=merged.count[merged.write];

		writed = remainingSeq + merged.count[merged.write] >= seqQuorum;
        if (writed)
        	merged.write++;


//      else
//        	DEBUG_LOG("Remove word : %s count : %d remainingSeq : %d total : %d Quorum : %d",
//        			  ecoUnhashWord(currentword,18),merged.count[merged.write],remainingSeq,maxcount+remainingSeq,seqQuorum);

	}  /* while loop */

//	DEBUG_LOG("r1 : %d r2 : %d qsize : %d nsame : %d tot : %d write : %s count : %d source : %d size : %d pop : %d push : %d  empty : %d",merged.read1,merged.read2-s1,qsize,nsame,qsize+nsame,ecoUnhashWord(currentword,18),currentcount,source,queue.size,queue.pop,queue.push,queue.empty)


	/*
	 * finish the merging with words not processed (AND apply the quorum property)
	 * they are stored in the second subtable (IF)
	 *              OR in the queue           (ELSE)
	 */

	if (merged.read2 < merged.size)
		{
			//DEBUG_LOG("end1 %d %d/%d  %d/%d",merged.write,merged.read1,s1,merged.read2,merged.size);
		for (;merged.read2 < merged.size;merged.read2++)
		{
			merged.words[merged.write]=merged.words[merged.read2];
			merged.count[merged.write]=merged.count[merged.read2];
	        if (remainingSeq + merged.count[merged.write] >= seqQuorum)
	        	merged.write++;

		}
		}
	else {
		//DEBUG_LOG("end2 %d %d/%d  %d/%d",merged.write,merged.read1,s1,merged.read2,merged.size);
		while (! queue.empty)
		{
//			DEBUG_LOG("write : %s count : %d write : %d size : %d pop : %d push : %d  empty : %d",ecoUnhashWord(queue.words[queue.pop],18),queue.count[queue.pop],merged.write,queue.size,queue.pop,queue.push,queue.empty)
			merged.words[merged.write]=queue.words[queue.pop];
			merged.count[merged.write]=queue.count[queue.pop];
			pop(&queue);
	        if (remainingSeq + merged.count[merged.write] >= seqQuorum)
	        	merged.write++;
		}
		}

	data->size = merged.write;

	cleanQueue(&queue);

//	DEBUG_LOG("Max count : %d remainingSeq : %d total : %d Quorum : %d",maxcount,remainingSeq,maxcount+remainingSeq,seqQuorum)
//	DEBUG_LOG("Second word : %s",ecoUnhashWord(data->words[1],18))
//	DEBUG_LOG("Last  word : %s",ecoUnhashWord(data->words[data->size-1],18))


}
