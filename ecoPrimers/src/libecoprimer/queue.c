/*
 * queue.c
 *
 *  Created on: 14 nov. 2008
 *      Author: coissac
 */

#include "ecoprimer.h"



pqueue_t newQueue(pqueue_t queue, uint32_t size)
{
	if (!queue)
		queue = ECOMALLOC(sizeof(queue_t),"Cannot allocate queue structure");

	queue->size=0;

	resizeQueue(queue,size);

	return queue;

}

pqueue_t resizeQueue(pqueue_t queue, uint32_t size)
{
	queue->pop=0;
	queue->push=0;
	queue->empty=TRUE;
	queue->full=FALSE;

	if (!queue->size)
	{
		queue->count=ECOMALLOC(size * sizeof(uint32_t),
				               "Cannot allocate count queue array"
				              );
		queue->words=ECOMALLOC(size * sizeof(word_t),
				               "Cannot allocate word queue array"
				              );
		queue->size=size;
	}
	else if (size > queue->size)
	{
		queue->count=ECOREALLOC(queue->count,
				                size * sizeof(uint32_t),
				               "Cannot allocate count queue array"
				              );
		queue->words=ECOREALLOC(queue->words,
				                size * sizeof(word_t),
				               "Cannot allocate word queue array"
				              );

		queue->size=size;
	}

	return queue;
}

pqueue_t cleanQueue(pqueue_t queue)
{
	if (queue->size)
	{
		if (queue->count)
			ECOFREE(queue->count,"Free count queue");
		if (queue->words)
			ECOFREE(queue->words,"Free words queue");
	}

	queue->size=0;

	return queue;
}

void push(pqueue_t queue, word_t word, uint32_t count)
{
	ECO_ASSERT(!queue->full,"Queue is full");

	queue->count[queue->push]=count;
	queue->words[queue->push]=word;

	queue->push++;

	if (queue->push==queue->size)
		queue->push=0;

	queue->full=queue->push==queue->pop;
	queue->empty=FALSE;
}

void pop(pqueue_t queue)
{
	ECO_ASSERT(!queue->empty,"Queue is empty");
	queue->pop++;

	if (queue->pop==queue->size)
		queue->pop=0;

	queue->empty=queue->push==queue->pop;
	queue->full=FALSE;
}
