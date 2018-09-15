/*
 * amplifiatree.c
 *
 *  Created on: 7 mars 2009
 *      Author: coissac
 */

#include  "ecoprimer.h"
#include <search.h>

static void cleanamplifia(pamplifia_t amplifia);
static void deleteamplifialist(pamplifialist_t list);
static int cmpamplifia(const void* p1,const void*p2);


static void cleanamplifiatlist(pamplifiacount_t list)
{
	if (list->amplifias)
		ECOFREE(list->amplifias,
				"Free amplifia list");
}

static void cleanamplifia(pamplifia_t amplifia)
{
   cleanamplifiatlist(&(amplifia->pcr));
}

static pamplifialist_t newamplifialist(pamplifialist_t parent, size_t size)
{
	pamplifialist_t tmp;

	tmp=ECOMALLOC(sizeof(amplifialist_t)+sizeof(amplifia_t)*(size-1),
			      "Cannot allocate new amplifia list");

	tmp->amplifiaslots=size;
	tmp->amplifiacount=0;
	tmp->next=NULL;

	if (parent)
		parent->next=(void*)tmp;

	return tmp;
}

static void deleteamplifialist(pamplifialist_t list)
{
	size_t i;

	if (list)
	{
		if (list->next)
		{
			deleteamplifialist(list->next);
			list->next=NULL;
		}
		for (i=0; i < list->amplifiacount; i++)
			cleanamplifia((list->amplifias)+i);

		ECOFREE(list,"Delete amplifia list");
	}

}

static int cmpamplifia(const void* p1,const void*p2)
{
	pamplifia_t pr1,pr2;

	pr1=(pamplifia_t)p1;
	pr2=(pamplifia_t)p2;

	if (pr1->p1 < pr2->p1) return -1;
	if (pr1->p1 > pr2->p1) return  1;

	if (pr1->asdirect1 < pr2->asdirect1) return -1;
	if (pr1->asdirect1 > pr2->asdirect1) return  1;

	if (pr1->p2 < pr2->p2) return -1;
	if (pr1->p2 > pr2->p2) return  1;

	if (pr1->asdirect2 < pr2->asdirect2) return -1;
	if (pr1->asdirect2 > pr2->asdirect2) return  1;

	return 0;
}

pamplifia_t amplifiaintree (amplifia_t key,
		            pamplifiatree_t amplifialist)
{
	if (!amplifialist->tree)
		return NULL;

	return *((pamplifia_t*)tsearch((const void *)(&key),
								&(amplifialist->tree),
						        cmpamplifia
							  ));
}

pamplifia_t insertamplifia(amplifia_t key,
		           pamplifiatree_t list)
{
	pamplifia_t current;
	pamplifia_t found;

	if (list->last->amplifiacount==list->last->amplifiaslots)
	{
		list->last->next=newamplifialist(list,100);
		list->last=list->last->next;
	}

	current = list->last->amplifias + list->last->amplifiacount;
	*current=key;

	found = *((pamplifia_t*)tsearch((const void *)current,
			                    &(list->tree),
			                    cmpamplifia));
	if (found==current)
		list->last->amplifiacount++;

	return found;
}

pamplifiatree_t initamplifiatree(pamplifiatree_t tree)
{
	if (!tree)
		tree = ECOMALLOC(sizeof(amplifiatree_t),"Cannot allocate amplifia tree");

	tree->first=newamplifialist(NULL,500);
	tree->last=tree->first;

	tree->tree=NULL;
}
