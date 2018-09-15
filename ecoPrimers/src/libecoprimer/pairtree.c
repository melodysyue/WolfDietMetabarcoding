/*
 * pairtree.c
 *
 *  Created on: 7 mars 2009
 *      Author: coissac
 */

#include  "ecoprimer.h"
#include <search.h>

static void cleanpair(ppair_t pair);
static void deletepairlist(ppairlist_t list);
static int cmppair(const void* p1,const void*p2);


static void cleanamplifiatlist(pamplifiacount_t list)
{
	if (list->amplifias)
		ECOFREE(list->amplifias,
				"Free amplifia list");
}

static void cleanpair(ppair_t pair)
{
   cleanamplifiatlist(&(pair->pcr));
}

static ppairlist_t newpairlist(ppairlist_t parent, size_t size)
{
	ppairlist_t tmp;

	tmp=ECOMALLOC(sizeof(pairlist_t)+sizeof(pair_t)*(size-1),
			      "Cannot allocate new pair list");

	tmp->pairslots=size;
	tmp->paircount=0;
	tmp->next=NULL;

	if (parent)
		parent->next=(void*)tmp;


	return tmp;
}

static void deletepairlist(ppairlist_t list)
{
	size_t i;

	if (list)
	{
		if (list->next)
		{
			deletepairlist(list->next);
			list->next=NULL;
		}
		for (i=0; i < list->paircount; i++)
			cleanpair((list->pairs)+i);

		ECOFREE(list,"Delete pair list");
	}

}

static int cmppair(const void* p1,const void*p2)
{
	ppair_t pr1,pr2;

	pr1=(ppair_t)p1;
	pr2=(ppair_t)p2;

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

ppair_t pairintree (pair_t key,
		            ppairtree_t pairlist)
{
	if (!pairlist->tree)
		return NULL;

	return *((ppair_t*)tsearch((const void *)(&key),
								&(pairlist->tree),
						        cmppair
							  ));
}

ppair_t insertpair(pair_t key,
		           ppairtree_t list)
{
	ppair_t current;
	ppair_t found;

	if (list->last->paircount==list->last->pairslots)
	{
		list->last->next=newpairlist(list->last,100);
		list->last=list->last->next;
	}

	current = list->last->pairs + list->last->paircount;
	*current=key;

	found = *((ppair_t*)tsearch((const void *)current,
			                    &(list->tree),
			                    cmppair));
	if (found==current)
		list->last->paircount++;

	return found;
}

ppairtree_t initpairtree(ppairtree_t tree)
{

	if (!tree)
		tree = ECOMALLOC(sizeof(pairtree_t),"Cannot allocate pair tree");

	tree->first=newpairlist(NULL,300);
	tree->last=tree->first;

	tree->tree=NULL;
	tree->count=0;

	return tree;
}
