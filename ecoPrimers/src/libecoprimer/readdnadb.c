/*
 * readdnadb.c
 *
 *  Created on: 7 nov. 2008
 *      Author: coissac
 */

#include "ecoprimer.h"

pecodnadb_t readdnadb(const char *name, ecotaxonomy_t *taxonomy, uint32_t *size,poptions_t options)
{
	ecoseq_t      *seq;
	uint32_t       buffsize=100;
	pecodnadb_t db;

	db = ECOMALLOC(buffsize*sizeof(ecoseq_t*),"I cannot allocate db memory");


	for(seq=ecoseq_iterator(name), *size=0;
	    seq;
	    seq=ecoseq_iterator(NULL)
	   )
	{
	  if (isExampleTaxon(taxonomy,seq->taxid,options) ||
		  isCounterExampleTaxon(taxonomy,seq->taxid,options))
	  {
			if (*size==buffsize)
			{
				buffsize*=2;
				db = ECOREALLOC(db,buffsize*sizeof(ecoseq_t*),"I cannot allocate db memory");
			}
			db[*size]=seq;
			(*size)++;
	  }
	  else
	  {
		  delete_ecoseq(seq);
	  }
	};

	db = ECOREALLOC(db,(*size)*sizeof(ecoseq_t*),"I cannot allocate db memory");

	return db;
}


void printSeqTest(pecodnadb_t seqdb,uint32_t seqdbsize)
{
	uint32_t i;
	char ch[11];
	ch [10] = '\0';
	
	for (i=0; i < seqdbsize; i++)
	{
		strncpy (ch, seqdb[i]->SQ, 10);
		fprintf (stderr, "seq %d = %s\n", i, ch);
	}
	exit (0);
}
