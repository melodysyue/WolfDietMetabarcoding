/*
 * taxstats.c
 *
 *  Created on: 12 mars 2009
 *      Author: coissac
 */

#include <search.h>
//void tdestroy (void *root, void (*free_node)(void *nodep));

#include  "ecoprimer.h"

static int cmptaxon(const void *t1, const void* t2);

void **tree_root = NULL;
int delete_passes = 0; 

void delete_twalkaction (const void *node, VISIT order, int level)
{
	switch (order)
	{
	case preorder:
		delete_passes++;
		break;
	case postorder:
		delete_passes++;
		break;
	case endorder:
		delete_passes++;
		break;
	case leaf:
		if (tree_root)
			tdelete (node, tree_root,cmptaxon);
		delete_passes++;
		break;
	}
}

void free_tree_nodes (void *tree)
{
	while (1)
	{
		delete_passes = 0;
		twalk (tree, delete_twalkaction);
		if (delete_passes <= 1) break;
	}
}

static int cmptaxon(const void *t1, const void* t2)
{
	const size_t taxid1=(size_t)t1;
	const size_t taxid2=(size_t)t2;

	// fprintf(stderr,"==> counted taxid1 : %d\n",taxid1);
	// fprintf(stderr,"==> counted taxid2 : %d\n",taxid2);

	if (taxid1 < taxid2)
		return -1;
	if (taxid1 > taxid2)
		return +1;
	return 0;
}

int32_t counttaxon(int32_t taxid)
{
	static void* taxontree=NULL;
	static int32_t taxoncount=0;

	// fprintf(stderr,"counted taxid : %d taxontree %p\n",taxid,taxontree);

	if (taxid==-1)
	{
		if (taxontree)
		{
			tree_root = (void **)&taxontree;
			//free_tree_nodes (taxontree);
			ECOFREE(taxontree,"Free taxon tree");
			tree_root = NULL;
		}
		taxontree=NULL;
		taxoncount=0;
		return 0;
	}


	if ((taxid > 0) && ((!taxontree) || (!tfind((void*)((size_t)taxid),&taxontree,cmptaxon))))
	{
		tsearch((void*)((size_t)taxid),&taxontree,cmptaxon);
		taxoncount++;
	}
	return taxoncount;
}

int32_t getrankdbstats(pecodnadb_t seqdb, uint32_t seqdbsize, ecotaxonomy_t *taxonomy,
					   poptions_t options)
{

	uint32_t i;
	ecotx_t  *taxon;
	ecotx_t  *tmptaxon;

	counttaxon(-1);
	options->intaxa = 0;

    for (i=0;i<seqdbsize;i++)
	{
    	taxon = &(taxonomy->taxons->taxon[seqdb[i]->taxid]);
	    seqdb[i]->isexample=isExampleTaxon(taxonomy,seqdb[i]->taxid,options);

		tmptaxon = eco_findtaxonatrank(taxon,
									   options->taxonrankidx);

	    // fprintf(stderr,"Taxid : %d %p\n",taxon->taxid,tmptaxon);

		if (tmptaxon)
		{
			// fprintf(stderr,"orig : %d trans : %d\n",taxon->taxid,
			// 		                                tmptaxon->taxid);

			seqdb[i]->ranktaxonid=tmptaxon->taxid;
			if (seqdb[i]->isexample)
				options->intaxa = counttaxon(tmptaxon->taxid);
		}
		else
			seqdb[i]->ranktaxonid=-1;
	}

	counttaxon(-1);
	options->outtaxa = 0;

    for (i=0;i<seqdbsize;i++)
		{
			if (seqdb[i]->ranktaxonid>=0 && !seqdb[i]->isexample)
				options->outtaxa = counttaxon(seqdb[i]->ranktaxonid);
		}

    return options->outtaxa + options->intaxa;
}


float taxonomycoverage(ppair_t pair, poptions_t options, pecodnadb_t seqdb,uint32_t seqdbsize)
{
	int32_t seqcount;
	int32_t i;
	int32_t incount=0;
	int32_t outcount=0;
	uint32_t j;


	memset (pair->coveredSeqs, 0, seqdbsize*sizeof (int));
	seqcount=pair->pcr.ampcount;

	counttaxon(-1);
	for (i=0; i < seqcount; i++)
		if (pair->pcr.amplifias[i].sequence->isexample
				&& pair->pcr.amplifias[i].sequence->ranktaxonid > 0 )
		{
			incount = counttaxon(pair->pcr.amplifias[i].sequence->ranktaxonid);

			for (j=0; j<seqdbsize; j++)
				if (pair->pcr.amplifias[i].sequence == seqdb[j])
					{pair->coveredSeqs[j] = 1; break;}
		}

	counttaxon(-1);
	for (i=0; i < seqcount; i++)
		if (!pair->pcr.amplifias[i].sequence->isexample
				&& pair->pcr.amplifias[i].sequence->ranktaxonid)
			outcount = counttaxon(pair->pcr.amplifias[i].sequence->ranktaxonid);


	pair->intaxa=incount;
	pair->outtaxa=outcount;
	pair->bc=(float)incount/options->intaxa;
	return pair->bc;
}

/*
static int cmpamp(const void *ampf1, const void* ampf2)
{
	int i;
	int j = 0;
	int incr = 1;
	char cd1;
	char cd2;
	int chd = 0;
	int len = 0;

	pamptotaxon_t pampf1 = (pamptotaxon_t) ampf1;
	pamptotaxon_t pampf2 = (pamptotaxon_t) ampf2;


	if (pampf1->strand != pampf2->strand)
	{
		incr = -1;
		j = pampf1->length - 1;
		
		if (pampf2->strand)
		{
			pampf1 = (pamptotaxon_t) ampf2;
			pampf2 = (pamptotaxon_t) ampf1;
			chd = 1;
		}
		//j = pampf2->length - 1; should have been here and pampf2 instead of pampf1?
	}

	len = (pampf1->length <= pampf2->length)? pampf1->length: pampf2->length;

	for (i = 0; i < len; i++, j += incr)
	{
		cd1 = pampf1->amplifia[i];
		if (incr == -1)
			cd2 = ecoComplementChar(pampf2->amplifia[j]);
		else
			cd2 = pampf2->amplifia[j];

		if (cd1 < cd2) return chd ? 1: -1;
		if (cd2 < cd1) return chd ? -1: 1;
	}

	if (pampf1->length > pampf2->length) return chd ? -1: 1;
	if (pampf2->length > pampf1->length) return chd ? 1: -1;

	return 0;
}*/


static int cmpamp(const void *ampf1, const void* ampf2)
{
	int i;
	char cd1;
	char cd2;
	int len = 0;
	char *ch1;
	char *ch2;
	int incr1;
	int incr2;

	pamptotaxon_t pampf1 = (pamptotaxon_t) ampf1;
	pamptotaxon_t pampf2 = (pamptotaxon_t) ampf2;

	ch1 = pampf1->amplifia;
	ch2 = pampf2->amplifia;
	
	incr1 = 1;
	incr2 = 1;
	
	if (!pampf1->strand)
		incr1 = -1;
	if (!pampf2->strand)
		incr2 = -1;
	
	len = (pampf1->length <= pampf2->length)? pampf1->length: pampf2->length;
	for (i = 0; i < len; i++)
	{
		cd1 = *ch1;
		if (incr1 == -1)
			cd1 = ecoComplementChar(*ch1);
		
		cd2 = *ch2;
		if (incr2 == -1)
			cd2 = ecoComplementChar(*ch2);
		
		if (cd1 < cd2) return -1;
		if (cd2 < cd1) return 1;
		
		ch1 += incr1;
		ch2 += incr2;
	}
	
	if (pampf1->length > pampf2->length) return 1;
	if (pampf2->length > pampf1->length) return -1;

	return 0;
}

void twalkaction (const void *node, VISIT order, int level)
{
	int32_t *taxid = (int32_t*)node;
	//const size_t taxid=(size_t)node;
	//printf ("\t%d:%p, ", *taxid, node);
	counttaxon(*taxid);
}

int32_t gtxid;
void twalkaction2 (const void *node, VISIT order, int level)
{
	int32_t *pt = (int32_t *) node;
	gtxid = *pt;
}

void taxonomyspecificity (ppair_t pair, pecodnadb_t seqdb,uint32_t seqdbsize)
{
	uint32_t i, j;
	uint32_t ampfindex = 0;
	int32_t taxid;
	uint32_t wellidentifiedcount;
	
	void *ampftree = NULL;
	pamptotaxon_t pcurrentampf;
	pamptotaxon_t *ptmp;

	pamptotaxon_t ampfwithtaxtree = ECOMALLOC(sizeof(amptotaxon_t) * pair->pcr.ampcount,"Cannot allocate amplifia tree");

	for (i = 0; i < pair->pcr.ampcount; i++)
	{
		/*populate taxon ids tree against each unique amplifia
		i.e set of taxon ids for each amplifia*/
		if (pair->pcr.amplifias[i].sequence->isexample)
		{
			ampfwithtaxtree[ampfindex].amplifia = pair->pcr.amplifias[i].amplifia;
			ampfwithtaxtree[ampfindex].strand = pair->pcr.amplifias[i].strand;
			ampfwithtaxtree[ampfindex].length = pair->pcr.amplifias[i].length;
			pcurrentampf = &ampfwithtaxtree[ampfindex];
			taxid = pair->pcr.amplifias[i].sequence->ranktaxonid;
			ptmp = tfind((const void*)pcurrentampf, &ampftree, cmpamp);
			if (ptmp == NULL)
			{
				pcurrentampf = &ampfwithtaxtree[ampfindex];
				tsearch((void*)pcurrentampf,&ampftree,cmpamp);
				ampfindex++;
			}
			else
				pcurrentampf = *ptmp;

			if (tfind((void*)((size_t)taxid), &(pcurrentampf->taxontree), cmptaxon) == NULL)
			{
				pcurrentampf->taxoncount++;
				tsearch((void*)((size_t)taxid),&(pcurrentampf->taxontree),cmptaxon);
			}
		}
	}

	memset (pair->wellIdentifiedSeqs, 0, seqdbsize*sizeof (int));
	//counttaxon(-1);
	for (i = 0; i < ampfindex; i++)
	{
		if (ampfwithtaxtree[i].taxoncount > 1)
		{
			//printf ("\nampfwithtaxtree[i].taxoncount: %d\n", ampfwithtaxtree[i].taxoncount);
			//twalk(ampfwithtaxtree[i].taxontree, twalkaction);
		}
		//TR 5/9/10 - added code for well identified seqs
		else if(ampfwithtaxtree[i].taxoncount == 1) /*well identified*/
		{
			gtxid = -1;
			twalk(ampfwithtaxtree[i].taxontree, twalkaction2);
			
			if (gtxid != -1)
			{
				for (j = 0; j < seqdbsize; j++)
					if (seqdb[j]->ranktaxonid == gtxid 
							&& seqdb[j]->isexample
							&&(pair->p1->directCount[j] > 0 
							|| pair->p1->reverseCount[j] > 0)
							&& (pair->p2->directCount[j] > 0 
							|| pair->p2->reverseCount[j] > 0))
					{
						pair->wellIdentifiedSeqs[j] = 1;
					}
			}
		}
	}
	//printf ("\n");
	counttaxon(-1);
	wellidentifiedcount = 0;
	for (j = 0; j < seqdbsize; j++)
		if (pair->wellIdentifiedSeqs[j] == 1)
			counttaxon(seqdb[j]->ranktaxonid);
	wellidentifiedcount = counttaxon(-2);
	//pair->notwellidentifiedtaxa = counttaxon(-2);
	pair->notwellidentifiedtaxa = (pair->intaxa-wellidentifiedcount); //counttaxon(-2);
	//pair->bs = ((float)pair->intaxa - (float)pair->notwellidentifiedtaxa) / pair->intaxa;
	pair->bs = ((float)wellidentifiedcount) / (float)pair->intaxa;
	
	ECOFREE (ampfwithtaxtree, "Free amplifia table");

}
