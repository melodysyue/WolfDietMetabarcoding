/*
 * pairs.c
 *
 *  Created on: 15 déc. 2008
 *      Author: coissac
 */

#include  "ecoprimer.h"
#include  <string.h>
#include  <stdlib.h>
#include "../libthermo/thermostats.h"

static void buildPrimerPairsForOneSeq(uint32_t seqid,
									  pecodnadb_t seqdb,
		 							  pprimercount_t primers,
								 	  ppairtree_t pairs,
									  poptions_t options);






/*************************************
 *
 *       pair collection management
 *
 *************************************/

#ifdef MASKEDCODE

char *addamplifiasetelem (ppair_t pair, char* amplifia, int32_t taxid)
{
	uint32_t i;
	uint32_t j;
	char *ampused = NULL;

	if(pair->ampsetcount == 0)
	{
		pair->ampsetcount = 500;
		pair->ampsetindex = 0;
		pair->ampset = ECOMALLOC(pair->ampsetcount * sizeof(ampseqset_t),"Cannot allocate amplifia set");
	}

	for (i = 0; i < pair->ampsetindex; i++)
	{
		if (strcmp (pair->ampset[i].amplifia, amplifia) == 0)
		{
			ampused = pair->ampset[i].amplifia;
			break;
		}
	}

	if (i == 0)
	{
		pair->ampset[i].seqidcount = 100;
		pair->ampset[i].seqidindex = 0;
		pair->ampset[i].taxonids = ECOMALLOC(pair->ampset[i].seqidcount * sizeof(uint32_t),"Cannot allocate amplifia sequence table");
	}

	if (pair->ampsetindex == pair->ampsetcount)
	{
		pair->ampsetcount += 500;
		pair->ampset = ECOREALLOC(pair->ampset, pair->ampsetcount * sizeof(ampseqset_t), "Cannot allocate amplifia set");
	}

	if (pair->ampset[i].seqidindex == pair->ampset[i].seqidcount)
	{
		pair->ampset[i].seqidcount += 100;
		pair->ampset[i].taxonids = ECOREALLOC(pair->ampset[i].taxonids, pair->ampset[i].seqidcount * sizeof(int32_t), "Cannot allocate amplifia sequence table");
	}

	if (pair->ampset[i].amplifia == NULL)
	{
		pair->ampset[i].amplifia = amplifia;
		pair->ampsetindex++;
	}

	for (j = 0; j < pair->ampset[i].seqidindex; j++)
	{
		if (pair->ampset[i].taxonids[j] == taxid) break;
	}

	if (j == pair->ampset[i].seqidindex)
		pair->ampset[i].taxonids[pair->ampset[i].seqidindex++] = taxid;
	return ampused;
}

void addtaxampsetelem (ppair_t pair, int32_t taxid, char *amplifia)
{
	uint32_t i;
	uint32_t j;

	if(pair->taxsetcount == 0)
	{
		pair->taxsetcount = 500;
		pair->taxsetindex = 0;
		pair->taxset = ECOMALLOC(pair->taxsetcount * sizeof(taxampset_t),"Cannot allocate taxon set");
	}

	for (i = 0; i < pair->taxsetindex; i++)
	{
		if (pair->taxset[i].taxonid == taxid) break;
	}

	if (i == 0)
	{
		pair->taxset[i].amplifiacount = 100;
		pair->taxset[i].amplifiaindex = 0;
		pair->taxset[i].amplifia = ECOMALLOC(pair->taxset[i].amplifiacount * sizeof(char *),"Cannot allocate amplifia table");
	}

	if (pair->taxsetindex == pair->taxsetcount)
	{
		pair->taxsetcount += 500;
		pair->taxset = ECOREALLOC(pair->taxset, pair->taxsetcount * sizeof(taxampset_t), "Cannot allocate taxon set");
	}

	if (pair->taxset[i].amplifiaindex == pair->taxset[i].amplifiacount)
	{
		pair->taxset[i].amplifiacount += 100;
		pair->taxset[i].amplifia = ECOREALLOC(pair->taxset[i].amplifia, pair->taxset[i].amplifiacount * sizeof(char *), "Cannot allocate amplifia table");
	}

	if (pair->taxset[i].taxonid == 0)
	{
		pair->taxset[i].taxonid = taxid;
		pair->taxsetindex++;
	}

	for (j = 0; j < pair->taxset[i].amplifiaindex; j++)
	{
		if (strcmp(pair->taxset[i].amplifia[j], amplifia) == 0) break;
	}

	if (j == pair->taxset[i].amplifiaindex)
	{
		pair->taxset[i].amplifia[j] = amplifia;
		pair->taxset[i].amplifiaindex++;
	}
}

char *getamplifia (pecoseq_t seq, uint32_t start, uint32_t len)
{
	fprintf(stderr,"start : %d length : %d\n",start,len);
	char *amplifia = ECOMALLOC((len + 1) * sizeof(char),"Cannot allocate amplifia");
	char *seqc = &seq->SQ[start];

	strncpy(amplifia, seqc, len);
	return amplifia;
}

#endif

/*TR: Added*/
ppairtree_t buildPrimerPairs(pecodnadb_t seqdb,uint32_t seqdbsize,pprimercount_t primers,poptions_t options)
{
	uint32_t i;
	ppairtree_t primerpairs;
	
	primerpairs = initpairtree(NULL);

	for (i=0; i < seqdbsize; i++)
	{
		buildPrimerPairsForOneSeq(i, seqdb, primers, primerpairs, options);
	}
	return primerpairs;
}

#define DMAX (2000000000)

static void buildPrimerPairsForOneSeq(uint32_t seqid,
									  pecodnadb_t seqdb,
		 							  pprimercount_t primers,
								 	  ppairtree_t pairs,
									  poptions_t options)
{
	static uint32_t    paircount=0;
	uint32_t           i,j,k;
	uint32_t           matchcount=0;
	pprimermatch_t     matches = NULL;
	//primermatchcount_t seqmatchcount;
	ppair_t            pcurrent;
	pair_t			   current;
	pprimer_t		   wswp;
	bool_t			   bswp;
	size_t			   distance;
	bool_t             strand;
	//char 			   prmr[50];
	//float			   mtemp;
	word_t 			   w1, w1a, omask = (0x1L << (options->strict_three_prime*2)) -1;
	word_t 			   w2, w2a;//, wtmp;
	uint32_t   		   bp1,bp2;
		
	//prmr[options->primer_length] = '\0';
	
	for (i=0;i < primers->size; i++)
	{
		matchcount+=primers->primers[i].directCount[seqid];
		matchcount+=primers->primers[i].reverseCount[seqid];
	}

	if (matchcount <= 0)
		return;

	matches = ECOMALLOC(matchcount * sizeof(primermatch_t),"Cannot allocate primers match table");

	for (i=0,j=0;i < primers->size; i++)
	{
		if (primers->primers[i].directCount[seqid])
		{
			if (primers->primers[i].directCount[seqid]==1)
			{
			matches[j].primer = primers->primers+i;
			matches[j].strand=TRUE;
				matches[j].position=primers->primers[i].directPos[seqid].value;
				j++;
			}
			else for (k=0; k < primers->primers[i].directCount[seqid]; k++,j++)
			{
				matches[j].primer = primers->primers+i;
				matches[j].strand=TRUE;
				matches[j].position=primers->primers[i].directPos[seqid].pointer[k];
			}
		}

		if (primers->primers[i].reverseCount[seqid])
		{
			if (primers->primers[i].reverseCount[seqid]==1)
			{
			matches[j].primer = primers->primers+i;
			matches[j].strand=FALSE;
				matches[j].position=primers->primers[i].reversePos[seqid].value;
				j++;
			}
			else for (k=0; k < primers->primers[i].reverseCount[seqid]; k++,j++)
			{
				matches[j].primer = primers->primers+i;
				matches[j].strand=FALSE;
				matches[j].position=primers->primers[i].reversePos[seqid].pointer[k];
			}
		}
	}

	if (matchcount>1)
	{
//		fprintf(stderr,"\n====================================\n");

		sortmatch(matches,matchcount); // sort in ascending order by position

		for (i=0; i < matchcount;i++)
		{
			// For all primers matching the sequence

			/*for(j=i+1;
			       (j<matchcount)
			    && ((distance=matches[j].position - matches[i].position - options->primer_length) < options->lmax);
			   j++
			   )//*/
			for (j=i+1; j<matchcount; j++)
			{
				if (matches[j].position - matches[i].position <= options->primer_length) continue;
				distance = matches[j].position - matches[i].position - options->primer_length;
				if (distance >= options->lmax) break;
				

			   // For all not too far primers

		       if ( (matches[i].primer->good || matches[j].primer->good)
		    		&& (distance > options->lmin)
		    		)
		       {
		    	   // If possible primer pair
		    	   current.p1 = matches[i].primer;
		    	   current.asdirect1=matches[i].strand;
		    	   current.p2 = matches[j].primer;
		    	   current.asdirect2= !matches[j].strand;
		    	   current.maxd=DMAX;
		    	   current.mind=DMAX;
		    	   current.sumd=0;
		    	   current.amplifiacount=0;
	    		   current.inexample=0;
	    		   current.outexample=0;
				   current.curseqid = 0;
				   current.refsequence=-1;
				   //current.p1temp = 100;
				   //current.p1mintemp = 100;
				   //current.p2temp = 100;
				   //current.p2mintemp = 100;

		    	   // Standardize the pair
	    		   strand = current.p2->word > current.p1->word;
		    	   if (!strand)
		    	   {
		    		   wswp = current.p1;
		    		   current.p1=current.p2;
		    		   current.p2=wswp;

		    		   bswp = current.asdirect1;
		    		   current.asdirect1=current.asdirect2;
		    		   current.asdirect2=bswp;
		    	   }

		    	   
		    	   //Code to make sure that if -3 option is given then
		    	   //3' end must match upto given number of base pairs
		    	   if (options->strict_three_prime > 0)
		    	   {
		    		   w1 = current.p1->word;
		    		   w2 = current.p2->word;
		    		   if (!current.asdirect1) //make sure that word is from 5' to 3'
		    			   w1=ecoComplementWord(w1,options->primer_length);

		    		   if (!current.asdirect2) //make sure that word is from 5' to 3'
		    			   w2=ecoComplementWord(w2,options->primer_length);
		    		   //now both w1 and w2 are from 5' to 3' end
		    		   bp1 = matches[i].position;
		    		   bp2 = matches[j].position;
		    		   if (!strand)
		    		   {
			    		   bp1 = matches[j].position;
			    		   bp2 = matches[i].position;
		    		   }
		    		   //get word of first approximate repeat
		    		   w1a = extractSite(seqdb[seqid]->SQ,bp1,options->primer_length,strand);
		    		   //get word of second approximate repeat
		    		   w2a = extractSite(seqdb[seqid]->SQ,bp2,options->primer_length,!strand);
		    		   		    		  		    		   
		    		   w1 = w1 & omask;   //keep only strict_three_prime bases on the right (3') end
		    		   w2 = w2 & omask;   //keep only strict_three_prime bases on the right (3') end
		    		   w1a = w1a & omask; //keep only strict_three_prime bases on the right (3') end
		    		   w2a = w2a & omask; //keep only strict_three_prime bases on the right (3') end
		    		
		    		   //now check that both words and primers of amplifia have same bases on 3' end
		    		   if ((w1 ^ w1a) != 0) continue; 
		    		   if ((w2 ^ w2a) != 0) continue;
		    	   }
		    	   
		    	   

		    	   // Look for the new pair in already seen pairs

		    	   pcurrent = insertpair(current,pairs);


		    	   if (seqdb[seqid]->isexample)

				   {
		    		   //pcurrent->inexample++;
		    		   pcurrent->sumd+=distance;
		    		   pcurrent->amplifiacount++;

		    		   if ((pcurrent->maxd==DMAX) || (distance > pcurrent->maxd))
			    		   pcurrent->maxd = distance;

			    	   if (distance < pcurrent->mind)
			    		   pcurrent->mind = distance;
		    	   }
		    	   //else
		    		 //  pcurrent->outexample++;

		    	   //for each pair we save current sequence id in the pair
		    	   //when we see this pair for the first time in currnet sequence
		    	   //because we want to increment inexample & outexample count
		    	   //only once for one sequence
				  if (pcurrent->curseqid != (seqid+1))
				   {
						if (seqdb[seqid]->isexample)
							pcurrent->inexample++;
		    	        else
		    		        pcurrent->outexample++;

						if (pcurrent->curseqid != 0)
							pcurrent->curseqid = seqid+1;
					}

					/*if ((pcurrent->outexample+pcurrent->inexample)==0)
					{
						fprintf(stderr,"pcurrent->outexample+pcurrent->inexample=0!\n");
						exit(0);
					}*/

		    	   if (pcurrent->curseqid == 0)//((pcurrent->outexample+pcurrent->inexample)==1)
		    	   {
					   pcurrent->curseqid = seqid+1;
					   paircount++;
					   pcurrent->pcr.ampslot=200;
					   pcurrent->pcr.ampcount=0;
					   pcurrent->pcr.amplifias = ECOMALLOC(sizeof(amplifia_t)*pcurrent->pcr.ampslot,
							                               "Cannot allocate amplifia table");
		    	   }
		    	   else
		    	   {
		    		   if (pcurrent->pcr.ampslot==pcurrent->pcr.ampcount)
		    		   {
		    			   pcurrent->pcr.ampslot+=200;
		    			   pcurrent->pcr.amplifias = ECOREALLOC(pcurrent->pcr.amplifias,
															    sizeof(amplifia_t)*pcurrent->pcr.ampslot,
		    			   							            "Cannot allocate amplifia table");
		    		   }
		    	   }

		    	   if (seqid==options->refseqid)
		    		   pcurrent->refsequence=seqid;
		    	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].length=distance;
		    	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].sequence=seqdb[seqid];
		    	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].strand=strand;
		    	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].begin=matches[i].position + options->primer_length;
		    	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].end=  matches[j].position - 1;

		    	   if (strand)
		    	   	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].amplifia=  seqdb[seqid]->SQ + matches[i].position + options->primer_length;
		    	   else
		    	   	   pcurrent->pcr.amplifias[pcurrent->pcr.ampcount].amplifia=  seqdb[seqid]->SQ + matches[j].position - 1 ;

		    	   
		    	   /*strncpy (prmr, seqdb[seqid]->SQ + matches[i].position, options->primer_length);
		    	   mtemp = nparam_CalcSelfTM (options->pnparm, prmr, options->primer_length) - 273.0;
		    	   if (mtemp < pcurrent->p1mintemp)
		    		   pcurrent->p1mintemp = mtemp;
		    	   //fprintf (stderr, "prmr1: %s\n", seqdb[seqid]->SQ);
		    	   strncpy (prmr, seqdb[seqid]->SQ + matches[j].position, options->primer_length);
		    	   mtemp = nparam_CalcSelfTM (options->pnparm, prmr, options->primer_length) - 273.0;
		    	   if (mtemp < pcurrent->p2mintemp)
		    		   pcurrent->p2mintemp = mtemp;
		    	   //fprintf (stderr, "prmr2: %s\n", prmr);

		    	   if (pcurrent->p1temp == 100)
		    		   pcurrent->p1temp = nparam_CalcSelfTM (options->pnparm, ecoUnhashWord(pcurrent->p1->word, options->primer_length), 0) - 273.0;
		    	   if (pcurrent->p2temp == 100)
		    		   pcurrent->p2temp = nparam_CalcSelfTM (options->pnparm, ecoUnhashWord(pcurrent->p2->word, options->primer_length), 0) - 273.0;
		    	   */
		    	   pcurrent->pcr.ampcount++;
//		    	   fprintf(stderr,"%c%c W1 : %s   direct : %c",
//		    			   "bG"[(int)pcurrent->p1->good],
//		    			   "bG"[(int)pcurrent->p2->good],
//		    			   ecoUnhashWord(pcurrent->p1->word, options->primer_length),
//		    			   "><"[(int)pcurrent->asdirect1]
//		    			   );
//
//		    	   fprintf(stderr,"   W2 : %s   direct : %c distance : %d (min/max/avg : %d/%d/%f) in/out: %d/%d %c (%d pairs)\n",
//		    			   ecoUnhashWord(pcurrent->p2->word, options->primer_length),
//		    			   "><"[(int)pcurrent->asdirect2],
//		    			   distance,
//		    			   pcurrent->mind,pcurrent->maxd,
//		    			   (pcurrent->inexample) ? (float)pcurrent->sumd/pcurrent->inexample:0.0,
//		    			   pcurrent->inexample,pcurrent->outexample,
//		    			   " N"[(pcurrent->outexample+pcurrent->inexample)==1],
//		    			   paircount
//
//		    			   );
//

		       }
		}
		 }
	}
   pairs->count=paircount;

}
