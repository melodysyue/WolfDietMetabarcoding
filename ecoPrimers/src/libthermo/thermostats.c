#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "thermostats.h"

word_t extractSite(char* sequence, size_t begin, size_t length, bool_t strand)
{
	char *c;
	char *start;
	uint32_t l;
	word_t site = 0;

	start=sequence+begin;
	if (!strand)
		start+=length-1;


	for (c=start,
	     l=0;
		 l<length;
		 l++,
		 c+=(strand)? 1:-1)
		site = (site << 2) | ((strand)? (*c):(~*c)&3);

	return site;
}

void getThermoProperties (ppair_t* pairs, size_t count, poptions_t options)
{
	size_t i, j,k,l;
	uint32_t   bp1,bp2;
	uint32_t   ep1,ep2;
	word_t w1;
	word_t w2;
	bool_t strand;

	char *sq,*sq1,*sq2,*c;
	char prmrd[50];
	char prmrr[50];
	char sqsite[50];
	double mtemp;

	for (i = 0; i < count; i++)
	{
		w1 = pairs[i]->p1->word;
		w2 = pairs[i]->p2->word;

		if (!pairs[i]->asdirect1)
			w1=ecoComplementWord(w1,options->primer_length);

		if (!pairs[i]->asdirect2)
			w2=ecoComplementWord(w2,options->primer_length);

		strncpy(prmrd,ecoUnhashWord(w1, options->primer_length),options->primer_length);
		strncpy(prmrr,ecoUnhashWord(w2, options->primer_length),options->primer_length);
		prmrd[options->primer_length]=0;
		prmrr[options->primer_length]=0;
		pairs[i]->p1temp = nparam_CalcSelfTM (options->pnparm, prmrd, options->primer_length) - 273.0;
		pairs[i]->p2temp = nparam_CalcSelfTM (options->pnparm, prmrr, options->primer_length) - 273.0;
		pairs[i]->p1mintemp = 100;
		pairs[i]->p2mintemp = 100;

		for (j = 0; j < pairs[i]->pcr.ampcount; j++)
			if (pairs[i]->pcr.amplifias[j].sequence->isexample)
		{

			sq = pairs[i]->pcr.amplifias[j].sequence->SQ;
			strand = pairs[i]->pcr.amplifias[j].strand;
			bp1 = pairs[i]->pcr.amplifias[j].begin - options->primer_length;
			bp2 = pairs[i]->pcr.amplifias[j].end + 1;

			if (!strand)
			{
				uint32_t tmp;
				tmp=bp1;
				bp1=bp2;
				bp2=tmp;
			}

//			printf("%s : %s, %c",prmrd,
//					               ecoUnhashWord(extractSite(sq,bp1,options->primer_length,strand),options->primer_length),
//					               "rd"[strand]);
		    mtemp = nparam_CalcTwoTM(options->pnparm,
		    		                 prmrd,
		    		                 ecoUnhashWord(extractSite(sq,bp1,options->primer_length,strand),options->primer_length),
		    		                 options->primer_length) - 273.0;
//		    printf(" %4.2f  %4.2f\n",pairs[i]->p1temp,mtemp);
		    if (mtemp < pairs[i]->p1mintemp)
				    		pairs[i]->p1mintemp = mtemp;

//			printf("%s : %s, %c\n",prmrr,ecoUnhashWord(extractSite(sq,bp2,options->primer_length,!strand),options->primer_length),
//					"rd"[strand]);
//
	    	mtemp = nparam_CalcTwoTM(options->pnparm,
	    			                 prmrr,
	    			                 ecoUnhashWord(extractSite(sq,bp2,options->primer_length,!strand),options->primer_length),
	    			                 options->primer_length) - 273.0;
	    	if (mtemp < pairs[i]->p2mintemp)
	    		pairs[i]->p2mintemp = mtemp;
		}

		if (w2 < w1)
		{
			mtemp = pairs[i]->p1temp;
			pairs[i]->p1temp = pairs[i]->p2temp;
			pairs[i]->p2temp = mtemp;

			mtemp = pairs[i]->p1mintemp;
			pairs[i]->p1mintemp = pairs[i]->p2mintemp;
			pairs[i]->p2mintemp = mtemp;
		}

	}
}
