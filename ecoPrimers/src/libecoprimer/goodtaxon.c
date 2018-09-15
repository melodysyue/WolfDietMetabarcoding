/*
 * goodtaxon.c
 *
 *  Created on: 7 nov. 2008
 *      Author: coissac
 */


#include "ecoprimer.h"

int isGoodTaxon(ecotaxonomy_t *taxonomy,int32_t taxon,poptions_t options)
{
	int result;

	result=((options->r == 0) || (eco_is_taxid_included(taxonomy,
			                            options->restricted_taxid,
			                            options->r,
			                            taxonomy->taxons->taxon[taxon].taxid)
		   )) &&
		   ((options->e == 0) || !(eco_is_taxid_included(taxonomy,
				                        options->exception_taxid,
		   		                        options->e,
		   		                        taxonomy->taxons->taxon[taxon].taxid)
		   ));

	return result;
}

int isExampleTaxon(ecotaxonomy_t *taxonomy,int32_t taxon,poptions_t options)
{
	int result;

	result=( (options->r == 0) || (eco_is_taxid_included(taxonomy,
			                            options->restricted_taxid,
			                            options->r,
			                            taxonomy->taxons->taxon[taxon].taxid)
		   )) &&
		   ((options->e == 0) || !(eco_is_taxid_included(taxonomy,
				                        options->exception_taxid,
		   		                        options->e,
		   		                        taxonomy->taxons->taxon[taxon].taxid)
		   ));

	return result;
}


int isCounterExampleTaxon(ecotaxonomy_t *taxonomy,int32_t taxon,poptions_t options)
{
	int result;

	result=((options->g != 0) && (eco_is_taxid_included(taxonomy,
								   options->ignored_taxid,
								   options->g,
								   taxonomy->taxons->taxon[taxon].taxid))
           ) || ((options->e != 0) && (eco_is_taxid_included(taxonomy,
				                        options->exception_taxid,
		   		                        options->e,
		   		                        taxonomy->taxons->taxon[taxon].taxid))
           );


	return result;
}
