/*
 * ahocorasick.h
 *
 *  Created on: 26 march 2011
 *      Author: tiayyba
 */
#include <inttypes.h>
#include "hashencoder.h"
#include "ahocorasick.h"

void ahoc_graphKeywordTree (aho_state *root);
aho_state *groot = NULL; //just for graph testing

#define BASEATINDEX(w, l, i) (uint8_t)((((w)&(0x3LLU<<(((l)-(i))*2)))>>(((l)-(i))*2)) & 0x3LLU)

void ahoc_addOutputElement (aho_state *node, bool_t isdirect, uint32_t idx)
{
	if (!node) return;
	if (node->output.count == 0)
		node->output.out_set = ECOMALLOC(sizeof(aho_output),
                                      "Cannot allocate memory for aho-corasick state output element");
	else
		node->output.out_set = ECOREALLOC(node->output.out_set, (node->output.count+1)*sizeof(aho_output),
                                      "Cannot allocate memory for aho-corasick state output element");
	node->output.out_set[node->output.count].wordidx = idx;
	node->output.out_set[node->output.count].isdirect = isdirect;
	node->output.count++;
}

//is the passed output element in the set
bool_t ahoc_isOutputIn (aho_state *node, aho_output ot)
{
	uint32_t i;

	for (i=0; i<node->output.count; i++)
		if (node->output.out_set[i].isdirect == ot.isdirect && node->output.out_set[i].wordidx == ot.wordidx) return TRUE;
	return FALSE;
}

//take union of output of the two nodes and put in node1
void ahoc_unionOutputElements (aho_state *node1, aho_state *node2)
{
	uint32_t i;

	for (i=0; i<node2->output.count; i++)
		if (ahoc_isOutputIn (node1, node2->output.out_set[i]) == FALSE)
			ahoc_addOutputElement (node1, node2->output.out_set[i].isdirect, node2->output.out_set[i].wordidx);
}

void ahoc_addKeyword (aho_state *root, word_t w, bool_t isdirect, uint32_t idx, poptions_t options)
{
	uint32_t i;
	aho_state *nextnode = root;
	uint8_t basecode;
	static uint32_t state_id = 0;
	
	//fprintf (stderr, "%s\n", ecoUnhashWord(w, options->primer_length));
	for (i=1; i<=options->primer_length; i++)
	{
		basecode = BASEATINDEX (w, options->primer_length, i);
		//fprintf (stderr, "%d", basecode);
		if (nextnode->next[basecode] == NULL)
		{
			//add new state
			nextnode->next[basecode] =  ECOMALLOC(sizeof(aho_state),
                                      "Cannot allocate memory for aho-corasick state");
			nextnode = nextnode->next[basecode];
			//initialize state
			nextnode->id = ++state_id;
			nextnode->next[0]=nextnode->next[1]=nextnode->next[2]=nextnode->next[3]=NULL;
			nextnode->fail = NULL;
			nextnode->output.count = 0;
		}
		else
			nextnode = nextnode->next[basecode];
	}
	//fprintf (stderr, "\n", basecode);
	//new pattern addess so add node ouptup element
	ahoc_addOutputElement (nextnode, isdirect, idx);
}

void ahoc_buildKeywordTree (aho_state *root, pwordcount_t words, poptions_t options)
{
	uint32_t i;
	if (!root) return;

	//init root
	root->id = 0;
	root->next[0]=root->next[1]=root->next[2]=root->next[3]=NULL;
	root->fail = NULL;
	root->output.count = 0;

	//now add each word as a pattern in the keyword tree
	for (i=0; i<words->size; i++)
	{
		//add direct word
		word_t w=WORD(words->words[i]);
		ahoc_addKeyword (root, w, TRUE, i, options);

		//add reverse word
		w=ecoComplementWord(w,options->primer_length);
		ahoc_addKeyword (root, w, FALSE, i, options);
	}

	//loop on root if some base has no out going edge from roots
	for (i=0; i<4; i++)
		if (root->next[i] == NULL)
			root->next[i] = root;
}

void ahoc_enqueue (aho_queue *ahoqueue, aho_state *node)
{
	queue_node *q;
	if (node == NULL) return;

	q = ECOMALLOC(sizeof(queue_node),
            "Cannot allocate memory for aho-corasick queue node");
	q->state_node = node;
	q->next = NULL;

	if (ahoqueue->first == NULL)
	{
		ahoqueue->first = q;
		ahoqueue->last = q;
	}
	else
	{
		ahoqueue->last->next = q;
		ahoqueue->last = q;
	}
}

aho_state *ahoc_dequeue (aho_queue *ahoqueue)
{
	aho_state *node = NULL;
	queue_node *q;

	if (ahoqueue->first == NULL) return node;
	q = ahoqueue->first;
	ahoqueue->first = q->next;

	node = q->state_node;
	ECOFREE (q, "Cannot free memory for aho-corasick queue node");
	return node;
}

//set fail links and output sets for the keyword tree
void ahoc_updateForFailAndOutput (aho_state *root)
{
	int32_t i;
	aho_queue Q;
	aho_state *node_r;
	aho_state *node_u;
	aho_state *node_v;

	//empty queue
	Q.first = NULL;
	Q.last = NULL;

	//for us alphabet has 4 elements, A=0, C=1, G=2 and T=3
	for (i=0; i<4; i++)
	{
		if (root->next[i] != root && root->next[i] != NULL)
		{
			root->next[i]->fail = root;
			ahoc_enqueue (&Q, root->next[i]);
		}
	}

	//while queue not empty
	while (Q.first != NULL)
	{
		node_r = ahoc_dequeue (&Q);
		for (i=0; i<4; i++)
		{
			if (node_r->next[i] != NULL)
			{
				node_u = node_r->next[i];
				ahoc_enqueue (&Q, node_u);
				node_v = node_r->fail;
				while (node_v->next[i] == NULL) 
					node_v = node_v->fail;
				node_u->fail = node_v->next[i];
				ahoc_unionOutputElements (node_u, node_u->fail);
			}
		}
	}
}

void ahoc_freeKeywordTree (aho_state *node)
{
	int i;
	for (i=0; i<4; i++)
		if (node->next[i])
			ahoc_freeKeywordTree (node->next[i]);
	if (node->output.count > 0)
		ECOFREE (node->output.out_set, "Free failed for node output");
	ECOFREE (node, "Free failed for node");
}

pprimercount_t ahoc_lookforStrictPrimers (pecodnadb_t database, uint32_t seqdbsize,uint32_t exampleCount,
		                     pwordcount_t words,poptions_t options)
{
	aho_state automaton_root;
	aho_state *curr_state;
	//uint32_t  inSequenceQuorum;
	uint32_t  outSequenceQuorum;
	pprimer_t  data;
	pprimercount_t primers;
	uint32_t   i, j, k;
	int32_t pos;
	uint32_t lmax;
	char *base;
	int8_t code;
	uint32_t   goodPrimers=0;
	static int iii=0;


	//inSequenceQuorum = (uint32_t)floor((float)exampleCount * options->sensitivity_quorum);
	outSequenceQuorum = (uint32_t)floor((float)(seqdbsize-exampleCount) * options->false_positive_quorum);

	//fprintf(stderr,"  Primers should be at least present in %d/%d example sequences\n",inSequenceQuorum,exampleCount);
	fprintf(stderr,"  Primers should not be present in more than %d/%d counterexample sequences\n",outSequenceQuorum,(seqdbsize-exampleCount));

	data = ECOMALLOC(words->size * sizeof(primer_t),
			         "Cannot allocate memory for fuzzy matching results");
	for (i=0; i < words->size; i++)
	{
		data[i].word=WORD(words->words[i]);
		data[i].inexample = 0;
		data[i].outexample= 0;

		data[i].directCount=ECOMALLOC(seqdbsize * sizeof(uint32_t),
		              	"Cannot allocate memory for primer position");
		data[i].directPos = ECOMALLOC(seqdbsize * sizeof(poslist_t),
				"Cannot allocate memory for primer position");
		data[i].reverseCount=ECOMALLOC(seqdbsize * sizeof(uint32_t),
		 		"Cannot allocate memory for primer position");
		data[i].reversePos = ECOMALLOC(seqdbsize * sizeof(poslist_t),
				"Cannot allocate memory for primer position");
	}

	//build keywords automaton
	ahoc_buildKeywordTree (&automaton_root, words, options);
	//set fail links and output sets
	ahoc_updateForFailAndOutput (&automaton_root);

	//debug; print keywordtree in a gv file
	//ahoc_graphKeywordTree (&automaton_root);

	//loop on each sequence for its each base and find words
	for (i=0; i < seqdbsize; i++)
	{
		if(database[i]->SQ_length <= options->primer_length) continue;

		lmax = database[i]->SQ_length;
		if (!options->circular)
			lmax += options->primer_length-1;
		curr_state = &automaton_root;

		for (j=0,base=database[i]->SQ; j<lmax; j++,base++)
		{
			if (i==(uint32_t)database[i]->SQ_length) base=database[i]->SQ;

			//code = encoder[(*base) - 'A'];
			code = *base;
			//if (iii++ < 30)
			//	fprintf (stderr, "%d:%d,", *base, code);
			if (code < 0 || code > 3)
			{
				//if error char, start from root for next character
				//+forget any incomplete words
				curr_state = &automaton_root;
				continue;
			}
			while (curr_state->next[code] == NULL) curr_state = curr_state->fail;
			curr_state = curr_state->next[code];

			//start position of primer is options->primer_length-1 chars back
			pos = j-options->primer_length+1;
			if (pos < 0) pos = database[i]->SQ_length+pos;

			//set output, if there is some output on this state then
			//+all words in the output set complete here, so increment their
			//+found properties for current sequence
			for (k=0; k<curr_state->output.count; k++)
			{
				if (curr_state->output.out_set[k].isdirect)
					data[curr_state->output.out_set[k].wordidx].directCount[i]++;
				else
					data[curr_state->output.out_set[k].wordidx].reverseCount[i]++;

				if (options->no_multi_match)
				{
					if ((data[curr_state->output.out_set[k].wordidx].directCount[i] + 
					  data[curr_state->output.out_set[k].wordidx].reverseCount[i]) > 1)
						//since multimach not allowd, set an indication on 1st seq position that
						//+ a multimatch was found, so that this word will be filtered out
						//+ and because of first postion we wont have to search the whole array
						//+ to find if it voilated nomultimatch constraint for some seq
						data[curr_state->output.out_set[k].wordidx].directCount[0] = 2;
					else
					{
						if (curr_state->output.out_set[k].isdirect)
							//direct word found on jth position of ith sequence
							data[curr_state->output.out_set[k].wordidx].directPos[i].value = (uint32_t)pos;
						else
							//reverse word found on jth position of ith sequence
							data[curr_state->output.out_set[k].wordidx].reversePos[i].value = (uint32_t)pos;
					}
				}
				else
				{
					//okay multi match allowed
					if (curr_state->output.out_set[k].isdirect)
					{
						if (data[curr_state->output.out_set[k].wordidx].directCount[i] == 1)
							data[curr_state->output.out_set[k].wordidx].directPos[i].value = (uint32_t)pos;
						else
						{
							//need to create or extend the positions list
							if (data[curr_state->output.out_set[k].wordidx].directCount[i] == 2)
							{
								//for second element, first was put in .value, so dont forget to copy that in the array too
								data[curr_state->output.out_set[k].wordidx].directPos[i].pointer = ECOMALLOC(2 * sizeof(uint32_t),
		                      													"Cannot allocate memory for primer position");
								data[curr_state->output.out_set[k].wordidx].directPos[i].pointer[0] = data[curr_state->output.out_set[k].wordidx].directPos[i].value;
								data[curr_state->output.out_set[k].wordidx].directPos[i].pointer[1] = (uint32_t)pos;
							}
							else
							{
								//for third or greater element
								data[curr_state->output.out_set[k].wordidx].directPos[i].pointer = ECOREALLOC(data[curr_state->output.out_set[k].wordidx].directPos[i].pointer, 
																	data[curr_state->output.out_set[k].wordidx].directCount[i] * sizeof(uint32_t),
		                      													"Cannot allocate memory for primer position");
								data[curr_state->output.out_set[k].wordidx].directPos[i].pointer[data[curr_state->output.out_set[k].wordidx].directCount[i]-1] = (uint32_t)pos;
							}
						}
					}
					else
					{
						if (data[curr_state->output.out_set[k].wordidx].reverseCount[i] == 1)
							data[curr_state->output.out_set[k].wordidx].reversePos[i].value = (uint32_t)pos;
						else
						{
							//need to create or extend the positions list
							if (data[curr_state->output.out_set[k].wordidx].reverseCount[i] == 2)
							{
								//for second element, first was put in .value, so dont forget to copy that in the array too
								data[curr_state->output.out_set[k].wordidx].reversePos[i].pointer = ECOMALLOC(2 * sizeof(uint32_t),
		                      													"Cannot allocate memory for primer position");
								data[curr_state->output.out_set[k].wordidx].reversePos[i].pointer[0] = data[curr_state->output.out_set[k].wordidx].reversePos[i].value;
								data[curr_state->output.out_set[k].wordidx].reversePos[i].pointer[1] = (uint32_t)pos;
							}
							else
							{
								//for third or greater element
								data[curr_state->output.out_set[k].wordidx].reversePos[i].pointer = ECOREALLOC(data[curr_state->output.out_set[k].wordidx].reversePos[i].pointer, 
																	data[curr_state->output.out_set[k].wordidx].reverseCount[i] * sizeof(uint32_t),
		                      													"Cannot allocate memory for primer position");
								data[curr_state->output.out_set[k].wordidx].reversePos[i].pointer[data[curr_state->output.out_set[k].wordidx].reverseCount[i]-1] = (uint32_t)pos;
							}
						}
					}
				}
				//dont forget to increment inexample or outexample count, but only once for a sequence
				if ((data[curr_state->output.out_set[k].wordidx].directCount[i] + 
				  data[curr_state->output.out_set[k].wordidx].reverseCount[i]) == 1)
				{
					if (database[i]->isexample)
						data[curr_state->output.out_set[k].wordidx].inexample++;
					else
						data[curr_state->output.out_set[k].wordidx].outexample++;
				}
			}
		}
	}

	//Only thing that remains is to remove the failed words
	for (i=0,j=0; i<words->size; i++)
	{
		fprintf(stderr,"Primers %5d/%lld analyzed  => sequence : %s in %d example and %d counterexample sequences       \r",
				i+1,words->size,ecoUnhashWord(data[i].word,options->primer_length),
				data[i].inexample,data[i].outexample);

		//if (data[i].inexample < inSequenceQuorum || (data[i].directCount[0] == 2 && options->no_multi_match))
		if (data[i].directCount[0] == 2 && options->no_multi_match)
		{
			//bad word, delete from the array
			for (k=0; k<seqdbsize; k++)
			{
				if (data[i].directCount[k] > 1)
					ECOFREE (data[i].directPos[k].pointer, "Cannot free position pointer.");
				if (data[i].reverseCount[k] > 1)
					ECOFREE (data[i].reversePos[k].pointer, "Cannot free position pointer.");
			}
			ECOFREE (data[i].directCount, "Cannot free position pointer.");
			ECOFREE (data[i].directPos, "Cannot free position pointer.");
			ECOFREE (data[i].reverseCount, "Cannot free position pointer.");
			ECOFREE (data[i].reversePos, "Cannot free position pointer.");
		}
		else
		{
			//data[i].good = data[i].inexample >= inSequenceQuorum && data[i].outexample <= outSequenceQuorum;
			data[i].good = data[i].outexample <= outSequenceQuorum;
			goodPrimers+=data[i].good? 1:0;
			if (j < i)
				data[j] = data[i];
			j++;
		}
	}
	fprintf(stderr,"\n\nOn %lld analyzed primers %d respect quorum conditions\n",words->size,goodPrimers);
	fprintf(stderr,"Conserved primers for further analysis : %d/%lld\n",j,words->size);

	primers = ECOMALLOC(sizeof(primercount_t),"Cannot allocate memory for primer table");
	primers->primers=ECOREALLOC(data,
			                    j * sizeof(primer_t),
								"Cannot reallocate memory for fuzzy matching results");
	primers->size=j;

	//free memory of keyword table
	for (i=0; i<4; i++)
		if (automaton_root.next[i] != &automaton_root)
			ahoc_freeKeywordTree (automaton_root.next[i]);

	return primers;
}

void ahoc_graphPrintNodesInfo (aho_state *node, FILE* gfile)
{
	uint32_t i;
	fprintf (gfile, "\"%d\"[\n", node->id);
	fprintf (gfile, "label=\"%d\\n", node->id);
	for (i=0; i<node->output.count; i++)
		fprintf (gfile, "%d%c,", node->output.out_set[i].wordidx, node->output.out_set[i].isdirect?'d':'r');
	fprintf (gfile, "\"\n];\n");

	for (i=0; i<4; i++)
		if (node->next[i] != NULL && node->next[i] != node)
			ahoc_graphPrintNodesInfo (node->next[i], gfile);
}

void ahoc_graphPrintNodesLinks (aho_state *node, FILE* gfile)
{
	uint32_t i;
	static int j=0;

	for (i=0; i<4; i++)
		if (node->next[i] != NULL && node->next[i] != node)
		{
			fprintf (gfile, "\"%d\" -> \"%d\" [\n", node->id, node->next[i]->id);
			fprintf (gfile, "label=\"%c\"\n];\n", "ACGT"[i]);
		}

	if (j++ < 40)
	if (node->fail != NULL && node->fail != groot)
	{
		fprintf (gfile, "\"%d\" -> \"%d\" [\n", node->id, node->fail->id);
		fprintf (gfile, "color= \"red\"\n];\n");
	}

	for (i=0; i<4; i++)
		if (node->next[i] != NULL && node->next[i] != node)
			ahoc_graphPrintNodesLinks (node->next[i], gfile);
}

void ahoc_graphKeywordTree (aho_state *root)
{
	FILE *gfile;

	groot=root;
	gfile = fopen ("keywordtree.gv", "w");
	fprintf (gfile, "digraph keywordtree {\n");
	ahoc_graphPrintNodesInfo (root, gfile);
	ahoc_graphPrintNodesLinks (root, gfile);
	fprintf (gfile, "}\n");
	fclose(gfile);
}

