#include "ecoPCR.h"
#include <stdlib.h>

static int eco_log_malloc = 0;
static size_t eco_amount_malloc=0;
static size_t eco_chunk_malloc=0;

void    eco_trace_memory_allocation()
{
	eco_log_malloc=1;
}

void    eco_untrace_memory_allocation()
{
	eco_log_malloc=0;
}

void ecoMallocedMemory()
{
	//eco_amount_malloc;
}

void   *eco_malloc(int64_t chunksize,
                   const char *error_message,
                   const char *filename,
                   int32_t    line)
{
	void * chunk;

	chunk = calloc(1,chunksize);


	if (!chunk)
		ecoError(ECO_MEM_ERROR,error_message,filename,line);

	eco_chunk_malloc++;

	if (eco_log_malloc)
		fprintf(stderr,
			    "Memory segment located at %p of size %d is allocated (file : %s [%d])",
			    chunk,
			    chunksize,
			    filename,
			    line);

	return chunk;
}

void   *eco_realloc(void *chunk,
                    int64_t newsize,
                    const char *error_message,
                    const char *filename,
                    int32_t    line)
{
	void *newchunk;

	if (newsize == 0)
	{
		if (chunk)
			free(chunk);
		return NULL;
	}

	newchunk = realloc(chunk,newsize);

	if (!newchunk)
           {
        fprintf(stderr,"Requested memory : %d\n",newsize);
		ecoError(ECO_MEM_ERROR,error_message,filename,line);
           }
	if (!chunk)
		eco_chunk_malloc++;

	if (eco_log_malloc)
		fprintf(stderr,
			    "Old memory segment %p is reallocated at %p with a size of %d (file : %s [%d])",
			    chunk,
			    newchunk,
			    newsize,
			    filename,
			    line);

	return newchunk;
}

void    eco_free(void *chunk,
                 const char *error_message,
                 const char *filename,
                 int32_t    line)
{
	free(chunk);

	if (eco_log_malloc)
		fprintf(stderr,
			    "Memory segment %p is released => %s (file : %s [%d])",
			    chunk,
			    error_message,
			    filename,
			    line);

	eco_chunk_malloc--;
}
