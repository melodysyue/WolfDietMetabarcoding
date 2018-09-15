/*
 * debug.h
 *
 *  Created on: 12 nov. 2008
 *      Author: coissac
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include <stdio.h>


#ifdef DEBUG

#define DEBUG_LOG(message,...) { \
	  char *text; \
	  (void)asprintf(&text,(message),##__VA_ARGS__); \
      fprintf(stderr,"DEBUG %s (line %d) : %s\n",__FILE__,__LINE__,(text)); \
      free(text); \
	}

#else

#define DEBUG_LOG(message, ...)

#endif

#endif /* DEBUG_H_ */
