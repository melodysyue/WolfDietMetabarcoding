/* ==================================================== */
/*      Copyright (c) Atelier de BioInformatique        */
/*      Dec. 94                                         */
/*      File: apat.h                                    */
/*      Purpose: pattern scan                           */
/*      History:                                        */
/*      28/12/94 : <Gloup> ascan first version          */
/*      14/05/99 : <Gloup> last revision                */
/* ==================================================== */


#ifndef H_apat
#define H_apat


#include "libstki.h"
#include "inttypes.h"
#include "../libecoPCR/ecoPCR.h"


/* ----------------------------------------------- */
/* constantes                                      */
/* ----------------------------------------------- */

#ifndef BUFSIZ
#define BUFSIZ          1024    /* io buffer size               */
#endif

#define MAX_NAME_LEN    BUFSIZ  /* max length of sequence name  */

#define ALPHA_LEN        4     /* alphabet length              */
                                /* *DO NOT* modify              */

#define MAX_PATTERN       4     /* max # of patterns            */
                                /* *DO NOT* modify              */

#define MAX_PAT_LEN      32     /* max pattern length           */
                                /* *DO NOT* modify              */

#define MAX_PAT_ERR      32     /* max # of errors              */
                                /* *DO NOT* modify              */

#define PATMASK 0x3ffffff       /* mask for 26 symbols          */
                                /* *DO NOT* modify              */

#define OBLIBIT 0x4000000       /* bit 27 to 1 -> oblig. pos    */
                                /* *DO NOT* modify              */

                                /* mask for position            */
#define ONEMASK 0x80000000      /* mask for highest position    */

                                /* masks for Levenhstein edit   */
#define OPER_IDT  0x00000000    /* identity                     */
#define OPER_INS  0x40000000    /* insertion                    */
#define OPER_DEL  0x80000000    /* deletion                     */
#define OPER_SUB  0xc0000000    /* substitution                 */

#define OPER_SHFT 30            /* <unused> shift               */

                                /* Levenhstein Opcodes          */
#define SOPER_IDT 0x0           /* identity                     */
#define SOPER_INS 0x1           /* insertion                    */
#define SOPER_DEL 0x2           /* deletion                     */
#define SOPER_SUB 0x3           /* substitution                 */

                                /* Levenhstein Opcodes masks    */
#define OPERMASK  0xc0000000    /* mask for Opcodes             */
#define NOPERMASK 0x3fffffff    /* negate of previous           */



/* ----------------------------------------------- */
/* data structures                                 */
/* ----------------------------------------------- */


typedef uint32_t pattern_t[ALPHA_LEN], *ppattern_t;

                                    /* -------------------- */
typedef struct {                    /* pattern              */
                                    /* -------------------- */
int      patlen;                    /* pattern length       */
int      maxerr;                    /* max # of errors      */
uint32_t omask;                     /* oblig. bits mask     */
bool_t   circular;                  /* is circular sequence */
} patternParam_t, *ppatternParam_t;


/* ----------------------------------------------- */
/* macros                                          */
/* ----------------------------------------------- */

#ifndef NEW
#define NEW(typ)                (typ*)malloc(sizeof(typ))
#define NEWN(typ, dim)          (typ*)malloc((unsigned long)(dim) * sizeof(typ))
#define REALLOC(typ, ptr, dim)  (typ*)realloc((void *) (ptr), (unsigned long)(dim) * sizeof(typ))
#define FREE(ptr)               free((void *) ptr)
#endif

/* ----------------------------------------------- */
/* prototypes                                      */
/* ----------------------------------------------- */

                                /* apat_search.c        */

int32_t ManberNoErr(pecoseq_t pseq,ppattern_t pat,
                    ppatternParam_t param,
                    StackiPtr stkpos);

int32_t ManberSub(pecoseq_t pseq,ppattern_t pat,
		          ppatternParam_t param,
		          StackiPtr stkpos);

int32_t ManberAll(pecoseq_t pseq,ppattern_t pat,
                  ppatternParam_t param,
                  StackiPtr stkpos);


#endif /* H_apat */

