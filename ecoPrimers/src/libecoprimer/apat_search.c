/* ==================================================== */
/*      Copyright (c) Atelier de BioInformatique        */
/*      Dec. 94                                         */
/*      File: apat_search.c                             */
/*      Purpose: recherche du pattern                   */
/*               algorithme de Baeza-Yates/Gonnet       */
/*                             Manber (agrep)           */
/*      History:                                        */
/*      07/12/94 : <MFS>   first version                */
/*      28/12/94 : <Gloup> revised version              */
/*      14/05/99 : <Gloup> last revision                */
/* ==================================================== */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libstki.h"
#include "apat.h"

#define POP             PopiOut
#define PUSH(s,v)       PushiIn(&(s),(v))
#define TOPCURS         CursiToTop
#define DOWNREAD        ReadiDown

#define KRONECK(x, msk) ((~x & msk) ? 0 : 1)
#define MIN(x, y)       ((x) < (y)  ? (x) : (y))


/* -------------------------------------------- */
/* Baeza-Yates/Manber algorithm                 */
/* NoError                                      */
/* -------------------------------------------- */
int32_t ManberNoErr(pecoseq_t pseq,ppattern_t pat,
        ppatternParam_t param,
        StackiPtr stkpos)
{
        int32_t     pos;
        uint32_t    smask, r;
        uint8_t     *data;
        int32_t    end;

        end = (size_t)(pseq->SQ_length);

        if (param->circular)
        	end+=param->patlen - 1;


                                        /* create local masks   */
        smask = r = 0x1L << param->patlen;
                                        /* init. scan           */
        data   = (uint8_t*)(pseq->SQ);

                                        /* loop on text data    */
        for (pos = 0 ; pos < end ; pos++,data++) {
        	if (pos==pseq->SQ_length)
        		data=(uint8_t*)(pseq->SQ);

        	if (*data < 4)
        		r = (r >> 1) & pat[*data];
        	else
        		r=0;

            if (r & 0x1L) {
                PUSH(stkpos, pos - param->patlen + 1);
            }

            r |= smask;
        }
        return stkpos->top;  /* aka # of hits        */
}

/* -------------------------------------------- */
/* Baeza-Yates/Manber algorithm                 */
/* Substitution only                            */
/*                                              */
/* Note : r array is stored as :                */
/*    0 0 r(0,j) r(0,j+1) r(1,j) r(1,j+1) ...   */
/*                                              */
/* -------------------------------------------- */
int32_t ManberSub(pecoseq_t pseq,ppattern_t pat,
        ppatternParam_t param,
        StackiPtr stkpos)
{
        int       e, found;
        int32_t     pos;
        uint32_t    smask, cmask, sindx;
        uint32_t    *pr, r[2 * MAX_PAT_ERR + 2];
        uint8_t     *data;
        int32_t    end;

        end = (size_t)(pseq->SQ_length);

        if (param->circular)
        	end+=param->patlen - 1;

                                        /* create local masks   */
        r[0] = r[1] = 0x0;

        cmask = smask = 0x1L << param->patlen;

        for (e = 0, pr = r + 3 ; e <= param->maxerr ; e++, pr += 2)
                *pr = cmask;

        cmask = ~ param->omask;
                                        /* init. scan           */
        data   = (uint8_t*)(pseq->SQ);

                                        /* loop on text data    */

        for (pos = 0 ; pos < end ; pos++,data++) {
        	if (pos==pseq->SQ_length)
        		data=(uint8_t*)(pseq->SQ);

            sindx  = (*data==4) ? 0:pat[*data];

            for (e = found = 0, pr = r ; e <= param->maxerr ; e++, pr += 2) {

                pr[2]  = pr[3] | smask;

                pr[3]  =   ((pr[0] >> 1) & cmask)       /* sub   */
                         | ((pr[2] >> 1) & sindx);      /* ident */

                if (pr[3] & 0x1L) {                     /* found */
                    if (! found)  {
                       PUSH(stkpos, pos - param->patlen + 1);
                    }
                    found++;
                }
            }
        }

        return stkpos->top;  /* aka # of hits        */
}


/* -------------------------------------------- */
/* Baeza-Yates/Manber algorithm                 */
/* API call to previous functions               */
/* -------------------------------------------- */
int32_t ManberAll(pecoseq_t pseq,ppattern_t pat,
        ppatternParam_t param,
        StackiPtr stkpos)
{
        if (param->maxerr == 0)
       	  return ManberNoErr(pseq,
		                     pat, param,
		                     stkpos);
        else
          return ManberSub(pseq,
                             pat, param,
                             stkpos);
}

