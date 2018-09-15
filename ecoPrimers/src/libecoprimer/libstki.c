/* ==================================================== */
/*      Copyright (c) Atelier de BioInformatique        */
/*      Mar. 92                                         */
/*      File: libstki.c                                 */
/*      Purpose: A library to deal with 'stacks' of     */
/*               integers                               */
/*      Note: 'stacks' are dynamic (i.e. size is        */
/*            automatically readjusted when needed)     */
/*      History:                                        */
/*      00/03/92 : <Gloup> first draft                  */
/*      15/08/93 : <Gloup> revised version              */
/*      14/05/99 : <Gloup> last revision                */
/* ==================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libstki.h"
#include "ecoprimer.h"


/* ============================ */
/* Constantes et Macros locales */
/* ============================ */

#define ExpandStack(stkh) ResizeStacki((stkh), (*stkh)->size << 1)

#define ShrinkStack(stkh) ResizeStacki((stkh), (*stkh)->size >> 1)


static int16_t sStkiLastError = kStkiNoErr;

/* -------------------------------------------- */
/* gestion des erreurs                          */
/* get/reset erreur flag                        */
/*                                              */
/* @function: StkiError                         */
/* -------------------------------------------- */

int16_t StkiError(bool_t reset)
{
        int16_t err;

        err = sStkiLastError;

        if (reset)
           sStkiLastError = kStkiNoErr;

        return err;

} /* end of StkiError */

/* -------------------------------------------- */
/* creation d'un stack                          */
/*                                              */
/* @function: NewStacki                         */
/* -------------------------------------------- */

StackiPtr NewStacki(int32_t size)
{
        StackiPtr stki;

        if (! (stki = NEW(Stacki)))
                return NULL;

        stki->size    = size;
        stki->top     = 0;
        stki->cursor  = 0;

        if ( ! (stki->val = NEWN(int32_t, size))) {
            sStkiLastError = kStkiMemErr;
            return FreeStacki(stki);
        }

        return stki;

} /* end of NewStacki */


/* -------------------------------------------- */
/* liberation d'un stack                        */
/*                                              */
/* @function: FreeStacki                        */
/* -------------------------------------------- */

StackiPtr FreeStacki(StackiPtr stki)
{
        if (stki) {
            if (stki->val)
                ECOFREE(stki->val,"Free stack values");
            ECOFREE(stki,"Free stack");
        }

        return NULL;

} /* end of FreeStacki */

/* -------------------------------------------- */
/* creation d'un vecteur de stacks              */
/*                                              */
/* @function: NewStackiVector                   */
/* -------------------------------------------- */

StackiHdle NewStackiVector(int32_t vectSize, int32_t stackSize)
{
        int32_t           i;
        StackiHdle      stkh;

        if (! (stkh = NEWN(StackiPtr, vectSize))) {
            sStkiLastError = kStkiMemErr;
            return NULL;
        }

        for (i = 0 ; i < vectSize ; i++)
            if (! (stkh[i] = NewStacki(stackSize)))
                return FreeStackiVector(stkh, i);

        return stkh;

} /* end of NewStackiVector */


/* -------------------------------------------- */
/* liberation d'un vecteur de stacks            */
/*                                              */
/* @function: FreeStackiVector                  */
/* -------------------------------------------- */

StackiHdle FreeStackiVector(StackiHdle stkh, int32_t vectSize)
{
        int32_t   i;

        if (stkh) {
            for (i = 0 ; i < vectSize ; i++)
                (void) FreeStacki(stkh[i]);
            ECOFREE(stkh,"Free stack vector");
        }

        return NULL;

} /* end of FreeStackiVector */

/* -------------------------------------------- */
/* resize d'un stack                            */
/*                                              */
/* @function: ResizeStacki                      */
/* -------------------------------------------- */

int32_t ResizeStacki(StackiHdle stkh, int32_t size)
{
        int32_t resize = 0;               /* assume error         */
        int32_t *val;

        if ((val = ECOREALLOC((*stkh)->val, size * sizeof(int32_t),"Cannot reallocate stack values"))) {
            (*stkh)->size = resize = size;
            (*stkh)->val = val;
        }

        if (! resize)
            sStkiLastError = kStkiMemErr;

        return resize;

} /* end of ResizeStacki */

/* -------------------------------------------- */
/* empilage(/lement)                            */
/*                                              */
/* @function: PushiIn                           */
/* -------------------------------------------- */

bool_t PushiIn(StackiHdle stkh, int32_t val)
{
        if (((*stkh)->top >= (*stkh)->size) && (! ExpandStack(stkh)))
            return FALSE;

        (*stkh)->val[((*stkh)->top)++] = val;

        return TRUE;

} /* end of PushiIn */

/* -------------------------------------------- */
/* depilage(/lement)                            */
/*                                              */
/* @function: PopiOut                           */
/* -------------------------------------------- */

bool_t PopiOut(StackiHdle stkh, int32_t *val)
{
        if ((*stkh)->top <= 0)
            return FALSE;

        *val = (*stkh)->val[--((*stkh)->top)];

        if (    ((*stkh)->top < ((*stkh)->size >> 1))
             && ((*stkh)->top > kMinStackiSize))

            (void) ShrinkStack(stkh);

        return TRUE;

} /* end of PopiOut */

/* -------------------------------------------- */
/* lecture descendante                          */
/*                                              */
/* @function: ReadiDown                         */
/* -------------------------------------------- */

bool_t ReadiDown(StackiPtr stki, int32_t *val)
{
        if (stki->cursor <= 0)
            return FALSE;

        *val = stki->val[--(stki->cursor)];

        return TRUE;

} /* end of ReadiDown */

/* -------------------------------------------- */
/* lecture ascendante                           */
/*                                              */
/* @function: ReadiUp                           */
/* -------------------------------------------- */

bool_t ReadiUp(StackiPtr stki, int32_t *val)
{
        if (stki->cursor >= stki->top)
            return FALSE;

        *val = stki->val[(stki->cursor)++];

        return TRUE;

} /* end of ReadiUp */

/* -------------------------------------------- */
/* remontee/descente du curseur                 */
/*                                              */
/* @function: CursiToTop                        */
/* @function: CursiToBottom                     */
/* -------------------------------------------- */

void CursiToTop(StackiPtr stki)
{
        stki->cursor = stki->top;

} /* end of CursiToTop */

void CursiToBottom(stki)
        StackiPtr stki;
{
        stki->cursor = 0;

} /* end of CursiToBottom */

/* -------------------------------------------- */
/* echange des valeurs cursor <-> (top - 1)     */
/*                                              */
/* @function: CursiSwap                         */
/* -------------------------------------------- */

void CursiSwap(StackiPtr stki)
{
        int32_t   tmp;

        if ((stki->top <= 0) || (stki->cursor < 0))
            return;

        tmp = stki->val[stki->cursor];
        stki->val[stki->cursor] = stki->val[stki->top - 1];
        stki->val[stki->top - 1] = tmp;

} /* end of CursiSwap */

/* -------------------------------------------- */
/* Recherche d'une valeur en stack a partir du  */
/* curseur courant en descendant.               */
/* on laisse le curseur a l'endroit trouve      */
/*                                              */
/* @function: SearchDownStacki                  */
/* -------------------------------------------- */

bool_t SearchDownStacki(StackiPtr stki, int32_t sval)
{
        int32_t   val;
        bool_t    more;

        while ((more = ReadiDown(stki, &val)))
            if (val == sval)
                break;

        return more;

} /* end of SearchDownStacki */

/* -------------------------------------------- */
/* Recherche dichotomique d'une valeur en stack */
/* le stack est suppose trie par valeurs        */
/* croissantes.                                 */
/* on place le curseur a l'endroit trouve       */
/*                                              */
/* @function: BinSearchStacki                   */
/* -------------------------------------------- */

bool_t BinSearchStacki(StackiPtr stki, int32_t sval)
{
        int32_t   midd, low, high, span;

        low  = 0;
        high = stki->top - 1;

        while (high >= low) {

            midd = (high + low) / 2;

            span = stki->val[midd] - sval;

            if (span == 0) {
                stki->cursor = midd;
                return TRUE;
            }

            if (span > 0)
                high = midd - 1;
            else
                low  = midd + 1;
        }

        return FALSE;

} /* end of BinSearchStacki */

/* -------------------------------------------- */
/* teste l'egalite *physique* de deux stacks    */
/*                                              */
/* @function: SameStacki                        */
/* -------------------------------------------- */

bool_t SameStacki(StackiPtr stki1, StackiPtr stki2)
{
        if (stki1->top != stki2->top)
            return FALSE;

        return ((memcmp(stki1->val, stki2->val,
                        stki1->top * sizeof(int32_t)) == 0) ? TRUE : FALSE);

} /* end of SameStacki */


/* -------------------------------------------- */
/* inverse l'ordre des elements dans un stack   */
/*                                              */
/* @function: ReverseStacki                     */
/* -------------------------------------------- */

bool_t ReverseStacki(StackiPtr stki)
{
        int32_t   *t, *b, swp;

        if (stki->top <= 0)
            return FALSE;

        b = stki->val;
        t = b + stki->top - 1;

        while (t > b) {
             swp  = *t;
             *t-- = *b;
             *b++ = swp;
        }

        return TRUE;

} /* end of ReverseStacki */

