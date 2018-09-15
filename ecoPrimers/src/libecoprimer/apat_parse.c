/* ==================================================== */
/*      Copyright (c) Atelier de BioInformatique        */
/*      Mar. 92                                         */
/*      File: apat_parse.c                              */
/*      Purpose: Codage du pattern                      */
/*      History:                                        */
/*      00/07/94 : <Gloup> first version (stanford)     */
/*      00/11/94 : <Gloup> revised for DNA/PROTEIN      */
/*      30/12/94 : <Gloup> modified EncodePattern       */
/*                         for manber search            */
/*      14/05/99 : <Gloup> indels added                 */
/* ==================================================== */

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "apat.h"
#include "ecoprimer.h"


 /* IUPAC Dna      */
static int32_t sDnaCode[]  =  {
		/* IUPAC */

		        0x00000001 /* A */, 0x0000000E /* B */, 0x00000002 /* C */,
		        0x0000000D /* D */, 0x00000000 /* E */, 0x00000000 /* F */,
		        0x00000004 /* G */, 0x0000000B /* H */, 0x00000000 /* I */,
		        0x00000000 /* J */, 0x0000000C /* K */, 0x00000000 /* L */,
		        0x00000003 /* M */, 0x0000000F /* N */, 0x00000000 /* O */,
		        0x00000000 /* P */, 0x00000000 /* Q */, 0x00000005 /* R */,
		        0x00000006 /* S */, 0x00000008 /* T */, 0x00000008 /* U */,
		        0x00000007 /* V */, 0x00000009 /* W */, 0x00000000 /* X */,
		        0x0000000A /* Y */, 0x00000000 /* Z */
};


/* -------------------------------------------- */
/* internal replacement of gets                 */
/* -------------------------------------------- */
static char *sGets(char *buffer, int size) {

        char *ebuf;

        if (! fgets(buffer, size-1, stdin))
           return NULL;

        /* remove trailing line feed */

        ebuf = buffer + strlen(buffer);

        while (--ebuf >= buffer) {
           if ((*ebuf == '\n') || (*ebuf == '\r'))
                *ebuf = '\000';
           else
                break;
        }

        return buffer;
}

/* -------------------------------------------- */
/* Interface                                    */
/* -------------------------------------------- */
