/**************************************************************************/
/*                                                                        */
/*  CINV Library / Shape Domain                                           */
/*                                                                        */
/*  Copyright (C) 2009-2011                                               */
/*    LIAFA (University of Paris Diderot and CNRS)                        */
/*                                                                        */
/*                                                                        */
/*  you can redistribute it and/or modify it under the terms of the GNU   */
/*  Lesser General Public License as published by the Free Software       */
/*  Foundation, version 3.                                                */
/*                                                                        */
/*  It is distributed in the hope that it will be useful,                 */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*  GNU Lesser General Public License for more details.                   */
/*                                                                        */
/*  See the GNU Lesser General Public License version 3.                  */
/*  for more details (enclosed in the file LICENSE).                      */
/*                                                                        */
/**************************************************************************/


#include "ap_passign0.h"
#include "ap_linexpr0.h"

  /* ================================================================== */
  /* Globals */
  /* ================================================================== */

passign0_t *passigns_ht; /* global hash table of ptr constraints */

void
ap_passign0_init (void)
{
  passigns_ht = NULL;
}

  /* ================================================================== */
  /* Constructors/Destructors */
  /* ================================================================== */


/* To use only when htable of passigns is freed */
void
shape_passign0_clear(passign0_t * c) {
    if (c->op == PA_ASSIGN_INT) {
        ap_linexpr0_free(c->info.data.expr);
        if (c->info.data.offsets) free(c->info.data.offsets);
    }
}

void
shape_passign0_array_clear(passign0_array_t * array) {
    size_t i;

    if (array->p != NULL) {
        for (i = 0; i < array->size; i++)
            if (array->p[i] != NULL) {
                //Warning:all passign0_t values are hashed, so not to be freed !
                //shape_passign0_clear(array->p[i]);
                //free(array->p[i]);
                array->p[i] = NULL;
            }
        free(array->p);
        array->p = NULL;
    }
}

/* ===================================================================== */
/* Global Set Manipulation */
/* ===================================================================== */

/* Search a value in the htable using keys lhs, lexpr, and texpr */
passign0_t *
shape_passign_search (ap_dim_t lhs,
                      size_t intdim, size_t ptrdim,
                      ap_linexpr0_t * lexpr, ap_texpr0_t * texpr)
{
  passign0_t *r, rr;
  unsigned keylen;
  /* search in the htable */
  keylen = offsetof (passign0_t, hh) - offsetof (passign0_t, lhs);
  r = NULL;
  rr.lhs = lhs;
  rr.intdim = intdim;
  rr.ptrdim = ptrdim;
  rr.lexpr = lexpr;
  rr.texpr = texpr;
  HASH_FIND (hh, passigns_ht, &rr.lhs, keylen, r);

  return r;
}

/* Add a value to the htable */
passign0_t *
shape_passign_add (passign0_t * a)
{
  passign0_t *r;
  unsigned keylen;
  /* search in the htable */
  keylen = offsetof (passign0_t, hh) - offsetof (passign0_t, lhs);
  r = NULL;
  HASH_FIND (hh, passigns_ht, &a->lhs, keylen, r);
  if (!r)
    {
      HASH_ADD (hh, passigns_ht, lhs, keylen, a);
      HASH_FIND (hh, passigns_ht, &a->lhs, keylen, r);
    }
  return r;
}


  /* ================================================================== */
  /* Printing */
  /* ================================================================== */

void
shape_passign_fdump(FILE * stream, passign0_t * a,
        size_t intdim, size_t ptrdim) {
    if (!a)
        fprintf(stream, "[NULL ptr assignement]");
    else if (a->intdim != intdim || a->ptrdim != ptrdim)
        fprintf(stream, "[Bad dimensions (%zu,%zu) in assignement of dimensions (%zu,%zu)]",
                intdim, ptrdim, a->intdim, a->ptrdim);
    else {
        shape_offset_fprint(stream, a->offx, intdim, a->x);

        fprintf(stream, " :=(kind=%d) ", (int) a->op);

        if (a->op == PA_ASSIGN_INT) {
            ap_linexpr0_fprint(stream, a->info.data.expr, NULL);
            shape_offsets_fprint(stream, a->info.data.offsets, intdim, ptrdim);
        }
        else if (a->op == PA_ASSIGN_PTR)
            shape_offset_fprint(stream, a->info.ptr.offy, intdim, a->info.ptr.y);
        else if (a->op == PA_ALLOC || a->op == PA_ALLOC_N)
            fprintf(stream, "new()"); /* TODO: add size of allocation */

        else if (a->op == PA_FREE)
            fprintf(stream, "free(%zu)", (size_t) a->x);

        else
            fprintf(stream, "unknown assignment");

    }
    fflush(stream);
}

void
shape_passign_array_fdump(FILE * stream, passign0_array_t * array,
        size_t intdim, size_t ptrdim) {
    if (!array || array->size == 0 || array->p == NULL)
        fprintf(stream, "[empty array]");
    else {
        size_t i;
        fprintf(stream, "[");
        for (i = 0; i < array->size; i++) {
            shape_passign_fdump(stream, array->p[i], intdim, ptrdim);
            if (i % 4 == 3)
                fprintf(stream, ",\n\t");
            else
                fprintf(stream, ", ");

        }
        fprintf(stream, "]");
    }
    fflush(stream);

}
