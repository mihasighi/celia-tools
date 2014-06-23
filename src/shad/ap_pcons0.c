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



#include "ap_pcons0.h"

/* To use only when htable of pcons is freed */
void
shape_pcons0_clear(pcons0_t * c) {
    if (c->type == DATA_CONS) {
        ap_lincons0_clear(&c->info.data.cons);
        if (c->info.data.offsets) free(c->info.data.offsets);
    }
}

void
shape_pcons0_array_clear(pcons0_array_t * array) {
    size_t i;

    if (array->p != NULL) {
        for (i = 0; i < array->size; i++)
            if (array->p[i] != NULL) {
                //Warning:all pcons0_t values are hashed, so not to be freed !
                //shape_pcons0_clear(array->p[i]);
                //free(array->p[i]);
                array->p[i] = NULL;
            }
        free(array->p);
        array->p = NULL;
    }
}

void
shape_dim_fprint(FILE * stream, size_t dim) {
    if (dim == AP_DIM_MAX) // NULL_DIM
        fprintf(stream, "NULL");
    else
        fprintf(stream, "x%zu", dim);
}

void
shape_offset_fprint(FILE * stream, int ofs, size_t intdim, size_t dim) {
    if (dim < intdim) {
        /* integer variable */
        shape_dim_fprint(stream, dim);
        return;
    }

    /* pointer variable */

    if (ofs >= 0 && ofs < ((int) intdim)) {
        /* ptr variable with index offset */
        fprintf(stream, "x%d[x%d]", dim, ofs);
    } else if (ofs >= ((int) intdim)) {
        /* ptr variable with field access */
        switch (ofs - intdim) {
            case OFFSET_DATA: fprintf(stream, " \\data(");
                break;
            case OFFSET_NEXT: fprintf(stream, " \\next(");
                break;
            case OFFSET_PREV: fprintf(stream, " \\prev(");
                break;
            default: fprintf(stream, " \\fld%d(", ofs - intdim);
                break;
        }
        shape_dim_fprint(stream, dim);
        fprintf(stream, ") ");
    } else if (ofs == OFFSET_NONE)
        shape_dim_fprint(stream, dim);
    else {
        /* ptr variable with function application, for logic */
        switch (ofs) {
            case OFFSET_LEN: fprintf(stream, " \\length(");
                break;
            case OFFSET_SUM: fprintf(stream, " \\sum(");
                break;
            case OFFSET_MSET: fprintf(stream, " \\mset(");
                break;
            case OFFSET_UCONS: fprintf(stream, " \\ucons(");
                break;
            default:
                fprintf(stream, " \\unknown(");
        }
        shape_dim_fprint(stream, dim);
        fprintf(stream, ") ");
    }
}

void
shape_offsets_fprint(FILE * stream, int* offsets, size_t intdim, size_t ptrdim) {
#ifndef NDEBUG
    fprintf(stdout, "\n====shape_offsets_fprint: dim=(%d)\n", ptrdim);
#endif
    if (offsets == NULL)
        fprintf(stream, " [NULL offsets]\n");
    else {
        size_t i;
        fprintf(stream, " offsets = [");
        for (i = 0; i < ptrdim; i++)
            shape_offset_fprint(stream, offsets[i], intdim, i + intdim);
        fprintf(stream, "]\n");
    }
}

void
shape_pcons_fdump(FILE * stream, pcons0_t * c) {
    if (!c)
        fprintf(stream, "[NULL ptr constraint]\n");
    else {
#ifndef NDEBUG
        fprintf(stdout, "\n====shape_pcons_fdump: dim=(%zu,%zu)\n", c->intdim, c->ptrdim);
        fflush(stdout);
#endif
        if (c->type == DATA_CONS) {
            ap_lincons0_fprint(stream, &c->info.data.cons, NULL);
            shape_offsets_fprint(stream, c->info.data.offsets, c->intdim, c->ptrdim);
        } 
	if (c->type == SL3_CONS) {
	  fprintf(stream, "pan/spec_");
	  ap_lincons0_fprint(stream, &c->info.data.cons, NULL);
	  fprintf(stream, ".smt");
	} else {
            shape_offset_fprint(stream, c->info.ptr.offx, c->intdim, c->info.ptr.x);
            switch (c->type) {
                case EQ_CONS:
                    fprintf(stream, " == ");
                    break;
                case NE_CONS:
                    fprintf(stream, " != ");
                    break;
                case REACH_CONS:
                    fprintf(stream, " --> ");
                    break;
                case ISO_CONS:
                    fprintf(stream, " iso ");
                    break;
                case SEP_CONS:
                    fprintf(stream, " * ");
                    break;
                default:
                    fprintf(stream, " Internal Error ");
                    break;
            }
            shape_offset_fprint(stream, c->info.ptr.offx, c->intdim, c->info.ptr.x);
        }
    }
    fflush(stream);
}

void
shape_pcons_array_fdump(FILE * stream, pcons0_array_t * array) {
    if (!array || array->size == 0 || array->p == NULL)
        fprintf(stream, "[empty array]");
    else {
        size_t i;
        fprintf(stream, "[ ");
        for (i = 0; i < array->size; i++) {
            shape_pcons_fdump(stream, array->p[i]);
            if (i % 4 == 3)
                fprintf(stream, ",\n\t");
            else
                fprintf(stream, ", ");
        }
        fprintf(stream, "]");
    }
    fflush(stream);

}
