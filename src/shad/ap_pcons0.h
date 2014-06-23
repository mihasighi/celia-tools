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


#ifndef AP_PCONS0_H_
#define AP_PCONS0_H_


#include "ap_dimension.h"
#include "ap_expr0.h"
#include "uthash.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C" {
#endif
    /* *INDENT-ON* */

    /* Offsets are applied to pointer variables and may be:
     * - field offsets given by (ofs-intdim)
     * - function offsets (used in logic)
     * - index offsets are in 0..intdim-1
     */
    typedef enum {
        OFFSET_DATA = 0, /* for data field */
        OFFSET_NEXT = 1, /* for next field */
        OFFSET_PREV = 2, /* for next field */
        OFFSET_LEN = -3, /* for length function */
        OFFSET_SUM = -4, /* for sum function */
        OFFSET_MSET = -5, /* for mset function */
        OFFSET_UCONS = -6, /* for ucons function */
        OFFSET_SL3 = -7, /* for sl3 formulas */
        OFFSET_NONE = -8, /* special use, take care! */
    } offset_t;

    /* Maximum number of fields supported, @see ptrfields_max in interp_eqn.ml */
#define OFFSET_FLD_MAX 10

    
    /**
     * Pointer constraints.
     *
     * Constraints allowed are:
     * - pointer constraints
     *     acyclic(x)    encoded by x > 0.0
     *     cyclic(x)     encoded by x->next > x
     *     x == NULL     encoded by x(->next) == 0.0
     *     x != NULL     encoded by x(->next) <> 0.0
     *     disjoint(x,y) encoded by x(->next) mod_1 y(->next)
     *     same(x,y)     encoded by x(->next) mod_0 y(->next)
     *     reach(x,y)    encoded by x(->next) >  y(->next)
     *
     * - integer constraint
     *     encoded by a linear constraint with dimension intdim+ptrdim,
     *     i.e., each ptr var (ptrdim) has an attached data and
     *           a vector of offsets representing the function
     *           to be applied to ptr dimensions
     */

    typedef enum {
        DATA_CONS = 0, /* only constraint on int vars or data fields */
        ACYCLIC_CONS, /* unary segment description */
        CYCLIC_CONS,
        ISO_CONS, /* segment isomorphism */
        SEP_CONS, /* segment separation */
        EQ_CONS, /* ptr equality */
        NE_CONS, /* ptr difference */
        REACH_CONS, /* binary reachability */
	SL3_CONS,  /* SL3 constraint */
        OTHER_CONS,
    } pcons0_typ_t;

    typedef struct {
        pcons0_typ_t type;

        size_t intdim;
        size_t ptrdim;

        union {

            struct {
                ap_dim_t x, y; /* in [intdim..intdim+ptrdim-1] or NULL_DIM */
                int offx, offy; /* offset of each variable */
            } ptr;

            /* TODO: support all the offsets in a constraint */
            struct {
                ap_lincons0_t cons; /* over dimensions intdim + ptrdim */
                int* offsets; /* of size ptrdim, only for ptr dimensions above */
            } data;

        } info;

        ap_lincons0_t *lcons;
        ap_tcons0_t *tcons;
        UT_hash_handle hh; /* make structure hashable */
        /* keys are lcons--tcons */
    } pcons0_t;

    /* Array of constraints */
    typedef struct pcons0_array_t {
        pcons0_t **p;
        size_t size;
    } pcons0_array_t;

    void shape_pcons0_clear(pcons0_t * c);
    /* Clear the data constraint in c */
    void shape_pcons0_array_clear(pcons0_array_t * array);
    /* Clear the constraints of the array, and then the array itself */


    void shape_offset_fprint(FILE * stream, int ofs, size_t intdim, size_t dim);

    void shape_offsets_fprint(FILE * stream, int* offsets, size_t intdim, size_t dim);

    void shape_pcons_fdump(FILE * stream, pcons0_t * c);
    void shape_pcons_array_fdump(FILE * stream, pcons0_array_t * a);
    /* Printing */

    /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* AP_PCONS0_H_ */
