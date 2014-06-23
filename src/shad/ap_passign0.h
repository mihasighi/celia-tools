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


#ifndef AP_PASSIGN0_H_
#define AP_PASSIGN0_H_

#include "ap_dimension.h"
#include "ap_expr0.h"
#include "uthash.h"
#include "ap_pcons0.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C" {
#endif
    /* *INDENT-ON* */

    /**
     * Internal representation of assignments. Assignments allowed are:
     *
     * - pointer           x[->f] = y[->f]   [f is ptr field]
     * [Note: before assigning a pointer, we suppose that it is first put to null]
     *
     * - integer content   x->data = expr
     * [expr is built from integer variables and integer content y->data]
     *
     * - integer variables v = expr
     * [expr is built from integer variables and integer content y->data]
     *
     * - allocation        x = alloc([n])   [n constant or integer var]
     *
     * - memory free       free(x)
     */
    typedef enum /* kind of assignment */ {
        PA_ASSIGN_INT = 0,
        PA_ASSIGN_PTR,
        PA_ALLOC,
        PA_ALLOC_N,
        PA_FREE,
        PA_OTHER
    } passign0_op_t;

    typedef struct {
        passign0_op_t op; /* type */

        /* left hand side */
        ap_dim_t x; /* lhs in [0..intdim+ptrdim-1] */
        int offx; /* offset of x (field or index) */
        /* fields are in negatives,
         * 0 is data field or index,
         * positive are index (integer vars) dimensions
         */

        /* right hand side */
        union {

            struct {
                ap_dim_t y; /* ptr dimension in [intdim..intdim+ptrdim-1] or NULL_DIM */
                int offy;   /* offset of y (index, field, function) */
            } ptr;

            struct {
                ap_linexpr0_t *expr; /* integer expr over dimensions [0..intdim+ptrdim-1] */
                int* offsets; /* only for ptr dimensions (data=0 or index dimension) */
            } data;

            struct {
                bool cst; /* alloc with a constant */

                union {
                    size_t c;
                    ap_dim_t l;
                } size;
            } alloc;

        } info;

        ap_dim_t lhs;
        size_t intdim;
        size_t ptrdim;
        ap_linexpr0_t *lexpr;
        ap_texpr0_t *texpr;
        UT_hash_handle hh; /* make structure hashable */
        /* keys are lhs--texpr */
    } passign0_t;

        /* Array of assignments */
    typedef struct passign0_array_t {
        passign0_t **p;
        size_t size;
    } passign0_array_t;

    void shape_passign0_clear(passign0_t * c);
    /* Clear the data constraint in c */
    void shape_passign0_array_clear(passign0_array_t * array);
    /* Clear the constraints of the array, and then the array itself */

    void shape_passign_fdump(FILE * stream, passign0_t * c,
            size_t intdim, size_t ptrdim);
    void shape_passign_array_fdump(FILE * stream, passign0_array_t * a,
            size_t intdim, size_t ptrdim);
    /* Printing */

    /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* AP_PASSIGN0_H_ */
