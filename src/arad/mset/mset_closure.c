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


#include "mset.h"
#include "mset_internal.h"
#include "apron2shape.h"
#include "ap_generic.h"
#include "ap_linexpr0.h"
#include "ap_lincons0.h"


/* All closures are in-place. */


/* ============================================================ */
/* Full Closure */
/* ============================================================ */

/*
 * TODO: priority 1 What it is here closure? - put in the set form - close
 * all the units
 */

bool
mset_close(mset_t * g, size_t dim) {
    return false;
}


/* ============================================================ */
/* Incremental Closure */
/* ============================================================ */

/*
 * TODO: priority 1
 *
 * ?? time. ?? space.
 */

bool
mset_close_incremental(mset_t * m, size_t dim, size_t v) {
    return false;
}


/* ============================================================ */
/* Sanity Check */
/* ============================================================ */

/* TODO: priority 1 */
bool
mset_check_closed(mset_t * m, size_t dim) {
    return false;
}


/* ============================================================ */
/* Topological closure operation */

/* ============================================================ */

mset_t *
mset_closure(ap_manager_t * man, bool destructive, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_CLOSURE, 0);
    if (!a)
        return NULL;
    mset_t *r = mset_alloc_internal(pr, a->datadim, a->segmdim);
    r->dcons = ap_abstract0_closure(pr->man_dcons, false, a->dcons);
    r->mscons = ap_abstract0_closure(pr->man_mscons, false, a->mscons);
    mset_strengthen_all(pr, r);
    if (destructive)
        mset_free_internal(pr, a);
    return r;
}


/* ============================================================ */
/* Strengthening operation */
/* ============================================================ */

/* Strengthen after changes done on dimension dim */
void
mset_strengthen_dim(mset_internal_t * pr, mset_t *a, ap_dim_t dim) {
    /* Heuristic 1: propagate all equalities on dim from a->dcons
     * to a->mscons
     */
    ap_linexpr0_t *lexpr = NULL;
    ap_lincons0_array_t arr = ap_lincons0_array_make(1);
    arr.p[0].constyp = AP_CONS_EQ;
    arr.p[0].scalar = NULL;
    size_t di;
    size_t alldim = a->datadim + a->segmdim;
    for (di = 0; di < alldim; di++)
        if (dim != di &&
                di != a->datadim // node NULL
                ) {
            // check on a->dcons the constraint data(dim) = di
            lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE,
                    a->datadim + 2 * a->segmdim);
            arr.p[0].linexpr0 = lexpr;
            ap_linexpr0_set_coeff_scalar_int(lexpr, dim, -1);
            ap_linexpr0_set_coeff_scalar_int(lexpr, di, +1);
            if (ap_abstract0_sat_lincons(pr->man_dcons, a->dcons, arr.p)) {
                // meet with the constraint
                ap_linexpr0_realloc(lexpr, DATA_DIM(a->datadim, 0) + 2 * a->segmdim);
                a->mscons =
                        ap_abstract0_meet_lincons_array(pr->man_mscons, true, a->mscons, &arr);
            }
            ap_linexpr0_free(lexpr);
        }
}

/* Strengthen for all dimensions */
void
mset_strengthen_all(mset_internal_t * pr, mset_t *a) {
    /* TODO */
}
