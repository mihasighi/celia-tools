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


#include <stdio.h>
#include "mset.h"
#include "mset_internal.h"
#include "shape_macros.h"
#include "apron2shape.h"
#include "ap_abstract0.h"
#include "ap_generic.h"
#include "pk.h"
#include "pkeq.h"
#include "ap_ppl.h"



/* ============================================================ */
/* Internal Management */
/* ============================================================ */

/* generic allocation routine */
inline mset_t *
mset_alloc_internal(mset_internal_t * pr, size_t intdim, size_t realdim) {
    mset_t *r;
    checked_malloc(r, mset_t, 1, sizeof (mset_t), return NULL;);
    r->datadim = intdim;
    r->segmdim = realdim;
    r->mscons = NULL;
    r->dcons = NULL;
    return r;
}

/* returns a top mset */
inline mset_t *
mset_alloc_top(mset_internal_t * pr, size_t intdim, size_t realdim) {
    mset_t *r = mset_alloc_internal(pr, intdim, realdim);
    r->dcons = ap_abstract0_top(pr->man_dcons, intdim + 2 * r->segmdim, 0);
    r->mscons =
            ap_abstract0_top(pr->man_mscons, DATA_DIM(intdim, 0) + 2 * r->segmdim,
            0);
    return r;
}

inline void
mset_free_internal(mset_internal_t * pr, mset_t * a) {
    if (a) {
        if (a->dcons) {
            ap_abstract0_free(pr->man_dcons, a->dcons);
            a->dcons = NULL;
        }
        if (a->mscons) {
            ap_abstract0_free(pr->man_mscons, a->mscons);
            a->mscons = NULL;
        }
        free(a);
    }
    return;
}

inline mset_t *
mset_copy_internal(mset_internal_t * pr, mset_t * a) {
    arg_assert(a, return NULL;);
    mset_t *r = mset_alloc_internal(pr, a->datadim, a->segmdim);
    r->dcons = ap_abstract0_copy(pr->man_dcons, a->dcons);
    r->mscons = ap_abstract0_copy(pr->man_mscons, a->mscons);
    return r;
}


/* ============================================================ */
/* Memory */

/* ============================================================ */

mset_t *
mset_copy(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_COPY, 0);
    return (a) ? mset_copy_internal(pr, a) : NULL;
}

void
mset_free(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_FREE, 0);
    mset_free_internal(pr, a);
}

size_t
mset_size(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_ASIZE, 0);
    if (!a)
        return 0;
    size_t s = ap_abstract0_size(pr->man_dcons, a->dcons);
    s += ap_abstract0_size(pr->man_mscons, a->mscons);
    return s;
}


/* ============================================================ */
/* Control of internal representation */
/* ============================================================ */

/* Return the set of minimal elements */
void
mset_minimize(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_MINIMIZE, 0);
    ap_abstract0_minimize(pr->man_dcons, a->dcons);
    ap_abstract0_minimize(pr->man_mscons, a->mscons);
    // TODO: put also saturation
}

/* Return the set of canonical elements */
void
mset_canonicalize(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_CANONICALIZE, 0);
    ap_abstract0_canonicalize(pr->man_dcons, a->dcons);
    ap_abstract0_canonicalize(pr->man_mscons, a->mscons);
    // TODO: put also saturation
}

/* TODO: priority 0 */
int
mset_hash(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_HASH, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return 0;
}

/* NOT IMPLEMENTED: do nothing */
void
mset_approximate(ap_manager_t * man, mset_t * a, int algorithm) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_APPROXIMATE, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
}

/* TODO: used ? */
bool
mset_is_minimal(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_CANONICALIZE, 0);
    return true;
}

/* TODO: used ? */
bool
mset_is_canonical(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_CANONICALIZE, 0);
    return true;
}

/* ============================================================ */
/* Basic Constructors */

/* ============================================================ */

mset_t *
mset_bottom(ap_manager_t * man, size_t intdim, size_t realdim) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_BOTTOM, 0);
    mset_t *r = mset_alloc_internal(pr, intdim, realdim);
    // all constraints are NULL, i.e., bottom
    return r;
}

mset_t *
mset_top(ap_manager_t * man, size_t intdim, size_t realdim) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_TOP, 0);
    mset_t *r = mset_alloc_top(pr, intdim, realdim);
    return r;
}

/* put constraints on data variables */
mset_t *
mset_of_box(ap_manager_t * man, size_t intdim, size_t realdim,
        ap_interval_t ** t) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_OF_BOX, 0);
    mset_t *r = mset_alloc_top(pr, intdim, realdim);
    ap_abstract0_free(pr->man_dcons, r->dcons);
    r->dcons = ap_abstract0_of_box(pr->man_dcons, intdim, 0, t); /* TODO */
    // TODO: put also saturation
    return r;
}

/* NOT IMPLEMENTED */
mset_t *
mset_of_generator_array(ap_manager_t * man, size_t intdim, size_t realdim,
        ap_generator0_array_t * ar) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_ADD_RAY_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}


/* ============================================================ */
/* Accessors */

/* ============================================================ */

ap_dimension_t
mset_dimension(ap_manager_t * man, mset_t * a) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_DIMENSION, 0);
    ap_dimension_t r;
    r.intdim = 0;
    r.realdim = 0;
    if (a) {
        r.intdim = a->datadim;
        r.realdim = a->segmdim;
    }
    return r;
}

/* ============================================================ */
/* Managers */

/* ============================================================ */

void
mset_internal_free(mset_internal_t * pr) {
    // ap_abstract0_manager_free(pr->man_dcons);
    pr->man_dcons = NULL;
    // ap_pkgrid_manager_free(pr->man_mscons);
    pr->man_mscons = NULL;
    free(pr);
}

ap_manager_t *
mset_manager_alloc(void) {
    ap_manager_t *man;
    mset_internal_t *pr;
    size_t i;

    pr = (mset_internal_t *) malloc(sizeof (mset_internal_t));
    assert(pr);

    // pr->man_dcons = pk_manager_alloc(false);
    pr->man_dcons = ap_ppl_poly_manager_alloc(false);
    pr->man_mscons = pkeq_manager_alloc();
    man = ap_manager_alloc("mset", "0.1 with (dcons=pk,mscons=pkeq)", pr,
            (void (*)(void *)) mset_internal_free);

    pr->man = man;

    man->funptr[AP_FUNID_COPY] = &mset_copy;
    man->funptr[AP_FUNID_FREE] = &mset_free;
    man->funptr[AP_FUNID_ASIZE] = &mset_size;
    man->funptr[AP_FUNID_MINIMIZE] = &mset_minimize;
    man->funptr[AP_FUNID_CANONICALIZE] = &mset_canonicalize;
    man->funptr[AP_FUNID_HASH] = &mset_hash;
    man->funptr[AP_FUNID_APPROXIMATE] = &mset_approximate;
    man->funptr[AP_FUNID_FPRINT] = &mset_fprint;
    man->funptr[AP_FUNID_FPRINTDIFF] = &mset_fprintdiff;
    man->funptr[AP_FUNID_FDUMP] = &mset_fdump;
    man->funptr[AP_FUNID_SERIALIZE_RAW] = &mset_serialize_raw;
    man->funptr[AP_FUNID_DESERIALIZE_RAW] = &mset_deserialize_raw;
    man->funptr[AP_FUNID_BOTTOM] = &mset_bottom;
    man->funptr[AP_FUNID_TOP] = &mset_top;
    man->funptr[AP_FUNID_OF_BOX] = &mset_of_box;
    man->funptr[AP_FUNID_DIMENSION] = &mset_dimension;
    man->funptr[AP_FUNID_IS_BOTTOM] = &mset_is_bottom;
    man->funptr[AP_FUNID_IS_TOP] = &mset_is_top;
    man->funptr[AP_FUNID_IS_LEQ] = &mset_is_leq;
    man->funptr[AP_FUNID_IS_EQ] = &mset_is_eq;
    man->funptr[AP_FUNID_IS_DIMENSION_UNCONSTRAINED] =
            &mset_is_dimension_unconstrained;
    man->funptr[AP_FUNID_SAT_INTERVAL] = &mset_sat_interval;
    man->funptr[AP_FUNID_SAT_LINCONS] = &mset_sat_lincons;
    man->funptr[AP_FUNID_SAT_TCONS] = &mset_sat_tcons;
    man->funptr[AP_FUNID_BOUND_DIMENSION] = &mset_bound_dimension;
    man->funptr[AP_FUNID_BOUND_LINEXPR] = &mset_bound_linexpr;
    man->funptr[AP_FUNID_BOUND_TEXPR] = &mset_bound_texpr;
    man->funptr[AP_FUNID_TO_BOX] = &mset_to_box;
    man->funptr[AP_FUNID_TO_LINCONS_ARRAY] = &mset_to_lincons_array;
    man->funptr[AP_FUNID_TO_TCONS_ARRAY] = &mset_to_tcons_array;
    man->funptr[AP_FUNID_TO_GENERATOR_ARRAY] = &mset_to_generator_array;
    man->funptr[AP_FUNID_MEET] = &mset_meet;
    man->funptr[AP_FUNID_MEET_ARRAY] = &mset_meet_array;
    man->funptr[AP_FUNID_MEET_LINCONS_ARRAY] = &mset_meet_lincons_array;
    man->funptr[AP_FUNID_MEET_TCONS_ARRAY] = &mset_meet_tcons_array;
    man->funptr[AP_FUNID_JOIN] = &mset_join;
    man->funptr[AP_FUNID_JOIN_ARRAY] = &mset_join_array;
    man->funptr[AP_FUNID_ADD_RAY_ARRAY] = &mset_add_ray_array;
    man->funptr[AP_FUNID_ASSIGN_LINEXPR_ARRAY] = &mset_assign_linexpr_array;
    man->funptr[AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY] =
            &mset_substitute_linexpr_array;
    man->funptr[AP_FUNID_ASSIGN_TEXPR_ARRAY] = &mset_assign_texpr_array;
    man->funptr[AP_FUNID_SUBSTITUTE_TEXPR_ARRAY] = &mset_substitute_texpr_array;
    man->funptr[AP_FUNID_ADD_DIMENSIONS] = &mset_add_dimensions;
    man->funptr[AP_FUNID_REMOVE_DIMENSIONS] = &mset_remove_dimensions;
    man->funptr[AP_FUNID_PERMUTE_DIMENSIONS] = &mset_permute_dimensions;
    man->funptr[AP_FUNID_FORGET_ARRAY] = &mset_forget_array;
    man->funptr[AP_FUNID_EXPAND] = &mset_expand;
    man->funptr[AP_FUNID_FOLD] = &mset_fold;
    man->funptr[AP_FUNID_WIDENING] = &mset_widening;
    man->funptr[AP_FUNID_CLOSURE] = &mset_closure;

    for (i = 0; i < AP_EXC_SIZE; i++)
        ap_manager_set_abort_if_exception(man, i, false);

    return man;
}

mset_t *
mset_of_abstract0(ap_abstract0_t * a) {
    return (mset_t *) a->value;
}

ap_abstract0_t *
abstract0_of_mset(ap_manager_t * man, mset_t * a) {
    ap_abstract0_t *r = malloc(sizeof (ap_abstract0_t));
    assert(r);
    r->value = a;
    r->man = ap_manager_copy(man);
    return r;
}
