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


#include "shape.h"
#include "shape_internal.h"
#include "shape_macros.h"
#include "apron2shape.h"
#include "ap_generic.h"



/* ============================================================ */
/* Tests */
/* ============================================================ */

/* The bottom value is the empty set. */
bool
shape_is_bottom(ap_manager_t * man, shape_t * a) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_IS_BOTTOM, 0);
    return (!a || a->m.size == 0 || a->msize == 0);
}

/* The top value is top ushape in the first position. */
bool
shape_is_top(ap_manager_t * man, shape_t * a) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_IS_TOP, 0);
    return (a && a->msize == 1 && ushape_is_top(man, a->m.p[0]));
}

bool
shape_is_leq(ap_manager_t * man, shape_t * a1, shape_t * a2) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_IS_LEQ, 0);
    size_t i, j;
    bool r;
    if (a1 == NULL || (a2 == NULL && shape_is_bottom(man, a1)))
        return true;
    if (a2 == NULL)
        return false;
#ifndef NDEBUG
    fprintf(stdout, "\n\tshape_is_leq: args ");
    shape_set_print_smt();
    shape_fdump_le(stdout, pr, a1, a2);
    shape_set_print_dot();
    fflush(stdout);
#endif
    /* for each i in a1 shall exists a j in a2 s.t. i <= j */
    r = true;
    for (i = 0; i < a1->m.size && r; i++)
        if (a1->m.p[i]) {
            for (j = 0; j < a2->m.size; j++)
                if (a2->m.p[j] && ushape_is_leq(man, a1->m.p[i], a2->m.p[j]))
                    break;
            if (j == a2->m.size) /* not element j s.t. i <= j */
                r = false;
        }
#ifndef NDEBUG
    fprintf(stdout, "\n\tshape_is_leq: returns %d\n", (int) r);
#endif
    return r;
}

/* TODO: better complexity?? */
bool
shape_is_eq(ap_manager_t * man, shape_t * a1, shape_t * a2) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_IS_EQ, 0);
    size_t i, j;
    bool *r;
    if (a1 == NULL && a2 == NULL)
        return true;
    if ((!a1 && a2) || (a1 && !a2) || (a1->msize != a2->msize))
        return false;
    checked_malloc(r, bool, sizeof (bool), a2->m.size, return false;
            );
    for (i = 0; i < a1->m.size; i++)
        if (a1->m.p[i]) {
            for (j = 0; j < a2->m.size; j++)
                if (a2->m.p[j] && ushape_is_eq(man, a1->m.p[i], a2->m.p[j]))
                    break;
            if (j == a2->m.size) {
                free(r);
                return false;
            }
            r[j] = true;
        }
    for (j = 0; j < a2->m.size; j++)
        if (a2->m.p[j] && !r[j]) {
            for (i = 0; i < a1->m.size; i++)
                if (a1->m.p[i] && ushape_is_eq(man, a1->m.p[i], a2->m.p[j]))
                    break;
            if (i == a1->m.size) {
                free(r);
                return false;
            }
        }
    free(r);
    return true;
}

/*
 * Possibly needed for checking linear constraints representing - aliasing
 * between pointer variables, e.g., x = y, - or constraints between the
 * program variables (lengths and data) in the assert statements.
 */
bool
shape_sat_lincons(ap_manager_t * man, shape_t * a, ap_lincons0_t * lincons) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_SAT_LINCONS, 0);
    size_t i;
    if (!a || !a->m.p)
        return false;
    for (i = 0; i < a->m.size; i++)
        if (a->m.p[i] && ushape_sat_lincons(man, a->m.p[i], lincons))
            return true;
    return false;
}

/*
 * Needed for checking constraints representing aliasing between pointer
 * variables, e.g., x*next = y, in the assert statements.
 */
bool
shape_sat_tcons(ap_manager_t * man, shape_t * a, ap_tcons0_t * cons) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_SAT_TCONS, 0);
    size_t i;
    if (!a || !a->m.p)
        return false;
    for (i = 0; i < a->m.size; i++)
        if (a->m.p[i] && ushape_sat_tcons(man, a->m.p[i], cons))
            return true;
    return false;
}

/* Interval constraints are only between non-pointer variables */
bool
shape_sat_interval(ap_manager_t * man, shape_t * a,
        ap_dim_t dim, ap_interval_t * itv) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_SAT_INTERVAL, 0);
    size_t i;
    if (!a || !a->m.p)
        return false;
    for (i = 0; i < a->m.size; i++)
        if (a->m.p[i] && ushape_sat_interval(man, a->m.p[i], dim, itv))
            return true;
    return false;
}

bool
shape_is_dimension_unconstrained(ap_manager_t * man, shape_t * a,
        ap_dim_t dim) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED, 0);
    size_t i;
    bool r = true;
    if (!a || !a->msize || !a->m.p)
        return false;
    for (i = 0; i < a->m.size && r; i++)
        r = ushape_is_dimension_unconstrained(man, a->m.p[i], dim);
    return r;
}

/* ============================================================ */
/* Extraction of properties */
/* ============================================================ */

/* NOT IMPLEMENTED */
ap_interval_t *
shape_bound_linexpr(ap_manager_t * man, shape_t * a, ap_linexpr0_t * expr) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_BOUND_LINEXPR, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
shape_bound_texpr(ap_manager_t * man, shape_t * a, ap_texpr0_t * expr) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_BOUND_TEXPR, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
shape_bound_dimension(ap_manager_t * man, shape_t * a, ap_dim_t dim) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_BOUND_DIMENSION, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/*
 * Used to print finally the abstract values at control points.
 * Generate lincons null = <file number>.
 */
ap_lincons0_array_t
shape_to_lincons_array(ap_manager_t * man, shape_t * a) {
    ap_lincons0_array_t ar;
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_TO_LINCONS_ARRAY, 0);
    ar = ap_lincons0_array_make(1);
    ap_linexpr0_t *expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, 0);
    ap_coeff_set_scalar_int(&expr->cst, pr->filenum);
    shape_fdump(stdout, man, a);
    ar.p[0] = ap_lincons0_make(AP_CONS_SUPEQ, expr, NULL);
    return ar;
}

/* NOT IMPLEMENTED */
ap_tcons0_array_t
shape_to_tcons_array(ap_manager_t * man, shape_t * a) {
    return ap_generic_to_tcons_array(man, a);
}

/* NOT IMPLEMENTED */
ap_interval_t **
shape_to_box(ap_manager_t * man, shape_t * a) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_TO_BOX, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/* NOT IMPLEMENTED */
ap_generator0_array_t
shape_to_generator_array(ap_manager_t * man, shape_t * a) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_TO_GENERATOR_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return ap_generator0_array_make(0);
}
