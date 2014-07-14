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


#include "ushape.h"
#include "ushape_internal.h"
#include "hgraph_internal.h"
#include "apron2shape.h"
#include "shape_macros.h"
#include "ap_abstract0.h"
#include "ap_generic.h"




/* ============================================================ */
/* Tests */
/* ============================================================ */

/* ********************************************************************** */
/* VI. ushape_t */
/* ********************************************************************** */

/* The bottom value is the NULL hgraph and bottom constraint. */
bool
ushape_is_bottom(ap_manager_t * man, ushape_t * a) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_IS_BOTTOM, 0);
    if (!a)
        return true;
    size_t i;
    for (i = 0; i < pr->size_scons; i++)
        if (!a->scons || !a->scons[i]
                || ap_abstract0_is_bottom(pr->man_scons[i], a->scons[i]))
            return true;
    return (!a->h || hgraph_is_bottom(man, a->h));
}

/* The top value is top hgraph and top linear constraint. */
bool
ushape_is_top(ap_manager_t * man, ushape_t * a) {
    ushape_internal_t *pr = ushape_init_from_manager(man, AP_FUNID_IS_TOP, 0);
    if (!a)
        return false;
    bool istop = hgraph_is_top(man, a->h);
    if (pr->size_scons > 0 && a->scons) {
        size_t i;
        for (i = 0; i < pr->size_scons && istop; i++)
            istop = ap_abstract0_is_top(pr->man_scons[i], a->scons[i]);
    }
    return istop;
}

/*
 * Comparison is done only for ushapes with the same graph. The test is
 * correct but not complete: true means really true, false means don't know
 * or false
 */
bool
ushape_is_leq(ap_manager_t * man, ushape_t * a1, ushape_t * a2) {
    ushape_internal_t *pr = ushape_init_from_manager(man, AP_FUNID_IS_LEQ, 0);
    size_t i;
    bool is_leq = true;

#ifndef NDEBUG
    fprintf(stdout, "\n+++++ushape_is_leq: on a1=(");
    ushape_fdump(stdout, man, a1);
    fprintf(stdout, ") <= a2=(");
    ushape_fdump(stdout, man, a2);
    fprintf(stdout, "\n)\n");
    fflush(stdout);
#endif
    if (ushape_is_bottom(man, a1) || ushape_is_top(man, a2))
        is_leq = true;
    else if (!a2 || !hgraph_is_eq(man, a1->h, a2->h))
        is_leq = false;
    else {
        arg_assert(a1->datadim == a2->datadim
                && a1->ptrdim == a2->ptrdim, return false;
                );
        /* hgraphs are the same, look at node constraints */
        assert(a1->scons && a2->scons);
        for (i = 0; i < pr->size_scons && is_leq; i++)
            if (!(a1->scons[i] != NULL && a2->scons[i] != NULL &&
                    ap_abstract0_is_leq(pr->man_scons[i], a1->scons[i], a2->scons[i])))
                is_leq = false;
    }

#ifndef NDEBUG
    fprintf(stdout, "\n+++++ushape_is_leq: returns %d\n", (int) is_leq);
    fflush(stdout);
#endif
    return is_leq;
}

bool
ushape_is_eq(ap_manager_t * man, ushape_t * a1, ushape_t * a2) {
    ushape_internal_t *pr = ushape_init_from_manager(man, AP_FUNID_IS_EQ, 0);
    size_t i;
    if (ushape_is_bottom(man, a1))
        return true;
    if (!a2 || !hgraph_is_eq(man, a1->h, a2->h))
        return false;
    arg_assert(a1->datadim == a2->datadim
            && a1->ptrdim == a2->ptrdim, return false;
            );
    /* hgraphs are the same, look and constraints */
    for (i = 0; i < pr->size_scons; i++)
        if (!ap_abstract0_is_eq
                (pr->man_scons[i], (a1->scons) ? a1->scons[i] : NULL,
                (a2->scons) ? a2->scons[i] : NULL))
            return false;
    return true;
}

/*
 * Possibly needed for checking linear constraints representing - aliasing
 * between pointer variables, e.g., x = y, - or constraints between the
 * program variables (lengths and data) in the assert statements.
 */
bool
ushape_sat_lincons(ap_manager_t * man, ushape_t * a, ap_lincons0_t * lincons) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_SAT_LINCONS, 0);
    size_t i;
    if (ushape_is_bottom(man, a))
        return false;

    /* transform the linear constraint into a pointer constraint */
    pcons0_t *pcons =
            shape_pcons_of_lincons(pr, lincons, a->datadim, a->ptrdim);
    /* check the constraint on the hgraph and the data/length constraints */
    bool sat = false;
    switch (pcons->type) {
        case EQ_CONS:
        case NE_CONS:
        case REACH_CONS:
            sat = hgraph_sat_pcons(pr, a->h, pcons);
            break;
            /*
          case REACHL_CONS:
            {
              sat = hgraph_sat_pcons (pr, a->h, pcons);
              if (sat && pr->size_scons > 0)
                {
                   // the domains on node/segments shall satisfy l >= 1 and l[node of x]
                   // - l[node of y] = l
                  ap_lincons0_t l1 =
                    shape_lincons_x_cst (AP_CONS_SUPEQ, pcons->info.ptr.l, 1,
                                         a->datadim, a->ptrdim);
                  ap_lincons0_t lxy =
                    shape_lincons_x_y_l (AP_CONS_EQ, 1, pcons->info.ptr.x, -1,
                                         pcons->info.ptr.y, -1, 0,
                                         a->datadim, a->ptrdim);
                  size_t *v2n = hgraph_get_var2node (a->h);
                  ap_lincons0_t nlxy =
                    shape_lincons_of_node (pr, &lxy, NULL, v2n, a->h->size,
                                           a->datadim, a->ptrdim);
                  free (v2n);
                  sat = false;
                  for (i = 0; i < pr->size_scons && !sat; i++)
                    sat =
                      ap_abstract0_sat_lincons (pr->man_scons[i], a->scons[i], &l1)
                      && ap_abstract0_sat_lincons (pr->man_scons[i], a->scons[i],
                                                   &nlxy);
                   // TODO: if the universal domain is not dealing with lengths return
                   // false??
                  ap_lincons0_clear (&l1);
                  ap_lincons0_clear (&lxy);
                  ap_lincons0_clear (&nlxy);
                }
              break;
            }
             */
        case DATA_CONS:
        { /* linear constraint may be only on len/data constraints */
            ap_lincons0_t ncons = pcons->info.data.cons; // TODO: update with offsets
            sat = false;
            for (i = 0; i < pr->size_scons && !sat; i++)
                sat =
                    ap_abstract0_sat_lincons(pr->man_scons[i],
                    (a->scons) ? a->scons[i] : NULL,
                    &ncons);
            break;
        }
        default:
            sat = false;
            break;
    }
    /* do not free pcons since they are hashed */
    return sat;
}

/*
 * Needed for checking constraints representing aliasing between pointer
 * variables, e.g., x*next = y, in the assert statements.
 */
bool
ushape_sat_tcons(ap_manager_t * man, ushape_t * a, ap_tcons0_t * cons) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_SAT_TCONS, 0);
    if (ushape_is_bottom(man, a))
        return false;
    /* transform the tree constraint into a pointer constraint */
    pcons0_t *pcons = shape_pcons_of_tcons(pr, cons, a->datadim, a->ptrdim);
    /* check the constraint on the hgraph and the data/length constraints */
    bool sat = false;
    size_t i;
    switch (pcons->type) {
        case EQ_CONS:
        case NE_CONS:
            sat = hgraph_sat_pcons(pr, a->h, pcons);
            break;
        case DATA_CONS:
        {
            size_t *v2n = hgraph_get_var2node(a->h);
            ap_lincons0_t ncons = // set kind (scalar) to data constraint
                    shape_lincons_of_node(pr, &(pcons->info.data.cons),
                    pcons->info.data.offsets, v2n,
                    a->h->size, a->datadim, a->ptrdim);
            free(v2n);
            sat = false;
            for (i = 0; i < pr->size_scons && !sat; i++)
                sat =
                    ap_abstract0_sat_lincons(pr->man_scons[i],
                    (a->scons) ? a->scons[i] : NULL,
                    &ncons);
            ap_lincons0_clear(&ncons);
            break;
        }
        default:
            /* cannot be reach constraints!! */
            sat = false;
            break;
    }
    /* do not free pcons since they are hashed */
    return sat;
}

/* Interval constraints are only between non-pointer variables */
bool
ushape_sat_interval(ap_manager_t * man, ushape_t * a,
        ap_dim_t dim, ap_interval_t * itv) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_SAT_INTERVAL, 0);
    arg_assert(dim < a->datadim, return false;
            );
    size_t i;
    bool sat = true;
    for (i = 0; i < pr->size_scons && sat; i++)
        sat = ap_abstract0_sat_interval(pr->man_scons[i], a->scons[i], dim, itv);
    return sat;
}

/* Check on the hgraph or the existential constraint */
bool
ushape_is_dimension_unconstrained(ap_manager_t * man, ushape_t * a,
        ap_dim_t dim) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_IS_DIMENSION_UNCONSTRAINED, 0);
    arg_assert(a, return false;
            );
    size_t ndim = (dim < a->datadim) ? dim : (a->datadim + VAR2NODE(a->h,
            DIM2PTR
            (dim,
            a->
            datadim)));
    bool r = true;
    size_t i;
    if (ndim >= a->datadim)
        r = hgraph_is_dimension_unconstrained(man, a->h, dim);
    for (i = 0; i < pr->size_scons && r; i++)
        r =
            ap_abstract0_is_dimension_unconstrained(pr->man_scons[i],
            (a->scons) ? a->
            scons[i] : NULL, ndim);
    return r;
}

/* ============================================================ */
/* Extraction of properties */
/* ============================================================ */

/* TODO: used? */
ap_interval_t *
ushape_bound_linexpr(ap_manager_t * man, ushape_t * a, ap_linexpr0_t * expr) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_BOUND_LINEXPR, 0);
    arg_assert(a && a->h, return NULL;
            );
    /* applied only for the existential part */
    /* the dimensions shall be changed */
    size_t *v2n = hgraph_get_var2node(a->h);
    ap_linexpr0_t *nexpr =
            shape_linexpr_of_node(pr, expr, v2n, a->h->size, a->datadim,
            a->ptrdim);
    free(v2n);
    ap_interval_t *itv = NULL;
    size_t i;
    for (i = 0; i < pr->size_scons && !itv; i++)
        itv = ap_abstract0_bound_linexpr(pr->man_scons[i], a->scons[i], nexpr);
    /* TODO: give the min interval? */
    ap_linexpr0_free(nexpr);
    return itv;
}

/* NOT IMPLEMENTED */
ap_interval_t *
ushape_bound_texpr(ap_manager_t * man, ushape_t * a, ap_texpr0_t * expr) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_BOUND_TEXPR, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/* NOT IMPLEMENTED */
ap_interval_t *
ushape_bound_dimension(ap_manager_t * man, ushape_t * a, ap_dim_t dim) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_BOUND_DIMENSION, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/* NOT IMPLEMENTED */
ap_lincons0_array_t
ushape_to_lincons_array(ap_manager_t * man, ushape_t * a) {
    ap_lincons0_array_t ar;
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_TO_LINCONS_ARRAY, 0);
    ar = ap_lincons0_array_make(1);
    ar.p[0] = ap_lincons0_make_unsat();
    return ar;
}

/* NOT IMPLEMENTED */
ap_tcons0_array_t
ushape_to_tcons_array(ap_manager_t * man, ushape_t * a) {
    return ap_generic_to_tcons_array(man, a);
}

/* NOT IMPLEMENTED */
ap_interval_t **
ushape_to_box(ap_manager_t * man, ushape_t * a) {
    ushape_internal_t *pr = ushape_init_from_manager(man, AP_FUNID_TO_BOX, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return NULL;
}

/* NOT IMPLEMENTED */
ap_generator0_array_t
ushape_to_generator_array(ap_manager_t * man, ushape_t * a) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_TO_GENERATOR_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return ap_generator0_array_make(0);
}
