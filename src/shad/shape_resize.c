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


/* ============================================================ */
/* Projections */
/* ============================================================ */

/* TODO: priority 0 */

/* not used */
shape_t *
shape_forget_array(ap_manager_t * man,
        bool destructive, shape_t * a,
        ap_dim_t * tdim, size_t size, bool project) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_FORGET_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}


/* ============================================================ */
/* Change and permutation of dimensions */

/* ============================================================ */

shape_t *
shape_add_dimensions(ap_manager_t * man,
        bool destructive, shape_t * a,
        ap_dimchange_t * dimchange, bool project) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_ADD_DIMENSIONS, 0);

#ifndef NDEBUG
    fprintf(stdout, "\n****shape_add_dimensions: dimchange=(");
    ap_dimchange_fprint(stdout, dimchange);
    fprintf(stdout, ") on ");
    shape_fdump(stdout, man, a);
    fprintf(stdout, "\n");
    fprintf(stdout, " with ");
    fflush(stdout);
#endif
    /* nothing to do */
    if (!a)
        return NULL;
    if (shape_is_bottom(man, a))
        if (destructive) {
            a->intdim += dimchange->intdim;
            //a->realdim += dimchange->realdim;
            a->realdim += dimchange->intdim;
            return a;
        } else
            return shape_bottom(man,
                a->intdim + dimchange->intdim,
                a->realdim + dimchange->realdim);
    else
        if (!dimchange || (dimchange->intdim + dimchange->realdim) == 0)
        return (destructive || !a) ? a : shape_copy(man, a);

    shape_t *b;
    /* non trivial case: add dimension to each ushapes */
    if (!destructive)
        b = shape_copy_internal(pr, a);
    else
        b = a;
    /* go */
    ushape_array_t arr;
    arr.p = NULL;
    arr.size = 0;
    ushape_array_init(pr, &arr, a->m.size);
    shape_t *r;
    size_t i, size;
    size = 0;
#ifndef NDEBUG
    fprintf(stdout, "\n****shape_add_dimensions: before scons\n");
    fflush(stdout);
#endif
    for (i = 0; i < a->msize; i++) {
        ushape_t *ur =
                ushape_add_dimensions(man, false, a->m.p[i], dimchange, project);
        if (ur && !ushape_is_bottom(man, ur))
            size += ushape_array_add(pr, true, &arr, size, false, false, ur); /* do not copy and
											 * destroy */
        else if (!ur)
            ushape_free(man, ur);
    }
#ifndef NDEBUG
    fprintf(stdout, "\n****shape_add_dimensions: before returns\n");
    fflush(stdout);
#endif
    if (size == 0)
        r = shape_bottom(man, a->intdim + dimchange->intdim, a->realdim + dimchange->realdim);
    else {
        ushape_array_t *rr;
        r = shape_alloc_internal(pr, size);
        rr = ushape_array_copy(pr, &arr, size);
        r->m = *rr;
        rr->p = NULL;
        free(rr);
        r->set = true;
        r->closed = false;
        r->msize = size;
        r->intdim = a->intdim + dimchange->intdim;
        r->realdim = a->realdim + dimchange->realdim;
    }
    if (destructive)
        shape_free_internal(pr, a);
#ifndef NDEBUG1
    fprintf(stdout, "\n****shape_add_dimensions returns: ");
    shape_fdump(stdout, man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return r;
}

shape_t *
shape_remove_dimensions(ap_manager_t * man,
        bool destructive, shape_t * a,
        ap_dimchange_t * dimchange) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_REMOVE_DIMENSIONS, 0);

#ifndef NDEBUG
    fprintf(stdout, "\n****shape_remove_dimensions:  dimchang=(");
    ap_dimchange_fprint(stdout, dimchange);
    fprintf(stdout, ") on ");
    shape_fdump(stdout, man, a);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    shape_t* r;
    /* nothing to do */
    if (!a) {
        r = NULL;
    } else if (shape_is_bottom(man, a)) {
        if (destructive) {
            a->intdim -= dimchange->intdim;
            a->realdim -= dimchange->realdim;
            r = a;
        } else
            r = shape_bottom(man, a->intdim - dimchange->intdim, a->realdim - dimchange->realdim);
    } else if (!dimchange || (dimchange->intdim + dimchange->realdim) == 0)
        r = (destructive || !a) ? a : shape_copy(man, a);
    else {
        shape_t *b;
        /* non trivial case: remove dimension to each ushapes */
        if (!destructive)
            b = shape_copy_internal(pr, a);
        else
            b = a;
        /* go */
        ushape_array_t arr;
        arr.p = NULL;
        arr.size = 0;
        ushape_array_init(pr, &arr, a->m.size);
        size_t i, size;
        size = 0;
        for (i = 0; i < a->msize; i++) {
            ushape_t *ur =
                    ushape_remove_dimensions(man, false, a->m.p[i], dimchange);
            if (ur && !ushape_is_bottom(man, ur))
                size += ushape_array_add(pr, true, &arr, size, false, false, ur); /* do not copy and
											 * destroy */
            else if (!ur)
                ushape_free(man, ur);
        }
#ifndef NDEBUG
        fprintf(stdout, "\n****shape_remove_dimensions: before returns\n");
        fflush(stdout);
#endif
        if (size == 0)
            r = shape_bottom(man, a->intdim - dimchange->intdim, a->realdim - dimchange->realdim);
        else {
            ushape_array_t *rr;
            r = shape_alloc_internal(pr, size);
            rr = ushape_array_copy(pr, &arr, size);
            r->m = *rr;
            rr->p = NULL;
            free(rr);
            r->set = true;
            r->closed = false;
            r->msize = size;
            r->intdim = a->intdim - dimchange->intdim;
            r->realdim = a->realdim - dimchange->realdim;
        }
        if (destructive)
            shape_free_internal(pr, a);
    }
#ifndef NDEBUG
    fprintf(stdout, "\n****shape_remove_dimensions returns: ");
    if (!r)
        fprintf(stdout, " NULL");
    else
        shape_fdump(stdout, man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return r;
}

shape_t *
shape_permute_dimensions(ap_manager_t * man,
        bool destructive, shape_t * a,
        ap_dimperm_t * permutation) {
    shape_internal_t *pr =
            shape_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);

#ifndef NDEBUG1
    fprintf(stdout, "\n****shape_permute_dimension: dimperm=(");
    ap_dimperm_fprint(stdout, permutation);
    fprintf(stdout, ") on ");
    shape_fdump(stdout, man, a);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif

    /* nothing to do */
    if (!a)
        return NULL;
    if (shape_is_bottom(man, a))
        if (destructive) {
            return a;
        } else
            return shape_bottom(man, a->intdim, a->realdim);
    else if (!permutation || permutation->size == 0 || shape_dimperm_is_id(permutation))
        return (destructive || !a) ? a : shape_copy(man, a);

    shape_t *b;
    shape_t *r;
    /* non trivial case: permute dimension to each ushapes */
    if (!destructive)
        b = shape_copy_internal(pr, a);
    else
        b = a;
    /* go */
    ushape_array_t arr;
    arr.p = NULL;
    arr.size = 0;
    ushape_array_init(pr, &arr, a->m.size);
    size_t i, size;
    size = 0;
#ifndef NDEBUG
    fprintf(stdout, "\n****shape_permute_dimensions: before scons\n");
    fflush(stdout);
#endif
    for (i = 0; i < a->msize; i++) {
        ushape_t *ur =
                ushape_permute_dimensions(man, false, a->m.p[i], permutation);
        if (ur && !ushape_is_bottom(man, ur))
            size += ushape_array_add(pr, true, &arr, size, false, false, ur); /* do not copy and
											 * destroy */
        else if (!ur)
            ushape_free(man, ur);
    }
#ifndef NDEBUG
    fprintf(stdout, "\n****shape_permute_dimensions: before returns\n");
    fflush(stdout);
#endif
    if (size == 0)
        r = shape_bottom(man, a->intdim, a->realdim);
    else {
        ushape_array_t *rr;
        r = shape_alloc_internal(pr, size);
        rr = ushape_array_copy(pr, &arr, size);
        r->m = *rr;
        rr->p = NULL;
        free(rr);
        r->set = true;
        r->closed = false;
        r->msize = size;
        r->intdim = a->intdim;
        r->realdim = a->realdim;
    }
    if (destructive)
        shape_free_internal(pr, a);
#ifndef NDEBUG1
    fprintf(stdout, "\n****shape_permute_dimensions returns: ");
    shape_fdump(stdout, man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return r;
}


/* ============================================================ */
/* Expansion and folding of dimensions */
/* ============================================================ */

/* TODO: priority 0 */

/* not used */
shape_t *
shape_expand(ap_manager_t * man,
        bool destructive, shape_t * a, ap_dim_t dim, size_t n) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_EXPAND, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}

/* TODO: priority 0 */
shape_t *
shape_fold(ap_manager_t * man,
        bool destructive, shape_t * a, ap_dim_t * tdim, size_t size) {
    shape_internal_t *pr = shape_init_from_manager(man, AP_FUNID_FOLD, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}
