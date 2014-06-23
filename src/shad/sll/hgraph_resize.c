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


#include "hgraph.h"
#include "hgraph_internal.h"
#include "shape_macros.h"


/* ============================================================ */
/* Projections */
/* ============================================================ */

/* TODO: priority 0 */

/*
 * not used because we suppose that all ptr variables are declared from the
 * beginning
 */
hgraph_t *
hgraph_forget_array(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t * tdim, size_t size, bool project) {
    hgraph_internal_t *pr =
            hgraph_init_from_manager(man, AP_FUNID_FORGET_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}


/* ============================================================ */
/* Change and permutation of dimensions */
/* ============================================================ */

/**
 * Clone all nodes except NODE_NULL and label them in dimchange.
 * @requires  dimchange add dimensions at the end of the slice
 * @requires  perm has the final size of the graph
 * @requires  eqrel has size the initial size of the graph
 */
hgraph_t *
hgraph_dup_dimensions(hgraph_internal_t *pr, hgraph_t *a,
        ap_dimchange_t* dimchange, ap_dimperm_t* perm, ap_dimperm_t* eqrel) {

#ifndef NDEBUG1
    fprintf(stdout, "!!!!hgraph_dup_dimensions: with dimchange=(");
    ap_dimchange_fprint(stdout, dimchange);
    fprintf(stdout, ") and a->ptrdim = %zu\n", a->ptrdim);
    fflush(stdout);
#endif

    if (!a)
        return a;

    /* Requirement on dimchange*/
    arg_assert(dimchange &&
            (dimchange->realdim == a->ptrdim) &&
            ((dimchange->realdim > 0) ?
            ((dimchange->dim[dimchange->intdim] == (a->datadim + a->ptrdim)) ||
            (dimchange->dim[dimchange->intdim] == a->datadim))
            : 1),
    return NULL;);

    /* Requirement on perm */
    arg_assert(perm &&
            (perm->size == (2 * a->size - 1)),
    return NULL;);

    /* Requirement on eqrel */
    arg_assert(eqrel &&
            (eqrel->size == a->size),
    return NULL;);

    /* Step 0: special case, no pointer variables other than NULL */
    if (a->ptrdim == 0) {
        hgraph_t* r = hgraph_top(pr->man, 2 * a->datadim, 0);
        return r;
    }
    /* Step 1: alloc new graph with 2*a->size -1 nodes
     */
    size_t newsize = perm->size;
    size_t newdim = 2 * a->ptrdim;
    hgraph_t* r = hgraph_alloc_internal(pr, newsize, 2 * a->datadim, newdim);
    // copy ptrdim from a and clone them
    size_t n, v;
    for (v = 0; v < a->ptrdim; v++) {
        node_info_copy(&r->info[v], &a->info[v]);
        size_t dupv = a->ptrdim + v;
        VAR2NODE(r, dupv) = (VAR2NODE(a, v) == NODE_NULL) ? NODE_NULL
                : (VAR2NODE(a, v) + a->size - 1); // -1 because NODE_NULL not copied
    }
    // copy size from a and clone them
    eqrel->dim[0] = 0; // maps NODE_NULL to itself
    for (n = 1; n < a->size; n++) {
        node_info_copy(&r->info[n + newdim], &a->info[n + a->ptrdim]);
        size_t dupn = a->size - 1 + n;
        eqrel->dim[n] = dupn;
        node_info_copy(&r->info[dupn + newdim], &a->info[n + a->ptrdim]);
        NODE_VAR(r, dupn) = NODE_VAR(a, n) + a->ptrdim;
        NODE_NEXT(r, dupn) = (NODE_NEXT(a,n) == NODE_NULL) ? NODE_NULL :
                                (NODE_NEXT(a, n) + a->size - 1);
    }

    /* Step 3: normally, the graph generated is already sorted,
     * but to be sure call sorting and test identity.
     */
    hgraph_node_sort(r, 1, perm);
    // arg_assert (shape_dimperm_is_id(perm), return NULL;);

#ifndef NDEBUG1
    fprintf(stdout, "!!!!hgraph_dup_dimensions returns: (");
    hgraph_fdump(stdout, pr->man, r);
    fprintf(stdout, ")\n perm= ");
    ap_dimperm_fprint(stdout, perm);
    for (n = 0; n < a->size; n++)
        fprintf(stdout, "eqrel[%zu] = %zu\n", n, eqrel->dim[n]);
    fflush(stdout);
#endif
    return r;
}

/**
 * Add pointer dimensions du the hgraph. and initialize them to NODE_NULL.
 */
hgraph_t *
hgraph_add_dimensions(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dimchange_t * dimchange, bool project) {

    hgraph_internal_t *pr =
            hgraph_init_from_manager(man, AP_FUNID_ADD_DIMENSIONS, 0);

    if (!a)
        return a;

    hgraph_t *r;
    size_t nptrdim = a->ptrdim + dimchange->realdim;
    size_t ndatadim = a->datadim + dimchange->intdim;
    if (dimchange->realdim == 0) {
        r = hgraph_copy_mem(pr, a);
        r->datadim = ndatadim;
    } else {
        r = hgraph_alloc_internal(pr, a->size, ndatadim, nptrdim);
        size_t * newvars; // maps old vars in new vars
        checked_malloc(newvars, size_t, sizeof (size_t), a->ptrdim, return NULL;
                );
        /* start copying the ptrdim part of info */
        size_t i; // iterates over the new info
        size_t j; // index in the dimchange->dim
        size_t k; // the offset of i wrt the a->info
        for (i = 0, j = 0, k = 0; j < dimchange->realdim || i < nptrdim; i++)
            if (j < dimchange->realdim
                    && ((dimchange->dim[j + dimchange->intdim] - a->datadim) == (i - k))) {
                // introduce a new dimension in i mapped to NULL if not unify
                if (pr->meet_algo == 0)
                    VAR2NODE(r, i) = NODE_NULL;
                else
                    VAR2NODE(r, i) = NODE_T_TOP;
                k++; // increment offset
                j++;
            } else {
                // copy the old dimension
                node_info_copy(&r->info[i], &a->info[i - k]);
                newvars[i - k] = i;
            }
        // go through the nodes and update old variables (except NODE_NULL)
        node_info_copy(&r->info[nptrdim], &a->info[a->ptrdim]);
        for (i = 1; i < a->size; i++) {
            node_info_copy(&r->info[i + nptrdim], &a->info[i + a->ptrdim]);
            NODE_VAR(r, i) = newvars[NODE_VAR(a, i)];
        }
        free(newvars);
    }
    if (destructive)
        hgraph_free_internal(pr, a);

    // NO need to sort because all added ptr variables label NODE_NULL
    hgraph_t* rr = hgraph_copy_internal(pr, r);
    hgraph_free_internal(pr, r);
    return rr;
}

/**
 * Remove real dimensions from dimchange and reorder nodes.
 * @requires  No anonymous nodes are generated
 *            which is true because this function is used
 *            only in procedure call/return.
 * @return    Cutpoints generated are signaled...
 */
hgraph_t *
hgraph_remove_dimensions_internal(hgraph_internal_t *pr, bool destructive,
        hgraph_t* a, ap_dimchange_t* dimchange, ap_dimperm_t* perm) {
    if (!a)
        return a;

    hgraph_t *r = hgraph_copy_mem(pr, a);
    hgraph_t *rr = NULL;

    r->datadim = a->datadim - dimchange->intdim;
    if (dimchange->realdim > 0) {
        size_t i;
        // Step 1: set all removed dimensions to NULL
        for (i = 0; i < dimchange->realdim; i++) {
            size_t v = dimchange->dim[dimchange->intdim + i] - a->datadim;
            VAR2NODE(r, v) = NODE_NULL;
        }

        // Step 2: call closure with parameter true to signal cutpoints
        ap_dim_array2_t anon;
        anon.size = 0;
        anon.p = NULL;
        rr = hgraph_close(pr, r, perm, &anon, true);

        // more anonymous cannot be generated
        arg_assert(anon.size == 0, return NULL;);

        hgraph_node_anon_array_clear(pr, &anon);
        // some nodes (!= NODE_NULL) are mapped to 0 in perm!
        hgraph_free_internal(pr, r);

        // Step 3: effectively remove ptrdim from info
        // No reordering is needed since eliminated ptr are put to NODE_NULL
        size_t nptrdim;
        nptrdim = rr->ptrdim - dimchange->realdim;
        r = hgraph_alloc_internal(pr, rr->size, rr->datadim, nptrdim);
        size_t * newvars; // maps old vars in new vars
        checked_malloc(newvars, size_t, sizeof (size_t), rr->ptrdim, return NULL;
                );
        /* start copying the ptrdim part of info */
        // i iterates over the new info
        size_t j; // index in the dimchange->dim
        size_t k; // the offset of i wrt the a->info
        for (i = 0, j = 0, k = 0; j < dimchange->realdim || i < nptrdim;)
            if (j < dimchange->realdim
                    && ((dimchange->dim[j + dimchange->intdim] - a->datadim) == (i + k))) {
                // ignore this dimension
                k++; // increment the offset
                j++;
            } else {
                // copy the old dimension
                node_info_copy(&r->info[i], &rr->info[i + k]);
                newvars[i + k] = i;
                i++;
            }
        // go through the nodes (except NODE_NULL) and update old variables
        node_info_copy(&r->info[nptrdim], &rr->info[rr->ptrdim]);
        for (i = 1; i < rr->size; i++) {
            node_info_copy(&r->info[i + nptrdim], &rr->info[i + rr->ptrdim]);
            NODE_VAR(r, i) = newvars[NODE_VAR(rr, i)];
        }
        free(newvars);
        hgraph_free_internal(pr, rr);
        if (destructive)
            hgraph_free_internal(pr, a);
    }
    rr = hgraph_copy_internal(pr, r);
    hgraph_free_internal(pr, r);
    return rr;
}

hgraph_t *
hgraph_remove_dimensions(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dimchange_t * dimchange) {
    hgraph_internal_t *pr =
            hgraph_init_from_manager(man, AP_FUNID_REMOVE_DIMENSIONS, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}

/* Permutes pointer variables and reorder nodes
 * inperm->size==a->ptrdim and inperm->dim contains only ptrdim dimensions
 */
hgraph_t *
hgraph_permute_dimensions_internal(hgraph_internal_t *pr, bool destructive,
        hgraph_t* a, ap_dimperm_t* inperm, ap_dimperm_t* nodeperm) {
    if (!a)
        return a;

    hgraph_t *r = hgraph_copy_mem(pr, a);

#ifndef NDEBUG1
    fprintf(stdout, "!!!!hgraph_permute_dimensions_internal: on a=(");
    hgraph_fdump(stdout, pr->man, a);
    fprintf(stdout, ") and inperm=(");
    ap_dimperm_fprint(stdout, inperm);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif
    // Step 1: permute dimensions in info
    size_t v;
    for (v = 0; v < a->ptrdim; v++)
        // the only thing to copy is the node labeled, i.e.,
        // the new dimension receives the node of the old (permuted dimension)
        VAR2NODE(r, inperm->dim[v]) = VAR2NODE(a, v);

#ifndef NDEBUG1
    fprintf(stdout, "!!!!hgraph_permute_dimensions_internal: after permuting dimensions r=(");
    hgraph_fdump(stdout, pr->man, r);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    // Step 2: compute new info for nodes
    ap_dim_array2_t anon;
    anon.size = 0;
    anon.p = NULL;
    hgraph_t* rr = hgraph_close(pr, r, nodeperm, &anon, false);
    if (anon.size > 0) {
        ERROR("Bad dimension permutation in hgraph!", return NULL;);
        free(anon.p);
    }
    // free the intermediate memory
    hgraph_free_internal(pr, r);
    if (destructive)
        hgraph_free_internal(pr, a);
#ifndef NDEBUG1
    fprintf(stdout, "!!!!hgraph_permute_dimensions_internal: returns r=(");
    hgraph_fdump(stdout, pr->man, rr);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif
    return rr; /* not recorded in hash table! */
}

hgraph_t *
hgraph_permute_dimensions(ap_manager_t * man, /* TODO: priority 0 */
        bool destructive, hgraph_t * a,
        ap_dimperm_t * permutation) {
    hgraph_internal_t *pr =
            hgraph_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}


/* ============================================================ */
/* Expansion and folding of dimensions */

/* ============================================================ */

hgraph_t *
hgraph_expand(ap_manager_t * man, /* TODO: priority 0 */
        bool destructive, hgraph_t * a, ap_dim_t dim, size_t n) {
    hgraph_internal_t *pr = hgraph_init_from_manager(man, AP_FUNID_EXPAND, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}

/* TODO: priority 0 */
hgraph_t *
hgraph_fold(ap_manager_t * man,
        bool destructive, hgraph_t * a, ap_dim_t * tdim, size_t size) {
    hgraph_internal_t *pr = hgraph_init_from_manager(man, AP_FUNID_FOLD, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}
