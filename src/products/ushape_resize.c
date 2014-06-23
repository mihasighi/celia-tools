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


#include "shape_macros.h"
#include "ushape.h"
#include "ushape_internal.h"


/* ============================================================ */
/* Projections */
/* ============================================================ */

/**
 * Not implemented because we do forward analysis.
 */
ushape_t *
ushape_forget_array(ap_manager_t * man, /* TODO: priority 0 */
        bool destructive, ushape_t * a,
        ap_dim_t * tdim, size_t size, bool project) {
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_FORGET_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}


/* ============================================================ */
/* Change and permutation of dimensions */
/* ============================================================ */

/**
 * Clone existing nodes and dimensions and
 * introduce equalities between clones.
 * @requires: the dimchange contains dimensions which are all 
 *            at the beginning / end of each slice
 * @requires: the cloned dimensions are all equal to the first/last one
 * @returns:  the permutation is an identity
 * @returns:  the eqrel array has a->size and contains cloned nodes
 */
ushape_t*
ushape_dup_dimensions(ushape_internal_t* pr,
        ushape_t * a,
        ap_dimchange_t * dimchange) {

#ifndef NDEBUG1
    fprintf(stdout, "\n====ushape_dup_dimensions: on=(");
    ushape_fdump(stdout, pr->man, a);
    fprintf(stdout, ") with dimchange=(");
    ap_dimchange_fprint(stdout, dimchange);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    arg_assert(dimchange &&
            (dimchange->intdim == a->datadim) &&
            (dimchange->realdim == a->ptrdim) &&
            ((dimchange->intdim > 0) ?
            ((dimchange->dim[0] == a->datadim) ||
            (dimchange->dim[0] == 0)) : 1) &&
            ((dimchange->realdim > 0) ?
            ((dimchange->dim[dimchange->intdim] == a->datadim) ||
            (dimchange->dim[dimchange->intdim] == (a->datadim + a->ptrdim))) : 1),
    return NULL;);

    /* Step 1: apply clone on the hgraph.
     *    Generate a permutation of nodes labeled by
     *    the new dimensions.
     */
    ap_dimperm_t perm; // permuation generated, identity
    ap_dimperm_t eqrel; // equality between nodes (already permuted)
    hgraph_t *h;
    ap_dimperm_init(&perm, 2 * a->h->size - 1);
    ap_dimperm_set_id(&perm);
    ap_dimperm_init(&eqrel, a->h->size);
    ap_dimperm_set_id(&eqrel);
    h = hgraph_dup_dimensions(pr, a->h, dimchange, &perm, &eqrel);
    // initial nodes shall be not permuted

    ushape_t* b = ushape_alloc_internal(pr, a->datadim + dimchange->intdim,
            a->ptrdim + dimchange->realdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    b->closed = NULL;

    /* Step 2: apply dimchange where realdim are the cloned nodes.
     */
    size_t nrealdim = a->h->size - 1; // number of nodes added
    ap_dimchange_t newchange;
    ap_dimchange_init(&newchange, dimchange->intdim, nrealdim);
    size_t i;
    for (i = 0; i < dimchange->intdim; i++)
        newchange.dim[i] = dimchange->dim[i];
    for (i = 0; i < nrealdim; i++)
        newchange.dim[dimchange->intdim + i] = dimchange->intdim + a->h->size;
    // new dimensions are added after the existing ones
    ushape_apply_dimperm_dimchange_n(pr, a, b, &perm, &newchange);
    // normally, duplicated dimensions are a block at the end / begining ...

    /* Step 3: introduce equality predicates on constraints on data words
     *    i.e. dimension i equals dimension a->datadim+i for i < intdim
     *    node perm->dim[i] equals node perm.dim[eqrel.dim[i]]
     */
#ifndef NDEBUG2
    fprintf(stdout, "\nstructural equality with (datadim=%zu,segmdim=%zu)\n", b->datadim, b->h->size);
    fflush(stdout);
#endif
    ap_lincons0_array_t eqcons = ap_lincons0_array_make(a->datadim + a->h->size - 1);
    // NODE_NULL not duplicated
    for (i = 0; i < a->datadim; i++) {
        // assigned to (i)-(i+a->datadim) = 0
        ap_scalar_t * zero = ap_scalar_alloc();
        ap_scalar_set_int(zero, OFFSET_DATA);
        ap_linexpr0_t* lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, b->datadim + b->h->size);
        ap_linexpr0_set_coeff_scalar_int(lexpr, i, 1);
        ap_linexpr0_set_coeff_scalar_int(lexpr, i + a->datadim, -1);
#ifndef NDEBUG2
        fprintf(stdout, "\n\t- (%d = %d) builds ", i, a->datadim + i);
        ap_linexpr0_fprint(stdout, lexpr, NULL);
        fprintf(stdout, "\n");
        fflush(stdout);
#endif
        eqcons.p[i] = ap_lincons0_make(AP_CONS_EQ, lexpr, zero);
    }
    for (i = 1; i < a->h->size; i++) {
        // assigned to i-eqrel.dim[i] mod(OFFSET_OTHER) 0
        ap_scalar_t * zero = ap_scalar_alloc();
        ap_scalar_set_int(zero, OFFSET_NONE);
        ap_linexpr0_t* lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, b->datadim + b->h->size);
        ap_linexpr0_set_coeff_scalar_int(lexpr, b->datadim + i, -1); // PB i-1
        ap_linexpr0_set_coeff_scalar_int(lexpr, b->datadim + eqrel.dim[i], 1); // PB -1
#ifndef NDEBUG2
        fprintf(stdout, "\n\t- (%d = %d) builds ", b->datadim + i, b->datadim + eqrel.dim[i]);
        ap_linexpr0_fprint(stdout, lexpr, NULL);
        fprintf(stdout, "\n");
        fflush(stdout);
#endif
        eqcons.p[a->datadim + i - 1] = ap_lincons0_make(AP_CONS_EQ, lexpr, zero);
    }

    ushape_t* r = ushape_alloc_internal(pr, b->datadim, b->ptrdim);
    r->h = hgraph_copy_internal(pr, b->h);
    r->closed = NULL;
    for (i = 0; i < pr->size_scons; i++)
        r->scons[i] = ap_abstract0_meet_lincons_array(pr->man_scons[i], false, b->scons[i],
            &eqcons);

    // free allocated resources
    ap_dimperm_clear(&perm);
    ap_dimperm_clear(&eqrel);
    ap_dimchange_clear(&newchange);
    ap_lincons0_array_clear(&eqcons);
    ushape_free_internal(pr, b);

#ifndef NDEBUG1
    fprintf(stdout, "\n====ushape_dup_dimensions returns: (");
    ushape_fdump(stdout, pr->man, r);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif
    return r;
}

/* Add the dimensions in constraints using @code{dimchange}. */
void
ushape_addrem_dimensions_scons(ushape_internal_t* pr,
        ushape_t* a,
        ushape_t* b,
        ap_dimchange_t* dimchange, bool isadd) {
    size_t i;
    ap_dimchange_t dchange;
    assert(a && a->h && b->h && (dimchange->realdim < b->h->size));
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_addrem_dimensions_scons: with dimchange=%c(",
            (isadd) ? '+' : '-');
    ap_dimchange_fprint(stdout, dimchange);
    fprintf(stdout, ")\n");
#endif
    ap_dimchange_init(&dchange, dimchange->intdim, dimchange->realdim);
    // copy data dimensions
    for (i = 0; i < dimchange->intdim; i++)
        dchange.dim[i] = dimchange->dim[i];
    // copy node dimensions
    for (i = 0; i < dimchange->realdim; i++) {
        // TODO: remove test below after validation
        arg_assert(a->datadim <= dimchange->dim[i + dimchange->intdim] &&
                dimchange->dim[i + dimchange->intdim] < a->datadim + a->h->size,
                ;);
        // always add dimensions in scons at the end of the existing dimensions???
        dchange.dim[dimchange->intdim + i] = dimchange->dim[i + dimchange->intdim];
    }
    for (i = 0; i < pr->size_scons; i++) {
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_addrem_dimensions_scons: with dchange=(\n");
        ap_dimchange_fprint(stdout, &dchange);
        fprintf(stdout, ") on scons=(");
        ap_abstract0_fprint(stdout, pr->man_scons[i], a->scons[i], NULL);
        fprintf(stdout, ")\n");
#endif
        if (isadd)
            b->scons[i] =
                ap_abstract0_add_dimensions(pr->man_scons[i],
                false, a->scons[i], &dchange, false);
        else
            b->scons[i] =
                ap_abstract0_remove_dimensions(pr->man_scons[i],
                false, a->scons[i], &dchange);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_addrem_dimensions_scons: after dchange scons=(\n");
        ap_abstract0_fprint(stdout, pr->man_scons[i], b->scons[i], NULL);
        fprintf(stdout, ")\n");
#endif
    }
    ap_dimchange_clear(&dchange);
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_addrem_dimensions_scons returns\n");
    fflush(stdout);
#endif
}

/* Depending on @code{project}, this functions has two behaviours:
 * - if project and dimchange->intdim=a->datadim and dimchange->realdim=a->ptrdim
 *   then duplicate all dimensions
 * - if not(project) simply add dimensions in dimchange
 */
ushape_t *
ushape_add_dimensions(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dimchange_t * dimchange, bool project) {
    if (!a)
        return NULL;
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_ADD_DIMENSIONS, 0);
    /* nothing to do */
    if (dimchange->intdim + dimchange->realdim == 0)
        return (destructive) ? a : ushape_copy_internal(pr, a);

    ushape_t *b;

    if (project)
        b = ushape_dup_dimensions(pr, a, dimchange);
    else {
        /* normal procedure */
        /* apply dimchange on hgraph */
        hgraph_t *h = hgraph_add_dimensions(man, false, a->h, dimchange, false);
        // nodes are not permuted because all new dimensions are set to NODE_NULL

        b = ushape_alloc_internal(pr, a->datadim + dimchange->intdim, a->ptrdim + dimchange->realdim);
        b->h = hgraph_copy_internal(pr, h);
        hgraph_free_internal(pr, h);
        b->closed = NULL;
        if (dimchange->intdim > 0) {
            // since no new nodes are added, dimchange shall be reset to 0 for realdim
            size_t ptrdim = dimchange->realdim;
            dimchange->realdim = 0;
            ushape_addrem_dimensions_scons(pr, a, b, dimchange, true);
            dimchange->realdim = ptrdim;
        } else {
            // copy constraints from a
            if (a->scons) {
                size_t i;
                for (i = 0; i < pr->size_scons; i++)
                    if (a->scons[i])
                        b->scons[i] = ap_abstract0_copy(pr->man_scons[i], a->scons[i]);
            }
        }
    }
    if (destructive)
        ushape_free_internal(pr, a);
    return b;
}

ushape_t *
ushape_remove_dimensions(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dimchange_t * dimchange) {
    if (!a)
        return NULL;
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_REMOVE_DIMENSIONS, 0);
    /* nothing to do */
    if (dimchange->intdim + dimchange->realdim == 0)
        return (destructive) ? a : ushape_copy_internal(pr, a);

    hgraph_t* h;
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    /* apply dimchange on hgraph
     * ==> signals cutpoints, updates datadim
     * ==> generates a permutation of nodes, with some nodes mapped to 0
     * ==> but not anonymous since deleted ptr dimensions are supposed to
     */
    h = hgraph_remove_dimensions_internal(pr, false, a->h, dimchange, &perm);

    ushape_t *b = ushape_alloc_internal(pr,
            a->datadim - dimchange->intdim,
            a->ptrdim - dimchange->realdim);

    b->h = hgraph_copy_internal(pr, h);
    b->closed = NULL;

    // first remove data dimensions
    if (dimchange->intdim > 0) {
        size_t ptrdim = dimchange->realdim;
        dimchange->realdim = 0;
        ushape_addrem_dimensions_scons(pr, a, b, dimchange, false);
        dimchange->realdim = ptrdim;
    } else {
        ushape_copy_internal_scons(pr, a, b);
    }

    // then, if something changed for nodes, apply permutation procedure (includes removing)
    bool toperm = false;
    size_t i;
    for (i = 1; i < perm.size && !toperm; i++)
        if (perm.dim[i] == 0 || perm.dim[i] != i)
            toperm = true;
    if (toperm)
        ushape_apply_dimperm(pr, b, b, &perm);

    ap_dimperm_clear(&perm);

    if (destructive)
        ushape_free_internal(pr, a);

    return b;
}

ushape_t *
ushape_permute_dimensions(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dimperm_t * permutation) {
    if (!a)
        return NULL;
    ushape_internal_t *pr =
            ushape_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
    /* nothing to do */
    if (!permutation->size)
        return (destructive) ? a : ushape_copy_internal(pr, a);

    /* apply permutation on hgraph ==> generates a permutation of nodes */
    ap_dimperm_t hperm; /* permutation for dimensions of hgraph */
    ap_dimperm_init(&hperm, a->ptrdim);
    ap_dimperm_set_id(&hperm);
    size_t i; // translate real dimension in ptr dimensions */
    for (i = a->datadim; i < permutation->size; i++)
        hperm.dim[i - a->datadim] = permutation->dim[i] - a->datadim;
    hgraph_t *h;
    ap_dimperm_t perm; /* node permutation */
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    h = hgraph_permute_dimensions_internal(pr, false, a->h, &hperm, &perm);

    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    b->closed = NULL;
    ap_dimperm_t dimperm;
    ap_dimperm_init(&dimperm, a->datadim + a->h->size);
    ap_dimperm_set_id(&dimperm);
    shape_dimperm_copy(&dimperm, a->datadim, &perm);
    for (i = 0; i < a->datadim; i++)
        dimperm.dim[i] = permutation->dim[i];

    ushape_apply_dimperm_n(pr, a, b, &dimperm);
    if (destructive)
        ushape_free_internal(pr, a);
    return b;
}

/* Apply only permutation of nodes, i.e., perm->size == a->h->size */
void
ushape_apply_dimperm(ushape_internal_t * pr, ushape_t * a, ushape_t * b,
        ap_dimperm_t * perm) {
    size_t i;
    assert(a && a->h && (perm->size >= a->h->size));
#ifndef NDEBUG1
    fprintf(stdout, "===== ushape_apply_dimperm: with perm=(");
    ap_dimperm_fprint(stdout, perm);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif
    if (shape_dimperm_is_id(perm)) {
        if (a != b) ushape_copy_internal_scons(pr, a, b);
    } else {
        ap_dimperm_t dimperm;
        ap_dimperm_init(&dimperm, a->datadim + perm->size);
        ap_dimperm_set_id(&dimperm);
        shape_dimperm_copy(&dimperm, a->datadim, perm);
        for (i = 0; i < pr->size_scons; i++)
            b->scons[i] =
                ap_abstract0_permute_dimensions(pr->man_scons[i], (a == b),
                a->scons[i], &dimperm);
        ap_dimperm_clear(&dimperm);
    }
    ap_dimperm_clear(perm);
}

/* Apply permutation of data variables AND nodes, i.e., perm->size == a->h->size+a->datadim */
void
ushape_apply_dimperm_n(ushape_internal_t * pr, ushape_t * a, ushape_t * b,
        ap_dimperm_t * perm) {
    size_t i;
    assert(a && a->h && (perm->size >= (a->h->size + a->datadim)));
    if (shape_dimperm_is_id(perm)) {
        if (a != b) ushape_copy_internal_scons(pr, a, b);
    } else {
        for (i = 0; i < pr->size_scons; i++)
            b->scons[i] =
                ap_abstract0_permute_dimensions(pr->man_scons[i], (a == b),
                a->scons[i], perm);
    }
    ap_dimperm_clear(perm);
}

/* Add ONE dimension and map it to n, then apply permutation of nodes,
 * i.e., perm->size == a->h->size */
void
ushape_apply_dimperm_dimchange(ushape_internal_t * pr, ushape_t * a,
        ushape_t * b, ap_dimperm_t * perm, size_t n) {
    size_t i;
    ap_dimperm_t dimperm;
    ap_dimchange_t dimchange;
    assert(a && a->h && (perm->size == a->h->size));
    ap_dimchange_init(&dimchange, 0, 1);
    dimchange.dim[0] = a->datadim + a->h->size;
    ap_dimperm_init(&dimperm, a->datadim + a->h->size + 1);
    ap_dimperm_set_id(&dimperm);
    shape_dimperm_copy(&dimperm, a->datadim, perm);
    ap_dimperm_clear(perm);
    dimperm.dim[a->datadim + a->h->size] = a->datadim + n;
    for (i = 0; i < pr->size_scons; i++) {
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_apply_dimperm_dimchange: with dimchange=(");
        ap_dimchange_fprint(stdout, &dimchange);
        fprintf(stdout, ") on scons=(\n");
        ap_abstract0_fprint(stdout, pr->man_scons[i], a->scons[i], NULL);
        fprintf(stdout, ")\n");
#endif
        b->scons[i] =
                ap_abstract0_add_dimensions(pr->man_scons[i],
                false, a->scons[i], &dimchange, false);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_apply_dimperm_dimchange: after dimchange scons=(\n");
        ap_abstract0_fprint(stdout, pr->man_scons[i], b->scons[i], NULL);
        fprintf(stdout, ")\n");
#endif
        // dimperm is clearly not an identity
        b->scons[i] =
                ap_abstract0_permute_dimensions(pr->man_scons[i],
                true, b->scons[i], &dimperm);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_apply_dimperm_dimchange: after dimperm scons=(\n");
        ap_abstract0_fprint(stdout, pr->man_scons[i], b->scons[i], NULL);
        fprintf(stdout, ")\n");
#endif
    }
    ap_dimperm_clear(&dimperm);
    ap_dimchange_clear(&dimchange);
}

/* Add dimensions in @code{change}, then apply permutation of nodes,
 * i.e., perm->size == b->h->size == a->hsize + change->realdim */
void
ushape_apply_dimperm_dimchange_n(ushape_internal_t * pr, ushape_t * a,
        ushape_t * b, ap_dimperm_t * perm, ap_dimchange_t* change) {
    size_t i;
    ap_dimperm_t dimperm;
    ap_dimchange_t dimchange;
    assert(a && a->h && b->h && (perm->size == b->h->size));
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_dimperm_dimchange_n: with perm=(\n");
    ap_dimperm_fprint(stdout, perm);
    fprintf(stdout, ") and change=(");

    ap_dimchange_fprint(stdout, change);
    fprintf(stdout, ")\n");
#endif
    /*
    ap_dimchange_init(&dimchange, change->intdim, change->realdim);
    for (i = 0; i < change->intdim; i++)
        dimchange.dim[i] = change->dim[i]; // or b->h->size?
    for (i = 0; i < change->realdim; i++)
        // always add dimensions in scons at the end of the existing dimensions???
        dimchange.dim[change->intdim + i] = change->dim[change->intdim + i]; // or b->h->size?
     */
    // always add dimensions in scons at the end of the existing dimensions???
    ap_dimperm_init(&dimperm, a->datadim + change->intdim + b->h->size);
    ap_dimperm_set_id(&dimperm);
    shape_dimperm_copy(&dimperm, a->datadim + change->intdim, perm);
    //ap_dimperm_clear (perm);
    //dimperm.dim[a->h->size] = n;
    for (i = 0; i < pr->size_scons; i++) {
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_dimperm_dimchange_n: with dimchange_n=(\n");
        ap_dimchange_fprint(stdout, change); // &dimchange);
        fprintf(stdout, ") on scons=(");
        ap_abstract0_fprint(stdout, pr->man_scons[i], a->scons[i], NULL);
        fprintf(stdout, ")\n");
#endif
        b->scons[i] =
                ap_abstract0_add_dimensions(pr->man_scons[i],
                false, a->scons[i], change, false); // &dimchange, false);

#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_dimperm_dimchange_n: after dimchange_n scons=(\n");
        ap_abstract0_fprint(stdout, pr->man_scons[i], b->scons[i], NULL);
        fprintf(stdout, ") apply dimperm=(");
        ap_dimperm_fprint(stdout, &dimperm);
        fprintf(stdout, ")\n");
#endif
        if (!shape_dimperm_is_id(&dimperm)) {
            b->scons[i] =
                    ap_abstract0_permute_dimensions(pr->man_scons[i],
                    true, b->scons[i], &dimperm);
#ifndef NDEBUG
            fprintf(stdout, "\n++++ushape_dimperm_dimchange_n: after dimperm_n scons=(");
            ap_abstract0_fprint(stdout, pr->man_scons[i], b->scons[i], NULL);
            fprintf(stdout, ")\n");
#endif
        }
    }
    ap_dimperm_clear(&dimperm);
    //ap_dimchange_clear(&dimchange);
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_dimperm_dimchange_n returns\n");
    fflush(stdout);
#endif
}

/* ============================================================ */
/* Expansion and folding of dimensions */
/* ============================================================ */

/* TODO: priority 0 */

/* not used: expand a dimension is done by add_dimension with project */
ushape_t *
ushape_expand(ap_manager_t * man,
        bool destructive, ushape_t * a, ap_dim_t dim, size_t n) {
    ushape_internal_t *pr = ushape_init_from_manager(man, AP_FUNID_EXPAND, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}

void
ushape_expand_internal(ushape_internal_t * pr, ushape_t * a, ushape_t * b,
        ap_dim_t dim, size_t n) {
    if (!a || !b || !b->scons || !dim || !n)
        return;
    size_t i;
    for (i = 0; i < pr->size_scons; i++)
        b->scons[i] = ap_abstract0_expand(pr->man_scons[i], false, a->scons[i],
            dim, n);
}

ushape_t *
ushape_fold(ap_manager_t * man,
        bool destructive, ushape_t * a, ap_dim_t * tdim, size_t size) {
    if (!a)
        return NULL;
    ushape_internal_t *pr = ushape_init_from_manager(man, AP_FUNID_FOLD, 0);
    if (!size)
        return (destructive) ? a : ushape_copy_internal(pr, a);
    size_t i;
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, a->h);
    b->closed = hgraph_copy_internal(pr, a->closed);
    for (i = 0; i < pr->size_scons; i++)
        b->scons[i] =
            ap_abstract0_fold(pr->man_scons[i], false, a->scons[i], tdim, size);
    if (destructive)
        ushape_free_internal(pr, a);
    return b;
}

void
ushape_fold_internal(ushape_internal_t * pr,
        ushape_t * a, ushape_t * b, ap_dim_array2_t * anon) {
    if (!a || !b || !b->scons || !anon->size)
        return;
    size_t i, s, changed;
    ap_dimperm_t perm;
    ap_dim_array2_t fanon; // array for folded dimensions
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_fold_internal: with anon [");
    for (i = 0; i < anon->size; i++)
        fprintf(stdout, "%zu(%d),", anon->p[i].size, anon->p[i].p[0]);
    fprintf(stdout, "]\n");
#endif
    // prepare the array of folded dimension to pass the ap_abstract0_fold, i.e.,
    // (1) all dimensions shall be real xor integer
    // (2) all dimensions shall be sorted
    // GO
    fanon.size = anon->size;
    fanon.p = (ap_dim_array_t *) malloc(fanon.size * sizeof (ap_dim_array_t));
    // (1) add a->datadim to all anon in order to obtain only real dimensions
    for (s = 0; s < anon->size; s++) {
        fanon.p[s].size = anon->p[s].size;
        fanon.p[s].p =
                (ap_dim_t *) malloc(anon->p[s].size * sizeof (ap_dim_t));
        for (i = 0; i < anon->p[s].size; i++)
            fanon.p[s].p[i] = anon->p[s].p[i] + a->datadim;
    }
    // (2) sort each entry of fanon and generate a permutation of dimensions folded
    ap_dimperm_init(&perm, a->datadim + a->h->size);
    ap_dimperm_set_id(&perm);
    changed = 0;
    for (s = 0; s < fanon.size; s++) { // insert sort entry s
        for (i = 1; i < fanon.p[s].size; i++) {
            ap_dim_t m, n;
            size_t j = 0;
            while (j < i && fanon.p[s].p[j] <= fanon.p[s].p[i])
                j++;
            m = fanon.p[s].p[i];
            while (j < i) {
                n = fanon.p[s].p[j];
                fanon.p[s].p[j] = m;
                m = n;
                j++;
            }
            fanon.p[s].p[i] = m;
        }
        // generate permutation: note that anon->p[s] are pairwise disjoint
        for (i = 0; i < anon->p[s].size; i++)
            if (fanon.p[s].p[i] != (a->datadim + anon->p[s].p[i])) {
                perm.dim[a->datadim + anon->p[s].p[i]] = fanon.p[s].p[i];
                changed = 1;
            }
    }
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_fold_internal: permutation [");
    ap_dimperm_fprint(stdout, &perm);
    fprintf(stdout, "]\n");
#endif
    // apply permutation before folding
    if (changed) {
        for (i = 0; i < pr->size_scons; i++)
            b->scons[i] =
                ap_abstract0_permute_dimensions(pr->man_scons[i], (a == b),
                a->scons[i], &perm);
    }
    // apply folding
    for (i = 0; i < pr->size_scons; i++) {
        if (changed)
            b->scons[i] =
                ap_abstract0_fold(pr->man_scons[i], true,
                b->scons[i], (ap_dim_t *) fanon.p[0].p,
                fanon.p[0].size);
        else
            b->scons[i] =
                ap_abstract0_fold(pr->man_scons[i], (a == b),
                a->scons[i], (ap_dim_t *) fanon.p[0].p,
                fanon.p[0].size);
        for (s = 1; s < anon->size; s++)
            b->scons[i] = ap_abstract0_fold(pr->man_scons[i], true, b->scons[i],
                (ap_dim_t *) fanon.p[s].p,
                fanon.p[s].size);
    }
    // apply the reverse permutation to obtain the final result
    if (changed) {
        ap_dimperm_t nperm;
        ap_dimperm_init(&nperm, perm.size);
        ap_dimperm_invert(&nperm, &perm);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_fold_internal: after folding ");
        ushape_fdump(stdout, pr->man, b);
        fprintf(stdout, "\n++++ushape_fold_internal: reverse permutation [");
        ap_dimperm_fprint(stdout, &nperm);
        fprintf(stdout, "]\n");
#endif
        for (i = 0; i < pr->size_scons; i++)
            b->scons[i] =
                ap_abstract0_permute_dimensions(pr->man_scons[i], true,
                b->scons[i], &nperm); 
        ap_dimperm_clear(&nperm);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_fold_internal: after permute ");
        ushape_fdump(stdout, pr->man, b);
        fprintf(stdout, "\n");
#endif
    }
    ap_dimperm_clear(&perm);
    // free all anon
    hgraph_node_anon_array_clear(pr, &fanon);
    hgraph_node_anon_array_clear(pr, anon);
}
