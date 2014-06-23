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
#include "ap_passign0.h"
#include "ap_pcons0.h"
#include "shape_macros.h"
#include "ap_generic.h"


/* ============================================================ */
/* Assignement and Substitutions */

/* ============================================================ */

ushape_t *
ushape_assign_data(ushape_internal_t * pr, ushape_t * a, size_t lhs,
        ap_linexpr0_t * lexpr, int* offsets) {
    /* TODO: check that lhs is not ptr null */

    ap_dim_t d = (lhs < a->datadim) ? lhs : NODE_DATADIM(a->datadim,
            a->h->size,
            VAR2NODE(a->h,
            DIM2PTR
            (lhs,
            a->
            datadim)));
    size_t *v2n = hgraph_get_var2node(a->h);
    ap_linexpr0_t *exprn = shape_linexpr_of_node(pr, lexpr, v2n, a->h->size,
            a->datadim, a->ptrdim);
#ifndef NDEBUG2
    fprintf(stdout, "%d:=", d);
    ap_linexpr0_fprint(stdout, exprn, NULL);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    free(v2n);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, a->h);
    b->closed = hgraph_copy_internal(pr, a->closed);
    bool isbot = false;
    size_t i;
    for (i = 0; i < pr->size_scons; i++) {
        b->scons[i] =
                ap_abstract0_assign_linexpr(pr->man_scons[i], false, a->scons[i],
                d, exprn, NULL);
        isbot = ap_abstract0_is_bottom(pr->man_scons[i], b->scons[i]);
    }
    ap_linexpr0_free(exprn);
    if (isbot)
        ushape_free_internal(pr, b);
    return b;
}

ushape_t *
ushape_substitute_data(ushape_internal_t * pr, ushape_t * a, size_t lhs,
        ap_linexpr0_t * lexpr, int* offsets) {
    /* TODO: check that lhs is not ptr null */

    ap_dim_t d = (lhs < a->datadim) ? lhs : NODE_DATADIM(a->datadim,
            a->h->size,
            VAR2NODE(a->h,
            DIM2PTR
            (lhs,
            a->
            datadim)));
    size_t *v2n = hgraph_get_var2node(a->h);
    ap_linexpr0_t *exprn = shape_linexpr_of_node(pr, lexpr, v2n, a->h->size,
            a->datadim, a->ptrdim);
#ifndef NDEBUG2
    fprintf(stdout, "%d:=", d);
    ap_linexpr0_fprint(stdout, exprn, NULL);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    free(v2n);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, a->h);
    b->closed = hgraph_copy_internal(pr, a->closed);
    bool isbot = false;
    size_t i;
    for (i = 0; i < pr->size_scons; i++) {
        b->scons[i] =
                ap_abstract0_substitute_linexpr(pr->man_scons[i], false, a->scons[i],
                d, exprn, NULL);
        isbot = ap_abstract0_is_bottom(pr->man_scons[i], b->scons[i]);
    }
    ap_linexpr0_free(exprn);
    if (isbot)
        ushape_free_internal(pr, b);
    return b;
}

ushape_t *
ushape_assign_x_y(ushape_internal_t * pr, ushape_t * a, size_t vx, size_t vy) {
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    hgraph_t *h = hgraph_assign_x_y(pr, a->h, vx, vy, &perm);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_apply_dimperm(pr, a, b, &perm);
    return b;
}

ushape_array_t *
ushape_substitute_x_null(ushape_internal_t * pr, bool isfree, ushape_t * a,
        size_t vx) {

    ushape_array_t * array_r;
    size_t i, rsize;
    ap_dimperm_t perm;
    array_r = ushape_array_make(pr, 1);
    rsize = 0;

    if (DIM2NODE(a->h, vx) != NODE_NULL) {
        ERROR("Bad substitution: vx is not null pointer variable!",;
                );
    }

    if (!isfree) {
        // x was at node i (i == NODE_NULL is added at the end)
        for (i = 1; i < a->h->size; i++) {
#ifndef NDEBUG
            fprintf(stdout, "\n++++ushape_substitute_null: map to node ");
            fprintf(stdout, "%d)\n", i);
            fflush(stdout);
#endif
            ap_dimperm_init(&perm, a->h->size);
            ap_dimperm_set_id(&perm);
            hgraph_t *h = hgraph_copy_mem(pr, a->h);
            hgraph_node_set_var(h, i, vx, &perm);
            ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
            b->h = hgraph_copy_internal(pr, h);
            hgraph_free_internal(pr, h);

            ushape_apply_dimperm(pr, a, b, &perm); // clear perm

            rsize += ushape_array_add(pr, true, array_r, 0, false, false, b);
        }
    }
    //x became garbage and its successor is any node from the graph
    for (i = 0; i < a->h->size; i++) {
        ap_dimperm_init(&perm, a->h->size);
        ap_dimperm_set_id(&perm);
        size_t nx; // new position of the node added
        hgraph_t *h = hgraph_node_add(pr, a->h, i, vx, 0, &nx, &perm);
        ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
        b->h = hgraph_copy_internal(pr, h);
        hgraph_free_internal(pr, h);

        ushape_apply_dimperm_dimchange(pr, a, b, &perm, nx); // nx is added

        rsize += ushape_array_add(pr, true, array_r, 0, false, false, b);
    }

    // x was pointing to a node that
    // (1) is now inside a word attached to some node != # (not for free)
    // (2) was garbage and pointing into the word as in (1)
    // The two cases are put together to avoid too much tests for exapansion.
    for (i = 1; i < a->h->size; i++) {
        // try to see if the node can be created in the middle of the word of i
        ushape_t *bex = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
        bex->h = hgraph_copy_internal(pr, a->h);
        ushape_expand_internal(pr, a, bex, i, 3); // does the split of i of a only for data
        // the node added has dimension a->size
        if (!ushape_is_bottom(pr->man, bex)) {
            // the expansion of i can be done!!
            // graph of bex is not more needed
            hgraph_free_internal(pr, bex->h);

            hgraph_t *h;
            ushape_t *b;

            // Case (1): the new node is labeled by x
            // add the new node in the graph
            if (!isfree) {
                ap_dimperm_init(&perm, a->h->size + 1);
                ap_dimperm_set_id(&perm);
                h = hgraph_node_expand(pr, a->h, i, vx, 0, &perm);
                b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
                b->h = hgraph_copy_internal(pr, h); // of size a->size + 1
                bex->h = h; // because data is already over a->size + 1
                ushape_apply_dimperm(pr, bex, b, &perm);
                // add to result
                rsize += ushape_array_add(pr, true, array_r, 0, false, false, b);
                hgraph_free_internal(pr, h);
            }

            // Case (2) two new nodes are added: the one labeled by x and
            // which points to the node splitting the word of i
            // do expansion
            ap_dimperm_init(&perm, a->h->size + 1);
            ap_dimperm_set_id(&perm);
            h = hgraph_node_expand(pr, a->h, i, NODE_VAR(a->h, i), NODE_VAR_NEXT(a->h, i) + 1, &perm);
            b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
            b->h = h; // of size a->size + 1
            bex->h = h; // because data is already over a->size + 1
            size_t nex = (size_t) perm.dim[a->h->size];
            ushape_apply_dimperm(pr, bex, b, &perm); // clear perm!

            // add the new node labeled by x which successor is the expanded node
            ap_dimperm_init(&perm, a->h->size + 1);
            ap_dimperm_set_id(&perm);
            size_t nx;
            hgraph_t *hr = hgraph_node_add(pr, h, nex, vx, 0, &nx, &perm);
            ushape_t *br = ushape_alloc_internal(pr, b->datadim, b->ptrdim);
            br->h = hgraph_copy_internal(pr, hr);
            ushape_apply_dimperm_dimchange(pr, b, br, &perm, nx);
            // add to result
            rsize += ushape_array_add(pr, true, array_r, 0, false, false, br);
            ushape_free_internal(pr, b);
            bex->h = NULL; // freed above
            hgraph_free_internal(pr, hr);
        }
        ushape_free_internal(pr, bex);
    }

    // x was already null
    if (!isfree)
        rsize += ushape_array_add(pr, true, array_r, 0, false, false, a);

    return array_r;
}

ushape_t *
ushape_substitute_x_y(ushape_internal_t * pr, ushape_t * a, size_t vx, size_t vy) {
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);

    ap_dim_array2_t anon;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    anon.size = 0;
    anon.p = NULL;

    hgraph_t *h = hgraph_substitute_x_y(pr, a->h, vx, vy, &perm, &anon);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_apply_dimperm(pr, a, b, &perm);
    return b;
}

ushape_t *
ushape_assign_x_next_y(ushape_internal_t * pr, ushape_t * a, size_t vx,
        size_t vy) {
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    hgraph_t *h = hgraph_assign_x_next_y(pr, a->h, vx, vy, &perm);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_apply_dimperm(pr, a, b, &perm);
    return b;
}

ushape_t *
ushape_substitute_x_next_y(ushape_internal_t * pr, ushape_t * a, size_t vx,
        size_t vy) {
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    hgraph_t *h = hgraph_substitute_x_next_y(pr, a->h, vx, vy, &perm);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_apply_dimperm(pr, a, b, &perm);
    return b;
}

ushape_t *
ushape_assign_x_y_next(ushape_internal_t * pr, bool conservative,
        ushape_t * a, size_t vx, size_t vy) {
    ap_dimperm_t perm;
    if (conservative)
        ap_dimperm_init(&perm, a->h->size);
    else
        ap_dimperm_init(&perm, a->h->size + 1);
    ap_dimperm_set_id(&perm);
    hgraph_t *h =
            hgraph_assign_x_y_next(pr, conservative, a->h, vx, vy, &perm);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_expand_internal(pr, a, b, (ap_dim_t) VAR2NODE(a->h, vy),
            (conservative) ? 1 : 2);
    ushape_apply_dimperm(pr, b, b, &perm);
    return b;
}

ushape_t *
ushape_substitute_x_y_next(ushape_internal_t * pr, bool conservative,
        ushape_t * a, size_t vx, size_t vy) {
    ap_dimperm_t perm;
    ap_dim_array2_t anon;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    anon.size = 0;
    anon.p = NULL;
    hgraph_t *h =
            hgraph_substitute_x_y_next(pr, conservative, a->h, vx, vy, &perm, &anon);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_expand_internal(pr, a, b, (ap_dim_t) VAR2NODE(a->h, vy),
            (conservative) ? 1 : 2);
    ushape_apply_dimperm(pr, b, b, &perm);
    return b;
}

/**
 * Vx shall be at NULL, vy is not NULL and vi is an integer variable
 */
ushape_t*
ushape_assign_x_y_at_i (ushape_internal_t * pr,
			ushape_t * a,
			size_t vx, size_t vy, size_t vi) 
{
  // 
  size_t px = DIM2PTR(vx,a->datadim);
  size_t py = DIM2PTR(vy,a->datadim);
  size_t nx = DIM2NODE(a->h, px);
  size_t ny = DIM2NODE(a->h, py);
  // assert (nx==0)
  // try to see if the node can be created 
  // in the middle of the word of vy at distance vi
  ushape_t *bex = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
  bex->h = hgraph_copy_internal(pr, a->h);
  ushape_expand_internal(pr, a, bex, ny, 4+vi); // does the split of vy of a only for data
  // the node added has dimension a->size
  if (!ushape_is_bottom(pr->man, bex)) {
    // the expansion of i can be done!!
    // graph of bex is not more needed
    hgraph_free_internal(pr, bex->h);
    
    hgraph_t *h;
    ushape_t *b;
    ap_dimperm_t perm;

    // Case (1): the new node is labeled by vx
    // add the new node in the graph
    ap_dimperm_init(&perm, a->h->size + 1);
    ap_dimperm_set_id(&perm);
    h = hgraph_node_expand(pr, a->h, ny, px, 0, &perm);
    b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h); // of size a->size + 1
    bex->h = h; // because data is already over a->size + 1
    ushape_apply_dimperm(pr, bex, b, &perm);
    return b;
  }
  return bex;
}

ushape_t *
ushape_assign_x_null(ushape_internal_t * pr, bool isfree, ushape_t * a,
        size_t vx) {
    ap_dimperm_t perm;
    ap_dim_array2_t anon;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    anon.size = 0;
    anon.p = NULL;
    hgraph_t *h;
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_assign_null: with a=(");
    ushape_fdump(stdout, pr->man, a);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif
    if (isfree)
        h = hgraph_assign_x_free(pr, a->h, vx, &perm, &anon);
    else
        h = hgraph_assign_x_null(pr, a->h, vx, &perm, &anon);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_assign_null: apply perm=[");
    ap_dimperm_fprint(stdout, &perm);
    fprintf(stdout, "] and anon=[");
    size_t i;
    for (i = 0; i < anon.size; i++)
        fprintf(stdout, "%zu(%d),", anon.p[i].size, anon.p[i].p[0]);
    fprintf(stdout, "]\n");
    fflush(stdout);
#endif
    ushape_fold_internal(pr, a, b, &anon);
#ifndef NDEBUG2
    fprintf(stdout, "\n++++ushape_assign_null: after fold: ");
    ushape_fdump(stdout, pr->man, b);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    ushape_apply_dimperm(pr, (anon.size == 0) ? a : b, b, &perm);
#ifndef NDEBUG
    fprintf(stdout, "\n++++ushape_assign_null returns: after dimperm: ");
    ushape_fdump(stdout, pr->man, b);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return b;
}

ushape_t *
ushape_assign_x_next_null(ushape_internal_t * pr, bool isfree, ushape_t * a,
        size_t vx) {
    ap_dimperm_t perm;
    ap_dim_array2_t anon;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    anon.size = 0;
    anon.p = NULL;
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    hgraph_t *h;
    if (isfree)
        h = hgraph_assign_x_next_free(pr, a->h, vx, &perm, &anon);
    else
        h = hgraph_assign_x_next_null(pr, a->h, vx, &perm, &anon);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    ushape_expand_internal(pr, a, b, (ap_dim_t) VAR2NODE(a->h, vx), 1);
    ushape_fold_internal(pr, b, b, &anon);
    ushape_apply_dimperm(pr, b, b, &perm);
    return b;
}

ushape_array_t *
ushape_substitute_x_next_null(ushape_internal_t * pr, bool isfree, ushape_t * a,
        size_t vx) {

    ushape_array_t * array_r;
    size_t i, rsize, size_data, size_ptr;
    ushape_t * a_copy;
    ap_dimperm_t perm, perm1;

    array_r = ushape_array_make(pr, 1);
    rsize = 0;
    size_data = 0;
    size_ptr = 0;
    size_t vy;

    if (NODE_NEXT(a->h, VAR2NODE(a->h, vx)) != NODE_NULL)
        ERROR("Bad substitution: x -> next should be null!",;
            );

    // x->next was any node from the graph
    for (i = 0; i < a->h->size; i++) {
        a_copy = ushape_copy_internal(pr, a);

#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_substitute_x_next_null: map to node ");
        fprintf(stdout, "%d)\n", i);
        fflush(stdout);
#endif
        ap_dimperm_init(&perm, a_copy->h->size);
        ap_dimperm_set_id(&perm);
        hgraph_t *h = hgraph_copy_mem(pr, a_copy->h);
        NODE_NEXT(h, VAR2NODE(h, vx)) = i;
        ushape_t *b = ushape_alloc_internal(pr, a_copy->datadim, a_copy->ptrdim);
        b->h = hgraph_copy_internal(pr, h);
        hgraph_free_internal(pr, h);

        ushape_apply_dimperm(pr, a_copy, b, &perm); // clear perm

        rsize += ushape_array_add(pr, true, array_r, 0, false, false, b);

    }

    //x->next became garbage and its successor is any node from the graph
    for (i = 0; i < a->h->size; i++) {

        a_copy = ushape_copy_internal(pr, a);

        ap_dimperm_init(&perm, a_copy->h->size);
        ap_dimperm_set_id(&perm);
        size_t nx_succ, vx_succ; // new position of the node added
        hgraph_t *h = hgraph_node_add(pr, a_copy->h, i, vx, 1, &nx_succ, &perm);
        ushape_t *b = ushape_alloc_internal(pr, a_copy->datadim, a_copy->ptrdim);
        b->h = hgraph_copy_internal(pr, h);
        hgraph_free_internal(pr, h);

        ushape_apply_dimperm_dimchange(pr, a_copy, b, &perm, nx_succ); // nx_succ is added

        ap_dimperm_init(&perm, a_copy->h->size + 1);
        ap_dimperm_set_id(&perm);
        hgraph_t *h1 = hgraph_copy_mem(pr, b->h);
        NODE_NEXT(h1, VAR2NODE(h1, vx)) = nx_succ;
        ushape_t *b1 = ushape_alloc_internal(pr, b->datadim, b->ptrdim);
        b1->h = hgraph_copy_mem(pr, h1);
        hgraph_free_internal(pr, h1);

        ushape_apply_dimperm(pr, a_copy, b1, &perm); // clear perm


        rsize += ushape_array_add(pr, true, array_r, 0, false, false, b1);
    }

    //x->next was pointing to a node that
    // (1) is now inside a word attached to some node != #
    // (2) was garbage and pointing into the word as in (1)
    // The two cases are put together to avoid too much tests for expansion.
    for (i = 1; i < a->h->size; i++) {

        a_copy = ushape_copy_internal(pr, a);

        // try to see if the node can be created in the middle of the word of i
        ushape_t *bex = ushape_alloc_internal(pr, a_copy->datadim, a_copy->ptrdim);
        bex->h = hgraph_copy_internal(pr, a_copy->h);
        ushape_expand_internal(pr, a_copy, bex, i, 3); // does the split of i of a only for data
        // the node added has dimension a->size
        if (!ushape_is_bottom(pr->man, bex)) {
            // the expansion of i can be done!!
            // graph of bex is not more needed
            hgraph_free_internal(pr, bex->h);

            // Case (1): the new node is labeled by x->next
            // add the new node in the graph
            ap_dimperm_init(&perm, a_copy->h->size + 1);
            ap_dimperm_set_id(&perm);
            hgraph_t *h = hgraph_node_expand(pr, a_copy->h, i, vx, 1, &perm);
            ushape_t *b = ushape_alloc_internal(pr, a_copy->datadim, a_copy->ptrdim);
            b->h = hgraph_copy_internal(pr, h); // of size a->size + 1
            bex->h = h; // because data is already over a->size + 1
            size_t nx_succ;
            nx_succ = perm.dim[a->h->size];
            ushape_apply_dimperm(pr, bex, b, &perm);

            ap_dimperm_init(&perm, a_copy->h->size + 1);
            ap_dimperm_set_id(&perm);
            hgraph_t *h11 = hgraph_copy_mem(pr, b->h);
            NODE_NEXT(h11, VAR2NODE(h11, vx)) = nx_succ;
            ushape_t *b11 = ushape_alloc_internal(pr, b->datadim, b->ptrdim);
            b11->h = hgraph_copy_internal(pr, h11);
            hgraph_free_internal(pr, h11);

            ushape_apply_dimperm(pr, a_copy, b11, &perm);

            // add to result
            rsize += ushape_array_add(pr, true, array_r, 0, false, false, b11);
            hgraph_free_internal(pr, h);

            // Case (2) two new nodes are added: the one labeled by x->next and
            // which points to the node splitting the word of i
            // do expansion
            ap_dimperm_init(&perm, a_copy->h->size + 1);
            ap_dimperm_set_id(&perm);
            h = hgraph_node_expand(pr, a_copy->h, i, NODE_VAR(a_copy->h, i), NODE_VAR_NEXT(a_copy->h, i) + 1, &perm);
            b = ushape_alloc_internal(pr, a_copy->datadim, a_copy->ptrdim);
            b->h = h; // of size a->size + 1
            bex->h = h; // because data is already over a->size + 1
            size_t nex = (size_t) perm.dim[a_copy->h->size];
            ushape_apply_dimperm(pr, bex, b, &perm); // clear perm!

            // add the new node labeled by x->next which successor is the expanded node
            ap_dimperm_init(&perm, a_copy->h->size + 1);
            ap_dimperm_set_id(&perm);
            size_t nx;
            hgraph_t *hr = hgraph_node_add(pr, h, nex, vx, 1, &nx, &perm);
            ushape_t *br = ushape_alloc_internal(pr, b->datadim, b->ptrdim);
            br->h = hgraph_copy_internal(pr, hr);
            ushape_apply_dimperm_dimchange(pr, b, br, &perm, nx);

            ap_dimperm_init(&perm, a_copy->h->size + 1);
            ap_dimperm_set_id(&perm);
            hgraph_t *h1 = hgraph_copy_mem(pr, br->h);
            NODE_NEXT(h1, VAR2NODE(h1, vx)) = nx;
            ushape_t *b1 = ushape_alloc_internal(pr, br->datadim, br->ptrdim);
            b1->h = hgraph_copy_mem(pr, h1);
            hgraph_free_internal(pr, h1);

            ushape_apply_dimperm(pr, a_copy, b1, &perm);

            // add to result
            rsize += ushape_array_add(pr, true, array_r, 0, false, false, b1);
            ushape_free_internal(pr, b);
            bex->h = NULL; // freed above
            hgraph_free_internal(pr, hr);
        }
        ushape_free_internal(pr, bex);
    }

    return array_r;
}

ushape_t *
ushape_assign_x_new(ushape_internal_t * pr, ushape_t * a, size_t vx) {
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    hgraph_t *h = hgraph_assign_x_new(pr, a->h, vx, &perm);
    ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
    b->h = hgraph_copy_internal(pr, h);
    hgraph_free_internal(pr, h);
    node_t nn = VAR2NODE(b->h, vx);
    ushape_apply_dimperm_dimchange(pr, a, b, &perm, nn);
    // done above: ap_dimperm_clear(&perm);
    // add constraint l[new-node] == 1
    ushape_expand_internal(pr, b, b, nn, 1);
    return b;
}

ushape_t *
ushape_substitute_x_new(ushape_internal_t * pr, ushape_t * a, size_t vx) {
    ap_dimperm_t perm;
    ap_dimperm_init(&perm, a->h->size);
    ap_dimperm_set_id(&perm);
    hgraph_t *h = hgraph_substitute_x_new(pr, a->h, vx, &perm);
    // no annonymous can be generated, but the node allocated is mapped to 0
    if (!hgraph_is_bottom(pr->man, h)) {
        ushape_t *b = ushape_alloc_internal(pr, a->datadim, a->ptrdim);
        b->h = hgraph_copy_internal(pr, h);
        hgraph_free_internal(pr, h);
        ushape_apply_dimperm(pr, a, b, &perm);
        return b;
    }
    return NULL;
}

ushape_array_t *
ushape_assign_passign(ushape_internal_t * pr, ushape_t * a, passign0_t * op) {
    ushape_array_t *r;
    ushape_t *b;
    size_t rsize;

    if (ushape_is_bottom(pr->man, a))
      return NULL;
    arg_assert(a && op, return NULL;
	       );
#ifndef NDEBUG1
    fprintf(stdout, "\n====ushape_assign_passign: on a=(");
    ushape_fdump(stdout, pr->man, a);
    fprintf(stdout, ")\n with passign=(");
    shape_passign_fdump(stdout, op, a->datadim, a->ptrdim);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    r = ushape_array_make(pr, 1);
    rsize = 0;
    // TODO: set all variables mapped to NODE_T_TOP in a->h to NULL node
    if (op->x < a->datadim || 
	op->op == PA_ASSIGN_INT) 
      { // data assignment, only a->scons affected
        b = ushape_assign_data(pr, a, op->x, 
			       op->info.data.expr, op->info.data.offsets);
        rsize += ushape_array_add(pr, true, r, rsize, true, true, b);
        return r;
      } 
    else 
      { // pure pointer assignment
        if (!IS_PTRDIM(op->lhs, a->datadim, a->ptrdim)) 
	  {
	    ERROR("Bad assignment: not a pointer left value!", return NULL;
		  );
	  }

        size_t vx = DIM2PTR(op->x, a->datadim);
        node_t nx = DIM2NODE(a->h, vx);

        if (nx == 0 && 
	    op->offx != OFFSET_NONE) 
	  {
            //ERROR("Bad assignment: null pointer dereference!", return NULL;
            //        );
            FEEDBACK("Bad assignment: null pointer dereference!", 
		     return NULL;
		     );
	  }

        if (op->offx == OFFSET_NONE &&
	    (op->op == PA_ALLOC || op->op == PA_ALLOC_N)) {
            // x := new
            b = ushape_assign_x_new(pr, a, vx);
            rsize += ushape_array_add(pr, true, r, 0, false, false, b);
            /* put b as built in r */
        } 
	else if (op->offx == OFFSET_NONE && 
		 op->op == PA_FREE) 
	  {
	    // x := free
	      b = ushape_assign_x_null(pr, true, a, vx);
	      rsize += ushape_array_add(pr, true, r, 0, false, false, b);
	  } 
	else if (op->offx == OFFSET_NONE && 
		 op->op == PA_ASSIGN_PTR) 
	  {
            // x := [null | y->next | y[i] ]
            if (!IS_NULLDIM(op->info.ptr.y) &&
		!IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim)) 
	      {
		ERROR("Bad assignment: not a pointer right value!",
		      return NULL;
		      );
	      }

            size_t vy = DIM2PTR(op->info.ptr.y, a->datadim);
            node_t ny = DIM2NODE(a->h, vy);

            if (ny == 0 && 
		op->info.ptr.offy != OFFSET_NONE) 
	      {
		// y==NULL and offset
		FEEDBACK("Bad assignment: null pointer dereference!", 
			 return NULL;
			 );
	      }

            if (//IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
		!IS_NULLDIM(op->info.ptr.y) &&
		ny != 0)
	      { // y is not null
		if (op->info.ptr.offy == OFFSET_NONE) 
		  { // x := y
		    b = ushape_assign_x_y(pr, a, vx, vy);
		    rsize += ushape_array_add(pr, true, r, 0, false, false, b); 
		    /* put b as built in r */
		  } 
		else if (op->info.ptr.offy == ((int) (a->datadim + OFFSET_NEXT))) 
		  { // x := y->next
		    // case (1): no new node
		    b = ushape_assign_x_y_next(pr, true, a, vx, vy);
		    rsize += ushape_array_add(pr, true, r, 0, false, false, b);
		    // case (2): new node
		    b = ushape_assign_x_y_next(pr, false, a, vx, vy);
		    rsize += ushape_array_add(pr, true, r, 0, false, false, b);
		  } 
		else if (IS_INTDIM(op->info.ptr.offy, a->datadim, a->ptrdim))
		  {
		    // x := y+i
		    // here some code specific to array encoding
		    if (nx != 0)
		      {
			// nothing to do, add a copy of a 
			rsize += ushape_array_add(pr, true, r, 0, 
						  true, false, a);
		      }
		    else
		      {
			// always a node created 
			b = ushape_assign_x_y_at_i(pr, a, op->x, 
						   op->info.ptr.y, op->info.ptr.offy);
			rsize += ushape_array_add(pr, true, r, 0, 
						  true, false, b);
		      }
		  }
		else 
		  {
		    ERROR("Bad assignment: not an allowed right value!",
			  return NULL;
			  );
		  }
	      }
	    else if (IS_NULLDIM(op->info.ptr.y) ||
		     ny == 0) 
	      { 
		// x := null | y (with y==null)
		if (nx != 0) 
		  {
		    b = ushape_assign_x_null(pr, false, a, vx);
		    rsize += ushape_array_add(pr, true, r, 0, false, false, b);
		  } 
		else 
		  { // else nothing to do, put b (copy of a) in r
		    b = ushape_copy_internal(pr, a);
		    rsize += ushape_array_add(pr, true, r, 0, false, false, b); /* put b as built in r */
		  }
	      } 
	  } 
	else 
	  { 
	    // op->type == NEXT_MODIFIER && IS_PTRDIM(op->lhs,a->datadim,a->ptrdim)
	    // x->next = [null | y]
	    if (op->offx != OFFSET_NONE && 
		op->info.ptr.offy != OFFSET_NONE) 
	      {
		ERROR("Bad assignment: not an elementary assignment!",
		      ;
		      );
	      }
            size_t vy = DIM2PTR(op->info.ptr.y, a->datadim);
            node_t ny = DIM2NODE(a->h, vy);
            if (IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
		!IS_NULLDIM(op->info.ptr.y) &&
		ny != 0) 
	      { // x->next := y, x->next is null
		b = ushape_assign_x_next_y(pr, a, vx, vy);
                rsize += ushape_array_add(pr, true, r, 0, false, false, b); /* put b as built in r */
	      } 
	    else if (IS_NULLDIM(op->info.ptr.y) ||
		     ny == 0) 
	      { // x->next = null
                b = ushape_assign_x_next_null(pr, false, a, vx);
                rsize += ushape_array_add(pr, true, r, 0, false, false, b);
	      } 
	    else 
	      {
	      ERROR("Bad assignment: not an elementary assignment!", 
		    ;
		    );
	      }
	  }
      }
    if (rsize)
      ushape_array_resize(pr, r, rsize);
    else if (r) 
      {
        ushape_array_clear(pr, r, r->size);
        free(r->p);
        free(r);
        r = NULL;
      }
    return r;
}

ushape_array_t *
ushape_assign_passign_array(ushape_internal_t * pr, bool destructive,
        ushape_t * a, passign0_array_t * arr) {
    ushape_array_t *r = ushape_array_make(pr, 1);
    r->p[0] = ushape_copy_internal(pr, a);
    size_t i;
    for (i = 0; i < arr->size; i++) { // apply each operation of the list on all elements of r
        ushape_array_t *rr = ushape_array_make(pr, r->size);
        size_t j;
        for (j = 0; j < r->size; j++)
            if (r->p[j]) {
                ushape_array_t *rj =
                        ushape_assign_passign(pr, r->p[j], arr->p[i]);
                if (rj) {
                    ushape_array_add_array(pr, true, rr, rj);
                    ushape_array_init(pr, rj, rj->size);
                    free(rj);
                }
            }
        ushape_array_init(pr, r, r->size);
        r = rr;
    }
    if (destructive)
        ushape_free_internal(pr, a);
    return r;
}

ushape_t *
ushape_assign_linexpr(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dim_t d, ap_linexpr0_t * expr, ushape_t * dest) {
    if (ushape_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : ushape_copy(man, a);
    else {
        ushape_t *b;
        passign0_array_t *op;
        ushape_array_t *r;
        ushape_internal_t *pr =
                ushape_init_from_manager(man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
        ap_dimension_t dim = ushape_dimension(man, a);
        if (!destructive)
            b = ushape_copy_internal(pr, a);
        else
            b = a;
        op =
                shape_passign_of_linexpr_array(pr, &d, &expr, 1, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_assign_linexpr: with passign=(");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
        fprintf(stdout, ")\n");
#endif
        /* go */
        r = ushape_assign_passign_array(pr, false, b, op);
        ushape_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            ushape_t *rr = ushape_meet(pr->man, false, b, dest);
            ushape_free_internal(pr, b);
            return rr;
        }
        return b;
    }
}

/**
 * Strengthen the data constraints of @code{a}.
 * At the present time, @code{ucons} is strengthened only
 * for patterns concerning the dimensions @code{tdim}.
 */
ushape_t*
ushape_strengthen(ushape_internal_t* pr,
        ushape_t* a,
        ap_dim_t* tdim, size_t size) {
    ushape_t* r = ushape_copy_internal(pr, a);
    // UCONS is not present, no strengthening
    if (pr->man_ucons >= pr->size_scons ||
            pr->man_mset >= pr->size_scons) /* TODO: may be improved */
        return r;

#ifndef NDEBUG2
    fprintf(stdout, "\n++++ushape_strengthen: tdim=[");
    for (size_t i = 0; i < size; i++)
        fprintf(stdout, "%zu, \n", tdim[i]);
    fprintf(stdout, "] \n");
    ushape_fdump(stdout, pr->man, a);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    // UCONS is present, strengthen
    // Step 1: build mset constraints tdim = 0 to signal tdim
    size_t eqsize = size + a->h->size;
    size_t csize = 0; // number of NOT NULL entries in eqsize
    ap_lincons0_array_t eqcons = ap_lincons0_array_make(eqsize);
    size_t i, k; // index of dimensions in tdim
    size_t n; // nodes
    for (i = 0; i < size; i++) {
        // build constraint node(tdim[i]) =(mset) 0 for nodes not NULL
        n = VAR2NODE(a->h, tdim[i] - a->datadim);
        if (n != NODE_NULL) {
            ap_scalar_t * zero = ap_scalar_alloc();
            ap_scalar_set_int(zero, OFFSET_MSET);
            ap_linexpr0_t* lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + a->h->size);
            ap_linexpr0_set_coeff_scalar_int(lexpr, a->datadim + n, 1);
            eqcons.p[csize] = ap_lincons0_make(AP_CONS_EQ, lexpr, zero);
            csize++;
        }
    }

    ap_lincons0_array_resize(&eqcons, csize);
    eqsize = csize;
#ifndef NDEBUG1
    fprintf(stdout, "\n++++ushape_strengthen: ?RRRbuilt constraints eqcons=(");
    fflush(stdout);
    if (eqcons.p) ap_lincons0_array_fprint(stdout, &eqcons, NULL);
    else
        fprintf(stdout, "\nNULL");
    fprintf(stdout, ") of size %zu\n", csize);
    fflush(stdout);
#endif
    if (pr->man_mset < pr->size_scons) {
        // Step 2: build mset constraints tdim - n = 0 to signal mset equality
        // between tdim and node x (different from tdim)
        // build the list of nodes which are not labeled in tdim
        bool *isin_a = (bool *) malloc(a->h->size * sizeof (bool));
        memset(isin_a, 0, a->h->size * sizeof (bool));
        for (i = 0; i < size; i++) {
            n = VAR2NODE(a->h, tdim[i] - a->datadim);
            isin_a[n] = true;
        }
        // iterate over each node not in a
        for (n = 1; n < a->h->size; n++)
            if (!isin_a[n]) {
                for (i = 0; i < size; i++) {
                    ap_linexpr0_t* lexpr = ap_linexpr0_alloc(AP_LINEXPR_DENSE,
                            a->datadim + a->h->size);
                    ap_dim_t ti = VAR2NODE(a->h, tdim[i] - a->datadim);
                    ap_linexpr0_set_coeff_scalar_int(lexpr, a->datadim + ti, 1);
                    ap_linexpr0_set_coeff_scalar_int(lexpr, a->datadim + n, -1);
                    ap_scalar_t * zero = ap_scalar_alloc();
                    ap_scalar_set_int(zero, OFFSET_MSET);
                    ap_lincons0_t lcons = ap_lincons0_make(AP_CONS_EQ, lexpr, zero);
#ifndef NDEBUG1
                    fprintf(stdout, "\n ti= %zu \n ++++ushape_strengthen: test constraints of type2 lcons=(", ti);
                    ap_lincons0_fprint(stdout, &lcons, NULL);
                    fprintf(stdout, ") \n");
                    fflush(stdout);
#endif

                    if (ap_abstract0_sat_lincons(pr->man_scons[pr->man_mset],
                            a->scons[pr->man_mset], &lcons)) {
#ifndef NDEBUG1
                        fprintf(stdout, "\n lcons is_sat a->h->size = %zu csize= %zu eqcons.size= %zu\n ", a->h->size, csize, eqcons.size);
                        fflush(stdout);
#endif
                        if (csize == eqsize) {
                            eqsize += a->h->size;
                            ap_lincons0_array_resize(&eqcons, eqsize);
                        }
                        /* send to other domains as nodes equalitie where tail has coefficient 2 */

                        eqcons.p[csize] = lcons;
                        csize++;

                    } else {
#ifndef NDEBUG1
                        fprintf(stdout, "\n lcons type 2 is_not_sat \n");
#endif
                    }


                    //#ifndef NDEBUG1
                    //					fprintf(stdout, "\n test for constraits of size 3");
                    //					fflush(stdout);
                    //#endif
                    //					// trasnform lcons, there might be constraints over 2 multisets as well as over 3 multisets
                    //					ap_lincons0_t lcons3 = ap_lincons0_copy(&lcons);
                    //
                    //					ap_linexpr0_realloc(lcons3.linexpr0,a->datadim + 2*a->h->size);
                    //
                    //					for (k = 1; k < a->h->size && !isin_a[k] && k!=n ; k++) {
                    //
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->datadim + a->h->size + ti, 1);
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->datadim + n, 1);
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->h->size + a->datadim + n, 1);
                    //
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->h->size + a->datadim + k, -1);
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->datadim + k, -1);
                    //
                    //#ifndef NDEBUG1
                    //						fprintf(stdout, "\n i=%zu tdim[i]=%zu n=%zu k=%zu ", i, a->datadim + a->h->size + ti, n, k);
                    //						fprintf(stdout, "\n++++ushape_strengthen: test constraints of type 3 lcons=(");
                    //						ap_lincons0_fprint(stdout, &lcons3, NULL);
                    //						fprintf(stdout, ") \n");
                    //						fflush(stdout);
                    //#endif
                    //
                    //						if (ap_abstract0_sat_lincons(pr->man_scons[pr->man_mset],
                    //								a->scons[pr->man_mset], &lcons3)){
                    //#ifndef NDEBUG1
                    //							fprintf(stdout, "\n lcons type 3 is_sat \n");
                    //							fflush(stdout);
                    //#endif
                    //							ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, tdim[i], 0);
                    //							ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->h->size + a->datadim + n, 0);
                    //							ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->h->size + a->datadim + k, 0);
                    //							ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->datadim + k, 0);
                    //
                    //							/* send to other domains as nodes equalitie where tail has coefficient 2 */
                    //
                    //							ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->datadim + k, -2);
                    //
                    //							if (csize == eqsize) {
                    //								eqsize += a->h->size;
                    //								ap_lincons0_array_resize(&eqcons, eqsize);
                    //							}
                    //
                    //							eqcons.p[csize] = ap_lincons0_copy(&lcons3);
                    //							csize++;
                    //						}
                    //						else{
                    //#ifndef NDEBUG1
                    //							fprintf(stdout, "\n lcons type 3 is_not_sat \n");
                    //							fflush(stdout);
                    //#endif
                    //						}
                    //
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->h->size + a->datadim + k, 0);
                    //						ap_linexpr0_set_coeff_scalar_int(lcons3.linexpr0, a->datadim + k, 0);
                    //
                    //					}//end check cts over 3 multisets
                    //					//ap_lincons0_clear(&lcons); // free also scalar and lexpr

                }
            }
    }
    ap_lincons0_array_resize(&eqcons, csize);
#ifndef NDEBUG1
    fprintf(stdout, "\n++++ushape_strengthen: built constraints eqcons=(");
    fflush(stdout);
    if (eqcons.p) ap_lincons0_array_fprint(stdout, &eqcons, NULL);
    else
        fprintf(stdout, "\nNULL");
    fprintf(stdout, ") of size %zu\n", csize);
    fflush(stdout);
#endif
    // Step 3: apply the constraint to the ucons
    if (r->scons[pr->man_ucons] != NULL)
        r->scons[pr->man_ucons] =
            ap_abstract0_meet_lincons_array(pr->man_scons[pr->man_ucons], true,
            r->scons[pr->man_ucons],
            &eqcons);
    else {
#ifndef NDEBUG1
        fprintf(stdout, "\n++++ushape_strengthen: empty object \n");
        fflush(stdout);
#endif
    }

    // free allocated resources
    ap_lincons0_array_clear(&eqcons);
    return r;
}

/**
 * Substitute actuals in tdim by the ini parameters in tnew.
 */
ushape_array_t*
ushape_substitute_actuals(ushape_internal_t* pr,
        ushape_t* a, ap_dim_t* tdim, ap_dim_t* tnew, size_t size) {
    ushape_array_t* r;
    ushape_t *b;
    hgraph_t* h;
    size_t i, rsize, size_ptr, size_data;
    ap_dim_t *tdim_ptr, *tdim_data;
    ap_dim_t *tnew_ptr, *tnew_data;
    ap_dimperm_t perm;

    if (ushape_is_bottom(pr->man, a))
        return NULL;
    arg_assert(a && tdim && tnew && size >= 1, return NULL;
            );
#ifndef NDEBUG1
    fprintf(stdout, "\n====ushape_substitute_actuals: (%d) with subst=(", pr->meet_algo);
    for (i = 0; i < size; i++)
        fprintf(stdout, "x%d/x%d, ", tnew[i], tdim[i]);
    fprintf(stdout, ") on a=(");
    ushape_fdump(stdout, pr->man, a);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    r = ushape_array_make(pr, 1);
    rsize = 0; // There is only one possible result for this operation
    size_data = 0;
    size_ptr = 0;
    // Step 0: generate a substitution for datadim, a substitution for ptrdim
    for (i = 0; i < size; i++)
        if (tdim[i] < a->datadim && tnew[i] < a->datadim)
            size_data++;
        else if (tdim[i] >= a->datadim && tnew[i] >= a->datadim)
            size_ptr++;
        else
            ERROR("bad dimension in substitution!",;);
    if (size_data >= 1) {
        checked_malloc(tdim_data, ap_dim_t, sizeof (ap_dim_t), size_data, return NULL;);
        checked_malloc(tnew_data, ap_dim_t, sizeof (ap_dim_t), size_data, return NULL;);
    } else {
        tdim_data = NULL;
        tnew_data = NULL;
    }

    if (size_ptr >= 1) {
        checked_malloc(tdim_ptr, ap_dim_t, sizeof (ap_dim_t), size_ptr, return NULL;);
        checked_malloc(tnew_ptr, ap_dim_t, sizeof (ap_dim_t), size_ptr, return NULL;);
    } else {
        tdim_ptr = NULL;
        tnew_ptr = NULL;
    }
    {
        size_t j, k;
        for (i = 0, j = 0, k = 0; i < size; i++)
            if (tdim[i] < a->datadim) {
                tdim_data[j] = tdim[i];
                tnew_data[j] = tnew[i];
                j++;
            } else {
                tdim_ptr[k] = tdim[i];
                tnew_ptr[k] = tnew[i];
                k++;
            }
    }
    // Substitution for ptr vars
    if (size_ptr >= 1) {
        // Step 1: strengthening
#ifndef NDEBUG1
        fprintf(stdout, "\n====ushape_substitute_actuals with ushape_strengthen");
        ushape_fdump(stdout, pr->man, a);
        fprintf(stdout, ")\n");
        fflush(stdout);
#endif
        b = ushape_strengthen(pr, a, tdim_ptr, size_ptr);

        // Step 2: substitute in the hgraph
        // generates permutation of nodes (with some garbage inside!)
        ap_dimperm_init(&perm, a->h->size);
        ap_dimperm_set_id(&perm);
        h = hgraph_substitute_array(pr, a->h, tdim_ptr, tnew_ptr, size_ptr, &perm);
        hgraph_free_internal(pr, b->h);
        b->h = hgraph_copy_internal(pr, h);
        hgraph_free_internal(pr, h);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_substitute_actuals: returns hgraph=[");
        hgraph_fdump(stdout, pr->man, b->h);
        fprintf(stdout, "] and perm=[");
        ap_dimperm_fprint(stdout, &perm);
        fprintf(stdout, "]\n");
        fflush(stdout);
#endif
        // Step 3: apply permutation (+ elimination) of nodes
        ushape_apply_dimperm(pr, b, b, &perm);
        // free perm
        ap_dimperm_clear(&perm);
    } else {
        b = ushape_copy_internal(pr, a);
    }

    // Step 4: apply substitution for the scalar variables
    if (size_data >= 1) {
        bool isbot = false;
        // build texpr from tnew_data
        ap_texpr0_t** texpr_data;
        checked_malloc(texpr_data, ap_texpr0_t*, sizeof (ap_texpr0_t*), size_data, return NULL;);
        for (i = 0; i < size_data; i++)
            texpr_data[i] = ap_texpr0_dim(tnew_data[i]);

        // apply substitution on each scons
        for (i = 0; i < pr->size_scons && !isbot; i++) {
            b->scons[i] =
                    ap_abstract0_substitute_texpr_array(pr->man_scons[i], true, b->scons[i],
                    tdim_data, texpr_data, size, NULL);
            isbot = ap_abstract0_is_bottom(pr->man_scons[i], b->scons[i]);
        }
        // free allocated data
        for (i = 0; i < size_data; i++)
            ap_texpr0_free(texpr_data[i]);
        free(texpr_data);
        if (isbot)
            ushape_free_internal(pr, b);
        b = NULL;
    }

    // put the new result
    rsize += ushape_array_add(pr, true, r, rsize, true, true, b);

    // Free allocated memory
    free(tdim_data);
    free(tdim_ptr);
    free(tnew_data);
    free(tnew_ptr);

    return r;
}

/**
 * Pre computation for assignements in passign.
 */
ushape_array_t *
ushape_substitute_passign(ushape_internal_t* pr,
        ushape_t * a, passign0_t * op) {
    ushape_array_t *r;
    ushape_t *b;
    ushape_array_t *b_array;

    if (ushape_is_bottom(pr->man, a))
        return NULL;
    arg_assert(a && op, return NULL;
            );

#ifndef NDEBUG1
    fprintf(stdout, "\n====ushape_substitute_passign: with passign=(");
    shape_passign_fdump(stdout, op, a->datadim, a->ptrdim);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    r = ushape_array_make(pr, 1);
    // TODO: set all variables mapped to NODE_T_TOP in a->h to NULL node
    if (op->x < a->datadim || op->op == PA_ASSIGN_INT) {
        // data assignment, only a->scons affected
        b = ushape_substitute_data(pr, a, op->x,
                op->info.data.expr, op->info.data.offsets);
        ushape_array_add(pr, true, r, 0, true, true, b);
        return r;
    } else { // pure pointer substitution
        if (!IS_PTRDIM(op->lhs, a->datadim, a->ptrdim)) {
            ERROR("Bad subsitution: not a pointer left value!", return NULL;
                    );
        }

        size_t vx = DIM2PTR(op->x, a->datadim);
        node_t nx = DIM2NODE(a->h, vx);

        if (op->offx == OFFSET_NONE) { // x := [new | free | null | y | y->next]
            switch (op->op) {
                case PA_ALLOC:
                case PA_ALLOC_N:
                { // x := new
                    b = ushape_substitute_x_new(pr, a, vx);
                    ushape_array_add(pr, true, r, 0, false, false, b); /* put b as built in r */
                    break;
                }

                case PA_FREE:
                { // x := free
                    if (nx == 0) {
                        b_array = ushape_substitute_x_null(pr, true, a, vx);
                        ushape_array_add_array(pr, true, r, b_array);
                    } else {
                        ERROR("Bad substitution: not an allowed right value!",
                        return NULL;);
                    }
                    break;
                }

                default:
                { // other assignments
                    if (!IS_NULLDIM(op->info.ptr.y) &&
                            !IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim)) {
                        ERROR("Bad assignment: not a pointer right value!", return NULL;
                                );
                    }

                    size_t vy = DIM2PTR(op->info.ptr.y, a->datadim);
                    node_t ny = DIM2NODE(a->h, vy);

                    if (ny == 0 && op->info.ptr.offy != OFFSET_NONE) {
                        //ERROR("Bad assignment: null pointer dereference!", return NULL;
                        //        );
                        FEEDBACK("Bad assignment: null pointer dereference!", return NULL;)
                    }

                    if (//IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
                            IS_NULLDIM(op->info.ptr.y)) { // x := null
                        if (nx == 0) {
                            b_array = ushape_substitute_x_null(pr, !IS_NULLDIM(op->info.ptr.y), a, vx);
                            ushape_array_add_array(pr, true, r, b_array);
                        } else {
                            ERROR("Bad substitution: not an allowed right value!",
                            return NULL;);
                        }
                    } else if (//IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
                            op->info.ptr.offy == OFFSET_NONE) { // x := y
                        b = ushape_substitute_x_y(pr, a, vx, vy);
                        ushape_array_add(pr, true, r, 0, false, false, b); /* put b as built in r */
                    } else if (// IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
                            //!IS_NULLDIM(op->info.ptr.y) &&
                            op->info.ptr.offy == ((int) a->datadim + OFFSET_NEXT)) {
                        // x := y->next
                        // case (1): no new node
                        b = ushape_substitute_x_y_next(pr, true, a, vx, vy);
                        ushape_array_add(pr, true, r, 0, false, false, b);
                        // case (2): new node
                        b = ushape_substitute_x_y_next(pr, false, a, vx, vy);
                        ushape_array_add(pr, true, r, 0, false, false, b);
                    } else {
                        ERROR("Bad substitution: not an allowed right value!",
                        return NULL;);
                    }
                    break;
                }
            } // end switch

        } else if (op->info.ptr.offy != ((int) a->datadim + OFFSET_NEXT)) {
            // x->next = [null | y]

            if (nx == 0) {
                //ERROR("Bad substitution: null pointer dereference!", return NULL;
                //        );
                FEEDBACK("Bad substitution: null pointer dereference!", return NULL;)
            }

            size_t vy = DIM2PTR(op->info.ptr.y, a->datadim);
            if (IS_NULLDIM(op->info.ptr.y) &&
                    op->info.ptr.offy == OFFSET_NONE) { // x->next = null
                b_array = ushape_substitute_x_next_null(pr, false, a, vx);
                ushape_array_add_array(pr, true, r, b_array);
            } else if (IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
                    !IS_NULLDIM(op->info.ptr.y) &&
                    op->info.ptr.offy == OFFSET_NONE) { // x->next := y,
                // x->next is null, y is not NULL
                b = ushape_substitute_x_next_y(pr, a, vx, vy);
                ushape_array_add(pr, true, r, 0, false, false, b);
                /* put b as built in r */
            } else {
                ERROR("Bad substitution: not an elementary substitution!",;
                        );
            }
        } else { // TODO: other pointer fields than next
            ERROR("Bad substitution: not an elementary substitution!",;
                    );
        }
    }
    if (r && r->size)
        ushape_array_resize(pr, r, r->size);
    else if (r) {
        ushape_array_clear(pr, r, r->size);
        free(r->p);
        free(r);
        r = NULL;
    }
#ifndef NDEBUG1
    fprintf(stdout, "\n====ushape_substitute_passign: returns %zu graphs=[",
            (r) ? r->size : 0);
    size_t i;
    for (i = 0; r && i < r->size; i++) {
        fprintf(stdout, "graph %zu = (", i);
        ushape_fdump(stdout, pr->man, r->p[i]);
        fprintf(stdout, ")\n");
    }
    fprintf(stdout, "]\n");
    fflush(stdout);
#endif

    return r;
}

ushape_array_t *
ushape_substitute_passign_array(ushape_internal_t* pr,
        ushape_t * a, passign0_array_t * arr) {
    ushape_array_t *r = ushape_array_make(pr, 1);
    r->p[0] = ushape_copy_internal(pr, a);
    size_t i;
    for (i = 0; i < arr->size; i++) { // apply each operation of the list on all elements of r
        ushape_array_t *rr = ushape_array_make(pr, r->size);
        size_t j;
        for (j = 0; j < r->size; j++)
            if (r->p[j]) {
                ushape_array_t *rj =
                        ushape_substitute_passign(pr, r->p[j], arr->p[i]);
                if (rj) {

                    ushape_array_add_array(pr, true, rr, rj);
                    ushape_array_init(pr, rj, rj->size);
                    free(rj);
                }
            }
        ushape_array_init(pr, r, r->size);
        r = rr;
    }
    return r;
}

/* used for pre-image computation */
ushape_t *
ushape_substitute_linexpr(ap_manager_t * man, /* TODO proprity 0 */
        bool destructive, ushape_t * a,
        ap_dim_t d, ap_linexpr0_t * expr, ushape_t * dest) {
    if (ushape_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : ushape_copy(man, a);
    else {
        ushape_t *b;
        passign0_array_t *op;
        ushape_array_t *r;
        ushape_internal_t *pr =
                ushape_init_from_manager(man, AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY, 0);
        ap_dimension_t dim = ushape_dimension(man, a);
        if (!destructive)
            b = ushape_copy_internal(pr, a);
        else
            b = a;
        op =
                shape_passign_of_linexpr_array(pr, &d, &expr, 1, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_assign_linexpr: with passign=(");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
        fprintf(stdout, ")\n");
#endif
        /* go */
        r = ushape_substitute_passign_array(pr, b, op);
        ushape_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            ushape_t *rr = ushape_meet(pr->man, false, b, dest);
            ushape_free_internal(pr, b);

            return rr;
        }
        return b;
    }
}

ushape_t *
ushape_assign_linexpr_array(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dim_t * tdim,
        ap_linexpr0_t ** texpr,
        size_t size, ushape_t * dest) {
    if (ushape_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : ushape_copy(man, a);
    else {
        ushape_t *b;
        passign0_array_t *op;
        ushape_array_t *r;
        ushape_internal_t *pr =
                ushape_init_from_manager(man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
        ap_dimension_t dim = ushape_dimension(man, a);
        if (!destructive)
            b = ushape_copy_internal(pr, a);
        else
            b = a;
        op =
                shape_passign_of_linexpr_array(pr, tdim, texpr, size, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_assign_linexpr_array: with passign_array=(");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
        fprintf(stdout, ")\n");
#endif
        /* go */
        r = ushape_assign_passign_array(pr, false, b, op);
        ushape_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            ushape_t *rr = ushape_meet(pr->man, false, b, dest);
            ushape_free_internal(pr, b);

            return rr;
        }
        return b;
    }
}

/* TODO: priority 0 */

/* used for pre-image computation */
ushape_t *
ushape_substitute_linexpr_array(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dim_t * tdim,
        ap_linexpr0_t ** texpr,
        size_t size, ushape_t * dest) {
    if (size == 1)
        return ushape_substitute_linexpr(man, destructive, a, tdim[0], texpr[0],
            dest);

    ushape_internal_t * pr =
            ushape_init_from_manager(man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");

    return a;
}

ushape_t *
ushape_assign_texpr_array(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dim_t * tdim,
        ap_texpr0_t ** texpr, size_t size, ushape_t * dest) {
    if (ushape_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : ushape_copy(man, a);
    else {
        ushape_t *b;
        passign0_array_t *op;
        ushape_array_t *r;
        ushape_internal_t *pr =
                shape_init_from_manager(man, AP_FUNID_ASSIGN_TEXPR_ARRAY, 0);
        ap_dimension_t dim = ushape_dimension(man, a);
        if (!destructive)
            b = ushape_copy_internal(pr, a);
        else
            b = a;
        op =
                shape_passign_of_texpr_array(pr, tdim, texpr, size, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_assign_texpr_array: with passign_array=(\n");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
        fprintf(stdout, ")\n");
#endif
        /* go */
        r = ushape_assign_passign_array(pr, false, b, op);
        ushape_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            ushape_t *rr = ushape_meet(pr->man, false, b, dest);
            ushape_free_internal(pr, b);

            return rr;
        }
        return b;
    }
}

/* TODO: priority 0 */

/* used only for pre-image computation */
ushape_t *
ushape_substitute_texpr_array(ap_manager_t * man,
        bool destructive, ushape_t * a,
        ap_dim_t * tdim,
        ap_texpr0_t ** texpr,
        size_t size, ushape_t * dest) {

    if (ushape_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : ushape_copy(man, a);
    else {
        ushape_t *b;
        passign0_array_t *op;
        ushape_array_t *r;
        ushape_internal_t *pr =
                shape_init_from_manager(man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY, 0);
        ap_dimension_t dim = ushape_dimension(man, a);
        if (!destructive)
            b = ushape_copy_internal(pr, a);
        else
            b = a;
        op =
                shape_passign_of_texpr_array(pr, tdim, texpr, size, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        fprintf(stdout, "\n++++ushape_substitute_texpr_array: with passign_array=(\n");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
        fprintf(stdout, ")\n");
#endif
        /* go */
        r = ushape_substitute_passign_array(pr, b, op);
        ushape_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            ushape_t *rr = ushape_meet(pr->man, false, b, dest);
            ushape_free_internal(pr, b);
            return rr;
        }
        return b;
    }
}
