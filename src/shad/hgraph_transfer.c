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
#include "apron2shape.h"
#include "shape_macros.h"
#include "ap_generic.h"



/* ============================================================ */
/* Assignement and Substitutions */
/* ============================================================ */

/* Assignment of special kind */

void
hgraph_node_anon_array_clear(hgraph_internal_t * pr, ap_dim_array2_t * anon) {
    if (anon && anon->size) {
        node_t i;
        for (i = 0; i < anon->size; i++) {
            free(anon->p[i].p);
            anon->p[i].p = NULL;
        }
        free(anon->p);
    }
}

/* Assignement vx := vy with vy not labelling null
 * No garbage, no anonymous nodes.
 * Permutation of old nodes.
 */
hgraph_t *
hgraph_assign_x_y(hgraph_internal_t * pr, hgraph_t * a,
        size_t vx, size_t vy, ap_dimperm_t * perm) {
    // x shall label null or top
    if (!(DIM2NODE(a, vx) == 0 || DIM2NODE(a, vx) == NODE_T_TOP)) {
        ERROR("Bad assignment: assign a not null pointer variable!",;
                );
    }
    // y shall be defined
    if (DIM2NODE(a, vy) == NODE_T_TOP) {
        ERROR("Bad assignment: assign from an undefined variable!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    hgraph_node_set_var(r, VAR2NODE(a, vy), vx, perm);
    return r;
}

/* Substitution vx := vy with vy not labelling null
 * No garbage, no anonymous nodes.
 * Permutation of old nodes.
 */
hgraph_t *
hgraph_substitute_x_y(hgraph_internal_t * pr, hgraph_t * a,
        size_t vx, size_t vy, ap_dimperm_t * perm, ap_dim_array2_t * anon) {

    if (DIM2NODE(a, vx) != DIM2NODE(a, vy)) {
        ERROR("Bad substitution: x should be equal to y!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    r = hgraph_assign_x_null(pr, a, vx, perm, anon);

    return r;
}

/* Assignment vx := vy->next
 * New node is created if not conservative.
 * No garbage, no anonymous nodes.
 * Permutation perm have size a->size if conservative, a->size+1 otherwise.
 */
hgraph_t *
hgraph_assign_x_y_next(hgraph_internal_t * pr, bool conservative,
        hgraph_t * a, size_t vx, size_t vy,
        ap_dimperm_t * perm) {
    node_t nx = DIM2NODE(a, vx);
    node_t ny = DIM2NODE(a, vy);
    node_t nny = NODE_NEXT(a, ny);
    hgraph_t *r = NULL;

    // x shall label null
    if (nx != 0) {
        ERROR("Bad assignment: assign a not null pointer variable!",;
                );
    }
    // y is not null or undefined
    if (ny == 0 || ny == NODE_T_TOP) {
        ERROR("Bad dereference: null pointer!",;
                );
    }
    if (nny == NODE_T_TOP) {
        ERROR("Bad assignment: assign from an undefined next pointer!",;
                );
    }
    if (conservative) { // case y->next is exactly x
        arg_assert(perm->size == a->size, return NULL;
                );
        r = hgraph_copy_mem(pr, a);
        hgraph_node_set_var(r, nny, vx, perm);
    } else { // case y->next is a new node
        node_t i;
        size_t j;
        arg_assert(perm->size == (a->size + 1), return NULL;
                );
        ap_dimperm_t perm1, perm2;
        ap_dimperm_init(&perm1, a->size);
        ap_dimperm_set_id(&perm1);
        r = hgraph_node_add(pr, a, nny, vx, 0, &i, &perm1);
        ap_dimperm_init(&perm2, r->size);
        ap_dimperm_set_id(&perm2);
        hgraph_node_set_succ(r, perm1.dim[ny], i, &perm2);
        for (j = 0; j < a->size; j++)
            perm->dim[j] = perm2.dim[perm1.dim[j]];
        perm->dim[a->size] = perm2.dim[i];
    }
    return r;
}

hgraph_t *
hgraph_substitute_x_y_next(hgraph_internal_t * pr, bool conservative,
        hgraph_t * a, size_t vx, size_t vy,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {

    if (DIM2NODE(a, vx) != NODE_NEXT(a, DIM2NODE(a, vy))) {
        ERROR("Bad substitution: x should be equal to y->next!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    r = hgraph_assign_x_null(pr, a, vx, perm, anon);

    return r;
}

/* Assignment vx := null
 * Garbage possible, anonymous possible.
 * Permutation has size a->size
 */
hgraph_t *
hgraph_assign_x_null(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {
    arg_assert(perm->size == a->size, return NULL;
            );
    hgraph_t *r = hgraph_copy_mem(pr, a);
    VAR2NODE(r, vx) = NODE_NULL;
    hgraph_t *rr = hgraph_close(pr, r, perm, anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Substitution vx := null
 * Garbage possible, anonymous possible.
 * Permutation has size a->size
 */
hgraph_t *
hgraph_substitute_x_null(hgraph_internal_t * pr, hgraph_t * a, size_t vx, size_t vy,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {
    arg_assert(perm->size == a->size, return NULL;
            );
    hgraph_t *r = hgraph_copy_mem(pr, a);
    VAR2NODE(r, vx) = vy;
    hgraph_t *rr = hgraph_close(pr, r, perm, anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Assignment vx := new
 * No garbage, no anonymous.
 * Permutation perm has size a->size, node a->size is the one added.
 */
hgraph_t *
hgraph_assign_x_new(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        ap_dimperm_t * perm) {
    arg_assert(perm && (perm->size == a->size), return NULL;
            );
    // x shall label null
    if (DIM2NODE(a, vx) != 0) {
        ERROR("Bad assignment: assign a not null pointer variable!",;
                );
    }
    node_t i;
    hgraph_t *r = hgraph_node_add(pr, a, 0, vx, 0, &i, perm);
    // assert (i == a->size);
    return r;
}

/* Substitution vx := new
 * No garbage, no anonymous, but change of dimension!
 * Permutation perm has size a->size.
 * At the end, the node DIM2NODE(a,vx) ismapped to 0.
 */
hgraph_t*
hgraph_substitute_x_new(hgraph_internal_t* pr, hgraph_t* a, size_t vx,
        ap_dimperm_t* perm) {

    arg_assert(perm && (perm->size == a->size), return NULL;
            );
    // vx shall be not null but isolated,
    // i.e., nobody points to and all ptr fields point to NULL
    node_t nx = DIM2NODE(a, vx);
    if (nx == NODE_NULL || NODE_NEXT(a, nx) != NODE_NULL) {
        return NULL; // do not generate graph for this case
    }
    node_t i;
    for (i = 0; i < a->ptrdim; i++)
        if (i != vx && DIM2NODE(a, i) == nx)
            return NULL; // no variable other than vx labels nx
    for (i = 0; i < a->size; i++)
        if (NODE_NEXT(a, i) == nx)
            return NULL; // no other node points to nx
    // put vx to NULL and recompute the graph
    hgraph_t *r = hgraph_copy_mem(pr, a);
    VAR2NODE(r, vx) = NODE_NULL;
    ap_dim_array2_t anon;
    anon.size = 0;
    anon.p = NULL;
    hgraph_t *rr = hgraph_close(pr, r, perm, &anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Assignment vx := free
 * Garbage possible, anonymous possible.
 * Permutation has size a->size.
 */
hgraph_t *
hgraph_assign_x_free(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {
    arg_assert(perm->size == a->size, return NULL;
            );
    // x is not null
    if (DIM2NODE(a, vx) == 0) {
        ERROR("Bad free: null pointer!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    node_t nx = DIM2NODE(r, vx);
    size_t i;
    for (i = 0; i < r->size; i++)
        if (NODE_NEXT(r, i) == nx)
            NODE_NEXT(r, i) = 0;
    VAR2NODE(r, vx) = 0;
    hgraph_t *rr = hgraph_close(pr, r, perm, anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Assignment vx->next := null
 * vx is a not null pointer.
 * Garbage possible, anonymous possible.
 * Permutation has size a->size.
 */
hgraph_t *
hgraph_assign_x_next_null(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {
    arg_assert(perm->size == a->size, return NULL;
            );
    // x is not null
    if (DIM2NODE(a, vx) == 0) {
        ERROR("Bad assignment: dereference of a null pointer!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    NODE_NEXT(r, VAR2NODE(r, vx)) = 0;
    hgraph_t *rr = hgraph_close(pr, r, perm, anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Substitute vx->next := null
 * vx is a not null pointer.
 * Garbage possible, anonymous possible.
 * Permutation has size a->size.
 */
hgraph_t *
hgraph_substitute_x_next_null(hgraph_internal_t * pr, hgraph_t * a, size_t vx, size_t vy,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {
    arg_assert(perm->size == a->size, return NULL;
            );
    // x is not null
    if (DIM2NODE(a, vx) == 0) {
        ERROR("Bad substitution: dereference of a null pointer!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    NODE_NEXT(r, VAR2NODE(r, vx)) = vy;
    hgraph_t *rr = hgraph_close(pr, r, perm, anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Assignment vx->next := free
 * Garbage possible, anonymous possible.
 * Permutation has size a->size.
 */
hgraph_t *
hgraph_assign_x_next_free(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        ap_dimperm_t * perm, ap_dim_array2_t * anon) {
    arg_assert(perm->size == a->size, return NULL;
            );
    // x is not null
    if (DIM2NODE(a, vx) == 0) {
        ERROR("Bad free: null pointer!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    node_t nnx = NODE_NEXT(r, DIM2NODE(r, vx));
    size_t i;
    for (i = 0; i < r->size; i++)
        if (NODE_NEXT(r, i) == nnx)
            NODE_NEXT(r, i) = 0;
    VAR2NODE(r, vx) = 0;
    hgraph_t *rr = hgraph_close(pr, r, perm, anon, false);
    hgraph_free_internal(pr, r);
    return rr;
}

/* Assignment vx->next := y
 * vx is a not null pointer, but vx->next yes!.
 * No garbage, no anonymous.
 * Permutation has size a->size.
 */
hgraph_t *
hgraph_assign_x_next_y(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        size_t vy, ap_dimperm_t * perm) {
    arg_assert(perm->size == a->size, return NULL;
            );
    // x is not null
    if (DIM2NODE(a, vx) == 0) {
        ERROR("Bad dereference: null pointer!",;
                );
    }
    // x->next is null
    if (NODE_NEXT(a, DIM2NODE(a, vx)) != 0) {
        ERROR("Bad assignment: not an elementary assignment!",;
                );
    }
    hgraph_t *r = hgraph_copy_mem(pr, a);
    hgraph_node_set_succ(r, DIM2NODE(a, vx), DIM2NODE(a, vy), perm);
    return r;
}

/* Assignment vx->next := y
 * vx is a not null pointer, but vx->next yes!.
 * No garbage, no anonymous.
 * Permutation has size a->size.
 */
hgraph_t *
hgraph_substitute_x_next_y(hgraph_internal_t * pr, hgraph_t * a, size_t vx,
        size_t vy, ap_dimperm_t * perm) {
    arg_assert(perm->size == a->size, return NULL;
            );

    // x -> next equals y
    if (NODE_NEXT(a, DIM2NODE(a, vx)) != DIM2NODE(a, vy)) {
        ERROR("Bad substitution: x->next and y should be the same!",;
                );
    }

    hgraph_t *r = hgraph_copy_mem(pr, a);
    hgraph_node_set_succ(r, DIM2NODE(a, vx), NODE_NULL, perm);

    return r;
}

/* Do the assignment op on the input hgraph a using the same
 * (if conservative) or not number of nodes.
 * The change in numbering nodes due to assignment is given by perm:
 * - old nodes which are garbage are mapped to 0
 * - annon nodes are mapped to the new node
 * - old nodes are mapped to new nodes
 * If more than pr->k anonymous nodes are created,
 * return in anon->size the number of such segments,
 *        in anon->orig (array of anon->size) old nodes head of the segment
 *        in anon->nodes (matrix of [anon->size][pr->k] old anonymous nodes ordered
 */
hgraph_t *
hgraph_assign_passign(hgraph_internal_t * pr, bool conservative,
        hgraph_t * a, passign0_t * op) {
    if (op->lhs < a->datadim || op->op == PA_ASSIGN_INT) { // data assignment, do nothing
        return hgraph_copy_internal(pr, a);
    }
    if (!IS_PTRDIM(op->lhs, a->datadim, a->ptrdim)) {
        ERROR("Bad assignment: not a pointer left value!",;
                );
    }
    if (!IS_NULLDIM(op->info.ptr.y) &&
            !IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim)) {
        ERROR("Bad assignment: not a pointer right value!",;
                );
    }

    hgraph_t *r = NULL;
    ap_dimperm_t perm;
    ap_dim_array2_t anon;
    size_t vx = DIM2PTR(op->x, a->datadim);
    node_t nx = DIM2NODE(a, vx);
    size_t vy = DIM2PTR(op->info.ptr.y, a->datadim);
    node_t ny = DIM2NODE(a, vy);

    if (op->offx == OFFSET_NONE) {
        // x := [null | new | free | y | y->next]
        if (IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
                !IS_NULLDIM(op->info.ptr.y) &&
                ny != 0 && ((int) op->info.ptr.offy) == OFFSET_NONE) {
            // x := y
            ap_dimperm_init(&perm, a->size);
            ap_dimperm_set_id(&perm);
            r = hgraph_assign_x_y(pr, a, vx, vy, &perm);
            ap_dimperm_clear(&perm);
        } else if (IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim) &&
                !IS_NULLDIM(op->info.ptr.y) &&
                op->info.ptr.offy == ((int) (a->datadim + OFFSET_NEXT))) {
            // x := y->next
            if (conservative) {
                // case (1): no new node
                ap_dimperm_init(&perm, a->size);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_y_next(pr, true, a, vx, vy, &perm);
                ap_dimperm_clear(&perm);
            } else {
                // case (2): new node
                ap_dimperm_init(&perm, a->size + 1);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_y_next(pr, false, a, vx, vy, &perm);
                ap_dimperm_clear(&perm);
            }
        } else if (op->op == PA_ALLOC || op->op == PA_ALLOC_N) {
            // x := new
            if (!conservative) {
                ap_dimperm_init(&perm, a->size);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_new(pr, a, vx, &perm);
                ap_dimperm_clear(&perm);
            }
            // else null
        } else if (op->op == PA_FREE) {
            // x := free
            if (conservative) {
                ap_dimperm_init(&perm, a->size);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_free(pr, a, vx, &perm, &anon);
                ap_dimperm_clear(&perm);
                hgraph_node_anon_array_clear(pr, &anon);
            }
            // else null
        } else if (IS_NULLDIM(op->info.ptr.y) ||
                ny == 0) {
            // x := null
            if (nx != 0 && conservative) {
                ap_dimperm_init(&perm, a->size);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_null(pr, a, vx, &perm, &anon);
                ap_dimperm_clear(&perm);
                hgraph_node_anon_array_clear(pr, &anon);
            } else {
                r = hgraph_copy_internal(pr, a);
            }
        } else {
            ERROR("Bad assignment: not an allowed right value!", return NULL;
                    );
        }
    } else if (op->offx == ((int) (a->datadim+OFFSET_NEXT))) {
        // x->next = [null | y]
        if (op->info.ptr.offy != OFFSET_NONE) {
            ERROR("Bad assignment: not an elementary assignment!",;
                    );
        }
        if (IS_NULLDIM(op->info.ptr.y) || ny == 0) { // x->next = null
            if (conservative) {
                ap_dimperm_init(&perm, a->size);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_next_null(pr, a, vx, &perm, &anon);
                ap_dimperm_clear(&perm);
                hgraph_node_anon_array_clear(pr, &anon);
            }
            // else return null
        } else if (IS_PTRDIM(op->info.ptr.y, a->datadim, a->ptrdim)) { // x->next = y
            if (conservative) {
                ap_dimperm_init(&perm, a->size);
                ap_dimperm_set_id(&perm);
                r = hgraph_assign_x_next_y(pr, a, vx, vy, &perm);
                ap_dimperm_clear(&perm);
            }
            // else return null
        } else {
            ERROR("Bad assignment: not an elementary assignment!",;
                    );
        }
    }
    return r;
}

hgraph_array_t *
hgraph_assign_passign_array(hgraph_internal_t * pr, bool destructive,
        hgraph_t * a, passign0_array_t * arr) {
    hgraph_array_t *r = hgraph_array_make(pr, 1);
    r->p[0] = (destructive) ? a : hgraph_copy_internal(pr, a);
    size_t i;
    for (i = 0; i < arr->size; i++) { // apply each operation of the list on all elements of r
        hgraph_array_t *rr = hgraph_array_make(pr, r->size);
        size_t j;
        for (j = 0; j < r->size; j++)
            if (r->p[j]) {
                hgraph_t *r1 =
                        hgraph_assign_passign(pr, true, r->p[i], arr->p[i]);
                if (r1)
                    hgraph_array_add(pr, rr, true, true, r1); /* copy & destroy */
                r1 = hgraph_assign_passign(pr, false, r->p[i], arr->p[i]);
                if (r1)
                    hgraph_array_add(pr, rr, true, true, r1); /* copy & destroy */
            }
        hgraph_array_clear(pr, r);
        r = rr;
    }
    return r;
}

/* Only done if real dimension (ptr vars) are used in constraints */
hgraph_t *
hgraph_assign_linexpr(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t d, ap_linexpr0_t * expr, hgraph_t * dest) {
    if (hgraph_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : hgraph_copy(man, a);
    else {
        hgraph_t *b;
        passign0_array_t *op;
        hgraph_array_t *r;
        hgraph_internal_t *pr =
                hgraph_init_from_manager(man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
        ap_dimension_t dim = hgraph_dimension(man, a);
        if (!destructive)
            b = hgraph_copy_mem(pr, a);
        else
            b = a;
        op =
                shape_passign_of_linexpr_array(pr, &d, &expr, 1, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        printf("\n========== shape assign with passign:\n");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
#endif
        /* go */
        r = hgraph_assign_passign_array(pr, false, b, op);
        hgraph_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            hgraph_t *rr = hgraph_meet(pr->man, false, b, dest);
            hgraph_free_internal(pr, b);
            return rr;
        }
        return b;
    }
}

/**
 * Special substitution of nodes through pointer variables.
 * Used in return from procedure.
 * @param size length of both arrays tdim, tnew
 * @param perm to be computed
 * @return
 */
hgraph_t *
hgraph_substitute_array(hgraph_internal_t* pr, hgraph_t* a,
        ap_dim_t *tdim, ap_dim_t *tnew,
        size_t size, ap_dimperm_t *perm) {
    hgraph_t* r;
    if (!a) return a;
    arg_assert(a && tdim && tnew && size >= 1, return NULL;);

    // Step 0: new graph with the same dimensions as a
    r = hgraph_alloc_internal(pr, a->size, a->datadim, a->ptrdim);

    // Step 1: fill the info field of r such that:
    // a- all variables pointing to the same nodes that tdim point now to
    //    the nodes of tnew
    // b- all nodes having as successor a node pointed by tdim points to
    //    the correspoding node of tnew
    ap_dimperm_t nperm; // gives the translation form old to new nodes
    ap_dimperm_init(&nperm, a->size);
    ap_dimperm_set_id(&nperm);
    size_t v, n;
    for (v = 0; v < size; v++) {
        n = VAR2NODE(a, tdim[v] - a->datadim);
        if (n != NODE_NULL)
            nperm.dim[n] = VAR2NODE(a, tnew[v] - a->datadim);
    }
    // 1(a)
    for (v = 0; v < a->ptrdim; v++) {
        VAR2NODE(r, v) = nperm.dim[VAR2NODE(a, v)];
    }
    // 1(b)
    for (n = 1; n < a->size; n++) {
        NODE_NEXT(r, n) = nperm.dim[NODE_NEXT(a, n)];
        // no need to update variables for nodes because this is done
        // by hgraph_close from the 0..ptrdim entries of info filled above
    }
#ifndef NDEBUG
    fprintf(stdout, "\n????hgraph_substitute_array: with array=(\n");
    for (v = 0; v < size; v++)
        fprintf(stdout, "%zu/%zu ", tnew[v] - a->datadim, tdim[v] - a->datadim);
    fprintf(stdout, "\n) generates hgraph with info=(\n");
    for (v = 0; v < a->ptrdim; v++)
        fprintf(stdout, "[v=%zu,n=%zu] ", v, VAR2NODE(r, v));
    for (n = 1; n < a->size; n++)
        fprintf(stdout, "[n=%zu,v=%zu] ", n, NODE_VAR(r, n));
    fprintf(stdout, "\n) for nperm = (\n");
    ap_dimperm_fprint(stdout, &nperm);
#endif
    // Step 2: close the graph
    ap_dim_array2_t anon;
    anon.size = 0;
    anon.p = NULL;
    hgraph_t *rr = hgraph_close(pr, r, perm, &anon, false);
    if (anon.size > 0) {
        ERROR("Anonymous nodes generated by substitution!", return NULL;);
        free(anon.p);
    }
    // free the intermediate memory
    ap_dimperm_clear(&nperm);
    hgraph_free_internal(pr, r);

    return rr;
}

/* TODO: priority 0 */

/* used for pre-image computation */
hgraph_t *
hgraph_substitute_linexpr(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t d, ap_linexpr0_t * expr, hgraph_t * dest) {
    hgraph_internal_t *pr =
            hgraph_init_from_manager(man, AP_FUNID_SUBSTITUTE_LINEXPR_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}

/* Only done if real dimension (ptr vars) are used in constraints */
hgraph_t *
hgraph_assign_linexpr_array(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t * tdim,
        ap_linexpr0_t ** texpr,
        size_t size, hgraph_t * dest) {
    if (hgraph_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : hgraph_copy(man, a);
    else {
        hgraph_t *b;
        passign0_array_t *op;
        hgraph_array_t *r;
        hgraph_internal_t *pr =
                hgraph_init_from_manager(man, AP_FUNID_ASSIGN_LINEXPR_ARRAY, 0);
        ap_dimension_t dim = hgraph_dimension(man, a);
        if (!destructive)
            b = hgraph_copy_mem(pr, a);
        else
            b = a;
        op =
                shape_passign_of_linexpr_array(pr, tdim, texpr, size, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        printf("\n========== shape assign with passign:\n");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
#endif
        /* go */
        r = hgraph_assign_passign_array(pr, false, b, op);
        hgraph_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            hgraph_t *rr = hgraph_meet(pr->man, false, b, dest);
            hgraph_free_internal(pr, b);
            return rr;
        }
        return b;
    }
}

/* TODO: priority 0 */

/* used for pre-image computation */
hgraph_t *
hgraph_substitute_linexpr_array(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t * tdim,
        ap_linexpr0_t ** texpr,
        size_t size, hgraph_t * dest) {
    if (size == 1)
        return hgraph_substitute_linexpr(man, destructive, a, tdim[0], texpr[0],
            dest);

    hgraph_internal_t *pr =
            hgraph_init_from_manager(man, AP_FUNID_SUBSTITUTE_TEXPR_ARRAY, 0);
    ap_manager_raise_exception(man, AP_EXC_NOT_IMPLEMENTED, pr->funid,
            "not implemented");
    return a;
}

/* Only for constraints between realdim */
hgraph_t *
hgraph_assign_texpr_array(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t * tdim,
        ap_texpr0_t ** texpr, size_t size, hgraph_t * dest) {
    if (hgraph_is_bottom(man, a))
        /* nothing to do */
        return (destructive) ? a : hgraph_copy(man, a);
    else {
        hgraph_t *b;
        passign0_array_t *op;
        hgraph_array_t *r;
        hgraph_internal_t *pr =
                shape_init_from_manager(man, AP_FUNID_ASSIGN_TEXPR_ARRAY, 0);
        ap_dimension_t dim = hgraph_dimension(man, a);
        if (!destructive)
            b = hgraph_copy_mem(pr, a);
        else
            b = a;
        op =
                shape_passign_of_texpr_array(pr, tdim, texpr, size, dim.intdim,
                dim.realdim);
#ifndef NDEBUG
        fprintf(stdout, "\n========== hgraph assign with passign: assigns \n");
        shape_passign_array_fdump(stdout, op, dim.intdim, dim.realdim);
        fprintf(stdout, "\n\t on ");
        hgraph_fdump(stdout, man, a);
        fflush(stdout);
#endif
        /* go */
        r = hgraph_assign_passign_array(pr, false, b, op);
#ifndef NDEBUG
        fprintf(stdout, "\n========== hgraph assign with passign: returns ");
        hgraph_array_fdump(stdout, man, r);
        fflush(stdout);
#endif
        hgraph_free_internal(pr, b);
        b = NULL;
        if (r)
            b = r->p[0];
        if (dest && b) { /* intersect r with dest */
            hgraph_t *rr = hgraph_meet(pr->man, false, b, dest);
            hgraph_free_internal(pr, b);
            return rr;
        }
        return b;
    }
}

/* TODO: priority 0 */

/* used only for pre-image computation */
hgraph_t *
hgraph_substitute_texpr_array(ap_manager_t * man,
        bool destructive, hgraph_t * a,
        ap_dim_t * tdim,
        ap_texpr0_t ** texpr,
        size_t size, hgraph_t * dest) {
    return ap_generic_substitute_texpr_array(man, destructive, a, tdim, texpr,
            size, dest);
}
