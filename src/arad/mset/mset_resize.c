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
#include "shape_macros.h"
#include "apron2shape.h"


/* ============================================================ */
/* Projections */
/* ============================================================ */

/* Forget nodes in tdim.
 * Each tdim[i] is in [a->datadim,a->datadim+a->segmdim).
 */
mset_t *
mset_forget_array(ap_manager_t * man,
        bool destructive, mset_t * a,
        ap_dim_t * tdim, size_t size, bool project) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_FORGET_ARRAY, 0);
    arg_assert(a && tdim && 0 < size && size < a->segmdim, return NULL;);
    // the dimension is kept
    mset_t *r = mset_alloc_internal(pr, a->datadim, a->segmdim);
    ap_dim_t *ntdim;
    size_t i;
    // GO for data and length dimensions
    ntdim = (ap_dim_t *) malloc(2 * size * sizeof (ap_dim_t));
    for (i = 0; i < size; i++) {
        ntdim[i] = tdim[i];
        ntdim[size + i] = a->segmdim + tdim[i];
    }
    r->dcons =
            ap_abstract0_forget_array(pr->man_dcons, false, a->dcons, ntdim,
            2 * size, project);
    // GO for multiset dimensions
    free(ntdim);
    ntdim = (ap_dim_t *) malloc(2 * size * sizeof (ap_dim_t));
    for (i = 0; i < size; i++) {
        ntdim[i] = tdim[i];
        ntdim[size + i] = tdim[i] + a->segmdim;
    }
    r->mscons =
            ap_abstract0_forget_array(pr->man_mscons, false, a->mscons, ntdim,
            size, project);
    free(ntdim);
    if (destructive)
        mset_free_internal(pr, a);
    return r;
}


/* ============================================================ */
/* Change and permutation of dimensions */
/* ============================================================ */

/* The dimchange may contain positions in datadim and at the end of segmentdim.
 */
mset_t *
mset_add_dimensions(ap_manager_t * man,
        bool destructive, mset_t * a,
        ap_dimchange_t * dimchange, bool project) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_ADD_DIMENSIONS, 0);
    size_t addr = dimchange->realdim;
    size_t addi = dimchange->intdim;
    arg_assert(a && dimchange
            && dimchange->dim != NULL
            , return NULL;);
    // this add a new dimension for data, for length, and for other abstractions (sum, multiset).
    mset_t *r = mset_alloc_internal(pr, a->datadim + addi, a->segmdim + addr);
    ap_dimchange_t dimadd;
    dimadd.realdim = 0;
    // GO for data and length dimensions
    dimadd.intdim = addi + 2 * addr;
    dimadd.dim = (ap_dim_t *) malloc(dimadd.intdim * sizeof (ap_dim_t));
    size_t i;
    for (i = 0; i < addi; i++)
        dimadd.dim[i] = dimchange->dim[i];
    for (i = 0; i < addr; i++) {
        arg_assert(dimchange->dim[addi + i] == (a->datadim + a->segmdim), return NULL;);
        dimadd.dim[addi + i] = a->datadim + a->segmdim;
        dimadd.dim[addi + addr + i] = a->datadim + 2 * a->segmdim;
    }
    r->dcons =
            ap_abstract0_add_dimensions(pr->man_dcons, false, a->dcons, &dimadd,
            project);
    // GO for multiset dimensions
    dimadd.intdim = addi + 2 * addr;
    dimadd.dim = (ap_dim_t *) malloc(dimadd.intdim * sizeof (ap_dim_t));
    for (i = 0; i < addi; i++)
        // TODO: is not working if length dimensions
        dimadd.dim[i] = dimchange->dim[i];
    for (i = 0; i < addr; i++) {
        arg_assert(dimchange->dim[addi + i] == (a->datadim + a->segmdim), return NULL;);
        dimadd.dim[addi + i] = DATA_DIM(a->datadim, 0) + a->segmdim;
        dimadd.dim[addi + addr + i] = DATA_DIM(a->datadim, 0) + 2 * a->segmdim;
    }
    r->mscons =
            ap_abstract0_add_dimensions(pr->man_mscons, false, a->mscons, &dimadd,
            project);
    free(dimadd.dim);
    if (destructive)
        mset_free_internal(pr, a);
    return r;
}

/* Used only to remove node dimensions,
 * so dimchange contains only values >= a->datadim.
 */
mset_t *
mset_remove_dimensions(ap_manager_t * man,
        bool destructive, mset_t * a,
        ap_dimchange_t * dimchange) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_REMOVE_DIMENSIONS, 0);
    arg_assert(dimchange, return NULL;);
    if (!a)
        return NULL;
    // for one segment to remove, remove several dimensions
    mset_t *r =
            mset_alloc_internal(pr, a->datadim - dimchange->intdim, a->segmdim - dimchange->realdim);
    size_t i;
    ap_dimchange_t dimrm;
    dimrm.realdim = 0;
    // GO for data and length dimensions
    dimrm.intdim = dimchange->intdim + 2 * dimchange->realdim;
    dimrm.dim = (ap_dim_t *) malloc(dimrm.intdim * sizeof (ap_dim_t));
    for (i = 0; i < dimchange->intdim; i++)
        dimrm.dim[i] = dimchange->dim[i];
    for (i = 0; i < dimchange->realdim; i++) {
        dimrm.dim[dimchange->intdim + 2 * i + 0] = dimchange->dim[i];
        dimrm.dim[dimchange->intdim + 2 * i + 1] = a->segmdim + dimchange->dim[i];
    }
    // dimrm shall be sorted!
    shape_dimchange_sort(&dimrm);
    r->dcons =
            ap_abstract0_remove_dimensions(pr->man_dcons, false, a->dcons, &dimrm);
    free(dimrm.dim);
    // GO for multiset dimensions
    dimrm.intdim = dimchange->intdim + 2 * dimchange->realdim;
    dimrm.dim = (ap_dim_t *) malloc(dimrm.intdim * sizeof (ap_dim_t));
    for (i = 0; i < dimchange->intdim; i++)
        dimrm.dim[i] = dimchange->dim[i];
    for (i = 0; i < dimchange->realdim; i++) {
        dimrm.dim[dimchange->intdim + 2 * i + 0] = dimchange->dim[i];
        dimrm.dim[dimchange->intdim + 2 * i + 1] =
                a->segmdim + dimchange->dim[i];
    }
    shape_dimchange_sort(&dimrm);
    r->mscons =
            ap_abstract0_remove_dimensions(pr->man_mscons, false, a->mscons, &dimrm);
    free(dimrm.dim);

    if (destructive)
        mset_free_internal(pr, a);
    return r;
}

/* Intensively used to permute nodes.
 * perm->size == a->datadim+a->segmdim
 * Nodes mapped to a->datadim has to be removed */
mset_t *
mset_permute_dimensions(ap_manager_t * man,
        bool destructive, mset_t * a, ap_dimperm_t * perm) {
    mset_internal_t *pr =
            mset_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
    if (!a)
        return NULL;
    assert(perm && perm->size == (a->datadim + a->segmdim));
    mset_t *r;
#ifndef NDEBUG2
    fprintf(stdout, "\n====mset_permute_dimensions: with a=(");
    mset_fdump(stdout, man, a);
    fprintf(stdout, ") and perm=[");
    ap_dimperm_fprint(stdout, perm);
    fprintf(stdout, "]\n");
    fflush(stdout);
#endif
    /*
     * Remove node dimensions mapped to a->datadim (except a->datadim itself) and build a new permutation
     */
    ap_dimchange_t dimrm;
    ap_dimperm_t newperm;
    size_t i, j, rmsize;
    dimrm.intdim = 0;
    dimrm.realdim = perm->size;
    dimrm.dim = (ap_dim_t *) malloc(perm->size * sizeof (ap_dim_t));
    newperm.size = perm->size;
    newperm.dim = (ap_dim_t *) malloc(perm->size * sizeof (ap_dim_t));
    for (i = 0, j = 0, rmsize = 0; i < perm->size; i++)
        if (i > a->datadim && perm->dim[i] == a->datadim) {
            dimrm.dim[rmsize] = i;
            rmsize++;
        } else {
            newperm.dim[j] = perm->dim[i];
            j++;
        }
    if (rmsize > 0) {
        dimrm.intdim = 0;
        dimrm.realdim = rmsize;
        dimrm.dim =
                (ap_dim_t *) realloc(dimrm.dim, rmsize * sizeof (ap_dim_t));
        // Warning: forget is not needed, the following code produce an error!
        // r = mset_forget_array (man, false, a, dimrm.dim, rmsize, true);
        // r = mset_remove_dimensions (man, true, r, &dimrm);
        r = mset_remove_dimensions(man, false, a, &dimrm);
        pr = mset_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
        newperm.size = perm->size - rmsize;
        newperm.dim =
                (ap_dim_t *) realloc(newperm.dim, newperm.size * sizeof (ap_dim_t));
    } else
        r = mset_copy_internal(pr, a);
    free(dimrm.dim);
    /*
     * Build permutations corresponding to each set of constraints.
     */
    ap_dimperm_t consperm;
    size_t newsegmsize = (newperm.size - a->datadim);
    // GO for data and length dimensions
    consperm.size = a->datadim + 2 * newsegmsize;
    consperm.dim = (ap_dim_t *) malloc(consperm.size * sizeof (ap_dim_t));
    ap_dimperm_set_id(&consperm);
    for (i = 0; i < newsegmsize; i++) {
        consperm.dim[a->datadim + i] = newperm.dim[a->datadim + i];
        consperm.dim[a->datadim + newsegmsize + i] =
                newsegmsize + newperm.dim[a->datadim + i];
    }
    r->dcons =
            ap_abstract0_permute_dimensions(pr->man_dcons, true, r->dcons,
            &consperm);
    free(consperm.dim);
    // GO for multiset constraints
    consperm.size = DATA_DIM(a->datadim, 0) + 2 * newsegmsize;
    consperm.dim = (ap_dim_t *) malloc(consperm.size * sizeof (ap_dim_t));
    ap_dimperm_set_id(&consperm);
    for (i = 0; i < newsegmsize; i++) {
        consperm.dim[DATA_DIM(a->datadim, 0) + i] =
                newperm.dim[a->datadim + i];
        consperm.dim[DATA_DIM(a->datadim, 0) + newsegmsize + i] =
                newsegmsize + newperm.dim[a->datadim + i];
    }
    r->mscons =
            ap_abstract0_permute_dimensions(pr->man_mscons, true, r->mscons,
            &consperm);
    free(consperm.dim);
    free(newperm.dim);
    if (destructive)
        mset_free_internal(pr, a);
#ifndef NDEBUG2
    fprintf(stdout, "\n====mset_permute_dimensions: returns r=(");
    mset_fdump(stdout, man, r);
    fprintf(stdout, ")\n");
    fflush(stdout);
#endif

    return r;
}


/* ============================================================ */
/* Expansion and folding of dimensions */

/* ============================================================ */


mset_t *
mset_singleton(mset_internal_t * pr, bool destructive, mset_t * a,
        ap_dim_t dim) {
    /*
     * intersect with constraints l[dim] = 1 and mstl(dim) = 0
     */
    mset_t *r = mset_alloc_internal(pr, a->datadim, a->segmdim);
    ap_lincons0_array_t arr = ap_lincons0_array_make(1);
    // l[dim]-1 == 0
    arr.p[0].constyp = AP_CONS_EQ;
    arr.p[0].linexpr0 =
            ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + 2 * a->segmdim);
    arr.p[0].scalar = NULL;
    ap_linexpr0_set_coeff_scalar_int(arr.p[0].linexpr0,
            a->datadim + a->segmdim + dim, 1);
    ap_linexpr0_set_cst_scalar_int(arr.p[0].linexpr0, -1);
    r->dcons =
            ap_abstract0_meet_lincons_array(pr->man_dcons, false, a->dcons, &arr);
    ap_lincons0_array_clear(&arr);
    // M[dim]==0
    arr = ap_lincons0_array_make(1);
    arr.p[0].constyp = AP_CONS_EQ;
    arr.p[0].linexpr0 =
            ap_linexpr0_alloc(AP_LINEXPR_DENSE,
            DATA_DIM(a->datadim, 0) + 2 * a->segmdim);
    arr.p[0].scalar = NULL;
    ap_linexpr0_set_coeff_scalar_int(arr.p[0].linexpr0,
            DATA_DIM(a->datadim,
            0) + a->segmdim + dim, 1);
    ap_linexpr0_set_cst_scalar_int(arr.p[0].linexpr0, 0);
    r->mscons =
            ap_abstract0_meet_lincons_array(pr->man_mscons, false, a->mscons, &arr);
    ap_lincons0_array_clear(&arr);
    mset_strengthen_dim(pr, r, dim);

    return r;
}

mset_t *
mset_unfold(mset_internal_t * pr, bool destructive, mset_t * a, ap_dim_t dim) {
    /*
     * The node added by expansion (n2) has dimension a->segmdim
     * Then apply substitutions/assignments:
     * mstl(n1) == mstl(n2) + mshd[n2];
     * mstl(n1) := 0;
     * l[n1] == l[n2] + 1;
     * l[n1] := 1;
     * meet l[n2]>=1
     */
    // Add dimension a->segmdim
    ap_dimchange_t dimadd;
    ap_dimchange_init(&dimadd, 0, 1);
    dimadd.dim = (ap_dim_t *) malloc(1 * sizeof (ap_dim_t));
    dimadd.dim[0] = a->datadim + a->segmdim;
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_unfold: on\n");
    mset_fdump(stdout, pr->man, a);
    fflush(stdout);
#endif
    mset_t *r = mset_add_dimensions(pr->man, false, a, &dimadd, false);
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_unfold: after add dimension\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    ap_dimchange_clear(&dimadd);
    // Build statements and apply them
    ap_dim_t n1 = dim;
    ap_dim_t n2 = a->segmdim;
    ap_dim_t msn1, msn2, dn2, ln1, ln2;
    ap_linexpr0_t *expr;
    // GO for multiset constraints
    msn1 = DATA_DIM(r->datadim, 0) + r->segmdim + n1;
    msn2 = DATA_DIM(r->datadim, 0) + r->segmdim + n2;
    dn2 = DATA_DIM(r->datadim, 0) + n2;
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, DATA_DIM(r->datadim, 0) + 2 * r->segmdim);
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_unfold: before substitute mstl(n1)\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    // mstl(n1) == mstl(n2) + mshd[n2]
    ap_linexpr0_set_coeff_scalar_int(expr, msn2, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
    r->mscons =
            ap_abstract0_substitute_linexpr(pr->man_mscons, true, r->mscons, msn1, expr,
            NULL);
#ifndef NDEBUG
    fprintf(stdout,
            "\n====mset_unfold: after substitute mstl(n1), before assign mstl(n1)\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    // mstl(n1) := 0
    ap_linexpr0_set_coeff_scalar_int(expr, msn2, 0);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 0);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
    r->mscons =
            ap_abstract0_assign_linexpr(pr->man_mscons, true, r->mscons, msn1, expr,
            NULL);
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_unfold: after assign mstl(n1)\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif

    // GO for length constraints
    ln1 = r->datadim + r->segmdim + n1;
    ln2 = r->datadim + r->segmdim + n2;
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 2 * r->segmdim);
    // l[n1] == l[n2] + 1;
    ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, 1);
    r->dcons =
            ap_abstract0_substitute_linexpr(pr->man_dcons, true, r->dcons, ln1, expr,
            NULL);
    // l[n1] := 1;
    ap_linexpr0_set_coeff_scalar_int(expr, ln2, 0);
    r->dcons =
            ap_abstract0_assign_linexpr(pr->man_dcons, true, r->dcons, ln1, expr,
            NULL);
    // meet l[n2]-1 >= 0
    ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, -1);
    ap_lincons0_array_t arr = ap_lincons0_array_make(1);
    arr.p[0] = ap_lincons0_make(AP_CONS_SUPEQ, expr, NULL);
    r->dcons =
            ap_abstract0_meet_lincons_array(pr->man_dcons, true, r->dcons, &arr);
    ap_lincons0_array_clear(&arr); // free also expr
    mset_strengthen_dim(pr, r, r->datadim + n1);
    return r;
}

mset_t*
mset_split(mset_internal_t* pr, bool destructive, mset_t* a, ap_dim_t dim) {
    /*
     * The node (n2) added by expansion has dimension a->segmdim.
     * Then apply substitutions/assignments:
     * mstl(n1) == mstl(n1) + mshd[n2] + mstl(n2);
     * l[n1] == l[n1] + l[n2];
     * meet l[n1]>=1 && l[n2]>=1
     */
    // Add dimension a->segmdim
    ap_dimchange_t dimadd;
    ap_dimchange_init(&dimadd, 0, 1);
    dimadd.dim = (ap_dim_t *) malloc(1 * sizeof (ap_dim_t));
    dimadd.dim[0] = a->datadim + a->segmdim;
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_split: on\n");
    mset_fdump(stdout, pr->man, a);
    fflush(stdout);
#endif
    mset_t *r = mset_add_dimensions(pr->man, false, a, &dimadd, false);
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_split: after add dimension\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    ap_dimchange_clear(&dimadd);
    // Build statements and apply them
    ap_dim_t n1 = dim;
    ap_dim_t n2 = a->segmdim;
    ap_dim_t msn1, msn2, dn2, ln1, ln2;
    ap_linexpr0_t *expr;
    // GO for multiset constraints
    msn1 = DATA_DIM(r->datadim, 0) + r->segmdim + n1;
    msn2 = DATA_DIM(r->datadim, 0) + r->segmdim + n2;
    dn2 = DATA_DIM(r->datadim, 0) + n2;
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, DATA_DIM(r->datadim, 0) + 2 * r->segmdim);
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_split: before substitute mstl(n1)\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    // mstl(n1) == msttl(n1) + mstl(n2) + mshd[n2]
    ap_linexpr0_set_coeff_scalar_int(expr, msn1, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, msn2, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
    r->mscons =
            ap_abstract0_substitute_linexpr(pr->man_mscons, true, r->mscons, msn1, expr,
            NULL);
#ifndef NDEBUG
    fprintf(stdout,
            "\n====mset_split: after substitute mstl(n1)\n");
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    ap_linexpr0_free(expr);

    // GO for length constraints
    ln1 = r->datadim + r->segmdim + n1;
    ln2 = r->datadim + r->segmdim + n2;
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 2 * r->segmdim);
    // l[n1] == l[n1] + l[n2];
    ap_linexpr0_set_coeff_scalar_int(expr, ln1, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
    r->dcons =
            ap_abstract0_substitute_linexpr(pr->man_dcons, true, r->dcons, ln1, expr,
            NULL);
    // meet l[n1]-1>=0 && l[n2]-1>=0
    ap_lincons0_array_t arr = ap_lincons0_array_make(2);
    ap_linexpr0_set_coeff_scalar_int(expr, ln1, 0);
    ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, -1);
    arr.p[0] = ap_lincons0_make(AP_CONS_SUPEQ, expr, NULL);
    ap_linexpr0_t* expr1 = ap_linexpr0_alloc(AP_LINEXPR_DENSE,
            DATA_DIM(r->datadim, 0) + 2 * r->segmdim);
    ap_linexpr0_set_coeff_scalar_int(expr1, ln1, 1);
    ap_linexpr0_set_cst_scalar_int(expr1, -1);
    arr.p[1] = ap_lincons0_make(AP_CONS_SUPEQ, expr1, NULL);
    r->dcons =
            ap_abstract0_meet_lincons_array(pr->man_dcons, true, r->dcons, &arr);
    ap_lincons0_array_clear(&arr); // free also expr
    mset_strengthen_dim(pr, r, r->datadim + n1);
    return r;
}

/*
 * Expand dimension @code{dim} using code @code{n} such that
 * - n=1 simulate the singleton operation,
 * - n=2 simulate the unfold operation, and
 * - n=3 simulate the split operation.
 */
mset_t *
mset_expand(ap_manager_t * man,
        bool destructive, mset_t * a, ap_dim_t dim, size_t n) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_EXPAND, 0);
    if (!a)
        return NULL;
    arg_assert(n <= 3 && dim > 0 && dim < a->segmdim, return NULL;);
    mset_t *r;
    switch (n) {
        case 1:
            r = mset_singleton(pr, destructive, a, dim);
            break;
        case 2:
            r = mset_unfold(pr, destructive, a, dim);
            break;
        default:
            r = mset_split(pr, destructive, a, dim);
    }
    if (destructive)
        mset_free_internal(pr, a);
#ifndef NDEBUG
    fprintf(stdout, "\n====mset_expand: code=(%zu) on\n", n);
    mset_fdump(stdout, pr->man, r);
    fflush(stdout);
#endif
    return r;
}

/* Used for the merge operation:
 * the first node of tdim is the cut node, the other nodes are anonymous nodes.
 * The size of the abstract value is not changed (it will be changed by permute!).
 * Values of tdim are in [a->datadim, a->datadim+a->segmdim).
 * Warning: if size==1 (folding of several segments)
 *          then return the object
 */
mset_t *
mset_fold(ap_manager_t * man,
        bool destructive, mset_t * a, ap_dim_t * tdim, size_t size) {
    mset_internal_t *pr = mset_init_from_manager(man, AP_FUNID_FOLD, 0);
    if (!a)
        return NULL;
    arg_assert(tdim && 0 < size && size < a->segmdim, return NULL;);

    if (size == 1)
        /* mark the end of the set of anonymous nodes to be folded, NOT USED HERE */
        return a;

    /*
     * The operations done are:
     * M(tdim[0]) := M(tdim[0]) + d(tdim[1]) + M(tdim[1]) +... + d(tdim[size-1]) + M(tdim[size-1])
     * l[tdim[0]] := l[tdim[0]] + ... + l[tdim[size-1]] (not directly size because of summary nodes!)
     */
    ap_linexpr0_t *expr;
    mset_t *r = mset_alloc_internal(pr, a->datadim, a->segmdim);
    size_t i, basedim;

#ifndef NDEBUG
    fprintf(stdout, "\n mset_fold: \n\t dimensions ");
    for (i = 0; i < size; i++)
        fprintf(stdout, "%d, ", tdim[i]);
    fprintf(stdout, "\n\t on ");
    mset_fdump(stdout, man, a);
    fflush(stdout);
#endif

    /* M(tdim[0]) := M(tdim[0]) + d(tdim[1]) + M(tdim[1]) +... + d(tdim[size-1]) + M(tdim[size-1]) */
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE,
            DATA_DIM(a->datadim, 0) + 2 * a->segmdim);
    basedim = a->segmdim;
    ap_linexpr0_set_coeff_scalar_int(expr, basedim + tdim[0], 1);
    for (i = 1; i < size; i++) {
        ap_linexpr0_set_coeff_scalar_int(expr,
                tdim[i],
                1);
        ap_linexpr0_set_coeff_scalar_int(expr, basedim + tdim[i], 1);
    }
#ifndef NDEBUG
    fprintf(stdout, "\n\t built dexpr ");
    ap_linexpr0_fprint(stdout, expr, NULL);
    fflush(stdout);
#endif

    r->mscons =
            ap_abstract0_assign_linexpr(pr->man_mscons, false, a->mscons,
            basedim + tdim[0], expr, NULL);
    ap_linexpr0_free(expr);


    /* l[tdim[0]] := l[tdim[0]] + ... + l[tdim[size-1]] */
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + 2 * a->segmdim);
    basedim = a->segmdim;
    for (i = 0; i < size; i++)
        ap_linexpr0_set_coeff_scalar_int(expr, basedim + tdim[i], 1);
#ifndef NDEBUG
    fprintf(stdout, "\n\t built lexpr ");
    ap_linexpr0_fprint(stdout, expr, NULL);
    fflush(stdout);
#endif
    r->dcons =
            ap_abstract0_assign_linexpr(pr->man_dcons, false, a->dcons,
            basedim + tdim[0], expr, NULL);
    ap_linexpr0_free(expr);
    if (destructive)
        mset_free_internal(pr, a);
#ifndef NDEBUG
    fprintf(stdout, "\n\t returns ");
    mset_fdump(stdout, man, r);
    fflush(stdout);
#endif
    return r;
}
