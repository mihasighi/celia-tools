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


#include "lsum.h"
#include "lsum_internal.h"
#include "shape_macros.h"
#include "apron2shape.h"


/* ============================================================ */
/* Projections */
/* ============================================================ */

/* Forget nodes in tdim.
 * Each tdim[i] is in [a->datadim,a->datadim+a->segmdim).
 */
lsum_t *
lsum_forget_array(ap_manager_t * man,
        bool destructive, lsum_t * a,
        ap_dim_t * tdim, size_t size, bool project) {
    lsum_internal_t *pr =
            lsum_init_from_manager(man, AP_FUNID_FORGET_ARRAY, 0);
    arg_assert(a && tdim && 0 < size && size < a->segmdim, return NULL;);
    // the dimension is kept
    lsum_t *r = lsum_alloc_internal(pr, a->datadim, a->segmdim);
    ap_dim_t *ntdim;
    size_t i;
    ntdim = (ap_dim_t *) malloc(3 * size * sizeof (ap_dim_t));
    for (i = 0; i < size; i++) {
        ntdim[i] = tdim[i];
        ntdim[size + i] = a->segmdim + tdim[i];
        ntdim[2 * size + i] = 2 * a->segmdim + tdim[i];
    }
    r->dcons =
            ap_abstract0_forget_array(pr->man_dcons, false, a->dcons, ntdim,
            3 * size, project);
    free(ntdim);
    if (destructive)
        lsum_free_internal(pr, a);
    return r;
}


/* ============================================================ */
/* Change and permutation of dimensions */
/* ============================================================ */

/* The dimchange may contain positions in datadim and at the end of segmentdim.
 */
lsum_t *
lsum_add_dimensions(ap_manager_t * man,
        bool destructive, lsum_t * a,
        ap_dimchange_t * dimchange, bool project) {
    lsum_internal_t *pr =
            lsum_init_from_manager(man, AP_FUNID_ADD_DIMENSIONS, 0);
    size_t addr = dimchange->realdim;
    size_t addi = dimchange->intdim;
    arg_assert(a && dimchange
            && dimchange->dim != NULL
            , return NULL;);
    // this add a new dimension for data, for sum, and for length
    lsum_t *r = lsum_alloc_internal(pr, a->datadim + addi, a->segmdim + addr);
    ap_dimchange_t dimadd;
    dimadd.realdim = 0;
    dimadd.intdim = addi + 3 * addr;
    dimadd.dim = (ap_dim_t *) malloc(dimadd.intdim * sizeof (ap_dim_t));
    size_t i;
    for (i = 0; i < addi; i++)
        dimadd.dim[i] = dimchange->dim[i];
    for (i = 0; i < addr; i++) {
        arg_assert(dimchange->dim[addi + i] == (a->datadim + a->segmdim), return NULL;);
        dimadd.dim[addi + i] = a->datadim + a->segmdim;
        dimadd.dim[addi + addr + i] = a->datadim + 2 * a->segmdim;
        dimadd.dim[addi + 2 * addr + i] = a->datadim + 3 * a->segmdim;
    }
    r->dcons =
            ap_abstract0_add_dimensions(pr->man_dcons, false, a->dcons, &dimadd,
            project);
    if (destructive)
        lsum_free_internal(pr, a);
    return r;
}

/* Used to remove data or node dimensions.
 */
lsum_t *
lsum_remove_dimensions(ap_manager_t * man,
        bool destructive, lsum_t * a,
        ap_dimchange_t * dimchange) {
    lsum_internal_t *pr =
            lsum_init_from_manager(man, AP_FUNID_REMOVE_DIMENSIONS, 0);
    arg_assert(dimchange, return NULL;);
    if (!a)
        return NULL;
    // for one dimension to remove, remove 3 dimensions
    lsum_t *r =
            lsum_alloc_internal(pr, a->datadim - dimchange->intdim, a->segmdim - dimchange->realdim);
    size_t i;
    ap_dimchange_t dimrm;
    dimrm.realdim = 0;
    dimrm.intdim = dimchange->intdim + 3 * dimchange->realdim;
    dimrm.dim = (ap_dim_t *) malloc(dimrm.intdim * sizeof (ap_dim_t));
    for (i = 0; i < dimchange->intdim; i++)
        dimrm.dim[i] = dimchange->dim[i];
    for (i = 0; i < dimchange->realdim; i++) {
        dimrm.dim[dimchange->intdim + 3 * i + 0] = dimchange->dim[dimchange->intdim + i];
        dimrm.dim[dimchange->intdim + 3 * i + 1] = a->segmdim + dimchange->dim[dimchange->intdim + i];
        dimrm.dim[dimchange->intdim + 3 * i + 2] = 2 * a->segmdim + dimchange->dim[dimchange->intdim + i];
    }
    // dimrm shall be sorted!
    shape_dimchange_sort(&dimrm);
    r->dcons =
            ap_abstract0_remove_dimensions(pr->man_dcons, false, a->dcons, &dimrm);
    // ap_abstract0_forget_array(pr->man_dcons, a->dcons, &dimrm.dim, dimrm.intdim, true);
    if (destructive)
        lsum_free_internal(pr, a);
    return r;
}

/* Intensively used to permute nodes.
 * perm->size == a->datadim+a->segmdim
 * Nodes mapped to a->datadim has to be removed */
lsum_t *
lsum_permute_dimensions(ap_manager_t * man,
        bool destructive, lsum_t * a, ap_dimperm_t * perm) {
    lsum_internal_t *pr =
            lsum_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
    if (!a)
        return NULL;
    assert(perm && perm->size == (a->datadim + a->segmdim));
    lsum_t *r;
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_permute_dimensions: ");
#endif
    /*
     * Remove node dimensions mapped to a->datadim (except a->datadim itself) and build a new permutation
     */
    ap_dimchange_t dimrm;
    ap_dimperm_t newperm;
    size_t i, j, k, rmsize, newsize;
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
        // r = lsum_forget_array (man, false, a, dimrm.dim, rmsize, true);
        // r = lsum_remove_dimensions (man, true, r, &dimrm);
        r = lsum_remove_dimensions(man, false, a, &dimrm);
        pr = lsum_init_from_manager(man, AP_FUNID_PERMUTE_DIMENSIONS, 0);
        newperm.size = perm->size - rmsize;
        newperm.dim =
                (ap_dim_t *) realloc(newperm.dim, newperm.size * sizeof (ap_dim_t));
    } else
        r = lsum_copy_internal(pr, a);
    free(dimrm.dim);
    /*
     * Build permutations corresponding to each set of constraints.
     */
    ap_dimperm_t consperm;
    size_t newsegmsize = (newperm.size - a->datadim);
    consperm.size = a->datadim + 3 * newsegmsize;
    consperm.dim = (ap_dim_t *) malloc(consperm.size * sizeof (ap_dim_t));
    ap_dimperm_set_id(&consperm);
    for (i = 0; i < newsegmsize; i++) {
        consperm.dim[a->datadim + i] = newperm.dim[a->datadim + i];
        consperm.dim[a->datadim + newsegmsize + i] =
                newsegmsize + newperm.dim[a->datadim + i];
        consperm.dim[a->datadim + 2 * newsegmsize + i] =
                2 * newsegmsize + newperm.dim[a->datadim + i];
    }
#ifndef NDEBUG2
    fprintf(stdout, "\n\tperm=");
    ap_dimperm_fprint(stdout, perm);
    fprintf(stdout, "\n\tnewperm=");
    ap_dimperm_fprint(stdout, &newperm);
    fprintf(stdout, "\n\tconsperm=");
    ap_dimperm_fprint(stdout, &consperm);
    fprintf(stdout, "\n");
#endif
    r->dcons =
            ap_abstract0_permute_dimensions(pr->man_dcons, true, r->dcons,
            &consperm);
    free(consperm.dim);
    free(newperm.dim);
    if (destructive)
        lsum_free_internal(pr, a);
    return r;
}


/* ============================================================ */
/* Expansion and folding of dimensions */

/* ============================================================ */


lsum_t *
lsum_singleton(lsum_internal_t * pr, bool destructive, lsum_t * a,
        ap_dim_t dim) {
    /*
     * intersect with constraints l[dim] = 1 and S(dim) = 0
     */
    lsum_t *r = lsum_alloc_internal(pr, a->datadim, a->segmdim);
    ap_lincons0_array_t arr = ap_lincons0_array_make(2);
    // l[dim]-1 == 0
    arr.p[0].constyp = AP_CONS_EQ;
    arr.p[0].linexpr0 =
            ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + 3 * a->segmdim);
    arr.p[0].scalar = NULL;
    ap_linexpr0_set_coeff_scalar_int(arr.p[0].linexpr0,
            a->datadim + 2 * a->segmdim + dim, 1);
    ap_linexpr0_set_cst_scalar_int(arr.p[0].linexpr0, -1);
    // S(dim) == 0
    arr.p[1].constyp = AP_CONS_EQ;
    arr.p[1].linexpr0 =
            ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + 3 * a->segmdim);
    arr.p[1].scalar = NULL;
    ap_linexpr0_set_coeff_scalar_int(arr.p[1].linexpr0,
            a->datadim + a->segmdim + dim, 1);
    ap_linexpr0_set_cst_scalar_int(arr.p[1].linexpr0, 0);
    r->dcons =
            ap_abstract0_meet_lincons_array(pr->man_dcons, false, a->dcons, &arr);
    ap_lincons0_array_clear(&arr);
    return r;
}

lsum_t *
lsum_unfold(lsum_internal_t * pr, bool destructive, lsum_t * a, ap_dim_t dim) {
    /*
     * The node added by expansion has dimension a->segmdim
     * Then apply substitutions/assignments:
     * S(n1) == S(n2) + d[n2];
     * S(n1) := 0;
     * l[n1] == l[n2] + 1;
     * l[n1] := 1;
     * meet l[n2]>=1
     */
    // Add dimension a->segmdim
    ap_dimchange_t dimadd;
    ap_dimchange_init(&dimadd, 0, 1);
    dimadd.dim = (ap_dim_t *) malloc(1 * sizeof (ap_dim_t));
    dimadd.dim[0] = a->datadim + a->segmdim;
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_unfold: on\n");
    lsum_fdump(stdout, pr->man, a);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    lsum_t *r = lsum_add_dimensions(pr->man, false, a, &dimadd, false);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_unfold: after add dimension\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    ap_dimchange_clear(&dimadd);
    // Build statements and apply them
    ap_dim_t n1 = dim;
    ap_dim_t n2 = a->segmdim;
    ap_dim_t sn1, sn2, dn2, ln1, ln2;
    ap_linexpr0_t *expr;
    sn1 = r->datadim + r->segmdim + n1;
    sn2 = r->datadim + r->segmdim + n2;
    dn2 = r->datadim + n2;
    ln1 = r->datadim + 2 * r->segmdim + n1;
    ln2 = r->datadim + 2 * r->segmdim + n2;
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 3 * r->segmdim);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_unfold: before substitute S(n1)\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    // S(n1) == S(n2) + d[n2]
    ap_linexpr0_set_coeff_scalar_int(expr, sn2, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
    r->dcons =
            ap_abstract0_substitute_linexpr(pr->man_dcons, true, r->dcons, sn1, expr,
            NULL);
#ifndef NDEBUG2
    fprintf(stdout,
            "\n====lsum_unfold: - after substitute S(n1), before assign S(n1)\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    // S(n1) := 0
    ap_linexpr0_set_coeff_scalar_int(expr, sn2, 0);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 0);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
    r->dcons =
            ap_abstract0_assign_linexpr(pr->man_dcons, true, r->dcons, sn1, expr,
            NULL);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_unfold: after assign S(n1)\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
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
    ap_lincons0_array_clear(&arr);
    return r;
}

lsum_t*
lsum_split(lsum_internal_t* pr, bool destructive, lsum_t* a, ap_dim_t dim,
	   size_t idim) {
    /*
     * The node (n2) added by expansion has dimension a->segmdim.
     * Then apply substitutions/assignments:
     * S(n1) == S(n1) + d[n2] + S(n2);
     * l[n1] == l[n1] + l[n2];
     * if idim == 0 then meet l[n1]>=1 && l[n2]>=1
     * else meet l[n1]==(idim-1) && l[n2]>=1
     */
    // Add dimension a->segmdim
    ap_dimchange_t dimadd;
    ap_dimchange_init(&dimadd, 0, 1);
    dimadd.dim = (ap_dim_t *) malloc(1 * sizeof (ap_dim_t));
    dimadd.dim[0] = a->datadim + a->segmdim;
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_split: on\n");
    lsum_fdump(stdout, pr->man, a);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    lsum_t *r = lsum_add_dimensions(pr->man, false, a, &dimadd, false);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_split: after add dimension\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    ap_dimchange_clear(&dimadd);
    // Build statements and apply them
    ap_dim_t n1 = dim;
    ap_dim_t n2 = a->segmdim;
    ap_dim_t sn1, sn2, dn2, ln1, ln2;
    ap_linexpr0_t *expr;
    sn1 = r->datadim + r->segmdim + n1;
    sn2 = r->datadim + r->segmdim + n2;
    dn2 = r->datadim + n2;
    ln1 = r->datadim + 2 * r->segmdim + n1;
    ln2 = r->datadim + 2 * r->segmdim + n2;
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, r->datadim + 3 * r->segmdim);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_split: before substitute S(n1)\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    // S(n1) == S(n1) + d[n2] + S(n2)
    ap_linexpr0_set_coeff_scalar_int(expr, sn1, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, sn2, 1);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
    r->dcons =
            ap_abstract0_substitute_linexpr(pr->man_dcons, true, r->dcons, sn1, expr,
            NULL);
#ifndef NDEBUG2
    fprintf(stdout,
            "\n====lsum_split: - after substitute S(n1)\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    ap_linexpr0_set_coeff_scalar_int(expr, sn1, 0);
    ap_linexpr0_set_coeff_scalar_int(expr, sn2, 0);
    ap_linexpr0_set_coeff_scalar_int(expr, dn2, 0);
    ap_linexpr0_set_cst_scalar_int(expr, 0);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_split: - before substitute l[n1]\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    // l[n1] == l[n1] + l[n2];
    ap_linexpr0_set_coeff_scalar_int(expr, ln1, 1);
    ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
    r->dcons =
            ap_abstract0_substitute_linexpr(pr->man_dcons, true, r->dcons, ln1, expr,
            NULL);
#ifndef NDEBUG2
    fprintf(stdout,
            "\n====lsum_split: - after substitute l[n1], before meet\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    if (idim == 0)
      {
	// meet l[n1]-1>=0 && l[n2]-1>=0
	ap_lincons0_array_t arr = ap_lincons0_array_make(2);
	ap_linexpr0_set_coeff_scalar_int(expr, ln1, 0);
	ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
	ap_linexpr0_set_cst_scalar_int(expr, -1);
	arr.p[0] = ap_lincons0_make(AP_CONS_SUPEQ, expr, NULL);
	ap_linexpr0_t* expr1 = ap_linexpr0_alloc(AP_LINEXPR_DENSE,
						 DATA_DIM(r->datadim, 0) + 3 * r->segmdim);
	ap_linexpr0_set_coeff_scalar_int(expr1, ln1, 1);
	ap_linexpr0_set_cst_scalar_int(expr1, -1);
	arr.p[1] = ap_lincons0_make(AP_CONS_SUPEQ, expr1, NULL);
	r->dcons =
	  ap_abstract0_meet_lincons_array(pr->man_dcons, true, r->dcons, &arr);
	ap_lincons0_array_clear(&arr);
      }
    else
      {
	// meet l[n1]-(idim-1)==0 && l[n2]-1>=0
	ap_lincons0_array_t arr = ap_lincons0_array_make(2);
	ap_linexpr0_set_coeff_scalar_int(expr, ln1, 0);
	ap_linexpr0_set_coeff_scalar_int(expr, ln2, 1);
	ap_linexpr0_set_cst_scalar_int(expr, -1);
	arr.p[0] = ap_lincons0_make(AP_CONS_SUPEQ, expr, NULL);
	ap_linexpr0_t* expr1 = ap_linexpr0_alloc(AP_LINEXPR_DENSE,
						 DATA_DIM(r->datadim, 0) + 3 * r->segmdim);
	ap_linexpr0_set_coeff_scalar_int(expr1, ln1, 1);
	ap_linexpr0_set_coeff_scalar_int(expr1, (idim-1), -1);
	ap_linexpr0_set_cst_scalar_int(expr1, 0);
	arr.p[1] = ap_lincons0_make(AP_CONS_EQ, expr1, NULL);
	r->dcons =
	  ap_abstract0_meet_lincons_array(pr->man_dcons, true, r->dcons, &arr);
	ap_lincons0_array_clear(&arr);
      }
#ifndef NDEBUG2
    fprintf(stdout,
            "\n====lsum_split: - before return\n");
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return r;
}

/*
 * Expand dimension @code{dim} using code @code{n} such that
 * - n=1 simulate the singleton operation,
 * - n=2 simulate the unfold operation, and
 * - n=3 simulate the split operation.
 */
lsum_t *
lsum_expand(ap_manager_t * man,
        bool destructive, lsum_t * a, ap_dim_t dim, size_t n) {
    lsum_internal_t *pr = lsum_init_from_manager(man, AP_FUNID_EXPAND, 0);
    if (!a)
        return NULL;
    arg_assert(n <= (3+a->datadim) && dim > 0 && dim < a->segmdim, return NULL;);
    lsum_t *r;
    switch (n) {
        case 1:
            r = lsum_singleton(pr, destructive, a, dim);
            break;
        case 2:
            r = lsum_unfold(pr, destructive, a, dim);
            break;
        default:
	  r = lsum_split(pr, destructive, a, dim, n-3);
    }
    if (destructive)
        lsum_free_internal(pr, a);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_expand: code=(%zu) on\n", n);
    lsum_fdump(stdout, pr->man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return r;
}

/* Used for the merge operation:
 * the first node of tdim is the cut node, the other nodes are anonymous nodes.
 * The size of the abstract value is not changed (it will be changed by permute!).
 * Values of tdim are in [a->datadim, a->datadim+a->segmdim) except 0 to signal
 * the end of folding segments.
 * Warning: if size==1 (folding of several segments)
 *          then return the object
 */
lsum_t *
lsum_fold(ap_manager_t * man,
        bool destructive, lsum_t * a, ap_dim_t * tdim, size_t size) {
    lsum_internal_t *pr = lsum_init_from_manager(man, AP_FUNID_FOLD, 0);
    if (!a)
        return NULL;
    arg_assert(tdim && 0 < size && size < a->segmdim, return NULL;);

    /*
     * The operations done are:
     * S(tdim[0]) := S(tdim[0]) + d(tdim[1]) + S(tdim[1]) +... + d(tdim[size-1]) + S(tdim[size-1])
     * l[tdim[0]] := l[tdim[0]] + ... + l[tdim[size-1]] (not directly size because of summary nodes!)
     */
    ap_linexpr0_t *expr;
    lsum_t *r = lsum_alloc_internal(pr, a->datadim, a->segmdim);
    size_t i, basedim;

#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_fold: dimensions=[");
    for (i = 0; i < size; i++)
        fprintf(stdout, "%d, ", tdim[i]);
    fprintf(stdout, "] on ");
    lsum_fdump(stdout, man, a);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif

    /* S(tdim[0]) := S(tdim[0]) + d(tdim[1]) + S(tdim[1]) +... + d(tdim[size-1]) + S(tdim[size-1]) */
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + 3 * a->segmdim);
    basedim = a->segmdim;
    /* S(tdim[0]) := S(tdim[0]) + d(tdim[1]) + S(tdim[1]) +... + d(tdim[size-1]) + S(tdim[size-1]) */
    ap_linexpr0_set_coeff_scalar_int(expr, basedim + tdim[0], 1);
    for (i = 1; i < size; i++) {
        ap_linexpr0_set_coeff_scalar_int(expr, tdim[i], 1);
        ap_linexpr0_set_coeff_scalar_int(expr, basedim + tdim[i], 1);
    }
    r->dcons =
            ap_abstract0_assign_linexpr(pr->man_dcons, false, a->dcons,
            basedim + tdim[0], expr, NULL);
    ap_linexpr0_free(expr);
    /* l[tdim[0]] := l[tdim[0]] + ... + l[tdim[size-1]] */
    expr = ap_linexpr0_alloc(AP_LINEXPR_DENSE, a->datadim + 3 * a->segmdim);
    basedim = 2 * a->segmdim;
    for (i = 0; i < size; i++)
        ap_linexpr0_set_coeff_scalar_int(expr, basedim + tdim[i], 1);
    r->dcons =
            ap_abstract0_assign_linexpr(pr->man_dcons, true, r->dcons,
            basedim + tdim[0], expr, NULL);
    ap_linexpr0_free(expr);
    if (destructive)
        lsum_free_internal(pr, a);
#ifndef NDEBUG2
    fprintf(stdout, "\n====lsum_fold returns: ");
    lsum_fdump(stdout, man, r);
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return r;
}
