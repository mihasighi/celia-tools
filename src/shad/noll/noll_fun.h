/* SHAD - Library of shape abstract domains
 * Copyright (C) 2012-2013 LIAFA
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Please read the COPYING file packaged in the distribution.
 */

/* 
 * noll_fun.h: API for the NOLL domain
 */

#ifndef __NOLL_FUN_H
#define __NOLL_FUN_H

#include <stdlib.h>
#include <stdio.h>
#include <glib.h>

#include "sh_stmt.h"
#include "sh_manager.h"
#include "noll.h"

#ifdef __cplusplus
extern "C"
{
#endif

struct noll_val_s;
typedef struct noll_val_s noll_val_t;
/* Abstract data type of shapes (defined in noll_internal.h). */

/* ============================================================ */
/* I.1 Memory */
/* ============================================================ */

/* in noll_build.c */

noll_val_t *noll_copy(sh_manager_t * man, noll_val_t * a);
/*
 * Return a copy of an abstract value, on which destructive update
 * does not affect the initial value.
 */

void noll_free(sh_manager_t * man, noll_val_t * a);
/* Free all the memory used by the abstract value */
size_t noll_size(sh_manager_t * man, noll_val_t * a);
/* Return the abstract size, i.e., number of graphs of an
 * abstract value.
 */

/* ============================================================ */
/* I.2 Control of the internal representation */
/* ============================================================ */

/* in noll_transform.c */

void noll_minimize(sh_manager_t * man, noll_val_t * a);
/*
 * Minimize the size of the representation of a. This may reduce the
 * number of local variables used.
 */

void noll_canonicalize(sh_manager_t * man, noll_val_t * a);
/*
 * Put the abstract value in canonical form. (not yet clear
 * definition)
 */

int noll_hash(sh_manager_t * man, noll_val_t * a);
/*
 * Return an hash value for the abstract value.  Two abstract values
 * in canonical from (according to @code{ap_abstract0_canonicalize})
 * and considered as equal by the function ap_abstract0_is_eq should
 * be given the same hash value (this implies more or less a
 * canonical form).
 */

void noll_approximate(sh_manager_t * man, noll_val_t * a, int algorithm);
/*
 * Perform some transformation on the abstract value, guided by the
 * field algorithm. The transformation may lose information.
 */

gboolean noll_is_minimal(sh_manager_t * man, noll_val_t * a);
/* NOLL only: No isomorphic graphs in set,
 no graphs with anonymous nodes except nil */

gboolean noll_is_canonical(sh_manager_t * man, noll_val_t * a);
/* NOLL only: no garbage, consistent graphs */

/* ============================================================ */
/* I.3 Printing */
/* ============================================================ */

/* in noll_print.c */

void noll_fprint(FILE * stream, sh_manager_t * man, noll_val_t * a);
/* Print the abstract value in a pretty way. */

void noll_fdump(FILE * stream, sh_manager_t * man, noll_val_t * a);
/*
 * Dump the internal representation of an abstract value, for
 * debugging purposes
 */

/* ********************************************************************** */
/* II. Constructor, accessors, tests and property extraction */
/* ********************************************************************** */

/* ============================================================ */
/* II.1 Basic constructors */
/* ============================================================ */

/* in noll_build.c */

noll_val_t *noll_bottom(sh_manager_t * man, size_t fid);
/* Create a bottom (empty) value */

noll_val_t *noll_top(sh_manager_t * man, size_t fid);
/* Create a top (universe) value */

noll_val_t *noll_empty(sh_manager_t * man, size_t fid);
/* Create an value for the empty heap */

noll_val_t *noll_of_formula(sh_manager_t* man, size_t fid, void* f);
/* Build graph from NOLL formula from (temporary) code */

/* ============================================================ */
/* II.2 Getters */
/* ============================================================ */

/* in noll_query.c */
size_t noll_dimension(sh_manager_t * man, noll_val_t * a);
/* Return the frame of the abstract values */

/* ============================================================ */
/* II.3 Tests */
/* ============================================================ */

/*
 * In abstract tests,
 *
 * - true means that the predicate is certainly true.
 *
 * - false means by default don't know (an exception has occurred, or
 * the exact computation was considered too expensive to be
 * performed).
 *
 * However, if the flag exact in the manager is true, then false means
 * really that the predicate is false.
 */

/* in noll_query.c */

gboolean noll_is_bottom(sh_manager_t * man, noll_val_t * a);
gboolean noll_is_top(sh_manager_t * man, noll_val_t * a);

gboolean noll_is_le(sh_manager_t * man, noll_val_t * a1, noll_val_t * a2);
/* inclusion check */

gboolean noll_is_eq(sh_manager_t * man, noll_val_t * a1, noll_val_t * a2);
/* equality check */

gboolean noll_sat_constraint(sh_manager_t * man, noll_val_t * a, size_t x,
		size_t y, gboolean iseq);
/* does the abstract value satisfy the constraint x iseq y ? */

gboolean noll_is_dimension_unconstrained(sh_manager_t * man, noll_val_t * a,
		size_t x);
/* is the dimension x unconstrained? */

gboolean noll_is_dimension_shared(sh_manager_t * man, noll_val_t * a, size_t x);
/* is the dimension x shared? */

gboolean noll_is_dimension_nil(sh_manager_t * man, noll_val_t * a, size_t x);
/* is the dimension x at nil? */

/* ********************************************************************** */
/* III. Operations */
/* ********************************************************************** */

/* ============================================================ */
/* III.1 Meet and Join */
/* ============================================================ */

/* in noll_nary.c */

noll_val_t *noll_meet(sh_manager_t * man, gboolean destructive,
		noll_val_t * a1, noll_val_t * a2);

noll_val_t *noll_join(sh_manager_t * man, gboolean destructive,
		noll_val_t * a1, noll_val_t * a2);
/* Meet and Join of 2 abstract values */

/*
 noll_val_t *noll_meet_array (sh_manager_t * man, GPtrArray* tab);
 noll_val_t *noll_join_array (sh_manager_t * man, GPtrArray* tab);
 * Meet and Join of an array of abstract values. Raises an
 * [[exc_invalid_argument]] exception if empty tab.
 */

noll_val_t *noll_meet_cons_array(sh_manager_t * man, gboolean destructive,
		noll_val_t * a, GPtrArray* array);
/*
 * Meet of an abstract value with a set of constraints.
 */

/* ============================================================ */
/* III.2 Assignment */
/* ============================================================ */

/* in noll_transfer.c */
noll_val_t *noll_assign_post(sh_manager_t * man, gboolean destructive,
		noll_val_t * a, sh_stmt_t* stmt, noll_val_t * dest);

noll_val_t *noll_assign_pre(sh_manager_t * man, gboolean destructive,
		noll_val_t * a, sh_stmt_t* stmt, noll_val_t * dest);
/* Assignment computation, i.e., all procedure.
 Materialization is provided below, consistency check is done
 in noll_is_canonical. */

noll_val_t *noll_expand(sh_manager_t * man, gboolean destructive,
		noll_val_t * a, size_t x, size_t f);
/*
 * Materialization step for the dimension x for the field f.
 */

/* ============================================================ */
/* III.3 Change of frame */
/* ============================================================ */

/* in noll_resize.c */
noll_val_t *noll_change_dimensions(sh_manager_t * man, gboolean destructive,
		noll_val_t * a, size_t oldframe, size_t newframe, GArray substitution);
/* Change of frame using the substitution to map vars in
 * the newframe to vars in the onld frame
 */

noll_val_t *noll_project_dimensions(sh_manager_t * man, gboolean destructive,
		noll_val_t * a, GArray dimensions);
/* Compute the graphs with only dimensions not at nil */

/* ============================================================ */
/* III.4 Widening */
/* ============================================================ */

/* in noll_widen.c */

noll_val_t *noll_widening(sh_manager_t * man, noll_val_t * a1, noll_val_t * a2);
noll_val_t * noll_fold (sh_manager_t * man, gboolean destructive,
		noll_val_t * a, size_t annon, GArray *preds);

/* Standard widening: set unstable constraints to +oo */

/* ============================================================ */
/* III.5 Closure operation */
/* ============================================================ */

/* in noll_closure.c */

noll_val_t *noll_closure(sh_manager_t * man, gboolean destructive,
		noll_val_t * a);
/*
 * Get minimal value wrt comparison.
 */

#ifdef __cplusplus
}
#endif

#endif /* __NOLL_FUN_H */
