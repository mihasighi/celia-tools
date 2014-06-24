/**************************************************************************/
/*                                                                        */
/*  CELIA Tools / SLL Abstract Domain                                     */
/*                                                                        */
/*  Copyright (C) 2009-2014                                               */
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


#ifndef HGRAPH_FUN_H_
#define HGRAPH_FUN_H_

/*
 * API of HGRAPH as absract domain.
 * Implements the APRON API, for the moment; not all methods are useful.
 */

/* dependencies */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NUMFLT_PRINT_PREC ap_scalar_print_prec

/* for APRON API only */
#include "ap_coeff.h"
#include "ap_dimension.h"
#include "ap_expr0.h"

/* for SHAD API */
#include "hgraph.h"
#include "shad.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */

  /* ============================================================ */
  /* I.1 Memory */
  /* ============================================================ */

struct _hgraph_t;
typedef struct _hgraph_t hgraph_t;
  /* Abstract data type of hgraphs (defined in hgraph_internal.h). */

hgraph_t *hgraph_copy (sh_manager_t * man, hgraph_t * a);
  /*
   * Return a copy of an abstract value, on which destructive update
   * does not affect the initial value.
   */

void hgraph_free (sh_manager_t * man, hgraph_t * a);
  /* Free all the memory used by the abstract value */

size_t hgraph_size (sh_manager_t * man, hgraph_t * a);
  /* Return the abstract size of an abstract value (see sh_manager_t) */

struct _shape_internal_t;
typedef struct _shape_internal_t hgraph_internal_t;
  /* Abstract data type of library-specific manager options. */


hgraph_t *hgraph_of_abstract0 (sh_abstract0_t * a);
sh_abstract0_t *abstract0_of_hgraph (sh_manager_t * man, hgraph_t * oct);
  /* Wrapping / unwrapping of hgraph_t in sh_abstract0_t (no copy) */


  /* ============================================================ */
  /* I.2 Control of the internal representation */
  /* ============================================================ */

void hgraph_minimize (sh_manager_t * man, hgraph_t * a);
  /*
   * Minimize the size of the representation of a. This may result in a
   * later recomputation of internal information.
   */

void hgraph_canonicalize (sh_manager_t * man, hgraph_t * a);
  /*
   * Put the abstract value in canonical form. (not yet clear
   * definition)
   */

int hgraph_hash (sh_manager_t * man, hgraph_t * a);
  /*
   * Return an hash value for the abstract value.  Two abstract values
   * in canonical from (according to @code{ap_abstract0_canonicalize})
   * and considered as equal by the function ap_abstract0_is_eq should
   * be given the same hash value (this implies more or less a
   * canonical form).
   */

void hgraph_approximate (sh_manager_t * man, hgraph_t * a, int algorithm);
  /*
   * Perform some transformation on the abstract value, guided by the
   * field algorithm.
   *
   * The transformation may lose information.  The argument "algorithm"
   * overrides the field algorithm of the structure of type ap_funopt_t
   * associated to ap_abstract0_approximate (commodity feature).
   */

bool hgraph_is_minimal (sh_manager_t * man, hgraph_t * a);

bool hgraph_is_canonical (sh_manager_t * man, hgraph_t * a);


  /* ============================================================ */
  /* I.3 Printing */
  /* ============================================================ */

void hgraph_fprint (FILE * stream,
                    sh_manager_t * man, hgraph_t * a, char **name_of_dim);
  /*
   * Print the abstract value in a pretty way, using function
   * name_of_dim to name dimensions
   */

void hgraph_fprintdiff (FILE * stream,
                        sh_manager_t * man,
                        hgraph_t * a1, hgraph_t * a2, char **name_of_dim);
  /*
   * Print the difference between a1 (old value) and a2 (new value),
   * using function name_of_dim to name dimensions. The meaning of
   * difference is library dependent.
   */

void hgraph_fdump (FILE * stream, sh_manager_t * man, hgraph_t * a);
  /*
   * Dump the internal representation of an abstract value, for
   * debugging purposes
   */


  /* ============================================================ */
  /* I.4 Serialization */
  /* ============================================================ */

ap_membuf_t hgraph_serialize_raw (sh_manager_t * man, hgraph_t * a);
  /*
   * Allocate a memory buffer (with malloc), output the abstract value
   * in raw binary format to it and return a pointer on the memory
   * buffer and the size of bytes written.  It is the user
   * responsability to free the memory afterwards (with free).
   */

hgraph_t *hgraph_deserialize_raw (sh_manager_t * man, void *ptr,
                                  size_t * size);
  /*
   * Return the abstract value read in raw binary format from the input
   * stream and store in size the number of bytes read
   */


  /* ********************************************************************** */
  /* II. Constructor, accessors, tests and property extraction */
  /* ********************************************************************** */

  /* ============================================================ */
  /* II.1 Basic constructors */
  /* ============================================================ */

  /*
   * We assume that dimensions [0..intdim-1] correspond to data
   * variables, and dimensions [intdim..intdim+realdim-1] to ptr
   * variables
   */

hgraph_t *hgraph_bottom (sh_manager_t * man, size_t intdim, size_t realdim);
  /* Create a bottom (empty) value */

hgraph_t *hgraph_top (sh_manager_t * man, size_t intdim, size_t realdim);
  /* Create a top (universe) value */

hgraph_t *hgraph_of_box (sh_manager_t * man,
                         size_t intdim, size_t realdim,
                         ap_interval_t ** tinterval);
  /*
   * Abstract an hypercube defined by the array of intervals of size
   * intdim+realdim
   */

hgraph_t *hgraph_of_lincons_array (sh_manager_t * man,
                                   size_t intdim, size_t realdim,
                                   ap_lincons0_array_t * array);
  /* Abstract a conjunction of linear constraints */

hgraph_t *hgraph_of_tcons_array (sh_manager_t * man,
                                 size_t intdim, size_t realdim,
                                 ap_tcons0_array_t * array);
  /* Abstract a conjunction of tree expressions constraints */


hgraph_t *hgraph_of_formula (hgraph_internal_t * pr,
                             sh_formula_t * f, size_t disj,
                             ap_dimperm_t * perm);
  /* Abstract the graph of the shadformula and return the permutation */


  /* ============================================================ */
  /* II.2 Accessors */
  /* ============================================================ */

ap_dimension_t hgraph_dimension (sh_manager_t * man, hgraph_t * a);
  /* Return the total number of dimensions of the abstract values */


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

bool hgraph_is_bottom (sh_manager_t * man, hgraph_t * a);
bool hgraph_is_top (sh_manager_t * man, hgraph_t * a);

bool hgraph_is_leq (sh_manager_t * man, hgraph_t * a1, hgraph_t * a2);
  /* inclusion check */

bool hgraph_is_eq (sh_manager_t * man, hgraph_t * a1, hgraph_t * a2);
  /* equality check */

bool hgraph_sat_lincons (sh_manager_t * man, hgraph_t * a,
                         ap_lincons0_t * cons);
  /* does the abstract value satisfy the linear constraint ? */

bool hgraph_sat_tcons (sh_manager_t * man, hgraph_t * a, ap_tcons0_t * cons);
  /* does the abstract value satisfy the tree expression constraint ? */

bool hgraph_sat_interval (sh_manager_t * man, hgraph_t * a,
                          ap_dim_t dim, ap_interval_t * interval);
  /* is the dimension included in the interval in the abstract value ? */

bool hgraph_is_dimension_unconstrained (sh_manager_t * man, hgraph_t * a,
                                        ap_dim_t dim);
  /* is the dimension unconstrained ? */

  /* ============================================================ */
  /* II.4 Extraction of properties */
  /* ============================================================ */

ap_interval_t *hgraph_bound_linexpr (sh_manager_t * man,
                                     hgraph_t * a, ap_linexpr0_t * expr);
  /*
   * Returns the interval taken by a linear expression over the
   * abstract value
   */

ap_interval_t *hgraph_bound_texpr (sh_manager_t * man,
                                   hgraph_t * a, ap_texpr0_t * expr);
  /*
   * Returns the interval taken by a tree expression over the abstract
   * value
   */

ap_interval_t *hgraph_bound_dimension (sh_manager_t * man,
                                       hgraph_t * a, ap_dim_t dim);
  /*
   * Returns the interval taken by the dimension over the abstract
   * value
   */

ap_lincons0_array_t hgraph_to_lincons_array (sh_manager_t * man,
                                             hgraph_t * a);
  /*
   * Converts an abstract value to a polyhedra (conjunction of linear
   * constraints).
   */

ap_tcons0_array_t hgraph_to_tcons_array (sh_manager_t * man, hgraph_t * a);
  /*
   * Converts an abstract value to a conjunction of tree expressions
   * constraints
   */

ap_interval_t **hgraph_to_box (sh_manager_t * man, hgraph_t * a);
  /*
   * Converts an abstract value to an interval/hypercube. The size of
   * the resulting array is hgraph_dimension(man,a).  This function can
   * be reimplemented by using hgraph_bound_linexpr
   */

ap_generator0_array_t hgraph_to_generator_array (sh_manager_t * man,
                                                 hgraph_t * a);
  /* Converts an abstract value to a system of generators. */


  /* ********************************************************************** */
  /* III. Operations */
  /* ********************************************************************** */

  /* ============================================================ */
  /* III.1 Meet and Join */
  /* ============================================================ */

hgraph_t *hgraph_meet (sh_manager_t * man, bool destructive, hgraph_t * a1,
                       hgraph_t * a2);
hgraph_t *hgraph_join (sh_manager_t * man, bool destructive, hgraph_t * a1,
                       hgraph_t * a2);
  /* Meet and Join of 2 abstract values */

hgraph_t *hgraph_meet_array (sh_manager_t * man, hgraph_t ** tab,
                             size_t size);
hgraph_t *hgraph_join_array (sh_manager_t * man, hgraph_t ** tab,
                             size_t size);
  /*
   * Meet and Join of an array of abstract values. Raises an
   * [[exc_invalid_argument]] exception if [[size==0]] (no way to
   * define the dimensionality of the result in such a case
   */

hgraph_t *hgraph_meet_lincons_array (sh_manager_t * man,
                                     bool destructive, hgraph_t * a,
                                     ap_lincons0_array_t * array);
  /*
   * Meet of an abstract value with a set of constraints (generalize
   * hgraph_of_lincons_array)
   */

hgraph_t *hgraph_meet_tcons_array (sh_manager_t * man,
                                   bool destructive, hgraph_t * a,
                                   ap_tcons0_array_t * array);
  /*
   * Meet of an abstract value with a set of tree expressions
   * constraints. (generalize hgraph_of_tcons_array)
   */

hgraph_t *hgraph_add_ray_array (sh_manager_t * man,
                                bool destructive, hgraph_t * a,
                                ap_generator0_array_t * array);
  /*
   * Generalized time elapse operator. Note: this is not like adding
   * arbitrary generators because: - hgraph_add_ray_array is strict -
   * array can only contain rays and lines, not vertices
   */

  /* ============================================================ */
  /* III.2 Assignement and Substitutions */
  /* ============================================================ */

hgraph_t *hgraph_assign_linexpr_array (sh_manager_t * man,
                                       bool destructive, hgraph_t * a,
                                       ap_dim_t * tdim,
                                       ap_linexpr0_t ** texpr,
                                       size_t size, hgraph_t * dest);
hgraph_t *hgraph_substitute_linexpr_array (sh_manager_t * man,
                                           bool destructive, hgraph_t * a,
                                           ap_dim_t * tdim,
                                           ap_linexpr0_t ** texpr,
                                           size_t size, hgraph_t * dest);
hgraph_t *hgraph_assign_texpr_array (sh_manager_t * man,
                                     bool destructive, hgraph_t * a,
                                     ap_dim_t * tdim,
                                     ap_texpr0_t ** texpr,
                                     size_t size, hgraph_t * dest);
hgraph_t *hgraph_substitute_texpr_array (sh_manager_t * man,
                                         bool destructive, hgraph_t * a,
                                         ap_dim_t * tdim,
                                         ap_texpr0_t ** texpr,
                                         size_t size, hgraph_t * dest);
  /*
   * Parallel Assignement and Substitution of several dimensions by
   * expressions in abstract value a.
   *
   * dest is an optional argument. If not NULL, semantically speaking, the
   * result of the transformation is intersected with dest. This is
   * useful for precise backward transformations in lattices like
   * intervals or octagons.
   */


  /* ============================================================ */
  /* III.3 Projections */
  /* ============================================================ */

hgraph_t *hgraph_forget_array (sh_manager_t * man,
                               bool destructive, hgraph_t * a,
                               ap_dim_t * tdim, size_t size, bool project);


  /* ============================================================ */
  /* III.4 Change and permutation of dimensions */
  /* ============================================================ */

hgraph_t *hgraph_add_dimensions (sh_manager_t * man,
                                 bool destructive, hgraph_t * a,
                                 ap_dimchange_t * dimchange, bool project);

hgraph_t *hgraph_remove_dimensions (sh_manager_t * man,
                                    bool destructive, hgraph_t * a,
                                    ap_dimchange_t * dimchange);

hgraph_t *hgraph_permute_dimensions (sh_manager_t * man,
                                     bool destructive,
                                     hgraph_t * a,
                                     ap_dimperm_t * permutation);
  /*
   * Size of the permutation is supposed to be equal to the dimension
   * of the abstract value
   */


  /* ============================================================ */
  /* III.5 Expansion and folding of dimensions */
  /* ============================================================ */

hgraph_t *hgraph_expand (sh_manager_t * man,
                         bool destructive, hgraph_t * a,
                         ap_dim_t dim, size_t n);
  /*
   * Expand the dimension dim into itself + n additional dimensions. It
   * results in (n+1) unrelated dimensions having same relations with
   * other dimensions. The (n+1) dimensions are put as follows:
   *
   * - original dimension dim
   *
   * - if the dimension is integer, the n additional dimensions are put at
   * the end of integer dimensions; if it is real, at the end of the
   * real dimensions.
   */

hgraph_t *hgraph_fold (sh_manager_t * man,
                       bool destructive, hgraph_t * a,
                       ap_dim_t * tdim, size_t size);
  /*
   * Fold the dimensions in the array tdim of size n>=1 and put the
   * result in the first dimension in the array. The other dimensions
   * of the array are then removed.
   */


  /* ============================================================ */
  /* III.6 Widening, Narrowing */
  /* ============================================================ */

hgraph_t *hgraph_widening (sh_manager_t * man, hgraph_t * a1, hgraph_t * a2);
  /* Standard widening: set unstable constraints to +oo */

hgraph_t *hgraph_widening_thresholds (sh_manager_t * man,
                                      hgraph_t * a1, hgraph_t * a2,
                                      ap_scalar_t ** array, size_t nb);
  /*
   * Widening with threshold. array is assumed to contain nb
   * thresholds, sorted in increasing order.
   */

hgraph_t *hgraph_narrowing (sh_manager_t * man, hgraph_t * a1, hgraph_t * a2);
  /* Standard narrowing: refine only +oo constraint */


  /* ============================================================ */
  /* III.7 Topological closure operation */
  /* ============================================================ */

hgraph_t *hgraph_closure (sh_manager_t * man, bool destructive, hgraph_t * a);
  /*
   * Returns the topological closure of a possibly opened abstract
   * value
   */


  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* HGRAPH_FUN_H_ */
