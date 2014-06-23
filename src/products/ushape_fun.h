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


#ifndef __USHAPE_FUN_H
#define __USHAPE_FUN_H


/* dependencies */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NUMFLT_PRINT_PREC ap_scalar_print_prec

#include "ap_dimension.h"
#include "ap_expr0.h"
#include "ap_manager.h"
#include "ushape.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */


  /* ============================================================ */
  /* I.1 Memory */
  /* ============================================================ */

  struct _ushape_t;
  typedef struct _ushape_t ushape_t;
  /* Abstract data type of ushapes (defined in ushape_internal.h). */

  struct _ushape_array_t;
  typedef struct _ushape_array_t ushape_array_t;
  /* Array of ushapes (defined in ushape_internal.h). */

  ushape_t *ushape_copy (ap_manager_t * man, ushape_t * a);
  /*
   * Return a copy of an abstract value, on which destructive update
   * does not affect the initial value.
   */

  void ushape_free (ap_manager_t * man, ushape_t * a);
  /* Free all the memory used by the abstract value */

  size_t ushape_size (ap_manager_t * man, ushape_t * a);
  /* Return the abstract size of an abstract value (see ap_manager_t) */

  struct _shape_internal_t;
  typedef struct _shape_internal_t ushape_internal_t;
  /* Abstract data type of library-specific manager options. */


  ushape_t *ushape_of_abstract0 (ap_abstract0_t * a);
  ap_abstract0_t *abstract0_of_ushape (ap_manager_t * man, ushape_t * oct);
  /* Wrapping / unwrapping of ushape_t in ap_abstract0_t (no copy) */


  /* ============================================================ */
  /* I.2 Control of the internal representation */
  /* ============================================================ */

  void ushape_minimize (ap_manager_t * man, ushape_t * a);
  /*
   * Minimize the size of the representation of a. This may result in a
   * later recomputation of internal information.
   */

  void ushape_canonicalize (ap_manager_t * man, ushape_t * a);
  /*
   * Put the abstract value in canonical form. (not yet clear
   * definition)
   */

  int ushape_hash (ap_manager_t * man, ushape_t * a);
  /*
   * Return an hash value for the abstract value.  Two abstract values
   * in canonical from (according to @code{ap_abstract0_canonicalize})
   * and considered as equal by the function ap_abstract0_is_eq should
   * be given the same hash value (this implies more or less a
   * canonical form).
   */

  void ushape_approximate (ap_manager_t * man, ushape_t * a, int algorithm);
  /*
   * Perform some transformation on the abstract value, guided by the
   * field algorithm.
   *
   * The transformation may lose information.  The argument "algorithm"
   * overrides the field algorithm of the structure of type ap_funopt_t
   * associated to ap_abstract0_approximate (commodity feature).
   */

  bool ushape_is_minimal (ap_manager_t * man, ushape_t * a);

  bool ushape_is_canonical (ap_manager_t * man, ushape_t * a);


  /* ============================================================ */
  /* I.3 Printing */
  /* ============================================================ */

  void ushape_fprint (FILE * stream,
                      ap_manager_t * man, ushape_t * a, char **name_of_dim);
  /*
   * Print the abstract value in a pretty way, using function
   * name_of_dim to name dimensions
   */

  void ushape_fprintdiff (FILE * stream,
                          ap_manager_t * man,
                          ushape_t * a1, ushape_t * a2, char **name_of_dim);
  /*
   * Print the difference between a1 (old value) and a2 (new value),
   * using function name_of_dim to name dimensions. The meaning of
   * difference is library dependent.
   */

  void ushape_fdump (FILE * stream, ap_manager_t * man, ushape_t * a);
  /*
   * Dump the internal representation of an abstract value, for
   * debugging purposes
   */


  /* ============================================================ */
  /* I.4 Serialization */
  /* ============================================================ */

  ap_membuf_t ushape_serialize_raw (ap_manager_t * man, ushape_t * a);
  /*
   * Allocate a memory buffer (with malloc), output the abstract value
   * in raw binary format to it and return a pointer on the memory
   * buffer and the size of bytes written.  It is the user
   * responsability to free the memory afterwards (with free).
   */

  ushape_t *ushape_deserialize_raw (ap_manager_t * man, void *ptr,
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
   * We assume that dimensions [0..intdim-1] correspond to integer
   * variables, and dimensions [intdim..intdim+realdim-1] to real
   * variables
   */

  ushape_t *ushape_bottom (ap_manager_t * man, size_t intdim, size_t realdim);
  /* Create a bottom (empty) value */

  ushape_t *ushape_top (ap_manager_t * man, size_t intdim, size_t realdim);
  /* Create a top (universe) value */

  ushape_t *ushape_of_box (ap_manager_t * man,
                           size_t intdim, size_t realdim,
                           ap_interval_t ** tinterval);
  /*
   * Abstract an hypercube defined by the array of intervals of size
   * intdim+realdim
   */

  ushape_t *ushape_of_lincons_array (ap_manager_t * man,
                                     size_t intdim, size_t realdim,
                                     ap_lincons0_array_t * array);
  /* Abstract a conjunction of linear constraints */

  ushape_t *ushape_of_tcons_array (ap_manager_t * man,
                                   size_t intdim, size_t realdim,
                                   ap_tcons0_array_t * array);
  /* Abstract a conjunction of tree expressions constraints */


  /* ============================================================ */
  /* II.2 Accessors */
  /* ============================================================ */

  ap_dimension_t ushape_dimension (ap_manager_t * man, ushape_t * a);
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

  bool ushape_is_bottom (ap_manager_t * man, ushape_t * a);
  bool ushape_is_top (ap_manager_t * man, ushape_t * a);

  bool ushape_is_leq (ap_manager_t * man, ushape_t * a1, ushape_t * a2);
  /* inclusion check */

  bool ushape_is_eq (ap_manager_t * man, ushape_t * a1, ushape_t * a2);
  /* equality check */

  bool ushape_sat_lincons (ap_manager_t * man, ushape_t * a,
                           ap_lincons0_t * cons);
  /* does the abstract value satisfy the linear constraint ? */

  bool ushape_sat_tcons (ap_manager_t * man, ushape_t * a, ap_tcons0_t * cons);
  /* does the abstract value satisfy the tree expression constraint ? */

  bool ushape_sat_interval (ap_manager_t * man, ushape_t * a,
                            ap_dim_t dim, ap_interval_t * interval);
  /* is the dimension included in the interval in the abstract value ? */

  bool ushape_is_dimension_unconstrained (ap_manager_t * man, ushape_t * a,
                                          ap_dim_t dim);
  /* is the dimension unconstrained ? */

  /* ============================================================ */
  /* II.4 Extraction of properties */
  /* ============================================================ */

  ap_interval_t *ushape_bound_linexpr (ap_manager_t * man,
                                       ushape_t * a, ap_linexpr0_t * expr);
  /*
   * Returns the interval taken by a linear expression over the
   * abstract value
   */

  ap_interval_t *ushape_bound_texpr (ap_manager_t * man,
                                     ushape_t * a, ap_texpr0_t * expr);
  /*
   * Returns the interval taken by a tree expression over the abstract
   * value
   */

  ap_interval_t *ushape_bound_dimension (ap_manager_t * man,
                                         ushape_t * a, ap_dim_t dim);
  /*
   * Returns the interval taken by the dimension over the abstract
   * value
   */

  ap_lincons0_array_t ushape_to_lincons_array (ap_manager_t * man,
                                               ushape_t * a);
  /*
   * Converts an abstract value to a polyhedra (conjunction of linear
   * constraints).
   */

  ap_tcons0_array_t ushape_to_tcons_array (ap_manager_t * man, ushape_t * a);
  /*
   * Converts an abstract value to a conjunction of tree expressions
   * constraints
   */

  ap_interval_t **ushape_to_box (ap_manager_t * man, ushape_t * a);
  /*
   * Converts an abstract value to an interval/hypercube. The size of
   * the resulting array is ushape_dimension(man,a).  This function can
   * be reimplemented by using ushape_bound_linexpr
   */

  ap_generator0_array_t ushape_to_generator_array (ap_manager_t * man,
                                                   ushape_t * a);
  /* Converts an abstract value to a system of generators. */


  /* ********************************************************************** */
  /* III. Operations */
  /* ********************************************************************** */

  /* ============================================================ */
  /* III.1 Meet and Join */
  /* ============================================================ */

  ushape_t *ushape_meet (ap_manager_t * man, bool destructive, ushape_t * a1,
                         ushape_t * a2);
  ushape_t *ushape_join (ap_manager_t * man, bool destructive, ushape_t * a1,
                         ushape_t * a2);
  /* Meet and Join of 2 abstract values */

  ushape_t *ushape_meet_array (ap_manager_t * man, ushape_t ** tab,
                               size_t size);
  ushape_t *ushape_join_array (ap_manager_t * man, ushape_t ** tab,
                               size_t size);
  /*
   * Meet and Join of an array of abstract values. Raises an
   * [[exc_invalid_argument]] exception if [[size==0]] (no way to
   * define the dimensionality of the result in such a case
   */

  ushape_t *ushape_meet_lincons_array (ap_manager_t * man,
                                       bool destructive, ushape_t * a,
                                       ap_lincons0_array_t * array);
  /*
   * Meet of an abstract value with a set of constraints (generalize
   * ushape_of_lincons_array)
   */

  ushape_t *ushape_meet_tcons_array (ap_manager_t * man,
                                     bool destructive, ushape_t * a,
                                     ap_tcons0_array_t * array);
  /*
   * Meet of an abstract value with a set of tree expressions
   * constraints.
   */

  ushape_t *ushape_add_ray_array (ap_manager_t * man,
                                  bool destructive, ushape_t * a,
                                  ap_generator0_array_t * array);
  /*
   * Generalized time elapse operator. Note: this is not like adding
   * arbitrary generators because: - ushape_add_ray_array is strict -
   * array can only contain rays and lines, not vertices
   */

  /* ============================================================ */
  /* III.2 Assignement and Substitutions */
  /* ============================================================ */

  ushape_t *ushape_assign_linexpr_array (ap_manager_t * man,
                                         bool destructive, ushape_t * a,
                                         ap_dim_t * tdim,
                                         ap_linexpr0_t ** texpr,
                                         size_t size, ushape_t * dest);
  ushape_t *ushape_substitute_linexpr_array (ap_manager_t * man,
                                             bool destructive, ushape_t * a,
                                             ap_dim_t * tdim,
                                             ap_linexpr0_t ** texpr,
                                             size_t size, ushape_t * dest);
  ushape_t *ushape_assign_texpr_array (ap_manager_t * man,
                                       bool destructive, ushape_t * a,
                                       ap_dim_t * tdim,
                                       ap_texpr0_t ** texpr,
                                       size_t size, ushape_t * dest);
  ushape_t *ushape_substitute_texpr_array (ap_manager_t * man,
                                           bool destructive, ushape_t * a,
                                           ap_dim_t * tdim,
                                           ap_texpr0_t ** texpr,
                                           size_t size, ushape_t * dest);
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

  ushape_t *ushape_forget_array (ap_manager_t * man,
                                 bool destructive, ushape_t * a,
                                 ap_dim_t * tdim, size_t size, bool project);


  /* ============================================================ */
  /* III.4 Change and permutation of dimensions */
  /* ============================================================ */

  ushape_t *ushape_add_dimensions (ap_manager_t * man,
                                   bool destructive, ushape_t * a,
                                   ap_dimchange_t * dimchange, bool project);

  ushape_t *ushape_remove_dimensions (ap_manager_t * man,
                                      bool destructive, ushape_t * a,
                                      ap_dimchange_t * dimchange);

  ushape_t *ushape_permute_dimensions (ap_manager_t * man,
                                       bool destructive,
                                       ushape_t * a,
                                       ap_dimperm_t * permutation);
  /*
   * Size of the permutation is supposed to be equal to the dimension
   * of the abstract value
   */


  /* ============================================================ */
  /* III.5 Expansion and folding of dimensions */
  /* ============================================================ */

  ushape_t *ushape_expand (ap_manager_t * man,
                           bool destructive, ushape_t * a,
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

  ushape_t *ushape_fold (ap_manager_t * man,
                         bool destructive, ushape_t * a,
                         ap_dim_t * tdim, size_t size);
  /*
   * Fold the dimensions in the array tdim of size n>=1 and put the
   * result in the first dimension in the array. The other dimensions
   * of the array are then removed.
   */


  /* ============================================================ */
  /* III.6 Widening, Narrowing */
  /* ============================================================ */

  ushape_t *ushape_widening (ap_manager_t * man, ushape_t * a1, ushape_t * a2);
  /* Standard widening: set unstable constraints to +oo */

  ushape_t *ushape_widening_thresholds (ap_manager_t * man,
                                        ushape_t * a1, ushape_t * a2,
                                        ap_scalar_t ** array, size_t nb);
  /*
   * Widening with threshold. array is assumed to contain nb
   * thresholds, sorted in increasing order.
   */

  ushape_t *ushape_narrowing (ap_manager_t * man, ushape_t * a1, ushape_t * a2);
  /* Standard narrowing: refine only +oo constraint */


  /* ============================================================ */
  /* III.7 Topological closure operation */
  /* ============================================================ */

  ushape_t *ushape_closure (ap_manager_t * man, bool destructive, ushape_t * a);
  /*
   * Returns the topological closure of a possibly opened abstract
   * value
   */


  /* ============================================================ */
  /* III.8 Additional operations */
  /* ============================================================ */

  ushape_t* ushape_of_lincons_array (ap_manager_t* man,
                                     size_t intdim, size_t realdim,
                                     ap_lincons0_array_t* array);
  ushape_t* ushape_of_tcons_array (ap_manager_t* man,
                                   size_t intdim, size_t realdim,
                                   ap_tcons0_array_t* array);
  /* Abstract a conjunction of tree/linear constraints.
   * Based on meet_*_array. */

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */


#endif /* __USHAPE_FUN_H */
