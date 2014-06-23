/*
 * ucons_fun.h
 *
 * Direct access to high level heap graphs functions bypassing the manager.
 * The only file you need to include to use constrained uconss.
 * The complete semantics of these functions is given in ap_abstract0.h
 *
 * APRON Library / Shape Domain
 *
 * Copyright (C) LIAFA 2009
 *
 */

/*
 * This file is part of the APRON Library, released under LGPL license.
 * Please read the COPYING file packaged in the distribution.
 */

#ifndef __UCONS_FUN_H
#define __UCONS_FUN_H


/* dependencies */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define NUMFLT_PRINT_PREC ap_scalar_print_prec

#include "ap_dimension.h"
#include "ap_expr0.h"
#include "ap_manager.h"
#include "ucons.h"


/* *INDENT-OFF* */
#ifdef __cplusplus
extern          "C"
{
#endif
  /* *INDENT-ON* */



  /* ============================================================ */
  /* I.1 Memory */
  /* ============================================================ */

struct _ucons_t;
typedef struct _ucons_t ucons_t;
  /* Abstract data type of uconss (defined in ucons_internal.h). */

struct _ucons_internal_t;
typedef struct _ucons_internal_t ucons_internal_t;
  /* Abstract data type of library-specific manager options. */

ucons_t *ucons_copy (ap_manager_t * man, ucons_t * a);
  /*
   * Return a copy of an abstract value, on which destructive update
   * does not affect the initial value.
   */

void ucons_free (ap_manager_t * man, ucons_t * a);
  /* Free all the memory used by the abstract value */

size_t ucons_size (ap_manager_t * man, ucons_t * a);
  /* Return the abstract size of an abstract value (see ap_manager_t) */


ucons_t *ucons_of_abstract0 (ap_abstract0_t * a);
ap_abstract0_t *abstract0_of_ucons (ap_manager_t * man, ucons_t * oct);
  /* Wrapping / unwrapping of ucons_t in ap_abstract0_t (no copy) */


  /* ============================================================ */
  /* I.2 Control of the internal representation */
  /* ============================================================ */

void ucons_minimize (ap_manager_t * man, ucons_t * a);
  /*
   * Minimize the size of the representation of a. This may result in a
   * later recomputation of internal information.
   */

void ucons_canonicalize (ap_manager_t * man, ucons_t * a);
  /*
   * Put the abstract value in canonical form. (not yet clear
   * definition)
   */

int ucons_hash (ap_manager_t * man, ucons_t * a);
  /*
   * Return an hash value for the abstract value.  Two abstract values
   * in canonical from (according to @code{ap_abstract0_canonicalize})
   * and considered as equal by the function ap_abstract0_is_eq should
   * be given the same hash value (this implies more or less a
   * canonical form).
   */

void ucons_approximate (ap_manager_t * man, ucons_t * a, int algorithm);
  /*
   * Perform some transformation on the abstract value, guided by the
   * field algorithm.
   *
   * The transformation may lose information.  The argument "algorithm"
   * overrides the field algorithm of the structure of type ap_funopt_t
   * associated to ap_abstract0_approximate (commodity feature).
   */

bool ucons_is_minimal (ap_manager_t * man, ucons_t * a);

bool ucons_is_canonical (ap_manager_t * man, ucons_t * a);


  /* ============================================================ */
  /* I.3 Printing */
  /* ============================================================ */

void ucons_fprint (FILE * stream,
		   ap_manager_t * man, ucons_t * a, char **name_of_dim);
  /*
   * Print the abstract value in a pretty way, using function
   * name_of_dim to name dimensions
   */

void ucons_fprintdiff (FILE * stream,
		       ap_manager_t * man,
		       ucons_t * a1, ucons_t * a2, char **name_of_dim);
  /*
   * Print the difference between a1 (old value) and a2 (new value),
   * using function name_of_dim to name dimensions. The meaning of
   * difference is library dependent.
   */

void ucons_fdump (FILE * stream, ap_manager_t * man, ucons_t * a);
  /*
   * Dump the internal representation of an abstract value, for
   * debugging purposes
   */

//void
//pattern_key_fprint (FILE * stream, ucons_internal_t *pr, pattern_key_t * a,
//              char **name_of_dim);

  /* ============================================================ */
  /* I.4 Serialization */
  /* ============================================================ */

ap_membuf_t ucons_serialize_raw (ap_manager_t * man, ucons_t * a);
  /*
   * Allocate a memory buffer (with malloc), output the abstract value
   * in raw binary format to it and return a pointer on the memory
   * buffer and the size of bytes written.  It is the user
   * responsability to free the memory afterwards (with free).
   */

ucons_t *ucons_deserialize_raw (ap_manager_t * man, void *ptr, size_t * size);
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

ucons_t *ucons_bottom (ap_manager_t * man, size_t intdim, size_t realdim);
  /* Create a bottom (empty) value */

ucons_t *ucons_top (ap_manager_t * man, size_t intdim, size_t realdim);
  /* Create a top (universe) value */

ucons_t *ucons_of_box (ap_manager_t * man,
		       size_t intdim, size_t realdim,
		       ap_interval_t ** tinterval);
  /*
   * Abstract an hypercube defined by the array of intervals of size
   * intdim+realdim
   */

ucons_t *ucons_of_lincons_array (ap_manager_t * man,
				 size_t intdim, size_t realdim,
				 ap_lincons0_array_t * array);
  /* Abstract a conjunction of linear constraints */

ucons_t *ucons_of_tcons_array (ap_manager_t * man,
			       size_t intdim, size_t realdim,
			       ap_tcons0_array_t * array);
  /* Abstract a conjunction of tree expressions constraints */


  /* ============================================================ */
  /* II.2 Accessors */
  /* ============================================================ */

ap_dimension_t ucons_dimension (ap_manager_t * man, ucons_t * a);
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

bool ucons_is_bottom (ap_manager_t * man, ucons_t * a);
bool ucons_is_top (ap_manager_t * man, ucons_t * a);

bool ucons_is_leq (ap_manager_t * man, ucons_t * a1, ucons_t * a2);
  /* inclusion check */

bool ucons_is_eq (ap_manager_t * man, ucons_t * a1, ucons_t * a2);
  /* equality check */

bool ucons_sat_lincons (ap_manager_t * man, ucons_t * a,
			ap_lincons0_t * cons);
  /* does the abstract value satisfy the linear constraint ? */

bool ucons_sat_tcons (ap_manager_t * man, ucons_t * a, ap_tcons0_t * cons);
  /* does the abstract value satisfy the tree expression constraint ? */

bool ucons_sat_interval (ap_manager_t * man, ucons_t * a,
			 ap_dim_t dim, ap_interval_t * interval);
  /* is the dimension included in the interval in the abstract value ? */

bool ucons_is_dimension_unconstrained (ap_manager_t * man, ucons_t * a,
				       ap_dim_t dim);
  /* is the dimension unconstrained ? */

  /* ============================================================ */
  /* II.4 Extraction of properties */
  /* ============================================================ */

ap_interval_t *ucons_bound_linexpr (ap_manager_t * man,
				    ucons_t * a, ap_linexpr0_t * expr);
  /*
   * Returns the interval taken by a linear expression over the
   * abstract value
   */

ap_interval_t *ucons_bound_texpr (ap_manager_t * man,
				  ucons_t * a, ap_texpr0_t * expr);
  /*
   * Returns the interval taken by a tree expression over the abstract
   * value
   */

ap_interval_t *ucons_bound_dimension (ap_manager_t * man,
				      ucons_t * a, ap_dim_t dim);
  /*
   * Returns the interval taken by the dimension over the abstract
   * value
   */

ap_lincons0_array_t ucons_to_lincons_array (ap_manager_t * man, ucons_t * a);
  /*
   * Converts an abstract value to a polyhedra (conjunction of linear
   * constraints).
   */

ap_tcons0_array_t ucons_to_tcons_array (ap_manager_t * man, ucons_t * a);
  /*
   * Converts an abstract value to a conjunction of tree expressions
   * constraints
   */

ap_interval_t **ucons_to_box (ap_manager_t * man, ucons_t * a);
  /*
   * Converts an abstract value to an interval/hypercube. The size of
   * the resulting array is ucons_dimension(man,a).  This function can
   * be reimplemented by using ucons_bound_linexpr
   */

ap_generator0_array_t ucons_to_generator_array (ap_manager_t * man,
						ucons_t * a);
  /* Converts an abstract value to a system of generators. */


  /* ********************************************************************** */
  /* III. Operations */
  /* ********************************************************************** */

  /* ============================================================ */
  /* III.1 Meet and Join */
  /* ============================================================ */

ucons_t *ucons_meet (ap_manager_t * man, bool destructive, ucons_t * a1,
		     ucons_t * a2);
ucons_t *ucons_join (ap_manager_t * man, bool destructive, ucons_t * a1,
		     ucons_t * a2);
  /* Meet and Join of 2 abstract values */

ucons_t *ucons_meet_array (ap_manager_t * man, ucons_t ** tab, size_t size);
ucons_t *ucons_join_array (ap_manager_t * man, ucons_t ** tab, size_t size);
  /*
   * Meet and Join of an array of abstract values. Raises an
   * [[exc_invalid_argument]] exception if [[size==0]] (no way to
   * define the dimensionality of the result in such a case
   */

ucons_t *ucons_meet_lincons_array (ap_manager_t * man,
				   bool destructive, ucons_t * a,
				   ap_lincons0_array_t * array);
  /*
   * Meet of an abstract value with a set of constraints (generalize
   * ucons_of_lincons_array)
   */

ucons_t *ucons_meet_tcons_array (ap_manager_t * man,
				 bool destructive, ucons_t * a,
				 ap_tcons0_array_t * array);
  /*
   * Meet of an abstract value with a set of tree expressions
   * constraints. (generalize ucons_of_tcons_array)
   */

ucons_t *ucons_add_ray_array (ap_manager_t * man,
			      bool destructive, ucons_t * a,
			      ap_generator0_array_t * array);
  /*
   * Generalized time elapse operator. Note: this is not like adding
   * arbitrary generators because: - ucons_add_ray_array is strict -
   * array can only contain rays and lines, not vertices
   */

  /* ============================================================ */
  /* III.2 Assignement and Substitutions */
  /* ============================================================ */

ucons_t *ucons_assign_linexpr_array (ap_manager_t * man,
				     bool destructive, ucons_t * a,
				     ap_dim_t * tdim,
				     ap_linexpr0_t ** texpr,
				     size_t size, ucons_t * dest);
ucons_t *ucons_substitute_linexpr_array (ap_manager_t * man,
					 bool destructive, ucons_t * a,
					 ap_dim_t * tdim,
					 ap_linexpr0_t ** texpr,
					 size_t size, ucons_t * dest);
ucons_t *ucons_assign_texpr_array (ap_manager_t * man,
				   bool destructive, ucons_t * a,
				   ap_dim_t * tdim,
				   ap_texpr0_t ** texpr,
				   size_t size, ucons_t * dest);
ucons_t *ucons_substitute_texpr_array (ap_manager_t * man,
				       bool destructive, ucons_t * a,
				       ap_dim_t * tdim,
				       ap_texpr0_t ** texpr,
				       size_t size, ucons_t * dest);
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

ucons_t *ucons_forget_array (ap_manager_t * man,
			     bool destructive, ucons_t * a,
			     ap_dim_t * tdim, size_t size, bool project);


  /* ============================================================ */
  /* III.4 Change and permutation of dimensions */
  /* ============================================================ */


ucons_t *ucons_add_dimensions (ap_manager_t * man,
			       bool destructive, ucons_t * a,
			       ap_dimchange_t * dimchange, bool project);

ucons_t *add_pattern_1 (ucons_internal_t * pr, ucons_t * r, size_t dim);

ucons_t * add_pattern_1_lx( ucons_internal_t *pr, ucons_t *r, size_t dim);

ucons_t *add_pattern_2_1 (ucons_internal_t * pr, ucons_t * r, size_t dim);

ucons_t *add_pattern_1_2 (ucons_internal_t * pr, ucons_t * r, size_t dim);

ucons_t *add_pattern_3_1 (ucons_internal_t * pr, ucons_t * r, size_t dim);

ucons_t *ucons_remove_dimensions (ap_manager_t * man,
				  bool destructive, ucons_t * a,
				  ap_dimchange_t * dimchange);

ucons_t *ucons_permute_dimensions (ap_manager_t * man,
				   bool destructive,
				   ucons_t * a, ap_dimperm_t * permutation);
  /*
   * Size of the permutation is supposed to be equal to the dimension
   * of the abstract value
   */


  /* ============================================================ */
  /* III.5 Expansion and folding of dimensions */
  /* ============================================================ */

ucons_t *ucons_expand (ap_manager_t * man,
		       bool destructive, ucons_t * a, ap_dim_t dim, size_t n);
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

ucons_t *ucons_fold (ap_manager_t * man,
		     bool destructive, ucons_t * a,
		     ap_dim_t * tdim, size_t size);
  /*
   * Fold the dimensions in the array tdim of size n>=1 and put the
   * result in the first dimension in the array. The other dimensions
   * of the array are then removed.
   */

/* is called when foding with pattern P12 while unrolling only one segment */
/* called from ucons_fold then minus_dim <=0 and update_length = true
 *  called from fold_without_closoure_of_P21 then minus_dim >0 and update_length = false
 * */
ucons_t * fold_with_closoure_of_P21(ucons_internal_t * pr, ucons_t * r, ap_dim_t * sdim,
		size_t size, ap_dim_t minus_dim, bool update_length);

ucons_t * fold_with_closure_P11_or_P12(ucons_internal_t * pr,
		ucons_t * r, ap_dim_t * tdim, size_t size,bool update_lenght);

ucons_t * fold_with_closure_succ_P11(ucons_internal_t * pr,
		ucons_t * r, ap_dim_t * tdim, size_t size,bool update_lenght);

ucons_t * fold_without_closoure_of_P21(ucons_internal_t * pr, ucons_t * a, ap_dim_t * node_tdim,
		size_t size_node_tdim,  bool update_length);

//ucons_t * merge_succesors_P21(ucons_internal_t *pr, ucons_t * r, ap_dim_t * gdim, size_t size_gdim);
ucons_t *merge_succesors_1 (ucons_internal_t * pr, ucons_t * r,
			    ap_dim_t * gdim, size_t size_gdim,bool update_lenght);

ucons_t *merge_succesors_2 (ucons_internal_t * pr, ucons_t * r,
			    ap_dim_t * gdim, size_t size_gdim,bool update_lenght);

ucons_t *generate_pattern1_from_succesors (ucons_internal_t * pr, ucons_t * r,
					   ap_dim_t * gdim, size_t size_gdim);

ucons_t *concat_nodes_1 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
			 size_t size_cdim,bool update_lenght);

//ucons_t *concat_nodes_P21 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
//			 size_t size_cdim);

void update_lenghts (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
			 size_t size_cdim);

ap_abstract0_t *concat_P11 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
			 size_t size_cdim);

ap_abstract0_t *concat_P12 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
			 size_t size_cdim);

ucons_t *concat_nodes_2 (ucons_internal_t * pr, ucons_t * r, ap_dim_t * cdim,
			 size_t size_cdim,bool update_lenght);

/* (y1,y2) = (tdim[i],y=1 ) */
ap_abstract0_t * create4_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i);

/* (y1,y2) _substitute of succ_P12  */
ap_abstract0_t * create3_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i);

/*(y1,y2) == (last(node_tdim[i-1]), node_tdim[i])*/
ap_abstract0_t * create2_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i);

/*(y1,y2) == (node_tdim[i-1], node_tdim[i])*/
ap_abstract0_t * create1_succ_P12( ucons_internal_t *pr,
		ucons_t * r, ap_dim_t* node_tdim, size_t node_tdim_size, size_t i);



ap_lincons0_t cons_dy0_dy1(size_t intdim, size_t realdim, ap_dim_t dy0, ap_dim_t dy1);
ap_lincons0_t cons_dy0(size_t intdim, size_t realdim, ap_dim_t dy0, ap_dim_t n);
ap_lincons0_t cons_k_m(size_t intdim, size_t realdim, ap_dim_t k, ap_dim_t m);

  /* ============================================================ */
  /* III.6 Widening, Narrowing */
  /* ============================================================ */

ucons_t *ucons_widening (ap_manager_t * man, ucons_t * a1, ucons_t * a2);
  /* Standard widening: set unstable constraints to +oo */

ucons_t *ucons_widening_thresholds (ap_manager_t * man,
				    ucons_t * a1, ucons_t * a2,
				    ap_scalar_t ** array, size_t nb);
  /*
   * Widening with threshold. array is assumed to contain nb
   * thresholds, sorted in increasing order.
   */

ucons_t *ucons_narrowing (ap_manager_t * man, ucons_t * a1, ucons_t * a2);
  /* Standard narrowing: refine only +oo constraint */


  /* ============================================================ */
  /* III.7 Topological closure operation */
  /* ============================================================ */
void close_eq(size_t ** tdim, size_t size);

/* intersection between universal constraints and multiset constraints */
ucons_t*
ucons_strengthen(ucons_internal_t* pr, ucons_t* a,size_t ** eqdim);

/* saturate using projections */
ucons_t*
ucons_saturation(ucons_internal_t* pr, ucons_t* a);

ucons_t *ucons_closure (ap_manager_t * man, bool destructive, ucons_t * a);
  /*
   * Returns the topological closure of a possibly opened abstract
   * value
  */
bool test_l_leq_v(ap_manager_t * man, ap_abstract0_t * dcons,
		size_t intdim, size_t realdim, ap_dim_t n, int v);
bool test_g_leq_v(ap_manager_t * man, ap_abstract0_t * dcons,
		size_t intdim, size_t realdim, ap_dim_t n, int v);

bool test_singleton(ap_manager_t * man, ap_abstract0_t * dcons, size_t intdim,
		size_t realdim, ap_dim_t n);

bool test_equal_length(ap_manager_t * man, ap_abstract0_t * dcons, size_t intdim,
		size_t realdim, ap_dim_t n, ap_dim_t m);

ap_dim_t * equal_lengths(ap_manager_t * man, ap_abstract0_t * dcons, size_t intdim,
		size_t realdim, ap_dim_t n, ap_dim_t *tdim, size_t size_tdim, size_t* size);

ap_lincons0_t cons_dy0_dy1(size_t intdim, size_t realdim, ap_dim_t dy0, ap_dim_t dy1);
ap_lincons0_t cons_dy0(size_t intdim, size_t realdim, ap_dim_t dy0, ap_dim_t n);
ap_lincons0_t cons_k_m(size_t intdim, size_t realdim, ap_dim_t k, ap_dim_t m);



/* r = meet(a, with EQMOD constraints on useg */
void ucons_meet_eq_structural_constraint(ucons_internal_t *pr, ucons_t *a,
		ucons_t *r, ap_dim_t *useg, size_t size_useg );


  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */

#endif /* __SHAPE_FUN_H */
