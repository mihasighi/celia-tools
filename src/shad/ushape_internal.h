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


#ifndef __USHAPE_INTERNAL_H_
#define __USHAPE_INTERNAL_H_

#include "hgraph_fun.h"
#include "hgraph_internal.h"
#include "ushape_fun.h"
#include "shape_manager.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */


  /* ********************************************************************** */
  /* Constrained shapes */
  /* ********************************************************************** */

  /* ============================================================ */
  /* Internal Representation */
  /* ============================================================ */

  /** Ushapes are built from a heap graph and an array (of size
   * given by the manager) of ``segment'' abstract domains.
   * These segment abstract domains are speaking about
   * the nodes and the segments of the heap wrt to
   * some abstraction of theirs property (see ucons and lcons domains).
   */

  struct _ushape_t
  {
    hgraph_t *h; /* heap graph (or NULL) */
    hgraph_t *closed; /* closed version of h (or NULL if not
				 * available) */

    size_t datadim; /* the number of length and data prg
				 * variables */
    size_t ptrdim; /* the number of ptr prg variables */

    void **scons; /* Array of node and segment constraints
				 * given by different abstract domains for
				 * each node/segment of the heap-graph. The
				 * size of the array is fixed by the manager
				 * field scons_size. The dimension of each
				 * domain is intdim=datadim+2*h->size, realdim=h->size */
  };

  /** Array of ushapes */
  struct _ushape_array_t
  {
    size_t size;
    ushape_t **p;
  };


#define ushape_init_from_manager(man,id,size)  ((ushape_internal_t*) shape_init_from_manager(man,id,size))

  /* Translation between dimensions */
#define NODE_DATADIM(datadim,size,node) (datadim + node)
#define NODE_LENDIM(datadim,size,node) (datadim + size + node)

  /* ============================================================ */
  /* Internal Management */
  /* ============================================================ */

  /* ********************************************************************** */
  /* I. ushape_t */
  /* ********************************************************************** */

  /* see ushape_representation.c */

  ushape_t *ushape_alloc_internal (ushape_internal_t * pr,
                                   size_t intdim, size_t realdim);
  /* Builds a value with a graph and constraints empty (NULL) */

  ushape_array_t *
  ushape_of_pcons_array (ushape_internal_t * pr,
                         size_t intdim, size_t realdim,
                         pcons0_array_t * array);
  ushape_t *ushape_make (ushape_internal_t * pr, size_t code,
                         size_t datadim, size_t ptrdim);
  ushape_t *ushape_random (ushape_internal_t * pr, size_t size,
                           size_t datadim, size_t ptrdim);
  /* Builds radomly a ushape with size nodes and given variables. */

  ushape_array_t *ushape_of_formula (ushape_internal_t * pr,
                                     sh_formula_t* f, size_t* rsize);
  /* Builds a shape from an SL3 formula. */

  void ushape_set_bottom (ushape_internal_t * pr, ushape_t * a);
  /* put all fields to NULL after freeing */

  void ushape_free_internal (ushape_internal_t * pr, ushape_t * a);
  /* Free a ushape */

  char ushape_check (ushape_internal_t * pr, ushape_t * a);
  /* Check integrity of the canonical representation of ushapes */

  ushape_t *ushape_copy_internal (ushape_internal_t * pr, ushape_t * a);
  void ushape_copy_internal_scons (ushape_internal_t * pr, ushape_t * a, ushape_t * b);
  int ushape_cmp_internal (ushape_internal_t * pr, ushape_t * a, ushape_t * b);

  size_t ushape_get_size (ushape_t * a);
  ap_dimension_t ushape_dimension_var (ap_manager_t * man, ushape_t * a);
  ushape_array_t ushape_canonicalize_internal (ushape_internal_t * pr,
                                               ushape_t * a);

  /* see ushape_predicate.c */

  bool ushape_is_equal (ushape_t * a, ushape_t * b);
  bool ushape_is_lt (ushape_t * a, ushape_t * b);
  int ushape_cmp (ushape_t * a, ushape_t * b);


  /* ============================================================ */
  /* Meet and Join */
  /* ============================================================ */

  /* see ushape_nary.c */

  ushape_array_t *ushape_meet_pcons_array (ushape_internal_t * pr,
                                           bool destructive, ushape_t * a,
                                           pcons0_array_t * c);
  /* Version for meet returning a set of ushapes from shape constraints */

  /* ============================================================ */
  /* Assignement and Substitutions */
  /* ============================================================ */

  /* see ushape_transfer.c */

  ushape_array_t *ushape_assign_passign_array (ushape_internal_t * pr,
                                               bool destructive, ushape_t * a,
                                               passign0_array_t * op);
  /* Version for assign for shape assignment */

  ushape_t* ushape_strengthen (ushape_internal_t* pr,
                               ushape_t* a,
                               ap_dim_t* tdim,
                               size_t size);
  ushape_array_t *ushape_substitute_actuals (ushape_internal_t* pr,
                                             ushape_t* a,
                                             ap_dim_t* tdim,
                                             ap_dim_t* tnew,
                                             size_t size);
  /* Substitute actuals and strengthen */

  ushape_array_t *ushape_substitute_passign_array (ushape_internal_t * pr,
                                                   ushape_t * a, passign0_array_t * op);
  /* Pre of assignment */

  /* ============================================================ */
  /* Change and permutation of dimensions */
  /* ============================================================ */

  /* see ushape_resize.c */

  void
  ushape_apply_dimperm (ushape_internal_t * pr, ushape_t * a, ushape_t * b,
                        ap_dimperm_t * perm);
  void
  ushape_apply_dimperm_n (ushape_internal_t * pr, ushape_t * a, ushape_t * b,
                          ap_dimperm_t * perm);

  void
  ushape_apply_dimperm_dimchange (ushape_internal_t * pr, ushape_t * a,
                                  ushape_t * b, ap_dimperm_t * perm, size_t n);
  void
  ushape_apply_dimperm_dimchange_n (ushape_internal_t * pr, ushape_t * a,
                                    ushape_t * b, ap_dimperm_t * perm, ap_dimchange_t* dimchange);

  void
  ushape_expand_internal (ushape_internal_t * pr, ushape_t * a, ushape_t * b,
                          ap_dim_t dim, size_t n);

  void
  ushape_fold_internal (ushape_internal_t * pr, ushape_t * a, ushape_t * b,
                        ap_dim_array2_t * anon);


  /* ********************************************************************** */
  /* II. ushape_array_t */
  /* ********************************************************************** */

  /* see ushape_representation.c */

  void ushape_array_init (ushape_internal_t * pr, ushape_array_t * a,
                          size_t size);
  /* Allocate an array of SIZE  NULL pointers */

  ushape_array_t *ushape_array_make (ushape_internal_t * pr, size_t size);
  /* Allocate an array of SIZE  NULL pointers */

  int ushape_array_add (ushape_internal_t * pr, bool isset,
                        ushape_array_t * arr, size_t msize, bool docopy,
                        bool destructive, ushape_t * a);
  /* Add an element to the array after position msize */

  ushape_array_t *ushape_array_add_array (ushape_internal_t * pr, bool isset,
                                          ushape_array_t * a,
                                          ushape_array_t * b);
  /* Add to a the array b with property isset */

  int ushape_array_resize (ushape_internal_t * pr, ushape_array_t * a,
                           size_t size);
  /* Change size of a to size */

  void ushape_array_clear (ushape_internal_t * pr, ushape_array_t * a,
                           size_t size);
  /* Free the first size elements of the array and put them to NULL */

  ushape_array_t *ushape_array_copy (ushape_internal_t * pr,
                                     ushape_array_t * src, size_t size);
  /* Get a copy of the array for the first size elements */


  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */



#endif /* __USHAPE_INTERNAL_H_ */
