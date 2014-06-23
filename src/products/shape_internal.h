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


#ifndef __SHAPE_INTERNAL_H_
#define __SHAPE_INTERNAL_H_

#include "hgraph_fun.h"
#include "ushape_fun.h"
#include "ushape_internal.h"
#include "shape_fun.h"
#include "shape_manager.h"

/* *INDENT-OFF* */
#ifdef __cplusplus
extern "C"
{
#endif
  /* *INDENT-ON* */

  /* ********************************************************************** */
  /* Shapes */
  /* ********************************************************************** */

  /* ============================================================ */
  /* Internal Representation */
  /* ============================================================ */

    /*
   * Shapes is a powerset domain based on constrained shapes, so shapes
   * are collections of unit shapes.
   */

  struct _shape_t
  {
    ushape_array_t m; /* collection of size unit shapes */
    size_t msize; /* <= m.size, the number of NOT NULL elements
				 * in m invariant: the first msize elements
				 * of m are not NULL i.e., no gaps of NULLs
				 * inside m */
    bool set; /* true if the collection is a set */
    bool closed; /* true if the collection has only canonical
				 * elements */
    size_t intdim;
    size_t realdim;
  };


  /* ============================================================ */
  /* Internal Management */
  /* ============================================================ */

  /* see shape_representation.c */

  shape_t *shape_alloc_internal (shape_internal_t * pr, size_t size);
  shape_t *shape_make (shape_internal_t * pr, size_t code,
                       size_t datadim, size_t ptrdim);
  shape_t *shape_random (shape_internal_t * pr, size_t size,
                         size_t datadim, size_t ptrdim);
  /* Builds radomly a shape with size nodes and given variables. */

  void shape_free_internal (shape_internal_t * pr, shape_t * a);

  char shape_check (shape_internal_t * pr, shape_t * a);
  /* Check integrity of the canonical representation of ushapes */

  shape_t *shape_copy_internal (shape_internal_t * pr, shape_t * a);
  int shape_cmp_internal (shape_internal_t * pr, shape_t * a, shape_t * b);

  /* ============================================================ */
  /* Direct Management: Use with care! */
  /* ============================================================ */

  /* see shape_representation.c */

  /* ============================================================ */
  /* Meet and Join */
  /* ============================================================ */

  /* see shape_nary.c */

  shape_t *
  shape_meet_pcons_array (shape_internal_t * pr,
    bool destructive, shape_t * a,
    pcons0_array_t * c);
  /* Version for meet returning a shape from shape constraints */


  /* ============================================================ */
  /* Interface with sl3 */
  /* ============================================================ */

  /* see shape_nary.c */

  typedef struct shape_spec_t
  {
    int id; // unique id given from code, use to sort the list, -1 for dp input
    shape_t* shape_pos;
    shape_t* shape_neg;
    struct shape_spec_t* next;
    struct shape_spec_t* prev;
  } shape_spec_t;
  // stores specifications used during the analysis
  extern shape_spec_t* shape_specs;

  shape_t *shape_of_spec (shape_internal_t* pr, int spec,
                          bool version);
  /* Build shapes from sl3 spec, positive and negative versions. */

  /* ============================================================ */
  /* Closure Algorithms */
  /* ============================================================ */

  /* see shape_closure.c */

  bool shape_close (shape_internal_t * pr, shape_t * a);
  /* In place closure */

  /* ============================================================ */
  /* Printing internals */
  /* ============================================================ */

  void shape_fdump_internal (FILE * stream, shape_internal_t * pr, shape_t * a);
  void shape_fdump_le (FILE * stream, shape_internal_t * pr, shape_t * a1, shape_t * a2);

  /* *INDENT-OFF* */
#ifdef __cplusplus
}
#endif
/* *INDENT-ON* */


#endif /* __SHAPE_INTERNAL_H_ */
